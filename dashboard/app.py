"""
Prosculpt Dashboard - Flask backend.

Run with:
    python app.py
then open http://localhost:5000 in a browser.

The server needs read access to the log file and to the pipeline's output
directory - run it on the machine (e.g. the cluster login node) that can
see that filesystem, and access it from your laptop via an SSH tunnel if
needed, e.g.:
    ssh -L 5000:localhost:5000 you@cluster
"""

from __future__ import annotations

import csv
import io
import os
import zipfile

from flask import Flask, jsonify, request, send_file, abort
from werkzeug.exceptions import HTTPException

import parser as P
import run_targets as RT
from run_api import run_api

app = Flask(__name__, static_folder="static", template_folder="templates")
app.register_blueprint(run_api)


# Every route here is an API this dashboard's own frontend calls with
# fetch()/apiGet() and expects JSON back from - Flask's/Werkzeug's default
# abort() page is HTML, though, so every `abort(400, description=...)`
# elsewhere in this app (there are many, across app.py and run_api.py) was
# silently losing its specific `description` on the way to the browser:
# apiGet()'s `res.json()` parse just fails against an HTML body and falls
# back to the generic statusText ("Bad Request") instead. This turns any
# HTTPException (however it was raised) into a plain {"error": "..."} body
# with the same status code, so callers actually see the reason.
@app.errorhandler(HTTPException)
def handle_http_exception(e):
    response = jsonify({"error": e.description})
    response.status_code = e.code
    return response


# A network-mounted drive going into a bad state mid-request (this
# deployment has hit both a Windows "insufficient system resources" and a
# "volume does not contain a recognized file system" from a flaky
# SSHFS-Win/WinFsp mount - real incidents, not hypothetical) surfaces as a
# raw OSError from whatever os.*/open()/glob.glob() call happened to be
# touching it at the time - never something a client request itself did
# wrong. Without this handler, Flask's default handling turns that into a
# generic HTML 500 page: apiGet() in app.js can't parse HTML as JSON, so
# the frontend shows a bare "Internal Server Error" with none of the
# actual reason - visible only in this process's own console/log, which
# is how this class of failure has had to be diagnosed by hand so far.
# Registering this on OSError (not just for one endpoint) covers every
# route uniformly, since the trigger is never which endpoint was called,
# only whether the mount happened to be in a bad state at that moment.
@app.errorhandler(OSError)
def handle_os_error(e):
    # Not just f"{e}" - OSError.__str__() formats its filename with repr(),
    # which doubles up every backslash in a Windows path (a real path with
    # one '' between components shows up as '\' in the message). That's
    # meant for a Python traceback, not a path someone's about to copy out
    # of this error and paste into Explorer to go check on it themselves -
    # e.filename/e.filename2 are the plain, unescaped strings the OS call
    # actually failed on, so those are used directly instead.
    detail = e.strerror or str(e)
    for fname in (e.filename, e.filename2):
        if fname:
            detail = f"{detail}: {fname}"
    response = jsonify(
        {
            "error": f"Filesystem error ({e.__class__.__name__}): {detail}. This usually means "
            "a network-mounted drive this dashboard reads from is temporarily unreachable "
            "- try again in a moment; if it keeps happening, that drive's mount likely needs "
            "attention on the machine running this dashboard."
        }
    )
    response.status_code = 502
    return response

def _require_job_dir() -> str:
    """A tracked job is identified directly by its own base output
    directory now (the one containing logs/ and one numbered 01/02/...
    subdirectory per SLURM array task) - no more resolving it by parsing a
    log file (see get_job_status()/discover_tasks() in parser.py). This
    just reads the ?job_dir= param every route needs; it does not require
    the directory to exist yet (a freshly-staged job that hasn't started
    writing output is a valid, if empty, thing to track - get_job_status()
    itself reports that state rather than erroring)."""
    job_dir = request.args.get("job_dir")
    if not job_dir:
        abort(400, description="Missing ?job_dir= parameter")
    return job_dir


def _get_task_dir(job_dir: str, task_num_raw) -> str:
    """Resolves ?task=N to that task's own directory under job_dir - looked
    up via discover_tasks() (not a naive job_dir/f"{n:02d}" join) so this
    still works regardless of exactly how many digits the folder name
    uses."""
    if task_num_raw is None:
        abort(400, description="Missing ?task= parameter")
    try:
        task_num = int(task_num_raw)
    except ValueError:
        abort(400, description="?task= must be an integer")
    for num, task_dir in P.discover_tasks(job_dir):
        if num == task_num:
            return task_dir
    abort(404, description=f"Task {task_num} not found under {job_dir}")


def _is_within(candidate: str, base: str) -> bool:
    # realpath() resolves symlinks/bind-mounts, so two different-looking
    # absolute paths that point at the same real location compare equal -
    # HPC storage very commonly exposes one directory under more than one
    # absolute path (e.g. a project directory reachable both at its "real"
    # location and via a symlink into a user's home folder), and `1_rfdiff/`
    # in particular is a plausible candidate for living on a different,
    # relocated/archived storage tier than the rest of output_dir. A purely
    # lexical prefix check rejects a perfectly legitimate file the moment it
    # arrives expressed via a different alias than the one `base` resolved
    # to - realpath() doesn't require the path to exist to be resolved, so
    # this is safe to call even for a path that turns out not to be there.
    #
    # Falls back to a purely lexical comparison (normpath, no symlink
    # resolution) if realpath() itself raises - not a hypothetical: on a
    # WinFsp-mounted network drive (SSHFS-Win, this deployment's actual
    # setup), realpath()'s underlying GetFinalPathNameByHandle call is a
    # known-unreliable Win32 API (see winfsp/winfsp#427, winfsp/sshfs-win#243),
    # raising WinError 1005 "volume does not contain a recognized file
    # system" even though the exact same path opens/reads/lists fine
    # through every other API - not a sign the path is actually
    # unreachable. Losing symlink-awareness in that fallback is an
    # acceptable trade: the property this function exists for (candidate
    # can't escape outside base) still holds via a plain prefix check, it
    # just won't additionally catch a symlink-based escape on a filesystem
    # where realpath() itself can't be trusted to work at all right now.
    # Both paths are re-derived the same way on any failure (not just
    # whichever call happened to raise) so they're never compared in a
    # part-resolved, part-lexical mismatched state.
    try:
        base_real = os.path.realpath(base)
        full_real = os.path.realpath(candidate)
    except OSError:
        base_real = os.path.normpath(base)
        full_real = os.path.normpath(candidate)
    return full_real == base_real or full_real.startswith(base_real + os.sep)


def _safe_join(base: str, *parts: str) -> str:
    """Join and ensure the result stays inside `base` (avoid path traversal
    from the file-browser / pdb-fetch endpoints)."""
    full = os.path.normpath(os.path.join(base, *parts))
    if not _is_within(full, base):
        abort(400, description="Invalid path")
    return full


# ---------------------------------------------------------------------------
# Pages
# ---------------------------------------------------------------------------


@app.route("/")
def index():
    from flask import render_template

    return render_template("index.html")


# ---------------------------------------------------------------------------
# Status / overview
# ---------------------------------------------------------------------------


@app.route("/api/status")
def api_status():
    job_dir = _require_job_dir()
    return jsonify(P.get_job_status(job_dir))


# ---------------------------------------------------------------------------
# Backbones (rfdiffusion)
# ---------------------------------------------------------------------------


@app.route("/api/backbones")
def api_backbones():
    """Every task's backbones, concatenated and tagged with task_num - for
    the common num_tasks=1 case this is exactly one task's list (identical
    shape to before, plus the tag); for num_tasks>1 the Backbones tab shows
    every task's backbones together rather than needing a separate tracked
    job per task."""
    job_dir = _require_job_dir()
    job_status = P.get_job_status(job_dir)
    if job_status.get("error"):
        abort(400, description=job_status["error"])
    out = []
    for task in job_status["tasks"]:
        stage = task["stage"]
        filtering_done = stage not in ("setup", "rfdiffusion", "filtering")
        for b in P.list_backbones(task["task_dir"], filtering_done=filtering_done):
            b["task_num"] = task["task_num"]
            out.append(b)
    return jsonify(out)


@app.route("/api/backbone_pdb")
def api_backbone_pdb():
    job_dir = _require_job_dir()
    task_dir = _get_task_dir(job_dir, request.args.get("task"))
    name = request.args.get("name")
    status = request.args.get("status", "passed")
    subdir = "failed_filters" if status == "failed_filter" else ""
    path = _safe_join(task_dir, "1_rfdiff", subdir, f"{name}.pdb")
    if not os.path.isfile(path):
        abort(404, description="pdb not found")
    return send_file(path, mimetype="chemical/x-pdb")


# ---------------------------------------------------------------------------
# Sequences (mpnn)
# ---------------------------------------------------------------------------


@app.route("/api/sequences")
def api_sequences():
    """Every task's sequences, concatenated - each backbone/monomer entry
    tagged with task_num so num_tasks>1 shows all tasks together (see
    /api/backbones's docstring for the same pattern)."""
    job_dir = _require_job_dir()
    out = {"backbones": [], "monomers": []}
    for task_num, task_dir in P.discover_tasks(job_dir):
        mp = P.list_sequences(task_dir)
        for key in ("backbones", "monomers"):
            for entry in mp[key]:
                entry["task_num"] = task_num
                out[key].append(entry)
    return jsonify(out)


# ---------------------------------------------------------------------------
# Models (modeling stage)
# ---------------------------------------------------------------------------


@app.route("/api/models")
def api_models():
    """Every task's models, concatenated and tagged with task_num (same
    pattern as /api/backbones/api/sequences above)."""
    job_dir = _require_job_dir()
    out = []
    for task_num, task_dir in P.discover_tasks(job_dir):
        cfg = P.extract_config(None, task_dir)
        for m in P.list_models(task_dir, model_monomer=bool(cfg.get("model_monomer"))):
            m["task_num"] = task_num
            out.append(m)
    return jsonify(out)


def _resolve_within_output_dir(output_dir: str, path: str) -> str:
    """`path` is normally an absolute path already under output_dir (as
    returned by list_models/list_backbones); re-validate it stays inside
    output_dir regardless of whether it arrived absolute or relative."""
    rel = os.path.relpath(path, output_dir) if os.path.isabs(path) else path
    return _safe_join(output_dir, rel)


def _resolve_within_output_dir_or_none(output_dir: str, path: str):
    """Same as _resolve_within_output_dir but returns None instead of
    aborting the request - for callers (like the zip export) that want to
    just skip a bad/missing entry rather than fail the whole response."""
    try:
        rel = os.path.relpath(path, output_dir) if os.path.isabs(path) else path
        full = os.path.normpath(os.path.join(output_dir, rel))
    except Exception:
        return None
    return full if _is_within(full, output_dir) else None


@app.route("/api/model_pdb")
def api_model_pdb():
    job_dir = _require_job_dir()
    task_dir = _get_task_dir(job_dir, request.args.get("task"))
    path = request.args.get("path")
    full = _resolve_within_output_dir(task_dir, path)
    if not os.path.isfile(full):
        abort(404, description="structure file not found")
    mimetype = "chemical/x-cif" if full.lower().endswith(".cif") else "chemical/x-pdb"
    return send_file(full, mimetype=mimetype)


@app.route("/api/model_confidence")
def api_model_confidence():
    job_dir = _require_job_dir()
    task_dir = _get_task_dir(job_dir, request.args.get("task"))
    path = request.args.get("path")
    full = _resolve_within_output_dir(task_dir, path)
    if not os.path.isfile(full):
        abort(404, description="confidence file not found")
    return jsonify(P.load_confidence(full))


@app.route("/api/trb")
def api_trb():
    """Generic .trb reader - takes any path under a task's own directory,
    so it works for a backbone's own _N.trb (Backbones tab) and for the
    RFdiffusion .trb reachable from a final result row's path_rfdiff column
    (Results tab, path derived client-side by swapping .pdb -> .trb)."""
    job_dir = _require_job_dir()
    task_dir = _get_task_dir(job_dir, request.args.get("task"))
    path = request.args.get("path")
    # A Backbones/Models-tab trb_path is already local (built server-side
    # from the already-translated task dir), so try it exactly as given
    # first. Only fall back to translate_remote_path() - for the Results
    # tab's path_rfdiff, read verbatim from final_output.csv (see
    # load_final_csv()'s docstring) as recorded on whichever cluster the
    # job ran on - if the direct interpretation doesn't resolve. Applying
    # the translation unconditionally used to be assumed a safe no-op for
    # the already-local case, but it isn't: when this dashboard runs
    # directly on the target host (not on a machine that mounts it, e.g.
    # via the Q: SSHFS mount), an already-correct local path can itself
    # match a mounts: entry's "remote" prefix and get rewritten into that
    # mount's "local" side (a Windows path, on Linux) - a real bug hit in
    # production, not just a hypothetical.
    full = _resolve_within_output_dir_or_none(task_dir, path)
    if not full or not os.path.isfile(full):
        full = _resolve_within_output_dir(task_dir, RT.translate_remote_path(path))
    if not os.path.isfile(full):
        abort(404, description="trb file not found")
    return jsonify(P.load_trb_provenance(full))


# ---------------------------------------------------------------------------
# Results (final)
# ---------------------------------------------------------------------------


def _load_merged_final_csv(job_dir: str):
    """Merges every task's own final_output.csv into one table - union of
    columns across tasks (a task missing a given column just gets a blank
    for it), plus a leading "task" column so every row stays attributable
    - the same idea api_export_filtered_multi() already uses to merge
    across several *tracked jobs*, just one level down, across one job's
    own tasks. None if no task has a final_output.csv yet."""
    all_columns: list[str] = []
    per_task = []
    for task_num, task_dir in P.discover_tasks(job_dir):
        result = P.load_final_csv(task_dir)
        if result is None or not result["columns"]:
            continue
        columns, rows = result["columns"], result["rows"]
        for c in columns:
            if c not in all_columns:
                all_columns.append(c)
        per_task.append((task_num, columns, rows))
    if not per_task:
        return None
    merged_rows = []
    for task_num, columns, rows in per_task:
        col_index = {c: i for i, c in enumerate(columns)}
        for r in rows:
            merged_rows.append(
                [str(task_num), *(r[col_index[c]] if c in col_index and col_index[c] < len(r) else "" for c in all_columns)]
            )
    return {"columns": ["task", *all_columns], "rows": merged_rows}


@app.route("/api/results")
def api_results():
    """Aggregated across every task - final_pdbs entries are tagged with
    task_num (a plain name string wouldn't be unique across tasks)."""
    job_dir = _require_job_dir()
    final_pdbs = []
    any_csv = False
    all_finished = True
    any_task = False
    for task_num, task_dir in P.discover_tasks(job_dir):
        any_task = True
        summary = P.results_summary(task_dir)
        for name in summary["final_pdbs"]:
            final_pdbs.append({"name": name, "task_num": task_num})
        any_csv = any_csv or summary["csv_exists"]
        all_finished = all_finished and summary["finished"]
    return jsonify({"final_pdbs": final_pdbs, "csv_exists": any_csv, "finished": any_task and all_finished})


@app.route("/api/final_pdb")
def api_final_pdb():
    job_dir = _require_job_dir()
    task_dir = _get_task_dir(job_dir, request.args.get("task"))
    name = request.args.get("name")
    path = _safe_join(task_dir, "final_pdbs", name)
    if not os.path.isfile(path):
        abort(404, description="pdb not found")
    return send_file(path, mimetype="chemical/x-pdb")


@app.route("/api/final_csv")
def api_final_csv():
    job_dir = _require_job_dir()
    result = _load_merged_final_csv(job_dir)
    if result is None:
        abort(404, description="final_output.csv not found")
    return jsonify(result)


@app.route("/api/export_filtered", methods=["POST"])
def api_export_filtered():
    """Zips up the currently-filtered final_output.csv rows plus the
    model_path pdb for each of them. (The dashboard's Export button calls
    this a ".zip" - true .rar creation needs a proprietary external binary
    that isn't reliably available on a cluster, so .zip is used instead;
    it needs no extra dependency and every OS can open it natively.)"""
    job_dir = request.args.get("job_dir")
    data = request.get_json(silent=True) or {}
    job_dir = data.get("job_dir") or job_dir
    row_ids = data.get("row_ids")
    if not job_dir:
        abort(400, description="Missing job_dir")

    result = _load_merged_final_csv(job_dir)
    if result is None:
        abort(404, description="final_output.csv not found")
    columns, data_rows = result["columns"], result["rows"]

    if row_ids is not None:
        wanted = {int(i) for i in row_ids}
        selected_rows = [r for i, r in enumerate(data_rows) if i in wanted]
    else:
        selected_rows = data_rows

    if not selected_rows:
        abort(400, description="No rows selected for export")

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        csv_buf = io.StringIO()
        writer = csv.writer(csv_buf)
        writer.writerow(columns)
        writer.writerows(selected_rows)
        zf.writestr("filtered_output.csv", csv_buf.getvalue())

        # _load_merged_final_csv() always puts "task" first - which task a
        # given row's model_path lives under, needed since the same model
        # filename can exist in more than one task's own final_pdbs/.
        if "model_path" in columns and "task" in columns:
            model_idx = columns.index("model_path")
            task_idx = columns.index("task")
            # Only used to actually locate+read the file for zipping - the
            # CSV text written above (writer.writerows(selected_rows)) is
            # already done and keeps the original, untranslated values, the
            # same as final_output.csv itself (see load_final_csv()'s
            # docstring for why that's deliberate).
            translate = RT.make_path_translator()
            seen = set()
            for r in selected_rows:
                p = r[model_idx] if model_idx < len(r) else ""
                task_str = r[task_idx] if task_idx < len(r) else ""
                if not p or (task_str, p) in seen:
                    continue
                seen.add((task_str, p))
                try:
                    task_dir = _get_task_dir(job_dir, task_str)
                except HTTPException:
                    continue
                full = _resolve_within_output_dir_or_none(task_dir, translate(p))
                if full and os.path.isfile(full):
                    zf.write(full, arcname=os.path.join("models", f"task_{task_str}", os.path.basename(full)))

    buf.seek(0)
    return send_file(
        buf,
        mimetype="application/zip",
        as_attachment=True,
        download_name="prosculpt_filtered_export.zip",
    )


def _job_label(job_dir: str) -> str:
    """Same derivation the frontend uses for a job's short display name
    (see jobLabel() in app.js) - kept in lockstep so the "source_job"
    column and the zip's per-job model subfolders match what the user
    sees on screen."""
    return os.path.basename(os.path.normpath(job_dir))


@app.route("/api/export_filtered_multi", methods=["POST"])
def api_export_filtered_multi():
    """Same idea as /api/export_filtered, but across every tracked job at
    once (the "All jobs results" tab's Export button): one combined CSV
    (columns are the union across jobs - a job missing a given column
    just gets a blank for it - plus a "source_job" column), and each
    job's model pdbs kept in their own subfolder under models/ so two
    jobs' identically-named model files can't collide in the zip."""
    data = request.get_json(silent=True) or {}
    jobs = data.get("jobs")
    if not jobs:
        abort(400, description="Missing jobs")

    all_columns: list[str] = []
    per_job_selected: list[tuple[str, str, list[str], list[list[str]]]] = []  # (job_dir, label, columns, rows)

    for job in jobs:
        job_dir = job.get("job_dir")
        row_ids = job.get("row_ids")
        if not job_dir:
            continue
        result = _load_merged_final_csv(job_dir)
        if result is None:
            continue
        columns, data_rows = result["columns"], result["rows"]
        if row_ids is not None:
            wanted = {int(i) for i in row_ids}
            selected_rows = [r for i, r in enumerate(data_rows) if i in wanted]
        else:
            selected_rows = data_rows
        if not selected_rows:
            continue
        for c in columns:
            if c not in all_columns:
                all_columns.append(c)
        per_job_selected.append((job_dir, _job_label(job_dir), columns, selected_rows))

    if not per_job_selected:
        abort(400, description="No rows selected for export")

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        csv_buf = io.StringIO()
        writer = csv.writer(csv_buf)
        writer.writerow(["source_job", *all_columns])
        for job_dir, label, columns, selected_rows in per_job_selected:
            col_index = {c: i for i, c in enumerate(columns)}
            for r in selected_rows:
                writer.writerow([label, *(r[col_index[c]] if c in col_index and col_index[c] < len(r) else "" for c in all_columns)])
        zf.writestr("filtered_output.csv", csv_buf.getvalue())

        # One translator, built once and reused across every job in this
        # export - global (spans all configured targets), not job-specific,
        # so there's no need to rebuild it per job. Only used to actually
        # locate+read each file for zipping; the CSV text written above
        # keeps the original, untranslated values (see load_final_csv()'s
        # docstring for why that's deliberate).
        translate = RT.make_path_translator()
        for job_dir, label, columns, selected_rows in per_job_selected:
            if "model_path" not in columns or "task" not in columns:
                continue
            model_idx = columns.index("model_path")
            task_idx = columns.index("task")
            seen = set()
            for r in selected_rows:
                p = r[model_idx] if model_idx < len(r) else ""
                task_str = r[task_idx] if task_idx < len(r) else ""
                if not p or (task_str, p) in seen:
                    continue
                seen.add((task_str, p))
                try:
                    task_dir = _get_task_dir(job_dir, task_str)
                except HTTPException:
                    continue
                full = _resolve_within_output_dir_or_none(task_dir, translate(p))
                if full and os.path.isfile(full):
                    zf.write(full, arcname=os.path.join("models", label, os.path.basename(full)))

    buf.seek(0)
    return send_file(
        buf,
        mimetype="application/zip",
        as_attachment=True,
        download_name="prosculpt_all_jobs_filtered_export.zip",
    )


# ---------------------------------------------------------------------------
# Crash / error log
# ---------------------------------------------------------------------------


@app.route("/api/error_log")
def api_error_log():
    """Returns one task's slurm .err file content whenever it exists,
    regardless of whether a crash/cancellation has actually been detected -
    unconditionally looked up (see find_task_log()), not only once
    something looks wrong. Unlike before, a task may simply have no log at
    all (never captured, or the job was submitted outside this dashboard) -
    that's reported as err_exists: False, log_found: False, not an error."""
    job_dir = _require_job_dir()
    task_dir = _get_task_dir(job_dir, request.args.get("task"))
    task_num = int(request.args.get("task"))
    log_info = P.find_task_log(job_dir, task_num)
    err_path = log_info["err_path"]
    if not err_path or not os.path.isfile(err_path):
        return jsonify({"err_exists": False, "err_path": err_path, "log_found": log_info["out_path"] is not None})
    err = P.read_text_file(err_path)
    return jsonify({"err_exists": True, "err_path": err_path, "log_found": True, **err})


@app.route("/api/output_log")
def api_output_log():
    """Returns one task's own .out log content, for the Output log tab.
    Unlike before, this may genuinely not exist (a task tracked without a
    log at all) - reported as log_found: False, not a 404, since that's
    an expected state now, not a broken request."""
    job_dir = _require_job_dir()
    task_dir = _get_task_dir(job_dir, request.args.get("task"))
    task_num = int(request.args.get("task"))
    log_info = P.find_task_log(job_dir, task_num)
    out_path = log_info["out_path"]
    if not out_path or not os.path.isfile(out_path):
        return jsonify({"log_found": False})
    return jsonify({"log_found": True, **P.read_text_file(out_path)})


# ---------------------------------------------------------------------------
# File browser (helper for picking the log file in the UI)
# ---------------------------------------------------------------------------


@app.route("/api/browse_roots")
def api_browse_roots():
    """Every configured target's best locally-reachable root (see
    list_browse_roots() in run_targets.py) - powers the file browser's
    "Jump to..." dropdown, so a multi-cluster setup can switch between
    clusters without typing/pasting each one's root path by hand."""
    return jsonify(RT.list_browse_roots())


@app.route("/api/browse")
def api_browse():
    # No ?path= at all means "just opened the browser" (the frontend omits
    # it for exactly that case, see openBrowse() in app.js) - start at the
    # locally-mounted projects root from dashboard_config.yaml if one's
    # configured (where real jobs actually are, unlike the process's own
    # home directory), falling back to the home directory the way this
    # always worked before that file existed.
    path = request.args.get("path") or RT.get_default_browse_root() or os.path.expanduser("~")
    path = os.path.abspath(path)
    # A pasted path that turns out to be a file (rather than a directory -
    # see the "paste a path" input, goToBrowsePath() in app.js) is resolved
    # to its containing folder - there's no equivalent of the old
    # file-preselect here, since a job is tracked by directory now, not by
    # picking one specific file.
    if os.path.isfile(path):
        path = os.path.dirname(path)
    if not os.path.isdir(path):
        abort(400, description="Not a directory")
    entries = []
    try:
        # Directories only - a job is tracked by pointing at its own output
        # directory now, so a plain file (whatever it is) is never a valid
        # target to select here; only worth listing at all.
        for name in sorted(os.listdir(path)):
            full = os.path.join(path, name)
            if os.path.isdir(full):
                entries.append({"name": name, "path": full, "is_dir": True})
    except PermissionError:
        abort(403, description="Permission denied")
    parent = os.path.dirname(path) if path != os.path.dirname(path) else None
    return jsonify({"path": path, "parent": parent, "entries": entries})


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    debug = os.environ.get("FLASK_DEBUG", "0") == "1"
    # threaded=True: Werkzeug's dev server otherwise handles exactly one
    # request at a time - with several people pointed at the same
    # dashboard, one request that's slow for any reason (a stalled network
    # mount under a glob(), a big final_output.csv, a wedged console write
    # to this very terminal - see dashboard/README.md's Troubleshooting
    # section) blocks every other user's request too, not just the one
    # that triggered it, since there's no second worker to pick up new
    # ones. One thread per request means a single slow/stuck request only
    # ever stalls itself.
    app.run(host="0.0.0.0", port=port, debug=debug, threaded=True)
