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

import parser as P

app = Flask(__name__, static_folder="static", template_folder="templates")

# Simple in-memory cache of the last resolved output_dir per log path, so
# helper endpoints (pdb/json fetchers) don't need to re-parse the whole log
# on every click.
_output_dir_cache: dict[str, str] = {}


def _get_output_dir(log_path: str) -> str:
    if log_path in _output_dir_cache:
        return _output_dir_cache[log_path]
    loc = P.resolve_output_dir(log_path)
    if not loc.get("output_dir"):
        abort(400, description=loc.get("error") or "Could not resolve output_dir")
    _output_dir_cache[log_path] = loc["output_dir"]
    return loc["output_dir"]


def _safe_join(base: str, *parts: str) -> str:
    """Join and ensure the result stays inside `base` (avoid path traversal
    from the file-browser / pdb-fetch endpoints)."""
    full = os.path.normpath(os.path.join(base, *parts))
    base_norm = os.path.normpath(base)
    if not (full == base_norm or full.startswith(base_norm + os.sep)):
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
    log_path = request.args.get("log")
    if not log_path:
        return jsonify({"error": "Missing ?log= parameter"}), 400
    status = P.get_status(log_path)
    if status.get("location", {}).get("output_dir"):
        _output_dir_cache[log_path] = status["location"]["output_dir"]
    return jsonify(status)


# ---------------------------------------------------------------------------
# Backbones (rfdiffusion)
# ---------------------------------------------------------------------------


@app.route("/api/backbones")
def api_backbones():
    log_path = request.args.get("log")
    output_dir = _get_output_dir(log_path)
    with open(log_path, "r", errors="replace") as f:
        lines = f.readlines()
    stage = P.detect_stage(lines)["stage"]
    filtering_done = stage not in ("setup", "rfdiffusion", "filtering")
    return jsonify(P.list_backbones(output_dir, filtering_done=filtering_done))


@app.route("/api/backbone_pdb")
def api_backbone_pdb():
    log_path = request.args.get("log")
    name = request.args.get("name")
    status = request.args.get("status", "passed")
    output_dir = _get_output_dir(log_path)
    subdir = "failed_filters" if status == "failed_filter" else ""
    path = _safe_join(output_dir, "1_rfdiff", subdir, f"{name}.pdb")
    if not os.path.isfile(path):
        abort(404, description="pdb not found")
    return send_file(path, mimetype="chemical/x-pdb")


# ---------------------------------------------------------------------------
# Sequences (mpnn)
# ---------------------------------------------------------------------------


@app.route("/api/sequences")
def api_sequences():
    log_path = request.args.get("log")
    output_dir = _get_output_dir(log_path)
    return jsonify(P.list_sequences(output_dir))


# ---------------------------------------------------------------------------
# Models (modeling stage)
# ---------------------------------------------------------------------------


@app.route("/api/models")
def api_models():
    log_path = request.args.get("log")
    output_dir = _get_output_dir(log_path)
    cfg = P.extract_config(log_path, output_dir)
    return jsonify(P.list_models(output_dir, model_monomer=bool(cfg.get("model_monomer"))))


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
    base_norm = os.path.normpath(output_dir)
    if full == base_norm or full.startswith(base_norm + os.sep):
        return full
    return None


@app.route("/api/model_pdb")
def api_model_pdb():
    log_path = request.args.get("log")
    path = request.args.get("path")
    output_dir = _get_output_dir(log_path)
    full = _resolve_within_output_dir(output_dir, path)
    if not os.path.isfile(full):
        abort(404, description="structure file not found")
    mimetype = "chemical/x-cif" if full.lower().endswith(".cif") else "chemical/x-pdb"
    return send_file(full, mimetype=mimetype)


@app.route("/api/model_confidence")
def api_model_confidence():
    log_path = request.args.get("log")
    path = request.args.get("path")
    output_dir = _get_output_dir(log_path)
    full = _resolve_within_output_dir(output_dir, path)
    if not os.path.isfile(full):
        abort(404, description="confidence file not found")
    return jsonify(P.load_confidence(full))


@app.route("/api/trb")
def api_trb():
    """Generic .trb reader - takes any path under output_dir, so it works
    for a backbone's own _N.trb (Backbones tab) and for the RFdiffusion
    .trb reachable from a final result row's path_rfdiff column (Results
    tab, path derived client-side by swapping .pdb -> .trb)."""
    log_path = request.args.get("log")
    path = request.args.get("path")
    output_dir = _get_output_dir(log_path)
    full = _resolve_within_output_dir(output_dir, path)
    if not os.path.isfile(full):
        abort(404, description="trb file not found")
    return jsonify(P.load_trb_provenance(full))


# ---------------------------------------------------------------------------
# Results (final)
# ---------------------------------------------------------------------------


@app.route("/api/results")
def api_results():
    log_path = request.args.get("log")
    output_dir = _get_output_dir(log_path)
    return jsonify(P.results_summary(output_dir))


@app.route("/api/final_pdb")
def api_final_pdb():
    log_path = request.args.get("log")
    name = request.args.get("name")
    output_dir = _get_output_dir(log_path)
    path = _safe_join(output_dir, "final_pdbs", name)
    if not os.path.isfile(path):
        abort(404, description="pdb not found")
    return send_file(path, mimetype="chemical/x-pdb")


@app.route("/api/final_csv")
def api_final_csv():
    log_path = request.args.get("log")
    output_dir = _get_output_dir(log_path)
    csv_path = os.path.join(output_dir, "final_output.csv")
    if not os.path.isfile(csv_path):
        abort(404, description="final_output.csv not found")
    with open(csv_path, newline="", errors="replace") as f:
        reader = csv.reader(f)
        rows = list(reader)
    if not rows:
        return jsonify({"columns": [], "rows": []})
    columns, data_rows = rows[0], rows[1:]
    return jsonify({"columns": columns, "rows": data_rows})


@app.route("/api/export_filtered", methods=["POST"])
def api_export_filtered():
    """Zips up the currently-filtered final_output.csv rows plus the
    model_path pdb for each of them. (The dashboard's Export button calls
    this a ".zip" - true .rar creation needs a proprietary external binary
    that isn't reliably available on a cluster, so .zip is used instead;
    it needs no extra dependency and every OS can open it natively.)"""
    log_path = request.args.get("log")
    data = request.get_json(silent=True) or {}
    log_path = data.get("log") or log_path
    row_ids = data.get("row_ids")
    if not log_path:
        abort(400, description="Missing log")

    output_dir = _get_output_dir(log_path)
    csv_path = os.path.join(output_dir, "final_output.csv")
    if not os.path.isfile(csv_path):
        abort(404, description="final_output.csv not found")

    with open(csv_path, newline="", errors="replace") as f:
        reader = csv.reader(f)
        rows = list(reader)
    if not rows:
        abort(404, description="final_output.csv is empty")
    columns, data_rows = rows[0], rows[1:]

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

        if "model_path" in columns:
            idx = columns.index("model_path")
            seen = set()
            for r in selected_rows:
                p = r[idx] if idx < len(r) else ""
                if not p or p in seen:
                    continue
                seen.add(p)
                full = _resolve_within_output_dir_or_none(output_dir, p)
                if full and os.path.isfile(full):
                    zf.write(full, arcname=os.path.join("models", os.path.basename(full)))

    buf.seek(0)
    return send_file(
        buf,
        mimetype="application/zip",
        as_attachment=True,
        download_name="prosculpt_filtered_export.zip",
    )


def _job_label(log_path: str) -> str:
    """Same derivation the frontend uses for a job's short display name
    (see jobLabel() in app.js) - kept in lockstep so the "source_job"
    column and the zip's per-job model subfolders match what the user
    sees on screen."""
    base = os.path.basename(log_path)
    return os.path.splitext(base)[0]


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
    per_job_selected: list[tuple[str, str, list[str], list[list[str]]]] = []  # (log_path, label, columns, rows)

    for job in jobs:
        log_path = job.get("log")
        row_ids = job.get("row_ids")
        if not log_path:
            continue
        output_dir = _get_output_dir(log_path)
        csv_path = os.path.join(output_dir, "final_output.csv")
        if not os.path.isfile(csv_path):
            continue
        with open(csv_path, newline="", errors="replace") as f:
            reader = csv.reader(f)
            rows = list(reader)
        if not rows:
            continue
        columns, data_rows = rows[0], rows[1:]
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
        per_job_selected.append((log_path, _job_label(log_path), columns, selected_rows))

    if not per_job_selected:
        abort(400, description="No rows selected for export")

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        csv_buf = io.StringIO()
        writer = csv.writer(csv_buf)
        writer.writerow(["source_job", *all_columns])
        for log_path, label, columns, selected_rows in per_job_selected:
            col_index = {c: i for i, c in enumerate(columns)}
            for r in selected_rows:
                writer.writerow([label, *(r[col_index[c]] if c in col_index and col_index[c] < len(r) else "" for c in all_columns)])
        zf.writestr("filtered_output.csv", csv_buf.getvalue())

        for log_path, label, columns, selected_rows in per_job_selected:
            if "model_path" not in columns:
                continue
            output_dir = _get_output_dir(log_path)
            idx = columns.index("model_path")
            seen = set()
            for r in selected_rows:
                p = r[idx] if idx < len(r) else ""
                if not p or p in seen:
                    continue
                seen.add(p)
                full = _resolve_within_output_dir_or_none(output_dir, p)
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
    """Returns the slurm .err file's content whenever it exists, regardless
    of whether a crash/cancellation has actually been detected - the
    dashboard looks this file up unconditionally (see locate_err_file), not
    only once something looks wrong."""
    log_path = request.args.get("log")
    if not log_path:
        return jsonify({"error": "Missing ?log= parameter"}), 400
    err_loc = P.locate_err_file(log_path)
    if not err_loc["err_exists"]:
        return jsonify({"err_exists": False, "err_path": err_loc["err_path"]})
    err = P.read_text_file(err_loc["err_path"])
    return jsonify({"err_exists": True, "err_path": err_loc["err_path"], **err})


@app.route("/api/output_log")
def api_output_log():
    """Returns the job's own .out log content, for the Output log tab -
    always available (unlike the .err file) since it's the same file the
    dashboard already requires to operate."""
    log_path = request.args.get("log")
    if not log_path:
        return jsonify({"error": "Missing ?log= parameter"}), 400
    if not os.path.isfile(log_path):
        abort(404, description="log file not found")
    return jsonify(P.read_text_file(log_path))


# ---------------------------------------------------------------------------
# File browser (helper for picking the log file in the UI)
# ---------------------------------------------------------------------------


@app.route("/api/browse")
def api_browse():
    path = request.args.get("path") or os.path.expanduser("~")
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        abort(400, description="Not a directory")
    entries = []
    try:
        for name in sorted(os.listdir(path)):
            full = os.path.join(path, name)
            is_dir = os.path.isdir(full)
            if not is_dir and not (name.endswith(".out") or name.endswith(".log") or name.endswith(".txt")):
                continue
            entries.append({"name": name, "path": full, "is_dir": is_dir})
    except PermissionError:
        abort(403, description="Permission denied")
    parent = os.path.dirname(path) if path != os.path.dirname(path) else None
    return jsonify({"path": path, "parent": parent, "entries": entries})


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    debug = os.environ.get("FLASK_DEBUG", "0") == "1"
    app.run(host="0.0.0.0", port=port, debug=debug)
