"""
Run-job API routes - submitting new jobs (built via the GUI or pointed at
an already-prepared project directory), checking `squeue`, and polling for
a just-submitted job's log file so it can be handed off to the tracking
tab. Kept as its own blueprint (registered from app.py) rather than piled
into app.py directly, which is already sizable and entirely read-only
today - this is the one place that writes files and runs subprocesses.
"""

from __future__ import annotations

import glob
import os
import re
import shutil
import tempfile

from flask import Blueprint, abort, jsonify, request
from werkzeug.utils import secure_filename

import job_staging as JS
import run_targets as RT

run_api = Blueprint("run_api", __name__)

# Matches sbatch's own "Submitted batch job <id>" line (still printed to
# stderr/stdout by sbatch itself even when slurm_runner.py additionally
# captures --parsable output) as a fallback if the --parsable job id line
# slurm_runner.py prints isn't found for some reason.
_JOB_ID_LINE_RE = re.compile(r"submitted to slurm with id (\d+)", re.IGNORECASE)
_FALLBACK_JOB_ID_RE = re.compile(r"Submitted batch job (\d+)")


def _target_or_400(name):
    try:
        return RT.get_target(name)
    except RT.RunTargetError as e:
        abort(400, description=str(e))


@run_api.route("/api/run/targets")
def api_run_targets():
    # get_defaults() never raises (see its own docstring) - dashboard-wide
    # defaults should still reach the frontend even if no target is
    # configured yet, so it's fetched unconditionally rather than only in
    # the success branch below.
    defaults = RT.get_defaults()
    try:
        return jsonify({"targets": RT.list_targets(), "defaults": defaults})
    except RT.RunTargetError as e:
        # Not configured yet is a normal, expected state (e.g. right after
        # a fresh checkout) - report it as data, not a 500.
        return jsonify({"error": str(e), "targets": [], "defaults": defaults})


@run_api.route("/api/run/squeue", methods=["POST"])
def api_run_squeue():
    body = request.get_json(silent=True) or {}
    target = _target_or_400(body.get("target"))
    try:
        result = RT.run_remote_or_local(target, ["squeue", "--me"])
    except RT.RunTargetError as e:
        abort(502, description=str(e))
    return jsonify(result)


def _collect_uploaded_files() -> dict[str, bytes]:
    """Reconstructs the {relative_path: content} file map the frontend
    sends as two parallel multipart fields: repeated `file_paths` (text,
    the path this file should land at inside the job directory - decided
    client-side, since that's also where the config YAML referencing those
    same paths gets built) and repeated `file_contents` (the actual
    files), in matching order."""
    paths = request.form.getlist("file_paths")
    uploads = request.files.getlist("file_contents")
    if len(paths) != len(uploads):
        abort(400, description="file_paths and file_contents must have the same length")
    files: dict[str, bytes] = {}
    for rel_path, upload in zip(paths, uploads):
        # secure_filename() strips directory components too aggressively
        # for paths like "alignments/x_Chain_B.a3m" - sanitize each
        # segment individually instead, then let job_staging's own
        # containment check catch anything that still tries to escape.
        safe_parts = [secure_filename(part) for part in rel_path.split("/") if part not in ("", ".", "..")]
        if not safe_parts:
            abort(400, description=f"Invalid file path from client: {rel_path!r}")
        files["/".join(safe_parts)] = upload.read()
    return files


@run_api.route("/api/run/submit", methods=["POST"])
def api_run_submit():
    mode = request.form.get("mode")
    target_name = request.form.get("target")
    target = _target_or_400(target_name)

    dir_renamed = False
    original_name = None
    final_name = None

    if mode == "build":
        job_name = request.form.get("job_name", "")
        config_yaml = request.form.get("config_yaml", "")
        original_name = JS.sanitize_job_name(job_name)
        if not config_yaml.strip():
            abort(400, description="config_yaml is empty")
        files = _collect_uploaded_files()
        # For a local target, projects_path already IS a local path - stage
        # straight there, letting stage_job() collision-check against the
        # local filesystem. For an ssh target, projects_path is a *remote*
        # path (not reachable as-is from this process) - the directory that
        # actually needs collision-checking is the remote one, so that's
        # listed first and stage_job() checks against that instead; the
        # local side is just a throwaway temp staging copy, always fresh,
        # scp'd over and discarded once the copy succeeds.
        if target["kind"] == "local":
            staging_root = target["projects_path"]
            existing_names = None
        else:
            staging_root = tempfile.mkdtemp(prefix="prosculpt_dashboard_stage_")
            try:
                existing_names = RT.list_remote_dir_names(target, target["projects_path"])
            except RT.RunTargetError as e:
                shutil.rmtree(staging_root, ignore_errors=True)
                abort(502, description=str(e))

        # The config filename tracks whatever name the directory actually
        # ends up with (unless the caller explicitly asked for a different
        # one), so e.g. "my_job_2/my_job.yaml" never happens - the
        # directory and its own config always agree.
        explicit_config_filename = request.form.get("config_filename")

        try:
            staged = JS.stage_job(job_name, staging_root, files, existing_names=existing_names)
        except JS.JobStagingError as e:
            abort(400, description=str(e))
        job_dir = staged["job_dir"]
        final_name = staged["name"]
        dir_renamed = staged["renamed"]
        config_filename = explicit_config_filename or f"{final_name}.yaml"
        try:
            with open(os.path.join(job_dir, config_filename), "w") as f:
                f.write(config_yaml)
        except OSError as e:
            shutil.rmtree(job_dir, ignore_errors=True)
            abort(400, description=f"Could not write config file: {e}")

        submit_from = job_dir
        remote_job_dir = None
        reported_local_job_dir = job_dir
        if target["kind"] == "ssh":
            try:
                copy_result = RT.copy_to_remote(target, job_dir, target["projects_path"], dest_name=final_name)
            except RT.RunTargetError as e:
                shutil.rmtree(job_dir, ignore_errors=True)
                abort(502, description=f"Could not copy job directory to {target_name}: {e}")
            shutil.rmtree(job_dir, ignore_errors=True)
            if copy_result["returncode"] != 0:
                abort(502, description=f"scp failed: {copy_result['stderr'] or copy_result['stdout']}")
            remote_job_dir = f"{target['projects_path']}/{final_name}"
            submit_from = remote_job_dir
            # job_dir was just deleted (it was only ever a throwaway local
            # staging copy for the scp) - reporting it as "the local job
            # dir" would be reporting a path that no longer exists, and
            # (worse) the frontend's pending-runs poller would try to
            # watch it instead of the real remote location, which is
            # exactly what silently broke auto-tracking for ssh targets.
            reported_local_job_dir = None

    elif mode == "existing":
        existing_path = request.form.get("existing_project_path", "").strip()
        config_filename = request.form.get("config_filename", "").strip()
        if not existing_path or not config_filename:
            abort(400, description="existing_project_path and config_filename are required")
        if not os.path.isdir(existing_path):
            abort(400, description=f"Not a directory (or not reachable from this machine): {existing_path}")
        JS.ensure_logs_dir(existing_path)
        job_dir = existing_path
        reported_local_job_dir = job_dir
        if target["kind"] == "local":
            submit_from = job_dir
            remote_job_dir = None
        else:
            # A hand-assembled directory may reference files this
            # dashboard knows nothing about - it's copied as a whole,
            # verbatim, rather than reconstructed file-by-file. Still
            # collision-checked against the remote side first, same as
            # "build" mode, so a name already in use there doesn't get
            # silently merged into by scp.
            original_name = os.path.basename(os.path.normpath(existing_path))
            try:
                existing_names = RT.list_remote_dir_names(target, target["projects_path"])
            except RT.RunTargetError as e:
                abort(502, description=str(e))
            final_name, dir_renamed = JS.resolve_available_name(original_name, lambda n: n in existing_names)
            try:
                copy_result = RT.copy_to_remote(target, existing_path, target["projects_path"], dest_name=final_name)
            except RT.RunTargetError as e:
                abort(502, description=f"Could not copy project directory to {target_name}: {e}")
            if copy_result["returncode"] != 0:
                abort(502, description=f"scp failed: {copy_result['stderr'] or copy_result['stdout']}")
            remote_job_dir = f"{target['projects_path']}/{final_name}"
            submit_from = remote_job_dir
            # existing_path itself still exists (nothing deleted it, unlike
            # the "build" temp-staging case above) - but it's the source
            # directory the user pointed at, not where the job actually
            # runs/writes logs, so it shouldn't be reported as "the" local
            # job dir either, for the same reason as above.
            reported_local_job_dir = None

    else:
        abort(400, description="mode must be 'build' or 'existing'")

    dry_run = request.form.get("dry_run") in ("1", "true", "True")
    argv = [target["python_path"], f"{target['slurm_runner_path']}/slurm_runner.py", config_filename]
    if dry_run:
        argv.append("--dry-run")
    try:
        result = RT.run_remote_or_local(target, argv, cwd=submit_from)
    except RT.RunTargetError as e:
        abort(502, description=str(e))

    job_id = None
    if not dry_run:
        m = _JOB_ID_LINE_RE.search(result["stdout"]) or _FALLBACK_JOB_ID_RE.search(result["stdout"])
        if m:
            job_id = m.group(1)

    return jsonify(
        {
            "ok": result["returncode"] == 0,
            "dry_run": dry_run,
            "job_id": job_id,
            # Only set when the job's own live directory is actually
            # reachable at this path from this process - for an ssh
            # target this is None (the real directory is remote_job_dir;
            # any local path involved along the way was either a
            # throwaway staging copy already deleted, or the user's own
            # source directory, which isn't where the job's logs land).
            "local_job_dir": reported_local_job_dir,
            "remote_job_dir": remote_job_dir,
            "stdout": result["stdout"],
            "stderr": result["stderr"],
            # A directory named `original_name` already existed at the
            # destination, so `final_name` (with "_2"/"_3"/... appended)
            # was used instead - the frontend surfaces this prominently
            # rather than letting a silently-different directory name go
            # unnoticed.
            "dir_renamed": dir_renamed,
            "original_name": original_name,
            "final_name": final_name,
        }
    )


@run_api.route("/api/run/check_pending")
def api_run_check_pending():
    """Looks for a just-submitted job's log file, so the frontend can
    promote a pending submission into the tracking tab once it exists.
    `roots` (repeated query param) restricts the search to configured
    target projects_path/local_mount_path values only - never an arbitrary
    filesystem glob from the client."""
    pattern = request.args.get("glob", "")
    if not pattern:
        abort(400, description="glob is required")
    allowed_roots = []
    for t in RT.list_targets():
        try:
            full = RT.get_target(t["name"])
        except RT.RunTargetError:
            continue
        allowed_roots.append(full.get("projects_path"))
        if full.get("local_mount_path"):
            allowed_roots.append(full["local_mount_path"])
    pattern_norm = os.path.normpath(pattern)
    if not any(pattern_norm.startswith(os.path.normpath(root) + os.sep) for root in allowed_roots if root):
        abort(400, description="glob must be under a configured target's projects_path")
    matches = sorted(glob.glob(pattern))
    return jsonify({"found": bool(matches), "path": matches[0] if matches else None})
