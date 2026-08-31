"""
Assembling a job's own directory on disk - the "create a directory named
after the job, with the input pdb, an alignments/ subdirectory, and the
config file" step described for the "Build a new job" flow. Pure
filesystem logic, no subprocess/SSH calls (see run_targets.py for those) -
kept separate so it's testable/reusable independent of local-vs-remote.

Not used for the "submit an existing project directory" flow - there, the
user has already assembled their own layout (which may reference files
this module knows nothing about, e.g. a custom filter script), so nothing
here rewrites or reorganizes it.
"""

from __future__ import annotations

import os
import re

_JOB_NAME_RE = re.compile(r"[^A-Za-z0-9_.-]+")


class JobStagingError(Exception):
    """Raised for a problem assembling the job directory - callers should
    surface `str(exc)` to the user rather than a stack trace."""


def sanitize_job_name(job_name: str) -> str:
    name = job_name.strip()
    if not name:
        raise JobStagingError("Job name cannot be empty")
    name = _JOB_NAME_RE.sub("_", name)
    if name in (".", ".."):
        raise JobStagingError("Invalid job name")
    return name


def resolve_available_name(name: str, is_taken) -> tuple[str, bool]:
    """Finds the first of `name`, `name_2`, `name_3`, ... for which
    `is_taken(candidate)` is False. Returns (final_name, renamed)."""
    if not is_taken(name):
        return name, False
    counter = 2
    while True:
        candidate = f"{name}_{counter}"
        if not is_taken(candidate):
            return candidate, True
        counter += 1


def stage_job(
    job_name: str,
    projects_path: str,
    files: dict[str, bytes | str],
    existing_names: set[str] | None = None,
) -> dict:
    """Creates `<projects_path>/<name>/` and writes every entry of `files`
    into it (relative path -> content; str as text, bytes as binary; parent
    directories, e.g. "alignments/x.a3m", are created as needed), then
    ensures logs/ exists.

    If a directory named after the (sanitized) job name already exists,
    this does NOT overwrite it or error out - it appends "_2", "_3", ... to
    the name until it finds one that's free, exactly like a browser
    downloading "file.pdf" a second time gets "file (1).pdf" instead of
    clobbering the first. The caller is expected to surface `renamed` to
    the user so a silently-different directory name is never a surprise.

    `existing_names`, if given, is checked instead of the local filesystem -
    needed for the ssh-target case, where the directory actually being
    collision-checked is on the *remote* side (see run_api.py, which lists
    the remote projects_path first) and a local os.path.exists() check
    against a throwaway local staging temp dir would never find anything.
    """
    name = sanitize_job_name(job_name)
    if existing_names is None:
        if not os.path.isdir(projects_path):
            raise JobStagingError(f"Configured projects_path does not exist or isn't reachable: {projects_path}")
        is_taken = lambda n: os.path.exists(os.path.join(projects_path, n))  # noqa: E731
    else:
        is_taken = lambda n: n in existing_names  # noqa: E731

    final_name, renamed = resolve_available_name(name, is_taken)
    job_dir = os.path.join(projects_path, final_name)
    os.makedirs(job_dir)
    try:
        for rel_path, content in files.items():
            full_path = os.path.join(job_dir, rel_path)
            full_norm = os.path.normpath(full_path)
            job_dir_norm = os.path.normpath(job_dir)
            if not (full_norm == job_dir_norm or full_norm.startswith(job_dir_norm + os.sep)):
                raise JobStagingError(f"Refusing to write outside the job directory: {rel_path!r}")
            os.makedirs(os.path.dirname(full_path), exist_ok=True)
            mode = "w" if isinstance(content, str) else "wb"
            with open(full_path, mode) as f:
                f.write(content)
        # No ensure_logs_dir() here anymore - slurm_runner.py creates
        # <output_dir>/logs/ itself, right before it actually needs it (see
        # its own comment), which covers both submission modes uniformly
        # without this module needing to guess where logs will end up.
    except Exception:
        # Don't leave a half-written job directory behind - either the
        # whole thing is there correctly or none of it is.
        import shutil

        shutil.rmtree(job_dir, ignore_errors=True)
        raise
    return {"job_dir": job_dir, "name": final_name, "renamed": renamed}
