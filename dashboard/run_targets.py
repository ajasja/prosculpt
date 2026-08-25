"""
Run targets - where a job actually gets submitted from.

A "target" is one destination the dashboard can submit a job to: either
"local" (the dashboard's own process runs `slurm_runner.py` directly, e.g.
because the dashboard itself is running on the cluster's login node) or
"ssh" (the dashboard shells out to the system `ssh`/`scp` binaries to stage
files and submit on a remote cluster it's on the same network as, using
whatever key/agent trust already exists for that host - no credentials of
any kind are read, stored, or passed by this module).

Targets are configured in `dashboard_config.yaml` next to this file
(git-ignored, host-specific - see dashboard_config.yaml.example for the
shape). This mirrors how `config/installation.yaml` already works for the
pipeline itself: a small, hand-edited YAML file rather than a database or
a settings UI. The same file also holds dashboard-wide `defaults:` (e.g.
which prediction model to preselect) unrelated to any one target - see
get_defaults() - so there's one config file for "how this dashboard talks
to a cluster and what it should default to", not several.
"""

from __future__ import annotations

import os
import shlex
import subprocess
from typing import Any, Optional

try:
    import yaml
except ImportError:  # pragma: no cover - matches parser.py's own guard
    yaml = None

_CONFIG_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "dashboard_config.yaml")

# Commands are given a generous timeout since `sbatch`/`ssh` to a busy login
# node can be slow, but this still bounds a hung connection from blocking a
# request forever.
_TIMEOUT_SECONDS = 60


class RunTargetError(Exception):
    """Raised for a config or usage problem - callers should surface
    `str(exc)` to the user rather than a stack trace."""


def _load_raw() -> dict[str, Any]:
    """The config file's contents as-is, or {} if it's missing/unreadable -
    used both by _load_config() (which additionally requires `targets` to
    be present) and by get_defaults() (which doesn't - dashboard-wide
    defaults should still work even before any target has been set up)."""
    if yaml is None or not os.path.isfile(_CONFIG_PATH):
        return {}
    with open(_CONFIG_PATH, "r") as f:
        return yaml.safe_load(f) or {}


def _load_config() -> dict[str, Any]:
    if yaml is None:
        raise RunTargetError("PyYAML is not installed - cannot read dashboard_config.yaml")
    if not os.path.isfile(_CONFIG_PATH):
        raise RunTargetError(
            f"No dashboard_config.yaml found at {_CONFIG_PATH} - copy dashboard_config.yaml.example "
            "and fill in at least one target before submitting a job."
        )
    data = _load_raw()
    if "targets" not in data or not data["targets"]:
        raise RunTargetError("dashboard_config.yaml has no targets defined")
    return data


def get_defaults() -> dict[str, Any]:
    """Dashboard-wide defaults (e.g. which prediction model the Run Job GUI
    should preselect) - lives under a `defaults:` key in the same config
    file as targets, so this file can grow to hold other dashboard-level
    settings later without a second config file. Always {} rather than
    raising - dashboard defaults are a nice-to-have, never a hard
    requirement the way a target's own settings are."""
    return _load_raw().get("defaults") or {}


def get_default_browse_root() -> Optional[str]:
    """Best local starting point for the Track job tab's filesystem browser
    (see /api/browse in app.py): the default target's local_mount_path (an
    ssh target's remote output, reachable locally) if it has one, else its
    projects_path (directly meaningful for a local target, and still a
    reasonable starting guess for an ssh target with no mount configured).
    None if there's no configured target to make any of this out of, or if
    the resulting path doesn't actually exist on this machine (a stale/
    wrong config entry shouldn't break the browser - it should just fall
    back to the caller's own default, same as before dashboard_config.yaml
    existed at all)."""
    try:
        data = _load_config()
    except RunTargetError:
        return None
    targets = data["targets"]
    name = data.get("default_target")
    if not name or name not in targets:
        # No explicit default - same "good enough, not going to guess
        # harder" spirit as get_target()'s own fallback: just take
        # whichever target happens to be first.
        name = next(iter(targets))
    root = targets[name].get("local_mount_path") or targets[name].get("projects_path")
    return root if root and os.path.isdir(root) else None


def list_targets() -> list[dict[str, Any]]:
    """Public, UI-safe view of every configured target - never includes
    anything secret (there is nothing secret to include: no passwords/keys
    are ever stored, only hostnames and paths)."""
    data = _load_config()
    out = []
    for name, t in data["targets"].items():
        out.append(
            {
                "name": name,
                "label": t.get("label", name),
                "kind": t.get("kind", "local"),
                "is_default": name == data.get("default_target"),
                "projects_path": t.get("projects_path"),
                # Only meaningful for an ssh target whose output is *also*
                # reachable as a normal mounted path from this machine - lets
                # the frontend translate a remote job directory back into a
                # local path it can poll for a just-appeared log file
                # without any further SSH calls. None if not configured.
                "local_mount_path": t.get("local_mount_path"),
            }
        )
    return out


def get_target(name: Optional[str] = None) -> dict[str, Any]:
    data = _load_config()
    targets = data["targets"]
    if name is None:
        name = data.get("default_target")
    if not name or name not in targets:
        raise RunTargetError(f"Unknown run target: {name!r}")
    t = dict(targets[name])
    t["name"] = name
    t.setdefault("kind", "local")
    # The interpreter to invoke slurm_runner.py with. Defaults to a bare
    # "python" (whatever's first on PATH) if not set, but that's exactly
    # the thing that tends to silently break: the dashboard's own process
    # (local case) or a non-interactive SSH session (remote case, where
    # `conda activate` commonly isn't even available - see the README)
    # usually isn't running inside prosculpt's own conda env. Setting this
    # to that env's absolute interpreter path (the same value as
    # installation.yaml's own prosculpt_python_path) sidesteps shell
    # activation entirely - invoking a conda env's bin/python directly
    # gets you that env's packages regardless of shell/session type.
    t.setdefault("python_path", "python")
    if t["kind"] not in ("local", "ssh"):
        raise RunTargetError(f"Target {name!r} has unknown kind {t['kind']!r} (expected 'local' or 'ssh')")
    if t["kind"] == "local":
        for key in ("slurm_runner_path", "projects_path"):
            if not t.get(key):
                raise RunTargetError(f"Target {name!r} is missing required key {key!r}")
    else:
        for key in ("ssh_host", "slurm_runner_path", "projects_path"):
            if not t.get(key):
                raise RunTargetError(f"Target {name!r} is missing required key {key!r}")
    return t


def run_remote_or_local(target: dict[str, Any], argv: list[str], cwd: Optional[str] = None) -> dict[str, Any]:
    """Runs `argv` either as a local subprocess (cwd=`cwd`) or, for an ssh
    target, over `ssh <host> 'cd <cwd> && <argv...>'` - one shared code path
    so job-submission and squeue never diverge in how they decide
    local-vs-remote. Returns {"returncode", "stdout", "stderr"}; never
    raises for a nonzero exit (that's a normal "the command failed" outcome
    the caller should inspect), only for things that mean the command
    couldn't even be attempted (e.g. `ssh`/`scp` not installed).
    """
    try:
        if target["kind"] == "local":
            proc = subprocess.run(
                argv, cwd=cwd, capture_output=True, text=True, timeout=_TIMEOUT_SECONDS
            )
        else:
            remote_cmd = shlex.join(argv)
            if cwd:
                remote_cmd = f"cd {shlex.quote(cwd)} && {remote_cmd}"
            proc = subprocess.run(
                ["ssh", target["ssh_host"], remote_cmd],
                capture_output=True,
                text=True,
                timeout=_TIMEOUT_SECONDS,
            )
    except FileNotFoundError as e:
        # For a local target this is virtually always argv[0] itself
        # missing (e.g. a wrong/unreachable python_path) - a *remote*
        # command not being found comes back as a nonzero exit + stderr
        # over the ssh connection instead, not a local Python exception,
        # so this branch is only ever the `ssh` binary itself missing.
        if target["kind"] == "local":
            raise RunTargetError(f"Could not run command ({e}) - check this target's python_path/slurm_runner_path.") from e
        raise RunTargetError(f"Could not run command ({e}) - is ssh installed and on PATH?") from e
    except subprocess.TimeoutExpired as e:
        raise RunTargetError(f"Command timed out after {_TIMEOUT_SECONDS}s: {e}") from e
    return {"returncode": proc.returncode, "stdout": proc.stdout, "stderr": proc.stderr}


def list_remote_dir_names(target: dict[str, Any], remote_dir: str) -> set[str]:
    """Names of the entries directly inside `remote_dir` on an ssh target -
    used to pick a collision-free job directory name before copying
    anything there (see job_staging.resolve_available_name()). Raises
    RunTargetError if the listing itself couldn't be done (e.g. the
    directory doesn't exist or isn't reachable) - deliberately NOT treated
    as "empty directory", since that could let a real problem (a typo'd
    projects_path) through silently as if there were no collision."""
    result = run_remote_or_local(target, ["ls", "-1", remote_dir])
    if result["returncode"] != 0:
        raise RunTargetError(
            f"Could not list {remote_dir!r} on {target['ssh_host']} to check for a name collision: "
            f"{result['stderr'].strip() or result['stdout'].strip() or 'ls failed'}"
        )
    return {line.strip() for line in result["stdout"].splitlines() if line.strip()}


def copy_to_remote(
    target: dict[str, Any], local_dir: str, remote_parent_dir: str, dest_name: Optional[str] = None
) -> dict[str, Any]:
    """`scp -r` a whole local directory to `<remote_parent_dir>/<dest_name>`
    (defaulting to local_dir's own basename) on an ssh target. Only
    meaningful for kind == "ssh". Callers that care about not silently
    merging into an already-existing remote directory (scp's own behavior
    when the destination already exists) should resolve a free `dest_name`
    first via list_remote_dir_names() + job_staging.resolve_available_name()."""
    if target["kind"] != "ssh":
        raise RunTargetError("copy_to_remote() only applies to ssh targets")
    name = dest_name or os.path.basename(os.path.normpath(local_dir))
    try:
        proc = subprocess.run(
            # -s forces the SFTP-based transfer protocol instead of scp's
            # legacy one. The legacy protocol reads the very first bytes
            # of the remote shell's own stdout as its binary handshake -
            # any shell startup script that prints so much as one stray
            # line before that (a login banner, a conda-activation
            # message, a module-load echo, ... all common on HPC login
            # shells) corrupts it with a cryptic "Received message too
            # long" error that has nothing to do with the actual transfer.
            # SFTP mode runs over a proper subsystem channel instead, so
            # it isn't vulnerable to that whole class of failure - safe to
            # force unconditionally, not just as a fix for one target's
            # particular shell config.
            ["scp", "-s", "-r", local_dir, f"{target['ssh_host']}:{remote_parent_dir}/{name}"],
            capture_output=True,
            text=True,
            timeout=_TIMEOUT_SECONDS * 5,  # a whole directory (pdb/alignments) can be slower than a bare command
        )
    except FileNotFoundError as e:
        raise RunTargetError(f"Could not run scp ({e}) - is scp installed and on PATH?") from e
    except subprocess.TimeoutExpired as e:
        raise RunTargetError(f"Copy timed out: {e}") from e
    return {"returncode": proc.returncode, "stdout": proc.stdout, "stderr": proc.stderr}
