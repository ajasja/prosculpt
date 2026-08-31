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


def _target_mounts(t: dict[str, Any]) -> list[dict[str, Any]]:
    """A target's `mounts:` list, normalized to always be a list (never
    None/missing). Each entry is `{remote: ..., local: ..., label: ...}` -
    `remote` is an absolute path as it appears on the cluster this target
    runs jobs on, `local` is wherever that same directory is actually
    reachable from this machine (a mounted/mapped drive), and `label` is
    optional (see list_browse_roots() for how it's used/defaulted).

    A target isn't limited to exactly one mount, or to mounting only
    `projects_path`/`slurm_runner_path` specifically - see
    translate_remote_path()'s docstring for why more than one remote root
    can be in play for a single target, and dashboard_config.yaml.example
    for the common two-mount shape. Nothing stops a target from mounting
    some third, unrelated remote directory too (e.g. shared scratch
    storage) - it's just another entry in this same list, with no
    dashboard code change needed to support it."""
    return t.get("mounts") or []


def _target_root_pairs(t: dict[str, Any]) -> list[tuple[Optional[str], Optional[str]]]:
    """The (remote_root, local_root) pairs worth checking a path against
    for one target, straight out of its `mounts:` list. A pair with either
    half missing is still returned (not filtered out here) - the callers
    below already skip a pair with a missing half, and returning it
    uniformly keeps this helper a pure "what are the pairs" question,
    independent of which ones happen to be usable right now."""
    return [(m.get("remote"), m.get("local")) for m in _target_mounts(t)]


def _translate_with_pairs(remote_path: Optional[str], pairs: list[tuple[Optional[str], Optional[str]]]) -> Optional[str]:
    """Core of translate_remote_path(), factored out so it can be reused
    against a `pairs` list that's already been built once (see
    make_path_translator()) instead of reloading/reparsing
    dashboard_config.yaml on every call - the difference between one
    config read per API request (translate_remote_path(), used for a
    single path like output_dir) and one per *cell* of a potentially huge
    final_output.csv (make_path_translator(), see its own docstring)."""
    if not remote_path:
        return remote_path
    remote_norm = remote_path.replace("\\", "/").rstrip("/")
    for remote_root, local_root in pairs:
        if not remote_root or not local_root:
            continue
        remote_root_norm = remote_root.replace("\\", "/").rstrip("/")
        if remote_norm == remote_root_norm or remote_norm.startswith(remote_root_norm + "/"):
            return local_root.rstrip("/\\") + remote_norm[len(remote_root_norm):]
    return remote_path


def translate_remote_path(remote_path: Optional[str]) -> Optional[str]:
    """Translates an absolute path *as it appears inside a job's own log
    file* (Hydra's own startup config dump prints PWD:/output_dir: lines
    in terms of whatever filesystem the job actually ran on - the remote
    cluster's, not necessarily wherever this dashboard process happens to
    be reading from) into whichever local path it's actually reachable at
    here, via the first configured target's `mounts:` entry whose `remote`
    is a prefix of it.

    Checking every configured mount (not just one fixed root) matters
    because a job's output_dir isn't always staged under the same remote
    directory in the first place - that's only guaranteed for a job the
    Run job tab itself submitted (always under `projects_path`; see
    job_staging.py). A job launched directly against the prosculpt
    installation instead - `python slurm_runner.py some_config.yaml
    ++output_dir=Examples/...` run by hand from `slurm_runner_path`
    itself, or via run_tests.py, both documented usages in the main
    README - gets an output_dir resolved relative to *slurm_runner_path*
    instead, since that's the job's PWD at submission time. Without a
    mount covering that root too, that path fell through to the
    "unchanged" fallback below and silently resolved to nothing on any
    machine whose mount doesn't happen to mirror the remote's entire home
    directory (see the mount-scoping note further down) - every listing
    (backbones/sequences/models, the progress numbers on the All jobs
    overview cards) came back empty despite the job's log itself still
    parsing fine, since the log's own path needs no translation at all.
    `dashboard_config.yaml.example` shows both of these as separate
    `mounts:` entries on the same target - nothing stops a target having
    more, for some other remote root a job's data can land under.

    This is the Track job side's equivalent of computeWatchDir() in
    run_job.js (which does the same projects_path -> local mount
    substitution for Run Job's own pending-submission tracking, though
    only for that one root - a job Run job itself submits is always staged
    under projects_path) - the two never shared this logic before because
    Track job's own path resolution (resolve_output_dir() in parser.py)
    predates Run Job entirely, and had until now only ever been exercised
    in setups where it worked by accident: either the dashboard runs
    directly on the same machine the job ran on (no translation needed -
    the job's own absolute paths already ARE the right local ones), or
    happened to have its mounted drive mirror the remote's entire root
    filesystem rather than just its projects directory (so a leading "/"
    coincidentally still landed in the right place, via Windows' own "no
    drive letter = root of the *current* drive" path convention) - neither
    holds once a mount is deliberately scoped to just one subtree, which
    is when this stopped working by accident and needed to actually be
    handled.

    Not just output_dir: the exact same substitution is needed for any
    other absolute path Prosculpt records at run time and the dashboard
    later reads back - e.g. final_output.csv's path_rfdiff/model_path/
    af2_pdb/af3_pdb/... columns (see load_final_csv() in parser.py, which
    uses make_path_translator() below rather than calling this function
    once per cell).

    Falls back to the path unchanged if nothing matches - covers both "no
    dashboard_config.yaml at all" and "this path isn't under any
    configured mount" (e.g. a job that predates any target/mount being
    configured), where the original direct-filesystem-access behavior is
    still exactly correct."""
    if not remote_path:
        return remote_path
    try:
        data = _load_config()
    except RunTargetError:
        return remote_path
    pairs = [pair for t in data["targets"].values() for pair in _target_root_pairs(t)]
    return _translate_with_pairs(remote_path, pairs)


def make_path_translator():
    """Loads dashboard_config.yaml *once* and returns a plain function
    `translate(path) -> path` that reuses that single load for as many
    paths as the caller needs translated - for batch use against many
    paths at once (e.g. every path-like column of a final_output.csv that
    can run into the tens of thousands of rows, see load_final_csv() in
    parser.py) where calling translate_remote_path() per-value would mean
    one disk read + YAML parse per value instead of one for the whole
    batch. Degrades to a no-op translator (returns every path unchanged)
    if there's no config at all - same fallback translate_remote_path()
    itself has, just resolved once up front instead of on every call."""
    try:
        data = _load_config()
    except RunTargetError:
        return lambda p: p
    pairs = [pair for t in data["targets"].values() for pair in _target_root_pairs(t)]
    return lambda p: _translate_with_pairs(p, pairs)


def get_default_browse_root() -> Optional[str]:
    """Best local starting point for the Track job tab's filesystem browser
    (see /api/browse in app.py): the default target's mount covering
    projects_path if it has one (most likely where actual job directories
    live), else any other usable mount it has, else projects_path itself
    directly (meaningful for a kind: local target, where it's already a
    local path and needs no mount at all - and still a reasonable starting
    guess for an ssh target with no mounts configured). None if there's no
    configured target to make any of this out of, or if the resulting path
    doesn't actually exist on this machine (a stale/wrong config entry
    shouldn't break the browser - it should just fall back to the
    caller's own default, same as before dashboard_config.yaml existed at
    all)."""
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
    t = targets[name]
    projects_path = t.get("projects_path")
    ordered = sorted(_target_mounts(t), key=lambda m: m.get("remote") != projects_path)
    for m in ordered:
        root = m.get("local")
        if root and os.path.isdir(root):
            return root
    return projects_path if projects_path and os.path.isdir(projects_path) else None


def _mount_label(target_label: str, t: dict[str, Any], m: dict[str, Any]) -> str:
    """The "Jump to..." dropdown label for one of a target's mounts -
    `target_label` plus a suffix distinguishing which remote root this
    particular mount covers, so e.g. "federico_arc" (projects_path) and
    "federico_arc (installation)" (slurm_runner_path) don't look identical
    in a dropdown that lists both. An explicit `label` on the mount entry
    always wins; otherwise this guesses from whichever of
    projects_path/slurm_runner_path the mount's `remote` matches (the two
    roots the dashboard itself cares about), falling back to the remote
    directory's own basename for anything else (e.g. a third, unrelated
    mount like shared scratch storage)."""
    explicit = m.get("label")
    if explicit:
        return f"{target_label} ({explicit})"
    remote = m.get("remote") or ""
    if remote == t.get("projects_path"):
        return target_label
    if remote == t.get("slurm_runner_path"):
        return f"{target_label} (installation)"
    basename = remote.replace("\\", "/").rstrip("/").rsplit("/", 1)[-1]
    return f"{target_label} ({basename})" if basename else target_label


def list_browse_roots() -> list[dict[str, Any]]:
    """One entry per configured target per usable mount (see
    _target_mounts()) - whichever of those actually exist on this machine
    right now. Powers the Track job file browser's "Jump to..." dropdown,
    so a multi-cluster (or multi-mount) setup doesn't require typing/
    pasting each root's path by hand to get started - a target commonly
    has more than one mount (see translate_remote_path()'s docstring for
    why) because not every job lives under the same remote root: one
    launched directly against the installation (`slurm_runner.py` run by
    hand, or via run_tests.py, rather than through the Run job tab) writes
    its output under slurm_runner_path instead of projects_path. Silently
    omits a mount whose `local` half is missing or isn't actually
    reachable right now (e.g. a network mount that's temporarily down) -
    same "don't guess, just skip it" spirit as get_default_browse_root().
    A target with no mounts configured at all still gets one entry for
    projects_path directly, same fallback get_default_browse_root() uses."""
    try:
        data = _load_config()
    except RunTargetError:
        return []
    out = []
    for name, t in data["targets"].items():
        label = t.get("label", name)
        mounts = _target_mounts(t)
        if not mounts:
            projects_path = t.get("projects_path")
            if projects_path and os.path.isdir(projects_path):
                out.append({"name": name, "label": label, "path": projects_path})
            continue
        for m in mounts:
            root = m.get("local")
            if root and os.path.isdir(root):
                out.append({"name": name, "label": _mount_label(label, t, m), "path": root})
    return out


def list_targets() -> list[dict[str, Any]]:
    """Public, UI-safe view of every configured target - never includes
    anything secret (there is nothing secret to include: no passwords/keys
    are ever stored, only hostnames and paths)."""
    data = _load_config()
    out = []
    for name, t in data["targets"].items():
        projects_path = t.get("projects_path")
        # The mount covering projects_path specifically, if any - exposed
        # under this historical field name because it's what
        # computeWatchDir() in run_job.js already reads to auto-track a
        # just-submitted job's log file without further SSH calls. A job
        # Run job itself submits is always staged under projects_path
        # (never any other mounted root), so that's the only mount this
        # particular frontend consumer ever needs to know about - it has
        # no equivalent need for a target's other mounts (e.g. one
        # covering slurm_runner_path), which only matter for tracking a
        # job that wasn't submitted through the Run job tab in the first
        # place.
        local_mount_path = next(
            (m.get("local") for m in _target_mounts(t) if m.get("remote") == projects_path),
            None,
        )
        out.append(
            {
                "name": name,
                "label": t.get("label", name),
                "kind": t.get("kind", "local"),
                "is_default": name == data.get("default_target"),
                "projects_path": projects_path,
                "local_mount_path": local_mount_path,
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
