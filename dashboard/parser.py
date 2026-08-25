"""
parser.py

All the "understand what a Prosculpt run is doing" logic lives here, kept
separate from the Flask routes in app.py so it can be tested / reused on its
own.

Two data sources are combined, deliberately:

  1. The log file - tells us *when* things happened, what the current stage
     is, and lets us compute rates (time per backbone, time per model) that
     the filesystem alone can't give us.

  2. The output directory on disk - tells us *what has actually been
     produced so far* (ground truth counts / listings). This is more
     reliable than trying to count log lines, because filters, retries,
     crashes, etc. can all make the log's bookkeeping diverge from what's
     really on disk.

Nothing here hard-codes the exact schema of AF3's summary_confidences.json
or Boltz's confidence_*.json - those are parsed generically (every scalar
key/value pair is surfaced) so the dashboard keeps working even if a field
gets renamed/added/removed in a future tool version. The same applies to the
run configuration (see extract_config) - every key Prosculpt was actually
run with is surfaced, not just a hard-coded subset.
"""

from __future__ import annotations

import glob
import json
import os
import pickle
import re
from datetime import datetime
from typing import Any, Optional

# No Flask dependency here either (see run_targets.py's own docstring) -
# only used to translate a job's own remote-filesystem paths (see
# resolve_output_dir() below) into wherever they're actually reachable
# from this machine, when dashboard_config.yaml has a matching target
# configured. Degrades to a no-op with no config file at all, so this
# module's "test/reuse standalone" promise still holds.
import run_targets as RT

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


TIMESTAMP_RE = re.compile(r"\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}),(\d{3})\]")


def _parse_ts(line: str) -> Optional[datetime]:
    m = TIMESTAMP_RE.search(line)
    if not m:
        return None
    dt = datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S")
    return dt.replace(microsecond=int(m.group(2)) * 1000)


# ---------------------------------------------------------------------------
# Locating the job
# ---------------------------------------------------------------------------


def resolve_output_dir(log_path: str) -> dict:
    """Reconstruct the absolute output directory from the log's PWD: and
    output_dir: lines, per the convention described by the pipeline author.
    """
    with open(log_path, "r", errors="replace") as f:
        lines = f.readlines()

    pwd = None
    for i, line in enumerate(lines):
        if line.strip() == "PWD:":
            if i + 1 < len(lines):
                pwd = lines[i + 1].strip()
            break

    output_dir_rel = None
    for line in lines:
        m = re.match(r"^output_dir:\s*(\S.*)$", line.strip())
        if m:
            output_dir_rel = m.group(1).strip().strip("'\"")
            break

    if not pwd or not output_dir_rel:
        return {
            "pwd": pwd,
            "output_dir_raw": output_dir_rel,
            "output_dir": None,
            "error": "Could not find PWD: and/or output_dir: in the log file.",
        }

    if os.path.isabs(output_dir_rel):
        full = output_dir_rel
    else:
        full = os.path.normpath(os.path.join(pwd, output_dir_rel))

    # `full` at this point is expressed in terms of wherever the job
    # actually ran (PWD:/output_dir: come straight from the log's own
    # Hydra config dump) - not necessarily a path this dashboard process
    # can read directly. See translate_remote_path()'s own docstring for
    # why this was silently wrong (empty backbones/sequences/models tabs,
    # no progress on the All jobs overview cards, despite the files
    # genuinely being on disk) for any setup whose mounted drive doesn't
    # happen to mirror the remote's entire root filesystem.
    full = RT.translate_remote_path(full)

    return {"pwd": pwd, "output_dir_raw": output_dir_rel, "output_dir": full, "error": None}


# ---------------------------------------------------------------------------
# Config extraction
#
# Prosculpt logs the *fully resolved* Hydra config once, near the top of
# every run, as a plain OmegaConf YAML dump:
#
#   [2026-07-23 17:03:52,833][__main__][INFO] - The following configuration was passed:
#   task_name: bltz_binder
#   slurm:
#     output: logs/slurm-%A_%a_%x.out
#     ...
#
# That block is a strict superset of whatever the user's input yaml
# contained (defaults get merged in, output_dir gets resolved to an
# absolute path, etc.), so it's the primary source here. input.yaml (a raw
# copy of the file the user passed to Hydra, written to the output dir
# early in the run) is only used as a fallback for keys the log block
# didn't have - e.g. if the log got truncated before that point, or is in
# an older format.
#
# Rather than hard-coding which keys to look for, the whole block is
# yaml-parsed and *every* top-level key is surfaced, so task-specific
# options (e.g. boltz2_templates, symmetry, ppi) show up automatically.
# ---------------------------------------------------------------------------

_CONFIG_DUMP_MARKER_RE = re.compile(r"The following configuration was passed:\s*$")


def _extract_yaml_config_block(lines: list[str]) -> Optional[str]:
    start = None
    for i, line in enumerate(lines):
        if _CONFIG_DUMP_MARKER_RE.search(line):
            start = i + 1
            break
    if start is None:
        return None

    block_lines = []
    for line in lines[start:]:
        # The block ends where the next properly-prefixed log line begins.
        if TIMESTAMP_RE.match(line):
            break
        block_lines.append(line)
    return "".join(block_lines) if block_lines else None


# Keys the rest of the dashboard relies on existing (with sane fallbacks),
# on top of whatever else gets surfaced generically.
_DEFAULTS: dict[str, Any] = {
    "num_designs_rfdiff": None,
    "num_seq_per_target_mpnn": None,
    "num_models": None,
    "af2_mpnn_cycles": 1,
    "model_monomer": False,
    "prediction_model": None,
    "task_name": None,
    "contig": None,
}


def extract_config(log_path: str, output_dir: Optional[str]) -> dict:
    cfg: dict[str, Any] = {}

    with open(log_path, "r", errors="replace") as f:
        lines = f.readlines()

    block = _extract_yaml_config_block(lines)
    if block and yaml is not None:
        try:
            parsed = yaml.safe_load(block)
            if isinstance(parsed, dict):
                cfg.update(parsed)
        except Exception:
            pass

    # Fallback / fill-in: read input.yaml directly if the log block was
    # missing or incomplete (e.g. truncated log, or a key it doesn't have).
    if output_dir and yaml is not None:
        yaml_path = os.path.join(output_dir, "input.yaml")
        if os.path.isfile(yaml_path):
            try:
                with open(yaml_path) as f:
                    doc = yaml.safe_load(f) or {}
                if isinstance(doc, dict):
                    for key, val in doc.items():
                        cfg.setdefault(key, val)
            except Exception:
                pass

    for key, default in _DEFAULTS.items():
        cfg.setdefault(key, default)

    # model_monomer sometimes arrives as a string ("true"/"True") depending
    # on how it was set on the command line - normalize it.
    if isinstance(cfg.get("model_monomer"), str):
        cfg["model_monomer"] = cfg["model_monomer"].strip().lower() == "true"

    if isinstance(cfg.get("af2_mpnn_cycles"), str):
        try:
            cfg["af2_mpnn_cycles"] = int(cfg["af2_mpnn_cycles"])
        except ValueError:
            pass

    return cfg


# ---------------------------------------------------------------------------
# Job / timing info (job id, node, start time, per-step durations)
# ---------------------------------------------------------------------------

JOB_HELLO_RE = re.compile(r"Hello from job (\d+) on (\S+) at (.+)$")
STEP_DURATION_RE = re.compile(r"\* \* \* (.+?) lasted ([\d.]+) s\. Running (.+?) \* \* \*")
TOTAL_DURATION_RE = re.compile(r"Total job duration:\s*(\d+)d (\d+)h (\d+)m (\d+)s")

STEP_LABELS = {
    "Before app start": "Startup",
    "general_config_prep": "Config preparation",
    "pass_config_to_rfdiff": "Preparing RFdiffusion config",
    "run_rfdiff": "RFdiffusion",
    "plugin_filters": "Backbone filtering",
    "rechain_rfdiff_pdbs": "Rechaining backbones",
    "skipping RfDiff": "Skipped RFdiffusion (input PDB provided)",
    "do_cycling": "ProteinMPNN + structure modeling",
    "final_operations": "Final scoring",
    "only final_operations": "Final scoring (analysis-only run)",
}


def _friendly_step_label(name: str) -> str:
    if name in STEP_LABELS:
        return STEP_LABELS[name]
    return name.replace("_", " ").strip().capitalize()


def parse_job_info(lines: list[str]) -> dict:
    for line in lines:
        m = JOB_HELLO_RE.search(line)
        if m:
            return {
                "job_id": m.group(1),
                "node": m.group(2),
                "started_at": m.group(3).strip(),
            }
    return {"job_id": None, "node": None, "started_at": None}


def parse_step_durations(lines: list[str]) -> dict:
    order: list[str] = []
    seconds_by_name: dict[str, float] = {}
    total_seconds = None

    for line in lines:
        m = STEP_DURATION_RE.search(line)
        if m:
            name, secs = m.group(1), float(m.group(2))
            if name not in seconds_by_name:
                order.append(name)
            seconds_by_name[name] = secs
            continue
        m = TOTAL_DURATION_RE.search(line)
        if m:
            d, h, mnt, s = (int(x) for x in m.groups())
            total_seconds = d * 86400 + h * 3600 + mnt * 60 + s

    steps = [
        {"name": name, "label": _friendly_step_label(name), "seconds": seconds_by_name[name]}
        for name in order
    ]
    return {"steps": steps, "total_seconds": total_seconds}


# ---------------------------------------------------------------------------
# Cycle info (af2_mpnn_cycles > 1)
# ---------------------------------------------------------------------------

# Not a proper log.info() call (no timestamp prefix) - just a stray print()
# in do_cycling(), but it's the only place the current cycle index is
# recorded, so we use it.
CYCLE_RE = re.compile(r"^cycleeeeee\s+(\d+)\s*$")


def current_cycle_info(lines: list[str], cfg: dict) -> dict:
    total_cycles = cfg.get("af2_mpnn_cycles") or 1
    current = 0
    seen_any = False
    for line in lines:
        m = CYCLE_RE.match(line.strip())
        if m:
            current = int(m.group(1))
            seen_any = True
    return {
        "current_cycle": current,
        "current_cycle_display": current + 1,
        "total_cycles": total_cycles,
        "is_multi_cycle": total_cycles > 1,
        "seen": seen_any,
    }


# ---------------------------------------------------------------------------
# Crash / cancellation / error-file detection
#
# The slurm .err file (same directory as the log, same basename, .err
# extension) is looked up unconditionally, from the very first status
# fetch - not only once a crash is already suspected - so the dashboard can
# say plainly whether one exists at all, and so it can be tailed on every
# poll to catch a Slurm-level cancellation (which never touches the .out
# log the rest of this module reads).
# ---------------------------------------------------------------------------

CRASH_RE = re.compile(r"There was an error running the command\.")


def detect_crash(lines: list[str]) -> dict:
    return {"crashed": any(CRASH_RE.search(line) for line in lines)}


def locate_err_file(log_path: str) -> dict:
    base, _ext = os.path.splitext(log_path)
    err_path = base + ".err"
    return {"err_path": err_path, "err_exists": os.path.isfile(err_path)}


TEXT_FILE_MAX_BYTES = 5_000_000  # safety cap so a huge log/traceback dump can't hang the browser
_CANCEL_TAIL_BYTES = 200_000  # Slurm appends its cancellation notice at the very end
CANCELLED_RE = re.compile(r"\bCANCELLED\b", re.IGNORECASE)


def read_text_file(path: str) -> dict:
    """Generic capped text read, used for both the .err file (Error log
    tab) and the .out log itself (Output log tab)."""
    size = os.path.getsize(path)
    truncated = size > TEXT_FILE_MAX_BYTES
    with open(path, "r", errors="replace") as f:
        content = f.read(TEXT_FILE_MAX_BYTES) if truncated else f.read()
    return {"content": content, "truncated": truncated, "size": size}


def scan_err_file_for_cancellation(err_path: str) -> dict:
    """Slurm writes a line like '*** JOB 123 ON node CANCELLED AT ... ***'
    to stderr when a job is cancelled (by the user or by the scheduler,
    e.g. hitting a time limit) - this never shows up in the .out log the
    rest of the module works from, so the .err file has to be checked
    directly. Only the tail is read since the notice is always the last
    thing Slurm appends."""
    if not os.path.isfile(err_path):
        return {"cancelled": False, "cancelled_line": None}
    try:
        size = os.path.getsize(err_path)
        with open(err_path, "rb") as f:
            if size > _CANCEL_TAIL_BYTES:
                f.seek(-_CANCEL_TAIL_BYTES, os.SEEK_END)
            tail = f.read().decode("utf-8", errors="replace")
    except OSError:
        return {"cancelled": False, "cancelled_line": None}
    for line in tail.splitlines():
        if CANCELLED_RE.search(line):
            return {"cancelled": True, "cancelled_line": line.strip()}
    return {"cancelled": False, "cancelled_line": None}


# ---------------------------------------------------------------------------
# Stage detection - a simple linear state machine over the log lines
# ---------------------------------------------------------------------------

STAGES = [
    "setup",
    "rfdiffusion",
    "filtering",
    "mpnn",
    "modeling",
    "scoring",
    "finished",
]

_STAGE_TRIGGERS: list[tuple[re.Pattern, str]] = [
    (re.compile(r"Running run_rfdiff"), "rfdiffusion"),
    (re.compile(r"Running plugin_filters"), "filtering"),
    (re.compile(r"Running rechain_rfdiff_pdbs"), "filtering"),
    (re.compile(r"Running do_cycling"), "mpnn"),
    (re.compile(r"protein_mpnn_run\.py"), "mpnn"),
    (re.compile(r"Running AlphaFold 3"), "modeling"),
    (re.compile(r"run_alphafold\.py"), "modeling"),
    (re.compile(r"boltz predict "), "modeling"),
    (re.compile(r"Running final_operations"), "scoring"),
    (re.compile(r"Final operations"), "scoring"),
    (re.compile(r"Running Finished"), "finished"),
    (re.compile(r"Total job duration"), "finished"),
]

_ERROR_RE = re.compile(r"(Traceback \(most recent call last\)|Error:|CUDA out of memory|srun: error|slurmstepd: error)")


def detect_stage(lines: list[str]) -> dict:
    stage = "setup"
    last_line_no = -1
    error_lines: list[str] = []
    for i, line in enumerate(lines):
        for pattern, target in _STAGE_TRIGGERS:
            if pattern.search(line):
                stage = target
                last_line_no = i
                break
        if _ERROR_RE.search(line):
            error_lines.append(line.strip())

    return {
        "stage": stage,
        "stage_trigger_line": last_line_no,
        "possible_errors": error_lines[-5:],  # last few, so we don't flood the UI
    }


# ---------------------------------------------------------------------------
# RFdiffusion stage details
# ---------------------------------------------------------------------------

MAKING_DESIGN_RE = re.compile(r"Making design .*?_(\d+)\s*$")


def rfdiffusion_progress(lines: list[str], output_dir: str, total: Optional[int]) -> dict:
    events: list[tuple[int, datetime]] = []
    for line in lines:
        m = MAKING_DESIGN_RE.search(line)
        if m:
            ts = _parse_ts(line)
            if ts:
                events.append((int(m.group(1)), ts))
    events_by_index = dict(events)

    # Ground truth: pdb files actually written to disk so far.
    completed_files = glob.glob(os.path.join(output_dir, "1_rfdiff", "_*.pdb"))
    completed = len(completed_files)

    # Per-design duration = time between successive "Making design" log
    # lines (independent of file-write timing, which lags slightly).
    durations = []
    for (idx_a, ts_a), (idx_b, ts_b) in zip(events, events[1:]):
        durations.append((ts_b - ts_a).total_seconds())

    avg_duration = sum(durations) / len(durations) if durations else None

    now = datetime.now()
    current_index = None
    current_elapsed = None
    if total is not None and completed < total:
        current_index = completed  # 0-indexed design currently in flight
        start_ts = events_by_index.get(current_index)
        if start_ts is None and events:
            start_ts = events[-1][1]  # fallback: most recent known start
        if start_ts is not None:
            current_elapsed = (now - start_ts).total_seconds()

    eta_seconds = None
    if avg_duration is not None and total is not None:
        remaining = max(total - completed, 0)
        eta_seconds = avg_duration * remaining

    return {
        "completed": completed,
        "total": total,
        "current_index": current_index,
        "current_elapsed_seconds": current_elapsed,
        "avg_seconds_per_backbone": avg_duration,
        "eta_seconds": eta_seconds,
    }


def list_backbones(output_dir: str, filtering_done: bool) -> list[dict]:
    """`filtering_done` reflects whether plugin_filters has actually
    finished (the caller derives this from detect_stage() - stage has
    reached "mpnn" or later). Until then, a backbone still sitting in
    1_rfdiff/ hasn't necessarily passed anything - it just hasn't been
    evaluated yet - so it's tagged "pending" rather than "passed", to
    avoid implying a filter verdict that hasn't actually happened."""
    rfdiff_dir = os.path.join(output_dir, "1_rfdiff")
    out = []
    if not os.path.isdir(rfdiff_dir):
        return out
    not_yet_filtered_status = "passed" if filtering_done else "pending"

    def _trb_path(pdb_path: str) -> Optional[str]:
        candidate = os.path.splitext(pdb_path)[0] + ".trb"
        return candidate if os.path.isfile(candidate) else None

    for pdb in sorted(glob.glob(os.path.join(rfdiff_dir, "_*.pdb"))):
        out.append({
            "name": os.path.splitext(os.path.basename(pdb))[0],
            "status": not_yet_filtered_status,
            "path": pdb,
            "trb_path": _trb_path(pdb),
        })
    failed_dir = os.path.join(rfdiff_dir, "failed_filters")
    if os.path.isdir(failed_dir):
        for pdb in sorted(glob.glob(os.path.join(failed_dir, "_*.pdb"))):
            out.append({
                "name": os.path.splitext(os.path.basename(pdb))[0],
                "status": "failed_filter",
                "path": pdb,
                "trb_path": _trb_path(pdb),
            })
    out.sort(key=lambda r: _natural_key(r["name"]))
    return out


def _natural_key(s: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


# ---------------------------------------------------------------------------
# RFdiffusion .trb sidecar (residue provenance: motif / fixed chain / de
# novo) - used for the "color by RFdiffusion provenance" viewer option.
# ---------------------------------------------------------------------------


def load_trb_provenance(trb_path: str) -> dict:
    """A .trb file is a pickled dict RFdiffusion writes alongside each
    backbone .pdb. Two fields drive this, both lists of (chain, resnum)
    pairs in the *generated* structure's own numbering:

      con_hal_pdb_idx         - residues taken from the reference structure
                                 that ended up in a redesigned chain
                                 ("Motif").
      complex_con_hal_pdb_idx - every non-sculpted residue (motif AND fixed
                                 chains together), correctly chained and
                                 numbered. Only present at all when the run
                                 actually had fixed chains.

    "Fixed chains" is derived as complex_con_hal_pdb_idx minus
    con_hal_pdb_idx (set difference on (chain, resnum) pairs) rather than
    read from receptor_con_hal_pdb_idx directly - that field turns out to
    get both the resnums AND the chain letters wrong (it doesn't reset
    resnums per chain, and mislabels which chain a residue belongs to).
    complex_con_hal_pdb_idx, by contrast, already uses correct chain
    letters and correct (chain-local) residue numbers, so no renumbering
    is needed on the result - unlike the old receptor_con_hal_pdb_idx-based
    approach. If complex_con_hal_pdb_idx is absent from the pickle at all,
    the run had no fixed chains - every non-sculpted residue is already
    covered by con_hal_pdb_idx, so fixed_chain is simply empty.

    Any residue in neither list was generated de novo ("Sculpted") - that
    set isn't enumerated here since the caller already knows the full
    residue list from the structure itself.
    """
    try:
        with open(trb_path, "rb") as f:
            data = pickle.load(f)
    except Exception as e:
        return {
            "error": (
                f"Could not read this .trb file ({e}). RFdiffusion .trb files are "
                "Python pickles that can reference numpy/torch objects - make sure "
                "those packages are installed in the dashboard's own Python "
                "environment (see requirements.txt)."
            )
        }

    def normalize(entries):
        out = []
        for entry in entries or []:
            try:
                chain, resnum = entry[0], entry[1]
                out.append([str(chain), int(resnum)])
            except Exception:
                continue  # skip anything not shaped like (chain, resnum)
        return out

    motif = normalize(data.get("con_hal_pdb_idx"))

    fixed_chain = []
    if "complex_con_hal_pdb_idx" in data:
        motif_set = {tuple(entry) for entry in motif}
        fixed_chain = [
            entry for entry in normalize(data.get("complex_con_hal_pdb_idx"))
            if tuple(entry) not in motif_set
        ]

    return {
        "motif": motif,
        "fixed_chain": fixed_chain,
    }


# ---------------------------------------------------------------------------
# MPNN stage details
# ---------------------------------------------------------------------------

FASTA_HEADER_RE = re.compile(
    r"score=([\d.\-]+).*?global_score=([\d.\-]+).*?seq_recovery=([\d.\-]+)"
)


def _parse_fasta(path: str) -> list[dict]:
    """A ProteinMPNN fasta file: while the pipeline is running, the first
    record is the input/reference sequence (no 'sample=' in its header);
    it gets removed right before the MPNN stage ends and we don't care
    about it either way, so we filter it out. Chains within a sequence are
    colon-separated (confirmed from real pipeline output)."""
    records = []
    header = None
    seq_lines: list[str] = []

    def flush():
        if header is None:
            return
        seq = "".join(seq_lines)
        chains = seq.split(":")
        meta = {}
        m = FASTA_HEADER_RE.search(header)
        if m:
            meta = {
                "score": float(m.group(1)),
                "global_score": float(m.group(2)),
                "seq_recovery": float(m.group(3)),
            }
        is_input = "sample=" not in header
        records.append(
            {
                "header": header,
                "chains": chains,
                "sequence": seq,
                "is_input": is_input,
                **meta,
            }
        )

    if not os.path.isfile(path):
        return records
    with open(path, errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
    flush()
    return records


def list_sequences(output_dir: str) -> dict:
    seqs_dir = os.path.join(output_dir, "2_mpnn", "seqs")
    backbones = []
    if os.path.isdir(seqs_dir):
        # NOTE: filenames are '_<n>.fa' only in the *first* MPNN/AF cycle.
        # From cycle 1 onwards (af2_mpnn_cycles > 1) the pipeline re-feeds
        # the previous cycle's models back into MPNN under names like
        # 'rf_0__model_1__cycle_1__itr_0__.fa' - so we can't anchor on a
        # leading underscore, just take every .fa file in the directory.
        for fa in sorted(glob.glob(os.path.join(seqs_dir, "*.fa")), key=lambda p: _natural_key(os.path.basename(p))):
            name = os.path.splitext(os.path.basename(fa))[0]
            records = _parse_fasta(fa)
            samples = [r for r in records if not r["is_input"]]
            backbones.append({"backbone": name, "num_samples": len(samples), "samples": samples})

    monomers = []
    mono_dir = os.path.join(seqs_dir, "monomers")
    if os.path.isdir(mono_dir):
        for fa in sorted(glob.glob(os.path.join(mono_dir, "*.fa")), key=lambda p: _natural_key(os.path.basename(p))):
            name = os.path.splitext(os.path.basename(fa))[0]
            records = _parse_fasta(fa)
            samples = [r for r in records if not r["is_input"]]
            monomers.append({"backbone": name, "num_samples": len(samples), "samples": samples})

    return {"backbones": backbones, "monomers": monomers}


# ---------------------------------------------------------------------------
# Modeling stage details
# ---------------------------------------------------------------------------

AF3_JOB_START_RE = re.compile(r"Running fold job (\S+)\.\.\.")
AF3_JOB_DONE_RE = re.compile(r"Fold job (\S+) done")
AF3_TIMING_RE = re.compile(
    r"Running model inference and extracting output structures with \d+ seed\(s\) took ([\d.]+) seconds"
)

BOLTZ_INVOKE_RE = re.compile(r"boltz predict\s+(\S+)\s")
BOLTZ_PROGRESS_RE = re.compile(
    r"Predicting DataLoader 0: 100%\|.*?\|\s*(\d+)/(\d+)\s*\[(\d+):(\d+)<"
)


def _af3_modeling_stats(lines: list[str]) -> dict:
    started: list[tuple[str, int]] = []
    done: set[str] = set()
    timings: list[float] = []
    for i, line in enumerate(lines):
        m = AF3_JOB_START_RE.search(line)
        if m:
            started.append((m.group(1), i))
        m = AF3_JOB_DONE_RE.search(line)
        if m:
            done.add(m.group(1))
        m = AF3_TIMING_RE.search(line)
        if m:
            timings.append(float(m.group(1)))

    completed = len(timings)
    current_name = None
    if started and started[-1][0] not in done:
        current_name = started[-1][0]

    avg = sum(timings) / len(timings) if timings else None
    return {"completed": completed, "current_name": current_name, "avg_seconds": avg, "unit": "sequence"}


def _boltz_modeling_stats(lines: list[str]) -> dict:
    invocations: list[tuple[str, int]] = []
    completions_after: dict[int, float] = {}
    timings: list[float] = []

    last_invoke_idx = None
    for i, line in enumerate(lines):
        m = BOLTZ_INVOKE_RE.search(line)
        if m:
            # pull the model_X folder name out of the yaml_inputs path
            path = m.group(1)
            mm = re.search(r"(model_\d+)", path)
            name = mm.group(1) if mm else path
            invocations.append((name, i))
            last_invoke_idx = i
        m = BOLTZ_PROGRESS_RE.search(line)
        if m and m.group(1) == m.group(2):  # only count fully-completed (N/N) lines
            minutes, seconds = int(m.group(3)), int(m.group(4))
            elapsed = minutes * 60 + seconds
            timings.append(elapsed)
            if last_invoke_idx is not None:
                completions_after[last_invoke_idx] = elapsed

    completed = len(completions_after)
    current_name = None
    if invocations:
        last_name, last_idx = invocations[-1]
        if last_idx not in completions_after:
            current_name = last_name

    avg = sum(timings) / len(timings) if timings else None
    return {"completed": completed, "current_name": current_name, "avg_seconds": avg, "unit": "backbone batch"}


def modeling_progress(lines: list[str], output_dir: str, cfg: dict, mp: dict) -> dict:
    """`mp` is the already-computed list_sequences(output_dir) result (the
    caller has it anyway for the mpnn summary), reused here instead of
    re-scanning the same directory.
    """
    prediction_model = (cfg.get("prediction_model") or "").upper()
    is_boltz = prediction_model.startswith("BOLTZ")
    model_monomer = bool(cfg.get("model_monomer"))

    # Ground truth: count actual completed structures on disk (one row per
    # confidence file, from list_models() - already tags monomer rows
    # correctly) rather than counting "took N seconds" log lines. This
    # naturally includes monomer predictions without any separate
    # accounting, and - importantly - updates the instant each individual
    # structure's confidence file is written, rather than only once a
    # whole batched AF3/Boltz subprocess call finishes and its buffered
    # stdout gets flushed to the log (which could easily look like "only
    # updates after the whole fasta file's worth of sequences are done").
    completed = len(list_models(output_dir, model_monomer=model_monomer))

    # Ground truth: exactly how many sequences were actually designed this
    # cycle - deliberately not num_seq_per_target_mpnn * num backbones,
    # since that config value is only honored in the first cycle; from
    # cycle 1 onwards MPNN is forced to design exactly one sequence per
    # input regardless of it (see list_sequences()'s note on cycle > 0
    # filenames), so counting actual fasta records sidesteps that entirely.
    total_designed = sum(len(bb["samples"]) for bb in mp["backbones"])

    monomer_multiplier = 2 if model_monomer else 1
    # Boltz produces num_models structures per sequence within a single
    # invocation (num_models doubles as its --diffusion_samples argument);
    # AF3 produces exactly one structure per sequence, so it has no such
    # multiplier - confirmed against real output directories for both.
    num_models = cfg.get("num_models") or 1
    per_sequence_multiplier = num_models if is_boltz else 1

    expected_total = (
        total_designed * per_sequence_multiplier * monomer_multiplier if total_designed else None
    )
    # NOTE: this - like `completed` above - is progress *within the
    # current cycle* only: 2_mpnn/ and 3_models/ both get wiped at the
    # start of the next cycle (see list_sequences() / list_models()), so
    # there's no way to observe earlier cycles' totals from disk once a
    # new cycle has started. This matches the Models tab, which shows the
    # same current-cycle-only scope.

    stats = _boltz_modeling_stats(lines) if is_boltz else _af3_modeling_stats(lines)

    eta_seconds = None
    if stats["avg_seconds"] is not None and expected_total is not None:
        remaining = max(expected_total - completed, 0)
        eta_seconds = stats["avg_seconds"] * remaining

    return {
        "prediction_model": cfg.get("prediction_model"),
        "completed": completed,
        "expected_total": expected_total,
        "unit": "structure",
        "current_name": stats["current_name"],
        "avg_seconds": stats["avg_seconds"],
        "eta_seconds": eta_seconds,
    }


def _pick_matching_file(candidates: list[str], conf_filename: str) -> Optional[str]:
    """Given a list of same-extension candidate filenames in a results
    folder, pick the one that actually matches this confidence file (a
    folder can hold more than one, e.g. Boltz writing one pdb per
    diffusion sample / num_models)."""
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]

    # Boltz: confidence_<name>.json <-> <name>.<ext>
    if conf_filename.startswith("confidence_"):
        expected_base = os.path.splitext(conf_filename[len("confidence_"):])[0]
        for p in candidates:
            if os.path.splitext(p)[0] == expected_base:
                return p

    # AF3: <name>_summary_confidences.json <-> <name>*.<ext>
    if conf_filename.endswith("summary_confidences.json"):
        prefix = conf_filename[: -len("_summary_confidences.json")]
        for p in candidates:
            if p.startswith(prefix):
                return p

    return candidates[0]  # fallback: best effort


def _match_structure_for_confidence(directory: str, conf_filename: str) -> tuple[Optional[str], Optional[str]]:
    """Pick the structure file that matches this confidence file, returning
    (path, format). Prefers .cif over .pdb: AF3 writes each model's .cif as
    soon as that model is done, but only writes the .pdb once *all* models
    for the job have finished - so preferring .cif lets the Models tab show
    a structure much sooner. Boltz only ever writes .pdb here, so this is a
    no-op for Boltz runs."""
    try:
        files = os.listdir(directory)
    except OSError:
        return None, None

    cif_pick = _pick_matching_file(sorted(f for f in files if f.endswith(".cif")), conf_filename)
    if cif_pick:
        return os.path.join(directory, cif_pick), "cif"

    pdb_pick = _pick_matching_file(sorted(f for f in files if f.endswith(".pdb")), conf_filename)
    if pdb_pick:
        return os.path.join(directory, pdb_pick), "pdb"

    return None, None


_CONF_MODEL_IDX_RE = re.compile(r"_model_(\d+)")

# AF3 writes one subfolder per seed/sample it modeled internally
# (seed-<n>_sample-<m>/), each with its own *_summary_confidences.json but
# only a .cif (no .pdb) - the top-level folder already has the
# representative pdb + confidences, so these would otherwise show up in the
# Models tab as entries with no structure to view. Skip them.
_SEED_SAMPLE_DIR_RE = re.compile(r"^seed-\d+_sample-\d+$")
_SKIP_DIR_NAMES = {"json_inputs", "yaml_inputs", "alignment_inputs"}


def list_models(output_dir: str, model_monomer: bool = False) -> list[dict]:
    """Ground-truth scan of 3_models/ on disk. Rather than assuming a
    fixed directory depth (which differs between AF3 and Boltz, and again
    between the multimer and monomers/ variants), this walks the whole
    model_X tree and treats every confidence-like file it finds
    (*summary_confidences.json for AF3, confidence_*.json for Boltz) as one
    row - which also means multiple structures per sequence (one per
    num_models) each get their own row instead of being collapsed.

    `model_monomer` (from the run config) decides whether non-monomer rows
    get tagged "complex" at all: if the run was never modeling a monomer
    side by side, there's nothing to contrast "complex" against, so no
    variant tag is set.
    """
    models_dir = os.path.join(output_dir, "3_models")
    out = []
    if not os.path.isdir(models_dir):
        return out

    # "model_N" (the 3_models/ subdirectory a row comes from) always
    # traces back to RFdiffusion's "_N.pdb"/"_N.trb" in 1_rfdiff/ - same
    # index, regardless of which sequence/sample/monomer-vs-complex
    # variant the row itself is, so the Models tab's viewer can offer the
    # same RFdiffusion-provenance coloring the Backbones/Results tabs do.
    def _trb_path_for_model(model_x: str) -> Optional[str]:
        m = re.match(r"model_(\d+)$", model_x)
        if not m:
            return None
        candidate = os.path.join(output_dir, "1_rfdiff", f"_{m.group(1)}.trb")
        return candidate if os.path.isfile(candidate) else None

    for model_x in sorted(os.listdir(models_dir), key=_natural_key):
        model_x_path = os.path.join(models_dir, model_x)
        if not os.path.isdir(model_x_path):
            continue
        trb_path = _trb_path_for_model(model_x)

        for root, dirs, files in os.walk(model_x_path):
            dirs[:] = [d for d in dirs if d not in _SKIP_DIR_NAMES and not _SEED_SAMPLE_DIR_RE.match(d)]
            conf_files = [
                f
                for f in files
                if f.endswith("summary_confidences.json") or (f.startswith("confidence_") and f.endswith(".json"))
            ]
            for cf in sorted(conf_files):
                conf_path = os.path.join(root, cf)
                structure_path, structure_format = _match_structure_for_confidence(root, cf)
                seq_name = os.path.basename(root)
                idx_m = _CONF_MODEL_IDX_RE.search(cf)
                is_monomer = "monomer" in root.lower()
                if is_monomer:
                    variant = "monomer"
                elif model_monomer:
                    variant = "complex"
                else:
                    variant = None
                out.append(
                    {
                        "model": model_x,
                        "sequence_name": seq_name,
                        "model_index": idx_m.group(1) if idx_m else None,
                        "variant": variant,
                        "structure_path": structure_path,
                        "structure_format": structure_format,
                        "confidence_path": conf_path,
                        "trb_path": trb_path,
                    }
                )

    out.sort(key=lambda r: (_natural_key(r["model"]), _natural_key(r["sequence_name"]), r["model_index"] or ""))
    return out


def load_confidence(path: str) -> dict:
    """Generic: surface every scalar field, plus small flat dicts/lists of
    scalars (e.g. AF3's chain_ptm, Boltz's chains_ptm - per-chain
    breakdowns are genuinely useful in a summary table). Skip anything
    that's a matrix / nested-deeper-than-one-level (e.g. chain_pair_iptm,
    pair_chains_iptm, per-residue pae arrays) since those aren't useful
    outside a dedicated matrix view."""
    try:
        with open(path) as f:
            data = json.load(f)
    except Exception as e:
        return {"error": str(e)}

    def is_scalar(x):
        return isinstance(x, (int, float, str, bool)) or x is None

    out = {}
    for k, v in data.items():
        if is_scalar(v):
            out[k] = v
        elif isinstance(v, list) and len(v) <= 8 and all(is_scalar(x) for x in v):
            out[k] = v
        elif isinstance(v, dict) and len(v) <= 8 and all(is_scalar(x) for x in v.values()):
            for sub_k, sub_v in v.items():
                out[f"{k}.{sub_k}"] = sub_v
        # else: skip (matrices, nested structures, etc.)
    return out


# ---------------------------------------------------------------------------
# Final results
# ---------------------------------------------------------------------------


def results_summary(output_dir: str) -> dict:
    final_dir = os.path.join(output_dir, "final_pdbs")
    csv_path = os.path.join(output_dir, "final_output.csv")
    pdbs = sorted(glob.glob(os.path.join(final_dir, "*.pdb"))) if os.path.isdir(final_dir) else []
    return {
        "final_pdbs": [os.path.basename(p) for p in pdbs],
        "csv_exists": os.path.isfile(csv_path),
        "csv_path": csv_path if os.path.isfile(csv_path) else None,
        "finished": bool(pdbs) and os.path.isfile(csv_path),
    }


# ---------------------------------------------------------------------------
# Top level: build the whole status payload
# ---------------------------------------------------------------------------


def get_status(log_path: str) -> dict:
    if not os.path.isfile(log_path):
        return {"error": f"Log file not found: {log_path}"}

    loc = resolve_output_dir(log_path)
    if loc.get("error") or not loc.get("output_dir"):
        return {"error": loc.get("error") or "Could not resolve output_dir", "location": loc}

    output_dir = loc["output_dir"]

    with open(log_path, "r", errors="replace") as f:
        lines = f.readlines()

    cfg = extract_config(log_path, output_dir)
    stage_info = detect_stage(lines)
    stage = stage_info["stage"]

    err_loc = locate_err_file(log_path)
    cancel_info = scan_err_file_for_cancellation(err_loc["err_path"]) if err_loc["err_exists"] else {
        "cancelled": False,
        "cancelled_line": None,
    }

    payload: dict[str, Any] = {
        "location": loc,
        "config": cfg,
        "stage": stage,
        "possible_errors": stage_info["possible_errors"],
        "output_dir_exists": os.path.isdir(output_dir),
        "job_info": parse_job_info(lines),
        "timing": parse_step_durations(lines),
        "cycle": current_cycle_info(lines, cfg),
        "crash": detect_crash(lines),
        "err_file": err_loc,
        "cancelled": cancel_info["cancelled"],
        "cancelled_line": cancel_info["cancelled_line"],
    }

    if not payload["output_dir_exists"]:
        return payload

    # Always include lightweight counts so the stepper / nav can show
    # progress hints even when not the active stage.
    payload["rfdiffusion"] = rfdiffusion_progress(lines, output_dir, cfg.get("num_designs_rfdiff"))

    # Ground truth for how many backbones actually survived plugin_filters
    # (as opposed to num_designs_rfdiff, which is just the planned count).
    # Blank (None) until filtering has actually finished - see
    # list_backbones() for why a backbone isn't tagged "passed" before then.
    filtering_done = stage not in ("setup", "rfdiffusion", "filtering")
    backbones = list_backbones(output_dir, filtering_done=filtering_done)
    payload["backbones_summary"] = {
        "total": len(backbones),
        "accepted": sum(1 for b in backbones if b["status"] == "passed") if filtering_done else None,
    }

    mp = list_sequences(output_dir)
    payload["mpnn"] = {
        "num_backbones_with_sequences": len(mp["backbones"]),
        "num_monomer_backbones": len(mp["monomers"]),
    }

    payload["modeling"] = modeling_progress(lines, output_dir, cfg, mp)
    payload["scoring"] = results_summary(output_dir)

    return payload
