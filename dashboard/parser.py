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
gets renamed/added/removed in a future tool version.
"""

from __future__ import annotations

import glob
import json
import os
import re
from datetime import datetime
from typing import Any, Optional

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

    return {"pwd": pwd, "output_dir_raw": output_dir_rel, "output_dir": full, "error": None}


# ---------------------------------------------------------------------------
# Config extraction (from the yaml dump at the top of the log, with a
# fallback to reading input.yaml directly off disk once it exists)
# ---------------------------------------------------------------------------

CONFIG_KEYS = {
    "num_designs_rfdiff": int,
    "num_seq_per_target_mpnn": int,
    "num_models": int,
    "af2_mpnn_cycles": int,
    "model_monomer": lambda v: str(v).strip().lower() == "true",
    "prediction_model": str,
    "task_name": str,
    "contig": str,
}


def _coerce(key, raw):
    caster = CONFIG_KEYS.get(key, str)
    try:
        return caster(raw)
    except Exception:
        return raw


def extract_config(log_path: str, output_dir: Optional[str]) -> dict:
    cfg: dict[str, Any] = {}
    with open(log_path, "r", errors="replace") as f:
        for line in f:
            stripped = line.strip()
            for key in CONFIG_KEYS:
                if key in cfg:
                    continue
                m = re.match(rf"^{re.escape(key)}:\s*(\S.*)$", stripped)
                if m:
                    cfg[key] = _coerce(key, m.group(1).strip().strip("'\""))
            if len(cfg) == len(CONFIG_KEYS):
                break

    # Fallback / cross-check: read input.yaml directly if present, filling
    # in anything the log-scrape missed.
    if output_dir and yaml is not None:
        yaml_path = os.path.join(output_dir, "input.yaml")
        if os.path.isfile(yaml_path):
            try:
                with open(yaml_path) as f:
                    doc = yaml.safe_load(f) or {}
                for key in CONFIG_KEYS:
                    if key not in cfg and key in doc:
                        cfg[key] = doc[key]
            except Exception:
                pass

    cfg.setdefault("num_designs_rfdiff", None)
    cfg.setdefault("num_seq_per_target_mpnn", None)
    cfg.setdefault("num_models", None)
    cfg.setdefault("af2_mpnn_cycles", 1)
    cfg.setdefault("model_monomer", False)
    cfg.setdefault("prediction_model", None)
    return cfg


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


def list_backbones(output_dir: str) -> list[dict]:
    rfdiff_dir = os.path.join(output_dir, "1_rfdiff")
    out = []
    if not os.path.isdir(rfdiff_dir):
        return out
    for pdb in sorted(glob.glob(os.path.join(rfdiff_dir, "_*.pdb"))):
        out.append({"name": os.path.splitext(os.path.basename(pdb))[0], "status": "passed", "path": pdb})
    failed_dir = os.path.join(rfdiff_dir, "failed_filters")
    if os.path.isdir(failed_dir):
        for pdb in sorted(glob.glob(os.path.join(failed_dir, "_*.pdb"))):
            out.append({"name": os.path.splitext(os.path.basename(pdb))[0], "status": "failed_filter", "path": pdb})
    out.sort(key=lambda r: _natural_key(r["name"]))
    return out


def _natural_key(s: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


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
        for fa in sorted(glob.glob(os.path.join(seqs_dir, "_*.fa")), key=lambda p: _natural_key(os.path.basename(p))):
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


def modeling_progress(lines: list[str], output_dir: str, cfg: dict) -> dict:
    prediction_model = (cfg.get("prediction_model") or "").upper()
    is_boltz = prediction_model.startswith("BOLTZ")

    stats = _boltz_modeling_stats(lines) if is_boltz else _af3_modeling_stats(lines)

    # Ground truth: how many backbones actually survived to the mpnn stage.
    seqs_dir = os.path.join(output_dir, "2_mpnn", "seqs")
    num_backbones = len(glob.glob(os.path.join(seqs_dir, "_*.fa")))
    monomer_multiplier = 2 if cfg.get("model_monomer") else 1
    cycles = cfg.get("af2_mpnn_cycles") or 1
    num_seq_per_target = cfg.get("num_seq_per_target_mpnn") or 1

    if is_boltz:
        expected_total = num_backbones * monomer_multiplier * cycles if num_backbones else None
    else:
        expected_total = (
            num_backbones * num_seq_per_target * monomer_multiplier * cycles if num_backbones else None
        )

    eta_seconds = None
    if stats["avg_seconds"] is not None and expected_total is not None:
        remaining = max(expected_total - stats["completed"], 0)
        eta_seconds = stats["avg_seconds"] * remaining

    return {
        "prediction_model": cfg.get("prediction_model"),
        "completed": stats["completed"],
        "expected_total": expected_total,
        "unit": stats["unit"],
        "current_name": stats["current_name"],
        "avg_seconds": stats["avg_seconds"],
        "eta_seconds": eta_seconds,
    }


def _match_pdb_for_confidence(directory: str, conf_filename: str) -> Optional[str]:
    """A results folder can contain more than one pdb (Boltz writes one per
    diffusion sample / num_models). Try to pick the one that actually
    matches this particular confidence file rather than just grabbing the
    first pdb we see."""
    try:
        pdbs = sorted(f for f in os.listdir(directory) if f.endswith(".pdb"))
    except OSError:
        return None
    if not pdbs:
        return None
    if len(pdbs) == 1:
        return os.path.join(directory, pdbs[0])

    # Boltz: confidence_<name>.json <-> <name>.pdb
    if conf_filename.startswith("confidence_"):
        expected = conf_filename[len("confidence_"):]
        expected = os.path.splitext(expected)[0] + ".pdb"
        if expected in pdbs:
            return os.path.join(directory, expected)

    # AF3: <name>_summary_confidences.json <-> <name>*.pdb
    if conf_filename.endswith("summary_confidences.json"):
        prefix = conf_filename[: -len("_summary_confidences.json")]
        for p in pdbs:
            if p.startswith(prefix):
                return os.path.join(directory, p)

    return os.path.join(directory, pdbs[0])  # fallback: best effort


_CONF_MODEL_IDX_RE = re.compile(r"_model_(\d+)")


def list_models(output_dir: str) -> list[dict]:
    """Ground-truth scan of 3_models/ on disk. Rather than assuming a
    fixed directory depth (which differs between AF3 and Boltz, and again
    between the multimer and monomers/ variants), this walks the whole
    model_X tree and treats every confidence-like file it finds
    (*summary_confidences.json for AF3, confidence_*.json for Boltz) as one
    row - which also means multiple structures per sequence (one per
    num_models) each get their own row instead of being collapsed.
    """
    models_dir = os.path.join(output_dir, "3_models")
    out = []
    if not os.path.isdir(models_dir):
        return out

    for model_x in sorted(os.listdir(models_dir), key=_natural_key):
        model_x_path = os.path.join(models_dir, model_x)
        if not os.path.isdir(model_x_path):
            continue

        for root, dirs, files in os.walk(model_x_path):
            dirs[:] = [d for d in dirs if d not in ("json_inputs", "yaml_inputs")]
            conf_files = [
                f
                for f in files
                if f.endswith("summary_confidences.json") or (f.startswith("confidence_") and f.endswith(".json"))
            ]
            for cf in sorted(conf_files):
                conf_path = os.path.join(root, cf)
                pdb_path = _match_pdb_for_confidence(root, cf)
                seq_name = os.path.basename(root)
                idx_m = _CONF_MODEL_IDX_RE.search(cf)
                is_monomer = "monomer" in root.lower()
                out.append(
                    {
                        "model": model_x,
                        "sequence_name": seq_name,
                        "model_index": idx_m.group(1) if idx_m else None,
                        "variant": "monomer" if is_monomer else "complex",
                        "pdb_path": pdb_path,
                        "confidence_path": conf_path,
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

    payload: dict[str, Any] = {
        "location": loc,
        "config": cfg,
        "stage": stage,
        "possible_errors": stage_info["possible_errors"],
        "output_dir_exists": os.path.isdir(output_dir),
    }

    if not payload["output_dir_exists"]:
        return payload

    # Always include lightweight counts so the stepper / nav can show
    # progress hints even when not the active stage.
    payload["rfdiffusion"] = rfdiffusion_progress(lines, output_dir, cfg.get("num_designs_rfdiff"))

    mp = list_sequences(output_dir)
    payload["mpnn"] = {
        "num_backbones_with_sequences": len(mp["backbones"]),
        "num_monomer_backbones": len(mp["monomers"]),
    }

    payload["modeling"] = modeling_progress(lines, output_dir, cfg)
    payload["scoring"] = results_summary(output_dir)

    return payload
