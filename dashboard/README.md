# Prosculpt Dashboard

A small local web dashboard for monitoring a running (or finished) Prosculpt
job: which stage it's at, progress/ETA for RFdiffusion and modeling, plus
browsable views of the backbones, designed sequences, individual models,
and final results.

It works by reading the job's log file and its output directory directly
off disk, so **it needs to run somewhere that has filesystem access to
both** — typically the same machine (e.g. a cluster login node) where you
launched the Prosculpt job.

## 1. Install

Needs Python 3.9+.

```bash
cd prosculpt_dashboard
pip install -r requirements.txt
```

## 2. Run

```bash
python app.py
```

This starts a local server on port 5000. Open **http://localhost:5000** in
your browser.

If you're running this on a remote machine (e.g. a cluster login node) and
browsing from your laptop, forward the port over SSH instead of exposing it
publicly:

```bash
ssh -L 5000:localhost:5000 you@cluster
```

then open http://localhost:5000 on your laptop as usual.

To use a different port: `PORT=8080 python app.py`.

## 3. Point it at a job

In the top bar, paste the full path to the job's log file (the `.out`/`.log`
file Slurm or your scheduler wrote), or click **Browse…** to navigate the
filesystem from the dashboard itself. Click **Load**.

The dashboard reconstructs the job's output directory the same way you
described: it reads the `PWD:` line (working directory at launch) and the
`output_dir:` line from the log, and joins them if `output_dir` is relative.

The page polls every 5 seconds for updates (toggle off with the
"auto-refresh" checkbox next to the log path field). It's safe to point it
at a log file for a job that's still running — everything is computed fresh
from the log + output directory on every refresh.

## What each tab shows

- **Overview** — a stepper showing the current pipeline stage
  (RFdiffusion → Filtering → MPNN → Modeling → Scoring → Finished), plus
  stage-specific detail: for RFdiffusion, how many backbones are done, the
  average time per backbone, and an ETA; for Modeling, the same but per
  designed sequence (AlphaFold3) or per RFdiffusion-backbone batch (Boltz,
  since Boltz models all sequences for one backbone in a single call). Any
  lines in the log that look like tracebacks/errors are surfaced here too.
- **Backbones** — every RFdiffusion backbone produced so far (including
  ones removed by filters, tagged accordingly), with a 3D structure viewer.
- **Sequences** — the ProteinMPNN output, grouped by backbone. Sequences
  from the same backbone are the same length, so they line up naturally
  without needing a real alignment algorithm; residues are colored by
  amino acid type.
- **Models** — every modeled structure produced so far (one row per
  sequence **and** per `num_models` sample — Boltz especially can produce
  several structures per sequence), with a structure viewer and a metrics
  table read straight from AlphaFold3's `summary_confidences.json` or
  Boltz's `confidence_*.json`.
- **Results** — once `final_pdbs/` and `final_output.csv` exist: a
  structure viewer over the final models and a sortable table of
  `final_output.csv`.

## Design notes / known limitations

- **Filesystem is ground truth, the log gives rates.** Counts and listings
  (how many backbones exist, which sequences/models exist) are read
  straight off disk, not inferred from the log — this keeps things correct
  even if filtering, retries, or crashes make the log's bookkeeping
  diverge from reality. The log is only used to compute *timing*: how long
  each stage / backbone / model took, and hence the ETA.
- **Confidence metrics are read generically.** Rather than hard-coding
  AlphaFold3's or Boltz's exact JSON schema (which could change), the
  dashboard surfaces every scalar field it finds, plus small per-chain
  breakdowns (e.g. `chain_ptm`, `chains_ptm.0`), and skips large matrices
  (e.g. `chain_pair_iptm`, per-residue PAE) that aren't useful in a
  summary table. If a tool changes its field names, this keeps working;
  if you want a specific field pinned/highlighted, it's a small edit to
  `load_confidence()` in `parser.py`.
- **Multi-cycle runs (`af2_mpnn_cycles` > 1).** The example logs used to
  build this only exercised a single cycle. The stage detector and
  progress counters are written generically enough that they should
  continue to work (each new cycle just re-triggers the MPNN → Modeling
  detection), but the "expected total" estimate for the Modeling stage's
  ETA multiplies by `af2_mpnn_cycles` from the config rather than being
  derived from direct observation — worth double-checking against a real
  multi-cycle run.
- **Boltz vs AF3 detection** is based on the `prediction_model` config
  value (`AF3` / `Boltz2`) read from the log/`input.yaml`.

## Project layout

```
app.py            Flask routes / API
parser.py          All log-parsing + filesystem-scanning logic (no Flask
                    dependency, easy to test standalone)
templates/index.html
static/css/style.css
static/js/app.js   Frontend (vanilla JS, no build step; uses 3Dmol.js
                    from a CDN for structure viewing)
requirements.txt
```

## Troubleshooting

- **"Could not find PWD: and/or output_dir: in the log file"** — the
  dashboard expects the exact log format Prosculpt writes (a `PWD:` line
  followed by the working directory on the next line, and a bare
  `output_dir: <path>` line from the printed config). If your log format
  differs, check `resolve_output_dir()` in `parser.py`.
- **A tab shows "No X generated yet" for something you know exists** —
  the dashboard looks for very specific file patterns (e.g. `_*.pdb` in
  `1_rfdiff/`, `*summary_confidences.json` / `confidence_*.json` under
  `3_models/`). If Prosculpt's output layout has changed, these glob
  patterns in `parser.py` are the place to update.
