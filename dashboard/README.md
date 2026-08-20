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
from the log + output directory on every refresh. Auto-refresh stops itself
(and unchecks the box, with a small note explaining why) once the job
reaches a terminal state — finished, crashed, or cancelled — since there's
nothing left to watch for; loading a different job (or the same one again)
always re-enables it. Selections (the structure shown in a viewer, which
list item is highlighted, list scroll position, the "Show full
configuration" disclosure) are preserved across a refresh rather than being
reset by it - a poll only touches what actually changed on disk.

A theme toggle (🌙/☀️) next to auto-refresh switches between a dark and a
light theme; it defaults to the OS's `prefers-color-scheme` and remembers
an explicit choice in local storage. The 3D structure viewers, the sequence
alignment viewer and the Error log's `<pre>` intentionally stay dark in
both themes, the way a code/terminal panel commonly does in an
otherwise-light editor.

Every 3D structure viewer (Backbones, Models, Results) is clickable:
clicking a residue highlights it and everything within 5 Å of it as
licorice sticks (both colored by element - CPK convention, N blue/O
red/S yellow/etc - with the clicked residue drawn thicker than its
neighbors rather than in a different flat color, so the element coloring
stays meaningful on both), on top of the usual cartoon, and shows a small
overlay in the corner of the viewer naming the residue plus how many
neighbors it found - click "clear" on that overlay, or click another
residue, to change it. Hovering (without clicking) shows a lightweight
floating tooltip naming whatever residue is under the cursor, without
touching the highlight. Both are built on
[NGL Viewer](https://nglviewer.org/)'s own picking API
(`stage.signals.clicked`/`hovered`) and its `getAtomSetWithinSelection` /
`getAtomSetWithinGroup` distance query for the neighbor search - see
`setupResidueInteraction()` in `app.js`. (An earlier version of this used
3Dmol.js; its residue-picking turned out to be unreliable in practice, so
the dashboard moved to NGL, which is more purpose-built for exactly this
interaction.) NGL's own default hover tooltip is detached from the DOM
right after each structure loads, since the dashboard supplies its own
app-styled one instead.

Every structure viewer (Backbones, Models, Results) has its own small
color toolbar above it, all driving the same settings (remembered in
local storage) so they stay in sync with each other no matter which tab
you're on: a **Background** select (a few preset colors, top-right of the
toolbar, always available) applies instantly via a live stage parameter;
a **Chain colors** select, shown whenever a viewer is coloring by chain,
picks between a few curated palettes (Vivid, Pastel, a colorblind-safe
Okabe-Ito set) - chosen deliberately over NGL's built-in `chainid`
scheme, which hashes the chain letter into a color that can land on tones
too dark to read against the (also dark) default background - paired
with a small legend naming which color is which chain. Changing either
setting recolors every currently-loaded viewer in place (swap the
cartoon representation for one using a new NGL color scheme via
`setCartoonColor()`) rather than reloading the structure, so it doesn't
reset your camera zoom/rotation.

The Backbones and Results tabs additionally have a **Color residues by**
/ **Color structure by** dropdown, offering **RFdiffusion provenance** as
an alternative to chain coloring. RFdiffusion writes a `.trb` sidecar
file (a Python pickle, despite the `.trb` extension -
`load_trb_provenance()` in `parser.py` unpickles it, which needs `numpy`
since the pickle contains `numpy` array/scalar objects) next to each
backbone `.pdb`, recording which residues were copied in from the
reference structure versus generated de novo. The dashboard buckets every
residue into one of three categories, colored consistently across both
tabs and shown in a small legend next to the dropdown, each with its own
color picker so the palette can be retuned: **Motif** (`con_hal_pdb_idx`
- reference residues present in a *redesigned* chain), **Fixed chains**
(`receptor_con_hal_pdb_idx` - reference residues in a *non-designed*
chain), and **Sculpted** (everything else - generated de novo). Residue
identity for this is chain letter + residue number, but the two files
this needs to apply to don't agree on numbering convention: RFdiffusion's
own raw backbone `.pdb` (Backbones tab) keeps counting residue numbers up
across chain boundaries instead of resetting per chain, while
AlphaFold3/Boltz's output (Results tab, via that row's `path_rfdiff`
column) resets every chain back to 1 like a normally-numbered PDB file
does. `load_trb_provenance()` returns `fixed_chain` residue numbers
already shifted onto the chain-local convention (each chain's lowest
resnum becomes 1, then 2, ...) - `con_hal_pdb_idx`/`motif` doesn't need
the same treatment since the designed chain(s) it refers to always come
first and so are never offset either way - and the frontend
(`computeStructureChainInfo()` in `app.js`) works out each *loaded*
structure's own per-chain offset the same way and normalizes onto it, so
the same provenance data colors correctly regardless of which
convention the structure in front of it happens to use. If a backbone
has no `.trb` file, or it can't be read (e.g. `numpy` isn't installed in
the dashboard's own environment), the dropdown falls back to chain
coloring with a short explanation in place of the legend rather than
failing silently.

## What each tab shows

- **Overview** — a stepper showing the current pipeline stage
  (RFdiffusion → Filtering → MPNN → Modeling → Scoring → Finished), with a
  status card directly underneath it showing what's happening right now:
  for RFdiffusion, how many backbones are done, the average time per
  backbone, and an ETA; for Modeling, one "structure" = one predicted
  model (one designed sequence, or one of Boltz's diffusion samples for a
  sequence), counted straight off disk as each one appears - including
  monomer predictions, which double the total when `model_monomer: true` -
  rather than from "took N seconds" log lines, which can lag well behind
  what's actually finished if AF3/Boltz's own stdout only gets flushed in
  batches; a short status line for the quicker stages. Below that: the
  Slurm job id, node and start time (parsed from the scheduler's "Hello
  from job ..." banner line); a **run plan** card with the planned
  RFdiffusion backbone count, the count that actually survived filtering
  ("Accepted backbones" - blank until filtering has finished), sequences
  per backbone, and (only shown when not the trivial value of 1) models
  per sequence and cycle count, plus "Planned designs" (backbones ×
  sequences/backbone) and "Final designs" (accepted backbones ×
  sequences/backbone, also blank until filtering has finished); a **step
  timings** card that fills in as each pipeline step finishes (RFdiffusion,
  filtering, MPNN+modeling, final scoring, ...) plus the total job
  duration once it's known; and — for multi-cycle runs (`af2_mpnn_cycles`
  > 1) — which cycle is currently running (note that the Modeling
  progress numbers reflect the current cycle only, the same as the Models
  tab, since disk doesn't retain earlier cycles' output). Any lines in the
  log that look like tracebacks/errors are
  surfaced here too. If the log contains Prosculpt's fatal-error line
  ("There was an error running the command.") a crash banner appears, and
  separately, if the slurm `.err` file's tail contains a Slurm cancellation
  notice ("*** JOB ...
  CANCELLED ...") a cancellation banner appears — both link to the
  **Error log** tab.
- **Error log** — the dashboard looks for the job's `.err` file (same
  directory as the log, same name, `.err` extension — the Slurm stderr
  convention) unconditionally from the very first status fetch, not only
  once something looks wrong. If found, an "Error log" stat on the
  Overview card says so and this tab is available to show its full
  contents (highlighted red in the tab bar, with a banner recapping why,
  if a crash or cancellation was actually detected); if not found, the
  Overview stat just says "not found" and the tab stays hidden.
- **Output log** — the job's own `.out` log, shown in full (same size cap
  as the Error log tab). Always available, since it's the same file the
  dashboard already needs to operate. If you're at the bottom of it when a
  refresh comes in, it keeps following the end (like `tail -f`); scroll up
  to read something and refreshes stop yanking you back down.
- **Backbones** — every RFdiffusion backbone produced so far (including
  ones removed by filters, tagged accordingly), with a 3D structure viewer.
  A backbone isn't tagged "ok" until `plugin_filters` has actually finished
  (tracked via the log's stage detection) — before that it's tagged
  "pending" rather than implying a filter verdict that hasn't happened yet.
  A "Show filtered out" toggle above the list hides the filtered-out ones.
  The viewer's "Color residues by" dropdown can switch to RFdiffusion
  provenance coloring (Motif/Fixed chains/Sculpted, see above) when a
  `.trb` sidecar is present.
- **Sequences** — the ProteinMPNN output, grouped by backbone, rendered in
  a shared alignment viewer: all samples for a backbone scroll together
  under one horizontal scrollbar (with a sticky label column and a residue
  ruler), rather than each sequence having its own independent scrollbar.
  Residues are colored by amino acid type. A "Highlight conserved
  residues" toggle (with a % threshold field) outlines residues that
  match the per-position consensus in at least that fraction of a
  backbone's samples — a sample whose residue differs from the consensus
  at a given position is never highlighted there, even if that position
  is otherwise highly conserved.
- **Models** — every modeled structure produced so far (one row per
  sequence **and** per `num_models` sample — Boltz especially can produce
  several structures per sequence), with a structure viewer and a metrics
  table read straight from AlphaFold3's `summary_confidences.json` or
  Boltz's `confidence_*.json`. AF3 writes each model's `.cif` as soon as
  that model is done but only writes the `.pdb` once *all* models in the
  job have finished, so the viewer prefers `.cif` when present (labelled
  as such) so structures show up as soon as possible instead of waiting
  for the whole job; Boltz only ever writes `.pdb`, so this has no effect
  there. AF3's internal per-seed/per-sample subfolders
  (`seed-<n>_sample-<m>/`, which have confidences but no structure) are
  skipped. Rows are only tagged "complex" when the run actually also
  modeled a monomer side by side (`model_monomer: true`); otherwise they're
  left untagged rather than mislabeled. On a multi-cycle run this tab shows
  only the current cycle's models, since `3_models/` is wiped by Prosculpt
  at the start of each new cycle — a note above the list says which cycle.
- **Results** — once `final_pdbs/` and `final_output.csv` exist: a
  structure viewer, a sequence alignment view (grouped by sequence length),
  and a table of `final_output.csv`. The table pins the `id` column on the
  left, pushes path-like columns (`model_path`, `af3_json`, `af3_pdb`,
  `path_rfdiff`) to the right end, and stripes rows for readability. A
  **Filters** panel lets you set a min/max cutoff on any numeric column,
  which hides out-of-range rows from the list, the alignment view and the
  table together. A **Color list by** control tints the model list and
  table rows on a red→green scale for a chosen numeric column, with a
  toggle for whether higher or lower is better - kept as a separate
  control from the structure viewer's own "Color structure by" dropdown
  (chain vs. RFdiffusion provenance, see above), since the two color
  different things at once. An **Export filtered** button
  downloads a `.zip` (not `.rar` — that needs a proprietary external
  binary that isn't reliably available on a cluster, so `.zip` is used
  instead; no extra dependency, opens natively everywhere) containing the
  currently-filtered rows as `filtered_output.csv` plus each of their
  `model_path` pdb files.

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
- **Run configuration is read generically too.** `extract_config()` in
  `parser.py` doesn't hard-code which config keys to look for; it parses
  the fully-resolved Hydra config Prosculpt logs near the top of every run
  ("The following configuration was passed: ..." followed by a plain YAML
  dump) and surfaces every key from it, falling back to `input.yaml` in
  the output directory for anything the log block is missing. So
  task-specific options (e.g. `boltz2_templates`, `symmetry`, `ppi`) show
  up automatically without any dashboard changes.
- **Multi-cycle runs (`af2_mpnn_cycles` > 1)** have been exercised against
  real logs. The current cycle is read from the `cycleeeeee <n>` marker
  do_cycling() prints. The Modeling stage's "completed"/"expected total"
  are both scoped to the *current* cycle only (backbones and models from
  earlier cycles aren't observable on disk once the next cycle starts,
  since `2_mpnn/` and `3_models/` both get wiped between cycles) - matching
  what the Models tab already shows. "expected total" is derived from the
  actual fasta records MPNN wrote this cycle (`num_seq_per_target_mpnn` is
  only honored in cycle 0; from cycle 1 onwards MPNN is forced to one
  sequence per input, so counting real records sidesteps needing to know
  which applies) times a per-sequence structure count (1 for AF3, `num_models`
  for Boltz - it doubles as Boltz's `--diffusion_samples` argument, unlike
  AF3 which produces exactly one structure per sequence) times 2 if
  `model_monomer: true`. avg-time/ETA still come from the log's timing
  lines (blended across whatever cycles are in the log), since that rate
  info genuinely isn't available from the filesystem. Note that from
  cycle 1 onwards, MPNN output fasta files are no longer named `_<n>.fa`
  but `rf_<n>__model_<m>__cycle_<c>__itr_<i>__.fa`; the dashboard's globs
  account for this.
- **Boltz vs AF3 detection** is based on the `prediction_model` config
  value (`AF3` / `Boltz2`) read from the log/`input.yaml`.

## Project layout

```
app.py            Flask routes / API
parser.py          All log-parsing + filesystem-scanning logic (no Flask
                    dependency, easy to test standalone)
templates/index.html
static/css/style.css
static/js/app.js   Frontend (vanilla JS, no build step; uses NGL Viewer
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
- **Error log tab doesn't appear** — it's only shown once a `.err` file is
  actually found next to the log (same basename, `.err` extension); the
  Overview card's "Error log" stat says "not found" in that case. A job
  that's merely slow or stalled, with no crash/cancellation, will still
  show the tab (if the `.err` file exists) but without a warning banner.
