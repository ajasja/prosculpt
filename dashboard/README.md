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
filesystem from the dashboard itself. Click **Add job**.

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
licorice sticks (colored by element - CPK convention, N blue/O red/S
yellow/etc - except the clicked residue's own carbons, which turn a
vivid yellow so it's unmistakable which one is actually selected among
its neighbors), on top of the usual cartoon, and shows a small overlay in
the corner of the viewer naming the residue plus how many neighbors it
found - click "clear" on that overlay, or click another residue, to
change it. Click-dragging over a range, or Ctrl/Cmd+clicking several
residues one at a time, selects multiple residues at once - all of them
get the vivid-yellow-carbon treatment, plus their combined surroundings
within 5 Å, exactly like a single click. Hovering (without clicking)
shows a lightweight floating tooltip naming whatever residue is under
the cursor, without touching the highlight. All of this is built on
[NGL Viewer](https://nglviewer.org/)'s own picking API
(`stage.signals.clicked`/`hovered`) and its `getAtomSetWithinSelection` /
`getAtomSetWithinGroup` distance query for the neighbor search - see
`setupResidueInteraction()` in `app.js`. (An earlier version of this used
3Dmol.js; its residue-picking turned out to be unreliable in practice, so
the dashboard moved to NGL, which is more purpose-built for exactly this
interaction.) NGL's own default hover tooltip is detached from the DOM
right after each structure loads, since the dashboard supplies its own
app-styled one instead.

The Backbones, Models and Results tabs each also show the currently
selected structure's sequence in an interactive panel below the viewer,
read straight off the already-loaded NGL structure (not a separately
fetched FASTA string, so a residue's position here can never drift out
of sync with the same residue in the 3D view next to it). Clicking a
residue in the sequence panel selects it exactly as a direct 3D click
would; dragging across a range or Ctrl/Cmd+clicking several residues
selects all of them, the same as doing it in the 3D view - the sequence
panel and the 3D viewer's own click/drag selection are two ways of
driving the same underlying highlight state (`buildStructureSequence()`,
`renderSequenceResidues()`, `setupSequencePanelInteraction()` in
`app.js`).

Every structure viewer (Backbones, Models, Results) has its own small
color toolbar above it, all driving the same settings (remembered in
local storage) so they stay in sync with each other no matter which tab
you're on: a **Background** select (a few preset colors, top-right of the
toolbar, always available) applies instantly via a live stage parameter;
a **Palette** select, shown whenever a viewer is coloring by chain, picks
between several curated palettes (Vivid, Pastel, a colorblind-safe
Okabe-Ito set, Bright, Sunset, Ocean) - chosen deliberately over NGL's
built-in `chainid` scheme, which hashes the chain letter into a color
that can land on tones too dark to read against the (also dark) default
background - paired with a small legend naming which color is which
chain, where each chain's swatch is itself a color `<input>` for manually
overriding just that chain. Manually recoloring a chain flips the
Palette select to a disabled "Custom" option (so it's clear a palette no
longer describes what's on screen); picking an actual palette afterward
clears every manual override and goes back to driving the colors
directly. Changing background or palette recolors every currently-loaded
viewer in place (swap the cartoon representation for one using a new NGL
color scheme via `setCartoonColor()`) rather than reloading the
structure, so it doesn't reset your camera zoom/rotation.

A **Show sidechains** control (three mutually-exclusive options: "Only
selected", "Interface", "All") governs which residues' sidechains are
drawn as licorice sticks on top of the cartoon, independent of whichever
residue is currently click-highlighted. "Only selected" (the default)
shows nothing extra beyond whatever the click/sequence-panel highlight
already draws. "Interface" shows the sidechains of every residue that has
any atom within 5 Å of an atom belonging to a *different* chain -
computed per chain via `structure.getAtomSet()` intersected with
`getAtomSetWithinSelection("not :chain", 5)`, then expanded to whole
residues and OR'd across chains (`computeInterfaceSele()` in `app.js`) -
useful for eyeballing a binder's interface at a glance without having to
click through it residue by residue. "All" shows every sidechain in the
structure. This setting, like the background and palette ones, is shared
across every viewer and remembered in local storage.

The Backbones, Models and Results tabs additionally have a **Color
residues by** / **Color structure by** dropdown, offering **RFdiffusion
provenance** as an alternative to chain coloring. RFdiffusion writes a
`.trb` sidecar file (a Python pickle, despite the `.trb` extension -
`load_trb_provenance()` in `parser.py` unpickles it, which needs `numpy`
since the pickle contains `numpy` array/scalar objects) next to each
backbone `.pdb`, recording which residues were copied in from the
reference structure versus generated de novo. The dashboard buckets every
residue into one of three categories, colored consistently across all
three tabs and shown in a small legend next to the dropdown, each with
its own color picker so the palette can be retuned: **Motif**
(`con_hal_pdb_idx` - reference residues present in a *redesigned* chain),
**Fixed chains** (`receptor_con_hal_pdb_idx` - reference residues in a
*non-designed* chain), and **Sculpted** (everything else - generated de
novo). Residue identity for this is chain letter + residue number, but
the files this needs to apply to don't agree on numbering convention:
RFdiffusion's own raw backbone `.pdb` (Backbones tab) keeps counting
residue numbers up across chain boundaries instead of resetting per
chain, while AlphaFold3/Boltz's output (Models and Results tabs) resets
every chain back to 1 like a normally-numbered PDB file does. Every
`3_models/model_N` directory traces back to the same-indexed
`1_rfdiff/_N.trb` regardless of which sequence/sample/monomer-vs-complex
variant a given row is (`list_models()` in `parser.py` resolves this and
stamps each row with its `trb_path`), and a Results row finds its `.trb`
the same way the Backbones tab does - same basename as that row's
`path_rfdiff`, `.trb` extension. `load_trb_provenance()` returns
`fixed_chain` residue numbers already shifted onto the chain-local
convention (each chain's lowest resnum becomes 1, then 2, ...) -
`con_hal_pdb_idx`/`motif` doesn't need the same treatment since the
designed chain(s) it refers to always come first and so are never offset
either way - and the frontend (`computeStructureChainInfo()` in
`app.js`) works out each *loaded* structure's own per-chain offset the
same way and normalizes onto it, so the same provenance data colors
correctly regardless of which convention the structure in front of it
happens to use. If a backbone/model has no `.trb` file, or it can't be
read (e.g. `numpy` isn't installed in the dashboard's own environment),
the dropdown falls back to chain coloring with a short explanation in
place of the legend rather than failing silently.

## 4. Track multiple jobs at once

The **Add job** field (and the **Browse…** button, which now lets you tick
several files before committing with "Add selected jobs") aren't limited
to one job - paste several log paths at once (one per line - it's a
`<textarea>` styled to look like a single-line field, since a real
`<input type="text">` can't actually hold a pasted newline) to track them
all. A pill bar appears above the usual tab strip once you do: **All jobs
overview** and **All jobs results** on the left, then one chip per
tracked job - its short label (the log's filename, extension stripped),
a colored status dot (blue = running, green = finished, red = crashed or
cancelled), and a **✕** to stop tracking it. Clicking a chip shows that
job's familiar Overview/Backbones/Sequences/Models/Results/Error
log/Output log tabs - exactly the single-job experience described above,
now with a banner above the tab strip naming which job you're looking at
(so it's never ambiguous once more than one is being tracked). Every
tracked job's status is polled independently of whichever one is
currently active (and stops once a job reaches a terminal state, same
reasoning as auto-refresh) so the chips and the overview list stay
current in the background.

**All jobs overview** is a card per tracked job - click one to switch to
that job's own tabs. Each card summarizes exactly what the single job's
own Overview tab shows: a progress box for whichever stage's counters are
most relevant right now ("Backbones done" while in/before RFdiffusion,
"Structures done" once modeling has started or finished, so a card never
goes blank just because a job has moved on to a later stage), an ETA for
that box (only filled in while the job is *actively* in that stage - a
job sitting in a later stage like scoring shows "—" rather than a stale
ETA), how many backbones passed filtering (blank/"pending" until
filtering has actually finished), and a "Total designs" count (accepted
backbones × sequences/backbone, same math as the single job Overview's
"Final designs" box). A job the dashboard can't resolve at all
(e.g. a log with no `output_dir:`/`PWD:` line) shows a short error message
in place of the metrics instead of a broken card.

**All jobs results** is the Results tab, unioned across every tracked
job: one combined list/table/alignment view/structure viewer, with a
synthetic **job** column identifying which job each row came from (jobs
with different metric columns - e.g. an AF3 run next to a Boltz run -
still combine cleanly, since the table columns are the union of every
job's own columns, blank where a given job doesn't have one). Filters,
"Color list by", and "Color structure by" all work the same way they do
on a single job's Results tab, just across the merged set. **Export
filtered (.zip)** downloads one archive covering every tracked job's
filtered rows - `filtered_output.csv` gains a `source_job` column, and
each job's model pdbs sit in their own `models/<job>/` subfolder so two
jobs' identically-named models (`model_0.pdb` is a very common name)
can't collide.

Under the hood, the per-job Results tab and All jobs results aren't two
separate implementations - `createResultsView()` in `app.js` is a
factory that both are instances of, parameterized by which DOM ids they
draw into and how to resolve a row back to its owning job's API calls
(`state.logPath` for the single-job instance, that row's own stamped-on
`__jobLogPath` for the multi-job one). Duplicating ~600 lines of
filtering/coloring/CSV/alignment/export logic across two nearly-identical
copies felt like exactly the kind of thing that quietly drifts out of
sync over time, so both tabs share one implementation instead.

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
  The viewer here has the same "Color residues by" RFdiffusion provenance
  option and interactive sequence panel as the Backbones and Results tabs
  (see above) — every `3_models/model_N` directory traces back to the
  same-indexed backbone's `.trb` file regardless of which sequence/
  sample/monomer-vs-complex variant a given row is, so provenance coloring
  is available here too. What the Models tab still doesn't have is a
  metrics-driven list — no filtering or "Color list by" a numeric column,
  since (unlike the Results tab) there's no single per-row CSV backing
  every model row here to filter or sort by; its own metrics table just
  shows that one selected model's own confidence values.
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
  `model_path` pdb files. A production job's result set can run into the
  tens of thousands of rows, so the model list, the metrics table, and
  each alignment-view sequence-length group are all **paginated** at 200
  rows (100 for an alignment group) per page rather than building DOM for
  every row at once - a small **‹ Prev / Next ›** footer appears once a
  set has more than one page. Filters and export still act on the *whole*
  filtered set regardless of which page you're looking at (a filtered
  export downloads every matching row, not just the current page); only
  what's actually built into the DOM is capped. Applying/clearing a
  filter or changing the metrics table's sort column jumps back to page 1
  of the affected view, since "page 3" means something different once the
  underlying row set or order has changed.

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
