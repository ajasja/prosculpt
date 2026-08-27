# Prosculpt Dashboard

A small local web dashboard for both **submitting** and **monitoring**
Prosculpt jobs. Two top-level tabs:

- **Run job** — build a config (guided form or hand-written YAML, or both
  kept in sync), optionally upload an input PDB / alignment files, and
  submit it - either as a local `slurm_runner.py` call (if the dashboard is
  running on the cluster itself) or over SSH to a cluster on the same
  network. See [§3](#3-run-a-new-job).
- **Track job** — everything this dashboard originally did: which stage a
  job is at, progress/ETA for RFdiffusion and modeling, plus browsable
  views of the backbones, designed sequences, individual models, and final
  results. A job submitted from the Run job tab is added here
  automatically once its log file appears on disk.

Tracking works by reading a job's log file and output directory directly
off disk, so **the dashboard needs filesystem access to both** — typically
because it's running on the same machine (e.g. a cluster login node) the
job runs on, or that filesystem is otherwise mounted where the dashboard
runs. Submitting a job additionally needs either local shell access to
`slurm_runner.py` (same-machine case) or a working `ssh`/`scp` client with
already-trusted key/agent access to the target cluster (remote case) - see
[dashboard_config.yaml.example](dashboard_config.yaml.example).

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

### Run automatically at Windows logon

`start_dashboard.ps1` (next to this README) wraps the same `python app.py`
call for exactly this: waits for a mapped network drive to actually
reconnect before launching (a persistent drive mapping isn't always ready
the instant a logon session starts), and logs everything to
`%USERPROFILE%\prosculpt_dashboard.log` since there's no console attached
once it's started this way. Open it and adjust the variables at the top
(`$DashboardDir`, `$Python`, `$Port`, `$DriveLetter`) for this machine
first, then either:

- **Task Scheduler** (recommended - can run hidden, restart on failure):
  Task Scheduler → Create Task → General: name it, "Run only when user is
  logged on" → Triggers: New → "At log on" → this user → Actions: New →
  "Start a program" → Program: `powershell.exe`, arguments:
  `-WindowStyle Hidden -ExecutionPolicy Bypass -File "C:\path\to\start_dashboard.ps1"`.
  On the **Settings** tab, uncheck (or raise) "Stop the task if it runs
  longer than 3 days" - that limit is normally meant for one-off jobs, and
  would otherwise silently kill the dashboard after 3 days of uptime.
- **Startup folder** (simpler, no admin needed): `Win+R` → `shell:startup`
  → create a shortcut with target
  `powershell.exe -WindowStyle Hidden -ExecutionPolicy Bypass -File "C:\path\to\start_dashboard.ps1"`.

Either way, since it runs with no visible window, check
`prosculpt_dashboard.log` if `http://localhost:5000` doesn't come up after
logging in - most often either the drive mapping timed out (raise
`$MaxWaitSecs`) or `$DashboardDir`/`$Python` needs adjusting for this
machine.

## 3. Run a new job

The **Run job** tab needs at least one *run target* configured before it
can submit anything - copy `dashboard/dashboard_config.yaml.example` to
`dashboard/dashboard_config.yaml` (git-ignored, host-specific - not
committed) and fill in at least one:

- **`kind: local`** - the dashboard itself runs `slurm_runner.py` as a
  local subprocess. Use this if the dashboard is running on the cluster
  (or a node with a working prosculpt installation) directly.
- **`kind: ssh`** - the dashboard shells out to your system's `ssh`/`scp`
  binaries to stage files and submit on a cluster it's on the same network
  as. Needs a working `ssh <ssh_host>` connection already trusted (key or
  agent - **no password is ever asked for, stored, or sent**; set up host
  key trust and key/agent auth once, outside the dashboard, the normal
  way). Optionally set `mounts` - a list of `{remote, local}` pairs, one
  per remote directory that's *also* reachable as a mounted filesystem
  path from wherever the dashboard runs - without at least one covering
  `projects_path`, submission still works, but a submitted job won't
  auto-appear on the Track job tab (see below), **and** any job's Track
  job listings (backbones/sequences/models/results) will come up empty
  even once added by hand, since there's nothing to translate that job's
  own remote output_dir against (see the next paragraph and
  Troubleshooting below).

Each target also needs `slurm_runner_path` (where the prosculpt
installation's `slurm_runner.py` lives - it needs its own working
`config/installation.yaml`, singularity images, etc. already set up there)
and `projects_path` (where new job directories get created). Also set
`python_path` to that installation's own conda env interpreter (the same
value as its `config/installation.yaml`'s `prosculpt_python_path`) -
without it, submission falls back to a bare `python`, which depends on
whatever's first on `$PATH` for the dashboard's own process (local) or a
non-interactive SSH session (remote) - usually **not** prosculpt's env,
and for SSH in particular `conda activate` often isn't even available in
that kind of session at all. Pointing straight at the env's own
interpreter avoids relying on shell activation entirely.

A `mounts` entry covering `projects_path` only helps track jobs staged
there - i.e. ones this dashboard's own Run job tab submitted. A job
launched directly against the installation instead (`python
slurm_runner.py some_config.yaml ++output_dir=Examples/...` run by hand
from `slurm_runner_path`, or via `run_tests.py` - both documented usages
of Prosculpt itself, independent of this dashboard) writes its output
under `slurm_runner_path` instead, which needs its own separate `mounts`
entry to be trackable the same way. Add as many entries as you actually
have mounted - none, one, or several, in any combination, including for
some other remote directory entirely if a job's data can land there too -
see `dashboard_config.yaml.example` for the full explanation and worked
examples.

The same file also has an optional top-level `defaults:` section for
dashboard-wide settings unrelated to any one target - currently just
`prediction_model` (which model the Run Job GUI's "Prediction model" field
preselects; falls back to Colabfold if unset). A home for anything else
dashboard-level that comes up later, without needing a second config file.

Two ways to submit:

- **Build a new job** - fill in the core settings and any of the
  independent, addable **modules** (Symmetry, RFdiffusion backbone
  filters, Inpaint sequence, Redesign-only, Partial diffusion, Boltz2
  templates - the last one only shown when the prediction model is
  Boltz2), optionally upload an input PDB (viewable, with a chain color
  legend, plus an interactive sequence panel - grouped into chunks with
  the chunk's starting residue number printed above it, for reading off
  residue numbers without hovering each one - contig/hotspot fields are
  always typed in by hand, the viewer doesn't feed them automatically) and,
  if "Use custom MSA" is checked, one or more `.a3m`
  alignment files (each filename must contain `Chain_<letter>`, matching
  the same convention `a3m_dir` already expects elsewhere in Prosculpt).
  The raw YAML underneath is kept in sync with the form both ways - edit
  either one; a manual YAML edit is validated (parseable, has a `contig`)
  before it's allowed to overwrite the form's values, and the form's own
  values always win back once you touch a form field again. On submit,
  the dashboard creates a directory named after the job (with the PDB, an
  `alignments/` subdirectory if any were uploaded, and the config yaml),
  and - for an `ssh` target - copies that whole directory over before
  submitting remotely.
- **Submit an existing project directory** - point at a directory you've
  already prepared by hand (your own config yaml plus anything it
  references, e.g. a custom filter script) - it's copied/submitted
  exactly as-is; nothing inside it is rewritten. Useful for a config more
  complex than the guided form covers, or one you're reusing from outside
  the dashboard entirely.

If a directory with the job's name already exists at the destination (the
`projects_path` directory locally, or that same path on the remote
cluster for an `ssh` target), it's never overwritten or merged into -
`_2`, `_3`, ... is appended until an unused name is found, and the result
banner says so plainly ("A directory named ... already existed there -
used ... instead"), including in a **Preview**. The config filename
tracks whatever name actually got used, so the directory and its own
config file never disagree about what the job is called.

Either way, **Preview (dry run)** runs `slurm_runner.py --dry-run` first
(prints the command it *would* run without touching `sbatch`) - worth
using before **Submit**, which asks for confirmation and actually
launches the job.

**squeue --me** (with a manual Refresh button, not auto-polled, so it
doesn't add load to the scheduler on every tab's auto-refresh timer) shows
your queued/running jobs for whichever target is selected.

Several fields (contig, "Also predict monomer", the redesign module's
designable-residues field, and the inpaint-sequence/partial-diffusion
range fields) have a small **?** next to their label - hover or focus it
for syntax examples and other guidance too long to leave always-visible
under the field itself.

**Pending submissions**: after a real (non-dry-run) submit, the dashboard
watches for that job's log file to appear (`logs/slurm-<job id>_*.out`
under the job's directory) and, once it does, adds it to the Track job tab
automatically - the exact same `addJobs()` a manually-pasted log path
would go through, so there's nothing job-submission-specific about how a
promoted job is tracked afterward. This only works for a `local` target,
or an `ssh` target with a `mounts` entry covering `projects_path` (see above); without
that, a pending entry says so and you add it to Track job by hand once you
know its log path.

## 4. Point it at a job

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

A **Show sidechains** control (four mutually-exclusive options: "Only
selected", "Interface", "Hydrophobic", "All") governs which residues'
sidechains are drawn as licorice sticks on top of the cartoon, independent
of whichever residue is currently click-highlighted. "Only selected" (the
default) shows nothing extra beyond whatever the click/sequence-panel
highlight already draws. "Interface" shows the sidechains of every residue
that has any atom within 5 Å of an atom belonging to a *different* chain -
computed per chain via `structure.getAtomSet()` intersected with
`getAtomSetWithinSelection("not :chain", 5)`, then expanded to whole
residues and OR'd across chains (`computeInterfaceSele()` in `app.js`) -
useful for eyeballing a binder's interface at a glance without having to
click through it residue by residue. "Hydrophobic" shows the sidechains of
every hydrophobic residue (ALA/VAL/LEU/ILE/PRO/PHE/MET/TRP), via NGL's own
built-in `hydrophobic` selection keyword - no separate residue list needed.
"All" shows every sidechain in the structure. This setting, like the
background and palette ones, is shared across every viewer and remembered
in local storage.

The Backbones, Models and Results tabs additionally have a **Color
residues by** / **Color structure by** dropdown, offering **RFdiffusion
provenance** as an alternative to chain coloring. RFdiffusion writes a
`.trb` sidecar file (a Python pickle, despite the `.trb` extension -
`load_trb_provenance()` in `parser.py` unpickles it, which needs `numpy`
since the pickle contains `numpy` array/scalar objects) next to each
backbone `.pdb`, recording which residues were copied in from the
reference structure versus generated de novo. The dashboard buckets every
residue into one of four categories, colored consistently across all
three tabs and shown in a small legend next to the dropdown, each with
its own color picker so the palette can be retuned: **Motif** (residues
present in a *redesigned* chain but taken from the reference, both
structure and sequence), **Fixed chains** (residues in a *non-designed*
chain, structure and sequence both taken from the reference verbatim),
**Inpainted sequence** (residues whose *structure* is fixed to the
reference but whose *sequence* ProteinMPNN is still free to redesign -
e.g. a fixed-chain residue with an unknown/masked identity), and
**Sculpted** (everything else - generated de novo, neither structure nor
sequence from the reference).

Earlier versions of this keyed provenance off one of RFdiffusion's own
`.trb` fields that pair a residue with its (chain, resnum) - every one
tried (`con_hal_pdb_idx`+`receptor_con_hal_pdb_idx`, then
`con_hal_pdb_idx`+`complex_con_hal_pdb_idx`) turned out to get chain
letters and/or residue numbers wrong in some case. `load_trb_provenance()`
now reads three fields keyed purely by a residue's flat, 0-based
*position* in the generated structure instead (residue 0 is the first
residue of the first chain, residue 1 the next, and so on across every
chain in file order - no chain-letter/resnum bookkeeping at all):
`con_hal_idx0` (flat indices of motif residues), `inpaint_str` (per-residue
bool - True wherever the structure is fixed to the reference), and
`inpaint_seq` (per-residue bool - True wherever the sequence is fixed to
the reference). `inpaint_str` False is Sculpted regardless of anything
else; of the rest (`inpaint_str` True), `inpaint_seq` False is Inpainted
sequence - this takes priority over `con_hal_idx0` membership, since a
motif residue can still have its sequence freed up for redesign and the
point of this coloring is to show what's actually still changeable; of
what's left, `con_hal_idx0` membership is Motif, and everything else is
Fixed chains.
Because this is purely positional, it needs no reconciling of
RFdiffusion's own raw backbone `.pdb` (Backbones tab, which keeps counting
resnums up across chain boundaries) against AlphaFold3/Boltz's output
(Models and Results tabs, which resets every chain back to 1) - residue
order (and count) is preserved end-to-end through the pipeline regardless
of that numbering difference, or of `rechain_rfdiff_pdbs()` reassigning
chain letters after RFdiffusion runs. The frontend matches a loaded
structure's residues up against this flat category list via NGL's own
`atom.residueIndex` (`makeProvenanceColorScheme()` in `app.js`) - the
same flat, 0-based, file-order position, so no chain/resnum lookup is
needed there either.

Every `3_models/model_N` directory traces back to the same-indexed
`1_rfdiff/_N.trb` regardless of which sequence/sample/monomer-vs-complex
variant a given row is (`list_models()` in `parser.py` resolves this and
stamps each row with its `trb_path`), and a Results row finds its `.trb`
the same way the Backbones tab does - same basename as that row's
`path_rfdiff`, `.trb` extension. If a backbone/model has no `.trb` file,
or it can't be read (e.g. `numpy` isn't installed in the dashboard's own
environment), the dropdown falls back to chain coloring with a short
explanation in place of the legend rather than failing silently.

## 5. Track multiple jobs at once

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

Also on the per-job tab bar, pushed to the far right and styled with a red
outline: a **Cancel job** button, enabled only while the job still looks to
be running (disabled once it's finished, crashed, or already been
cancelled). It runs `scancel` against the job's own Slurm job ID (the same
one shown on the Overview tab, parsed from the log's "Hello from job ..."
line). Since Track job has no other way to know which configured target a
tracked job's cluster actually is - a job just added by pasting a log path
was never necessarily submitted through the Run job tab - clicking it opens
a small dialog asking which target owns the job, and requires an explicit
"I'm sure I want to cancel this job" confirmation before anything is
actually sent.

## Design notes / known limitations

- **File-serving paths are checked by *real* location, not by string.**
  Every endpoint that serves a file under `output_dir` (`.trb`, model pdbs,
  confidence json, ...) re-validates the requested path actually resolves
  inside `output_dir` before reading it (`_is_within()` in `app.py`), using
  `os.path.realpath()` rather than a plain string-prefix check. This
  matters on HPC storage in particular: it's common for a project
  directory to be reachable under more than one absolute path (e.g. a
  symlink or bind-mount from a user's home folder into shared/Ceph
  storage, or an older archived-elsewhere `1_rfdiff/`), and a lexical
  check would reject a perfectly legitimate file the moment it's expressed
  via a different alias than the one `output_dir` happened to resolve to -
  `realpath()` resolves symlinks first, so both aliases compare equal.
- **`final_output.csv` is returned exactly as written; paths inside it are
  translated only at the point something actually opens a file, not
  baked into the data itself.** `output_dir` itself is translated in
  `resolve_output_dir()` (see `mounts`
  in §3) - everything the Backbones/Sequences/Models tabs and the Results
  table list is built from that already-local `output_dir` server-side, so
  none of it needs any further translation. `final_output.csv` is
  different: `load_final_csv()` in `parser.py` returns it completely
  as-is, including any path-like column (`model_path`, `path_rfdiff`, and
  whichever others a given prediction backend adds), because Prosculpt
  itself records those as paths *on the cluster the job ran on*, and nothing
  reads that file to *display* rewritten data - the raw values are what
  populate the results table and get written into "Export filtered"'s own
  `filtered_output.csv`. The two places that actually need to *open* one of
  those files - `path_rfdiff`, fetched for the Results tab's "Color by
  RFdiffusion provenance" feature, and `model_path`, looked up to include
  the right pdb bytes in an export zip - translate it themselves, right
  where they use it (`RT.translate_remote_path()`/`RT.make_path_translator()`,
  called from `api_trb()`/`api_export_filtered()`/`api_export_filtered_multi()`
  in `app.py`), rather than the read itself rewriting anything. Anything
  that reads a path out of job data going forward and needs to actually
  open the file it points to (a future confidence-JSON field, a future
  `.trb` field, ...) should do the same: translate it right there, not in
  whatever function first reads the surrounding data.
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
app.py              Flask routes / API for the Track job side
run_api.py          Flask routes / API for job submission/monitoring
                    (submit/squeue/cancel/check_pending) - registered as a
                    blueprint from app.py; squeue/cancel are also used by
                    Track job's own "Cancel job" button, not just Run job
run_targets.py       Run-target config + dashboard-wide defaults (both read
                    from dashboard_config.yaml) + the shared local/ssh
                    execution primitive
job_staging.py       Assembles a new job's directory on disk (pdb,
                    alignments/, logs/, config yaml) for the "build a new
                    job" flow
dashboard_config.yaml.example   Copy to dashboard_config.yaml and fill in -
                    git-ignored
parser.py            All log-parsing + filesystem-scanning logic (no Flask
                    dependency, easy to test standalone)
templates/index.html      Track job markup (Run job's own markup lives in
                    run_job.html, included from here)
templates/run_job.html
static/css/style.css
static/js/app.js          Track job frontend (vanilla JS, no build step;
                    uses NGL Viewer from a CDN for structure viewing)
static/js/run_job_schema.js  Run job's field/module schema (CORE_FIELDS,
                    MODULES, CONTIG_HELP) - no logic, just data
static/js/run_job.js         Run job's form rendering, YAML two-way sync,
                    upload handling, submit/squeue/pending-runs logic
requirements.txt
start_dashboard.ps1  Windows: launch at logon via Task Scheduler/Startup -
                    see "Run automatically at Windows logon" above
```

## Troubleshooting

- **"Could not find PWD: and/or output_dir: in the log file"** — the
  dashboard expects the exact log format Prosculpt writes (a `PWD:` line
  followed by the working directory on the next line, and a bare
  `output_dir: <path>` line from the printed config). If your log format
  differs, check `resolve_output_dir()` in `parser.py`.
- **A tab shows "No X generated yet" for something you know exists** —
  two different causes, depending on whether this is a remote (`ssh`)
  target:
  - The dashboard looks for very specific file patterns (e.g. `_*.pdb` in
    `1_rfdiff/`, `*summary_confidences.json` / `confidence_*.json` under
    `3_models/`). If Prosculpt's output layout has changed, these glob
    patterns in `parser.py` are the place to update.
  - For an `ssh` target: this is the expected result of an *unmapped*
    remote path, not a bug - it means the dashboard resolved the job's
    output_dir correctly (as a path on the remote cluster's own
    filesystem) but has no `mounts` entry to translate it into a path
    actually reachable from here, so every listing comes up empty while
    the log itself still loads fine (the log's own path needs no
    translation - only the paths *printed inside* it do). Check that this
    job's output_dir actually falls under one of this target's `mounts`
    (each entry's `remote` needs to be an exact prefix of it) - typically
    one entry covering `projects_path` (jobs the Run job tab itself
    staged) and another covering `slurm_runner_path` (jobs launched
    directly against the installation, e.g. via `run_tests.py` or a
    hand-run `slurm_runner.py` call), though a target can have as many
    mounts as it actually needs - see the `mounts` notes in
    `dashboard_config.yaml.example` and §3 above. A job whose output_dir
    falls under none of them (e.g. one that predates any target being
    configured, or one run somewhere unrelated to all of them) can't be
    translated no matter what's configured - `translate_remote_path()` in
    `run_targets.py` is the place to look if you need something more
    involved than a list of mount prefixes already covers.
- **Error log tab doesn't appear** — it's only shown once a `.err` file is
  actually found next to the log (same basename, `.err` extension); the
  Overview card's "Error log" stat says "not found" in that case. A job
  that's merely slow or stalled, with no crash/cancellation, will still
  show the tab (if the `.err` file exists) but without a warning banner.
- **The whole dashboard freezes for every user at once, then resumes the
  instant someone presses a key in its console window** — this is Windows
  console "QuickEdit Mode": clicking, scrolling, or selecting text in a
  `cmd.exe`/`powershell.exe` window pauses that console entirely, so any
  process (this one included) blocks the moment it tries to write another
  line to it - including Werkzeug's own per-request access log line (the
  `"GET ... 200 -"` lines) - until a key is pressed to cancel the
  selection. Two things compound this into "the whole dashboard is dead
  for everyone", not just a paused terminal:
  - `app.run()` was single-threaded by default, so the one worker thread
    stuck on that blocked console write couldn't pick up any other
    request either - fixed by passing `threaded=True` (already done in
    `app.py`), so one stuck request (a frozen console write, or anything
    else slow - a stalled network mount under a glob(), a huge
    `final_output.csv`, ...) only ever stalls itself, not every other
    user's request.
  - The freeze itself still happens (threading doesn't stop the console
    from pausing, just stops it from taking the whole server down with
    it) - to stop it happening at all, either disable QuickEdit Mode on
    that console window (right-click the title bar → Properties →
    Options → untick "QuickEdit Mode"), or don't run the dashboard in a
    visible, clickable console window in the first place - launch it via
    `start_dashboard.ps1` instead (see "Run automatically at Windows
    logon" above), which redirects all output to a log file and runs with
    no console window to accidentally click into.
