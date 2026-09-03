# Prosculpt Dashboard

A local web dashboard for submitting and monitoring Prosculpt jobs.

- **Run job** — build a config (guided form or YAML, kept in sync) and submit it, locally or over SSH.
- **Track job** — monitor one or more jobs: stage/progress/ETA, backbones, sequences, models, and results, read straight off disk.

A job is tracked by pointing at its own **output directory** (the folder with `logs/` and one numbered `01/`, `02/`, ... subdirectory per SLURM array task). The dashboard needs filesystem access to that directory - typically because it's running on the same machine the job runs on (e.g. a cluster login node).

## Install

```bash
cd dashboard
pip install -r requirements.txt
```

Needs Python 3.9+.

## Run

```bash
python app.py
```

Opens on **http://localhost:5000**. Different port: `PORT=8080 python app.py`, or set `defaults: {port: 8080}` in `dashboard_config.yaml` (see below) to make it permanent - the env var wins if both are set.

It listens on all interfaces, so running it on a cluster login node makes it reachable at `http://<that machine>:5000` for anyone on the same network. If you'd rather it not be, tunnel in instead:

```bash
ssh -L 5000:localhost:5000 you@cluster
```

and open http://localhost:5000 on your laptop as usual.

Optional password gate: set `auth: {enabled: true, password: "..."}` in `dashboard_config.yaml` (see below) to require a password before anything else loads. Worth turning on if you're sharing it over a network - it's a plain-text password meant to keep casual traffic out, not to protect anything sensitive.

### Keep it running

- **Windows** (logon via Task Scheduler or the Startup folder): `start_dashboard.ps1` - adjust the variables at its own top, then see its header comment for the exact setup steps. Logs to `%USERPROFILE%\prosculpt_dashboard.log`.
- **Linux** (under `screen`): `start_dashboard.sh` - adjust the variables at its own top, then:
  ```bash
  screen -dmS prosculpt_dashboard bash /path/to/start_dashboard.sh   # start, detached
  screen -r prosculpt_dashboard                                     # reattach (Ctrl-A D to detach again)
  screen -X -S prosculpt_dashboard quit                             # stop
  ```
  Logs to `~/prosculpt_dashboard.log`, and stays visible live when reattached. If something on the machine keeps killing the process, run it through `watchdog.sh` instead (same commands, just swap in `watchdog.sh` for `start_dashboard.sh`) - it relaunches `start_dashboard.sh` automatically whenever it exits, logging each restart to the same file.

## Configure run targets (for Run job)

Copy `dashboard_config.yaml.example` to `dashboard_config.yaml` (git-ignored, machine-specific) and fill in at least one target:

- **`kind: local`** - the dashboard runs `slurm_runner.py` itself. Use when it's running on the cluster (or a node with a working prosculpt installation) directly.
- **`kind: ssh`** - the dashboard shells out to your system's `ssh`/`scp` to a cluster on the same network. Needs an already-trusted connection (key/agent - no password is ever asked for, stored, or sent).

Each target needs `slurm_runner_path` (the prosculpt installation), `projects_path` (where new job directories go), and `python_path` (that installation's own conda env interpreter - avoids depending on shell activation). An `ssh` target's optional `mounts` let a submitted job get auto-tracked and its results actually browsable - see the example file for the full explanation and worked examples.

## Run job

Two ways to submit:

- **Build a new job** - core settings plus optional, independently-addable modules (Symmetry, RFdiffusion backbone filters, Inpaint sequence, Redesign-only, Partial diffusion, Boltz2 templates, Filtering & post-filtering scoring), an optional PDB/alignment upload, and a pre-filtering scoring script override. The raw YAML underneath stays in sync with the form both ways.
- **Submit an existing project directory** - point at a directory you've already prepared by hand; submitted as-is (copied over first for an `ssh` target), nothing inside it rewritten. One project directory can be resubmitted more than once - if its `output_dir` is still there from a previous run, a fresh `_2`, `_3`, ... is picked automatically.

**Submit** asks for confirmation and launches the job; a successful submission is added to Track job automatically. **squeue --me** shows your queue for the selected target.

## Track job

Paste a job's output directory path (or **Browse…** to pick one, or several, from the filesystem) and click **Add job**. Track several jobs at once and, on top of each job's own tabs, you get three combined views across all of them: **All jobs overview**, **All jobs results**, and **All jobs filtered results**.

Per-job tabs:

- **Overview** - current stage, progress/ETA, run plan, step timings, crash/cancellation banners, and the post-filtering scoring job's own state if the job configured one.
- **Backbones / Sequences / Models** - RFdiffusion backbones, MPNN sequences (shared alignment viewer), and predicted structures, each with an interactive 3D viewer (click/hover residue picking, chain or RFdiffusion-provenance coloring, palettes).
- **All Results / Filtered Results** - every model, or the subset that passed prosculpt's own `filtering:` stage - a metrics table, numeric filters, sort/color-by-metric, and zip export.
- **Error log / Output log** - the job's `.err`/`.out`, if found.
- **Cancel job** - sends `scancel` for the job's own Slurm ID, behind a confirmation dialog.

Auto-refreshes every 5s (toggle off next to the job field); stops once a job reaches a terminal state (finished/crashed/cancelled).

## Project layout

```
app.py                Flask routes / API for the Track job side
run_api.py            Flask routes / API for job submission/monitoring -
                      registered as a blueprint from app.py
run_targets.py        Run-target config + dashboard-wide defaults + the
                      shared local/ssh execution primitive
job_staging.py        Assembles a job's directory on disk for the
                      "build a new job" flow
parser.py             Log-parsing + filesystem-scanning logic (no Flask
                      dependency, easy to test standalone)
dashboard_config.yaml.example   Copy to dashboard_config.yaml and fill in -
                      git-ignored
templates/index.html  Track job markup (Run job's own markup lives in
                      run_job.html, included from here)
templates/run_job.html
templates/login.html  The password screen, when auth is enabled
static/css/style.css
static/js/app.js      Track job frontend (vanilla JS, no build step)
static/js/run_job_schema.js  Run job's field/module schema - no logic, just data
static/js/run_job.js  Run job's form rendering, YAML two-way sync,
                      upload handling, submit/squeue logic
requirements.txt
start_dashboard.ps1   Windows: launch at logon via Task Scheduler/Startup
start_dashboard.sh    Linux: launch under `screen` on a headnode
watchdog.sh           Linux: wraps start_dashboard.sh in a restart loop
```

## Troubleshooting

- **"No output has been generated yet"** - the directory looks right (it has a `logs/`) but no `01/`, `02/`, ... task directory exists yet. Normal for a job that's queued or still in its first step.
- **"The tracked directory doesn't look like an output directory"** - no numbered task directories *and* no `logs/`. Usually a directory level off: pointed at the project directory containing the output directory, or at a single task's own `01/` instead of its parent.
- **A tab shows "No X generated yet" for something that exists** - the dashboard looks for specific file patterns per task (see `parser.py`); if Prosculpt's output layout changed, that's where to update it.
- **Dashboard frozen for everyone, resumes the instant someone presses a key in its console window** - Windows console "QuickEdit Mode" pauses the whole console (and so the process) on click/select. Disable it (console title bar → Properties → Options → untick QuickEdit Mode), or run via `start_dashboard.ps1` instead, which has no console window to click into.
