// Prosculpt Dashboard frontend
// Vanilla JS, no build step. Talks to the Flask API in app.py.

const STORAGE_KEY = "prosculpt_dashboard_log_path";
const THEME_KEY = "prosculpt_dashboard_theme";
const APP_MODE_KEY = "prosculpt_dashboard_app_mode";
const POLL_MS = 5000;

const STAGE_ORDER = ["rfdiffusion", "filtering", "mpnn", "modeling", "scoring", "finished"];
const STAGE_LABEL = {
  setup: "Starting up",
  rfdiffusion: "RFdiffusion",
  filtering: "Filtering",
  mpnn: "ProteinMPNN",
  modeling: "Modeling",
  scoring: "Scoring / final ops",
  finished: "Finished",
};

const AA_COLORS = {
  A: "#8fd19e", V: "#8fd19e", L: "#8fd19e", I: "#8fd19e", M: "#8fd19e",
  F: "#f4a261", W: "#f4a261", Y: "#f4a261",
  K: "#6ea8fe", R: "#6ea8fe", H: "#8ecae6",
  D: "#ef5a6f", E: "#ef5a6f",
  N: "#c77dff", Q: "#c77dff",
  S: "#ffd166", T: "#ffd166",
  C: "#ffb703", G: "#adb5bd", P: "#e5989b",
};

// Columns in final_output.csv that hold filesystem paths - cumbersome in the
// middle of the table, so they get pushed to the end (see D1 in the spec).
const PATH_COLUMNS = ["model_path", "af3_json", "af3_pdb", "path_rfdiff"];

// state.logPath is "the currently active job" - the one the per-job tabs
// (Overview/Backbones/.../Output log) show. Multi-job tracking adds a
// layer above that: state.jobs is every job being watched, and
// state.topView picks between showing the active job's tabs ("job") or
// one of the two aggregate tabs ("all-overview"/"all-results") - see
// switchTopView()/switchToJob() further down. The Results-tab-shaped
// state that used to live directly on `state` (results/csvSort/
// loadedResultPdbName/resultsColorMode) now lives inside createResultsView()'s
// own closure instead, once per view instance, so the per-job Results tab
// and the All jobs results tab don't share (or fight over) one copy of it.
let state = {
  // "run" vs "track" is a higher-level split than everything else below -
  // Track job is the entire existing tracking UI (topView/pill bar/per-job
  // tabs), unchanged; Run job is a wholly separate sibling. See
  // showAppMode()/initAppModeTabs().
  appMode: localStorage.getItem(APP_MODE_KEY) || "track",
  logPath: localStorage.getItem(STORAGE_KEY) || "",
  activeTab: "overview",
  lastStatus: null,
  pollTimer: null,
  browsePath: null,

  jobs: [], // ordered list of tracked log paths
  jobStatuses: {}, // logPath -> { status, terminal }
  topView: "job", // "job" | "all-overview" | "all-results"

  backbonesCache: [],
  selectedBackboneKey: null,
  loadedBackboneKey: null,
  showFilteredBackbones: true,
  backboneColorMode: "chain", // "chain" | "provenance"

  sequencesCache: null,

  modelsCache: [],
  selectedModelKey: null,
  loadedModelStructureKey: null,
  loadedModelConfKey: null,
  modelsColorMode: "chain", // "chain" | "provenance"
};

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

function qs(sel, root = document) { return root.querySelector(sel); }
function qsa(sel, root = document) { return Array.from(root.querySelectorAll(sel)); }

// params.log lets a caller target a specific job's API explicitly (the
// All jobs tabs, where a row/job isn't necessarily the active one) -
// otherwise this defaults to whichever job is currently active.
async function apiGet(path, params = {}) {
  const url = new URL(path, window.location.origin);
  params.log = params.log || state.logPath;
  Object.entries(params).forEach(([k, v]) => { if (v !== undefined && v !== null) url.searchParams.set(k, v); });
  const res = await fetch(url);
  if (!res.ok) {
    let msg = res.statusText;
    try { const j = await res.json(); msg = j.error || j.message || msg; } catch (e) {}
    const err = new Error(msg);
    // Lets a caller tell "this exact request is malformed/invalid and will
    // never succeed no matter how many times it's retried" (4xx) apart
    // from "probably transient, worth trying again" (network failure, a
    // 5xx) - see pollPendingRuns() in run_job.js for the caller that
    // actually needs this distinction.
    err.status = res.status;
    throw err;
  }
  return res.json();
}

// Mirrors apiGet()'s own error handling, for the handful of call sites
// that fetch() a PDB/text file directly (structure viewers) instead of
// going through apiGet() (JSON API calls) - a failure there gets rendered
// straight into the viewer's own "could not load" message rather than
// thrown up to a caller that already knows how to display an Error.
// Every error path in this app (including the OSError handler in app.py)
// returns {"error": "..."} on failure, so that's tried first; falls back
// to the response's statusText if the body isn't parseable JSON at all
// (e.g. a proxy/webserver-level failure that never reached Flask).
async function fetchErrorMessage(res) {
  try {
    const j = await res.json();
    return j.error || j.message || res.statusText || `HTTP ${res.status}`;
  } catch (e) {
    return res.statusText || `HTTP ${res.status}`;
  }
}

function fmtSeconds(s) {
  if (s === null || s === undefined || isNaN(s)) return "—";
  s = Math.round(s);
  if (s < 60) return `${s}s`;
  const h = Math.floor(s / 3600);
  const m = Math.floor((s % 3600) / 60);
  const sec = s % 60;
  if (h > 0) return `${h}h ${m}m`;
  return `${m}m ${sec}s`;
}

function escapeHtml(str) {
  return String(str).replace(/[&<>"']/g, (c) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
  }[c]));
}

function basename(p) {
  return String(p || "").split(/[\\/]/).pop();
}

// Rebuild a scrollable container's contents while keeping its scroll
// position - a plain innerHTML replace resets scrollTop to 0, which is
// disruptive on every 5s poll if the user is mid-scroll through a list.
function rerenderPreservingScroll(el, renderFn) {
  const top = el.scrollTop;
  renderFn();
  el.scrollTop = top;
}

// ---------------------------------------------------------------------
// Theme
// ---------------------------------------------------------------------

function applyTheme(theme) {
  document.documentElement.dataset.theme = theme;
  qs("#themeToggleBtn").textContent = theme === "light" ? "🌙" : "☀️";
  qs("#themeToggleBtn").title = theme === "light" ? "Switch to dark theme" : "Switch to light theme";
}

function initTheme() {
  const stored = localStorage.getItem(THEME_KEY);
  const systemPrefersLight = window.matchMedia && window.matchMedia("(prefers-color-scheme: light)").matches;
  applyTheme(stored || (systemPrefersLight ? "light" : "dark"));
  qs("#themeToggleBtn").addEventListener("click", () => {
    const next = document.documentElement.dataset.theme === "light" ? "dark" : "light";
    applyTheme(next);
    localStorage.setItem(THEME_KEY, next);
  });
}

// ---------------------------------------------------------------------
// Top bar / log loading
// ---------------------------------------------------------------------

function initTopbar() {
  const input = qs("#logPathInput");

  qs("#loadBtn").addEventListener("click", () => {
    addJobs(input.value);
    input.value = "";
  });
  // Plain Enter adds the job(s) (matching the old single-line input's
  // behavior); Shift+Enter inserts an actual newline instead, for typing
  // a second path by hand rather than pasting a multi-line list.
  input.addEventListener("keydown", (e) => {
    if (e.key === "Enter" && !e.shiftKey) {
      e.preventDefault();
      qs("#loadBtn").click();
    }
  });

  qs("#autoRefreshToggle").addEventListener("change", (e) => {
    qs("#autoRefreshNote").classList.add("hidden");
    if (e.target.checked) startPolling(); else stopPolling();
  });

  // No fallback path on the very first open (state.browsePath is still
  // null then) - openBrowse() omits ?path= entirely in that case, letting
  // the server pick a starting directory (the configured projects root,
  // if any - see get_default_browse_root() in run_targets.py) instead of
  // this always defaulting to "." (the dashboard process's own cwd,
  // rarely where anyone's actual jobs are).
  qs("#browseBtn").addEventListener("click", () => openBrowse(state.browsePath));
  qs("#browseCloseBtn").addEventListener("click", () => qs("#browseModal").classList.add("hidden"));
  qs("#browseAddSelectedBtn").addEventListener("click", () => {
    if (!browseSelectedPaths.size) return;
    addJobs(Array.from(browseSelectedPaths).join("\n"));
    qs("#browseModal").classList.add("hidden");
  });
  qs("#browseFilterInput").addEventListener("input", (e) => renderBrowseEntries(e.target.value));

  loadBrowseRoots();
  qs("#browseRootSelect").addEventListener("change", (e) => {
    const path = e.target.value;
    e.target.value = ""; // a jump, not a persistent selection - avoid it going stale as the user navigates further
    if (path) openBrowse(path);
  });
  qs("#browseGoToPathBtn").addEventListener("click", () => goToBrowsePath());
  qs("#browseGoToPathInput").addEventListener("keydown", (e) => {
    if (e.key === "Enter") { e.preventDefault(); goToBrowsePath(); }
  });
}

function goToBrowsePath() {
  const input = qs("#browseGoToPathInput");
  const path = input.value.trim();
  if (!path) return;
  input.value = "";
  openBrowse(path);
}

// One dropdown entry per configured target whose root is actually reachable
// from this machine right now (see list_browse_roots() in run_targets.py) -
// lets a multi-cluster setup jump straight to another cluster's project
// root instead of pasting/typing its path by hand. Silently leaves the
// dropdown at just its placeholder option if none are configured/reachable
// (e.g. no dashboard_config.yaml yet) - the existing folder-by-folder
// browsing and the "paste a path" input above still work either way, so
// this is a nice-to-have, never a hard requirement to open the browser at
// all.
async function loadBrowseRoots() {
  const sel = qs("#browseRootSelect");
  try {
    const roots = await fetch("/api/browse_roots").then((r) => r.json());
    if (!Array.isArray(roots) || !roots.length) return;
    roots.forEach((r) => {
      const opt = document.createElement("option");
      opt.value = r.path;
      opt.textContent = r.label;
      opt.title = r.path;
      sel.appendChild(opt);
    });
  } catch (e) {
    // Leave the dropdown at just its placeholder - not worth surfacing an
    // error for a nice-to-have shortcut when the rest of the browser still
    // works fine without it.
  }
}

// Loading a (possibly new, possibly still-running) job should always start
// from a clean slate: auto-refresh resumes even if the previous job we
// were watching had finished/crashed/been cancelled and auto-stopped
// polling, and every "don't refetch/redraw a viewer unless the underlying
// file actually changed" guard gets cleared - those guards compare
// filenames (backbone name, model path, result pdb name) which are reused
// across unrelated jobs all the time ("_0.pdb", "model_0", ...), so
// without this a genuinely different job's structures could silently fail
// to replace what an earlier job had already drawn into the viewers.
function prepareForNewJob() {
  qs("#autoRefreshToggle").checked = true;
  qs("#autoRefreshNote").classList.add("hidden");
  startPolling();
  resetPerJobViewState();
}

function resetPerJobViewState() {
  state.backbonesCache = [];
  state.selectedBackboneKey = null;
  state.loadedBackboneKey = null;

  state.sequencesCache = null;

  state.modelsCache = [];
  state.selectedModelKey = null;
  state.loadedModelStructureKey = null;
  state.loadedModelConfKey = null;

  resultsView.reset();

  // Clear anything already drawn immediately, rather than leaving the
  // previous job's structure on screen until the next tab visit's fetch
  // happens to overwrite it.
  ["#backboneViewer", "#modelViewer", "#resultViewer"].forEach((sel) => {
    const el = qs(sel);
    if (el) el.innerHTML = "";
  });
  const backboneLabel = qs("#backboneViewerLabel"); if (backboneLabel) backboneLabel.textContent = "";
  const modelLabel = qs("#modelViewerLabel"); if (modelLabel) modelLabel.textContent = "";
  const resultsSeq = qs("#resultsSequenceContent"); if (resultsSeq) resultsSeq.innerHTML = "";
}

// The current directory's raw listing (from the last successful
// /api/browse call) - kept around so the filter box (renderBrowseEntries)
// can re-render from it instantly on every keystroke without refetching.
let browseListing = { parent: null, entries: [] };
// Which file paths are checked, tracked independently of the DOM (not
// read back from ".browse-entry input:checked") - the filter box hides
// non-matching entries by not rendering them at all, which would silently
// forget a checked box's state the moment its row scrolled out of the
// (filtered) DOM otherwise. Reset per directory, same as before the
// filter box existed - navigating to a different directory already
// replaced the whole listing (and so, implicitly, any selection) even
// before this.
let browseSelectedPaths = new Set();

// File rows get a checkbox instead of loading immediately on click, so
// several jobs can be picked in one visit before committing with "Add
// selected jobs" (directories still navigate on click, same as before).
async function openBrowse(path) {
  qs("#browseModal").classList.remove("hidden");
  qs("#browseFilterInput").value = "";
  browseSelectedPaths = new Set();
  updateBrowseSelectedCount();
  try {
    // No `path` at all (the very first open, before state.browsePath is
    // ever set) omits ?path= entirely instead of sending an empty string -
    // the server then picks its own starting directory (see
    // get_default_browse_root() in run_targets.py) rather than this
    // needing to know or guess one itself.
    const url = path ? `/api/browse?path=${encodeURIComponent(path)}` : "/api/browse";
    const data = await fetch(url).then((r) => r.json());
    state.browsePath = data.path;
    qs("#browsePath").textContent = data.path;
    browseListing = { parent: data.parent, entries: data.entries };
    if (data.preselect) browseSelectedPaths.add(data.preselect);
    renderBrowseEntries("");
    updateBrowseSelectedCount();
  } catch (err) {
    browseListing = { parent: null, entries: [] };
    qs("#browseEntries").innerHTML = `<p class="muted">${escapeHtml(err.message)}</p>`;
  }
}

// Purely a client-side re-render of the already-fetched browseListing -
// case-insensitive substring match on the entry's own name. The parent
// ("..") row is never filtered out (it's navigation, not a search result)
// - only actual entries. Checkbox state comes from browseSelectedPaths,
// not a fresh `false` default, so re-filtering (which rebuilds these rows
// from scratch) never drops an existing selection.
function renderBrowseEntries(filterText) {
  const wrap = qs("#browseEntries");
  wrap.innerHTML = "";
  const { parent, entries } = browseListing;
  if (parent) {
    const up = document.createElement("div");
    up.className = "browse-entry";
    up.innerHTML = "⬆️ ..";
    up.addEventListener("click", () => openBrowse(parent));
    wrap.appendChild(up);
  }
  const needle = filterText.trim().toLowerCase();
  const filtered = needle ? entries.filter((e) => e.name.toLowerCase().includes(needle)) : entries;
  if (needle && !filtered.length) {
    wrap.insertAdjacentHTML("beforeend", `<p class="muted">No entries match "${escapeHtml(filterText)}".</p>`);
  }
  filtered.forEach((e) => {
    const div = document.createElement("div");
    div.className = "browse-entry";
    if (e.is_dir) {
      div.innerHTML = `📁 ${escapeHtml(e.name)}`;
      div.addEventListener("click", () => openBrowse(e.path));
    } else {
      const checked = browseSelectedPaths.has(e.path) ? "checked" : "";
      div.innerHTML = `<input type="checkbox" data-path="${escapeHtml(e.path)}" ${checked}> 📄 ${escapeHtml(e.name)}`;
      const cb = div.querySelector("input");
      cb.addEventListener("change", () => {
        if (cb.checked) browseSelectedPaths.add(e.path); else browseSelectedPaths.delete(e.path);
        updateBrowseSelectedCount();
      });
      div.addEventListener("click", (ev) => { if (ev.target.tagName !== "INPUT") cb.click(); });
    }
    wrap.appendChild(div);
  });
}

function updateBrowseSelectedCount() {
  const n = browseSelectedPaths.size;
  qs("#browseSelectedCount").textContent = n ? `${n} selected` : "";
  qs("#browseAddSelectedBtn").disabled = n === 0;
}

// ---------------------------------------------------------------------
// Multi-job tracking: the top-level selector (two "All jobs ..." tabs
// plus one chip per tracked job) above the per-job tab bar.
// ---------------------------------------------------------------------

const JOBS_STORAGE_KEY = "prosculpt_dashboard_jobs";

function jobLabel(logPath) {
  const base = basename(logPath);
  return base.replace(/\.[^.]+$/, "") || logPath;
}

function persistJobs() {
  localStorage.setItem(JOBS_STORAGE_KEY, JSON.stringify(state.jobs));
}

function loadJobsFromStorage() {
  try {
    const stored = JSON.parse(localStorage.getItem(JOBS_STORAGE_KEY) || "[]");
    if (Array.isArray(stored)) state.jobs = stored.filter((p) => typeof p === "string" && p);
  } catch (e) {}
  // Backward compatibility: a pre-multi-job session only ever persisted a
  // single active log path (STORAGE_KEY) - fold that in as a tracked job
  // too so upgrading doesn't silently drop what was already loaded.
  if (state.logPath && !state.jobs.includes(state.logPath)) state.jobs.unshift(state.logPath);
}

// Splits on newlines or commas so pasting several paths at once (or the
// Browse modal's multi-select "Add selected jobs") works the same way a
// single "paste one path, click Add job" does.
//
// switchToLast defaults to true (a human explicitly adding job(s) clearly
// wants to see the last one) but is passed false by Run Job's background
// auto-promotion (see pollPendingRuns() in run_job.js) - that call site
// can fire repeatedly in quick succession as a multi-task array job's
// tasks are discovered one at a time, and switchToJob() is not idempotent
// enough to call back-to-back like that: it resets per-job view state
// (backbones/sequences/models caches, clears those viewers' DOM) and
// kicks off a fresh refreshAll() fetch every time it's called, so several
// calls in a row race each other - a later switch's reset can land while
// an earlier switch's fetch is still in flight, leaving whichever job
// ends up "active" with its backbones/sequences view corrupted (empty)
// even though its own status/progress - a plain single-value overwrite,
// not staged cache data - still resolves fine. A background add has no
// business yanking the user's current view around at all, on top of that.
function addJobs(rawInput, switchToLast = true) {
  const paths = rawInput.split(/[\n,]+/).map((p) => p.trim()).filter(Boolean);
  if (!paths.length) return;
  let last = null;
  paths.forEach((p) => {
    if (!state.jobs.includes(p)) {
      state.jobs.push(p);
      fetchJobStatus(p); // don't wait for the next poll tick to know its stage
    }
    last = p;
  });
  persistJobs();
  renderJobTabsBar();
  if (last && switchToLast) switchToJob(last);
  else if (state.topView === "all-overview") renderAllOverview(); // reflect the new card even without switching to it
}

function removeJob(logPath) {
  const idx = state.jobs.indexOf(logPath);
  if (idx === -1) return;
  state.jobs.splice(idx, 1);
  delete state.jobStatuses[logPath];
  persistJobs();
  if (state.logPath === logPath && state.topView === "job") {
    if (state.jobs.length) {
      switchToJob(state.jobs[0]);
    } else {
      state.logPath = "";
      localStorage.removeItem(STORAGE_KEY);
      stopPolling();
      showTopView("all-overview");
    }
  } else {
    renderJobTabsBar();
    if (state.topView === "all-overview") renderAllOverview();
  }
}

function switchToJob(logPath) {
  state.logPath = logPath;
  localStorage.setItem(STORAGE_KEY, logPath);
  prepareForNewJob();
  refreshAll(true);
  // Jump back to the Overview tab for the newly-active job, rather than
  // leaving whatever tab was selected for a *previous* job showing that
  // job's now-irrelevant content until the user happens to click it again.
  const overviewBtn = qs('.tab-btn[data-tab="overview"]');
  if (overviewBtn) overviewBtn.click();
  showTopView("job");
}

function showTopView(view) {
  state.topView = view;
  qs("#jobView").classList.toggle("hidden", view !== "job");
  qs("#tab-all-overview").classList.toggle("hidden", view !== "all-overview");
  qs("#tab-all-results").classList.toggle("hidden", view !== "all-results");
  renderJobTabsBar();
  renderActiveJobBanner();
  if (view === "all-overview") renderAllOverview();
  else if (view === "all-results") allResultsView.load();
}

function statusClassFor(logPath) {
  const entry = state.jobStatuses[logPath];
  if (!entry || !entry.status || entry.status.error) return "";
  const status = entry.status;
  if (status.cancelled) return "cancelled";
  if (status.crash && status.crash.crashed) return "crashed";
  if (status.stage === "finished") return "finished";
  return "running";
}

function statusTextFor(logPath) {
  const entry = state.jobStatuses[logPath];
  if (!entry || !entry.status) return "loading…";
  const status = entry.status;
  if (status.error) return "error";
  if (status.cancelled) return "cancelled";
  if (status.crash && status.crash.crashed) return "crashed";
  if (status.stage === "finished") return "finished";
  return STAGE_LABEL[status.stage] || status.stage || "running";
}

function renderJobTabsBar() {
  qsa(".job-tab-btn.special").forEach((btn) => {
    btn.classList.toggle("active", state.topView === btn.dataset.top);
  });
  const wrap = qs("#jobChips");
  wrap.innerHTML = "";
  state.jobs.forEach((logPath) => {
    const chip = document.createElement("button");
    chip.className = "job-chip" + (state.topView === "job" && state.logPath === logPath ? " active" : "");
    chip.innerHTML = `<span class="job-chip-dot ${statusClassFor(logPath)}"></span>` +
      `<span title="${escapeHtml(logPath)}">${escapeHtml(jobLabel(logPath))}</span>` +
      `<span class="job-chip-remove" title="Stop tracking this job">✕</span>`;
    chip.addEventListener("click", (e) => {
      if (e.target.closest(".job-chip-remove")) { removeJob(logPath); return; }
      switchToJob(logPath);
    });
    wrap.appendChild(chip);
  });
}

function renderActiveJobBanner() {
  const el = qs("#activeJobBanner");
  if (state.topView !== "job" || !state.logPath) { el.classList.add("hidden"); return; }
  el.classList.remove("hidden");
  el.innerHTML = `Viewing job: <b>${escapeHtml(jobLabel(state.logPath))}</b> ` +
    `<span class="muted" title="${escapeHtml(state.logPath)}">(${escapeHtml(state.logPath)})</span>`;
}

// Picks which of RFdiffusion's or the modeling stage's own progress
// counters is the relevant "current stage" one to show for a job: while
// a job is actually in that stage its own ETA applies, but once it's
// moved on (mpnn/filtering sit between rfdiffusion and modeling; scoring/
// finished sit after modeling) the counters stay at their final value on
// disk, so showing them still reads as "done" rather than going blank.
function jobOverviewProgress(status) {
  const stage = status.stage;
  if ((stage === "modeling" || stage === "scoring" || stage === "finished") && status.modeling) {
    const m = status.modeling;
    const unit = m.unit ? m.unit.charAt(0).toUpperCase() + m.unit.slice(1) : "Model";
    return {
      label: `${unit}s done`,
      value: `${m.completed} / ${m.expected_total ?? "?"}`,
      pct: m.expected_total ? Math.min(100, Math.round((m.completed / m.expected_total) * 100)) : 0,
      eta: stage === "modeling" ? m.eta_seconds : null,
    };
  }
  if (status.rfdiffusion) {
    const rf = status.rfdiffusion;
    return {
      label: "Backbones done",
      value: `${rf.completed} / ${rf.total ?? "?"}`,
      pct: rf.total ? Math.min(100, Math.round((rf.completed / rf.total) * 100)) : 0,
      eta: stage === "rfdiffusion" ? rf.eta_seconds : null,
    };
  }
  return null;
}

// One card per tracked job, summarizing exactly the numbers the single
// job's own Overview tab computes (same status payload, same "accepted
// backbones × sequences/backbone" planned-designs math) - so switching
// to a job's own Overview tab never shows different math than this
// glanced-at-from-the-list version did.
function renderJobOverviewCard(logPath) {
  const cls = statusClassFor(logPath);
  const badgeClass = cls === "finished" ? "passed" : (cls === "crashed" || cls === "cancelled") ? "failed_filter" : "pending";
  const label = jobLabel(logPath);
  const head = `
    <div class="job-overview-head">
      <span class="job-chip-dot ${cls}"></span>
      <span class="job-overview-name" title="${escapeHtml(logPath)}">${escapeHtml(label)}</span>
      <span class="badge ${badgeClass}">${escapeHtml(statusTextFor(logPath))}</span>
    </div>`;

  const entry = state.jobStatuses[logPath];
  const status = entry && entry.status;
  if (!status || status.error) {
    const msg = status && status.error ? status.error : "Loading…";
    return `<div class="job-overview-card" data-log="${escapeHtml(logPath)}">${head}<p class="muted" style="margin:8px 0 0">${escapeHtml(msg)}</p></div>`;
  }

  const progress = jobOverviewProgress(status);
  const backboneSummary = status.backbones_summary || {};
  const accepted = backboneSummary.accepted; // null until filtering has finished
  const perBackbone = status.config ? status.config.num_seq_per_target_mpnn : null;
  const finalDesigns = (typeof accepted === "number" && typeof perBackbone === "number") ? accepted * perBackbone : null;

  const boxes = [
    progress ? { label: progress.label, value: progress.value } : { label: "Progress", value: "—" },
    { label: "ETA (current stage)", value: progress && progress.eta != null ? fmtSeconds(progress.eta) : "—" },
    { label: "Accepted backbones", value: accepted != null ? `${accepted} / ${backboneSummary.total ?? "?"}` : "pending" },
    { label: "Total designs", value: finalDesigns != null ? finalDesigns.toLocaleString() : "—" },
  ];
  const boxesHtml = boxes.map((b) => `<div class="stat-box"><div class="label">${escapeHtml(b.label)}</div><div class="value small">${escapeHtml(String(b.value))}</div></div>`).join("");
  const bar = progress ? `<div class="progress-bar-outer" style="margin:10px 0 0"><div class="progress-bar-inner" style="width:${progress.pct}%"></div></div>` : "";

  return `<div class="job-overview-card" data-log="${escapeHtml(logPath)}">${head}${bar}<div class="status-grid job-overview-grid">${boxesHtml}</div></div>`;
}

function renderAllOverview() {
  const el = qs("#allOverviewList");
  if (!state.jobs.length) {
    el.innerHTML = `<p class="muted">Add a job above to get started.</p>`;
    return;
  }
  el.innerHTML = state.jobs.map((logPath) => renderJobOverviewCard(logPath)).join("");
  qsa(".job-overview-card", el).forEach((card) => {
    card.addEventListener("click", () => switchToJob(card.dataset.log));
  });
}

function initJobTabsBar() {
  qsa(".job-tab-btn.special").forEach((btn) => {
    btn.addEventListener("click", () => showTopView(btn.dataset.top));
  });
}

// "Run job" vs "Track job" - a split above everything else in this file.
// Track job is the entire existing UI (state.topView/showTopView() and
// everything under it), left completely untouched; this only toggles
// which of the two top-level roots is visible.
function showAppMode(mode) {
  const previousMode = state.appMode;
  state.appMode = mode;
  localStorage.setItem(APP_MODE_KEY, mode);
  qs("#trackJobRoot").classList.toggle("hidden", mode !== "track");
  qs("#runJobRoot").classList.toggle("hidden", mode !== "run");
  qsa(".app-mode-btn").forEach((btn) => btn.classList.toggle("active", btn.dataset.mode === mode));
  if (mode === "run" && !state.runJobInitialized) {
    state.runJobInitialized = true;
    initRunJobTab();
  }
  // Coming back to Track job from Run job, default to the overview across
  // every tracked job rather than leaving whatever single job's tabs (or
  // All jobs results) happened to be showing before the user left for Run
  // job - that's rarely still what they want to see first, especially
  // right after submitting something new. Only on that specific
  // transition, not on every Track job visit (e.g. clicking between Track
  // job's own tabs, or a page reload landing back on "track" via
  // localStorage, shouldn't fight the user's own navigation).
  if (mode === "track" && previousMode === "run") {
    showTopView("all-overview");
  }
}

function initAppModeTabs() {
  qsa(".app-mode-btn").forEach((btn) => {
    btn.addEventListener("click", () => showAppMode(btn.dataset.mode));
  });
  showAppMode(state.appMode);
}

// A lightweight, always-on sweep (independent of the active job's own
// auto-refresh toggle) that keeps every tracked job's status - and hence
// its chip's color dot / the All jobs overview list - current. Skips a
// job once it's reached a terminal state, same reasoning as
// stopPollingForTerminalState() for the active job's own polling: no
// point re-fetching something that can't change anymore.
async function fetchJobStatus(logPath) {
  try {
    const status = await apiGet("/api/status", { log: logPath });
    state.jobStatuses[logPath] = { status, terminal: !status.error && isTerminalStatus(status) };
  } catch (e) {
    state.jobStatuses[logPath] = { status: { error: e.message }, terminal: false };
  }
  renderJobTabsBar();
  if (state.topView === "all-overview") renderAllOverview();
}

function refreshAllJobStatuses() {
  state.jobs.forEach((logPath) => {
    const entry = state.jobStatuses[logPath];
    if (entry && entry.terminal) return;
    fetchJobStatus(logPath);
  });
}

function startJobStatusPolling() {
  refreshAllJobStatuses();
  setInterval(refreshAllJobStatuses, POLL_MS);
}

// ---------------------------------------------------------------------
// Tabs
// ---------------------------------------------------------------------

function initTabs() {
  qsa(".tab-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      qsa(".tab-btn").forEach((b) => b.classList.remove("active"));
      qsa(".tab-panel").forEach((p) => p.classList.remove("active"));
      btn.classList.add("active");
      qs(`#tab-${btn.dataset.tab}`).classList.add("active");
      state.activeTab = btn.dataset.tab;
      refreshActiveTabData();
    });
  });
}

function refreshActiveTabData() {
  if (!state.logPath) return;
  if (state.activeTab === "backbones") loadBackbones();
  else if (state.activeTab === "sequences") loadSequences();
  else if (state.activeTab === "models") loadModels();
  else if (state.activeTab === "results") resultsView.load();
  else if (state.activeTab === "error") loadErrorTab();
  else if (state.activeTab === "outputlog") loadOutputLog();
}

// ---------------------------------------------------------------------
// Cancel job
// ---------------------------------------------------------------------
// Cancelling actually runs scancel against a real cluster, so this is
// deliberately more ceremonious than the rest of the tracking UI: it asks
// which configured target the job's cluster is (the dashboard has no other
// way to know - a job just added by pasting a log path was never
// necessarily submitted through the Run job tab, so there's no target on
// record for it) and requires an explicit "I'm sure" confirmation before
// the request is ever sent. jobId/logPath are captured once when the modal
// opens (rather than re-read from state.lastStatus at submit time) so a
// poll tick landing mid-dialog can't quietly swap out which job the button
// is about to cancel.
let cancelJobCtx = { jobId: null, logPath: null };

function initCancelJobModal() {
  qs("#cancelJobBtn").addEventListener("click", openCancelJobModal);
  qs("#cancelJobCloseBtn").addEventListener("click", closeCancelJobModal);
  qs("#cancelJobBackBtn").addEventListener("click", closeCancelJobModal);
  qs("#cancelJobConfirmCheck").addEventListener("change", refreshCancelJobConfirmEnabled);
  qs("#cancelJobTargetSelect").addEventListener("change", refreshCancelJobConfirmEnabled);
  qs("#cancelJobConfirmBtn").addEventListener("click", doCancelJob);
}

async function openCancelJobModal() {
  if (!state.logPath) return;
  const job = (state.lastStatus && state.lastStatus.job_info) || {};
  cancelJobCtx = { jobId: job.job_id || null, logPath: state.logPath };

  const resultEl = qs("#cancelJobResult");
  resultEl.classList.add("hidden");
  resultEl.textContent = "";
  qs("#cancelJobStatus").textContent = "";
  const checkbox = qs("#cancelJobConfirmCheck");
  checkbox.checked = false;
  checkbox.disabled = false;
  qs("#cancelJobBackBtn").disabled = false;

  const intro = qs("#cancelJobIntro");
  if (cancelJobCtx.jobId) {
    intro.innerHTML = `Cancel SLURM job <b>#${escapeHtml(cancelJobCtx.jobId)}</b> for <b>${escapeHtml(jobLabel(state.logPath))}</b>? This runs <code>scancel</code> on the cluster right away - it cannot be undone.`;
  } else {
    intro.innerHTML = `This job's SLURM job ID hasn't been found in its log yet (it may not have started running), so it can't be cancelled from here yet. Try again once the <b>Overview</b> tab shows a "Slurm job" number.`;
  }
  checkbox.disabled = !cancelJobCtx.jobId;

  const targetSelect = qs("#cancelJobTargetSelect");
  const targetsWarning = qs("#cancelJobTargetsWarning");
  targetSelect.innerHTML = `<option value="">Loading…</option>`;
  targetSelect.disabled = true;
  targetsWarning.classList.add("hidden");
  qs("#cancelJobModal").classList.remove("hidden");

  try {
    const data = await apiGet("/api/run/targets");
    const targets = data.targets || [];
    if (!targets.length) {
      targetSelect.innerHTML = `<option value="">(no targets configured)</option>`;
      targetsWarning.textContent = data.error || "No run targets configured - copy dashboard/dashboard_config.yaml.example to dashboard_config.yaml and fill in at least one target.";
      targetsWarning.classList.remove("hidden");
    } else {
      targetSelect.innerHTML = targets
        .map((t) => `<option value="${escapeHtml(t.name)}" ${t.is_default ? "selected" : ""}>${escapeHtml(t.label)} (${t.kind})</option>`)
        .join("");
      targetSelect.disabled = false;
    }
  } catch (e) {
    targetSelect.innerHTML = `<option value="">(could not load targets)</option>`;
    targetsWarning.textContent = String(e.message || e);
    targetsWarning.classList.remove("hidden");
  }
  refreshCancelJobConfirmEnabled();
}

function closeCancelJobModal() {
  qs("#cancelJobModal").classList.add("hidden");
}

function refreshCancelJobConfirmEnabled() {
  const ready = !!cancelJobCtx.jobId &&
    qs("#cancelJobConfirmCheck").checked &&
    !!qs("#cancelJobTargetSelect").value;
  qs("#cancelJobConfirmBtn").disabled = !ready;
}

async function doCancelJob() {
  const target = qs("#cancelJobTargetSelect").value;
  const jobId = cancelJobCtx.jobId;
  if (!jobId || !target) return;

  const confirmBtn = qs("#cancelJobConfirmBtn");
  const backBtn = qs("#cancelJobBackBtn");
  const statusEl = qs("#cancelJobStatus");
  const resultEl = qs("#cancelJobResult");
  confirmBtn.disabled = true;
  backBtn.disabled = true;
  qs("#cancelJobConfirmCheck").disabled = true;
  qs("#cancelJobTargetSelect").disabled = true;
  statusEl.innerHTML = `<span class="spinner"></span>Cancelling…`;
  resultEl.classList.add("hidden");

  try {
    const res = await fetch("/api/run/cancel", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ target, job_id: jobId }),
    });
    const data = await res.json().catch(() => ({ ok: false, stderr: `HTTP ${res.status}` }));
    statusEl.textContent = "";
    resultEl.classList.remove("hidden");
    if (data.ok) {
      resultEl.className = "ok";
      resultEl.textContent = `Sent scancel for job #${jobId}. It may take a few seconds for the cluster to actually stop the job - the status banners here will update on the next refresh.`;
      if (state.logPath === cancelJobCtx.logPath) refreshAll(false);
    } else {
      resultEl.className = "err";
      resultEl.textContent = data.error || data.stderr || data.stdout || "scancel failed for an unknown reason.";
    }
  } catch (e) {
    statusEl.textContent = "";
    resultEl.classList.remove("hidden");
    resultEl.className = "err";
    resultEl.textContent = String(e.message || e);
  } finally {
    backBtn.disabled = false;
    qs("#cancelJobConfirmCheck").disabled = !cancelJobCtx.jobId;
    qs("#cancelJobTargetSelect").disabled = false;
    refreshCancelJobConfirmEnabled();
  }
}

// ---------------------------------------------------------------------
// Polling
// ---------------------------------------------------------------------

function startPolling() {
  stopPolling();
  state.pollTimer = setInterval(() => refreshAll(false), POLL_MS);
}
function stopPolling() {
  if (state.pollTimer) clearInterval(state.pollTimer);
  state.pollTimer = null;
}

function isTerminalStatus(status) {
  return status.stage === "finished" || !!(status.crash && status.crash.crashed) || !!status.cancelled;
}

function terminalReason(status) {
  if (status.cancelled) return "job was cancelled";
  if (status.crash && status.crash.crashed) return "job crashed";
  return "job finished";
}

function stopPollingForTerminalState(status) {
  if (!state.pollTimer) return; // already stopped (or user stopped it manually)
  stopPolling();
  qs("#autoRefreshToggle").checked = false;
  const note = qs("#autoRefreshNote");
  note.textContent = `(stopped — ${terminalReason(status)})`;
  note.classList.remove("hidden");
}

async function refreshAll(showLoadingState) {
  if (!state.logPath) return;
  await loadStatus(showLoadingState);
  refreshActiveTabData();
}

// ---------------------------------------------------------------------
// Overview / status
// ---------------------------------------------------------------------

async function loadStatus(showLoadingState) {
  const errBanner = qs("#errorBanner");
  try {
    const status = await apiGet("/api/status");
    state.lastStatus = status;
    if (status.error) {
      errBanner.textContent = status.error;
      errBanner.classList.remove("hidden");
      qs("#stepper").innerHTML = "";
      qs("#stageCard").innerHTML = "";
      qs("#overviewCard").innerHTML = `<p class="muted">Fix the log path above and click Load again.</p>`;
      qs("#timingCard").classList.add("hidden");
      return;
    }
    errBanner.classList.add("hidden");
    renderStepper(status);
    renderOverview(status);
    renderModelsCycleNote(status);
    if (isTerminalStatus(status)) stopPollingForTerminalState(status);
  } catch (err) {
    errBanner.textContent = err.message;
    errBanner.classList.remove("hidden");
  }
}

function renderStepper(status) {
  const stage = status.stage;
  const idx = STAGE_ORDER.indexOf(stage);
  const html = STAGE_ORDER.map((s, i) => {
    let cls = "step";
    if (i < idx) cls += " done";
    else if (i === idx) cls += " active";
    const arrow = i < STAGE_ORDER.length - 1 ? '<span class="step-arrow">→</span>' : "";
    return `<div class="${cls}">${STAGE_LABEL[s]}</div>${arrow}`;
  }).join("");
  qs("#stepper").innerHTML = html;
}

function cyclePillHtml(cycle) {
  if (!cycle || !cycle.is_multi_cycle) return "";
  return ` <span class="cycle-pill">Cycle ${cycle.current_cycle_display} / ${cycle.total_cycles}</span>`;
}

function simpleStageMessage(icon, html) {
  return `<div class="stage-message"><span class="stage-message-icon">${icon}</span><p>${html}</p></div>`;
}

// What's happening right now, shown directly under the stepper: either the
// full progress-bar block (RFdiffusion / Modeling) or a short status line
// with an icon for the quicker/simpler stages.
function renderStageBlock(status, cycleSuffix) {
  if (status.stage === "rfdiffusion" && status.rfdiffusion) {
    return renderRfdiffStage(status.rfdiffusion);
  }
  if (status.stage === "modeling" && status.modeling) {
    return renderModelingStage(status.modeling, cycleSuffix);
  }
  if (status.stage === "mpnn") {
    return simpleStageMessage("🧵", `Generating sequences with ProteinMPNN${cycleSuffix} — this step is fast and usually finishes in seconds. Check the <b>Sequences</b> tab once it's done.`);
  }
  if (status.stage === "filtering") {
    return simpleStageMessage("🔍", `Filtering backbones and rebuilding chains before ProteinMPNN — a quick housekeeping step.`);
  }
  if (status.stage === "scoring") {
    return simpleStageMessage(`<span class="spinner"></span>`, `Running final scoring (final_operations)… this produces <code>final_pdbs/</code> and <code>final_output.csv</code>. Check the <b>Results</b> tab once it's done.`);
  }
  if (status.stage === "finished") {
    const fp = (status.scoring && status.scoring.final_pdbs) || [];
    return simpleStageMessage("✅", `Job finished with <b>${fp.length}</b> final model(s). See the <b>Results</b> tab.`);
  }
  return simpleStageMessage("⚙️", `Setting up the run (reading config, preparing RFdiffusion)…`);
}

function renderOverview(status) {
  const card = qs("#overviewCard");
  const cfg = status.config || {};
  const loc = status.location || {};
  const job = status.job_info || {};
  const cycle = status.cycle || {};
  const errFile = status.err_file || {};
  const cycleSuffix = cyclePillHtml(cycle);

  // The "Show full configuration" <details> gets rebuilt below - remember
  // whether it was open so re-rendering (every poll) doesn't collapse it.
  const prevDetails = qs("details", card);
  const wasConfigOpen = prevDetails ? prevDetails.open : false;

  qs("#stageCard").innerHTML = renderStageBlock(status, cycleSuffix);

  const errBlock = (status.possible_errors && status.possible_errors.length)
    ? `<div class="card" style="border-color:var(--danger)"><h4 style="margin-top:0;color:var(--danger)">⚠ Possible errors detected in log</h4><pre style="white-space:pre-wrap;font-size:12px;margin:0">${escapeHtml(status.possible_errors.join("\n"))}</pre></div>`
    : "";

  const errLogValue = errFile.err_exists
    ? `<a class="crash-link" id="gotoErrorTabFromStat">available</a>`
    : `<span class="muted">not found</span>`;

  card.innerHTML = `
    <div class="status-grid">
      <div class="stat-box"><div class="label">Job</div><div class="value small">${escapeHtml(cfg.task_name || "—")}</div></div>
      <div class="stat-box"><div class="label">Prediction model</div><div class="value small">${escapeHtml(cfg.prediction_model || "—")}</div></div>
      <div class="stat-box"><div class="label">Output dir</div><div class="value small" style="word-break:break-all">${escapeHtml(loc.output_dir || "—")}</div></div>
      <div class="stat-box"><div class="label">Current stage</div><div class="value">${STAGE_LABEL[status.stage] || status.stage}${cycleSuffix}</div></div>
      <div class="stat-box"><div class="label">Slurm job</div><div class="value small">${job.job_id ? "#" + escapeHtml(job.job_id) : "—"}</div></div>
      <div class="stat-box"><div class="label">Node</div><div class="value small">${escapeHtml(job.node || "—")}</div></div>
      <div class="stat-box"><div class="label">Started</div><div class="value small">${escapeHtml(job.started_at || "—")}</div></div>
      <div class="stat-box"><div class="label">Error log</div><div class="value small">${errLogValue}</div></div>
    </div>
    ${renderPlanSummary(status)}
    ${errBlock}
    <details style="margin-top:14px">
      <summary style="cursor:pointer;color:var(--muted)">Show full configuration (${Object.keys(cfg).length} keys)</summary>
      <table class="config-table">
        ${renderConfigRows(cfg)}
      </table>
    </details>
  `;

  const newDetails = qs("details", card);
  if (newDetails) newDetails.open = wasConfigOpen;

  const gotoLink = qs("#gotoErrorTabFromStat", card);
  if (gotoLink) gotoLink.addEventListener("click", () => qs('.tab-btn[data-tab="error"]').click());

  renderTimingCard(status.timing);
  renderStatusBanners(status);
}

function renderPlanSummary(status) {
  const cfg = status.config || {};
  const cycle = status.cycle || {};
  const backboneSummary = status.backbones_summary || {};
  const backbones = cfg.num_designs_rfdiff;
  const accepted = backboneSummary.accepted; // null until filtering has actually finished
  const perBackbone = cfg.num_seq_per_target_mpnn;
  const perSeq = cfg.num_models;
  const cycles = cycle.total_cycles || cfg.af2_mpnn_cycles || 1;

  let plannedDesigns = null;
  if (typeof backbones === "number" && typeof perBackbone === "number") {
    plannedDesigns = backbones * perBackbone;
  }
  let finalDesigns = null;
  if (typeof accepted === "number" && typeof perBackbone === "number") {
    finalDesigns = accepted * perBackbone;
  }

  const boxes = [
    { label: "RFdiffusion backbones", value: backbones ?? "—" },
    { label: "Accepted backbones", value: accepted ?? "—" },
    { label: "Sequences / backbone", value: perBackbone ?? "—" },
  ];
  // Only worth a box when the run actually varies from the trivial case.
  if (perSeq != null && perSeq !== 1) boxes.push({ label: "Models / sequence", value: perSeq });
  if (cycles !== 1) boxes.push({ label: "Cycles", value: cycles });
  boxes.push({ label: "Planned designs", value: plannedDesigns != null ? plannedDesigns.toLocaleString() : "—" });
  boxes.push({ label: "Final designs", value: finalDesigns != null ? finalDesigns.toLocaleString() : "—" });

  const boxesHtml = boxes.map((b) => `<div class="stat-box"><div class="label">${escapeHtml(b.label)}</div><div class="value">${escapeHtml(String(b.value))}</div></div>`).join("");
  return `<div class="status-grid" style="margin-top:14px">${boxesHtml}</div>`;
}

function renderConfigRows(cfg) {
  return Object.entries(cfg).map(([k, v]) => {
    if (v !== null && typeof v === "object") {
      return `<tr><td>${escapeHtml(k)}</td><td class="config-nested">${escapeHtml(JSON.stringify(v, null, 2))}</td></tr>`;
    }
    return `<tr><td>${escapeHtml(k)}</td><td>${escapeHtml(String(v))}</td></tr>`;
  }).join("");
}

function renderTimingCard(timing) {
  const card = qs("#timingCard");
  const content = qs("#timingContent");
  if (!timing || !timing.steps || !timing.steps.length) {
    card.classList.add("hidden");
    return;
  }
  card.classList.remove("hidden");
  const rows = timing.steps.map((s) => `<tr><td>${escapeHtml(s.label)}</td><td>${fmtSeconds(s.seconds)}</td></tr>`).join("");
  const totalRow = timing.total_seconds != null
    ? `<tr class="total-row"><td>Total job duration</td><td>${fmtSeconds(timing.total_seconds)}</td></tr>`
    : "";
  content.innerHTML = `<table>${rows}${totalRow}</table>`;
}

// Crash banner, cancellation banner, and the Error log tab's visibility/
// alert state are all driven by the same status fields, so they're kept
// in sync from one place.
function renderStatusBanners(status) {
  const crash = status.crash || {};
  const errFile = status.err_file || {};

  // output_dir_exists is only False once the log itself has told us where
  // output_dir is - a job that's still in "setup" and hasn't reached that
  // config-dump line yet has output_dir_exists === undefined here, not
  // false, so this deliberately checks `=== false` rather than `!`.
  const outputDirBanner = qs("#outputDirBanner");
  if (status.output_dir_exists === false) {
    outputDirBanner.classList.remove("hidden");
    const resolvedPath = (status.location && status.location.output_dir) || "(unknown)";
    outputDirBanner.innerHTML = `<span class="crash-icon">📁</span><div><b>Output directory not reachable from this machine.</b> The resolved path is <code>${escapeHtml(resolvedPath)}</code>, but it doesn't exist here, so Backbones/Sequences/Models/Results will stay empty until it does. If that same path opens fine when you paste it directly into File Explorer, the process running this dashboard likely can't see it even though you can: a mapped network drive is only visible within the Windows session it was connected in, so a dashboard started via Task Scheduler/as a background service may not see a drive mapped in your own interactive login.</div>`;
  } else {
    outputDirBanner.classList.add("hidden");
  }

  const crashBanner = qs("#crashBanner");
  if (crash.crashed) {
    crashBanner.classList.remove("hidden");
    const note = errFile.err_exists
      ? `See the <a class="crash-link" id="gotoErrorTabCrash">Error log</a> tab for the full traceback.`
      : `No matching slurm error (.err) file was found next to the log.`;
    crashBanner.innerHTML = `<span class="crash-icon">⚠</span><div><b>This job appears to have crashed.</b> The log contains "There was an error running the command." ${note}</div>`;
    const link = qs("#gotoErrorTabCrash", crashBanner);
    if (link) link.addEventListener("click", () => qs('.tab-btn[data-tab="error"]').click());
  } else {
    crashBanner.classList.add("hidden");
  }

  const cancelBanner = qs("#cancelBanner");
  if (status.cancelled) {
    cancelBanner.classList.remove("hidden");
    cancelBanner.innerHTML = `<span class="crash-icon">🛑</span><div><b>This job appears to have been cancelled.</b>${status.cancelled_line ? ` <code>${escapeHtml(status.cancelled_line)}</code>` : ""} See the <a class="crash-link" id="gotoErrorTabCancel">Error log</a> tab.</div>`;
    const link = qs("#gotoErrorTabCancel", cancelBanner);
    if (link) link.addEventListener("click", () => qs('.tab-btn[data-tab="error"]').click());
  } else {
    cancelBanner.classList.add("hidden");
  }

  const tabBtn = qs("#errorTabBtn");
  tabBtn.classList.toggle("visible", !!errFile.err_exists);
  tabBtn.classList.toggle("alert", !!crash.crashed || !!status.cancelled);

  updateCancelJobBtn(status);
}

// Cancel job is only meaningful for a job that's actually still running -
// once it's reached a terminal state there's nothing left for scancel to
// do, so the button is disabled (same "is this job done" check already
// used to stop polling - see isTerminalStatus()) rather than left clickable
// but pointless.
function updateCancelJobBtn(status) {
  const btn = qs("#cancelJobBtn");
  if (!btn) return;
  const terminal = isTerminalStatus(status);
  btn.disabled = terminal;
  btn.title = terminal
    ? "This job has already finished, crashed, or been cancelled - nothing to cancel."
    : "Cancel this SLURM job";
}

function renderRfdiffStage(rf) {
  const pct = rf.total ? Math.min(100, Math.round((rf.completed / rf.total) * 100)) : 0;
  return `
    <h4 style="margin-top:0">RFdiffusion backbone generation</h4>
    <div class="progress-bar-outer"><div class="progress-bar-inner" style="width:${pct}%"></div></div>
    <div class="status-grid">
      <div class="stat-box"><div class="label">Backbones</div><div class="value">${rf.completed} / ${rf.total ?? "?"}</div></div>
      <div class="stat-box"><div class="label">Avg time / backbone</div><div class="value small">${fmtSeconds(rf.avg_seconds_per_backbone)}</div></div>
      <div class="stat-box"><div class="label">Currently on</div><div class="value small">${rf.current_index !== null ? `backbone #${rf.current_index}` : "—"}</div></div>
      <div class="stat-box"><div class="label">Est. time remaining</div><div class="value small">${fmtSeconds(rf.eta_seconds)}</div></div>
    </div>
  `;
}

function renderModelingStage(m, cycleSuffix) {
  const pct = m.expected_total ? Math.min(100, Math.round((m.completed / m.expected_total) * 100)) : 0;
  return `
    <h4 style="margin-top:0">Structure modeling (${escapeHtml(m.prediction_model || "")})${cycleSuffix || ""}</h4>
    <div class="progress-bar-outer"><div class="progress-bar-inner" style="width:${pct}%"></div></div>
    <div class="status-grid">
      <div class="stat-box"><div class="label">Completed (${escapeHtml(m.unit)})</div><div class="value">${m.completed} / ${m.expected_total ?? "?"}</div></div>
      <div class="stat-box"><div class="label">Avg time / ${escapeHtml(m.unit)}</div><div class="value small">${fmtSeconds(m.avg_seconds)}</div></div>
      <div class="stat-box"><div class="label">Currently on</div><div class="value small" style="word-break:break-all">${m.current_name ? escapeHtml(m.current_name) : "—"}</div></div>
      <div class="stat-box"><div class="label">Est. time remaining</div><div class="value small">${fmtSeconds(m.eta_seconds)}</div></div>
    </div>
    <p class="muted" style="margin-top:10px;font-size:12px">
      Each unit is one predicted structure (one designed sequence, or - for Boltz - one of its diffusion samples), counted straight off disk as soon as it's written, including monomer predictions when the run models those too.
      ${cycleSuffix ? " Completed/total reflect the current cycle only (earlier cycles' models are cleared from disk once the next cycle starts); the average time and ETA are estimated from timing across the whole log." : ""}
    </p>
  `;
}

// ---------------------------------------------------------------------
// Error log tab
// ---------------------------------------------------------------------

async function loadErrorTab() {
  const el = qs("#errorTabContent");
  const status = state.lastStatus;
  const errFile = status && status.err_file;
  if (!errFile || !errFile.err_exists) {
    el.innerHTML = `<p class="muted">No slurm error file found${errFile && errFile.err_path ? ` (looked for <code>${escapeHtml(errFile.err_path)}</code>)` : ""}.</p>`;
    return;
  }
  el.innerHTML = `<p class="muted">Loading…</p>`;
  try {
    const data = await apiGet("/api/error_log");
    if (!data.err_exists) {
      el.innerHTML = `<p class="muted">No slurm error file found at <code>${escapeHtml(data.err_path || "?")}</code>.</p>`;
      return;
    }
    let banners = "";
    if (status.crash && status.crash.crashed) {
      banners += `<div class="alert-banner"><span class="crash-icon">⚠</span><div><b>This job crashed.</b> The log contains "There was an error running the command."</div></div>`;
    }
    if (status.cancelled) {
      banners += `<div class="alert-banner"><span class="crash-icon">🛑</span><div><b>This job was cancelled.</b>${status.cancelled_line ? ` <code>${escapeHtml(status.cancelled_line)}</code>` : ""}</div></div>`;
    }
    const trunc = data.truncated
      ? `<p class="muted">(truncated — showing the first ${data.content.length.toLocaleString()} of ${data.size.toLocaleString()} bytes)</p>`
      : "";
    el.innerHTML = `
      ${banners}
      <p class="muted">Contents of <code>${escapeHtml(data.err_path)}</code>:</p>
      ${trunc}
      <pre class="log-pre">${escapeHtml(data.content)}</pre>
    `;
  } catch (err) {
    el.innerHTML = `<p class="muted">${escapeHtml(err.message)}</p>`;
  }
}

// ---------------------------------------------------------------------
// Output log tab
// ---------------------------------------------------------------------

async function loadOutputLog() {
  const el = qs("#outputLogContent");
  const prevPre = qs("#outputLogPre", el);
  // "Smart tail": if the user was scrolled to the bottom (or this is the
  // first load), keep following the end of the log on every refresh -
  // otherwise leave their scroll position alone, the same way `tail -f`
  // stops auto-scrolling once you scroll up to read something.
  const wasAtBottom = !prevPre || (prevPre.scrollTop + prevPre.clientHeight >= prevPre.scrollHeight - 20);
  const prevScrollTop = prevPre ? prevPre.scrollTop : 0;
  try {
    const data = await apiGet("/api/output_log");
    if (data.error) throw new Error(data.error);
    const trunc = data.truncated
      ? `<p class="muted">(truncated — showing the first ${data.content.length.toLocaleString()} of ${data.size.toLocaleString()} bytes)</p>`
      : "";
    el.innerHTML = `
      <p class="muted">Contents of <code>${escapeHtml(state.logPath)}</code>:</p>
      ${trunc}
      <pre class="log-pre" id="outputLogPre">${escapeHtml(data.content)}</pre>
    `;
    const pre = qs("#outputLogPre", el);
    pre.scrollTop = wasAtBottom ? pre.scrollHeight : prevScrollTop;
  } catch (err) {
    el.innerHTML = `<p class="muted">${escapeHtml(err.message)}</p>`;
  }
}

// ---------------------------------------------------------------------
// Backbones tab
// ---------------------------------------------------------------------

function applyBackboneSelectionHighlight() {
  qsa("#backbonesList .list-item").forEach((el) => {
    el.classList.toggle("selected", el.dataset.key === state.selectedBackboneKey);
  });
}

const BACKBONE_BADGE_LABEL = { passed: "ok", pending: "pending", failed_filter: "filtered" };

function renderBackbonesFilterNote() {
  const el = qs("#backbonesFilterNote");
  const anyPending = (state.backbonesCache || []).some((b) => b.status === "pending");
  el.classList.toggle("hidden", !anyPending);
  if (anyPending) {
    el.textContent = `Filtering hasn't finished yet, so these backbones haven't been evaluated — none are tagged "ok" until it does.`;
  }
}

function renderBackbonesList() {
  const listEl = qs("#backbonesList");
  renderBackbonesFilterNote();
  rerenderPreservingScroll(listEl, () => {
    const all = state.backbonesCache || [];
    const backbones = state.showFilteredBackbones ? all : all.filter((b) => b.status !== "failed_filter");
    if (!backbones.length) {
      listEl.innerHTML = all.length
        ? `<p class="muted">All backbones were filtered out — toggle "Show filtered out" to see them.</p>`
        : `<p class="muted">No backbones generated yet.</p>`;
      return;
    }
    listEl.innerHTML = "";
    backbones.forEach((b) => {
      // A backbone only ever lives in one place at a time (1_rfdiff/ or
      // 1_rfdiff/failed_filters/), so its name alone is a stable key -
      // deliberately *not* including status, since status can legitimately
      // change (pending -> passed) for the same backbone between polls,
      // and that shouldn't look like a different backbone got selected.
      const key = b.name;
      const div = document.createElement("div");
      div.className = "list-item" + (state.selectedBackboneKey === key ? " selected" : "");
      div.dataset.key = key;
      div.innerHTML = `<span>${escapeHtml(b.name)}</span><span class="badge ${b.status}">${BACKBONE_BADGE_LABEL[b.status] || b.status}</span>`;
      div.addEventListener("click", () => selectBackbone(b));
      listEl.appendChild(div);
    });
    const current = backbones.find((b) => b.name === state.selectedBackboneKey) || backbones[0];
    if (current) selectBackbone(current);
  });
}

async function loadBackbones() {
  const listEl = qs("#backbonesList");
  try {
    state.backbonesCache = await apiGet("/api/backbones");
    renderBackbonesList();
  } catch (err) {
    listEl.innerHTML = `<p class="muted">${escapeHtml(err.message)}</p>`;
  }
}

async function selectBackbone(b) {
  const key = b.name;
  state.selectedBackboneKey = key;
  applyBackboneSelectionHighlight();
  // Same structure already showing - don't recreate the 3D viewer (that
  // would discard the user's current zoom/rotation) just because a poll
  // tick re-rendered the list.
  if (state.loadedBackboneKey === key) return;
  state.loadedBackboneKey = key;
  qs("#backboneViewerLabel").textContent = `${b.name} (${b.status})`;
  const url = `/api/backbone_pdb?log=${encodeURIComponent(state.logPath)}&name=${encodeURIComponent(b.name)}&status=${b.status}`;
  try {
    const res = await fetch(url);
    if (!res.ok) throw new Error(await fetchErrorMessage(res));
    const pdbText = await res.text();
    await renderMol("#backboneViewer", pdbText);
    if (state.backboneColorMode === "provenance") {
      await applyProvenanceColoring("#backboneViewer", b.trb_path, qs("#backboneColorLegend"));
    }
  } catch (e) {
    qs("#backboneViewer").innerHTML = `<p class="muted" style="padding:20px">Could not load structure: ${escapeHtml(e.message)}</p>`;
  }
}

const RESIDUE_HIGHLIGHT_RADIUS = 5; // Angstroms - "surrounding residues" cutoff

// ---------------------------------------------------------------------
// Viewer settings (background color, chain color palette, provenance
// category colors) - global, shared across all three structure viewers,
// adjustable from each viewer's own color toolbar (Backbones/Models/
// Results tabs) and remembered in local storage.
// ---------------------------------------------------------------------

const VIEWER_SETTINGS_KEY = "prosculpt_dashboard_viewer_settings";

const BACKGROUND_PRESETS = { black: "#05070c", navy: "#0f1420", white: "#ffffff", light: "#e8ebf2" };

// Chain colors: NGL's built-in "chainid" scheme hashes the chain letter
// into a color that can land on muddy, hard-to-see-on-black tones, so
// these are curated instead. "cvd" is the Okabe-Ito colorblind-safe set.
// Each palette is ordered so its first two entries - by far the most
// common case, a 2-chain binder/receptor design - contrast strongly with
// each other, not just technically differ.
const CHAIN_PALETTES = {
  vivid: ["#5b9bff", "#ff6b6b", "#4ecdc4", "#ffd166", "#c77dff", "#f4a261", "#06d6a0", "#f72585"],
  pastel: ["#8ecae6", "#ffb4a2", "#b8f2e6", "#fff3b0", "#d0bfff", "#ffd6a5", "#a0e7b5", "#ffc2e2"],
  cvd: ["#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"],
  bright: ["#00e5ff", "#ff1744", "#76ff03", "#ffea00", "#d500f9", "#ff9100", "#00e676", "#f50057"],
  sunset: ["#ff5e5b", "#fcbf49", "#e63946", "#ff9f1c", "#ffbf69", "#d62828", "#f77f00", "#ffa5ab"],
  ocean: ["#00b4d8", "#7209b7", "#48cae4", "#3a86ff", "#80ffdb", "#4895ef", "#90e0ef", "#4cc9f0"],
};

const DEFAULT_PROVENANCE_COLORS = { motif: "#ffd166", fixed_chain: "#8b93a7", inpainted_seq: "#ef476f", sculpted: "#06d6a0" };
const HEX_COLOR_RE = /^#[0-9a-fA-F]{6}$/;

let viewerSettings = {
  background: "black",
  chainPalette: "vivid",
  // Per-chain-index color overrides on top of the chosen palette - set by
  // clicking a chain's own swatch in the legend, see
  // wireChainLegendColorInputs(). Keyed by chainIndex (the same number
  // chainColorScheme() indexes the palette array with), not by chain
  // letter, since letters aren't guaranteed to mean the same thing across
  // different structures the way "first chain in the file" reliably is.
  chainColorOverrides: {},
  provenanceColors: { ...DEFAULT_PROVENANCE_COLORS },
  // "selected" - no bulk sidechain representation, just whatever a direct
  // residue click/sequence-panel selection already shows via
  // setupResidueInteraction(). "interface" and "all" are mutually
  // exclusive with each other (and with "selected") by construction, one
  // <input type="radio"> group per viewer rather than two independent
  // checkboxes - showing every sidechain and showing only the interface
  // ones don't make sense to have on at once, so there was never a real
  // reason to let them be.
  sidechainMode: "selected", // "selected" | "interface" | "all" | "hydrophobic"
};
const SIDECHAIN_MODES = ["selected", "interface", "all", "hydrophobic"];

function loadViewerSettings() {
  try {
    const stored = JSON.parse(localStorage.getItem(VIEWER_SETTINGS_KEY) || "{}");
    if (BACKGROUND_PRESETS[stored.background]) viewerSettings.background = stored.background;
    if (CHAIN_PALETTES[stored.chainPalette]) viewerSettings.chainPalette = stored.chainPalette;
    if (stored.chainColorOverrides && typeof stored.chainColorOverrides === "object") {
      Object.entries(stored.chainColorOverrides).forEach(([idx, v]) => {
        if (/^\d+$/.test(idx) && typeof v === "string" && HEX_COLOR_RE.test(v)) viewerSettings.chainColorOverrides[idx] = v;
      });
    }
    if (stored.provenanceColors && typeof stored.provenanceColors === "object") {
      Object.keys(DEFAULT_PROVENANCE_COLORS).forEach((cat) => {
        const v = stored.provenanceColors[cat];
        if (typeof v === "string" && HEX_COLOR_RE.test(v)) viewerSettings.provenanceColors[cat] = v;
      });
    }
    if (SIDECHAIN_MODES.includes(stored.sidechainMode)) viewerSettings.sidechainMode = stored.sidechainMode;
  } catch (e) {}
}
function saveViewerSettings() {
  localStorage.setItem(VIEWER_SETTINGS_KEY, JSON.stringify(viewerSettings));
}

// Registered once; reads the *current* palette (and any per-chain
// overrides) live via closure each time it colors an atom, so switching
// palettes or picking a custom chain color never needs re-registering -
// just re-coloring (see applyChainPaletteSetting()).
let chainColorSchemeId = null;
function chainColorScheme() {
  if (chainColorSchemeId) return chainColorSchemeId;
  chainColorSchemeId = NGL.ColormakerRegistry.addScheme(function () {
    this.atomColor = function (atom) {
      const override = viewerSettings.chainColorOverrides[atom.chainIndex];
      if (override) return parseInt(override.slice(1), 16);
      const palette = CHAIN_PALETTES[viewerSettings.chainPalette] || CHAIN_PALETTES.vivid;
      return parseInt(palette[atom.chainIndex % palette.length].slice(1), 16);
    };
  });
  return chainColorSchemeId;
}

// Which chains are actually present in a loaded structure, and where each
// one starts - used both to render a "which color is which chain" legend
// (chainIndex, matching what chainColorScheme() keys off) and to normalize
// RFdiffusion provenance residue numbers onto whatever numbering
// convention the loaded structure happens to use (see
// makeProvenanceColorScheme() below).
function computeStructureChainInfo(structure) {
  const info = {};
  structure.eachAtom((atom) => {
    const c = atom.chainname;
    if (!(c in info)) info[c] = { minResno: atom.resno, chainIndex: atom.chainIndex };
    else if (atom.resno < info[c].minResno) info[c].minResno = atom.resno;
  });
  return info;
}

// Interface residues: any residue with at least one atom within
// INTERFACE_RADIUS Å of an atom belonging to a *different* chain - one
// chain at a time, since NGL's selection language has no "different
// chain" concept to hand it directly. For each chain: everything within
// the radius of "not this chain" (getAtomSetWithinSelection) necessarily
// includes some of the other chains' own atoms too, which
// makeIntersection() with "this chain's atoms" (getAtomSet) discards,
// leaving just this chain's side of the interface; getAtomSetWithinGroup()
// then expands that atom-level set to whole residues, same as the
// residue-click highlight's neighbor search uses elsewhere in this file.
// Returns an NGL selection string (all chains' interface residues ORed
// together), or null for a single-chain structure, where "interface"
// isn't a meaningful concept.
const INTERFACE_RADIUS = 5; // Angstroms

function computeInterfaceSele(structure) {
  const chainNames = Object.keys(computeStructureChainInfo(structure));
  if (chainNames.length < 2) return null;
  const parts = [];
  chainNames.forEach((chainName) => {
    try {
      const chainAtoms = structure.getAtomSet(new NGL.Selection(`:${chainName}`));
      const nearOtherChains = structure.getAtomSetWithinSelection(new NGL.Selection(`not :${chainName}`), INTERFACE_RADIUS);
      const interfaceAtoms = chainAtoms.makeIntersection(nearOtherChains);
      if (interfaceAtoms.getSize() > 0) {
        parts.push(structure.getAtomSetWithinGroup(interfaceAtoms).toSeleString());
      }
    } catch (e) {
      // Best-effort, same reasoning as the residue-click highlight.
    }
  });
  return parts.length ? parts.join(" or ") : null;
}

// Each swatch is itself an <input type="color"> (not just a static color
// chip) so a chain's color can be overridden directly, same idea as the
// provenance legend's per-category pickers - see wireChainLegendColorInputs().
function renderChainLegend(structure) {
  const info = computeStructureChainInfo(structure);
  const palette = CHAIN_PALETTES[viewerSettings.chainPalette] || CHAIN_PALETTES.vivid;
  const chains = Object.entries(info).sort((a, b) => a[1].chainIndex - b[1].chainIndex);
  return chains.map(([name, { chainIndex }]) => {
    const color = viewerSettings.chainColorOverrides[chainIndex] || palette[chainIndex % palette.length];
    return `<span class="legend-item"><input type="color" class="legend-color-input" data-chain-index="${chainIndex}" value="${color}" title="Change chain ${escapeHtml(name)}'s color">Chain ${escapeHtml(name)}</span>`;
  }).join("");
}

function wireChainLegendColorInputs(legendEl) {
  qsa(".legend-color-input", legendEl).forEach((input) => {
    input.addEventListener("input", (e) => {
      viewerSettings.chainColorOverrides[e.target.dataset.chainIndex] = e.target.value;
      saveViewerSettings();
      applyChainColorLive(legendEl);
    });
  });
}

// Recolors every chain-mode viewer after a per-chain override changes,
// without rebuilding `skipLegendEl`'s own DOM. Rebuilding it (the old
// behavior, via applyChainPaletteSetting()) would destroy the very
// <input type="color"> the user is still dragging in, which abruptly
// closes the native color picker on every pixel of movement - the
// RFdiffusion provenance legend's picker never has this problem, since
// its own color-change handler never rebuilds its legend's DOM either.
// Every *other* currently-open chain-mode viewer's legend still gets
// rebuilt, since a shared per-chain-index override can affect more than
// one of them at once.
function applyChainColorLive(skipLegendEl) {
  Object.keys(viewerRegistry).forEach((selector) => {
    const entry = viewerRegistry[selector];
    if (!entry || !entry.component || entry.colorMode === "provenance") return;
    setCartoonColor(selector, chainColorScheme());
    const legendSel = VIEWER_LEGEND_IDS[selector];
    const legendEl = legendSel && qs(legendSel);
    if (legendEl && legendEl !== skipLegendEl) {
      legendEl.classList.remove("hidden");
      legendEl.innerHTML = renderChainLegend(entry.component.structure);
      wireChainLegendColorInputs(legendEl);
    }
  });
  syncColorControlSelects(); // flips the palette select(s) to "Custom"
}

// Maps each structure viewer to the color-legend element next to its
// color-mode toolbar, so a shared helper (refreshChainLegendFor) can
// update whichever one applies without each caller needing to know.
const VIEWER_LEGEND_IDS = {
  "#backboneViewer": "#backboneColorLegend",
  "#resultViewer": "#resultsColorLegend",
  "#modelViewer": "#modelColorLegend",
  "#allResultViewer": "#allResultsColorLegend",
  "#runJobPdbViewer": "#runJobPdbColorLegend",
};

function refreshChainLegendFor(selector) {
  const entry = viewerRegistry[selector];
  const legendSel = VIEWER_LEGEND_IDS[selector];
  if (!entry || !entry.component || !legendSel) return;
  const legendEl = qs(legendSel);
  if (!legendEl) return;
  legendEl.classList.remove("hidden");
  legendEl.innerHTML = renderChainLegend(entry.component.structure);
  wireChainLegendColorInputs(legendEl);
}

// RFdiffusion residue provenance coloring ("color by .trb"): every residue
// in a generated structure is either part of a grafted motif (present in a
// redesigned chain but taken from the reference), part of an entirely
// fixed/non-designed chain (further split into "fully fixed" vs. "fixed
// structure, redesigned sequence"), or generated de novo. See
// parser.load_trb_provenance() for where con_hal_idx0 / inpaint_str /
// inpaint_seq come from, and why categories arrive as a flat array indexed
// by residue position rather than by (chain, resnum).
const PROVENANCE_LABELS = { motif: "Motif", fixed_chain: "Fixed chains", inpainted_seq: "Inpainted sequence", sculpted: "Sculpted" };

// Unlike the chain palette (one global scheme that reads shared settings
// live), each structure has its own provenance data, so a fresh scheme is
// registered per use rather than trying to share/mutate one - keeps two
// viewers showing different structures' provenance from stepping on each
// other. `categories[atom.residueIndex]` works directly, with no
// chain-letter/resnum lookup at all: parser.load_trb_provenance() already
// returns one category per residue in the generated structure's own file
// order, and NGL's atom.residueIndex is exactly that same flat, 0-based,
// file-order position within whichever structure is currently loaded -
// this is what sidesteps RFdiffusion's raw backbone .pdb (continuous
// resnums across chains) vs. AlphaFold3/Boltz's output (resnums reset per
// chain) ever needing to be reconciled, and also survives
// rechain_rfdiff_pdbs() reassigning chain letters after RFdiffusion runs.
function makeProvenanceColorScheme(categories) {
  return NGL.ColormakerRegistry.addScheme(function () {
    this.atomColor = function (atom) {
      const category = (categories && categories[atom.residueIndex]) || "sculpted";
      return parseInt(viewerSettings.provenanceColors[category].slice(1), 16);
    };
  });
}

// Color swatches double as <input type="color"> pickers here so a user
// can retune the provenance palette; reads viewerSettings.provenanceColors
// live for its initial value each time it's (re-)rendered, and the picked
// color takes effect via wireProvenanceLegendColorInputs()/
// applyProvenanceColorSetting() rather than by re-rendering this markup.
function renderProvenanceLegend() {
  return Object.keys(PROVENANCE_LABELS).map((cat) =>
    `<span class="legend-item"><input type="color" class="legend-color-input" data-category="${cat}" value="${viewerSettings.provenanceColors[cat]}" title="Change ${escapeHtml(PROVENANCE_LABELS[cat])} color">${PROVENANCE_LABELS[cat]}</span>`
  ).join("");
}

function wireProvenanceLegendColorInputs(legendEl) {
  qsa(".legend-color-input", legendEl).forEach((input) => {
    input.addEventListener("input", (e) => {
      viewerSettings.provenanceColors[e.target.dataset.category] = e.target.value;
      saveViewerSettings();
      applyProvenanceColorSetting();
    });
  });
}

// One NGL Stage per viewer container (#backboneViewer / #modelViewer /
// #resultViewer), tracked so a later renderMol() call on the same
// container disposes the old one first - each Stage owns its own WebGL
// context, and browsers cap how many can be alive at once, so leaving old
// ones around would eventually break rendering entirely. `generation`
// guards against a slower, now-stale load (e.g. from rapid-fire clicking
// through the list) finishing after a newer one and clobbering it.
// `lastText`/`lastFormat` cache what's currently loaded so a viewer
// setting change (e.g. chain palette) can re-render without re-fetching.
// `colorMode` ("chain" | "provenance") tracks which scheme is currently
// applied, so a global chain-palette change doesn't clobber a viewer
// that's deliberately showing provenance coloring.
const viewerRegistry = {};

async function renderMol(selector, structureText, format = "pdb") {
  const el = qs(selector);
  if (!viewerRegistry[selector]) viewerRegistry[selector] = { stage: null, generation: 0 };
  const entry = viewerRegistry[selector];
  const myGeneration = ++entry.generation;

  if (entry.stage) {
    entry.stage.dispose();
    entry.stage = null;
  }
  el.innerHTML = "";

  const stage = new NGL.Stage(el, { backgroundColor: BACKGROUND_PRESETS[viewerSettings.background] || BACKGROUND_PRESETS.black });
  // NGL wraps its canvas in a plain position:relative div sized with
  // literal pixel width/height (not percentages) - as an ordinary in-flow
  // child, that div's own size normally feeds right back into el's size
  // whenever el doesn't have a fixed height of its own, which .mol-viewer
  // (flex-basis:auto, so its own resize handle can work - see style.css)
  // doesn't until the user has actually dragged it. Left alone, that's a
  // runaway loop: el sizes to fit the wrapper, ensureResizeObserver()
  // below sees el's new size and resizes the canvas to match, the
  // (still in-flow) wrapper grows to fit *that*, and so on - which is
  // exactly what made a viewer keep growing on its own after a tab
  // became visible again, with no drag involved. Taking the wrapper out
  // of normal flow entirely (absolute + inset:0) breaks the loop for
  // good, regardless of whatever pixel size NGL sets it to internally.
  if (el.firstElementChild) {
    el.firstElementChild.style.position = "absolute";
    el.firstElementChild.style.inset = "0";
  }
  const blob = new Blob([structureText], { type: "text/plain" });

  let component;
  try {
    component = await stage.loadFile(blob, { ext: format });
  } catch (e) {
    if (myGeneration === entry.generation) {
      el.innerHTML = `<p class="muted" style="padding:20px">Could not render structure.</p>`;
    }
    stage.dispose();
    return null;
  }
  if (myGeneration !== entry.generation) {
    stage.dispose(); // a newer selection started before this one finished loading
    return null;
  }

  entry.stage = stage;
  entry.component = component;
  entry.lastText = structureText;
  entry.lastFormat = format;
  entry.colorMode = "chain";
  entry.provenanceSchemeId = null;
  // NGL's own default tooltip re-shows itself (resets display:block) on
  // every hover, so a one-time style change doesn't stick - detach it
  // from the DOM instead. It's still safe for NGL to keep writing to it
  // internally; it just has nothing to render. Our own app-styled
  // tooltip (.residue-tooltip, wired up in setupResidueInteraction)
  // replaces it.
  if (stage.tooltip && stage.tooltip.parentNode) stage.tooltip.parentNode.removeChild(stage.tooltip);
  entry.cartoonRepr = component.addRepresentation("cartoon", { color: chainColorScheme() });
  entry.sidechainRepr = null;
  setSidechainMode(selector, viewerSettings.sidechainMode);
  stage.autoView();
  entry.interaction = setupResidueInteraction(stage, component, el);
  refreshChainLegendFor(selector);
  ensureResizeObserver(selector, el);
  return stage;
}

// .mol-viewer has a native CSS "resize: vertical" handle (drag the
// bottom-right corner) so the user can make a viewer bigger/smaller, but
// NGL only re-fits its canvas on the browser *window* resizing, not on
// its own container being resized some other way - without this, a
// manually-resized viewer would keep rendering at its old size, stretched
// or clipped to fill the new box. Attached once per selector (not once
// per renderMol() call, which would happen every time a new structure
// loads) since the callback always looks up the *current* stage by
// selector rather than closing over one particular Stage instance, so it
// keeps working across every structure that selector ever loads.
function ensureResizeObserver(selector, el) {
  const entry = viewerRegistry[selector];
  if (entry.resizeObserverAttached) return;
  entry.resizeObserverAttached = true;
  // ResizeObserver always fires once immediately after observe() starts,
  // reporting the element's *current* size even though nothing has
  // actually changed yet - reacting to that first call would run
  // handleResize() right on top of the stage.autoView() that already
  // just fit the camera to the freshly-loaded structure, which is
  // redundant at best. Only genuine subsequent size changes matter here.
  let firstCall = true;
  const ro = new ResizeObserver(() => {
    if (firstCall) { firstCall = false; return; }
    const current = viewerRegistry[selector];
    if (current && current.stage) current.stage.handleResize();
  });
  ro.observe(el);
}

// Swaps the cartoon's color scheme on an already-loaded viewer in place -
// cheap (no reload/re-fetch), and doesn't reset the user's camera
// position the way rebuilding the whole viewer would.
function setCartoonColor(selector, colorSchemeId) {
  const entry = viewerRegistry[selector];
  if (!entry || !entry.component) return;
  if (entry.cartoonRepr) entry.component.removeRepresentation(entry.cartoonRepr);
  entry.cartoonRepr = entry.component.addRepresentation("cartoon", { color: colorSchemeId });
}

// "Show sidechains: All" adds every residue's sidechain as thin,
// element-colored (CPK) licorice sticks alongside the cartoon - a lighter
// radiusScale than the click-highlight's sticks (see
// setupResidueInteraction) so this stays a background detail layer rather
// than competing with a deliberate residue highlight/selection. The
// selection includes the alpha carbon (NGL's ".CA" atom-name syntax)
// alongside "sidechain" itself - NGL's plain "sidechain" keyword excludes
// CA, which would draw each sidechain as a stick floating just short of
// the cartoon ribbon instead of actually connecting to it.
// "Interface" restricts that same treatment to computeInterfaceSele()'s
// residues; a single-chain structure has no interface there, so it's a
// silent no-op rather than an error for e.g. a monomer prediction.
// "Hydrophobic" restricts it to NGL's own built-in "hydrophobic" residue
// class (ALA/VAL/LEU/ILE/PRO/PHE/MET/TRP) - no separate residue list
// needed here, NGL's selection language already defines this keyword.
// "Selected" removes the representation entirely - a direct residue
// click/sequence-panel selection already shows its own sidechains via
// setupResidueInteraction(), independent of this.
function setSidechainMode(selector, mode) {
  const entry = viewerRegistry[selector];
  if (!entry || !entry.component) return;
  if (entry.sidechainRepr) {
    entry.component.removeRepresentation(entry.sidechainRepr);
    entry.sidechainRepr = null;
  }
  if (mode === "all") {
    entry.sidechainRepr = entry.component.addRepresentation("licorice", {
      sele: "sidechain or .CA", colorScheme: "element", radiusScale: 0.7,
    });
  } else if (mode === "interface") {
    const interfaceSele = computeInterfaceSele(entry.component.structure);
    if (interfaceSele) {
      entry.sidechainRepr = entry.component.addRepresentation("licorice", {
        sele: `(${interfaceSele}) and (sidechain or .CA)`, colorScheme: "element", radiusScale: 0.7,
      });
    }
  } else if (mode === "hydrophobic") {
    entry.sidechainRepr = entry.component.addRepresentation("licorice", {
      sele: "hydrophobic and (sidechain or .CA)", colorScheme: "element", radiusScale: 0.7,
    });
  }
}

// Background is a cheap, genuinely live stage parameter - applies
// instantly to every currently-open viewer, no reload needed.
function applyBackgroundSetting() {
  const hex = BACKGROUND_PRESETS[viewerSettings.background] || BACKGROUND_PRESETS.black;
  Object.values(viewerRegistry).forEach((entry) => {
    if (entry.stage) entry.stage.setParameters({ backgroundColor: hex });
  });
}

// Shared by every viewer's "Show sidechains" radio group (kept in sync
// via the same class-based delegated-listener pattern as the background
// and chain-palette selects - see initColorControls()).
function applySidechainModeSetting() {
  Object.keys(viewerRegistry).forEach((selector) => setSidechainMode(selector, viewerSettings.sidechainMode));
}

// chainColorScheme() reads viewerSettings.chainPalette live, so any
// viewer currently colored "by chain" just needs its cartoon
// representation (and legend) rebuilt to re-evaluate colors - viewers
// showing provenance coloring are left alone.
function applyChainPaletteSetting() {
  Object.keys(viewerRegistry).forEach((selector) => {
    const entry = viewerRegistry[selector];
    if (entry && entry.colorMode !== "provenance") {
      setCartoonColor(selector, chainColorScheme());
      refreshChainLegendFor(selector);
    }
  });
}

// makeProvenanceColorScheme() reads viewerSettings.provenanceColors live,
// so a viewer already showing provenance coloring just needs its cartoon
// representation rebuilt (via its previously-registered scheme id) to
// pick up a retuned category color - no need to re-fetch the .trb or
// rebuild the provenance map itself.
function applyProvenanceColorSetting() {
  Object.keys(viewerRegistry).forEach((selector) => {
    const entry = viewerRegistry[selector];
    if (entry && entry.colorMode === "provenance" && entry.provenanceSchemeId) {
      setCartoonColor(selector, entry.provenanceSchemeId);
    }
  });
}

// --- Viewer color controls (background / chain palette selects) ---
// Each structure-viewing tab (Backbones, Models, Results, All jobs
// results) has its own copy of the background and chain-palette selects,
// but they all drive the same global `viewerSettings` - one delegated
// listener handles every copy via a shared class, and
// syncColorControlSelects() keeps whichever ones aren't currently visible
// in sync so they don't show a stale value when their tab is switched
// back to.

function syncColorControlSelects() {
  qsa(".bg-color-select").forEach((el) => { el.value = viewerSettings.background; });
  // Any per-chain override in play means the palette select can no longer
  // honestly claim to be "Vivid"/"Pastel"/etc - it shows "Custom" instead
  // until the user picks a real palette again (which clears the overrides
  // - see the change handler below).
  const hasOverrides = Object.keys(viewerSettings.chainColorOverrides).length > 0;
  qsa(".chain-palette-select").forEach((el) => { el.value = hasOverrides ? "custom" : viewerSettings.chainPalette; });
  // Each radio's own value ("selected"/"interface"/"all") is what's
  // authoritative, so this just checks whichever one matches the current
  // mode - one shared value across however many of these radio groups
  // (one per viewer tab) are on the page at once.
  qsa(".sidechain-mode-radio").forEach((el) => { el.checked = el.value === viewerSettings.sidechainMode; });
}

function initColorControls() {
  syncColorControlSelects();
  document.addEventListener("change", (e) => {
    if (e.target.matches(".bg-color-select")) {
      viewerSettings.background = e.target.value;
      saveViewerSettings();
      syncColorControlSelects();
      applyBackgroundSetting();
    } else if (e.target.matches(".chain-palette-select")) {
      // "custom" is a disabled option (see index.html) so the user can
      // never actually pick it - this only ever fires with a real
      // palette name, which overwrites any per-chain overrides, same
      // reasoning as picking a whole new palette from scratch.
      viewerSettings.chainPalette = e.target.value;
      viewerSettings.chainColorOverrides = {};
      saveViewerSettings();
      syncColorControlSelects();
      applyChainPaletteSetting();
    } else if (e.target.matches(".sidechain-mode-radio")) {
      viewerSettings.sidechainMode = e.target.value;
      saveViewerSettings();
      syncColorControlSelects();
      applySidechainModeSetting();
    }
  });
}

// Fetches a .trb file's provenance data and colors `selector`'s cartoon by
// it, updating the legend; falls back to chain coloring (with an
// explanatory note in the legend slot) if there's no .trb path or the
// fetch fails. Shared by the Backbones and Results tabs' "color residues
// by" controls.
async function applyProvenanceColoring(selector, trbPath, legendEl, logOverride) {
  const entry = viewerRegistry[selector];
  if (!entry) return;
  if (!trbPath) {
    entry.colorMode = "chain";
    entry.provenanceSchemeId = null;
    setCartoonColor(selector, chainColorScheme());
    if (legendEl) legendEl.innerHTML = `<span class="muted">No .trb file found for this backbone.</span>`;
    return;
  }
  try {
    const data = await apiGet("/api/trb", { path: trbPath, log: logOverride });
    if (data.error) throw new Error(data.error);
    const schemeId = makeProvenanceColorScheme(data.categories || []);
    entry.colorMode = "provenance";
    entry.provenanceSchemeId = schemeId;
    setCartoonColor(selector, schemeId);
    if (legendEl) {
      legendEl.innerHTML = renderProvenanceLegend();
      wireProvenanceLegendColorInputs(legendEl);
    }
  } catch (e) {
    entry.colorMode = "chain";
    entry.provenanceSchemeId = null;
    setCartoonColor(selector, chainColorScheme());
    if (legendEl) legendEl.innerHTML = `<span class="muted">Could not load .trb data: ${escapeHtml(e.message)}</span>`;
  }
}

// Color scheme for a click/selection highlight's own residue(s), as
// opposed to its neighbors: every element keeps its normal CPK color
// (reused from NGL's own built-in "element" scheme via getScheme(),
// rather than reimplementing that lookup table) except carbon, which
// becomes a vivid yellow - a clearly visible marker of "this is the
// selected residue" without needing a bigger radiusScale the way the
// neighbors-vs-selected distinction used to work.
let highlightColorSchemeId = null;
function highlightColorScheme() {
  if (highlightColorSchemeId) return highlightColorSchemeId;
  const elementScheme = NGL.ColormakerRegistry.getScheme({ scheme: "element" });
  highlightColorSchemeId = NGL.ColormakerRegistry.addScheme(function () {
    this.atomColor = function (atom) {
      if (atom.element === "C") return 0xffee00;
      return elementScheme.atomColor(atom);
    };
  });
  return highlightColorSchemeId;
}

// Click a residue to highlight it plus everything within
// RESIDUE_HIGHLIGHT_RADIUS Å of it, with a small overlay naming both.
// Hovering shows a lightweight floating tooltip for quick scanning
// without disturbing whatever's currently highlighted. Both use NGL's own
// picking API (stage.signals.clicked/hovered) and its selection language
// (structure.getAtomSetWithinSelection + getAtomSetWithinGroup for the
// "everything within N Å, expanded to whole residues" query) - no
// separate distance/highlight logic needed.
//
// The actual highlight-building logic (highlightResidue/highlightMultiple/
// clearHighlight) is exposed on the returned object rather than kept
// private, so something outside this function - the Results tab's
// sequence panel - can trigger the exact same highlight a 3D click would,
// instead of reimplementing it. setOnChange() lets that same panel keep
// its own residue selection in sync when a highlight is triggered the
// other way (a direct click in the 3D view).
function setupResidueInteraction(stage, component, container) {
  const structure = component.structure;

  const info = document.createElement("div");
  info.className = "residue-info hidden";
  container.appendChild(info);

  const tooltip = document.createElement("div");
  tooltip.className = "residue-tooltip hidden";
  stage.viewer.container.appendChild(tooltip);

  let highlightReprs = [];
  let onChange = null;

  function clearHighlight() {
    highlightReprs.forEach((r) => component.removeRepresentation(r));
    highlightReprs = [];
    info.classList.add("hidden");
    if (onChange) onChange([]);
  }

  // Single-residue highlight: the given residue plus everything within
  // RESIDUE_HIGHLIGHT_RADIUS Å of it, both as element-colored licorice
  // sticks (the residue itself drawn thicker than its neighbors instead
  // of in a different flat color, so element coloring stays meaningful
  // on both). Resname is looked up from the structure itself rather than
  // requiring the caller to know it, since a sequence-panel click only
  // has chain+resno to give us.
  function highlightResidue(chainname, resno) {
    highlightReprs.forEach((r) => component.removeRepresentation(r));
    highlightReprs = [];

    const residueSele = `${resno} and :${chainname}`;
    let resname = "?";
    let neighborCount = 0;
    try {
      structure.eachAtom((ap) => { resname = ap.resname; }, new NGL.Selection(residueSele));
      const nearAtoms = structure.getAtomSetWithinSelection(new NGL.Selection(residueSele), RESIDUE_HIGHLIGHT_RADIUS);
      const nearResidues = structure.getAtomSetWithinGroup(nearAtoms); // whole residues, not just the atoms caught by the radius

      const seen = new Set();
      structure.eachAtom((ap) => seen.add(`${ap.chainname}:${ap.resno}`), new NGL.Selection(nearResidues.toSeleString()));
      seen.delete(`${chainname}:${resno}`);
      neighborCount = seen.size;

      const neighborSele = `( ${nearResidues.toSeleString()} ) and not ( ${residueSele} )`;
      highlightReprs.push(component.addRepresentation("licorice", { sele: neighborSele, colorScheme: "element", radiusScale: 1.2 }));
    } catch (e) {
      // Highlight is best-effort - the info line below still works even
      // if the neighbor query fails for some edge-case selection.
    }
    // Same radiusScale as the neighbors - highlightColorScheme()'s vivid
    // yellow carbons are what marks this one out now, not extra size.
    highlightReprs.push(component.addRepresentation("licorice", { sele: residueSele, colorScheme: highlightColorScheme(), radiusScale: 1.2 }));

    info.innerHTML = `<b>${escapeHtml(resname || "?")} ${escapeHtml(String(resno))}</b> · chain ${escapeHtml(chainname || "?")}
      <span class="muted">— ${neighborCount} nearby residue${neighborCount === 1 ? "" : "s"} within ${RESIDUE_HIGHLIGHT_RADIUS} Å</span>
      <a class="clear-link">clear</a>`;
    info.classList.remove("hidden");
    info.querySelector(".clear-link").addEventListener("click", clearHighlight);

    if (onChange) onChange([{ chain: chainname, resno }]);
  }

  // Multi-residue highlight: the same idea as highlightResidue() above,
  // just seeded from every selected residue at once instead of one - the
  // whole selection's neighborhood (within RESIDUE_HIGHLIGHT_RADIUS Å) is
  // shown as normal-CPK licorice, and the selected residues themselves
  // use highlightColorScheme()'s vivid-yellow-carbon treatment so they
  // still read as "the selection" among their own neighbors.
  function highlightMultiple(residues) {
    highlightReprs.forEach((r) => component.removeRepresentation(r));
    highlightReprs = [];
    if (!residues.length) {
      info.classList.add("hidden");
      if (onChange) onChange([]);
      return;
    }

    const residueSele = residues.map((r) => `(${r.resno} and :${r.chain})`).join(" or ");
    let neighborCount = 0;
    try {
      const nearAtoms = structure.getAtomSetWithinSelection(new NGL.Selection(residueSele), RESIDUE_HIGHLIGHT_RADIUS);
      const nearResidues = structure.getAtomSetWithinGroup(nearAtoms); // whole residues, not just the atoms caught by the radius

      const seen = new Set();
      structure.eachAtom((ap) => seen.add(`${ap.chainname}:${ap.resno}`), new NGL.Selection(nearResidues.toSeleString()));
      residues.forEach((r) => seen.delete(`${r.chain}:${r.resno}`));
      neighborCount = seen.size;

      const neighborSele = `( ${nearResidues.toSeleString()} ) and not ( ${residueSele} )`;
      highlightReprs.push(component.addRepresentation("licorice", { sele: neighborSele, colorScheme: "element", radiusScale: 1.2 }));
    } catch (e) {
      // Best-effort, same reasoning as highlightResidue().
    }
    highlightReprs.push(component.addRepresentation("licorice", { sele: residueSele, colorScheme: highlightColorScheme(), radiusScale: 1.2 }));

    info.innerHTML = `<b>${residues.length} residues selected</b>
      <span class="muted">— ${neighborCount} nearby residue${neighborCount === 1 ? "" : "s"} within ${RESIDUE_HIGHLIGHT_RADIUS} Å</span>
      <a class="clear-link">clear</a>`;
    info.classList.remove("hidden");
    info.querySelector(".clear-link").addEventListener("click", clearHighlight);

    if (onChange) onChange(residues);
  }

  stage.signals.clicked.add((pickingProxy) => {
    if (!pickingProxy || !pickingProxy.atom) return;
    highlightResidue(pickingProxy.atom.chainname, pickingProxy.atom.resno);
  });

  stage.signals.hovered.add((pickingProxy) => {
    if (pickingProxy && pickingProxy.atom) {
      const atom = pickingProxy.atom;
      const cp = pickingProxy.canvasPosition;
      tooltip.textContent = `${atom.resname || "?"} ${atom.resno}`;
      tooltip.style.left = `${cp.x + 8}px`;
      tooltip.style.bottom = `${cp.y + 8}px`;
      tooltip.classList.remove("hidden");
    } else {
      tooltip.classList.add("hidden");
    }
  });

  return {
    highlightResidue,
    highlightMultiple,
    clearHighlight,
    setOnChange(cb) { onChange = cb; },
  };
}

// ---------------------------------------------------------------------
// Shared alignment viewer (Sequences tab + Results tab) - A / D2
//
// One shared horizontal scrollbar per group of same-length sequences,
// rather than one per sequence: all rows sit in a single overflow-x:auto
// container so scrolling any of them scrolls the whole block, and the
// label column stays pinned via position:sticky.
// ---------------------------------------------------------------------

function coloredSeq(seq, conservedChain) {
  return Array.from(seq).map((c, i) => {
    const color = AA_COLORS[c] || "#3a4560";
    const info = conservedChain && conservedChain[i];
    const highlight = info && info.isConserved && info.char === c;
    const cls = "aa" + (highlight ? " conserved" : "");
    // Conserved: solid background in the amino acid's own color with a
    // fixed dark, high-contrast text color (readable at 1-char width -
    // an outline/glow was getting lost against the neighboring residues).
    // Not conserved: the usual faint tint.
    const style = highlight ? `background:${color};color:#0b0f18` : `background:${color}1a;color:${color}`;
    return `<span class="${cls}" style="${style}">${c}</span>`;
  }).join("");
}

function renderRulerRow(chains) {
  let out = '<div class="align-ruler"><div class="align-label"></div><div class="align-seq">';
  let pos = 0;
  chains.forEach((chain, ci) => {
    for (let i = 0; i < chain.length; i++) {
      pos++;
      out += pos % 10 === 0 ? `<span class="tick">${pos}</span>` : `<span class="tick"></span>`;
    }
    if (ci < chains.length - 1) out += '<span class="chain-sep"> </span>';
  });
  out += "</div></div>";
  return out;
}

function renderAlignmentBlock(rows, options = {}) {
  // rows: [{label, chains: [seq, ...], swatch?: cssColor}]
  // options.conservation: per-chain array of {char, isConserved} (Sequences tab only)
  if (!rows.length) return `<p class="muted">Nothing to show.</p>`;
  const conservation = options.conservation || null;
  const ruler = renderRulerRow(rows[0].chains);
  const body = rows.map((r) => {
    const chainsHtml = r.chains
      .map((c, ci) => coloredSeq(c, conservation ? conservation[ci] : null))
      .join('<span class="chain-sep">/</span>');
    const swatch = r.swatch ? `<span class="swatch" style="background:${r.swatch}"></span>` : "";
    return `<div class="align-row"><div class="align-label" title="${escapeHtml(r.label)}">${swatch}${escapeHtml(r.label)}</div><div class="align-seq">${chainsHtml}</div></div>`;
  }).join("");
  return `<div class="align-scroll">${ruler}${body}</div>`;
}

// ---------------------------------------------------------------------
// Sequences tab
// ---------------------------------------------------------------------

// Per-position consensus across a backbone's samples: samples are
// guaranteed the same length/chain layout (same backbone), so residues
// line up index-for-index without needing a real alignment algorithm.
function computeConservation(samples, thresholdPct) {
  if (!samples.length) return null;
  const nChains = samples[0].chains.length;
  const perChain = [];
  for (let ci = 0; ci < nChains; ci++) {
    const len = samples[0].chains[ci] ? samples[0].chains[ci].length : 0;
    const positions = [];
    for (let pos = 0; pos < len; pos++) {
      const counts = {};
      let total = 0;
      samples.forEach((s) => {
        const chain = s.chains[ci];
        if (!chain || pos >= chain.length) return;
        const c = chain[pos];
        counts[c] = (counts[c] || 0) + 1;
        total++;
      });
      let bestChar = null, bestCount = 0;
      Object.entries(counts).forEach(([c, n]) => { if (n > bestCount) { bestCount = n; bestChar = c; } });
      const pct = total ? (bestCount / total) * 100 : 0;
      positions.push({ char: bestChar, isConserved: pct >= thresholdPct });
    }
    perChain.push(positions);
  }
  return perChain;
}

function renderSeqGroup(title, backbones) {
  if (!backbones.length) return "";
  const highlightOn = qs("#conservationToggle").checked;
  const thresholdPct = Math.max(1, Math.min(100, parseFloat(qs("#conservationThreshold").value) || 90));
  let html = `<h3>${escapeHtml(title)}</h3>`;
  backbones.forEach((bb) => {
    const rows = bb.samples.map((s, i) => {
      const meta = (s.score !== undefined)
        ? ` — score ${s.score} · global ${s.global_score} · recovery ${s.seq_recovery}`
        : "";
      return { label: `sample ${i + 1}${meta}`, chains: s.chains };
    });
    const conservation = highlightOn ? computeConservation(bb.samples, thresholdPct) : null;
    html += `<div class="backbone-block"><h4>${escapeHtml(bb.backbone)} <span class="muted" style="font-size:12px;font-weight:400">(${bb.num_samples} sequence${bb.num_samples === 1 ? "" : "s"})</span></h4>`;
    html += renderAlignmentBlock(rows, { conservation });
    html += `</div>`;
  });
  return html;
}

function renderSequencesContent() {
  const el = qs("#sequencesContent");
  const data = state.sequencesCache;
  if (!data) {
    el.innerHTML = `<p class="muted">Load a job first.</p>`;
    return;
  }
  if (!data.backbones.length && !data.monomers.length) {
    el.innerHTML = `<p class="muted">No sequences generated yet.</p>`;
    return;
  }
  el.innerHTML = renderSeqGroup("Designed sequences", data.backbones) + renderSeqGroup("Monomer sequences", data.monomers);
}

async function loadSequences() {
  const el = qs("#sequencesContent");
  try {
    state.sequencesCache = await apiGet("/api/sequences");
    renderSequencesContent();
  } catch (err) {
    el.innerHTML = `<p class="muted">${escapeHtml(err.message)}</p>`;
  }
}

// ---------------------------------------------------------------------
// Models tab
// ---------------------------------------------------------------------

function renderModelsCycleNote(status) {
  const el = qs("#modelsCycleNote");
  const cycle = status && status.cycle;
  if (!cycle || !cycle.is_multi_cycle) {
    el.classList.add("hidden");
    return;
  }
  el.classList.remove("hidden");
  el.innerHTML = `${cyclePillHtml(cycle)} Showing models from the current cycle only — earlier cycles' models are cleared from disk when the next cycle starts.`;
}

function modelKey(m) {
  return `${m.model}|${m.sequence_name}|${m.model_index}|${m.variant || ""}`;
}

function applyModelsSelectionHighlight() {
  qsa("#modelsList .list-item").forEach((el) => {
    el.classList.toggle("selected", el.dataset.key === state.selectedModelKey);
  });
}

function renderModelsList() {
  const listEl = qs("#modelsList");
  const filterVal = (qs("#modelsFilter").value || "").toLowerCase();
  const models = state.modelsCache;
  if (!models.length) {
    listEl.innerHTML = `<p class="muted">No models generated yet.</p>`;
    return;
  }
  rerenderPreservingScroll(listEl, () => {
    listEl.innerHTML = "";
    const visible = [];
    models.forEach((m) => {
      const label = `${m.model} / ${m.sequence_name}` + (m.model_index !== null && m.model_index !== undefined ? ` (model ${m.model_index})` : "");
      if (filterVal && !label.toLowerCase().includes(filterVal)) return;
      visible.push(m);
      const key = modelKey(m);
      const div = document.createElement("div");
      div.className = "list-item" + (state.selectedModelKey === key ? " selected" : "");
      div.dataset.key = key;
      const badge = m.variant ? `<span class="badge ${m.variant}">${m.variant}</span>` : "";
      div.innerHTML = `<span style="overflow:hidden;text-overflow:ellipsis">${escapeHtml(label)}</span>${badge}`;
      div.addEventListener("click", () => selectModel(m));
      listEl.appendChild(div);
    });
    const current = visible.find((m) => modelKey(m) === state.selectedModelKey) || visible[0];
    if (current) selectModel(current);
  });
}

async function loadModels() {
  const listEl = qs("#modelsList");
  try {
    state.modelsCache = await apiGet("/api/models");
    renderModelsList();
  } catch (err) {
    listEl.innerHTML = `<p class="muted">${escapeHtml(err.message)}</p>`;
  }
}

async function selectModel(m) {
  const key = modelKey(m);
  state.selectedModelKey = key;
  applyModelsSelectionHighlight();

  // Structure and confidence data are fetched/redrawn independently, and
  // only when the underlying file actually changed (e.g. AF3's .cif
  // upgrading to a final .pdb, or confidences appearing) - so an
  // already-loaded 3D view never gets reset just because a poll re-ran
  // this with the same selection.
  const structureKey = m.structure_path ? `${m.structure_path}|${m.structure_format}` : null;
  if (state.loadedModelStructureKey !== structureKey) {
    state.loadedModelStructureKey = structureKey;
    const label = qs("#modelViewerLabel");
    if (structureKey) {
      const url = `/api/model_pdb?log=${encodeURIComponent(state.logPath)}&path=${encodeURIComponent(m.structure_path)}`;
      try {
        const res = await fetch(url);
        if (!res.ok) throw new Error(await fetchErrorMessage(res));
        const text = await res.text();
        await renderMol("#modelViewer", text, m.structure_format || "pdb");
        label.textContent = m.structure_format === "cif" ? "Preview from .cif — final .pdb not written yet" : "";
        if (state.modelsColorMode === "provenance") {
          await applyProvenanceColoring("#modelViewer", m.trb_path, qs("#modelColorLegend"));
        }
        renderModelsSequencePanel();
      } catch (e) {
        qs("#modelViewer").innerHTML = `<p class="muted" style="padding:20px">Could not load structure: ${escapeHtml(e.message)}</p>`;
        label.textContent = "";
        renderModelsSequencePanel();
      }
    } else {
      qs("#modelViewer").innerHTML = `<p class="muted" style="padding:20px">No structure file found yet for this model.</p>`;
      label.textContent = "";
      renderModelsSequencePanel();
    }
  }

  const confKey = m.confidence_path || null;
  if (state.loadedModelConfKey !== confKey) {
    state.loadedModelConfKey = confKey;
    const metricsEl = qs("#modelMetrics");
    if (confKey) {
      const url = `/api/model_confidence?log=${encodeURIComponent(state.logPath)}&path=${encodeURIComponent(confKey)}`;
      const metrics = await fetch(url).then((r) => r.json());
      const rows = Object.entries(metrics).map(([k, v]) => `<tr><td>${escapeHtml(k)}</td><td>${escapeHtml(String(v))}</td></tr>`).join("");
      metricsEl.innerHTML = `<table>${rows}</table>`;
    } else {
      metricsEl.innerHTML = `<p class="muted">No confidence metrics found yet.</p>`;
    }
  }
}

// Interactive sequence panel for the Models tab - same read-straight-off-
// the-loaded-structure approach and click/drag/ctrl-click interaction as
// the Results tab's sequence panel (buildStructureSequence/
// renderSequenceResidues/setupSequencePanelInteraction are shared code).
function renderModelsSequencePanel() {
  const el = qs("#modelsSequenceContent");
  const entry = viewerRegistry["#modelViewer"];
  if (!entry || !entry.component) {
    el.innerHTML = "";
    return;
  }
  const chains = buildStructureSequence(entry.component.structure);
  el.innerHTML = renderSequenceResidues(chains);
  setupSequencePanelInteraction(el, "#modelViewer");
}

// ---------------------------------------------------------------------
// Results tab
// ---------------------------------------------------------------------

// The .trb sits next to the RFdiffusion pdb that final_output.csv's
// path_rfdiff column already points to - same basename, .trb extension,
// same convention as the Backbones tab.
function trbPathForResultRow(row) {
  return row.path_rfdiff ? row.path_rfdiff.replace(/\.pdb$/i, ".trb") : null;
}

// Interactive sequence panel (Results tab only): shows the currently
// selected row's sequence read straight off the already-loaded NGL
// structure, not a separately-fetched FASTA string - so a residue's
// position in this panel can never drift out of sync with the same
// residue in the 3D view next to it. Clicking a residue here calls the
// exact same highlightResidue() a direct 3D click would; selecting a
// range (drag) or several residues (ctrl/cmd+click) calls
// highlightMultiple() instead, which outlines the whole set.

function buildStructureSequence(structure) {
  const byChain = new Map();
  structure.eachResidue((rp) => {
    if (!rp.isProtein()) return; // skip waters/ligands/other hetero groups
    if (!byChain.has(rp.chainname)) byChain.set(rp.chainname, []);
    byChain.get(rp.chainname).push({ resno: rp.resno, code: rp.getResname1(), resname: rp.resname });
  });
  return Array.from(byChain, ([chain, residues]) => ({ chain, residues }));
}

// Residues are grouped into fixed-size chunks, each with its own starting
// residue number printed above it - a plain ruler-free run of letters
// left the residue number only reachable one-at-a-time by hovering, which
// is exactly the friction this panel exists to remove (its callers point
// users here specifically to "read off residue numbers" for contigs/
// hotspots). Chunking (rather than one number per residue) keeps it
// legible; chunks are inline-flex columns (number over residues) so the
// browser still wraps between them like plain inline content - no chunk
// itself splits across a line. qsa(".seq-residue", ...) still finds every
// residue span in document order regardless of this extra nesting, so
// setupSequencePanelInteraction()'s drag-select is unaffected.
const SEQ_CHUNK_SIZE = 10;

function renderSequenceResidues(chains) {
  return chains.map(({ chain, residues }) => {
    const chunks = [];
    for (let i = 0; i < residues.length; i += SEQ_CHUNK_SIZE) chunks.push(residues.slice(i, i + SEQ_CHUNK_SIZE));
    const chunkHtml = chunks.map((chunk) => {
      const spans = chunk.map((r) => {
        const color = AA_COLORS[r.code] || "#3a4560";
        return `<span class="seq-residue" data-chain="${escapeHtml(chain)}" data-resno="${r.resno}" title="${escapeHtml(r.resname)} ${r.resno}" style="background:${color}1a;color:${color}">${escapeHtml(r.code)}</span>`;
      }).join("");
      return `<span class="seq-chunk"><span class="seq-chunk-num">${chunk[0].resno}</span><span class="seq-chunk-residues">${spans}</span></span>`;
    }).join("");
    return `<div class="seq-chain"><span class="seq-chain-label">Chain ${escapeHtml(chain)} · ${residues.length} residues</span><div class="seq-residues">${chunkHtml}</div></div>`;
  }).join("");
}

// Plain click = select just this residue (mirrors a direct 3D click).
// Ctrl/Cmd+click = toggle this residue in the current multi-selection.
// Click-and-drag = select the dragged range. The 3D outline is only
// rebuilt on mouseup rather than on every mousemove - rebuilding an NGL
// surface representation per pixel of drag would be needlessly heavy -
// but the sequence panel's own selection still repaints live during the
// drag so the drag itself still feels immediate.
function setupSequencePanelInteraction(container, selector) {
  const residueEls = qsa(".seq-residue", container);
  let selected = new Set();
  let dragAnchor = null;
  let isDragging = false;

  function keyOf(el) { return `${el.dataset.chain}:${el.dataset.resno}`; }
  function paintSelection() {
    residueEls.forEach((el) => el.classList.toggle("selected", selected.has(keyOf(el))));
  }
  function applyHighlight() {
    const entry = viewerRegistry[selector];
    if (!entry || !entry.interaction) return;
    const list = residueEls
      .filter((el) => selected.has(keyOf(el)))
      .map((el) => ({ chain: el.dataset.chain, resno: parseInt(el.dataset.resno, 10) }));
    if (list.length === 0) entry.interaction.clearHighlight();
    else if (list.length === 1) entry.interaction.highlightResidue(list[0].chain, list[0].resno);
    else entry.interaction.highlightMultiple(list);
  }

  residueEls.forEach((el, idx) => {
    el.addEventListener("mousedown", (e) => {
      if (e.ctrlKey || e.metaKey) {
        const k = keyOf(el);
        if (selected.has(k)) selected.delete(k); else selected.add(k);
        paintSelection();
        applyHighlight();
        return; // a ctrl/cmd-click toggles one residue; it doesn't start a drag
      }
      selected = new Set([keyOf(el)]);
      dragAnchor = idx;
      isDragging = true;
      paintSelection();
    });
    el.addEventListener("mouseenter", () => {
      if (!isDragging || dragAnchor === null) return;
      const lo = Math.min(dragAnchor, idx);
      const hi = Math.max(dragAnchor, idx);
      selected = new Set(residueEls.slice(lo, hi + 1).map(keyOf));
      paintSelection();
    });
  });

  // A drag can end with the mouse anywhere on the page, not just over a
  // residue span, so this listener lives on the document rather than the
  // panel itself - and gets replaced (not stacked) on every re-render, by
  // stashing the handler reference on the container, which persists
  // across re-renders even though its contents get replaced each time.
  if (container._seqMouseupHandler) document.removeEventListener("mouseup", container._seqMouseupHandler);
  const mouseupHandler = () => {
    if (!isDragging) return;
    isDragging = false;
    dragAnchor = null;
    applyHighlight();
  };
  container._seqMouseupHandler = mouseupHandler;
  document.addEventListener("mouseup", mouseupHandler);

  // Keeps the panel in sync when a highlight is triggered the other way -
  // clicking directly in the 3D view - so that doesn't leave the sequence
  // panel still showing an earlier, now-stale selection.
  const entry = viewerRegistry[selector];
  if (entry && entry.interaction) {
    entry.interaction.setOnChange((residues) => {
      selected = new Set(residues.map((r) => `${r.chain}:${r.resno}`));
      paintSelection();
    });
  }
}

// ---------------------------------------------------------------------
// Results view - a factory rather than one hardwired set of functions,
// so the per-job Results tab and the All jobs results tab share the
// exact same filtering/coloring/CSV/alignment/export logic instead of
// two copies that could quietly drift apart. `cfg` supplies which DOM
// ids and viewer selector a given instance draws into, how to fetch its
// data (one job vs every tracked job merged), and how to resolve which
// job's API a given row belongs to.
// ---------------------------------------------------------------------

// A production job's result set can be tens of thousands of rows. Both
// the model list and the metrics table render one row's worth of DOM per
// row currently in view - at that scale, building it for the whole set
// at once (as this used to do) freezes the tab for minutes. Bounding how
// many rows are actually in the DOM at a time keeps each render fast
// regardless of how big the underlying job is.
const RESULTS_PAGE_SIZE = 200;
// Per sequence-length group in the alignment view - see renderAlignment().
const ALIGN_GROUP_ROW_CAP = 100;

// Small "‹ Prev  1–200 of 15,000  Next ›" footer, shared by the results
// list and the metrics table. Renders nothing when everything already
// fits on one page, so it's a no-op for the (very common) small-job case.
function renderPagerControls(el, page, pageSize, totalCount, onChange) {
  if (totalCount <= pageSize) { el.innerHTML = ""; return; }
  const totalPages = Math.max(1, Math.ceil(totalCount / pageSize));
  const start = page * pageSize + 1;
  const end = Math.min(totalCount, start + pageSize - 1);
  el.innerHTML = `
    <button class="pager-btn" ${page === 0 ? "disabled" : ""} data-dir="-1">‹ Prev</button>
    <span class="pager-info">${start.toLocaleString()}–${end.toLocaleString()} of ${totalCount.toLocaleString()} · page ${page + 1}/${totalPages}</span>
    <button class="pager-btn" ${page >= totalPages - 1 ? "disabled" : ""} data-dir="1">Next ›</button>
  `;
  qsa(".pager-btn", el).forEach((btn) => {
    if (btn.disabled) return;
    btn.addEventListener("click", () => onChange(page + parseInt(btn.dataset.dir, 10)));
  });
}

function createResultsView(cfg) {
  let results = {
    columns: [], rows: [], numericColumns: [], displayColumns: [], idCol: null,
    filters: {}, colorField: "", colorDirection: "higher", selectedRowId: null,
    // sortField/sortDir (1 = ascending, -1 = descending) is the ONE sort
    // state shared by the list, the metrics table, and the "Sort list by"
    // toolbar control - set from any of the three (the dropdown, or
    // clicking a table header) and every one of them reflects it, rather
    // than the table having its own independent click-to-sort state the
    // list and dropdown never knew about (which is how this used to work).
    sortField: "", sortDir: 1,
    listPage: 0, tablePage: 0,
  };
  let loadedPdbKey = null;
  let colorMode = "chain"; // "chain" | "provenance"
  // Per sequence-length group in the alignment view (keyed by the group's
  // length) - which ALIGN_GROUP_ROW_CAP-sized page is currently shown.
  let expandedAlignGroups = new Map();

  function ingestCsv(csv) {
    const { columns, rows, rowMeta } = csv;
    const objRows = rows.map((r, i) => {
      const obj = { __rowId: i };
      columns.forEach((c, ci) => { obj[c] = r[ci]; });
      if (rowMeta && rowMeta[i]) { obj.__jobLogPath = rowMeta[i].log; obj.__origRowId = rowMeta[i].origRowId; }
      return obj;
    });

    const numericColumns = columns.filter((c) => {
      const vals = objRows.map((r) => r[c]).filter((v) => v !== "" && v !== undefined && v !== null);
      if (!vals.length) return false;
      return vals.every((v) => !isNaN(parseFloat(v)) && isFinite(v));
    });

    const idCol = columns.includes("id") ? "id" : null;
    const jobCol = columns.includes("job") ? "job" : null;
    const pathCols = columns.filter((c) => PATH_COLUMNS.includes(c));
    const rest = columns.filter((c) => c !== idCol && c !== jobCol && !pathCols.includes(c));
    const displayColumns = [...(idCol ? [idCol] : []), ...(jobCol ? [jobCol] : []), ...rest, ...pathCols];

    const prev = results;
    results = {
      columns, rows: objRows, numericColumns, displayColumns, idCol,
      filters: prev.filters || {}, colorField: prev.colorField || "",
      colorDirection: prev.colorDirection || "higher", selectedRowId: prev.selectedRowId,
      sortField: prev.sortField || "", sortDir: prev.sortDir || 1,
      // Reset to page 1 on every fresh load - the row set may have changed
      // size/order entirely (new poll tick, new filters upstream), so
      // whatever page the user was on no longer means anything reliable.
      listPage: 0, tablePage: 0,
    };
    expandedAlignGroups = new Map();
    // drop filters/colorField/sortField that no longer refer to a real column
    Object.keys(results.filters).forEach((c) => { if (!numericColumns.includes(c)) delete results.filters[c]; });
    if (results.colorField && !numericColumns.includes(results.colorField)) results.colorField = "";
    if (results.sortField && !displayColumns.includes(results.sortField)) results.sortField = "";
  }

  function getFilteredRows() {
    const { rows, filters } = results;
    const cols = Object.keys(filters);
    if (!cols.length) return rows;
    return rows.filter((r) => {
      for (const col of cols) {
        const f = filters[col];
        if (f.min == null && f.max == null) continue;
        const val = parseFloat(r[col]);
        if (isNaN(val)) return false;
        if (f.min != null && val < f.min) return false;
        if (f.max != null && val > f.max) return false;
      }
      return true;
    });
  }

  // Shared by the list, the metrics table, and (indirectly, since it
  // reads results.sortField/sortDir) the "Sort list by" toolbar control -
  // one sort order, computed once per render rather than the table having
  // its own separate copy of this logic the list never saw. Numeric
  // columns compare numerically; anything else (including a numeric
  // column with the odd blank/non-numeric cell) falls back to a plain
  // string compare, same as the table's click-to-sort always did.
  function getSortedRows(rows) {
    const { sortField, sortDir } = results;
    if (!sortField) return rows;
    return [...rows].sort((a, b) => {
      const av = a[sortField], bv = b[sortField];
      const an = parseFloat(av), bn = parseFloat(bv);
      const cmp = !isNaN(an) && !isNaN(bn) ? an - bn : String(av ?? "").localeCompare(String(bv ?? ""));
      return cmp * sortDir;
    });
  }

  function computeColorScale(rows) {
    const { colorField, colorDirection } = results;
    if (!colorField) return null;
    const vals = rows.map((r) => parseFloat(r[colorField])).filter((v) => !isNaN(v));
    if (!vals.length) return null;
    return { min: Math.min(...vals), max: Math.max(...vals), field: colorField, direction: colorDirection };
  }

  function colorForRow(row, scale) {
    if (!scale) return null;
    const val = parseFloat(row[scale.field]);
    if (isNaN(val)) return null;
    let t = scale.max === scale.min ? 1 : (val - scale.min) / (scale.max - scale.min);
    if (scale.direction === "lower") t = 1 - t;
    t = Math.max(0, Math.min(1, t));
    const hue = t * 145; // 0 = red (worst) .. 145 = green (best)
    // This background is deliberately dark regardless of theme (so the hue
    // stays legible), so it needs a fixed light text color paired with it -
    // the ambient theme text turns dark in light mode, which would be
    // unreadable against it.
    return { bg: `hsl(${hue}, 42%, 21%)`, accent: `hsl(${hue}, 70%, 52%)`, text: "#eef1f8" };
  }

  async function load() {
    const listEl = qs(cfg.ids.list);
    try {
      const data = await cfg.fetchData();
      if (!data.csv_exists) {
        results.rows = [];
        if (!data.final_pdbs.length) {
          listEl.innerHTML = `<p class="muted">${cfg.emptyMessage}</p>`;
        } else {
          // Fallback (no csv yet, e.g. mid-scoring): plain filename list.
          listEl.innerHTML = "";
          data.final_pdbs.forEach((entry) => {
            const label = entry.jobLabel ? `${entry.jobLabel} / ${entry.name}` : entry.name;
            const div = document.createElement("div");
            div.className = "list-item";
            div.innerHTML = `<span style="overflow:hidden;text-overflow:ellipsis">${escapeHtml(label)}</span>`;
            div.addEventListener("click", async () => {
              qsa(".list-item", listEl).forEach((i) => i.classList.remove("selected"));
              div.classList.add("selected");
              const url = `/api/final_pdb?log=${encodeURIComponent(entry.log)}&name=${encodeURIComponent(entry.name)}`;
              const pdbText = await fetch(url).then((r) => r.text());
              await renderMol(cfg.viewerSelector, pdbText);
            });
            listEl.appendChild(div);
          });
        }
        qs(cfg.ids.alignContent).innerHTML = `<p class="muted">Not available yet.</p>`;
        qs(cfg.ids.csvWrap).innerHTML = `<p class="muted">${cfg.emptyCsvMessage}</p>`;
        qs(cfg.ids.count).textContent = "";
        populateColorBySelect();
        populateSortBySelect();
        renderFiltersList();
        return;
      }
      ingestCsv(data);
      populateColorBySelect();
      populateSortBySelect();
      renderFiltersList();
      updateFiltersBadge();
      renderAll();
    } catch (err) {
      listEl.innerHTML = `<p class="muted">${escapeHtml(err.message)}</p>`;
    }
  }

  function renderAll() {
    const total = results.rows.length;
    const shown = getFilteredRows().length;
    qs(cfg.ids.count).textContent = total ? `Showing ${shown} / ${total} model${total === 1 ? "" : "s"}` : "";
    renderList();
    renderAlignment();
    renderCsvTableFn();
  }

  // A production job's result set can run into the tens of thousands of
  // rows - building one DOM element per row (times however many per-row
  // elements a given view needs, e.g. one <span> per residue in the
  // alignment view) for the *entire* set at once is what made the list,
  // table and alignment views freeze the tab for minutes on a large job.
  // Paging bounds every one of them to RESULTS_PAGE_SIZE rows actually in
  // the DOM at a time, same idea a normal search results page uses.
  function pageOf(rows, page, pageSize) {
    const totalPages = Math.max(1, Math.ceil(rows.length / pageSize));
    const clamped = Math.max(0, Math.min(page, totalPages - 1));
    const start = clamped * pageSize;
    return { pageRows: rows.slice(start, start + pageSize), clamped, totalPages, start };
  }

  function renderList() {
    const listEl = qs(cfg.ids.list);
    const pagerEl = qs(cfg.ids.listPager);
    if (!results.rows.length) return; // handled by load()'s fallback path
    const rows = getSortedRows(getFilteredRows());
    if (!rows.length) {
      listEl.innerHTML = `<p class="muted">No models match the current filters.</p>`;
      if (pagerEl) pagerEl.innerHTML = "";
      return;
    }
    const { pageRows, clamped } = pageOf(rows, results.listPage, RESULTS_PAGE_SIZE);
    results.listPage = clamped;
    if (pagerEl) {
      renderPagerControls(pagerEl, clamped, RESULTS_PAGE_SIZE, rows.length, (newPage) => {
        results.listPage = newPage;
        renderList();
      });
    }
    rerenderPreservingScroll(listEl, () => {
      // The color scale (min/max for the red→green gradient) is computed
      // from the *whole* filtered set, not just this page - otherwise the
      // same score would render a different color depending on which page
      // it happened to land on.
      const scale = computeColorScale(rows);
      listEl.innerHTML = "";
      pageRows.forEach((row) => {
        const label = row[results.idCol] || basename(row.model_path) || `row ${row.__rowId}`;
        const div = document.createElement("div");
        div.className = "list-item" + (results.selectedRowId === row.__rowId ? " selected" : "");
        const color = colorForRow(row, scale);
        if (color) {
          div.style.background = color.bg;
          div.style.borderLeftColor = color.accent;
          div.style.color = color.text;
        }
        const jobTag = cfg.isMulti && row.job ? ` <span class="badge pending">${escapeHtml(row.job)}</span>` : "";
        div.innerHTML = `<span style="overflow:hidden;text-overflow:ellipsis">${escapeHtml(String(label))}${jobTag}</span>`;
        div.addEventListener("click", () => selectRow(row));
        listEl.appendChild(div);
      });
      const current = rows.find((r) => r.__rowId === results.selectedRowId) || rows[0];
      if (current) selectRow(current);
    });
  }

  function renderSequencePanel() {
    const el = qs(cfg.ids.seqContent);
    const entry = viewerRegistry[cfg.viewerSelector];
    if (!entry || !entry.component) {
      el.innerHTML = "";
      return;
    }
    const chains = buildStructureSequence(entry.component.structure);
    el.innerHTML = renderSequenceResidues(chains);
    setupSequencePanelInteraction(el, cfg.viewerSelector);
  }

  async function selectRow(row) {
    results.selectedRowId = row.__rowId;
    renderListSelectionOnly();
    renderTableSelectionOnly();
    const pdbName = basename(row.model_path);
    const rowLog = cfg.getLog(row);
    // Same structure already showing - skip refetching so the viewer's
    // camera (zoom/rotation) isn't reset on every poll tick. Keyed by job
    // too, not just the pdb name - two different jobs can both have a
    // model literally named "model_0.pdb".
    const pdbKey = `${rowLog}|${pdbName}`;
    if (loadedPdbKey === pdbKey) return;
    loadedPdbKey = pdbKey;
    if (!pdbName) return;
    const url = `/api/final_pdb?log=${encodeURIComponent(rowLog)}&name=${encodeURIComponent(pdbName)}`;
    try {
      const res = await fetch(url);
      if (!res.ok) throw new Error(await fetchErrorMessage(res));
      const pdbText = await res.text();
      const stage = await renderMol(cfg.viewerSelector, pdbText);
      if (colorMode === "provenance") {
        await applyProvenanceColoring(cfg.viewerSelector, trbPathForResultRow(row), qs(cfg.ids.viewerLegend), rowLog);
      }
      // renderMol() returns null (without throwing) both on a genuine load
      // failure and when a newer selection raced past this one, so only
      // rebuild the sequence panel on an actual successful load -
      // otherwise it could rebuild from a *previous* row's stale
      // structure still sitting in the registry.
      if (stage) renderSequencePanel(); else qs(cfg.ids.seqContent).innerHTML = "";
    } catch (e) {
      qs(cfg.viewerSelector).innerHTML = `<p class="muted" style="padding:20px">Could not load structure for this row: ${escapeHtml(e.message)}</p>`;
      qs(cfg.ids.seqContent).innerHTML = "";
    }
  }

  function renderListSelectionOnly() {
    // Cheap re-render (only this page's rows are in the DOM) so the
    // highlighted row + viewer stay in sync. If the selected row isn't on
    // the page currently shown (e.g. a poll re-selected a row that's now
    // on a different page), nothing in the visible list gets highlighted -
    // simpler and safer than silently jumping the user to another page.
    const { pageRows } = pageOf(getSortedRows(getFilteredRows()), results.listPage, RESULTS_PAGE_SIZE);
    const listEl = qs(cfg.ids.list);
    qsa(".list-item", listEl).forEach((el, i) => {
      const row = pageRows[i];
      if (!row) return;
      el.classList.toggle("selected", row.__rowId === results.selectedRowId);
    });
  }

  // Same idea as renderListSelectionOnly() - a cheap toggle over whatever
  // rows the metrics table currently has in the DOM, rather than
  // rebuilding the whole table (which would also reset its scroll
  // position and re-run every color-by/id-column computation) just to
  // move a selection highlight. Also does nothing if the selected row
  // isn't on the table's own currently-shown page - the table has its own
  // independent page from the list (results.tablePage vs listPage), so
  // "selected" not being visible right now is expected, not a bug.
  function renderTableSelectionOnly() {
    const wrap = qs(cfg.ids.csvWrap);
    qsa("tbody tr", wrap).forEach((tr) => {
      const rid = parseInt(tr.dataset.rowId, 10);
      tr.classList.toggle("selected", rid === results.selectedRowId);
    });
  }

  function renderAlignment() {
    const el = qs(cfg.ids.alignContent);
    if (!results.rows.length) {
      el.innerHTML = `<p class="muted">Not available yet.</p>`;
      return;
    }
    const rows = getFilteredRows();
    if (!rows.length) {
      el.innerHTML = `<p class="muted">No models match the current filters.</p>`;
      return;
    }
    if (!results.columns.includes("sequence")) {
      el.innerHTML = `<p class="muted">final_output.csv has no "sequence" column.</p>`;
      return;
    }
    const scale = computeColorScale(rows);

    // Group by sequence length: unrelated designs won't be the same length,
    // but this keeps each group's residues lined up in one shared scrollbar.
    const groups = new Map();
    rows.forEach((row) => {
      const seq = row.sequence || "";
      const key = seq.length;
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key).push(row);
    });

    // Each group already only shows 5 rows at a time (CSS-scrolled), but
    // that scroll cap alone doesn't stop the browser from having to build
    // one <span> per residue for *every* row in an oversized group up
    // front - a group with thousands of same-length designs was exactly
    // what made this view freeze the tab. Each group gets its own pager
    // (same ALIGN_GROUP_ROW_CAP-sized pages as the list/table use) rather
    // than an unbounded "show everything" escape hatch - a "show all"
    // button on a 15,000-row group would just recreate the same freeze
    // for whoever clicked it.
    const html = [...groups.keys()].sort((a, b) => b - a).map((len) => {
      const groupRows = groups.get(len);
      const groupPage = expandedAlignGroups.get(len) || 0;
      const { pageRows: shown, clamped } = pageOf(groupRows, groupPage, ALIGN_GROUP_ROW_CAP);
      expandedAlignGroups.set(len, clamped);
      const alignRows = shown.map((row) => {
        const color = colorForRow(row, scale);
        let label = String(row[results.idCol] || basename(row.model_path) || `row ${row.__rowId}`);
        if (cfg.isMulti && row.job) label = `${row.job} / ${label}`;
        return { label, chains: String(row.sequence || "").split(":"), swatch: color ? color.accent : null };
      });
      return `<div class="backbone-block"><h4>${len} residues <span class="muted" style="font-size:12px;font-weight:400">(${groupRows.length} model${groupRows.length === 1 ? "" : "s"})</span></h4>${renderAlignmentBlock(alignRows)}<div class="pager-bar align-group-pager" data-len="${len}"></div></div>`;
    }).join("");
    el.innerHTML = html;
    qsa(".align-group-pager", el).forEach((pagerEl) => {
      const len = parseInt(pagerEl.dataset.len, 10);
      const groupRows = groups.get(len);
      renderPagerControls(pagerEl, expandedAlignGroups.get(len) || 0, ALIGN_GROUP_ROW_CAP, groupRows.length, (newPage) => {
        expandedAlignGroups.set(len, newPage);
        renderAlignment();
      });
    });
  }

  function renderCsvTableFn() {
    const wrap = qs(cfg.ids.csvWrap);
    const pagerEl = qs(cfg.ids.tablePager);
    if (!results.rows.length) {
      wrap.innerHTML = `<p class="muted">Not available yet.</p>`;
      if (pagerEl) pagerEl.innerHTML = "";
      return;
    }
    const rows = getFilteredRows();
    if (!rows.length) {
      wrap.innerHTML = `<p class="muted">No rows match the current filters.</p>`;
      if (pagerEl) pagerEl.innerHTML = "";
      return;
    }
    const { displayColumns, idCol } = results;
    const scale = computeColorScale(rows);
    // Same sort (results.sortField/sortDir) the list and the "Sort list
    // by" toolbar control use - see getSortedRows(). Clicking a column
    // header below updates that shared state too, rather than the table
    // keeping its own separate sort no other view agreed with.
    const sortedRows = getSortedRows(rows);

    const { pageRows, clamped } = pageOf(sortedRows, results.tablePage, RESULTS_PAGE_SIZE);
    results.tablePage = clamped;
    if (pagerEl) {
      renderPagerControls(pagerEl, clamped, RESULTS_PAGE_SIZE, sortedRows.length, (newPage) => {
        results.tablePage = newPage;
        renderCsvTableFn();
      });
    }

    const thead = displayColumns.map((c) => {
      const cls = c === idCol ? ' class="id-col"' : (PATH_COLUMNS.includes(c) ? ' class="path-col"' : "");
      const arrow = results.sortField === c ? (results.sortDir === 1 ? " ▲" : " ▼") : "";
      return `<th data-col="${escapeHtml(c)}"${cls}>${escapeHtml(c)}${arrow}</th>`;
    }).join("");

    const tbody = pageRows.map((r, idx) => {
      const color = colorForRow(r, scale);
      const rowBg = color ? color.bg : (idx % 2 === 1 ? "var(--panel-2)" : "transparent");
      const idCellBg = color ? color.bg : (idx % 2 === 1 ? "var(--panel-2)" : "var(--panel)");
      // The row's own `color` (inherited by every <td>) is enough to fix
      // text against the tint everywhere *except* .path-col, which sets its
      // own `color: var(--muted)` in CSS - an inline style is needed there
      // to actually win over that class rule.
      const rowStyle = color ? `background:${rowBg};color:${color.text}` : `background:${rowBg}`;
      const tds = displayColumns.map((c) => {
        const isId = c === idCol;
        const isPath = PATH_COLUMNS.includes(c);
        const cls = isId ? ' class="id-col"' : (isPath ? ' class="path-col"' : "");
        const style = isId ? ` style="background:${idCellBg}"` : (isPath && color ? ` style="color:${color.text}"` : "");
        const val = r[c] ?? "";
        return `<td${cls}${style} title="${escapeHtml(String(val))}">${escapeHtml(String(val))}</td>`;
      }).join("");
      // .selected uses borders (not background/outline - background is
      // already spoken for by color-by, and outline on a <tr> renders
      // inconsistently across browsers with border-collapse) so it stays
      // visible layered on top of any color-by tint - see .data-table
      // tbody tr.selected in style.css.
      const selectedCls = r.__rowId === results.selectedRowId ? " selected" : "";
      return `<tr class="${selectedCls}" style="${rowStyle}" data-row-id="${r.__rowId}">${tds}</tr>`;
    }).join("");

    const wrapTop = wrap.scrollTop;
    wrap.innerHTML = `<table class="data-table"><thead><tr>${thead}</tr></thead><tbody>${tbody}</tbody></table>`;
    wrap.scrollTop = wrapTop;

    qsa("th", wrap).forEach((th) => {
      th.addEventListener("click", () => {
        const c = th.dataset.col;
        // Clicking the column already being sorted flips direction;
        // clicking a different one starts it fresh at ascending.
        results.sortDir = results.sortField === c ? -results.sortDir : 1;
        results.sortField = c;
        results.listPage = 0; results.tablePage = 0; // a new sort order means "page 3" means something different now
        syncSortControls();
        renderAll(); // the list needs to re-sort too now, not just this table
      });
    });
    qsa("tbody tr", wrap).forEach((tr) => {
      tr.addEventListener("click", () => {
        const rid = parseInt(tr.dataset.rowId, 10);
        const row = results.rows.find((r) => r.__rowId === rid);
        if (row) selectRow(row);
      });
    });
  }

  // --- Color-by controls (D4) ---

  function populateColorBySelect() {
    const sel = qs(cfg.ids.colorByField);
    const dirSel = qs(cfg.ids.colorByDirection);
    const prev = results.colorField;
    sel.innerHTML = `<option value="">None</option>` + results.numericColumns
      .map((c) => `<option value="${escapeHtml(c)}">${escapeHtml(c)}</option>`).join("");
    sel.value = results.numericColumns.includes(prev) ? prev : "";
    results.colorField = sel.value;
    dirSel.value = results.colorDirection;
  }

  // --- Sort-by controls ---
  // Every displayColumn (not just numericColumns, unlike color-by - a
  // plain string compare is a perfectly fine sort even if it wouldn't
  // make sense as a color gradient) is offered here, and this one state
  // (results.sortField/sortDir) drives the list, the metrics table, AND
  // this dropdown - whichever of the three last changed it, the other two
  // pick it up next render (see getSortedRows(), and the table's own
  // header-click handler, which calls syncSortControls() too).

  function populateSortBySelect() {
    const sel = qs(cfg.ids.sortByField);
    const dirSel = qs(cfg.ids.sortByDirection);
    const prev = results.sortField;
    sel.innerHTML = `<option value="">None</option>` + results.displayColumns
      .map((c) => `<option value="${escapeHtml(c)}">${escapeHtml(c)}</option>`).join("");
    sel.value = results.displayColumns.includes(prev) ? prev : "";
    results.sortField = sel.value;
    dirSel.value = String(results.sortDir);
  }

  // Cheap sync (just .value, not a full rebuild) for when sortField/sortDir
  // changed from somewhere other than these controls themselves (a table
  // header click) - rebuilding the <select>'s options here would be
  // needless work and would risk fighting an in-progress click on it.
  function syncSortControls() {
    qs(cfg.ids.sortByField).value = results.sortField;
    qs(cfg.ids.sortByDirection).value = String(results.sortDir);
  }

  // --- Filters panel (D3) ---

  function initFiltersPanel() {
    const panel = qs(cfg.ids.filtersPanel);
    qs(cfg.ids.filtersBtn).addEventListener("click", (e) => {
      e.stopPropagation();
      panel.classList.toggle("hidden");
    });
    panel.addEventListener("click", (e) => e.stopPropagation());
    document.addEventListener("click", () => panel.classList.add("hidden"));
    qs(cfg.ids.clearFiltersBtn).addEventListener("click", () => {
      results.filters = {};
      results.listPage = 0; results.tablePage = 0; // the filtered set is changing
      renderFiltersList();
      updateFiltersBadge();
      renderAll();
    });
  }

  function renderFiltersList() {
    const el = qs(cfg.ids.filtersList);
    const cols = results.numericColumns || [];
    if (!cols.length) {
      el.innerHTML = `<p class="muted" style="margin:0">No numeric columns found yet.</p>`;
      return;
    }
    el.innerHTML = cols.map((c) => {
      const f = results.filters[c] || {};
      const active = f.min != null || f.max != null;
      return `
        <div class="filter-row${active ? " filter-active" : ""}">
          <span class="filter-col-name" title="${escapeHtml(c)}">${escapeHtml(c)}</span>
          <input type="number" class="filter-min" data-col="${escapeHtml(c)}" placeholder="min" value="${f.min ?? ""}">
          <span class="filter-dash">–</span>
          <input type="number" class="filter-max" data-col="${escapeHtml(c)}" placeholder="max" value="${f.max ?? ""}">
        </div>`;
    }).join("");
    qsa(".filter-min, .filter-max", el).forEach((inp) => {
      inp.addEventListener("input", () => {
        const col = inp.dataset.col;
        const isMin = inp.classList.contains("filter-min");
        const val = inp.value === "" ? null : parseFloat(inp.value);
        const f = results.filters[col] || {};
        if (isMin) f.min = val; else f.max = val;
        if (f.min == null && f.max == null) delete results.filters[col];
        else results.filters[col] = f;
        inp.closest(".filter-row").classList.toggle("filter-active", f.min != null || f.max != null);
        updateFiltersBadge();
        results.listPage = 0; results.tablePage = 0; // the filtered set is changing
        renderAll();
      });
    });
  }

  function updateFiltersBadge() {
    const n = Object.keys(results.filters).length;
    const badge = qs(cfg.ids.filtersBadge);
    badge.textContent = String(n);
    badge.classList.toggle("hidden", n === 0);
  }

  // --- "Color structure by" (chain / RFdiffusion provenance) ---

  function initModeDropdown() {
    qs(cfg.ids.structColorMode).addEventListener("change", async (e) => {
      colorMode = e.target.value;
      const legendEl = qs(cfg.ids.viewerLegend);
      const paletteWrap = qs(cfg.ids.chainPaletteWrap);
      if (colorMode === "provenance") {
        paletteWrap.classList.add("hidden");
        const row = results.rows.find((r) => r.__rowId === results.selectedRowId);
        await applyProvenanceColoring(cfg.viewerSelector, row ? trbPathForResultRow(row) : null, legendEl, row ? cfg.getLog(row) : undefined);
      } else {
        paletteWrap.classList.remove("hidden");
        setCartoonColor(cfg.viewerSelector, chainColorScheme());
        const entry = viewerRegistry[cfg.viewerSelector];
        if (entry) { entry.colorMode = "chain"; entry.provenanceSchemeId = null; }
        refreshChainLegendFor(cfg.viewerSelector);
      }
    });
  }

  // --- Export filtered (F) ---

  // Rows are tagged with which job they came from (getLog) - for a
  // multi-job export that means grouping by job and sending each job its
  // own row indices, since a merged row's __rowId is an index into the
  // *combined* list, not that job's own final_output.csv.
  function groupRowsByJob(rows) {
    const byJob = new Map();
    rows.forEach((r) => {
      const log = cfg.getLog(r);
      if (!byJob.has(log)) byJob.set(log, []);
      byJob.get(log).push(r.__origRowId);
    });
    return Array.from(byJob, ([log, row_ids]) => ({ log, row_ids }));
  }

  async function doExport() {
    if (!results.rows.length) return;
    const rows = getFilteredRows();
    if (!rows.length) {
      alert("No rows match the current filters — nothing to export.");
      return;
    }
    const btn = qs(cfg.ids.exportBtn);
    const originalText = btn.textContent;
    btn.textContent = "Exporting…";
    btn.disabled = true;
    try {
      const body = cfg.isMulti
        ? { jobs: groupRowsByJob(rows) }
        : { log: state.logPath, row_ids: rows.map((r) => r.__rowId) };
      const res = await fetch(cfg.exportUrl, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
      });
      if (!res.ok) {
        let msg = res.statusText;
        try { const j = await res.json(); msg = j.error || j.message || msg; } catch (e) {}
        throw new Error(msg);
      }
      const blob = await res.blob();
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = cfg.exportFilename;
      document.body.appendChild(a);
      a.click();
      a.remove();
      URL.revokeObjectURL(url);
    } catch (err) {
      alert("Export failed: " + err.message);
    } finally {
      btn.textContent = originalText;
      btn.disabled = false;
    }
  }

  function initExportButton() {
    qs(cfg.ids.exportBtn).addEventListener("click", doExport);
  }

  function initColorByControls() {
    qs(cfg.ids.colorByField).addEventListener("change", (e) => { results.colorField = e.target.value; renderAll(); });
    qs(cfg.ids.colorByDirection).addEventListener("change", (e) => { results.colorDirection = e.target.value; renderAll(); });
  }

  function initSortByControls() {
    qs(cfg.ids.sortByField).addEventListener("change", (e) => {
      results.sortField = e.target.value;
      results.listPage = 0; results.tablePage = 0; // a new sort order means "page 3" means something different now
      renderAll();
    });
    qs(cfg.ids.sortByDirection).addEventListener("change", (e) => {
      results.sortDir = parseInt(e.target.value, 10);
      results.listPage = 0; results.tablePage = 0;
      renderAll();
    });
  }

  function initAll() {
    initFiltersPanel();
    initModeDropdown();
    initColorByControls();
    initSortByControls();
    initExportButton();
  }

  function reset() {
    results = {
      columns: [], rows: [], numericColumns: [], displayColumns: [], idCol: null,
      filters: {}, colorField: "", colorDirection: "higher", selectedRowId: null,
      sortField: "", sortDir: 1,
      listPage: 0, tablePage: 0,
    };
    expandedAlignGroups = new Map();
    loadedPdbKey = null;
    colorMode = "chain";
  }

  return { load, renderAll, selectRow, initAll, reset };
}

// Single-job fetch: what the per-job Results tab's createResultsView()
// instance uses to load its data (its own job's /api/results +
// /api/final_csv, no merging needed).
async function fetchSingleJobResults() {
  const res = await apiGet("/api/results");
  if (!res.csv_exists) {
    return { csv_exists: false, final_pdbs: (res.final_pdbs || []).map((name) => ({ name, log: state.logPath })) };
  }
  const csv = await apiGet("/api/final_csv");
  return { csv_exists: true, columns: csv.columns, rows: csv.rows };
}

// Multi-job fetch: fetches every tracked job's results in parallel and
// merges them into one table - a synthetic leading "job" column (each
// job's short label, see jobLabel()) plus the union of every job's own
// columns, so jobs with different metric sets (e.g. different task
// types) still combine into one sensible table rather than erroring or
// silently dropping columns. rowMeta keeps each merged row's original
// job + row index around (read by ingestCsv() above) since that's needed
// to resolve a row back to the right job's own API/export calls.
async function fetchAllJobsResults() {
  const jobs = state.jobs.slice();
  const perJob = await Promise.all(jobs.map(async (log) => {
    try {
      const res = await apiGet("/api/results", { log });
      if (!res.csv_exists) return { log, csv_exists: false, final_pdbs: res.final_pdbs || [] };
      const csv = await apiGet("/api/final_csv", { log });
      return { log, csv_exists: true, columns: csv.columns, rows: csv.rows };
    } catch (e) {
      return { log, csv_exists: false, final_pdbs: [] };
    }
  }));

  const withCsv = perJob.filter((j) => j.csv_exists);
  if (!withCsv.length) {
    return {
      csv_exists: false,
      final_pdbs: perJob.flatMap((j) => (j.final_pdbs || []).map((name) => ({ name, log: j.log, jobLabel: jobLabel(j.log) }))),
    };
  }

  const unionColumns = [];
  withCsv.forEach((j) => j.columns.forEach((c) => { if (!unionColumns.includes(c)) unionColumns.push(c); }));
  const columns = ["job", ...unionColumns];
  const rows = [];
  const rowMeta = [];
  withCsv.forEach((j) => {
    const colIndex = {};
    j.columns.forEach((c, i) => { colIndex[c] = i; });
    j.rows.forEach((r, origIdx) => {
      const row = [jobLabel(j.log)];
      unionColumns.forEach((c) => { row.push(c in colIndex ? (r[colIndex[c]] ?? "") : ""); });
      rows.push(row);
      rowMeta.push({ log: j.log, origRowId: origIdx });
    });
  });
  return { csv_exists: true, columns, rows, rowMeta };
}

const resultsView = createResultsView({
  ids: {
    list: "#finalPdbsList", viewer: "#resultViewer", viewerLegend: "#resultsColorLegend",
    seqContent: "#resultsSequenceContent", csvWrap: "#csvTableWrap", alignContent: "#resultsAlignContent",
    count: "#resultsCount", colorByField: "#colorByField", colorByDirection: "#colorByDirection",
    sortByField: "#sortByField", sortByDirection: "#sortByDirection",
    filtersBtn: "#filtersBtn", filtersPanel: "#filtersPanel", filtersList: "#filtersList",
    filtersBadge: "#filtersBadge", clearFiltersBtn: "#clearFiltersBtn", exportBtn: "#exportFilteredBtn",
    structColorMode: "#resultsStructColorMode", chainPaletteWrap: "#resultsChainPaletteWrap",
    listPager: "#resultsListPager", tablePager: "#resultsTablePager",
  },
  viewerSelector: "#resultViewer",
  isMulti: false,
  getLog: () => state.logPath,
  fetchData: fetchSingleJobResults,
  exportUrl: "/api/export_filtered",
  exportFilename: "prosculpt_filtered_export.zip",
  emptyMessage: "No final models yet — this appears once final_operations completes.",
  emptyCsvMessage: "Not available yet — created once final_operations finishes.",
});

const allResultsView = createResultsView({
  ids: {
    list: "#allFinalPdbsList", viewer: "#allResultViewer", viewerLegend: "#allResultsColorLegend",
    seqContent: "#allResultsSequenceContent", csvWrap: "#allCsvTableWrap", alignContent: "#allResultsAlignContent",
    count: "#allResultsCount", colorByField: "#allColorByField", colorByDirection: "#allColorByDirection",
    sortByField: "#allSortByField", sortByDirection: "#allSortByDirection",
    filtersBtn: "#allFiltersBtn", filtersPanel: "#allFiltersPanel", filtersList: "#allFiltersList",
    filtersBadge: "#allFiltersBadge", clearFiltersBtn: "#allClearFiltersBtn", exportBtn: "#allExportFilteredBtn",
    structColorMode: "#allResultsStructColorMode", chainPaletteWrap: "#allResultsChainPaletteWrap",
    listPager: "#allResultsListPager", tablePager: "#allResultsTablePager",
  },
  viewerSelector: "#allResultViewer",
  isMulti: true,
  getLog: (row) => row.__jobLogPath,
  fetchData: fetchAllJobsResults,
  exportUrl: "/api/export_filtered_multi",
  exportFilename: "prosculpt_all_jobs_filtered_export.zip",
  emptyMessage: "No final models yet across any tracked job.",
  emptyCsvMessage: "Not available yet.",
});

// ---------------------------------------------------------------------
// Init
// ---------------------------------------------------------------------

document.addEventListener("DOMContentLoaded", () => {
  initTheme();
  loadViewerSettings();
  loadJobsFromStorage();
  initColorControls();
  initTopbar();
  initAppModeTabs();
  initJobTabsBar();
  initTabs();
  initCancelJobModal();
  resultsView.initAll();
  allResultsView.initAll();
  qs("#modelsFilter").addEventListener("input", renderModelsList);
  qs("#showFilteredToggle").addEventListener("change", (e) => {
    state.showFilteredBackbones = e.target.checked;
    renderBackbonesList();
  });
  qs("#backboneColorMode").addEventListener("change", async (e) => {
    state.backboneColorMode = e.target.value;
    const legendEl = qs("#backboneColorLegend");
    const paletteWrap = qs("#backboneChainPaletteWrap");
    if (state.backboneColorMode === "provenance") {
      paletteWrap.classList.add("hidden");
      const b = (state.backbonesCache || []).find((x) => x.name === state.selectedBackboneKey);
      await applyProvenanceColoring("#backboneViewer", b ? b.trb_path : null, legendEl);
    } else {
      paletteWrap.classList.remove("hidden");
      setCartoonColor("#backboneViewer", chainColorScheme());
      const entry = viewerRegistry["#backboneViewer"];
      if (entry) { entry.colorMode = "chain"; entry.provenanceSchemeId = null; }
      refreshChainLegendFor("#backboneViewer");
    }
  });
  qs("#modelColorMode").addEventListener("change", async (e) => {
    state.modelsColorMode = e.target.value;
    const legendEl = qs("#modelColorLegend");
    const paletteWrap = qs("#modelChainPaletteWrap");
    if (state.modelsColorMode === "provenance") {
      paletteWrap.classList.add("hidden");
      const m = (state.modelsCache || []).find((x) => modelKey(x) === state.selectedModelKey);
      await applyProvenanceColoring("#modelViewer", m ? m.trb_path : null, legendEl);
    } else {
      paletteWrap.classList.remove("hidden");
      setCartoonColor("#modelViewer", chainColorScheme());
      const entry = viewerRegistry["#modelViewer"];
      if (entry) { entry.colorMode = "chain"; entry.provenanceSchemeId = null; }
      refreshChainLegendFor("#modelViewer");
    }
  });
  qs("#conservationToggle").addEventListener("change", renderSequencesContent);
  qs("#conservationThreshold").addEventListener("input", () => {
    if (qs("#conservationToggle").checked) renderSequencesContent();
  });

  startJobStatusPolling();
  if (state.logPath && state.jobs.includes(state.logPath)) {
    showTopView("job");
    refreshAll(true);
    startPolling();
  } else if (state.jobs.length) {
    showTopView("all-overview");
  } else {
    renderJobTabsBar();
  }
});
