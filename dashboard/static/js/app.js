// Prosculpt Dashboard frontend
// Vanilla JS, no build step. Talks to the Flask API in app.py.

const STORAGE_KEY = "prosculpt_dashboard_log_path";
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

let state = {
  logPath: localStorage.getItem(STORAGE_KEY) || "",
  activeTab: "overview",
  lastStatus: null,
  pollTimer: null,
  browsePath: null,
  selectedModelIdx: null,
  modelsCache: [],
  csvCache: null,
  csvSort: { col: null, dir: 1 },
};

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

function qs(sel, root = document) { return root.querySelector(sel); }
function qsa(sel, root = document) { return Array.from(root.querySelectorAll(sel)); }

async function apiGet(path, params = {}) {
  const url = new URL(path, window.location.origin);
  params.log = state.logPath;
  Object.entries(params).forEach(([k, v]) => { if (v !== undefined && v !== null) url.searchParams.set(k, v); });
  const res = await fetch(url);
  if (!res.ok) {
    let msg = res.statusText;
    try { const j = await res.json(); msg = j.error || j.message || msg; } catch (e) {}
    throw new Error(msg);
  }
  return res.json();
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

// ---------------------------------------------------------------------
// Top bar / log loading
// ---------------------------------------------------------------------

function initTopbar() {
  const input = qs("#logPathInput");
  input.value = state.logPath;

  qs("#loadBtn").addEventListener("click", () => {
    state.logPath = input.value.trim();
    localStorage.setItem(STORAGE_KEY, state.logPath);
    refreshAll(true);
  });
  input.addEventListener("keydown", (e) => { if (e.key === "Enter") qs("#loadBtn").click(); });

  qs("#autoRefreshToggle").addEventListener("change", (e) => {
    if (e.target.checked) startPolling(); else stopPolling();
  });

  qs("#browseBtn").addEventListener("click", () => openBrowse(state.browsePath || "."));
  qs("#browseCloseBtn").addEventListener("click", () => qs("#browseModal").classList.add("hidden"));
}

async function openBrowse(path) {
  qs("#browseModal").classList.remove("hidden");
  try {
    const data = await fetch(`/api/browse?path=${encodeURIComponent(path)}`).then((r) => r.json());
    state.browsePath = data.path;
    qs("#browsePath").textContent = data.path;
    const wrap = qs("#browseEntries");
    wrap.innerHTML = "";
    if (data.parent) {
      const up = document.createElement("div");
      up.className = "browse-entry";
      up.innerHTML = "⬆️ ..";
      up.addEventListener("click", () => openBrowse(data.parent));
      wrap.appendChild(up);
    }
    data.entries.forEach((e) => {
      const div = document.createElement("div");
      div.className = "browse-entry";
      div.innerHTML = `${e.is_dir ? "📁" : "📄"} ${escapeHtml(e.name)}`;
      div.addEventListener("click", () => {
        if (e.is_dir) openBrowse(e.path);
        else {
          qs("#logPathInput").value = e.path;
          state.logPath = e.path;
          localStorage.setItem(STORAGE_KEY, state.logPath);
          qs("#browseModal").classList.add("hidden");
          refreshAll(true);
        }
      });
      wrap.appendChild(div);
    });
  } catch (err) {
    qs("#browseEntries").innerHTML = `<p class="muted">${escapeHtml(err.message)}</p>`;
  }
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
  else if (state.activeTab === "results") loadResults();
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
      qs("#overviewCard").innerHTML = `<p class="muted">Fix the log path above and click Load again.</p>`;
      return;
    }
    errBanner.classList.add("hidden");
    renderStepper(status);
    renderOverview(status);
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

function renderOverview(status) {
  const card = qs("#overviewCard");
  const cfg = status.config || {};
  const loc = status.location || {};

  let stageBlock = "";
  if (status.stage === "rfdiffusion" && status.rfdiffusion) {
    stageBlock = renderRfdiffStage(status.rfdiffusion);
  } else if (status.stage === "modeling" && status.modeling) {
    stageBlock = renderModelingStage(status.modeling);
  } else if (status.stage === "mpnn") {
    stageBlock = `<p>Generating sequences with ProteinMPNN — this step is fast and usually finishes in seconds. Check the <b>Sequences</b> tab once it's done.</p>`;
  } else if (status.stage === "filtering") {
    stageBlock = `<p>Filtering backbones and rebuilding chains before ProteinMPNN — a quick housekeeping step.</p>`;
  } else if (status.stage === "scoring") {
    stageBlock = `<p><span class="spinner"></span>Running final scoring (final_operations)… this produces <code>final_pdbs/</code> and <code>final_output.csv</code>. Check the <b>Results</b> tab once it's done.</p>`;
  } else if (status.stage === "finished") {
    const fp = (status.scoring && status.scoring.final_pdbs) || [];
    stageBlock = `<p>✅ Job finished with <b>${fp.length}</b> final model(s). See the <b>Results</b> tab.</p>`;
  } else {
    stageBlock = `<p>Setting up the run (reading config, preparing RFdiffusion)…</p>`;
  }

  const errBlock = (status.possible_errors && status.possible_errors.length)
    ? `<div class="card" style="border-color:var(--danger)"><h4 style="margin-top:0;color:var(--danger)">⚠ Possible errors detected in log</h4><pre style="white-space:pre-wrap;font-size:12px;margin:0">${escapeHtml(status.possible_errors.join("\n"))}</pre></div>`
    : "";

  card.innerHTML = `
    <div class="status-grid">
      <div class="stat-box"><div class="label">Job</div><div class="value small">${escapeHtml(cfg.task_name || "—")}</div></div>
      <div class="stat-box"><div class="label">Prediction model</div><div class="value small">${escapeHtml(cfg.prediction_model || "—")}</div></div>
      <div class="stat-box"><div class="label">Output dir</div><div class="value small" style="word-break:break-all">${escapeHtml(loc.output_dir || "—")}</div></div>
      <div class="stat-box"><div class="label">Current stage</div><div class="value">${STAGE_LABEL[status.stage] || status.stage}</div></div>
    </div>
    <div style="margin-top:18px">${stageBlock}</div>
    ${errBlock}
    <details style="margin-top:14px">
      <summary style="cursor:pointer;color:var(--muted)">Run configuration</summary>
      <table class="config-table">
        ${Object.entries(cfg).map(([k, v]) => `<tr><td>${escapeHtml(k)}</td><td>${escapeHtml(String(v))}</td></tr>`).join("")}
      </table>
    </details>
  `;
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

function renderModelingStage(m) {
  const pct = m.expected_total ? Math.min(100, Math.round((m.completed / m.expected_total) * 100)) : 0;
  return `
    <h4 style="margin-top:0">Structure modeling (${escapeHtml(m.prediction_model || "")})</h4>
    <div class="progress-bar-outer"><div class="progress-bar-inner" style="width:${pct}%"></div></div>
    <div class="status-grid">
      <div class="stat-box"><div class="label">Completed (${escapeHtml(m.unit)})</div><div class="value">${m.completed} / ${m.expected_total ?? "?"}</div></div>
      <div class="stat-box"><div class="label">Avg time / ${escapeHtml(m.unit)}</div><div class="value small">${fmtSeconds(m.avg_seconds)}</div></div>
      <div class="stat-box"><div class="label">Currently on</div><div class="value small" style="word-break:break-all">${m.current_name ? escapeHtml(m.current_name) : "—"}</div></div>
      <div class="stat-box"><div class="label">Est. time remaining</div><div class="value small">${fmtSeconds(m.eta_seconds)}</div></div>
    </div>
    <p class="muted" style="margin-top:10px;font-size:12px">
      ${m.prediction_model && m.prediction_model.toUpperCase().startsWith("BOLTZ")
        ? "Boltz models all sequences for one RFdiffusion backbone in a single batch, so progress is tracked per backbone batch."
        : "AlphaFold3 models one designed sequence at a time, so progress is tracked per sequence."}
    </p>
  `;
}

// ---------------------------------------------------------------------
// Backbones tab
// ---------------------------------------------------------------------

let currentBackboneViewer = null;

async function loadBackbones() {
  const listEl = qs("#backbonesList");
  try {
    const backbones = await apiGet("/api/backbones");
    if (!backbones.length) {
      listEl.innerHTML = `<p class="muted">No backbones generated yet.</p>`;
      return;
    }
    listEl.innerHTML = "";
    backbones.forEach((b) => {
      const div = document.createElement("div");
      div.className = "list-item";
      div.innerHTML = `<span>${escapeHtml(b.name)}</span><span class="badge ${b.status}">${b.status === "passed" ? "ok" : "filtered"}</span>`;
      div.addEventListener("click", () => selectBackbone(b, div));
      listEl.appendChild(div);
    });
    // auto-select first if nothing selected
    if (!qs(".list-item.selected", listEl) && backbones.length) {
      selectBackbone(backbones[0], listEl.children[0]);
    }
  } catch (err) {
    listEl.innerHTML = `<p class="muted">${escapeHtml(err.message)}</p>`;
  }
}

async function selectBackbone(b, el) {
  qsa("#backbonesList .list-item").forEach((i) => i.classList.remove("selected"));
  if (el) el.classList.add("selected");
  qs("#backboneViewerLabel").textContent = `${b.name} (${b.status})`;
  const url = `/api/backbone_pdb?log=${encodeURIComponent(state.logPath)}&name=${encodeURIComponent(b.name)}&status=${b.status}`;
  const pdbText = await fetch(url).then((r) => r.text());
  renderMol("#backboneViewer", pdbText);
}

function renderMol(selector, pdbText) {
  const el = qs(selector);
  el.innerHTML = "";
  const viewer = $3Dmol.createViewer(el, { backgroundColor: "#05070c" });
  viewer.addModel(pdbText, "pdb");
  viewer.setStyle({}, { cartoon: { colorscheme: "chainHetatm" } });
  viewer.zoomTo();
  viewer.render();
  return viewer;
}

// ---------------------------------------------------------------------
// Sequences tab
// ---------------------------------------------------------------------

function coloredSeq(seq) {
  return Array.from(seq).map((c) => {
    const color = AA_COLORS[c] || "#3a4560";
    return `<span class="aa" style="background:${color}1a;color:${color}">${c}</span>`;
  }).join("");
}

function renderSeqGroup(title, backbones) {
  if (!backbones.length) return "";
  let html = `<h3>${escapeHtml(title)}</h3>`;
  backbones.forEach((bb) => {
    html += `<div class="backbone-block"><h4>${escapeHtml(bb.backbone)} <span class="muted" style="font-size:12px;font-weight:400">(${bb.num_samples} sequence${bb.num_samples === 1 ? "" : "s"})</span></h4>`;
    bb.samples.forEach((s, i) => {
      const chainsHtml = s.chains.map((c) => coloredSeq(c)).join('<span class="chain-sep">/</span>');
      const meta = (s.score !== undefined)
        ? `score ${s.score} · global ${s.global_score} · recovery ${s.seq_recovery}`
        : "";
      html += `<div class="seq-row"><span class="seq-label">sample ${i + 1}${meta ? " — " + escapeHtml(meta) : ""}</span><span>${chainsHtml}</span></div>`;
    });
    html += `</div>`;
  });
  return html;
}

async function loadSequences() {
  const el = qs("#sequencesContent");
  try {
    const data = await apiGet("/api/sequences");
    if (!data.backbones.length && !data.monomers.length) {
      el.innerHTML = `<p class="muted">No sequences generated yet.</p>`;
      return;
    }
    el.innerHTML = renderSeqGroup("Designed sequences", data.backbones) + renderSeqGroup("Monomer sequences", data.monomers);
  } catch (err) {
    el.innerHTML = `<p class="muted">${escapeHtml(err.message)}</p>`;
  }
}

// ---------------------------------------------------------------------
// Models tab
// ---------------------------------------------------------------------

async function loadModels() {
  const listEl = qs("#modelsList");
  try {
    const models = await apiGet("/api/models");
    state.modelsCache = models;
    renderModelsList();
    if (!qsa("#modelsList .list-item.selected").length && models.length) {
      selectModel(0);
    }
  } catch (err) {
    listEl.innerHTML = `<p class="muted">${escapeHtml(err.message)}</p>`;
  }
}

function renderModelsList() {
  const listEl = qs("#modelsList");
  const filterVal = (qs("#modelsFilter").value || "").toLowerCase();
  const models = state.modelsCache;
  if (!models.length) {
    listEl.innerHTML = `<p class="muted">No models generated yet.</p>`;
    return;
  }
  listEl.innerHTML = "";
  models.forEach((m, idx) => {
    const label = `${m.model} / ${m.sequence_name}` + (m.model_index !== null && m.model_index !== undefined ? ` (model ${m.model_index})` : "");
    if (filterVal && !label.toLowerCase().includes(filterVal)) return;
    const div = document.createElement("div");
    div.className = "list-item" + (idx === state.selectedModelIdx ? " selected" : "");
    div.innerHTML = `<span style="overflow:hidden;text-overflow:ellipsis">${escapeHtml(label)}</span><span class="badge ${m.variant}">${m.variant}</span>`;
    div.addEventListener("click", () => selectModel(idx));
    listEl.appendChild(div);
  });
}

async function selectModel(idx) {
  state.selectedModelIdx = idx;
  renderModelsList();
  const m = state.modelsCache[idx];
  if (!m) return;

  const metricsEl = qs("#modelMetrics");
  metricsEl.innerHTML = `<p class="muted">Loading…</p>`;

  if (m.pdb_path) {
    const url = `/api/model_pdb?log=${encodeURIComponent(state.logPath)}&path=${encodeURIComponent(m.pdb_path)}`;
    const pdbText = await fetch(url).then((r) => r.text());
    renderMol("#modelViewer", pdbText);
  } else {
    qs("#modelViewer").innerHTML = `<p class="muted" style="padding:20px">No pdb file found yet for this model.</p>`;
  }

  if (m.confidence_path) {
    const url = `/api/model_confidence?log=${encodeURIComponent(state.logPath)}&path=${encodeURIComponent(m.confidence_path)}`;
    const metrics = await fetch(url).then((r) => r.json());
    const rows = Object.entries(metrics).map(([k, v]) => `<tr><td>${escapeHtml(k)}</td><td>${escapeHtml(String(v))}</td></tr>`).join("");
    metricsEl.innerHTML = `<table>${rows}</table>`;
  } else {
    metricsEl.innerHTML = `<p class="muted">No confidence metrics found yet.</p>`;
  }
}

// ---------------------------------------------------------------------
// Results tab
// ---------------------------------------------------------------------

async function loadResults() {
  const listEl = qs("#finalPdbsList");
  const csvWrap = qs("#csvTableWrap");
  try {
    const res = await apiGet("/api/results");
    if (!res.final_pdbs.length) {
      listEl.innerHTML = `<p class="muted">No final models yet — this appears once final_operations completes.</p>`;
    } else {
      listEl.innerHTML = "";
      res.final_pdbs.forEach((name) => {
        const div = document.createElement("div");
        div.className = "list-item";
        div.innerHTML = `<span style="overflow:hidden;text-overflow:ellipsis">${escapeHtml(name)}</span>`;
        div.addEventListener("click", async () => {
          qsa("#finalPdbsList .list-item").forEach((i) => i.classList.remove("selected"));
          div.classList.add("selected");
          const url = `/api/final_pdb?log=${encodeURIComponent(state.logPath)}&name=${encodeURIComponent(name)}`;
          const pdbText = await fetch(url).then((r) => r.text());
          renderMol("#resultViewer", pdbText);
        });
        listEl.appendChild(div);
      });
    }

    if (res.csv_exists) {
      const csv = await apiGet("/api/final_csv");
      state.csvCache = csv;
      renderCsvTable();
    } else {
      csvWrap.innerHTML = `<p class="muted">Not available yet — created once final_operations finishes.</p>`;
    }
  } catch (err) {
    listEl.innerHTML = `<p class="muted">${escapeHtml(err.message)}</p>`;
  }
}

function renderCsvTable() {
  const wrap = qs("#csvTableWrap");
  const { columns, rows } = state.csvCache;
  if (!columns.length) { wrap.innerHTML = `<p class="muted">Empty file.</p>`; return; }

  let sortedRows = rows;
  const { col, dir } = state.csvSort;
  if (col !== null) {
    sortedRows = [...rows].sort((a, b) => {
      const av = a[col], bv = b[col];
      const an = parseFloat(av), bn = parseFloat(bv);
      let cmp;
      if (!isNaN(an) && !isNaN(bn)) cmp = an - bn;
      else cmp = String(av).localeCompare(String(bv));
      return cmp * dir;
    });
  }

  const thead = columns.map((c, i) => `<th data-col="${i}">${escapeHtml(c)}${col === i ? (dir === 1 ? " ▲" : " ▼") : ""}</th>`).join("");
  const tbody = sortedRows.map((r) => `<tr>${r.map((v) => `<td>${escapeHtml(v)}</td>`).join("")}</tr>`).join("");
  wrap.innerHTML = `<table class="data-table"><thead><tr>${thead}</tr></thead><tbody>${tbody}</tbody></table>`;

  qsa("th", wrap).forEach((th) => {
    th.addEventListener("click", () => {
      const c = parseInt(th.dataset.col, 10);
      state.csvSort = { col: c, dir: state.csvSort.col === c ? -state.csvSort.dir : 1 };
      renderCsvTable();
    });
  });
}

// ---------------------------------------------------------------------
// Init
// ---------------------------------------------------------------------

document.addEventListener("DOMContentLoaded", () => {
  initTopbar();
  initTabs();
  qs("#modelsFilter").addEventListener("input", renderModelsList);
  if (state.logPath) refreshAll(true);
  startPolling();
});
