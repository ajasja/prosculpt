// ---------------------------------------------------------------------
// Run Job tab - build/submit a new job, or point at an already-prepared
// project directory, either locally or (over SSH) on a remote cluster.
// Depends on run_job_schema.js (CORE_FIELDS/MODULES/CONTIG_HELP) and on
// several things already defined in app.js (qs/qsa/escapeHtml/apiGet,
// renderMol/setupResidueInteraction/buildStructureSequence/
// renderSequenceResidues/setupSequencePanelInteraction, addJobs) - loaded
// after both of those, sharing their top-level scope (no build step, no
// module system, matching the rest of this app).
// ---------------------------------------------------------------------

let runJob = {
  submitMode: "build", // "build" | "existing"
  core: {},
  modules: {}, // moduleKey -> { enabled, values: {}, items: [] }
  pdbUpload: null, // { file, name, text }
  a3mUploads: [], // [{ file, name }]
  targets: [],
  suppressYamlSync: false,
  yamlIsValid: true,
};

function initRunJobDefaults() {
  runJob.core = {};
  CORE_FIELDS.forEach((f) => {
    if (f.default !== undefined) runJob.core[f.key] = f.default;
  });
  runJob.modules = {};
  Object.keys(MODULES).forEach((key) => {
    const spec = MODULES[key];
    runJob.modules[key] = { enabled: false, values: {}, items: [] };
    (spec.fields || []).forEach((f) => {
      if (f.default !== undefined) runJob.modules[key].values[f.key] = f.default;
    });
  });
}

// Reads a nested path (e.g. "ppi.hotspot_res") out of a real nested
// object - used only for pulling values back out of a *parsed YAML*
// object (which is genuinely nested), never for runJob's own form state,
// which is stored flat (see the note in renderCoreFields()).
function dottedGet(obj, key) {
  return key.split(".").reduce((o, k) => (o == null ? undefined : o[k]), obj);
}

// Small "?" hover/focus tooltip for field-level guidance (a syntax
// example, a numbering convention) that's worth having but too long to
// leave permanently visible under every field - see .info-tip/
// .info-tip-bubble in style.css. `html` is raw markup, not escaped: every
// call site passes a fixed string written in this file or in
// run_job_schema.js, never anything from user input or a fetched
// response, so there's nothing to sanitize. tabindex="0" lets it be
// reached and shown via keyboard, not just a mouse - see the
// showInfoTip()/hideInfoTip() delegation below for how visibility and
// position are actually driven.
function infoTip(html) {
  return `<span class="info-tip" tabindex="0">?<span class="info-tip-bubble">${html}</span></span>`;
}

// Bubble visibility/position is JS-driven (a .visible class + inline
// top/left), not pure CSS :hover/:focus - the Contig field's tooltip in
// particular is long enough, and close enough to the top of .list-pane
// (which scrolls, `overflow-y: auto`), that a CSS-only bubble anchored
// above the icon got its own top half silently clipped by that ancestor's
// scroll viewport. `position: fixed` (set in style.css) escapes that
// clipping entirely - it's computed against the browser viewport, not
// any scrolling ancestor - but fixed coordinates can't be expressed in
// static CSS, so this measures the icon and bubble on every show and
// picks whichever of above/below actually has room, clamped to the
// viewport on both axes.
let activeInfoTipIcon = null;

function positionInfoTipBubble(icon, bubble) {
  const iconRect = icon.getBoundingClientRect();
  const bubbleW = bubble.offsetWidth;
  const bubbleH = bubble.offsetHeight;
  const margin = 8;
  const spaceBelow = window.innerHeight - iconRect.bottom;
  const spaceAbove = iconRect.top;
  const openAbove = spaceBelow < bubbleH + margin && spaceAbove > spaceBelow;
  const top = openAbove
    ? Math.max(margin, iconRect.top - bubbleH - margin)
    : Math.min(Math.max(margin, window.innerHeight - bubbleH - margin), iconRect.bottom + margin);
  const left = Math.max(margin, Math.min(iconRect.left, window.innerWidth - bubbleW - margin));
  bubble.style.top = `${top}px`;
  bubble.style.left = `${left}px`;
}

function showInfoTip(icon) {
  if (activeInfoTipIcon === icon) return;
  hideInfoTip();
  const bubble = qs(".info-tip-bubble", icon);
  if (!bubble) return;
  activeInfoTipIcon = icon;
  // Made visible off-screen first so offsetWidth/offsetHeight reflect its
  // real rendered size (a display:none element measures as 0x0) before
  // it's actually placed - both DOM writes happen synchronously in this
  // same tick, before the browser's next paint, so there's no visible
  // jump from the off-screen position to the real one.
  bubble.style.top = "-9999px";
  bubble.style.left = "-9999px";
  bubble.classList.add("visible");
  positionInfoTipBubble(icon, bubble);
}

function hideInfoTip() {
  if (!activeInfoTipIcon) return;
  const bubble = qs(".info-tip-bubble", activeInfoTipIcon);
  if (bubble) bubble.classList.remove("visible");
  activeInfoTipIcon = null;
}

// Delegated on document (rather than wired per-icon) since .info-tip
// elements are recreated on nearly every re-render (renderCoreFields(),
// renderModules(), ...) - delegation means nothing ever needs re-binding.
// mouseover/mouseout (not mouseenter/mouseleave) because they bubble,
// which delegation requires; contains(relatedTarget) skips the
// move-within-the-same-icon case (e.g. onto the bubble's own text).
document.addEventListener("mouseover", (e) => {
  const icon = e.target.closest(".info-tip");
  if (icon) showInfoTip(icon);
});
document.addEventListener("mouseout", (e) => {
  const icon = e.target.closest(".info-tip");
  if (icon && !icon.contains(e.relatedTarget)) hideInfoTip();
});
document.addEventListener("focusin", (e) => {
  const icon = e.target.closest(".info-tip");
  if (icon) showInfoTip(icon);
});
document.addEventListener("focusout", (e) => {
  const icon = e.target.closest(".info-tip");
  if (icon && !icon.contains(e.relatedTarget)) hideInfoTip();
});

// The contig field's tooltip is assembled from CONTIG_HELP (rather than a
// static string in run_job_schema.js) so there's exactly one place this
// guidance is written.
function contigHelpTooltipHtml() {
  return CONTIG_HELP.map(
    (h) => `<div class="info-tip-section"><b>${escapeHtml(h.title)}</b><p>${escapeHtml(h.body)}</p><code>${escapeHtml(h.example)}</code></div>`
  ).join("");
}

// ---------------------------------------------------------------------
// Field rendering - one function per field.type, all sharing the same
// (value, onChange) contract so CORE_FIELDS/module fields/repeatable
// item fields can all go through the same renderer.
// ---------------------------------------------------------------------
function renderFieldControl(field, value, onChange) {
  const id = "rjf_" + Math.random().toString(36).slice(2);
  const helpHtml = field.help ? `<div class="rj-field-help">${escapeHtml(field.help)}</div>` : "";
  const tooltipHtml = field.tooltip ? infoTip(field.tooltip) : "";
  let controlHtml = "";
  switch (field.type) {
    case "text":
      controlHtml = `<input type="text" id="${id}" value="${escapeHtml(value ?? "")}">`;
      break;
    case "textarea":
      controlHtml = `<textarea id="${id}" rows="2">${escapeHtml(value ?? "")}</textarea>`;
      break;
    case "number":
      controlHtml = `<input type="number" id="${id}" value="${value ?? ""}" ${field.min != null ? `min="${field.min}"` : ""} ${field.step != null ? `step="${field.step}"` : ""}>`;
      break;
    case "boolean":
      controlHtml = `<label class="toggle-label"><input type="checkbox" id="${id}" ${value ? "checked" : ""}> ${escapeHtml(field.label)}${tooltipHtml}</label>`;
      break;
    case "select":
      controlHtml = `<select id="${id}">${field.options.map((o) => `<option value="${escapeHtml(o)}" ${o === value ? "selected" : ""}>${escapeHtml(o)}</option>`).join("")}</select>`;
      break;
    case "list":
      controlHtml = `<input type="text" id="${id}" placeholder="comma-separated" value="${escapeHtml((value || []).join(", "))}">`;
      break;
    default:
      controlHtml = `<span class="muted">(unsupported field type: ${escapeHtml(field.type)})</span>`;
  }
  const labelHtml = field.type === "boolean" ? "" : `<label for="${id}">${escapeHtml(field.label)}${field.required ? " *" : ""}${tooltipHtml}</label>`;
  const wrap = document.createElement("div");
  wrap.className = "rj-field";
  wrap.innerHTML = `${labelHtml}${controlHtml}${helpHtml}`;
  const input = qs(`#${id}`, wrap);
  const readValue = () => {
    if (field.type === "boolean") return input.checked;
    if (field.type === "number") return input.value === "" ? null : parseFloat(input.value);
    if (field.type === "list") return input.value.split(",").map((s) => s.trim()).filter(Boolean);
    return input.value;
  };
  input.addEventListener(field.type === "select" || field.type === "boolean" ? "change" : "input", () => onChange(readValue()));
  return wrap;
}

// ---------------------------------------------------------------------
// Core fields + PDB upload/viewer + a3m upload
// ---------------------------------------------------------------------
// Field keys are stored FLAT (the literal string "ppi.hotspot_res" as one
// object key, not nested {ppi: {hotspot_res: ...}}) - buildConfigObject()
// already does its own explicit per-field YAML-nesting translation below,
// so form state itself never needs to be a nested structure; a dotted key
// here is just a naming convention (matching the eventual YAML path) for
// readability, not a real path to walk.
function renderCoreFields() {
  const el = qs("#runJobCoreFields");
  el.innerHTML = "";
  CORE_FIELDS.forEach((field) => {
    const value = runJob.core[field.key];
    // The contig field's tooltip is built from CONTIG_HELP at render time
    // rather than stored on the schema entry itself - see
    // contigHelpTooltipHtml().
    const fieldForRender = field.key === "contig" ? { ...field, tooltip: contigHelpTooltipHtml() } : field;
    const ctrl = renderFieldControl(fieldForRender, value, (v) => {
      runJob.core[field.key] = v;
      if (field.key === "prediction_model") renderModules(); // boltz2_templates visibility depends on this
      if (field.key === "use_a3m") qs("#runJobA3mSection").classList.toggle("hidden", !v);
      syncFormToYaml();
    });
    el.appendChild(ctrl);
  });
  qs("#runJobA3mSection").classList.toggle("hidden", !runJob.core.use_a3m);
}

function handlePdbUpload(file) {
  const reader = new FileReader();
  reader.onload = async () => {
    const text = reader.result;
    runJob.pdbUpload = { file, name: file.name, text };
    qs("#runJobPdbFilename").textContent = file.name;
    await renderMol("#runJobPdbViewer", text, "pdb");
    const entry = viewerRegistry["#runJobPdbViewer"];
    if (entry && entry.component) {
      const chains = buildStructureSequence(entry.component.structure);
      qs("#runJobPdbSequencePanel").innerHTML = renderSequenceResidues(chains);
      setupSequencePanelInteraction(qs("#runJobPdbSequencePanel"), "#runJobPdbViewer");
    }
    syncFormToYaml();
  };
  reader.readAsText(file);
}

function clearPdbUpload() {
  runJob.pdbUpload = null;
  qs("#runJobPdbFilename").textContent = "";
  qs("#runJobPdbViewer").innerHTML = "";
  qs("#runJobPdbSequencePanel").innerHTML = "";
  const legendEl = qs("#runJobPdbColorLegend");
  legendEl.innerHTML = "";
  legendEl.classList.add("hidden");
  syncFormToYaml();
}

function renderA3mList() {
  const el = qs("#runJobA3mList");
  if (!runJob.a3mUploads.length) {
    el.innerHTML = `<p class="muted">No alignment files added yet.</p>`;
    return;
  }
  el.innerHTML = runJob.a3mUploads
    .map((u, i) => {
      const warn = /Chain_[A-Za-z]/.test(u.name) ? "" : ` <span class="rj-warn">— doesn't look like it contains "Chain_&lt;letter&gt;"</span>`;
      return `<div class="rj-a3m-row"><span>${escapeHtml(u.name)}</span>${warn}<button type="button" class="rj-remove-btn" data-i="${i}">✕</button></div>`;
    })
    .join("");
  qsa(".rj-remove-btn", el).forEach((btn) => {
    btn.addEventListener("click", () => {
      runJob.a3mUploads.splice(parseInt(btn.dataset.i, 10), 1);
      renderA3mList();
      syncFormToYaml();
    });
  });
}

// ---------------------------------------------------------------------
// Modules
// ---------------------------------------------------------------------
function renderModules() {
  const el = qs("#runJobModules");
  el.innerHTML = "";
  Object.keys(MODULES).forEach((key) => {
    const spec = MODULES[key];
    if (spec.onlyWhen && !spec.onlyWhen(runJob.core)) return;
    const state = runJob.modules[key];
    const box = document.createElement("div");
    box.className = "rj-module" + (state.enabled ? " rj-module-enabled" : "");
    const bodyId = `rjmod_body_${key}`;
    box.innerHTML = `
      <label class="rj-module-toggle">
        <input type="checkbox" id="rjmod_${key}" ${state.enabled ? "checked" : ""}>
        <b>${escapeHtml(spec.label)}</b>
      </label>
      <div class="rj-field-help">${escapeHtml(spec.help || "")}</div>
      <div id="${bodyId}" class="rj-module-body ${state.enabled ? "" : "hidden"}"></div>
    `;
    el.appendChild(box);
    qs(`#rjmod_${key}`, box).addEventListener("change", (e) => {
      state.enabled = e.target.checked;
      qs(`#${bodyId}`, box).classList.toggle("hidden", !state.enabled);
      box.classList.toggle("rj-module-enabled", state.enabled);
      syncFormToYaml();
    });
    const bodyEl = qs(`#${bodyId}`, box);
    if (spec.repeatable) {
      renderRepeatableModuleBody(bodyEl, key, spec, state);
    } else {
      (spec.fields || []).forEach((field) => {
        const value = state.values[field.key];
        const ctrl = renderFieldControl(field, value, (v) => {
          state.values[field.key] = v;
          syncFormToYaml();
        });
        bodyEl.appendChild(ctrl);
      });
    }
  });
}

function renderRepeatableModuleBody(bodyEl, moduleKey, spec, state) {
  const listEl = document.createElement("div");
  listEl.className = "rj-repeat-list";
  bodyEl.appendChild(listEl);
  const addBtn = document.createElement("button");
  addBtn.type = "button";
  addBtn.className = "secondary";
  addBtn.textContent = "+ Add";
  addBtn.addEventListener("click", () => {
    const item = {};
    (spec.itemFields || []).forEach((f) => { if (f.default !== undefined) item[f.key] = f.default; });
    state.items.push(item);
    renderRepeatItems();
    syncFormToYaml();
  });
  bodyEl.appendChild(addBtn);

  function renderRepeatItems() {
    listEl.innerHTML = "";
    state.items.forEach((item, i) => {
      const card = document.createElement("div");
      card.className = "rj-repeat-item";
      const header = document.createElement("div");
      header.className = "rj-repeat-item-header";
      header.innerHTML = `<span class="muted">#${i + 1}</span>`;
      const removeBtn = document.createElement("button");
      removeBtn.type = "button";
      removeBtn.className = "rj-remove-btn";
      removeBtn.textContent = "✕";
      removeBtn.addEventListener("click", () => {
        state.items.splice(i, 1);
        renderRepeatItems();
        syncFormToYaml();
      });
      header.appendChild(removeBtn);
      card.appendChild(header);
      (spec.itemFields || []).forEach((field) => {
        if (field.type === "pdb_upload_small") {
          const wrap = document.createElement("div");
          wrap.className = "rj-field";
          wrap.innerHTML = `<label>${escapeHtml(field.label)}${field.required ? " *" : ""}</label>
            <input type="file" accept=".pdb">
            <span class="muted rj-small-filename">${escapeHtml(item.__pdbFilename || "no file chosen")}</span>`;
          qs("input[type=file]", wrap).addEventListener("change", (e) => {
            const file = e.target.files[0];
            if (!file) return;
            item.__pdbFile = file;
            item.__pdbFilename = file.name;
            qs(".rj-small-filename", wrap).textContent = file.name;
            syncFormToYaml();
          });
          card.appendChild(wrap);
          return;
        }
        // Generic version of the same pattern, for any other repeatable
        // module that needs a small per-item file upload (currently just
        // the filters module's filter_script) - keyed by field.key (not
        // hardcoded like pdb_upload_small above) on item.__uploads so two
        // different file_upload_small fields on the same item type, if
        // that's ever needed, wouldn't clobber each other's stored file.
        if (field.type === "file_upload_small") {
          item.__uploads = item.__uploads || {};
          const current = item.__uploads[field.key];
          const wrap = document.createElement("div");
          wrap.className = "rj-field";
          const helpHtml = field.help ? `<div class="rj-field-help">${escapeHtml(field.help)}</div>` : "";
          wrap.innerHTML = `<label>${escapeHtml(field.label)}${field.required ? " *" : ""}</label>
            <input type="file" ${field.accept ? `accept="${escapeHtml(field.accept)}"` : ""}>
            <span class="muted rj-small-filename">${escapeHtml(current ? current.name : "no file chosen")}</span>
            ${helpHtml}`;
          qs("input[type=file]", wrap).addEventListener("change", (e) => {
            const file = e.target.files[0];
            if (!file) return;
            item.__uploads[field.key] = { file, name: file.name };
            qs(".rj-small-filename", wrap).textContent = file.name;
            syncFormToYaml();
          });
          card.appendChild(wrap);
          return;
        }
        const value = item[field.key];
        const ctrl = renderFieldControl(field, value, (v) => {
          item[field.key] = v;
          syncFormToYaml();
        });
        card.appendChild(ctrl);
      });
      listEl.appendChild(card);
    });
  }
  renderRepeatItems();
}

// Assigns each filters-module item's uploaded script a target path under
// filters/ in the job directory - the original filename by default (so
// the project directory reads naturally if someone browses it later,
// rather than a synthetic name), disambiguated with a _2, _3, ... suffix
// before the extension if two items happen to upload same-named files.
// Without that, two same-named uploads would silently collide into one
// file server-side (both would stage to the same "filters/<name>" path,
// see _collect_uploaded_files() in run_api.py, which is keyed by path) -
// each filter's config entry would still get a filter_script value, but
// one of them would silently be pointed at the wrong script's content.
// Returns items.length paths in the same order as `items`, with null for
// any item that hasn't had a file uploaded yet.
function resolveFilterScriptPaths(items) {
  const used = new Set();
  return items.map((item) => {
    const upload = item.__uploads && item.__uploads.filter_script;
    if (!upload) return null;
    const dot = upload.name.lastIndexOf(".");
    const stem = dot > 0 ? upload.name.slice(0, dot) : upload.name;
    const ext = dot > 0 ? upload.name.slice(dot) : "";
    let candidate = `filters/${upload.name}`;
    let n = 2;
    while (used.has(candidate)) {
      candidate = `filters/${stem}_${n}${ext}`;
      n += 1;
    }
    used.add(candidate);
    return candidate;
  });
}

// ---------------------------------------------------------------------
// Config object -> YAML, and back
// ---------------------------------------------------------------------
function buildConfigObject() {
  const core = runJob.core;
  const mods = runJob.modules;
  const cfg = {};
  cfg.task_name = (core.task_name || core.job_name || "").trim();
  // Empty, not omitted - slurm_runner.py does yaml_data["slurm"] directly
  // (no .get() fallback), so the key has to exist even with nothing under
  // it. Deliberately not setting output/error here: slurm_runner.py forces
  // sbatch's own -o/-e under output_dir/logs/ by default now regardless of
  // what this says, so writing a value here would just be misleading -
  // looks like it's choosing where the log goes when it no longer is.
  cfg.slurm = {};
  cfg.num_tasks = core.num_tasks || 1;
  if (runJob.pdbUpload) cfg.pdb_path = "input.pdb";
  cfg.output_dir = "output";
  cfg.contig = core.contig || "";
  cfg.num_designs_rfdiff = core.num_designs_rfdiff;
  cfg.num_seq_per_target_mpnn = core.num_seq_per_target_mpnn;
  cfg.af2_mpnn_cycles = core.af2_mpnn_cycles;
  cfg.num_models = core.num_models;
  cfg.chain_break_cutoff_A = core.chain_break_cutoff_A;
  cfg.sampling_temp = core.sampling_temp;
  cfg.backbone_noise = core.backbone_noise;
  cfg.use_soluble_model = !!core.use_soluble_model;
  if (core.use_a3m) {
    cfg.use_a3m = true;
    cfg.a3m_dir = "alignments";
  }
  cfg.model_monomer = !!core.model_monomer;
  cfg.prediction_model = core.prediction_model;

  const passToRfdiff = new Set(["inference", "potentials"]);
  cfg.inference = { generatedWithProsculpt: true };
  cfg.potentials = { generatedWithProsculpt: true };

  if (mods.symmetry && mods.symmetry.enabled && mods.symmetry.values["inference.symmetry"]) {
    cfg.inference.symmetry = mods.symmetry.values["inference.symmetry"];
  }

  let needsDenoiser = false;
  const hotspots = (core["ppi.hotspot_res"] || []).filter(Boolean);
  if (hotspots.length) {
    cfg.ppi = { hotspot_res: hotspots };
    passToRfdiff.add("ppi");
    needsDenoiser = true;
  }

  if (mods.redesign && mods.redesign.enabled) {
    cfg.skipRfDiff = true;
    needsDenoiser = true;
    if (mods.redesign.values.designable_residues) cfg.designable_residues = mods.redesign.values.designable_residues;
  }

  if (mods.partial_diffusion && mods.partial_diffusion.enabled) {
    cfg.partial_diffusion = true;
    cfg.diffuser = { partial_T: mods.partial_diffusion.values["diffuser.partial_T"] };
    passToRfdiff.add("contigmap");
    passToRfdiff.add("diffuser");
    needsDenoiser = true;
    const provideSeq = mods.partial_diffusion.values["contigmap.provide_seq"];
    if (provideSeq) {
      cfg.contigmap = cfg.contigmap || {};
      cfg.contigmap.provide_seq = provideSeq;
    }
  }

  if (mods.inpaint_seq && mods.inpaint_seq.enabled) {
    cfg.contigmap = cfg.contigmap || {};
    cfg.contigmap.inpaint_seq = mods.inpaint_seq.values["contigmap.inpaint_seq"];
    passToRfdiff.add("contigmap");
  }

  if (needsDenoiser) {
    cfg.denoiser = { noise_scale_ca: 0, noise_scale_frame: 0 };
    passToRfdiff.add("denoiser");
  }

  cfg.pass_to_rfdiff = PASS_TO_RFDIFF_ORDER.filter((k) => passToRfdiff.has(k));

  if (mods.filters && mods.filters.enabled && mods.filters.items.length) {
    // Resolved (and stashed back onto each item as __resolvedScriptPath)
    // here rather than separately in collectUploadFiles() at submit time -
    // this function already runs on every relevant form change (see
    // syncFormToYaml()), so by submit time it's already up to date, and
    // there's exactly one place computing it instead of two that could
    // disagree about which item got which disambiguated path.
    const scriptPaths = resolveFilterScriptPaths(mods.filters.items);
    mods.filters.items.forEach((it, i) => { it.__resolvedScriptPath = scriptPaths[i]; });
    cfg.rfdiff_backbone_filters = mods.filters.items.map((it, i) => ({
      filter_name: it.filter_name,
      filter_script: scriptPaths[i] || "",
      delete_failed: !!it.delete_failed,
    }));
  }

  if (
    mods.boltz2_templates &&
    mods.boltz2_templates.enabled &&
    mods.boltz2_templates.items.length &&
    core.prediction_model === "Boltz2"
  ) {
    cfg.boltz2_templates = mods.boltz2_templates.items.map((it, i) => ({
      pdb: `boltz_templates/template_${i}.pdb`,
      chain_id: (it.chain_id || []).filter(Boolean),
      force: it.force !== false,
      threshold: it.threshold != null ? it.threshold : 1,
    }));
  }

  cfg.defaults = ["installation", "_self_"];
  return cfg;
}

function syncFormToYaml() {
  if (runJob.suppressYamlSync) return;
  const cfg = buildConfigObject();
  const text = jsyaml.dump(cfg, { lineWidth: -1 });
  runJob.suppressYamlSync = true;
  qs("#runJobYamlEditor").value = text;
  runJob.suppressYamlSync = false;
  setYamlError(null);
}

function setYamlError(message) {
  const el = qs("#runJobYamlError");
  runJob.yamlIsValid = !message;
  if (message) {
    el.textContent = message;
    el.classList.remove("hidden");
  } else {
    el.classList.add("hidden");
  }
  updateSubmitButtonsState();
}

function updateSubmitButtonsState() {
  const disabled = !runJob.yamlIsValid;
  qsa(".rj-submit-btn").forEach((b) => { b.disabled = disabled; });
}

let yamlSyncDebounceTimer = null;
function onYamlEditorInput() {
  if (runJob.suppressYamlSync) return;
  clearTimeout(yamlSyncDebounceTimer);
  yamlSyncDebounceTimer = setTimeout(() => {
    const text = qs("#runJobYamlEditor").value;
    let parsed;
    try {
      parsed = jsyaml.load(text);
    } catch (e) {
      setYamlError("Could not parse YAML: " + e.message);
      return;
    }
    if (!parsed || typeof parsed !== "object") {
      setYamlError("YAML must define a mapping (key: value pairs) at the top level.");
      return;
    }
    if (!parsed.contig) {
      setYamlError('Missing required field: "contig"');
      return;
    }
    // Best-effort form sync from the parsed object - only touches
    // fields that exist in the parsed YAML, leaving everything else
    // (and the raw text itself) untouched rather than guessing.
    runJob.suppressYamlSync = true;
    CORE_FIELDS.forEach((f) => {
      // parsed is the actual (nested) YAML object, so dottedGet's real
      // path-walking is correct here - only the flat runJob.core storage
      // convention needs a plain assignment.
      const v = dottedGet(parsed, f.key);
      if (v !== undefined) runJob.core[f.key] = v;
    });
    renderCoreFields();
    runJob.suppressYamlSync = false;
    setYamlError(null);
  }, 300);
}

// ---------------------------------------------------------------------
// Submission
// ---------------------------------------------------------------------
async function loadRunTargets() {
  const data = await apiGet("/api/run/targets");
  runJob.targets = data.targets || [];
  applyDashboardDefaults(data.defaults || {});
  const sel = qs("#runJobTarget");
  if (!runJob.targets.length) {
    sel.innerHTML = `<option value="">(no targets configured)</option>`;
    qs("#runJobTargetsWarning").classList.remove("hidden");
    qs("#runJobTargetsWarning").textContent =
      data.error || "No run targets configured - copy dashboard/dashboard_config.yaml.example to dashboard_config.yaml and fill in at least one target.";
    return;
  }
  qs("#runJobTargetsWarning").classList.add("hidden");
  sel.innerHTML = runJob.targets
    .map((t) => `<option value="${escapeHtml(t.name)}" ${t.is_default ? "selected" : ""}>${escapeHtml(t.label)} (${t.kind})</option>`)
    .join("");
}

// Applied once, right after the tab's own static defaults have already
// been rendered (initRunJobDefaults()/renderCoreFields() in
// initRunJobTab() run synchronously, well before this fetch resolves) -
// only overrides fields dashboard_config.yaml's `defaults:` actually sets,
// and only if the user hasn't already changed that field away from the
// schema's own static default in the brief window before this resolves.
function applyDashboardDefaults(defaults) {
  if (defaults.prediction_model && runJob.core.prediction_model === "Colabfold") {
    runJob.core.prediction_model = defaults.prediction_model;
    renderCoreFields();
    renderModules(); // boltz2_templates visibility depends on prediction_model
    syncFormToYaml();
  }
}

function collectUploadFiles() {
  const paths = [];
  const files = [];
  if (runJob.pdbUpload) {
    paths.push("input.pdb");
    files.push(runJob.pdbUpload.file);
  }
  runJob.a3mUploads.forEach((u) => {
    paths.push(`alignments/${u.name}`);
    files.push(u.file);
  });
  const boltzMod = runJob.modules.boltz2_templates;
  if (boltzMod && boltzMod.enabled) {
    boltzMod.items.forEach((item, i) => {
      if (item.__pdbFile) {
        paths.push(`boltz_templates/template_${i}.pdb`);
        files.push(item.__pdbFile);
      }
    });
  }
  const filtersMod = runJob.modules.filters;
  if (filtersMod && filtersMod.enabled) {
    // __resolvedScriptPath is set by buildConfigObject() (via
    // resolveFilterScriptPaths()), which syncFormToYaml() already keeps
    // current on every relevant change - reused here rather than
    // recomputed so the path a filter's config entry actually references
    // and the path its script gets staged at can never drift apart.
    filtersMod.items.forEach((item) => {
      const upload = item.__uploads && item.__uploads.filter_script;
      if (upload && item.__resolvedScriptPath) {
        paths.push(item.__resolvedScriptPath);
        files.push(upload.file);
      }
    });
  }
  return { paths, files };
}

function validateBeforeSubmit() {
  if (!runJob.yamlIsValid) return "Fix the YAML error above first.";
  if (runJob.submitMode === "build") {
    if (!(runJob.core.job_name || "").trim()) return "Job name is required.";
    if (!(runJob.core.contig || "").trim()) return "Contig is required.";
    for (const key of Object.keys(runJob.modules)) {
      const spec = MODULES[key];
      const state = runJob.modules[key];
      if (spec.requiresPdb && state.enabled && !runJob.pdbUpload) {
        return `The "${spec.label}" module requires an input PDB - upload one above.`;
      }
    }
  } else {
    if (!(qs("#runJobExistingPath").value || "").trim()) return "Project directory path is required.";
    if (!(qs("#runJobExistingConfigFilename").value || "").trim()) return "Config filename is required.";
  }
  const target = qs("#runJobTarget").value;
  if (!target) return "Choose a run target.";
  return null;
}

async function submitRunJob(dryRun) {
  const err = validateBeforeSubmit();
  if (err) {
    showRunJobResult({ ok: false, stderr: err }, dryRun);
    return;
  }
  const target = qs("#runJobTarget").value;
  const fd = new FormData();
  fd.append("target", target);
  fd.append("dry_run", dryRun ? "1" : "0");
  fd.append("mode", runJob.submitMode);

  if (runJob.submitMode === "build") {
    const jobName = (runJob.core.job_name || "").trim();
    fd.append("job_name", jobName);
    fd.append("config_yaml", qs("#runJobYamlEditor").value);
    const { paths, files } = collectUploadFiles();
    paths.forEach((p) => fd.append("file_paths", p));
    files.forEach((f) => fd.append("file_contents", f));
  } else {
    fd.append("existing_project_path", qs("#runJobExistingPath").value.trim());
    fd.append("config_filename", qs("#runJobExistingConfigFilename").value.trim());
  }

  const btns = qsa(".rj-submit-btn");
  btns.forEach((b) => (b.disabled = true));
  showRunJobWaiting(dryRun);
  try {
    const res = await fetch("/api/run/submit", { method: "POST", body: fd });
    const data = await res.json().catch(() => ({ ok: false, stderr: `HTTP ${res.status}` }));
    showRunJobResult(data, dryRun);
    if (data.ok && !dryRun && data.job_id) {
      addSubmittedJobToTracking(data);
    }
  } catch (e) {
    showRunJobResult({ ok: false, stderr: String(e) }, dryRun);
  } finally {
    btns.forEach((b) => (b.disabled = false));
    stopRunJobWaiting();
  }
}

// The submit/preview round trip runs sbatch (or a dry-run of it) over SSH
// on a real cluster - genuinely slow enough (multi-second, sometimes much
// more) that leaving the button just greyed out with no other feedback
// reads as "did this hang?" rather than "this is working". Shown in the
// same #runJobResult slot showRunJobResult() will overwrite once the real
// response lands - an elapsed counter (not a countdown) since there's no
// reliable estimate of how long a given cluster's login node will take to
// respond, only the hard upper bound (the request timeout) worth naming.
let runJobWaitingTimer = null;

function showRunJobWaiting(dryRun) {
  const el = qs("#runJobResult");
  el.classList.remove("hidden");
  el.className = "rj-result rj-result-waiting";
  const label = dryRun ? "Running preview" : "Submitting";
  const targetLabel = qs("#runJobTarget").selectedOptions[0]?.textContent || "target";
  let elapsed = 0;
  const render = () => {
    el.innerHTML = `<h4>${escapeHtml(label)}…</h4><p class="muted">Waiting for a response from ${escapeHtml(targetLabel)} - this runs over SSH and can take a while (${elapsed}s elapsed, times out at 60s).</p>`;
  };
  render();
  clearInterval(runJobWaitingTimer);
  runJobWaitingTimer = setInterval(() => {
    elapsed += 1;
    render();
  }, 1000);
}

function stopRunJobWaiting() {
  clearInterval(runJobWaitingTimer);
  runJobWaitingTimer = null;
}

function showRunJobResult(data, dryRun) {
  const el = qs("#runJobResult");
  el.classList.remove("hidden");
  const title = dryRun ? "Preview" : data.ok ? "Submitted" : "Submission failed";
  const cls = data.ok ? "rj-result-ok" : "rj-result-error";
  el.className = "rj-result " + cls;
  const renameWarning = data.dir_renamed
    ? `<div class="rj-warn rj-rename-warning">⚠ A directory named "${escapeHtml(data.original_name)}" already existed there - used "${escapeHtml(data.final_name)}" instead.</div>`
    : "";
  const lines = [];
  if (data.job_id) lines.push(`Job ID: ${escapeHtml(data.job_id)}`);
  if (data.remote_job_dir) lines.push(`Remote directory: ${escapeHtml(data.remote_job_dir)}`);
  else if (data.local_job_dir) lines.push(`Directory: ${escapeHtml(data.local_job_dir)}`);
  const body = [data.stdout, data.stderr].filter(Boolean).join("\n");
  el.innerHTML = `<h4>${title}</h4>${renameWarning}${lines.map((l) => `<div>${l}</div>`).join("")}<pre>${escapeHtml(body)}</pre>`;
}

// For a "local" target, the job dir the server just reported IS a path
// this same server can read directly. For an "ssh" target, only a
// *remote* path is known - translatable to a locally-reachable path only
// if that target configured local_mount_path (see
// dashboard_config.yaml.example); if not, there's nothing to add
// automatically and the user is told to add it to Track job by hand once
// they know its path.
//
// remote_job_dir is checked FIRST, not local_job_dir - the backend only
// ever sets both at once as a historical accident to guard against (it
// shouldn't, post-fix, but checking remote first means this function is
// correct even if that ever regresses): whenever a job actually ran on an
// ssh target, remote_job_dir is where its output really lands, never
// whatever local_job_dir happens to say.
function computeLocalJobDir(data) {
  if (data.remote_job_dir) {
    const target = runJob.targets.find((t) => t.name === qs("#runJobTarget").value);
    if (!target || !target.local_mount_path || !target.projects_path) return null;
    if (!data.remote_job_dir.startsWith(target.projects_path)) return null;
    return target.local_mount_path + data.remote_job_dir.slice(target.projects_path.length);
  }
  return data.local_job_dir || null;
}

// Track job now tracks a job by its own output directory directly (see
// get_job_status() in parser.py) rather than by discovering a log file
// that might not exist yet - stage_job() already creates the job
// directory synchronously before sbatch ever runs, so there's nothing to
// wait for: a successful submission can be added to Track job right away.
// This replaces the old poll-for-the-log-file-to-appear "pending runs"
// mechanism entirely - a genuine simplification this whole change enabled,
// not just new complexity to replace old complexity with.
function addSubmittedJobToTracking(data) {
  const jobDir = computeLocalJobDir(data);
  const jobName = data.final_name || (runJob.core.job_name || "").trim() || "(unnamed)";
  if (!jobDir) {
    showRunJobToast(`"${jobName}" submitted (job ${data.job_id}) - this target has no local_mount_path configured, so it can't be added to Track job automatically; add it by hand once you know its directory.`);
    return;
  }
  // false: this is an automatic hand-off, not something the user just
  // asked to see - don't yank them over to the Track job tab.
  addJobs(jobDir, false);
  showRunJobToast(`"${jobName}" submitted (job ${data.job_id}) - added to Track job.`);
}

function showRunJobToast(message) {
  const el = qs("#runJobToast");
  if (!el) return;
  el.textContent = message;
  el.classList.remove("hidden");
  clearTimeout(el._hideTimer);
  el._hideTimer = setTimeout(() => el.classList.add("hidden"), 8000);
}

// ---------------------------------------------------------------------
// squeue
// ---------------------------------------------------------------------
async function refreshSqueue() {
  const target = qs("#runJobTarget").value;
  const el = qs("#runJobSqueueOutput");
  if (!target) {
    el.textContent = "Choose a run target first.";
    return;
  }
  el.textContent = "Loading…";
  try {
    const res = await fetch("/api/run/squeue", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ target }),
    });
    const data = await res.json().catch(() => ({ stdout: "", stderr: `HTTP ${res.status}` }));
    el.textContent = data.stdout || data.stderr || "(no output)";
  } catch (e) {
    el.textContent = String(e);
  }
}

// ---------------------------------------------------------------------
// Build pane split-view resizer (drag the handle between the form and the
// PDB viewer). The viewer itself doesn't need any explicit resize call
// here - .mol-viewer already stretches to fill .viewer-pane's width (flex
// column, default align-items:stretch) and ensureResizeObserver() (see
// app.js) is already watching #runJobPdbViewer's own size and calls
// stage.handleResize() whenever it changes, so dragging the column
// automatically keeps NGL's canvas in sync for free.
// ---------------------------------------------------------------------
const RJ_SPLIT_WIDTH_KEY = "prosculpt_dashboard_run_job_split_width";
const RJ_SPLIT_DEFAULT_WIDTH = 680;
const RJ_SPLIT_MIN_WIDTH = 340; // never narrower than the pre-resize default
const RJ_SPLIT_HANDLE_WIDTH = 6;
const RJ_SPLIT_VIEWER_MIN_WIDTH = 320; // leave the viewer pane at least this wide

function applyRunJobSplitWidth(px) {
  const splitEl = qs("#runJobSplitView");
  if (splitEl) splitEl.style.gridTemplateColumns = `${px}px ${RJ_SPLIT_HANDLE_WIDTH}px 1fr`;
}

// Clamps against the split-view's OWN current width, not the window's -
// .list-pane/.viewer-pane live inside whatever the page's layout gives
// .split-view, which isn't full window width. Returns `px` unchanged if
// the element can't be measured right now (e.g. Run job tab isn't the
// active one, so it's display:none and reports 0 width) rather than
// clamping against a bogus 0, since the stored width is still valid for
// whenever it next becomes visible.
function clampRunJobSplitWidth(px) {
  const splitEl = qs("#runJobSplitView");
  const totalWidth = splitEl ? splitEl.getBoundingClientRect().width : 0;
  if (!totalWidth) return Math.max(RJ_SPLIT_MIN_WIDTH, px);
  const maxWidth = Math.max(RJ_SPLIT_MIN_WIDTH, totalWidth - RJ_SPLIT_HANDLE_WIDTH - RJ_SPLIT_VIEWER_MIN_WIDTH);
  return Math.min(maxWidth, Math.max(RJ_SPLIT_MIN_WIDTH, px));
}

function loadRunJobSplitWidth() {
  const stored = parseInt(localStorage.getItem(RJ_SPLIT_WIDTH_KEY), 10);
  return Number.isFinite(stored) && stored > 0 ? stored : RJ_SPLIT_DEFAULT_WIDTH;
}

function initRunJobSplitResizer() {
  const splitEl = qs("#runJobSplitView");
  const handle = qs("#runJobSplitResizer");
  if (!splitEl || !handle) return;
  applyRunJobSplitWidth(clampRunJobSplitWidth(loadRunJobSplitWidth()));

  handle.addEventListener("mousedown", (e) => {
    e.preventDefault();
    const startX = e.clientX;
    const startWidth = qs(".list-pane", splitEl).getBoundingClientRect().width;
    handle.classList.add("dragging");
    document.body.classList.add("rj-resizing");

    function onMove(ev) {
      const next = clampRunJobSplitWidth(startWidth + (ev.clientX - startX));
      applyRunJobSplitWidth(next);
    }
    function onUp() {
      document.removeEventListener("mousemove", onMove);
      document.removeEventListener("mouseup", onUp);
      handle.classList.remove("dragging");
      document.body.classList.remove("rj-resizing");
      const finalWidth = Math.round(qs(".list-pane", splitEl).getBoundingClientRect().width);
      localStorage.setItem(RJ_SPLIT_WIDTH_KEY, String(finalWidth));
    }
    document.addEventListener("mousemove", onMove);
    document.addEventListener("mouseup", onUp);
  });

  // Double-click resets back to the original default width.
  handle.addEventListener("dblclick", () => {
    localStorage.removeItem(RJ_SPLIT_WIDTH_KEY);
    applyRunJobSplitWidth(RJ_SPLIT_DEFAULT_WIDTH);
  });

  // Re-clamp on window resize so a previously-valid stored width can't
  // leave the viewer pane crushed (or the handle pushed off-screen) after
  // the browser window itself gets narrower.
  window.addEventListener("resize", () => {
    applyRunJobSplitWidth(clampRunJobSplitWidth(loadRunJobSplitWidth()));
  });
}

// ---------------------------------------------------------------------
// Mode switch (build vs existing) and init
// ---------------------------------------------------------------------
function setSubmitMode(mode) {
  runJob.submitMode = mode;
  qs("#runJobBuildPane").classList.toggle("hidden", mode !== "build");
  qs("#runJobExistingPane").classList.toggle("hidden", mode !== "existing");
  qsa(".rj-mode-btn").forEach((b) => b.classList.toggle("active", b.dataset.mode === mode));
}

function initRunJobTab() {
  initRunJobDefaults();
  renderCoreFields();
  renderModules();
  renderA3mList();
  syncFormToYaml();
  initRunJobSplitResizer();

  qs("#runJobPdbUpload").addEventListener("change", (e) => {
    if (e.target.files[0]) handlePdbUpload(e.target.files[0]);
  });
  qs("#runJobPdbClear").addEventListener("click", clearPdbUpload);
  qs("#runJobA3mUpload").addEventListener("change", (e) => {
    Array.from(e.target.files).forEach((f) => runJob.a3mUploads.push({ file: f, name: f.name }));
    e.target.value = "";
    renderA3mList();
    syncFormToYaml();
  });
  qs("#runJobYamlEditor").addEventListener("input", onYamlEditorInput);

  qsa(".rj-mode-btn").forEach((btn) => btn.addEventListener("click", () => setSubmitMode(btn.dataset.mode)));
  setSubmitMode("build");

  qs("#runJobPreviewBtn").addEventListener("click", () => submitRunJob(true));
  qs("#runJobSubmitBtn").addEventListener("click", () => {
    if (!confirm("Submit this job to the cluster now? (Use Preview first if you're not sure.)")) return;
    submitRunJob(false);
  });
  qs("#runJobSqueueRefreshBtn").addEventListener("click", refreshSqueue);

  loadRunTargets();
}
