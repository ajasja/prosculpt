// ---------------------------------------------------------------------
// Run Job - config schema
//
// Deliberately NOT a set of rigid per-task-type forms. There's one flat
// set of core fields (present in nearly every prosculpt config) plus a
// checklist of independently-toggleable MODULES that can be layered onto
// any config in any combination (symmetry + filters + inpaint_seq all at
// once is a perfectly valid job) - matching how the actual config schema
// works (see Examples/*.yaml) far better than forcing "binder" vs
// "motif-scaffolding" vs ... into separate mutually-exclusive forms.
//
// contig itself is always free text - RFdiffusion's contig syntax has too
// much expressive range to usefully constrain with form fields, so this
// only offers written guidance (CONTIG_HELP) rather than trying to build
// it from structured inputs.
// ---------------------------------------------------------------------

const CONTIG_HELP = [
  {
    title: "Binder design",
    body:
      "Include the target's own chain/residue range (from the input PDB) unmasked, a 0-length chain break, " +
      "then the binder chain as a length range. Example: a 90-150 residue binder against target chain B, " +
      "residues 1-150:",
    example: "[B1-150/0 90-150]",
  },
  {
    title: "Motif scaffolding",
    body:
      "Mix a fixed motif segment (residues copied from the input PDB) with flexible scaffold ranges around it - " +
      "nothing else is required beyond the contig itself. Example: a fixed motif at chain A residues 163-181, " +
      "flanked by 5-30 residue flexible scaffold on each side:",
    example: "[5-30/A163-181/5-30]",
  },
  {
    title: "Unconditional generation",
    body: "A bare length (or length range) with no chain reference - no input PDB needed. Example: exactly 150 residues:",
    example: "[150-150]",
  },
  {
    title: "Symmetric oligomers",
    body:
      "Usually just a length range for one subunit - RFdiffusion replicates it per the Symmetry module's group. Example:",
    example: "[60-80]",
  },
];

// Core fields - present (as a form control, not necessarily a required
// value) regardless of which modules are active.
const CORE_FIELDS = [
  { key: "job_name", label: "Job name", type: "text", required: true,
    help: "Used as the submitted directory name and as task_name unless overridden below." },
  { key: "task_name", label: "Task name (optional override)", type: "text",
    help: "Defaults to the job name above." },
  { key: "num_tasks", label: "Number of tasks", type: "number", default: 1, min: 1, step: 1 },
  { key: "contig", label: "Contig", type: "textarea", required: true,
    help: "Hover the ? for syntax examples per scenario." },
  { key: "prediction_model", label: "Prediction model", type: "select",
    options: ["Colabfold", "AF3", "Boltz2"], default: "Colabfold" },
  { key: "num_designs_rfdiff", label: "# RFdiffusion designs", type: "number", default: 2, min: 1, step: 1 },
  { key: "num_seq_per_target_mpnn", label: "# MPNN sequences / backbone", type: "number", default: 2, min: 1, step: 1 },
  { key: "af2_mpnn_cycles", label: "MPNN/prediction cycles", type: "number", default: 1, min: 1, step: 1 },
  { key: "num_models", label: "# models / sequence", type: "number", default: 1, min: 1, step: 1 },
  { key: "chain_break_cutoff_A", label: "Chain break cutoff (Å)", type: "number", default: 2, step: 0.1 },
  { key: "sampling_temp", label: "MPNN sampling temperature", type: "number", default: 0.1, step: 0.01 },
  { key: "backbone_noise", label: "Backbone noise", type: "number", default: 0.0, step: 0.01 },
  { key: "use_soluble_model", label: "Use soluble MPNN model", type: "boolean", default: true },
  { key: "model_monomer", label: "Also predict monomer", type: "boolean", default: false,
    tooltip: "Use only when designing binders or symmetry. Do not activate if using templates." },
  { key: "use_a3m", label: "Use custom MSA (.a3m) alignments", type: "boolean", default: false },
  { key: "ppi.hotspot_res", label: "PPI hotspot residues (optional)", type: "list", itemType: "text",
    help: "e.g. A59, A83, A91 - residues on the target used as RFdiffusion binder-design hotspots." },
];

// Independently-toggleable modules. Each is off by default; enabling one
// reveals its own fields (or, for `repeatable` modules, an "Add" button
// building a list of entries from `itemFields`). `configKey` writes a
// top-level list key directly; `passToRfdiff`/`setsConfig` merge into the
// generated config's inference/pass_to_rfdiff machinery instead - see
// buildConfigObject() in run_job.js.
const MODULES = {
  symmetry: {
    label: "Symmetry",
    help: "Diffuse a symmetric oligomer.",
    fields: [
      { key: "inference.symmetry", label: "Symmetry group", type: "select", required: true,
        options: ["c2", "c3", "c4", "c5", "c6", "d2", "d3", "d4", "tetrahedral", "octahedral", "icosahedral"],
        help: "Cyclic (c_n), dihedral (d_n), tetrahedral, octahedral, or icosahedral." },
    ],
  },
  filters: {
    label: "RFdiffusion backbone filters",
    help: "Post-filter generated backbones with a plugin script before they proceed to MPNN. Addable to any job.",
    repeatable: true,
    configKey: "rfdiff_backbone_filters",
    itemFields: [
      { key: "filter_name", label: "Filter name", type: "text", required: true },
      // Uploaded (not typed) - see the "file_upload_small" handling in
      // renderRepeatItems()/collectUploadFiles() in run_job.js. Staged
      // into a filters/ subdirectory of the job's own project directory,
      // with filter_script set to that path automatically - matches how
      // boltz2_templates' own "pdb" field uploads a file rather than
      // asking for a path to one that would have to already exist
      // somewhere the submitting machine (not necessarily this one) can
      // see.
      { key: "filter_script", label: "Filter script", type: "file_upload_small", accept: ".py", required: true,
        help: "Uploaded into this job's filters/ subdirectory; the path below is set automatically." },
      { key: "delete_failed", label: "Delete failed structures", type: "boolean", default: false },
    ],
  },
  inpaint_seq: {
    label: "Inpaint sequence",
    help: "Mask sequence identity in specific residue ranges while keeping the backbone shape there. Addable to any job.",
    fields: [
      { key: "contigmap.inpaint_seq", label: "Inpaint ranges", type: "text", required: true,
        help: "The backbone stays fixed in these ranges but MPNN is free to redesign the sequence.",
        tooltip: "Example: [A10-15/A20-25]" },
    ],
    passToRfdiff: ["contigmap"],
  },
  redesign: {
    label: "Redesign only (skip RFdiffusion)",
    help: "Skip RFdiffusion entirely and redesign sequence directly on the input PDB.",
    requiresPdb: true,
    setsConfig: { skipRfDiff: true },
    usesDenoiser: true,
    fields: [
      { key: "designable_residues", label: "Designable residues", type: "text",
        help: "Residues MPNN may redesign; everything else stays fixed.",
        tooltip: "Example: A1-50. Any chain that should stay entirely fixed (not redesigned) must still be " +
          "listed here too, given as just its chain letter (e.g. B) - not a residue range." },
    ],
  },
  partial_diffusion: {
    label: "Partial diffusion",
    help: "Partially noise and re-diffuse the input structure instead of generating from scratch.",
    requiresPdb: true,
    setsConfig: { partial_diffusion: true },
    passToRfdiff: ["contigmap", "diffuser"],
    usesDenoiser: true,
    fields: [
      { key: "diffuser.partial_T", label: "Partial T", type: "number", required: true, min: 1, step: 1,
        help: "How many forward-diffusion noising steps to apply before reverse diffusion - higher redesigns more of the original structure, lower preserves more of it. Typical range: 5-20." },
      { key: "contigmap.provide_seq", label: "Fixed-sequence ranges (optional)", type: "text",
        help: "Ranges whose sequence identity stays fixed during partial diffusion.",
        tooltip: "Example: [10-20/45-60]. Residue numbers here have no chain letter - they start from 0 and " +
          "keep counting across chain boundaries (they don't reset at the start of each chain)." },
    ],
  },
  boltz2_templates: {
    label: "Boltz2 structural templates",
    help: "Guide Boltz2's prediction with one or more known template structures - only available with prediction model = Boltz2.",
    onlyWhen: (core) => core.prediction_model === "Boltz2",
    repeatable: true,
    configKey: "boltz2_templates",
    itemFields: [
      { key: "pdb", label: "Template PDB", type: "pdb_upload_small", required: true },
      { key: "chain_id", label: "Applies to chain ID(s)", type: "list", itemType: "text", required: true,
        help: "Chain letter(s) of the DESIGNED sequence this template applies to, e.g. B, C, D." },
      { key: "force", label: "Force use of this template", type: "boolean", default: true },
      { key: "threshold", label: "Threshold", type: "number", default: 1, step: 0.1 },
    ],
  },
};

// Fixed order pass_to_rfdiff entries are written in, matching the
// convention every Examples/*.yaml already follows.
const PASS_TO_RFDIFF_ORDER = ["inference", "potentials", "ppi", "denoiser", "contigmap", "diffuser"];
