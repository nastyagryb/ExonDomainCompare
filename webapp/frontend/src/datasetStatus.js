// One canonical dataset-status vocabulary shared by every generic page/badge so
// page-level and dataset-level status can never disagree (Part B). Gene symbols
// are never referenced. The single source of truth is the validated
// protein-coordinate model (its per-model status + boundary_dashboard.page_mode).

export const DATASET_STATUS_META = {
  pre_cluster_ready: ["neutral", "Pre-cluster ready"],
  cluster_processing: ["minor", "Cluster processing"],
  post_cluster_partial: ["minor", "Post-cluster partial"],
  results_ready: ["accepted", "Results ready"],
  failed: ["excluded", "Failed"],
  unavailable: ["neutral", "Unavailable"],
};

// Resolve the canonical status from the coordinate model (+ optional dataset info
// for the mid-cluster "processing" and "failed" transients, which the model alone
// cannot express). Never returns results_ready while the domain layer is pending.
export function datasetStatusFromModel(coordModel, datasetInfo) {
  const models = coordModel?.models || [];
  const runStatus = String(datasetInfo?.status || "").toLowerCase();
  if (runStatus === "failed") return "failed";
  if (!models.length) {
    if (["cluster_submitted", "cluster_running", "cluster_processing"].some((s) => runStatus.includes(s)))
      return "cluster_processing";
    return "unavailable";
  }
  const pageMode = coordModel?.boundary_dashboard?.page_mode;
  const allAvailable = models.every((m) => m.status === "available");
  if (pageMode && pageMode.includes("results_ready")) {
    return runStatus === "results_partial" ? "post_cluster_partial" : "results_ready";
  }
  if (allAvailable) return runStatus === "results_partial" ? "post_cluster_partial" : "results_ready";
  if (["cluster_submitted", "cluster_running", "cluster_processing"].some((s) => runStatus.includes(s)))
    return "cluster_processing";
  return "pre_cluster_ready";
}

export function datasetStatusLabel(status) {
  return DATASET_STATUS_META[status]?.[1] || "Unavailable";
}

// The run API still speaks its own pipeline vocabulary ("cluster_required",
// "results_partial", "created", …). Map it onto the canonical six states so a
// dataset badge can never display a token no other view knows about.
const RUN_STATUS_ALIASES = {
  cluster_required: "pre_cluster_ready",
  precluster_ready: "pre_cluster_ready",
  pre_interpro_ready: "pre_cluster_ready",
  models_ready: "pre_cluster_ready",
  created: "pre_cluster_ready",
  cluster_submitted: "cluster_processing",
  cluster_running: "cluster_processing",
  running: "cluster_processing",
  results_partial: "post_cluster_partial",
  post_interpro_partial: "post_cluster_partial",
  complete: "results_ready",
  completed: "results_ready",
  error: "failed",
  missing: "unavailable",
};

export function normalizeDatasetStatus(status) {
  const s = String(status || "").toLowerCase();
  if (!s) return "unavailable";
  if (DATASET_STATUS_META[s]) return s;
  if (RUN_STATUS_ALIASES[s]) return RUN_STATUS_ALIASES[s];
  if (s.includes("fail") || s.includes("error")) return "failed";
  if (s.includes("running") || s.includes("processing")) return "cluster_processing";
  if (s.includes("partial")) return "post_cluster_partial";
  if (s.includes("cluster")) return "pre_cluster_ready";
  return "unavailable";
}

// For post_cluster_partial / pre_cluster_ready: name the exact post-cluster
// products that are still missing, per species, instead of a generic "pending".
export function missingPostClusterItems(coordModel) {
  const out = [];
  for (const m of coordModel?.models || []) {
    if (m.status === "available") continue;
    const layers = [];
    if (!(m.representative_domains || []).length) layers.push("representative InterPro domains");
    if (m.tm_analysis?.pending || !(m.tm_analysis && "n_tm_regions" in m.tm_analysis))
      layers.push("pyTMHMM topology");
    if (!(m.exon_boundaries || []).some((b) => b.signed_distance != null))
      layers.push("exon–domain boundary classification");
    out.push({
      species: m.species_id || m.species || m.protein_id || "species",
      missing: layers.length ? layers : ["post-cluster analysis outputs"],
    });
  }
  return out;
}
