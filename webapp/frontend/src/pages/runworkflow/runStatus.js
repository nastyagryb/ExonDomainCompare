// Shared reading of a run's status, so the list card and the detail header
// cannot disagree about what a status means.

export const RUNNING_STATUSES = new Set([
  "running", "pre_interpro_running", "post_interpro_running", "cluster_running",
]);

export const FAILED_STATUSES = new Set([
  "failed", "core_model_collection_failed", "incomplete",
]);

export function statusBadgeCls(status) {
  if (status === "results_ready") return "accepted";
  if (FAILED_STATUSES.has(status)) return "excluded";
  if (RUNNING_STATUSES.has(status)) return "info";
  return "neutral";
}
