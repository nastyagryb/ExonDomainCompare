// Run-level lifecycle states, shared by the poller and by the tests that assert when
// polling must stop. Deliberately separate from the analysis-level availability states in
// pages/viewers/common.js: a run finishing and an analysis being possible are different
// questions, and merging the two vocabularies is what made a correct single-exon run look
// unfinished.

/** How often an active run's status is polled. */
export const POLL_INTERVAL_MS = 10000;

/**
 * States at which polling stops because nothing further will change on its own.
 *
 * `results_ready` belongs here even when optional analyses resolved as `not_applicable`:
 * that is a finished run, not a partial one.
 */
export const TERMINAL_RUN_STATES = new Set([
  "results_ready",
  "results_partial",
  "partial",
  "failed",
  "cancelled",
  "deleted",
  "scientifically_unavailable",
]);

/** States in which a run is still progressing and is therefore polled. */
export const ACTIVE_RUN_STATES = new Set([
  "created",
  "precluster_running",
  "cluster_input_ready",
  "cluster_required",
  "cluster_processing",
  "cluster_running",
  "postcluster_running",
  "website_indices_building",
  "running",
]);

/** Whether a run in this state should still be polled. */
export const isActiveRunState = (status) => !TERMINAL_RUN_STATES.has(String(status || ""));

/**
 * Newest run first, by the timestamp the run id starts with.
 *
 * The backend already returns this order; sorting on the client as well means a run that
 * was just created and inserted optimistically lands in the right place before the
 * registry has been asked again.
 */
export function sortNewestFirst(list) {
  return [...(list || [])].sort((a, b) =>
    String(b?.run_id || "").localeCompare(String(a?.run_id || "")));
}
