/**
 * One display representation for exploratory candidate regions.
 *
 * The coordinate model (src/exondomaincompare/shared_gene_analysis/protein_coordinate_model.py)
 * already groups raw candidate rows into display clusters and assigns each cluster a
 * deterministic rank and lane: rows are merged only when they describe the same
 * primary-protein interval in the same alignment block, so an insertion and a deletion
 * reported over one interval stay two clusters. This module is the single consumer of
 * that layout, so the interactive architecture track and the exported SVG / PDF / PNG
 * cannot disagree about which clusters are shown or which lane a box sits in.
 *
 * Nothing here changes candidate detection, ranking, evidence or interpretation — the
 * raw table stays complete and unaggregated; this only decides what is drawable.
 */

/** Lanes shown before the reader asks for all candidates. */
export const DEFAULT_CANDIDATE_LANES = 3;
/** A box narrower than this many px cannot hold its own label legibly. */
export const MIN_LABEL_PX = 22;
/** Boxes narrower than this are counted into the density strip instead of drawn wide. */
export const DENSE_PX = 3;

export const candidateKey = (c) => c?.display_cluster_id || c?.candidate_id || c?.id;

/**
 * Clusters to draw, plus their lane geometry.
 *
 * @param {Array}  clusters          candidate_regions from the coordinate model
 * @param {Object} opts
 * @param {string} opts.selectedId   cluster that must stay visible whatever the mode
 * @param {boolean} opts.showAll     "Show all candidates"
 * @returns {{visible: Array, laneCount: number, hiddenCount: number, total: number,
 *            laneOf: Function, byLane: Array}}
 */
export function candidateDisplayLayout(clusters, { selectedId = null, showAll = false } = {}) {
  const all = (clusters || []).filter((c) => c && c.start != null && c.end != null);
  const laneOf = (c) => (Number.isInteger(c.display_lane) ? c.display_lane : 0);
  const defaultLanes = all[0]?.default_display_lanes ?? DEFAULT_CANDIDATE_LANES;

  // The selected cluster is never dropped, so selecting a deeply ranked candidate
  // opens exactly as many lanes as it needs and no more.
  const selected = all.find((c) => candidateKey(c) === selectedId) || null;
  const lanesNeeded = showAll
    ? Math.max(0, ...all.map(laneOf)) + 1
    : Math.max(defaultLanes, selected ? laneOf(selected) + 1 : 0);

  const visible = all.filter((c) => laneOf(c) < lanesNeeded);
  const byLane = Array.from({ length: lanesNeeded }, (_, i) => visible.filter((c) => laneOf(c) === i));
  return {
    visible,
    byLane,
    laneOf,
    laneCount: lanesNeeded,
    total: all.length,
    hiddenCount: all.length - visible.length,
  };
}

/** A cluster's label is drawn only when its own box can hold it. */
export const candidateLabelFits = (widthPx) => widthPx >= MIN_LABEL_PX;

/**
 * Compact density indication for a lane: how many of its clusters are too narrow to
 * read at the current scale. Reported rather than silently drawn wider, so a dense
 * stretch of short intervals is visible as a count instead of a smear of boxes.
 */
export function laneDensity(laneClusters, widthOf) {
  let narrow = 0;
  let minStart = null;
  let maxEnd = null;
  for (const c of laneClusters || []) {
    if (widthOf(c) < DENSE_PX) {
      narrow += 1;
      minStart = minStart == null ? c.start : Math.min(minStart, c.start);
      maxEnd = maxEnd == null ? c.end : Math.max(maxEnd, c.end);
    }
  }
  return { narrow, start: minStart, end: maxEnd };
}

/** Full tooltip text: the label plus the support the cluster stands for. */
export function candidateTooltip(c) {
  const parts = [
    `${c.label || candidateKey(c)} · ${c.candidate_type || "candidate"} (exploratory, not validated)`,
    `alignment block ${c.alignment_block || "—"} · confidence ${c.confidence || "—"}`,
    `${c.n_member_candidates ?? 1} raw candidate row(s) · ${c.n_comparisons ?? 1} isoform comparison(s)`
    + ` · ${c.n_supporting_isoforms ?? 0} supporting isoform(s)`,
    `rank ${c.display_rank ?? "—"} · lane ${c.display_lane ?? 0} · source ${c.source_file || "—"}`,
  ];
  if (c.supporting_isoforms?.length) parts.push(`isoforms: ${c.supporting_isoforms.join(", ")}`);
  return parts.join("\n");
}
