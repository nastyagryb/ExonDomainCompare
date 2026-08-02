// The shared synteny display model.
//
// One normaliser feeds the interactive viewer and the SVG/PDF/PNG exports, so
// what a reader sees on screen and what ends up in the paper are the same set
// of loci in the same order with the same styles. The backend contract
// (`shared_synteny_v1`) is the authority; the legacy shape that older indices
// still carry is projected onto it here rather than being rendered separately.
//
// Two invariants hold for every gene and every dataset:
//   * the target locus occupies its own central slot and is never counted as an
//     upstream or downstream neighbour;
//   * displayed counts are stated, never implied — a species with four real
//     downstream genes reads as four, and no fifth gene is invented to fill the
//     grid.

import { PALETTE } from "./figureSpec.js";
import { FEATURE_STYLES } from "./semanticStyles.js";

export const TARGET_SLOT = 0;
export const DEFAULT_PER_SIDE = 5;

/** Readable label and exact definition per orthology class (mirrors the backend). */
export const ORTHOLOGY_CLASSES = {
  target: {
    label: "Target gene",
    definition: "The gene this analysis is about, shown in its own central slot.",
    fill: FEATURE_STYLES.synteny_target.fill,
    stroke: FEATURE_STYLES.synteny_target.stroke,
    text: FEATURE_STYLES.synteny_target.text,
  },
  exact: {
    label: "Resolved ortholog",
    definition: "The neighbouring locus carries a curated gene symbol that matches "
      + "the reference symbol exactly.",
    fill: "#E3F0E8", stroke: "#2F9E6F", text: "#1C6B4B",
  },
  curated: {
    label: "Curated ortholog",
    definition: "The neighbouring locus was assigned by a curated orthology resource "
      + "rather than by symbol matching alone.",
    fill: "#E3F0E8", stroke: "#2F9E6F", text: "#1C6B4B",
  },
  rbh: {
    label: "Best reciprocal hit",
    definition: "The neighbouring locus is the reciprocal best protein hit against "
      + "the reference proteome.",
    fill: "#EAF1FA", stroke: "#3F7FC0", text: "#285A8F",
  },
  weak: {
    label: "Weak hit",
    definition: "A protein hit exists but falls below the identity or coverage "
      + "threshold used for a confident assignment.",
    fill: "#FDF6E9", stroke: "#C98A1F", text: "#8A5C11", dashed: true,
  },
  ambiguous: {
    label: "Ambiguous paralog",
    definition: "Several loci hit the same reference gene, so the assignment cannot "
      + "be resolved to one ortholog.",
    fill: "#FDF1DE", stroke: "#C97A16", text: "#8A5C11",
  },
  placeholder: {
    label: "Placeholder locus",
    definition: "Placeholder locus label; curated gene symbol unavailable. "
      + "The genomic position is known.",
    fill: "#F5F1E6", stroke: "#B8A272", text: "#7A5A1C", italic: true,
  },
  unresolved: {
    label: "Unresolved locus",
    definition: "The locus is annotated in the assembly but no orthology assignment "
      + "was attempted or succeeded.",
    fill: PALETTE.paper, stroke: PALETTE.axis, text: PALETTE.muted,
  },
};

export const orthologyStyle = (cls) =>
  ORTHOLOGY_CLASSES[cls] || ORTHOLOGY_CLASSES.unresolved;

const PLACEHOLDER = /^(?:LOC\d+|GENE\d+|ENS\w*G\d+)$/i;

export const isPlaceholderLocus = (symbol) => {
  const text = String(symbol || "").trim();
  return !text || PLACEHOLDER.test(text);
};

const readableStatus = (status) => String(status || "")
  .replace(/_/g, " ").replace(/^\w/, (c) => c.toUpperCase());

/** Project one legacy neighbour node onto the canonical locus shape. */
function legacyLocus(node) {
  const isTarget = Boolean(node.is_target ?? node.is_anchor);
  const symbol = String(node.symbol || node.raw_symbol || "").trim();
  const placeholder = !isTarget && isPlaceholderLocus(symbol);
  const legacyClass = node.method_class === "anchor" ? "target" : node.method_class;
  const cls = node.orthology_class
    || (isTarget ? "target" : placeholder ? "placeholder" : legacyClass) || "unresolved";
  const slot = node.slot_x ?? 0;
  return {
    ...node,
    slot_x: slot,
    side: node.side || (isTarget ? "target" : slot < 0 ? "upstream" : "downstream"),
    rank: node.rank ?? Math.abs(slot),
    is_target: isTarget,
    symbol: symbol || node.gene_id || "",
    source_symbol: node.source_symbol ?? node.raw_symbol ?? symbol,
    resolved_symbol: node.resolved_symbol ?? (placeholder ? "" : symbol),
    placeholder,
    orthology_class: cls,
    orthology_label: node.orthology_label || orthologyStyle(cls).label,
    orthology_definition: node.orthology_definition || orthologyStyle(cls).definition,
  };
}

/**
 * One species row in the shape both renderers consume.
 *
 * Accepts a canonical contract row, a legacy `neighbors5`/`neighbors10` row, or
 * the flat single-species index, and always returns explicit counts so the view
 * can state what it is showing instead of implying a full grid.
 */
export function normaliseSpeciesRow(row, { gene = "", perSide = DEFAULT_PER_SIDE } = {}) {
  if (!row) return null;

  let loci = Array.isArray(row.loci) && row.loci.length ? row.loci.map(legacyLocus) : null;
  if (!loci) {
    const legacy = row.neighbors10?.length ? row.neighbors10
      : (row.neighbors5?.length ? row.neighbors5 : (row.neighbors || row.neighbours || []));
    loci = legacy.map(legacyLocus);
  }

  const flanking = loci.filter((n) => !n.is_target);
  const upAll = flanking.filter((n) => n.side === "upstream")
    .sort((a, b) => a.rank - b.rank);
  const downAll = flanking.filter((n) => n.side === "downstream")
    .sort((a, b) => a.rank - b.rank);
  const limit = Math.max(0, Number(row.requested_neighbour_count) || perSide)
    || Math.max(upAll.length, downAll.length);
  const upstream = upAll.slice(0, limit).reverse();
  const downstream = downAll.slice(0, limit);

  const target = loci.find((n) => n.is_target) || legacyLocus({
    is_target: true, slot_x: TARGET_SLOT, side: "target",
    symbol: row.target_symbol || gene || "target gene",
    strand: row.target_strand || "",
  });

  const displayedUp = upstream.length;
  const displayedDown = downstream.length;
  const total = displayedUp + displayedDown;
  const countsLabel = row.counts_label || (total
    ? `${total} flanking ${total === 1 ? "locus" : "loci"} shown · `
      + `${displayedUp} upstream · ${displayedDown} downstream`
    : "No flanking loci available");

  return {
    speciesId: row.species_id || row.species || "",
    displayName: row.display_species_name || row.species_id || row.species || "",
    gene: row.gene_symbol || gene,
    target,
    upstream,
    downstream,
    loci: [...upstream, target, ...downstream],
    counts: {
      upstreamAvailable: row.upstream_count_available ?? upAll.length,
      downstreamAvailable: row.downstream_count_available ?? downAll.length,
      displayedUpstream: displayedUp,
      displayedDownstream: displayedDown,
      requested: limit,
    },
    countsLabel,
    truncationStatus: row.truncation_status
      || (displayedUp < upAll.length || displayedDown < downAll.length
        ? "truncated_to_request"
        : (displayedUp < limit || displayedDown < limit ? "fewer_available" : "complete")),
    omissionReason: row.omission_reason || "",
    targetPosition: row.target_position || "",
    targetCoordinateSource: row.target_coordinate_source || "annotation",
    statusLabel: row.synteny_status_label
      || readableStatus(row.synteny_status || row.synteny_class),
    statusDefinition: row.synteny_status_definition || "",
    comparisonAvailable: Boolean(row.comparison_available),
    isReview: Boolean(row.is_review),
    isHumanReference: Boolean(row.is_human_reference_control),
    classesPresent: [...new Set([...upstream, target, ...downstream]
      .map((n) => n.orthology_class))],
  };
}

/** Every species row of an index, in index order, already normalised. */
export function normaliseSyntenyIndex(data, { gene = "" } = {}) {
  const g = data?.gene_symbol || data?.target_symbol || gene;
  const rows = (data?.species || []).map((r) => normaliseSpeciesRow(r, { gene: g }));
  if (rows.length) return rows.filter(Boolean);
  // Flat single-species index (older generic runs).
  const flat = data?.neighbours || data?.neighbors;
  if (!flat?.length) return [];
  const one = normaliseSpeciesRow({
    species_id: data.species_id,
    display_species_name: data.species_id,
    gene_symbol: g,
    neighbors: flat.map((n) => ({
      ...n,
      symbol: n.symbol || n.neighbor_symbol,
      slot_x: (Number(n.order) || 0) * (n.side === "upstream" ? -1 : 1),
      strand: n.orientation || n.strand,
      method_class: (n.resolved || n.status === "resolved") ? "exact" : "unresolved",
    })),
  }, { gene: g });
  return one ? [one] : [];
}

/**
 * Column grid shared by every rendered row.
 *
 * The widest displayed side decides the half-width, so the target sits at the
 * exact centre and a human reference row lines up slot-for-slot with the
 * species row even when one of them has fewer real neighbours.
 */
export function slotGrid(rows) {
  const perSide = Math.max(1, ...rows.filter(Boolean).flatMap((r) => [
    r.counts.displayedUpstream, r.counts.displayedDownstream,
  ]));
  return {
    perSide,
    columns: perSide * 2 + 1,
    targetColumn: perSide,
    columnOf: (locus) => perSide + (locus.is_target ? 0 : locus.slot_x),
  };
}

/**
 * Per-species neighbourhood coverage of a rendered set of rows.
 *
 * A comparative figure has to say which species it actually covers, so a reader
 * never mistakes a silently dropped species for one without neighbours. Rows
 * that carry no flanking locus at all count as unavailable rather than complete.
 */
export function syntenyCoverage(rows) {
  const present = (rows || []).filter(Boolean);
  const label = (r) => r.displayName || r.speciesId;
  const complete = [];
  const partial = [];
  const unavailable = [];
  for (const r of present) {
    const flanking = r.counts.displayedUpstream + r.counts.displayedDownstream;
    if (!flanking) unavailable.push(label(r));
    else if (r.truncationStatus === "complete"
      || r.truncationStatus === "truncated_to_request") complete.push(label(r));
    else partial.push(label(r));
  }
  return {
    requested: present.length,
    shown: present.length,
    complete,
    partial,
    unavailable,
  };
}

/**
 * An explicit "no neighbourhood resolved" row for a dataset species that the
 * synteny index does not cover, so the species stays visible instead of being
 * dropped from a comparative figure.
 */
export function unresolvedSpeciesRow({ speciesId, displayName, gene, reason = "" }) {
  return normaliseSpeciesRow({
    species_id: speciesId,
    display_species_name: displayName || speciesId,
    gene_symbol: gene,
    target_symbol: gene,
    neighbors: [],
    truncation_status: "unavailable",
    omission_reason: reason || "no local gene neighbourhood resolved for this species",
    counts_label: "No flanking loci available",
    synteny_status_label: "neighbourhood unavailable",
  }, { gene });
}

/** Only the classes actually drawn, in a stable reading order. */
export function legendEntries(rows) {
  const order = ["target", "exact", "curated", "rbh", "weak", "ambiguous",
    "placeholder", "unresolved"];
  const present = new Set(rows.filter(Boolean).flatMap((r) => r.classesPresent));
  return order.filter((c) => present.has(c))
    .map((c) => ({ cls: c, ...orthologyStyle(c) }));
}

/** The exact rows behind a rendered view, for the source TSV. */
export function syntenyRowsTsv(rows) {
  const cols = ["species_id", "slot_x", "side", "rank", "is_target", "symbol",
    "source_symbol", "resolved_symbol", "gene_id", "placeholder", "orthology_class",
    "mapping_confidence", "strand", "distance", "seqid", "genomic_start", "genomic_end"];
  const out = [cols.join("\t")];
  for (const r of rows.filter(Boolean)) {
    for (const n of r.loci) {
      out.push(cols.map((c) => {
        if (c === "species_id") return r.speciesId;
        const v = n[c];
        return v === null || v === undefined ? "" : String(v);
      }).join("\t"));
    }
  }
  return `${out.join("\n")}\n`;
}
