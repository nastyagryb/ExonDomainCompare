// The shared scientific visual specification.
//
// One semantic style per scientific feature kind, used by BOTH the interactive
// React viewers and the SVG/PDF/PNG publication renderers. Before this module the
// two sides carried independent colour tables — the interactive exon was #c9d3e2
// while the exported exon was #A9BED4, and the boundary classes disagreed
// entirely — so "the figure on screen and the figure in the paper" were only
// nominally the same figure.
//
// Rules for using it:
//   * A scientific colour, stroke or marker is defined here and nowhere else.
//     Component CSS may still own layout, spacing, cursors, transitions and hover.
//   * Interactive marks pass these values as explicit SVG attributes, so a figure
//     stays legible when the component stylesheet is absent — which is exactly the
//     situation of an exported or standalone SVG.
//   * `labelPriority` orders label placement when space is short: higher wins.
//   * `printFallback` is the grey a monochrome print keeps, for reviewers who
//     print in black and white.
//
// Colours come from the publication palette in figureSpec.js, which is
// colour-blind-safe (Okabe–Ito derived). The interactive side adopts them rather
// than the reverse, so the accepted exported figures stay byte-identical.

import { PALETTE } from "./figureSpec.js";

/** Semantic keys. Using a symbol that is not listed here is a programming error. */
export const FEATURE_KEYS = [
  "coding_exon", "alternative_exon", "shared_exon", "shifted_boundary",
  "representative_domain", "family_superfamily", "tm_helix", "candidate_region",
  "boundary_exact", "boundary_near", "boundary_inside", "boundary_outside",
  "boundary_uncertain", "selected_feature", "validated_event",
  "primary_sequence", "alternative_sequence", "gap",
  "variable_region", "conserved_region",
  "protein_backbone", "exon_boundary_tick", "member_signature",
  "functional_site", "disorder_region",
  "synteny_target",
];

const F = (fill, stroke, extra = {}) => ({
  fill,
  stroke,
  strokeWidth: 0.8,
  opacity: 1,
  text: PALETTE.ink,
  marker: null,
  labelPriority: 5,
  printFallback: "#8a8a8a",
  ...extra,
});

export const FEATURE_STYLES = {
  // --- exon structure ------------------------------------------------------ //
  coding_exon: F(PALETTE.exon, PALETTE.exonEdge, {
    text: PALETTE.ink, labelPriority: 6, printFallback: "#c4c4c4",
  }),
  shared_exon: F(PALETTE.exon, PALETTE.exonEdge, {
    text: PALETTE.ink, labelPriority: 5, printFallback: "#c4c4c4",
  }),
  alternative_exon: F(PALETTE.exonAlt, "#B8801A", {
    text: PALETTE.ink, labelPriority: 8, printFallback: "#6f6f6f",
  }),
  shifted_boundary: F(PALETTE.exon, "#B8480A", {
    strokeWidth: 1.2, labelPriority: 7, printFallback: "#9a9a9a",
    marker: "edge",
  }),

  // --- domain architecture ------------------------------------------------- //
  representative_domain: F(PALETTE.domain, "#00456B", {
    text: PALETTE.paper, labelPriority: 9, printFallback: "#6b6b6b",
  }),
  family_superfamily: F(PALETTE.family, "#8F979F", {
    text: PALETTE.ink, labelPriority: 3, printFallback: "#d0d0d0",
  }),
  tm_helix: F(PALETTE.tm, "#9C4F7C", {
    text: PALETTE.paper, labelPriority: 7, printFallback: "#8f8f8f",
  }),

  // --- genomic context ------------------------------------------------------ //
  // The gene a synteny neighbourhood is centred on. It is deliberately the
  // strongest mark in the track: its own slot, never a flanking neighbour.
  synteny_target: F(PALETTE.domain, "#00456B", {
    text: PALETTE.paper, strokeWidth: 1.3, labelPriority: 10,
    printFallback: "#6b6b6b",
  }),

  // --- exploratory candidates ---------------------------------------------- //
  candidate_region: F(PALETTE.candidate, PALETTE.candidateEdge, {
    opacity: 0.55, text: "#8A5008", labelPriority: 4, printFallback: "#e0e0e0",
  }),

  // --- validated events ---------------------------------------------------- //
  // A validated event — FGFR2's IIIb/IIIc cassette exon is the one this project
  // has — must not share the exploratory candidate's style. A reader who cannot
  // tell an established alternative-splicing event from a positional guess has
  // been misled by the figure, so the two get visibly different marks: an opaque
  // fill and a full outline against the candidate's translucent wash.
  validated_event: F(PALETTE.exonAlt, "#8A4A00", {
    text: PALETTE.ink, strokeWidth: 1.1, labelPriority: 9,
    printFallback: "#7a7a7a",
  }),

  // --- boundary classes ---------------------------------------------------- //
  // Marker shapes carry the class a second time, so the classes stay separable
  // in monochrome print and for colour-blind readers.
  boundary_exact: F(PALETTE.exact_domain_edge, PALETTE.exact_domain_edge, {
    marker: "diamond", labelPriority: 9, printFallback: "#4a4a4a",
  }),
  boundary_near: F(PALETTE.near_domain_edge, PALETTE.near_domain_edge, {
    marker: "circle", labelPriority: 8, printFallback: "#6b6b6b",
  }),
  boundary_inside: F(PALETTE.inside_domain, PALETTE.inside_domain, {
    marker: "square", labelPriority: 7, printFallback: "#8a8a8a",
  }),
  boundary_outside: F(PALETTE.outside_annotated_domains, PALETTE.outside_annotated_domains, {
    marker: "triangle", labelPriority: 6, printFallback: "#a8a8a8",
  }),
  boundary_uncertain: F(PALETTE.unavailable_or_uncertain, PALETTE.unavailable_or_uncertain, {
    marker: "cross", labelPriority: 2, printFallback: "#c8c8c8",
  }),

  // --- selection ----------------------------------------------------------- //
  // A selection is an outline, never a translucent slab over several tracks: it
  // must not recolour the feature whose class the reader is trying to judge.
  selected_feature: F("none", PALETTE.ink, {
    strokeWidth: 1.8, labelPriority: 10, printFallback: "#000000", marker: "outline",
  }),

  // --- further annotation tracks ------------------------------------------- //
  // Not in the original nineteen, but present in the viewers: leaving them as
  // literals would have reopened exactly the drift this module closes.
  protein_backbone: F("#526476", "#3E4C5A", {
    strokeWidth: 0.5, labelPriority: 1, printFallback: "#9a9a9a",
  }),
  exon_boundary_tick: F("none", "#6B7A8C", {
    strokeWidth: 0.6, labelPriority: 3, printFallback: "#a8a8a8", marker: "tick",
  }),
  member_signature: F("#6F7F92", "#55636F", {
    text: PALETTE.paper, labelPriority: 3, printFallback: "#b8b8b8",
  }),
  functional_site: F("#B25A00", "#8A4500", {
    marker: "tick", labelPriority: 5, printFallback: "#8f8f8f",
  }),
  disorder_region: F("#B79A57", "#8F7539", {
    opacity: 0.75, labelPriority: 2, printFallback: "#d4d4d4",
  }),

  // --- alignment ----------------------------------------------------------- //
  primary_sequence: F(PALETTE.exonPrimary, "#5E7F9F", {
    text: PALETTE.paper, labelPriority: 9, printFallback: "#9a9a9a",
  }),
  alternative_sequence: F(PALETTE.exon, PALETTE.exonEdge, {
    labelPriority: 5, printFallback: "#c4c4c4",
  }),
  gap: F("none", PALETTE.axis, {
    strokeWidth: 0.6, opacity: 0.7, labelPriority: 1, printFallback: "#e8e8e8",
  }),
  variable_region: F(PALETTE.exonAlt, "#B8801A", {
    opacity: 0.9, labelPriority: 8, printFallback: "#6f6f6f",
  }),
  conserved_region: F(PALETTE.identity, "#0B5A26", {
    text: PALETTE.paper, labelPriority: 6, printFallback: "#7a7a7a",
  }),
};

/** Text roles. Kept apart from features: these are ink, not scientific colour. */
export const TEXT_ROLES = {
  axis: { fill: "#33404F", fontSize: 10 },
  axisEmphasis: { fill: "#33404F", fontSize: 10.5, fontWeight: 600 },
  trackLabel: { fill: PALETTE.muted, fontSize: 9.5 },
  featureLabel: { fill: PALETTE.ink, fontSize: 9.5, fontWeight: 600 },
  onFeatureLabel: { fill: PALETTE.paper, fontSize: 10, fontWeight: 600 },
  candidateLabel: { fill: "#8A5008", fontSize: 10, fontWeight: 600 },
  muted: { fill: PALETTE.muted, fontSize: 9 },
  empty: { fill: PALETTE.muted, fontSize: 10, fontStyle: "italic" },
};

/**
 * Non-feature furniture: paper, grid, rules, hairlines.
 *
 * These carry no scientific meaning, but they still have to be explicit, because a
 * grid line that vanishes in an exported SVG costs the reader the coordinate frame.
 */
export const CHROME = {
  paper: PALETTE.paper,
  grid: PALETTE.grid,
  rule: "#D2D9E6",
  axisLine: PALETTE.axis,
};

/** Repeated domain instances stay separable; index is the 1-based display order. */
export const DOMAIN_INSTANCE_COLOURS = [
  "#0072B2", "#56B4E9", "#3B7EA1", "#7570B3", "#1B7837", "#8C6D31",
];

export const domainInstanceFill = (displayOrder) =>
  DOMAIN_INSTANCE_COLOURS[((displayOrder || 1) - 1) % DOMAIN_INSTANCE_COLOURS.length];

/** Canonical boundary class -> semantic key. One mapping for both renderers. */
export const BOUNDARY_CLASS_STYLE = {
  exact_domain_edge: "boundary_exact",
  near_domain_edge: "boundary_near",
  inside_domain: "boundary_inside",
  outside_annotated_domains: "boundary_outside",
  unavailable_or_uncertain: "boundary_uncertain",
};

/** The legend vocabulary. Both sides must name a class the same way. */
export const BOUNDARY_CLASS_LABEL = {
  exact_domain_edge: "Exact domain edge",
  near_domain_edge: "Near domain edge",
  inside_domain: "Inside domain",
  outside_annotated_domains: "Outside annotated domains",
  unavailable_or_uncertain: "Uncertain / unavailable",
};

/**
 * The style for a semantic key.
 *
 * Throws on an unknown key rather than falling back to a default: a silent
 * default is how a feature ends up painted black in an exported SVG, which is the
 * class of defect this module exists to prevent.
 */
export function featureStyle(key) {
  const style = FEATURE_STYLES[key];
  if (!style) throw new Error(`unknown semantic feature key: ${key}`);
  return style;
}

/**
 * Explicit SVG paint attributes for a semantic key, ready to spread onto a React
 * element. `selected` adds the selection outline without touching the fill, so
 * the feature keeps its scientific colour while selected.
 */
export function featureProps(key, { selected = false, opacity, faint = false } = {}) {
  const s = featureStyle(key);
  const sel = FEATURE_STYLES.selected_feature;
  return {
    fill: s.fill,
    stroke: selected ? sel.stroke : s.stroke,
    strokeWidth: selected ? sel.strokeWidth : s.strokeWidth,
    fillOpacity: opacity ?? (faint ? Math.min(s.opacity, 0.35) : s.opacity),
  };
}

/** Explicit text paint for a role. */
export function textProps(role) {
  const r = TEXT_ROLES[role];
  if (!r) throw new Error(`unknown text role: ${role}`);
  const out = { fill: r.fill, fontSize: r.fontSize };
  if (r.fontWeight) out.fontWeight = r.fontWeight;
  if (r.fontStyle) out.fontStyle = r.fontStyle;
  return out;
}

/** Paint for a boundary record's canonical class. */
export const boundaryProps = (canonicalClass, opts) =>
  featureProps(BOUNDARY_CLASS_STYLE[canonicalClass] || "boundary_uncertain", opts);

export const boundaryStyleKey = (canonicalClass) =>
  BOUNDARY_CLASS_STYLE[canonicalClass] || "boundary_uncertain";
