// Publication figures for the generic Comparative Exon–Domain Boundary Explorer.
//
// Every figure is built from the canonical comparative index that
// scripts/shared_gene_analysis/boundary_dashboard.py publishes. Nothing here
// re-derives which boundaries are comparable, which species were mapped, or how
// confident a mapping is — a second implementation of that logic in the browser
// would be free to disagree with the backend, and a reader could not tell which
// answer the paper was based on.
//
// Rendering reuses the same vector primitives, presets and semantic palette as the
// single-species main figures, so a comparative panel and a per-species panel can
// sit in one figure without a visible style seam. There is no HTML, no CSS and no
// screenshot step in this path: SVG, PDF and PNG all come from the same marks.

import {
  createFigure, preset, PALETTE, textWidth, axisTicks,
} from "./figureSpec.js";
import {
  headerBlock, finalise, speciesLabel, boundaryClassColour, boundaryClassLabel, tsv,
} from "./mainFigures.js";
import { canonClass } from "./boundaryClasses.js";
import { speciesCompare } from "./speciesOrder.js";

// The band within which a boundary counts as sitting on a domain edge. Published by
// the backend as near_edge_band_aa; this is only the fallback for older indices.
export const NEAR_EDGE_BAND_AA = 5;

// Cell fills for the two non-biological matrix states. They are deliberately pale and
// unsaturated so they cannot be mistaken for a boundary class.
export const MATRIX_STATE_FILL = {
  boundary_absent_or_unmapped: "#eef1f4",
  result_pending: "#f6efe0",
};

export function matrixCellFill(state) {
  const cls = canonClass(state);
  // canonClass() maps anything unknown onto "unavailable_or_uncertain", so the
  // non-biological states have to be checked first or they would be painted as a
  // boundary class.
  if (MATRIX_STATE_FILL[state]) return MATRIX_STATE_FILL[state];
  return boundaryClassColour(cls);
}

/** Short species tag ("Gallus gallus" -> "G. gallus") for dense plot annotation. */
export function speciesTag(speciesId, scientificName) {
  const name = scientificName || speciesLabel(speciesId || "");
  const parts = String(name).trim().split(/\s+/);
  if (parts.length < 2) return name;
  return `${parts[0][0]}. ${parts.slice(1).join(" ")}`;
}

const signedLabel = (v) => (v > 0 ? `+${v}` : String(v));

/**
 * A mapping is only drawn as a confirmed cross-species pair when the backend called
 * it supported. Tentative mappings are shown but never connected by a solid line,
 * because a connector reads as "this is the same boundary" and that is exactly the
 * claim a tentative mapping does not support.
 */
export const isSupported = (mappingStatus) => mappingStatus === "supported_comparable"
  || mappingStatus === "high_confidence_comparable";

// --------------------------------------------------------------------------- //
// 1. Comparative matrix
// --------------------------------------------------------------------------- //
/**
 * Species × comparable-boundary-group matrix.
 *
 * @param mode "signed" (default), "absolute" or "class" — the same three toggles the
 *             interactive matrix offers, so the exported figure matches the screen.
 */
export function comparativeMatrixFigureSpec({
  gene, matrix = [], groups = [], selectedGroupId = null, mode = "signed",
  nearEdgeThreshold = NEAR_EDGE_BAND_AA, presetName = "double",
}) {
  const P = preset(presetName);
  const rowLabels = matrix.map((r) => speciesTag(r.species_id, r.scientific_name));
  const gutter = Math.max(70, ...rowLabels.map((s) => textWidth(s, P.font.label) + 10));
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right;
  const nCols = groups.length || 1;
  const colW = Math.max(18, (x1 - x0) / nCols);
  const rowH = Math.max(16, P.font.label + 8);
  const headH = P.font.small * 2 + 8;

  const bodyH = matrix.length * rowH;
  const fig = createFigure({ preset: presetName, height: bodyH + headH + 200 });

  const modeText = mode === "class" ? "boundary class only"
    : mode === "absolute" ? "absolute distance to nearest representative-domain edge"
      : "signed distance to nearest representative-domain edge";
  let y = headerBlock(fig, {
    title: gene,
    subtitle: `Comparative exon–domain boundary matrix · ${modeText}`,
    question: "Which exon–domain boundaries are comparable across species, and does "
      + "each species place that boundary at the same distance from a domain edge?",
  });
  y += 6;

  // Column headers: group id plus a mapping-status marker, so a tentative column is
  // identifiable in print without the interactive tooltip.
  const top = y + headH;
  groups.forEach((g, i) => {
    const cx = x0 + i * colW + colW / 2;
    const sel = selectedGroupId && g.comparable_boundary_group_id === selectedGroupId;
    fig.text(cx, y + P.font.small, String(g.comparable_boundary_group_id || "").replace(/^CBG/, ""),
      { size: "small", anchor: "middle", fill: PALETTE.ink,
        weight: sel ? "bold" : "normal" });
    if (!isSupported(g.mapping_status)) {
      fig.text(cx, y + P.font.small * 2 + 2, "t",
        { size: "small", anchor: "middle", fill: PALETTE.muted, italic: true });
    }
    if (sel) {
      fig.rect(x0 + i * colW, top - 2, colW, bodyH + 4,
        { fill: "none", stroke: PALETTE.ink, lw: P.lw.outline });
    }
  });

  matrix.forEach((row, ri) => {
    const ry = top + ri * rowH;
    fig.text(x0 - 6, ry + rowH / 2 + P.font.label * 0.35,
      speciesTag(row.species_id, row.scientific_name),
      { size: "label", anchor: "end", italic: true });
    (row.cells || []).forEach((c, ci) => {
      const cx = x0 + ci * colW;
      fig.rect(cx, ry, colW - 1, rowH - 1,
        { fill: matrixCellFill(c.state), stroke: PALETTE.paper, lw: P.lw.thin });
      if (mode === "class" || !c.observed) return;
      const v = mode === "absolute" ? c.absolute_distance : c.signed_distance;
      if (v == null) return;
      // Label colour follows cell luminance so the number stays legible on both the
      // dark edge classes and the pale outside-domain class.
      const dark = ["exact_domain_edge", "near_domain_edge", "inside_domain"]
        .includes(canonClass(c.state));
      fig.text(cx + (colW - 1) / 2, ry + rowH / 2 + P.font.small * 0.35,
        mode === "absolute" ? String(v) : signedLabel(v),
        { size: "small", anchor: "middle", fill: dark ? PALETTE.paper : PALETTE.ink });
    });
  });

  y = top + bodyH + 10;
  fig.text(x0, y, `Columns: comparable-boundary groups (CBG) · "t" marks a tentative `
    + `alignment mapping · near-edge threshold ±${nearEdgeThreshold} aa`,
    { size: "small", fill: PALETTE.muted });
  y += P.font.small + 6;

  const classes = [...new Set((matrix.flatMap((r) => r.cells || []))
    .filter((c) => c.observed).map((c) => canonClass(c.state)))];
  y = fig.legend(x0, y, [
    ...classes.map((c) => [boundaryClassColour(c), boundaryClassLabel(c)]),
    [MATRIX_STATE_FILL.boundary_absent_or_unmapped, "No comparable boundary mapped"],
    [MATRIX_STATE_FILL.result_pending, "Species analysis pending"],
  ]);
  return finalise(fig, y + 2);
}

// --------------------------------------------------------------------------- //
// 2. Paired (dumbbell) signed-distance plot
// --------------------------------------------------------------------------- //
/**
 * One row per comparable-boundary group; one dot per species observation.
 *
 * With the two species of the current real datasets this is a paired dot / dumbbell
 * plot: the two raw observations and the gap between them are the whole result, and a
 * boxplot of n = 2 would draw a distribution that does not exist.
 *
 * The layout scales without a rewrite. ``summary`` selects the extra annotation:
 *   - "pair"  (2 species)    connector between the two observations
 *   - "range" (3–4 species)  range bar plus a median tick, raw points kept
 *   - "box"   (5+ species)   quartile box plus median, raw points kept
 * Raw points are never replaced by a summary, at any species count.
 */
export function pairedSignedDistanceFigureSpec({
  gene, groups = [], stats = [], selectedGroupId = null,
  nearEdgeThreshold = NEAR_EDGE_BAND_AA, summary = "auto", presetName = "double",
}) {
  const P = preset(presetName);
  const statById = new Map(stats.map((s) => [s.comparable_boundary_group_id, s]));

  const rows = groups.map((g) => {
    const obs = (g.per_species_native_positions || [])
      .filter((o) => o.signed_distance != null)
      .slice()
      .sort((a, b) => speciesCompare(a.species_id, b.species_id));
    return {
      id: g.comparable_boundary_group_id,
      label: String(g.comparable_boundary_group_id || "").replace(/^CBG/, "CBG "),
      mappingStatus: g.mapping_status,
      supported: isSupported(g.mapping_status),
      msaColumn: g.msa_column,
      obs,
      stat: statById.get(g.comparable_boundary_group_id) || {},
    };
  }).filter((r) => r.obs.length);

  const nSpeciesMax = Math.max(1, ...rows.map((r) => r.obs.length));
  const mode = summary !== "auto" ? summary
    : nSpeciesMax <= 2 ? "pair" : nSpeciesMax <= 4 ? "range" : "box";

  // Species are identified by a fixed lane inside each row (first species on top),
  // not by a label next to every dot: with sixteen groups and two species those
  // thirty-two labels collide with the markers they are meant to explain, and the
  // near-zero differences in real data put both dots in the same place.
  const speciesOrder = [...new Set(rows.flatMap((r) => r.obs.map((o) => o.species_id)))]
    .sort((a, b) => String(a).localeCompare(String(b)));
  const speciesName = new Map();
  for (const r of rows) {
    for (const o of r.obs) speciesName.set(o.species_id, o.scientific_name);
  }

  const gutter = Math.max(56, ...rows.map((r) => textWidth(r.label, P.font.label) + 10));
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right - 18;
  // Enough height for one clearly separated lane per species.
  const rowH = Math.max(14, P.font.label + 7, speciesOrder.length * 6 + 6);

  const all = rows.flatMap((r) => r.obs.map((o) => Number(o.signed_distance)));
  const maxAbs = Math.max(10, ...all.map((v) => Math.abs(v)));
  const lim = Math.ceil(maxAbs / 10) * 10;
  const scale = (v) => x0 + ((v + lim) / (2 * lim)) * (x1 - x0);

  const bodyH = rows.length * rowH;
  const fig = createFigure({ preset: presetName, height: bodyH + 210 });

  let y = headerBlock(fig, {
    title: gene,
    subtitle: "Signed distance from each comparable exon–domain boundary to the nearest "
      + "representative-domain edge, per species",
    question: "Do the species place the same comparable boundary at the same distance, "
      + "and on the same side, of the nearest domain edge?",
  });
  y += 8;

  const top = y;
  const bottom = top + bodyH;

  // Near-edge band and the domain-edge zero line, drawn behind the data.
  fig.rect(scale(-nearEdgeThreshold), top - 3,
    scale(nearEdgeThreshold) - scale(-nearEdgeThreshold), bodyH + 6,
    { fill: PALETTE.grid, stroke: "none", opacity: 0.9 });
  fig.text(scale(0), top - 6, `±${nearEdgeThreshold} aa`,
    { size: "small", anchor: "middle", fill: PALETTE.muted });
  fig.line(scale(0), top - 3, scale(0), bottom + 3,
    { stroke: PALETTE.ink, lw: P.lw.rule });

  rows.forEach((r, ri) => {
    const cy = top + ri * rowH + rowH / 2;
    const sel = selectedGroupId && r.id === selectedGroupId;
    if (sel) {
      fig.rect(P.margin.left, cy - rowH / 2, x1 - P.margin.left + 14, rowH,
        { fill: PALETTE.grid, stroke: "none", opacity: 0.55 });
    }
    fig.text(x0 - 6, cy + P.font.label * 0.35, r.label,
      { size: "label", anchor: "end", fill: sel ? PALETTE.ink : PALETTE.muted,
        weight: sel ? "bold" : "normal" });

    // Each observation's plotted point, so the connector joins the markers that are
    // actually drawn. A horizontal-only connector vanishes whenever the two species
    // agree exactly, which is common in real data and is precisely the case a reader
    // most needs to see as a connected pair.
    const laneGap = Math.min(6, rowH / Math.max(1, speciesOrder.length));
    const pointOf = (o) => {
      const li = Math.max(0, speciesOrder.indexOf(o.species_id));
      const dy = speciesOrder.length > 1
        ? (li - (speciesOrder.length - 1) / 2) * laneGap : 0;
      return [scale(Number(o.signed_distance)), cy + dy];
    };
    const xs = r.obs.map((o) => scale(Number(o.signed_distance)));
    const lo = Math.min(...xs);
    const hi = Math.max(...xs);

    if (mode === "pair" && r.obs.length === 2) {
      const [ax, ay] = pointOf(r.obs[0]);
      const [bx, by] = pointOf(r.obs[1]);
      if (r.supported) {
        // A solid connector asserts "the same boundary in both species". It is drawn
        // only for supported mappings.
        fig.line(ax, ay, bx, by,
          { stroke: PALETTE.boundary, lw: P.lw.rule, opacity: 0.8 });
      } else {
        // Tentative: the positions are close, but equivalence is not established, so
        // the pair is hinted at rather than asserted.
        fig.line(ax, ay, bx, by,
          { stroke: PALETTE.muted, lw: P.lw.thin, dash: "1,2", opacity: 0.9 });
      }
    } else if (mode === "range" || mode === "box") {
      const sorted = [...r.obs.map((o) => Number(o.signed_distance))].sort((a, b) => a - b);
      if (mode === "box" && sorted.length >= 5) {
        const q = (p) => {
          const idx = (sorted.length - 1) * p;
          const f = Math.floor(idx);
          return sorted[f] + (sorted[Math.min(f + 1, sorted.length - 1)] - sorted[f]) * (idx - f);
        };
        fig.rect(scale(q(0.25)), cy - rowH * 0.26,
          scale(q(0.75)) - scale(q(0.25)), rowH * 0.52,
          { fill: PALETTE.grid, stroke: PALETTE.axis, lw: P.lw.thin, opacity: 0.9 });
      } else {
        fig.line(lo, cy, hi, cy, { stroke: PALETTE.axis, lw: P.lw.thin, opacity: 0.9 });
      }
      const med = r.stat.median_signed_distance;
      if (med != null) {
        fig.line(scale(med), cy - rowH * 0.28, scale(med), cy + rowH * 0.28,
          { stroke: PALETTE.ink, lw: P.lw.rule });
      }
    }

    // Raw observations. Class colour carries the boundary class; the marker outline
    // carries the domain edge the distance was measured to; a small species tag
    // keeps the two dots of a pair distinguishable in print.
    r.obs.forEach((o) => {
      const cx = scale(Number(o.signed_distance));
      // Two species reporting the same distance would overplot exactly, so each
      // species keeps its own lane. The lane is the species identity, which is why it
      // is derived from the run-wide species order rather than from this row.
      const li = Math.max(0, speciesOrder.indexOf(o.species_id));
      const dy = speciesOrder.length > 1
        ? (li - (speciesOrder.length - 1) / 2) * Math.min(6, rowH / speciesOrder.length)
        : 0;
      const colour = boundaryClassColour(canonClass(o.boundary_class));
      const rad = P.marker * (sel ? 1.3 : 1);
      if (o.nearest_edge === "start") {
        fig.circle(cx, cy + dy, rad,
          { fill: PALETTE.paper, stroke: colour, lw: P.lw.outline });
      } else {
        fig.circle(cx, cy + dy, rad,
          { fill: colour, stroke: sel ? PALETTE.ink : "none", lw: P.lw.thin });
      }
    });
  });

  y = bottom + 10;
  const { major } = axisTicks(-lim, lim, 8);
  fig.line(x0, y, x1, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of major) {
    if (t < -lim || t > lim) continue;
    fig.line(scale(t), y, scale(t), y + 3.5, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(scale(t), y + 4 + P.font.tick, signedLabel(t),
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  }
  y += 4 + P.font.tick * 2 + 4;
  fig.text((x0 + x1) / 2, y,
    "Signed distance to nearest representative-domain edge (aa) · 0 = domain edge",
    { size: "label", anchor: "middle" });
  y += P.font.label + 6;

  const classes = [...new Set(rows.flatMap((r) => r.obs.map(
    (o) => canonClass(o.boundary_class))))];
  y = fig.legend(x0, y, classes.map((c) => [boundaryClassColour(c), boundaryClassLabel(c)]));

  // Which lane belongs to which species, with the full scientific name. Without
  // this the vertical offset is a decoration rather than an identity, and the
  // reader cannot tell which point is which organism.
  if (speciesOrder.length > 1) {
    y += 3;
    fig.text(x0, y + P.font.small, "Marker lane per species",
      { size: "small", fill: PALETTE.muted });
    let lx = x0 + 108;
    speciesOrder.forEach((sid, i) => {
      const lane = i === 0 ? "upper"
        : i === speciesOrder.length - 1 ? "lower" : `lane ${i + 1}`;
      const name = speciesName.get(sid) || sid;
      fig.circle(lx + 3, y + P.font.small - 3, P.marker,
        { fill: PALETTE.paper, stroke: PALETTE.ink, lw: P.lw.thin });
      fig.text(lx + 11, y + P.font.small, `${lane} · `,
        { size: "small", fill: PALETTE.muted });
      fig.text(lx + 11 + (lane.length + 3) * P.font.small * 0.5, y + P.font.small,
        name, { size: "small", style: "italic", fill: PALETTE.ink });
      lx += 32 + (lane.length + 3 + name.length) * P.font.small * 0.52;
    });
    y += P.font.small + 4;
  }
  const connectorNote = mode === "pair"
    ? "Solid connector: supported cross-species mapping · dotted connector: tentative "
      + "mapping, equivalence not established"
    : "Bar: observed range · vertical tick: median · raw observations always shown";
  fig.text(x0, y + P.font.small + 1, connectorNote,
    { size: "small", fill: PALETTE.muted });
  y += P.font.small + 3;
  fig.text(x0, y + P.font.small + 1,
    "Open marker: distance to a domain start edge · filled marker: to a domain end edge "
    + "· negative: upstream of that edge",
    { size: "small", fill: PALETTE.muted });
  return finalise(fig, y + P.font.small + 4);
}

// --------------------------------------------------------------------------- //
// 3. Boundary-position consistency summary
// --------------------------------------------------------------------------- //
/**
 * Per-group consistency, drawn as the cross-species difference with the species
 * coverage and mapping status beside it.
 *
 * For two species the bar is literally the gap between the two raw observations, which
 * is why it is labelled "cross-species difference" and not "variance" or "spread": at
 * n = 2 those words would dress up a single subtraction as a distribution.
 */
export function consistencySummaryFigureSpec({
  gene, stats = [], groups = [], selectedGroupId = null, presetName = "single",
}) {
  const P = preset(presetName);
  const groupById = new Map(groups.map((g) => [g.comparable_boundary_group_id, g]));
  const rows = stats
    .filter((s) => s.cross_species_difference != null)
    .map((s) => ({
      id: s.comparable_boundary_group_id,
      label: String(s.comparable_boundary_group_id || "").replace(/^CBG/, "CBG "),
      diff: Number(s.cross_species_difference),
      coverage: s.species_coverage ?? s.mapping_coverage,
      cls: canonClass(s.dominant_class),
      classesDiffer: Boolean(s.classes_differ),
      supported: isSupported(
        s.mapping_status || groupById.get(s.comparable_boundary_group_id)?.mapping_status),
      nRaw: (s.raw_signed_distances || []).length,
    }));

  const gutter = Math.max(56, ...rows.map((r) => textWidth(r.label, P.font.label) + 10));
  const x0 = P.margin.left + gutter;
  // The value label sits outside the bar, so the plot area has to give way to the
  // longest one — otherwise the widest bar's annotation is clipped off the page.
  const valueLabel = (r) => `${r.diff} aa${r.classesDiffer ? " · classes differ" : ""}`;
  const labelRoom = Math.max(24, ...rows.map((r) => textWidth(valueLabel(r), P.font.small) + 6));
  const x1 = P.widthPt - P.margin.right - labelRoom;
  const rowH = Math.max(13, P.font.label + 6);
  const maxDiff = Math.max(1, ...rows.map((r) => r.diff));
  const lim = Math.ceil(maxDiff / 5) * 5;
  const scale = (v) => x0 + (v / lim) * (x1 - x0);

  const bodyH = rows.length * rowH;
  const fig = createFigure({ preset: presetName, height: bodyH + 200 });

  const nRaw = Math.max(0, ...rows.map((r) => r.nRaw));
  let y = headerBlock(fig, {
    title: gene,
    subtitle: nRaw <= 2
      ? "Boundary-position consistency · difference between the two species' raw "
        + "signed distances"
      : "Boundary-position consistency · range of the species' signed distances",
    question: "For each comparable boundary, how much do the species disagree about "
      + "its distance to the nearest domain edge?",
  });
  y += 8;
  const top = y;

  rows.forEach((r, ri) => {
    const ry = top + ri * rowH;
    const sel = selectedGroupId && r.id === selectedGroupId;
    fig.text(x0 - 6, ry + rowH / 2 + P.font.label * 0.35, r.label,
      { size: "label", anchor: "end", fill: sel ? PALETTE.ink : PALETTE.muted,
        weight: sel ? "bold" : "normal" });
    const w = Math.max(0.6, scale(r.diff) - x0);
    // A tentative mapping gets an outlined bar: the number is real, but the claim
    // that the two observations describe the same boundary is not established.
    if (r.supported) {
      fig.rect(x0, ry + rowH * 0.2, w, rowH * 0.6,
        { fill: boundaryClassColour(r.cls), stroke: "none" });
    } else {
      fig.rect(x0, ry + rowH * 0.2, w, rowH * 0.6,
        { fill: PALETTE.paper, stroke: boundaryClassColour(r.cls), lw: P.lw.outline });
    }
    fig.text(x0 + w + 3, ry + rowH / 2 + P.font.small * 0.35, valueLabel(r),
      { size: "small", fill: PALETTE.muted });
  });

  y = top + bodyH + 8;
  const { major } = axisTicks(0, lim, 6);
  fig.line(x0, y, x1, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of major) {
    if (t < 0 || t > lim) continue;
    fig.line(scale(t), y, scale(t), y + 3.5, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(scale(t), y + 4 + P.font.tick, String(t),
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  }
  y += 4 + P.font.tick * 2 + 4;
  fig.text((x0 + x1) / 2, y, "Cross-species difference in signed distance (aa)",
    { size: "label", anchor: "middle" });
  y += P.font.label + 6;
  // Two short lines rather than one long one: a single line at this width runs off the
  // page in the narrow single-column preset.
  fig.text(P.margin.left, y,
    "Bar colour: dominant boundary class · outlined bar: tentative mapping, "
    + "comparability not established.",
    { size: "small", fill: PALETTE.muted });
  y += P.font.small + 3;
  fig.text(P.margin.left, y,
    "A small difference is evidence of a consistent boundary position, not of "
    + "evolutionary conservation.",
    { size: "small", fill: PALETTE.muted });
  return finalise(fig, y + P.font.small + 4);
}

// --------------------------------------------------------------------------- //
// 4. Comparative local architecture for the selected group
// --------------------------------------------------------------------------- //
/**
 * The selected comparable boundary in each species' own local domain context.
 *
 * Each species gets one track drawn on its own native amino-acid axis but with the
 * tracks aligned on the boundary itself, so the domain edges either line up or
 * visibly do not. A shared MSA column axis is used for the tick labels when the
 * backend mapped every observation to the same column; otherwise the panel says so
 * rather than implying a common coordinate system that does not hold.
 */
export function comparativeArchitectureFigureSpec({
  gene, group, models = [], windowAa = 140, presetName = "double",
}) {
  const P = preset(presetName);
  const obs = (group?.per_species_native_positions || []);
  const modelBySpecies = new Map(models.map((m) => [m.species_id, m]));

  const trackH = 16;
  const laneH = trackH + P.font.small + 16;
  const fig = createFigure({ preset: presetName, height: obs.length * laneH + 210 });

  const supported = isSupported(group?.mapping_status);
  const columns = [...new Set(obs.map((o) => o.msa_column).filter((c) => c != null))];
  const sharedColumn = columns.length === 1 ? columns[0] : null;

  let y = headerBlock(fig, {
    title: gene,
    subtitle: `Comparative local architecture around ${group?.comparable_boundary_group_id || ""}`
      + (sharedColumn != null ? ` · alignment column ${sharedColumn}` : ""),
    question: "Does the same comparable boundary sit in the same local domain context "
      + "in every species?",
  });
  y += 4;

  if (!supported) {
    fig.text(P.margin.left, y + P.font.small,
      "Mapping is tentative: the tracks are aligned on each species' own boundary, "
      + "which does not establish that they are the same junction.",
      { size: "small", fill: PALETTE.muted, italic: true });
    y += P.font.small + 6;
  } else if (sharedColumn == null) {
    fig.text(P.margin.left, y + P.font.small,
      "Observations map to different alignment columns; tracks are aligned on each "
      + "species' own boundary rather than on a shared coordinate axis.",
      { size: "small", fill: PALETTE.muted, italic: true });
    y += P.font.small + 6;
  }
  y += 6;

  const labels = obs.map((o) => speciesTag(o.species_id, o.scientific_name));
  const gutter = Math.max(78, ...labels.map((s) => textWidth(s, P.font.label) + 12));
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right;
  const half = windowAa / 2;
  // Offset coordinates relative to the boundary: 0 is the boundary in every lane.
  const scale = (delta) => x0 + ((delta + half) / windowAa) * (x1 - x0);

  obs.slice()
    .sort((a, b) => speciesCompare(a.species_id, b.species_id))
    .forEach((o, i) => {
      const ly = y + i * laneH;
      const pos = Number(o.native_position);
      fig.text(x0 - 8, ly + trackH / 2 + P.font.label * 0.35,
        speciesTag(o.species_id, o.scientific_name),
        { size: "label", anchor: "end", italic: true });

      // Protein backbone for the visible window.
      fig.rect(x0, ly + trackH / 2 - 1.2, x1 - x0, 2.4,
        { fill: PALETTE.grid, stroke: "none" });

      const model = modelBySpecies.get(o.species_id);
      const domains = (model?.representative_domains || []);
      for (const d of domains) {
        const ds = Number(d.start) - pos;
        const de = Number(d.end) - pos;
        if (de < -half || ds > half) continue;
        const cs = Math.max(ds, -half);
        const ce = Math.min(de, half);
        fig.rect(scale(cs), ly, scale(ce) - scale(cs), trackH,
          { fill: PALETTE.domain, stroke: PALETTE.ink, lw: P.lw.thin, opacity: 0.85, rx: 1.5 });
        const w = scale(ce) - scale(cs);
        const label = d.short_label || d.label || "";
        if (w > textWidth(label, P.font.small) + 6) {
          fig.text(scale(cs) + w / 2, ly + trackH / 2 + P.font.small * 0.35, label,
            { size: "small", anchor: "middle", fill: PALETTE.paper });
        }
      }

      // The boundary itself, and the domain edge the distance was measured to.
      const colour = boundaryClassColour(canonClass(o.boundary_class));
      fig.line(scale(0), ly - 3, scale(0), ly + trackH + 3,
        { stroke: colour, lw: P.lw.rule });
      if (o.nearest_edge_position != null) {
        const edgeDelta = Number(o.nearest_edge_position) - pos;
        if (edgeDelta >= -half && edgeDelta <= half) {
          fig.line(scale(edgeDelta), ly - 1, scale(edgeDelta), ly + trackH + 1,
            { stroke: PALETTE.ink, lw: P.lw.thin, dash: "2,2" });
        }
      }
      fig.text(x0, ly + trackH + P.font.small + 4,
        `${o.exon_transition || ""} · native aa ${o.native_position}`
        + (o.msa_column != null ? ` · column ${o.msa_column}` : "")
        + ` · ${signedLabel(Number(o.signed_distance))} aa to `
        + `${o.nearest_edge || "?"} edge of ${o.nearest_domain_label || "no annotated domain"}`,
        { size: "small", fill: PALETTE.muted });
    });

  y += obs.length * laneH + 4;
  const { major } = axisTicks(-half, half, 8);
  fig.line(x0, y, x1, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of major) {
    if (t < -half || t > half) continue;
    fig.line(scale(t), y, scale(t), y + 3.5, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(scale(t), y + 4 + P.font.tick, signedLabel(t),
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  }
  y += 4 + P.font.tick * 2 + 4;
  fig.text((x0 + x1) / 2, y,
    "Amino-acid offset from the comparable boundary (aa) · 0 = boundary in each species",
    { size: "label", anchor: "middle" });
  y += P.font.label + 6;
  fig.text(x0, y,
    "Blue blocks: representative domains · coloured vertical line: the exon boundary "
    + "· dashed line: the domain edge the distance was measured to",
    { size: "small", fill: PALETTE.muted });
  return finalise(fig, y + P.font.small + 4);
}

// --------------------------------------------------------------------------- //
// Source tables
// --------------------------------------------------------------------------- //
const LONG_COLUMNS = [
  "comparable_boundary_group_id", "mapping_method", "mapping_status",
  "mapping_confidence", "msa_column", "species_id", "scientific_name",
  "taxonomic_group", "protein_id", "transcript_id", "boundary_id", "exon_transition",
  "native_position", "nearest_domain_instance_id", "nearest_domain_accession",
  "nearest_domain_label", "nearest_domain_start", "nearest_domain_end",
  "nearest_edge", "nearest_edge_position", "signed_distance", "absolute_distance",
  "boundary_class", "domain_annotation_available",
];

/** Long-format table: one row per species observation of a comparable boundary. */
export function comparativeLongTsv(groups = []) {
  const rows = groups.flatMap((g) => (g.per_species_native_positions || []).map((o) => ({
    ...o,
    comparable_boundary_group_id: g.comparable_boundary_group_id,
    mapping_method: o.mapping_method || g.mapping_method,
    mapping_status: o.mapping_status || g.mapping_status,
    mapping_confidence: o.mapping_confidence ?? g.confidence,
    msa_column: o.msa_column ?? g.msa_column,
  })));
  return tsv(rows, LONG_COLUMNS);
}

/** Wide matrix table, in the value mode currently shown on screen. */
export function comparativeMatrixTsv(matrix = [], groups = [], mode = "signed") {
  const ids = groups.map((g) => g.comparable_boundary_group_id);
  const rows = matrix.map((r) => {
    const out = { species_id: r.species_id, scientific_name: r.scientific_name };
    (r.cells || []).forEach((c) => {
      const key = c.comparable_boundary_group_id;
      if (!c.observed) {
        out[key] = c.state;
      } else if (mode === "class") {
        out[key] = c.state;
      } else if (mode === "absolute") {
        out[key] = c.absolute_distance;
      } else {
        out[key] = c.signed_distance;
      }
    });
    return out;
  });
  return tsv(rows, ["species_id", "scientific_name", ...ids]);
}

const MAPPING_COLUMNS = [
  "comparable_boundary_group_id", "mapping_method", "mapping_status",
  "mapping_confidence", "msa_column", "n_species_mapped", "species_ids",
  "species_coverage", "cross_species_difference", "distance_range_min",
  "distance_range_max", "dominant_class", "classes_differ",
  "domain_annotation_available_in_all", "primary_statistic",
];

/** How each comparable-boundary group was established, and how consistent it is. */
export function comparableMappingTsv(groups = [], stats = []) {
  const statById = new Map(stats.map((s) => [s.comparable_boundary_group_id, s]));
  const rows = groups.map((g) => {
    const s = statById.get(g.comparable_boundary_group_id) || {};
    const range = s.distance_range || [];
    const members = g.per_species_native_positions || [];
    return {
      comparable_boundary_group_id: g.comparable_boundary_group_id,
      mapping_method: g.mapping_method,
      mapping_status: g.mapping_status,
      mapping_confidence: g.confidence,
      msa_column: g.msa_column,
      n_species_mapped: members.length,
      species_ids: members.map((m) => m.species_id).join(","),
      species_coverage: s.species_coverage ?? s.mapping_coverage,
      cross_species_difference: s.cross_species_difference,
      distance_range_min: range[0],
      distance_range_max: range[1],
      dominant_class: s.dominant_class,
      classes_differ: s.classes_differ,
      domain_annotation_available_in_all: s.domain_annotation_available_in_all,
      primary_statistic: s.primary_statistic,
    };
  });
  return tsv(rows, MAPPING_COLUMNS);
}
