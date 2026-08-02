// Publication figure builders for the within-species protein isoform alignment.
//
// Three figures share one data model and one visual language:
//
//   1. alignmentOverviewFigureSpec  — the complete alignment at column
//      resolution: one row per protein model, difference and gap structure,
//      conservation, variable blocks and the candidate intervals.
//   2. wrappedAlignmentFigureSpecs  — the same alignment at residue resolution,
//      wrapped Jalview-style into blocks of 60–100 columns across as many pages
//      as it takes.
//   3. candidateAlignmentFigureSpec — one candidate interval at residue
//      resolution, with exon association and both coordinate systems.
//
// Every figure is assembled with figureSpec.js instead of by concatenating SVG
// strings, so each one renders both to a standalone SVG and to a real vector PDF
// whose text stays selectable. Nothing here rasterises anything.
//
// Colour carries meaning and nothing else: the primary sequence is dark, a
// matching residue is neutral, a residue differing from the primary is
// highlighted, a gap is light, and a candidate interval is a pale amber overlay.
// No isoform ever receives a categorical colour of its own — a colour per model
// would encode identity, which the row label already carries, and would leave no
// colour free for the alignment content itself.

import {
  PALETTE, axisTicks, createFigure, placeBlockLabels, preset, textWidth,
  // The explicit extension keeps this module importable by plain Node, so the
  // figures can be rendered and validated without a browser or bundler.
} from "./figureSpec.js";
import { speciesLabel, tsv } from "./mainFigures.js";

export { tsv };

/**
 * Alignment colour semantics, shared with the interactive alignment view so the
 * two cannot drift apart. Derived from the figure palette wherever a matching
 * role already exists there.
 */
export const ALN_COLOURS = {
  paper: PALETTE.paper,
  ink: PALETTE.ink,
  muted: PALETTE.muted,
  axis: PALETTE.axis,
  grid: PALETTE.grid,
  // primary protein: the darkest mark in the figure
  primary: "#39536E",
  // an alternative model's aligned residues
  match: PALETTE.exon,
  matchEdge: PALETTE.exonEdge,
  // residue cells at letter resolution
  primaryCell: "#DCE4EC",
  matchCell: "#F1F4F7",
  // the same cells inside a candidate interval, pre-blended with the candidate
  // amber: an opaque colour renders identically in the SVG and in the PDF, while
  // a translucent overlay would not, because a minimal PDF has no soft mask.
  primaryBandCell: "#E7DCC8",
  matchBandCell: "#FBF1DE",
  // a residue that differs from the primary at the same alignment column
  diff: "#D55E00",
  diffCell: "#F6DCC9",
  diffInk: "#6E2F00",
  // gap in this model, and the edge of a contiguous gap block
  gap: "#EEF1F5",
  gapEdge: "#B9C0C8",
  gapDensity: PALETTE.muted,
  // per-column conservation across the shown models
  conservation: PALETTE.identity,
  // structural annotation: contiguous variable blocks
  block: PALETTE.boundary,
  exon: PALETTE.exon,
  exonEdge: PALETTE.exonEdge,
  // exploratory candidate interval
  candidate: PALETTE.candidate,
  candidateEdge: PALETTE.candidateEdge,
  // the primary protein's own coordinate system
  reference: PALETTE.domain,
  trackBase: "#FAFBFC",
};

const GAP = "-";
const isGap = (ch) => !ch || ch === GAP;

// --------------------------------------------------------------------------- //
// Per-column statistics
// --------------------------------------------------------------------------- //

/**
 * Per-column state across the shown rows.
 *
 * `conservation` is the fraction of *all* shown rows carrying the column's most
 * frequent residue, so a column that is well conserved among the models which
 * have residues there still scores low when half of them are gapped. That is the
 * property an isoform alignment needs: a deletion is not conservation.
 *
 * `diffFraction` only counts alternative models that have a residue where the
 * primary also has one, so it measures substitution rather than indel structure;
 * `gapFraction` measures that separately.
 */
export function alignmentProfile(rows, nCols, primarySeq = "") {
  const nRows = Math.max(1, rows.length);
  const variable = new Uint8Array(nCols);
  const gapped = new Uint8Array(nCols);
  const substitution = new Uint8Array(nCols);
  const conserved = new Uint8Array(nCols);
  const conservation = new Float64Array(nCols);
  const gapFraction = new Float64Array(nCols);
  const diffFraction = new Float64Array(nCols);

  for (let c = 0; c < nCols; c++) {
    const counts = new Map();
    let gaps = 0;
    let compared = 0;
    let differing = 0;
    const pch = primarySeq[c] || GAP;
    for (const r of rows) {
      const ch = (r.seq || "")[c] || GAP;
      if (isGap(ch)) gaps += 1;
      else counts.set(ch, (counts.get(ch) || 0) + 1);
      if (r.is_primary || isGap(pch) || isGap(ch)) continue;
      compared += 1;
      if (ch !== pch) differing += 1;
    }
    let modal = 0;
    for (const v of counts.values()) if (v > modal) modal = v;
    const distinct = counts.size + (gaps > 0 ? 1 : 0);
    variable[c] = distinct > 1 ? 1 : 0;
    gapped[c] = gaps > 0 ? 1 : 0;
    substitution[c] = counts.size > 1 ? 1 : 0;
    conserved[c] = counts.size === 1 && gaps === 0 ? 1 : 0;
    conservation[c] = modal / nRows;
    gapFraction[c] = gaps / nRows;
    diffFraction[c] = compared ? differing / compared : 0;
  }
  return { variable, gapped, substitution, conserved, conservation, gapFraction,
    diffFraction, nRows: rows.length, nCols };
}

/** Column statistics for the interactive alignment, from the same computation. */
export function columnStats(rows, nCols) {
  const primary = rows.find((r) => r.is_primary) || rows[0];
  return alignmentProfile(rows, nCols, primary?.seq || "");
}

/** Contiguous runs where `flags[c]` is set, at least `minLen` columns long. */
export function runsOf(flags, minLen = 1) {
  const out = [];
  let s = -1;
  for (let c = 0; c < flags.length; c++) {
    if (flags[c]) { if (s < 0) s = c; }
    else if (s >= 0) { if (c - s >= minLen) out.push([s, c - 1]); s = -1; }
  }
  if (s >= 0 && flags.length - s >= minLen) out.push([s, flags.length - 1]);
  return out;
}

/** Identical and compared residue pairs against the primary over [c0, c1]. */
function identityCounts(seq, primary, c0 = 0, c1 = null) {
  const end = c1 == null ? (primary || "").length - 1 : c1;
  let same = 0;
  let compared = 0;
  for (let c = c0; c <= end; c++) {
    const p = (primary || "")[c] || GAP;
    const s = (seq || "")[c] || GAP;
    if (isGap(p) || isGap(s)) continue;
    compared += 1;
    if (p === s) same += 1;
  }
  return { same, compared };
}

/** Percentage identity to the primary over the columns where both have a residue. */
export function identityPct(seq, primary, c0 = 0, c1 = null) {
  const { same, compared } = identityCounts(seq, primary, c0, c1);
  return compared ? Math.round((same / compared) * 100) : null;
}

/**
 * Identity as a figure label.
 *
 * A single substitution in a 370-residue protein is 99.7% identity, which rounds
 * to "100%" and would then read as "identical". One decimal is kept just below
 * 100% so a real difference is never reported as none.
 */
export function identityLabel(seq, primary, c0 = 0, c1 = null) {
  const { same, compared } = identityCounts(seq, primary, c0, c1);
  if (!compared) return "no aligned residues";
  const pct = (same / compared) * 100;
  if (same === compared) return "100%";
  return `${pct >= 99.5 ? pct.toFixed(1) : Math.round(pct)}%`;
}

/** Map a 1-based primary-protein residue to its 0-based alignment column. */
export function aaToColumn(primaryAligned, aa) {
  if (!primaryAligned || aa == null) return null;
  let residue = 0;
  for (let c = 0; c < primaryAligned.length; c++) {
    if (!isGap(primaryAligned[c])) {
      residue += 1;
      if (residue === aa) return c;
    }
  }
  return null;
}

/** Column → 1-based primary residue index, for every column the primary occupies. */
export function primaryResidueByColumn(primaryAligned) {
  const map = new Map();
  let residue = 0;
  for (let c = 0; c < (primaryAligned || "").length; c++) {
    if (!isGap(primaryAligned[c])) {
      residue += 1;
      map.set(c, residue);
    }
  }
  return map;
}

/** Which alternative models differ from the primary anywhere inside [c0, c1]. */
export function affectedProteins(rows, primarySeq, c0, c1) {
  const out = new Set();
  for (const r of rows) {
    if (r.is_primary) continue;
    for (let c = c0; c <= c1; c++) {
      if (((r.seq || "")[c] || GAP) !== ((primarySeq || "")[c] || GAP)) {
        out.add(r.protein_id);
        break;
      }
    }
  }
  return out;
}

/** Row ordering used by every alignment view: primary first, then accession. */
export function orderRows(rows) {
  return [...(rows || [])].sort(
    (a, b) => (b.is_primary ? 1 : 0) - (a.is_primary ? 1 : 0)
      || String(a.protein_id).localeCompare(String(b.protein_id)));
}

/** Curation state of a protein model, defaulting to "predicted". */
export const curationMark = (row) =>
  (row?.curation_status === "curated" ? "curated" : "predicted");

// --------------------------------------------------------------------------- //
// Shared scaffolding — mirrors the single-species main figures
// --------------------------------------------------------------------------- //

// The same header block as mainFigures.js: bold title, italic species, a factual
// subtitle and the scientific question the figure answers.
function headerBlock(fig, { title, species, subtitle, question }) {
  const P = fig.preset;
  let y = P.margin.top + P.font.title;
  if (title) {
    fig.text(P.margin.left, y, title, { size: "title", weight: "bold" });
    if (species) {
      fig.text(P.margin.left + textWidth(title, P.font.title) + 4, y,
        `· ${speciesLabel(species)}`, { size: "title", italic: true });
    }
    y += P.font.subtitle + 4;
  }
  if (subtitle) {
    fig.text(P.margin.left, y, subtitle, { size: "subtitle", fill: ALN_COLOURS.muted });
    y += P.font.subtitle + 3;
  }
  if (question) {
    fig.text(P.margin.left, y, question,
      { size: "small", fill: ALN_COLOURS.muted, italic: true });
    y += P.font.small + 4;
  }
  return y;
}

/** Left-hand track label in the reserved gutter, so it can never overlap a mark. */
function trackLabel(fig, x, yCentre, text, { bold = false, fill = ALN_COLOURS.ink } = {}) {
  fig.text(x, yCentre + fig.preset.font.small * 0.35, text,
    { size: "small", anchor: "end", fill, weight: bold ? "bold" : "normal" });
}

/**
 * Bin alignment columns down to the plot resolution.
 *
 * At 823 columns in a page-wide figure one column is well under a hairline wide,
 * so a per-column mark would alias away. Binning makes every mark at least
 * `minPt` wide; the figure footnote states how many columns a mark stands for,
 * because that is the difference between a drawing convention and a lie.
 */
function columnBins(nCols, x, minPt = 1.3) {
  const perColumn = (x(nCols) - x(0)) / Math.max(1, nCols);
  const size = Math.max(1, Math.ceil(minPt / Math.max(1e-6, perColumn)));
  const bins = [];
  for (let c0 = 0; c0 < nCols; c0 += size) {
    const c1 = Math.min(nCols - 1, c0 + size - 1);
    bins.push({ c0, c1, x: x(c0), w: Math.max(0.6, x(c1 + 1) - x(c0)) });
  }
  return { bins, size };
}

/** Mean of `values` over [c0, c1]. */
function binMean(values, c0, c1) {
  let sum = 0;
  for (let c = c0; c <= c1; c++) sum += values[c];
  return sum / (c1 - c0 + 1);
}

/** True when any flag in [c0, c1] is set. */
function binAny(flags, c0, c1) {
  for (let c = c0; c <= c1; c++) if (flags[c]) return true;
  return false;
}

/**
 * A density track: one bar per bin, drawn upwards on a fixed 0–1 scale.
 *
 * The scale is never stretched to the data, because a gap density of 0.25 and one
 * of 1.0 must not look the same. A dashed rule at 0.5 and the observed maximum in
 * the gutter make a short bar readable anyway.
 */
function densityTrack(fig, { bins, values, y, h, colour, x0, x1 }) {
  const P = fig.preset;
  fig.rect(x0, y, x1 - x0, h,
    { fill: ALN_COLOURS.trackBase, stroke: ALN_COLOURS.grid, lw: P.lw.thin });
  fig.line(x0, y + h / 2, x1, y + h / 2,
    { stroke: ALN_COLOURS.axis, lw: P.lw.thin, dash: "1.5 2.5" });
  for (const b of bins) {
    const v = Math.max(0, Math.min(1, binMean(values, b.c0, b.c1)));
    if (v <= 0) continue;
    const bh = Math.max(0.6, v * h);
    fig.rect(b.x, y + h - bh, b.w, bh, { fill: colour, stroke: "none" });
  }
  return y + h;
}

/** An opaque amber ribbon marking a candidate interval outside the plot content. */
function candidateRibbon(fig, { intervals, y, h, x, x0, x1 }) {
  const P = fig.preset;
  for (const c of intervals) {
    const cx = Math.max(x0, x(c.col_start));
    const cw = Math.max(1.2, Math.min(x1, x(c.col_end + 1)) - cx);
    fig.rect(cx, y, cw, h,
      { fill: ALN_COLOURS.candidate, stroke: ALN_COLOURS.candidateEdge, lw: P.lw.thin });
  }
  return y + h;
}

/** Alignment-column axis; the first and last column are always labelled. */
function columnAxis(fig, { x0, x1, y, nCols, x, label }) {
  const P = fig.preset;
  const { major, minor } = axisTicks(1, nCols, Math.max(5, Math.round((x1 - x0) / 80)));
  fig.line(x0, y, x1, y, { stroke: ALN_COLOURS.axis, lw: P.lw.rule });
  for (const t of minor) {
    if (t < 1 || t > nCols) continue;
    fig.line(x(t - 1), y, x(t - 1), y + 2,
      { stroke: ALN_COLOURS.axis, lw: P.lw.thin, opacity: 0.7 });
  }
  const drawn = [];
  const place = (value) => {
    const text = String(value);
    const px = x(value - 1);
    const w = textWidth(text, P.font.tick);
    if (drawn.some(([a, b]) => px - w / 2 < b + 3 && px + w / 2 > a - 3)) return;
    drawn.push([px - w / 2, px + w / 2]);
    fig.line(px, y, px, y + 3.5, { stroke: ALN_COLOURS.axis, lw: P.lw.rule });
    fig.text(px, y + 4 + P.font.tick, text,
      { size: "tick", anchor: "middle", fill: ALN_COLOURS.muted });
  };
  // Termini first, so a crowded interior tick yields to them.
  place(1);
  place(nCols);
  for (const t of major) if (t > 1 && t < nCols) place(t);
  let bottom = y + 4 + P.font.tick + 2;
  if (label) {
    fig.text((x0 + x1) / 2, bottom + P.font.label, label,
      { size: "label", anchor: "middle" });
    bottom += P.font.label + 3;
  }
  return bottom;
}

/**
 * Second axis carrying the primary protein's own residue numbering, so a column
 * position can be read back as a coordinate on the reference protein.
 */
function primaryAaAxis(fig, { x0, x1, y, colByAa, proteinLength, x, label }) {
  const P = fig.preset;
  const { major } = axisTicks(1, proteinLength, Math.max(4, Math.round((x1 - x0) / 110)));
  const values = [1, ...major.filter((t) => t > 1 && t < proteinLength), proteinLength];
  fig.line(x0, y, x1, y, { stroke: ALN_COLOURS.reference, lw: P.lw.thin, opacity: 0.8 });
  const drawn = [];
  for (const aa of values) {
    const col = colByAa(aa);
    if (col == null) continue;
    const px = x(col);
    const text = String(aa);
    const w = textWidth(text, P.font.tick);
    if (drawn.some(([a, b]) => px - w / 2 < b + 3 && px + w / 2 > a - 3)) continue;
    drawn.push([px - w / 2, px + w / 2]);
    fig.line(px, y, px, y + 3, { stroke: ALN_COLOURS.reference, lw: P.lw.thin });
    fig.text(px, y + 3.5 + P.font.tick, text,
      { size: "tick", anchor: "middle", fill: ALN_COLOURS.reference });
  }
  let bottom = y + 3.5 + P.font.tick + 2;
  if (label) {
    fig.text((x0 + x1) / 2, bottom + P.font.small, label,
      { size: "small", anchor: "middle", fill: ALN_COLOURS.reference });
    bottom += P.font.small + 2;
  }
  return bottom;
}

const EXPLORATORY_NOTE = "Exploratory analysis; no isoform difference shown here is a "
  + "validated splicing event.";

const CONSERVATION_NOTE = "Conservation is the fraction of the shown models carrying the "
  + "column's most frequent residue, so a gapped column scores low.";

/** Map candidate intervals given in primary residues onto alignment columns. */
function mapCandidates(candidates, primarySeq) {
  return (candidates || [])
    .map((c) => ({
      ...c,
      label: c.label || c.rank_label || c.candidate_id || "candidate",
      col_start: c.col_start ?? aaToColumn(primarySeq, c.aa_start),
      col_end: c.col_end ?? aaToColumn(primarySeq, c.aa_end),
    }))
    .filter((c) => c.col_start != null && c.col_end != null)
    .sort((a, b) => a.col_start - b.col_start);
}

// --------------------------------------------------------------------------- //
// Figure 1 — Full isoform alignment overview
// --------------------------------------------------------------------------- //

/**
 * What is the difference, gap and conservation structure of the complete
 * within-species isoform alignment?
 *
 * One row per protein model over the whole column axis. Coverage shows which
 * columns a model occupies, marks inside the row show where it differs from the
 * primary, and three aggregate tracks characterise the columns themselves.
 */
export function alignmentOverviewFigureSpec({
  rows: rawRows, nCols, gene, species, primaryId, transcriptId,
  candidates = [], tool = "MAFFT", presetName = "full",
}) {
  const P = preset(presetName);
  const rows = orderRows(rawRows);
  const primary = rows.find((r) => r.is_primary) || rows[0];
  const primarySeq = primary?.seq || "";
  const profile = alignmentProfile(rows, nCols, primarySeq);
  const identity = new Map(rows.map((r) => [r.protein_id,
    r.is_primary ? "100%" : identityLabel(r.seq || "", primarySeq)]));

  // The gutter is sized from the text that actually goes into it, so a label can
  // never reach into the plot area.
  const idLines = rows.map((r) => `${r.protein_id}${r.is_primary ? " ★" : ""}`);
  const metaLines = rows.map((r) => [r.transcript_id,
    r.protein_length != null ? `${r.protein_length} aa` : null,
    curationMark(r)].filter(Boolean).join(" · "));
  const identLines = rows.map((r) => identity.get(r.protein_id) || "n/a");
  const gutter = Math.min(216, Math.max(126, Math.max(
    ...idLines.map((s, i) => textWidth(s, P.font.small) + 10
      + textWidth(identLines[i], P.font.small)),
    ...metaLines.map((s) => textWidth(s, P.font.small)),
  ) + 22));

  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right - 6;
  const x = (c) => x0 + (Math.max(0, Math.min(nCols, c)) / Math.max(1, nCols)) * (x1 - x0);
  const { bins, size: binSize } = columnBins(nCols, x);

  const ROW_H = 12;
  const ROW_PITCH = 21;
  const fig = createFigure({ preset: presetName, height: 560 });

  const nCurated = rows.filter((r) => curationMark(r) === "curated").length;
  let y = headerBlock(fig, {
    title: `${gene || "Gene"} · isoform alignment`,
    species,
    subtitle: `${rows.length} protein models (${nCurated} curated) · ${nCols} alignment `
      + `columns · ${tool} · primary ${primaryId || primary?.protein_id || "n/a"}`
      + (transcriptId ? ` / ${transcriptId}` : ""),
    question: "Where along the alignment do the protein isoforms of this gene differ "
      + "from the primary protein?",
  });
  y += 4;

  // --- candidate labels, above the body -------------------------------------
  const mapped = mapCandidates(candidates, primarySeq);
  if (mapped.length) {
    const placed = placeBlockLabels(mapped.map((c) => ({
      x0: x(c.col_start), x1: x(c.col_end + 1),
      label: `${c.label} · aa ${c.aa_start}–${c.aa_end}`,
    })), { size: P.font.small, rows: 2 });
    const rowsUsed = Math.max(1, ...placed.map((p) => (p.row ?? 0) + 1));
    trackLabel(fig, x0 - 8, y + P.font.small / 2, "Candidates",
      { fill: ALN_COLOURS.candidateEdge });
    for (const p of placed) {
      if (p.mode === "none") continue;
      const ly = y + (p.row ?? 0) * (P.font.small + 2) + P.font.small;
      const lx = Math.min(Math.max(p.labelX, x0 + p.width / 2), x1 - p.width / 2);
      fig.text(lx, ly, p.label,
        { size: "small", anchor: "middle", fill: ALN_COLOURS.candidateEdge });
    }
    y += rowsUsed * (P.font.small + 2) + 3;
  }

  // --- candidate interval ribbon, above the body ----------------------------
  const RIBBON_H = 4.5;
  if (mapped.length) {
    trackLabel(fig, x0 - 8, y + RIBBON_H / 2, "interval",
      { fill: ALN_COLOURS.candidateEdge });
    candidateRibbon(fig, { intervals: mapped, y, h: RIBBON_H, x, x0, x1 });
    y += RIBBON_H + 3;
  }

  // --- reserved-gutter column headings --------------------------------------
  const bodyTop = y + P.font.small + 2;
  fig.text(P.margin.left + 8, y + P.font.small, "protein model · transcript · length",
    { size: "small", fill: ALN_COLOURS.muted });
  fig.text(x0 - 6, y + P.font.small, "identity",
    { size: "small", anchor: "end", fill: ALN_COLOURS.muted });
  const bodyBottom = bodyTop + (rows.length - 1) * ROW_PITCH + ROW_H + 2;
  const colByAa = (aa) => aaToColumn(primarySeq, aa);
  const proteinLength = primary?.protein_length
    || (primarySeq ? primarySeq.replace(/-/g, "").length : nCols);

  // --- one row per protein model --------------------------------------------
  rows.forEach((r, i) => {
    const ry = bodyTop + i * ROW_PITCH;
    const seq = r.seq || "";
    const isPri = Boolean(r.is_primary);
    const centre = ry + ROW_H / 2 + P.font.small * 0.35 - 1;

    // A filled marker is a curated model, an open one a predicted model.
    fig.circle(P.margin.left + 2.6, ry + ROW_H / 2 - 1, 2.2, {
      fill: curationMark(r) === "curated" ? ALN_COLOURS.primary : ALN_COLOURS.paper,
      stroke: ALN_COLOURS.primary, lw: P.lw.thin,
    });
    fig.text(P.margin.left + 8, centre, idLines[i],
      { size: "small", weight: isPri ? "bold" : "normal" });
    fig.text(x0 - 6, centre, identLines[i],
      { size: "small", anchor: "end",
        fill: isPri ? ALN_COLOURS.ink : ALN_COLOURS.conservation });
    fig.text(P.margin.left + 8, ry + ROW_H + P.font.small + 0.5, metaLines[i],
      { size: "small", fill: ALN_COLOURS.muted });

    // Gap background first, then contiguous runs of aligned residues over it.
    fig.rect(x0, ry, x1 - x0, ROW_H, { fill: ALN_COLOURS.gap, stroke: "none" });
    const gapFlags = new Uint8Array(nCols);
    for (let c = 0; c < nCols; c++) gapFlags[c] = isGap(seq[c]) ? 1 : 0;
    let runStart = -1;
    for (let c = 0; c <= nCols; c++) {
      const present = c < nCols && !gapFlags[c];
      if (present && runStart < 0) runStart = c;
      if (!present && runStart >= 0) {
        fig.rect(x(runStart), ry, Math.max(0.6, x(c) - x(runStart)), ROW_H,
          { fill: isPri ? ALN_COLOURS.primary : ALN_COLOURS.match, stroke: "none" });
        runStart = -1;
      }
    }
    // Gap-block edges: the ends of an insertion or a deletion.
    for (const [gs, ge] of runsOf(gapFlags, 1)) {
      fig.line(x(gs), ry, x(gs), ry + ROW_H,
        { stroke: ALN_COLOURS.gapEdge, lw: P.lw.thin, opacity: 0.9 });
      fig.line(x(ge + 1), ry, x(ge + 1), ry + ROW_H,
        { stroke: ALN_COLOURS.gapEdge, lw: P.lw.thin, opacity: 0.9 });
    }
    // Differences from the primary, binned to the plot resolution.
    if (!isPri) {
      const differs = new Uint8Array(nCols);
      for (let c = 0; c < nCols; c++) {
        const ch = seq[c] || GAP;
        const pch = primarySeq[c] || GAP;
        differs[c] = !isGap(ch) && !isGap(pch) && ch !== pch ? 1 : 0;
      }
      for (const b of bins) {
        if (!binAny(differs, b.c0, b.c1)) continue;
        fig.rect(b.x, ry, b.w, ROW_H, { fill: ALN_COLOURS.diff, stroke: "none" });
      }
    }
    fig.rect(x0, ry, x1 - x0, ROW_H,
      { fill: "none", stroke: ALN_COLOURS.matchEdge, lw: P.lw.thin });
  });

  // Candidate interval edges through the body — opaque dashes, so the interval is
  // locatable in every row without an overlay that the PDF cannot reproduce.
  for (const c of mapped) {
    for (const cx of [x(c.col_start), x(c.col_end + 1)]) {
      fig.line(cx, bodyTop - 2, cx, bodyBottom,
        { stroke: ALN_COLOURS.candidateEdge, lw: P.lw.thin, dash: "2.5 2" });
    }
  }
  y = bodyTop + rows.length * ROW_PITCH + 4;

  // --- primary-protein coordinate mapping -----------------------------------
  // Which alignment columns carry a residue of the primary protein: an empty
  // stretch here is an insertion in some other model relative to the primary, and
  // is the reason a column index cannot simply be read as a residue number.
  const mapH = 6;
  trackLabel(fig, x0 - 8, y + mapH / 2, "Primary aa mapping");
  fig.rect(x0, y, x1 - x0, mapH,
    { fill: ALN_COLOURS.trackBase, stroke: ALN_COLOURS.grid, lw: P.lw.thin });
  const primaryOccupied = new Uint8Array(nCols);
  for (let c = 0; c < nCols; c++) primaryOccupied[c] = isGap(primarySeq[c]) ? 0 : 1;
  for (const [s, e] of runsOf(primaryOccupied, 1)) {
    fig.rect(x(s), y, Math.max(0.6, x(e + 1) - x(s)), mapH,
      { fill: ALN_COLOURS.reference, stroke: "none" });
  }
  for (const aa of axisTicks(1, proteinLength, 8).major) {
    const col = colByAa(aa);
    if (col == null) continue;
    fig.line(x(col), y - 1.5, x(col), y + mapH + 1.5,
      { stroke: ALN_COLOURS.ink, lw: P.lw.thin });
  }
  y += mapH + 7;

  // --- aggregate column tracks ---------------------------------------------
  const majorVar = runsOf(profile.variable, 3);
  const majorGap = runsOf(profile.gapped, 3);

  // 1. variable columns — where the shown models disagree at all
  const varH = 8;
  trackLabel(fig, x0 - 8, y + varH / 2, "Variable columns");
  fig.rect(x0, y, x1 - x0, varH,
    { fill: ALN_COLOURS.trackBase, stroke: ALN_COLOURS.grid, lw: P.lw.thin });
  for (const b of bins) {
    if (!binAny(profile.variable, b.c0, b.c1)) continue;
    fig.rect(b.x, y, b.w, varH, { fill: ALN_COLOURS.block, stroke: "none", opacity: 0.85 });
  }
  y += varH;
  // The three longest blocks are labelled with their column range; the rest are
  // counted in the footnote rather than crowding the track.
  const longest = [...majorVar].sort((a, b) => (b[1] - b[0]) - (a[1] - a[0])).slice(0, 3)
    .sort((a, b) => a[0] - b[0]);
  if (longest.length) {
    const placed = placeBlockLabels(longest.map(([s, e]) => ({
      x0: x(s), x1: x(e + 1), label: `${s + 1}–${e + 1}`,
    })), { size: P.font.small, rows: 1 });
    let used = false;
    for (const p of placed) {
      if (p.mode === "none") continue;
      used = true;
      const lx = Math.min(Math.max(p.labelX, x0 + p.width / 2), x1 - p.width / 2);
      fig.text(lx, y + P.font.small + 1.5, p.label,
        { size: "small", anchor: "middle", fill: ALN_COLOURS.muted });
    }
    if (used) y += P.font.small + 3;
  }
  y += 5;

  // 2. substitution / difference density, 3. gap density, 4. conservation
  const TRACK_H = 16;
  for (const [label, values, colour] of [
    ["Difference density", profile.diffFraction, ALN_COLOURS.diff],
    ["Gap density", profile.gapFraction, ALN_COLOURS.gapDensity],
    ["Conservation", profile.conservation, ALN_COLOURS.conservation],
  ]) {
    let max = 0;
    for (let c = 0; c < nCols; c++) if (values[c] > max) max = values[c];
    // Scale 0–1 with a rule at 0.5; the observed maximum keeps a short bar readable.
    trackLabel(fig, x0 - 8, y + TRACK_H / 2 - P.font.small * 0.6, label);
    trackLabel(fig, x0 - 8, y + TRACK_H / 2 + P.font.small * 0.6,
      `0–1, max ${max.toFixed(2)}`, { fill: ALN_COLOURS.muted });
    densityTrack(fig, { bins, values, y, h: TRACK_H, colour, x0, x1 });
    y += TRACK_H + 5;
  }
  y += 3;

  // --- both coordinate axes -------------------------------------------------
  y = columnAxis(fig, { x0, x1, y, nCols, x, label: "Alignment column" });
  y += 2;
  y = primaryAaAxis(fig, { x0, x1, y, colByAa, proteinLength, x,
    label: `${primaryId || primary?.protein_id || "primary"} residue position (aa)` });
  y += 6;

  // --- legend and footnotes -------------------------------------------------
  const legend = [
    [ALN_COLOURS.primary, "primary protein (aligned residues)"],
    [ALN_COLOURS.match, "alternative model (aligned residues)"],
    [ALN_COLOURS.diff, "residue differing from primary"],
    [ALN_COLOURS.gap, "gap in this isoform"],
    [ALN_COLOURS.block, "variable column (models disagree)"],
    [ALN_COLOURS.gapDensity, "gap density across models"],
    [ALN_COLOURS.conservation, "conservation (consensus fraction)"],
  ];
  if (mapped.length) legend.push([ALN_COLOURS.candidate, "exploratory candidate interval"]);
  y = fig.legend(P.margin.left, y + P.font.legend, legend);

  let varCols = 0;
  let gapCols = 0;
  for (let c = 0; c < nCols; c++) {
    varCols += profile.variable[c];
    gapCols += profile.gapped[c];
  }
  const notes = [
    `Major variable blocks: ${majorVar.length} · major gap blocks: ${majorGap.length} · `
    + `${varCols} variable and ${gapCols} gap-containing columns of ${nCols}.`,
    CONSERVATION_NOTE,
    binSize > 1
      ? `One plotted mark spans ${binSize} alignment columns at this width; a difference `
        + "or gap mark is drawn where any column within it qualifies."
      : null,
    "Filled marker: curated model · open marker: predicted model · identity: percentage "
    + "of aligned residue pairs identical to the primary protein.",
    EXPLORATORY_NOTE,
  ].filter(Boolean);
  for (const note of notes) {
    fig.text(P.margin.left, y + P.font.small, note,
      { size: "small", fill: ALN_COLOURS.muted, italic: note === EXPLORATORY_NOTE });
    y += P.font.small + 2.5;
  }

  return finalise(fig, y);
}

// --------------------------------------------------------------------------- //
// Residue-level block, shared by the wrapped export and the candidate detail
// --------------------------------------------------------------------------- //

/**
 * One wrapped block of the alignment at residue resolution.
 *
 * Draws, in this fixed order: the candidate band label, the alignment-column
 * ruler, one labelled row per model with residue letters, the gap-block
 * boundaries, the variable-column markers, the conservation bars, an optional
 * exon-association track and the primary-protein coordinate row. Returns the y
 * below the block.
 */
function residueBlock(fig, {
  rows, primarySeq, c0, c1, x0, cellW, y, rowH, profile, aaByColumn,
  labelFor, gutterRight, showResidues = true, bands = [], exons = null, tickEvery = 10,
}) {
  const P = fig.preset;
  const span = c1 - c0 + 1;
  const x = (c) => x0 + (c - c0) * cellW;
  const cellFill = Math.max(0.6, cellW - 0.4);
  const blockH = rows.length * rowH;
  const inBlock = bands
    .map((b) => ({ ...b, b0: Math.max(c0, b.col_start), b1: Math.min(c1, b.col_end) }))
    .filter((b) => b.b1 >= b.b0);
  let top = y;

  // --- candidate band labels and ribbon -------------------------------------
  if (bands.length) {
    const placed = placeBlockLabels(
      inBlock.map((b) => ({ x0: x(b.b0), x1: x(b.b1 + 1), label: b.label })),
      { size: P.font.small, rows: 1 });
    for (const p of placed) {
      if (p.mode === "none") continue;
      const lx = Math.min(Math.max(p.labelX, x0 + p.width / 2),
        x0 + span * cellW - p.width / 2);
      fig.text(lx, top + P.font.small, p.label,
        { size: "small", anchor: "middle", fill: ALN_COLOURS.candidateEdge, weight: "bold" });
    }
    top += P.font.small + 1;
    candidateRibbon(fig, { intervals: inBlock.map((b) => ({ col_start: b.b0, col_end: b.b1 })),
      y: top, h: 2.6, x, x0, x1: x0 + span * cellW });
    top += 3.4;
  }

  // --- alignment-column ruler ----------------------------------------------
  const rulerY = top + P.font.tick;
  for (let c = c0; c <= c1; c++) {
    if ((c + 1) % tickEvery) continue;
    fig.text(x(c) + cellW / 2, rulerY, String(c + 1),
      { size: "tick", anchor: "middle", fill: ALN_COLOURS.muted });
    fig.line(x(c) + cellW / 2, rulerY + 1.5, x(c) + cellW / 2, rulerY + 3.5,
      { stroke: ALN_COLOURS.axis, lw: P.lw.thin });
  }
  fig.text(gutterRight, rulerY, `cols ${c0 + 1}–${c1 + 1}`,
    { size: "tick", anchor: "end", fill: ALN_COLOURS.muted });
  const gridTop = rulerY + 4;

  // Which columns of this block lie inside a candidate interval. The interval is
  // shown by warming the residue cells rather than by covering them, so the
  // residue letters stay fully legible.
  const inBand = new Uint8Array(span);
  for (const b of inBlock) {
    for (let c = b.b0; c <= b.b1; c++) inBand[c - c0] = 1;
  }

  // --- residue rows ---------------------------------------------------------
  const letterPt = Math.max(5.5, Math.min(P.font.label, cellW * 0.94));
  rows.forEach((r, i) => {
    const ry = gridTop + i * rowH;
    const seq = r.seq || "";
    const isPri = Boolean(r.is_primary);
    fig.text(gutterRight, ry + rowH / 2 + P.font.small * 0.34, labelFor(r),
      { size: "small", anchor: "end", weight: isPri ? "bold" : "normal" });
    for (let c = c0; c <= c1; c++) {
      const ch = seq[c] || GAP;
      const pch = primarySeq[c] || GAP;
      const gap = isGap(ch);
      const differs = !gap && !isGap(pch) && ch !== pch;
      const cx = x(c);
      const band = inBand[c - c0] === 1;
      fig.rect(cx, ry, cellFill, rowH - 0.8, {
        fill: gap ? ALN_COLOURS.gap
          : differs ? ALN_COLOURS.diffCell
            : isPri ? (band ? ALN_COLOURS.primaryBandCell : ALN_COLOURS.primaryCell)
              : (band ? ALN_COLOURS.matchBandCell : ALN_COLOURS.matchCell),
        stroke: "none",
      });
      if (gap) {
        // An explicit dash keeps a deletion readable in monochrome print.
        fig.line(cx + cellFill * 0.2, ry + (rowH - 0.8) / 2,
          cx + cellFill * 0.8, ry + (rowH - 0.8) / 2,
          { stroke: ALN_COLOURS.gapEdge, lw: P.lw.rule });
      } else if (showResidues) {
        fig.text(cx + cellFill / 2, ry + rowH / 2 + letterPt * 0.34, ch, {
          size: letterPt, anchor: "middle",
          fill: differs ? ALN_COLOURS.diffInk : isPri ? ALN_COLOURS.primary : ALN_COLOURS.ink,
          weight: differs || isPri ? "bold" : "normal",
        });
      }
    }
  });
  let bottom = gridTop + blockH;

  // --- candidate interval edges --------------------------------------------
  for (const b of inBlock) {
    if (b.col_start >= c0) {
      fig.line(x(b.b0), gridTop - 1.5, x(b.b0), bottom + 1.5,
        { stroke: ALN_COLOURS.candidateEdge, lw: P.lw.rule, dash: "2.5 2" });
    }
    if (b.col_end <= c1) {
      fig.line(x(b.b1 + 1), gridTop - 1.5, x(b.b1 + 1), bottom + 1.5,
        { stroke: ALN_COLOURS.candidateEdge, lw: P.lw.rule, dash: "2.5 2" });
    }
  }

  // --- gap-block boundaries across the whole block --------------------------
  for (let c = c0; c <= c1; c++) {
    if (!profile.gapped[c]) continue;
    if (c === 0 || !profile.gapped[c - 1]) {
      fig.line(x(c), gridTop - 1, x(c), bottom + 1,
        { stroke: ALN_COLOURS.gapEdge, lw: P.lw.rule, opacity: 0.9 });
    }
    if (c === profile.nCols - 1 || !profile.gapped[c + 1]) {
      fig.line(x(c) + cellW, gridTop - 1, x(c) + cellW, bottom + 1,
        { stroke: ALN_COLOURS.gapEdge, lw: P.lw.rule, opacity: 0.9 });
    }
  }

  // --- variable-column markers ---------------------------------------------
  const markH = 3.5;
  bottom += 2;
  trackLabel(fig, gutterRight, bottom + markH / 2, "variable", { fill: ALN_COLOURS.muted });
  for (let c = c0; c <= c1; c++) {
    if (!profile.variable[c]) continue;
    fig.rect(x(c), bottom, cellFill, markH, { fill: ALN_COLOURS.block, stroke: "none" });
  }
  bottom += markH + 2;

  // --- conservation bars ----------------------------------------------------
  const consH = 8;
  trackLabel(fig, gutterRight, bottom + consH / 2, "conservation", { fill: ALN_COLOURS.muted });
  fig.rect(x0, bottom, span * cellW, consH,
    { fill: ALN_COLOURS.trackBase, stroke: ALN_COLOURS.grid, lw: P.lw.thin });
  for (let c = c0; c <= c1; c++) {
    const v = Math.max(0, Math.min(1, profile.conservation[c]));
    if (v <= 0) continue;
    const h = Math.max(0.5, v * consH);
    fig.rect(x(c), bottom + consH - h, cellFill, h,
      { fill: ALN_COLOURS.conservation, stroke: "none" });
  }
  bottom += consH + 2;

  // --- exon association (candidate detail only) -----------------------------
  if (exons && exons.length) {
    const exonH = 8;
    trackLabel(fig, gutterRight, bottom + exonH / 2, "primary exons",
      { fill: ALN_COLOURS.muted });
    for (const e of exons) {
      const e0 = Math.max(c0, e.col_start);
      const e1 = Math.min(c1, e.col_end);
      if (e1 < e0) continue;
      const ex = x(e0);
      const ew = Math.max(1, x(e1 + 1) - ex);
      fig.rect(ex, bottom, ew, exonH,
        { fill: ALN_COLOURS.exon, stroke: ALN_COLOURS.exonEdge, lw: P.lw.thin });
      if (ew >= textWidth(e.label, P.font.small) + 3) {
        fig.text(ex + ew / 2, bottom + exonH / 2 + P.font.small * 0.34, e.label,
          { size: "small", anchor: "middle", fill: ALN_COLOURS.ink });
      }
    }
    bottom += exonH + 2;
  }

  // --- primary-protein coordinate row --------------------------------------
  const aaY = bottom + P.font.tick;
  trackLabel(fig, gutterRight, bottom + P.font.tick / 2 + 1, "primary aa",
    { fill: ALN_COLOURS.reference });
  let lastX = -Infinity;
  for (let c = c0; c <= c1; c++) {
    const aa = aaByColumn.get(c);
    if (aa == null || (aa !== 1 && aa % tickEvery)) continue;
    const px = x(c) + cellW / 2;
    const w = textWidth(String(aa), P.font.tick);
    if (px - w / 2 < lastX + 2) continue;
    lastX = px + w / 2;
    fig.text(px, aaY, String(aa),
      { size: "tick", anchor: "middle", fill: ALN_COLOURS.reference });
  }
  return aaY + 2;
}

/** Height of one residue block, so pages can be filled before they are drawn. */
function residueBlockHeight(P, nRows, rowH, { bands = false, exons = false } = {}) {
  return (bands ? P.font.small + 4.4 : 0)   // candidate label and ribbon
    + P.font.tick + 4                       // column ruler
    + nRows * rowH                          // residue rows
    + 2 + 3.5 + 2                           // variable-column markers
    + 8 + 2                                 // conservation
    + (exons ? 10 : 0)                      // exon association
    + P.font.tick + 2;                      // primary-aa row
}

// --------------------------------------------------------------------------- //
// Figure 2 — Wrapped residue-level alignment
// --------------------------------------------------------------------------- //

const WRAP_HEADER_PT = 62;
const WRAP_FOOTER_PT = 58;
const WRAP_BLOCK_GAP = 8;

const wrapRowHeight = (P) => Math.max(8.2, P.font.small + 2.6);

/**
 * Block and page geometry of the wrapped export, without building the figures.
 *
 * Blocks are distributed evenly over the pages they need, so the last page is
 * not left with a single block while the others are full.
 */
export function wrappedAlignmentLayout({ nRows, nCols, colsPerBlock = 80,
  presetName = "full", hasCandidates = true }) {
  const P = preset(presetName);
  const perBlock = Math.max(60, Math.min(100, Math.round(colsPerBlock)));
  const blockH = residueBlockHeight(P, nRows, wrapRowHeight(P), { bands: hasCandidates })
    + WRAP_BLOCK_GAP;
  const usable = P.maxHeightPt - WRAP_HEADER_PT - WRAP_FOOTER_PT - P.margin.bottom;
  const capacity = Math.max(1, Math.floor(usable / blockH));
  const nBlocks = Math.max(1, Math.ceil(nCols / perBlock));
  const nPages = Math.max(1, Math.ceil(nBlocks / capacity));
  return { colsPerBlock: perBlock, nBlocks, nPages,
    blocksPerPage: Math.ceil(nBlocks / nPages), pageCapacity: capacity };
}

/**
 * The complete alignment at residue resolution, wrapped Jalview-style.
 *
 * An 823-column alignment cannot be printed legibly on one line, so it is cut
 * into blocks of `colsPerBlock` columns and the blocks are spread over as many
 * pages as the page height allows. Every block repeats the sequence labels, the
 * conservation row, the variable-column markers and both coordinate systems, so
 * a page can be read on its own.
 *
 * Returns one figure specification per page.
 */
export function wrappedAlignmentFigureSpecs({
  rows: rawRows, nCols, gene, species, primaryId, transcriptId,
  candidates = [], tool = "MAFFT", colsPerBlock = 80, presetName = "full",
}) {
  const P = preset(presetName);
  const rows = orderRows(rawRows);
  const primary = rows.find((r) => r.is_primary) || rows[0];
  const primarySeq = primary?.seq || "";
  const profile = alignmentProfile(rows, nCols, primarySeq);
  const aaByColumn = primaryResidueByColumn(primarySeq);
  const bands = mapCandidates(candidates, primarySeq);

  const labelFor = (r) => `${r.protein_id}${r.is_primary ? " ★" : ""}`;
  const gutter = Math.min(154, Math.max(80, Math.max(
    ...rows.map((r) => textWidth(labelFor(r), P.font.small)),
    textWidth("conservation", P.font.small)) + 12));
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right - 4;
  const rowH = wrapRowHeight(P);

  const layout = wrappedAlignmentLayout({ nRows: rows.length, nCols, colsPerBlock,
    presetName, hasCandidates: bands.length > 0 });
  const cellW = (x1 - x0) / layout.colsPerBlock;

  const blocks = [];
  for (let c0 = 0; c0 < nCols; c0 += layout.colsPerBlock) {
    blocks.push([c0, Math.min(nCols - 1, c0 + layout.colsPerBlock - 1)]);
  }
  const pages = [];
  for (let i = 0; i < blocks.length; i += layout.blocksPerPage) {
    pages.push(blocks.slice(i, i + layout.blocksPerPage));
  }

  return pages.map((pageBlocks, pageIndex) => {
    const fig = createFigure({ preset: presetName, height: P.maxHeightPt });
    let y = headerBlock(fig, {
      title: `${gene || "Gene"} · isoform alignment (residue level)`,
      species,
      subtitle: `${rows.length} protein models · columns ${pageBlocks[0][0] + 1}–`
        + `${pageBlocks[pageBlocks.length - 1][1] + 1} of ${nCols} · page ${pageIndex + 1} of `
        + `${pages.length} · ${layout.colsPerBlock} alignment columns per block · ${tool} · `
        + `primary ${primaryId || primary?.protein_id || "n/a"}`
        + (transcriptId ? ` / ${transcriptId}` : ""),
      question: pageIndex === 0
        ? "Which residues actually differ between the protein isoforms of this gene?"
        : "Continued: which residues differ between the protein isoforms of this gene?",
    });
    y += 2;

    for (const [c0, c1] of pageBlocks) {
      y = residueBlock(fig, {
        rows, primarySeq, c0, c1, x0, cellW, y, rowH, profile, aaByColumn,
        labelFor, gutterRight: x0 - 5,
        bands: bands.filter((b) => b.col_end >= c0 && b.col_start <= c1),
      }) + WRAP_BLOCK_GAP;
    }

    y = fig.legend(P.margin.left, y + P.font.legend, [
      [ALN_COLOURS.primaryCell, "primary protein residue"],
      [ALN_COLOURS.matchCell, "residue matching the primary"],
      [ALN_COLOURS.diffCell, "residue differing from primary"],
      [ALN_COLOURS.gap, "gap in this isoform"],
      [ALN_COLOURS.block, "variable column"],
      [ALN_COLOURS.conservation, "conservation"],
      ...(bands.length ? [[ALN_COLOURS.candidate, "exploratory candidate interval"]] : []),
    ]);
    fig.text(P.margin.left, y + P.font.small,
      `Above each block: alignment column · below it: ${primaryId || primary?.protein_id
      || "primary"} residue position (aa). ${CONSERVATION_NOTE}`,
      { size: "small", fill: ALN_COLOURS.muted });
    y += P.font.small + 2.5;
    fig.text(P.margin.left, y + P.font.small, EXPLORATORY_NOTE,
      { size: "small", fill: ALN_COLOURS.muted, italic: true });
    y += P.font.small + 2;

    return finalise(fig, y);
  });
}

// --------------------------------------------------------------------------- //
// Figure 3 — Candidate-focused alignment detail
// --------------------------------------------------------------------------- //

/**
 * What exactly happens to each protein model inside one candidate interval?
 *
 * Every model is shown, whether or not it is affected, because "unaffected" is a
 * result rather than a reason to omit a row. Residue letters, insertions,
 * deletions, substitutions, variable columns, gap boundaries, the coding exons of
 * the primary protein the interval falls in, per-model identity inside the
 * interval and both coordinate systems are all present.
 *
 * This is a companion to the full overview, never a replacement for it.
 */
export function candidateAlignmentFigureSpec({
  rows: rawRows, nCols, gene, species, primaryId, candidate,
  affected = null, exons = [], flank = 12, maxColsPerBlock = 76, presetName = "full",
}) {
  const P = preset(presetName);
  const rows = orderRows(rawRows);
  const primary = rows.find((r) => r.is_primary) || rows[0];
  const primarySeq = primary?.seq || "";
  const profile = alignmentProfile(rows, nCols, primarySeq);
  const aaByColumn = primaryResidueByColumn(primarySeq);

  const bandStart = candidate?.col_start ?? aaToColumn(primarySeq, candidate?.aa_start) ?? 0;
  const bandEnd = candidate?.col_end ?? aaToColumn(primarySeq, candidate?.aa_end) ?? nCols - 1;
  const c0 = Math.max(0, bandStart - flank);
  const c1 = Math.min(nCols - 1, bandEnd + flank);
  const span = c1 - c0 + 1;
  const affectedSet = affected instanceof Set
    ? affected
    : new Set(affected || affectedProteins(rows, primarySeq, bandStart, bandEnd));

  // Identity inside the candidate interval only — the number the candidate is
  // actually about, which a whole-alignment identity averages away. A model that
  // is entirely gapped there has no residue pair to compare, which is a different
  // statement from "0% identical".
  const localId = new Map(rows.map((r) => [r.protein_id,
    r.is_primary ? "100%" : identityLabel(r.seq || "", primarySeq, bandStart, bandEnd)]));
  const roleOf = (r) => (r.is_primary ? "reference"
    : affectedSet.has(r.protein_id) ? "affected" : "unaffected");
  const labelFor = (r) => `${r.protein_id}${r.is_primary ? " ★" : ""} · ${roleOf(r)} · `
    + `${localId.get(r.protein_id)}`;

  const gutter = Math.min(206, Math.max(126, Math.max(
    ...rows.map((r) => textWidth(labelFor(r), P.font.small)),
    textWidth("conservation", P.font.small)) + 12));
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right - 4;

  // The interval is wrapped rather than squeezed, so the residue letters stay
  // legible however wide the candidate is.
  const nBlocks = Math.max(1, Math.ceil(span / maxColsPerBlock));
  const perBlock = Math.ceil(span / nBlocks);
  const cellW = (x1 - x0) / perBlock;
  const rowH = Math.max(9, P.font.label + 2.4);

  // Coding exons of the primary protein, projected onto alignment columns.
  const exonBands = (exons || [])
    .map((e) => {
      const s = Number(e.start ?? e.start_aa);
      const t = Number(e.end ?? e.end_aa);
      const cs = aaToColumn(primarySeq, s);
      const ce = aaToColumn(primarySeq, t);
      if (!Number.isFinite(s) || !Number.isFinite(t) || cs == null || ce == null) return null;
      return { label: e.label || "", col_start: cs, col_end: ce, start: s, end: t };
    })
    .filter((e) => e && e.col_end >= c0 && e.col_start <= c1);
  const band = { label: candidate?.label || "candidate",
    col_start: bandStart, col_end: bandEnd };
  const insideExons = exonBands
    .filter((e) => e.col_end >= bandStart && e.col_start <= bandEnd)
    .map((e) => e.label).filter(Boolean);

  const fig = createFigure({ preset: presetName, height: P.maxHeightPt });
  let y = headerBlock(fig, {
    title: `${gene || "Gene"} · ${candidate?.label || "candidate"} alignment detail`,
    species,
    subtitle: `primary aa ${candidate?.aa_start ?? "?"}–${candidate?.aa_end ?? "?"} on `
      + `${primaryId || primary?.protein_id || "n/a"} = alignment columns ${bandStart + 1}–`
      + `${bandEnd + 1} · ${affectedSet.size} of ${Math.max(0, rows.length - 1)} alternative `
      + `models affected · ±${flank} columns of flanking context`
      + (insideExons.length ? ` · primary exons ${insideExons.join(", ")}` : ""),
    question: "Which residues does this candidate interval change, and in which protein "
      + "models?",
  });
  y += 2;

  for (let i = 0; i < nBlocks; i++) {
    const b0 = c0 + i * perBlock;
    const b1 = Math.min(c1, b0 + perBlock - 1);
    if (b1 < b0) break;
    y = residueBlock(fig, {
      rows, primarySeq, c0: b0, c1: b1, x0, cellW, y, rowH, profile, aaByColumn,
      labelFor, gutterRight: x0 - 5, exons: exonBands,
      bands: (band.col_end >= b0 && band.col_start <= b1) ? [band] : [],
    }) + 9;
  }

  let nVar = 0;
  let nGapCols = 0;
  for (let c = bandStart; c <= bandEnd; c++) {
    nVar += profile.variable[c];
    nGapCols += profile.gapped[c];
  }

  y = fig.legend(P.margin.left, y + P.font.legend, [
    [ALN_COLOURS.primaryCell, "primary protein residue (reference)"],
    [ALN_COLOURS.matchCell, "residue matching the primary"],
    [ALN_COLOURS.diffCell, "residue differing from primary (substitution)"],
    [ALN_COLOURS.gap, "gap in this isoform (insertion / deletion)"],
    [ALN_COLOURS.block, "variable column"],
    [ALN_COLOURS.conservation, "conservation"],
    [ALN_COLOURS.exon, "coding exon of the primary protein"],
    [ALN_COLOURS.candidate, "candidate interval"],
  ]);
  for (const note of [
    `${nVar} variable and ${nGapCols} gap-containing columns inside the candidate interval `
    + "· row label: protein model, role, identity within the interval.",
    `Above each block: alignment column · below it: ${primaryId || primary?.protein_id
    || "primary"} residue position (aa).`,
    "Exploratory candidate region — not a biologically validated splicing event.",
  ]) {
    fig.text(P.margin.left, y + P.font.small, note,
      { size: "small", fill: ALN_COLOURS.muted, italic: note.startsWith("Exploratory") });
    y += P.font.small + 2.5;
  }

  return finalise(fig, y);
}

// --------------------------------------------------------------------------- //
// Source data that ships with the figures
// --------------------------------------------------------------------------- //

/** The exact alignment the figures were drawn from, as FASTA. */
export function alignmentFasta(rows, { width = 60 } = {}) {
  const out = [];
  for (const r of orderRows(rows)) {
    out.push(`>${[r.protein_id, r.transcript_id, r.is_primary ? "primary" : "alternative",
      r.protein_length != null ? `${r.protein_length}aa` : null,
      curationMark(r)].filter(Boolean).join(" ")}`);
    const seq = r.seq || "";
    for (let i = 0; i < seq.length; i += width) out.push(seq.slice(i, i + width));
  }
  return `${out.join("\n")}\n`;
}

export const ALIGNMENT_SUMMARY_COLUMNS = [
  "protein_id", "transcript_id", "is_primary", "curation_status", "protein_length",
  "alignment_columns", "identity_to_primary_pct", "gap_columns", "differing_columns",
  "candidate_id", "candidate_alignment_columns", "candidate_identity_pct", "candidate_role",
];

/** Per-model source table for the alignment figures. */
export function alignmentSummaryRows({ rows: rawRows, nCols, candidates = [] }) {
  const rows = orderRows(rawRows);
  const primary = rows.find((r) => r.is_primary) || rows[0];
  const primarySeq = primary?.seq || "";
  const first = mapCandidates(candidates, primarySeq)[0];
  const affectedFirst = first
    ? affectedProteins(rows, primarySeq, first.col_start, first.col_end) : new Set();

  return rows.map((r) => {
    const seq = r.seq || "";
    let gaps = 0;
    let diff = 0;
    for (let c = 0; c < nCols; c++) {
      const ch = seq[c] || GAP;
      if (isGap(ch)) gaps += 1;
      else if (!isGap(primarySeq[c]) && ch !== primarySeq[c]) diff += 1;
    }
    return {
      protein_id: r.protein_id,
      transcript_id: r.transcript_id ?? "",
      is_primary: r.is_primary ? "true" : "false",
      curation_status: curationMark(r),
      protein_length: r.protein_length ?? seq.replace(/-/g, "").length,
      alignment_columns: nCols,
      identity_to_primary_pct: r.is_primary ? 100 : (identityPct(seq, primarySeq) ?? ""),
      gap_columns: gaps,
      differing_columns: diff,
      candidate_id: first?.label ?? "",
      candidate_alignment_columns: first ? `${first.col_start + 1}-${first.col_end + 1}` : "",
      candidate_identity_pct: first
        ? (r.is_primary ? 100
          : (identityPct(seq, primarySeq, first.col_start, first.col_end) ?? ""))
        : "",
      candidate_role: first
        ? (r.is_primary ? "reference"
          : (affectedFirst.has(r.protein_id) ? "affected" : "unaffected"))
        : "",
    };
  });
}

/** Trim the canvas to what was drawn, as the main figures do. */
function finalise(fig, contentBottom) {
  fig.resize(contentBottom + fig.preset.margin.bottom);
  return fig;
}
