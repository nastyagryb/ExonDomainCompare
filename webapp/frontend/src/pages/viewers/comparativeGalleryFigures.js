// Publication figures for the Comparative Figure Gallery (Parts 3–5).
//
// Consumes the shared comparative_dataset_index (and, for Boundary figures, the
// canonical boundary_dashboard multi-species contract). No independent
// comparable-boundary or comparable-domain algorithm lives here.
//
// Every builder in this module draws a scientific visualisation. Numbers that
// would otherwise become a prose page (pairwise identity, boundary-alignment
// coverage, isoform counts) are shown as marks inside a figure, and their exact
// values ship in the card's source table — never as a card whose page is text.

import {
  createFigure, preset, PALETTE, textWidth, axisTicks, placeBlockLabels,
} from "./figureSpec.js";
import {
  headerBlock, finalise, speciesLabel, domainInstances,
  boundaryClassColour, boundaryClassLabel,
} from "./mainFigures.js";
import { domainInstanceFill, FEATURE_STYLES } from "./semanticStyles.js";
import { normaliseSyntenyIndex, unresolvedSpeciesRow } from "./syntenyModel.js";
import {
  syntenyNeighbourhoodFigureSpec, neighbourConservationMatrixFigureSpec,
} from "./syntenyFigures.js";
import { speciesTag, isSupported } from "./comparativeFigures.js";
import { canonClass } from "./boundaryClasses.js";
import { orderSpeciesIds, orderSpeciesRows } from "./speciesOrder.js";

// Annotation-matrix cell states. Pending and unavailable are deliberately pale
// and unsaturated so they cannot be read as a detection result.
const STATE_FILL = {
  detected: PALETTE.domain,
  "not detected": PALETTE.paper,
  pending: "#f6efe0",
  unavailable: "#f0e6e6",
  "uncertain mapping": PALETTE.exonAlt,
};
const STATE_ORDER = ["detected", "not detected", "pending", "unavailable",
  "uncertain mapping"];

// The canonical taxonomic order, shared with every other comparative view.
// Alphabetical order used to put Danio rerio between Callithrix and Equus,
// which reads as though a fish belonged among the primates and mammals.
function speciesOrder(rows) {
  return orderSpeciesRows(rows);
}

/**
 * Species gutter: the binomial in italics with the primary protein underneath.
 *
 * Returns the width the caller has to reserve, so no lane label can ever be
 * clipped by the plotting area.
 */
function speciesGutter(P, rows, extra = 0) {
  const names = rows.map((r) => speciesLabel(r.scientific_name || r.species_id));
  const ids = rows.map((r) => String(r.protein_id || ""));
  return Math.max(
    92 + extra,
    ...names.map((s) => textWidth(s, P.font.label) + 10 + extra),
    ...ids.map((s) => textWidth(s, P.font.small) + 10 + extra),
  );
}

function drawSpeciesLabel(fig, x, cy, row) {
  const P = fig.preset;
  fig.text(x, cy - 1, speciesLabel(row.scientific_name || row.species_id),
    { size: "label", anchor: "end", italic: true });
  if (row.protein_id) {
    fig.text(x, cy + P.font.small + 1, row.protein_id,
      { size: "small", anchor: "end", fill: PALETTE.muted });
  }
}

/** Parse an aligned FASTA into [{ header, id, seq }]. */
export function parseAlignedFasta(text) {
  const records = [];
  let cur = null;
  for (const line of String(text || "").split(/\r?\n/)) {
    if (line.startsWith(">")) {
      const header = line.slice(1).trim();
      cur = { header, id: header.split(/[\s|]+/)[0], seq: "" };
      records.push(cur);
    } else if (cur) {
      cur.seq += line.trim();
    }
  }
  return records.filter((r) => r.seq.length);
}

// --------------------------------------------------------------------------- //
// 4A. MSA-aligned exon architecture (principal comparative figure)
// --------------------------------------------------------------------------- //
/**
 * Coding exons of one primary protein per species on the shared MSA-column axis.
 *
 * Exon fills carry meaning rather than identity: an exon whose aligned interval
 * is matched in every species is drawn in the shared-exon colour, an exon
 * without such a counterpart in the alternative-exon colour. Consecutive shared
 * exons alternate only in outline weight, so the reader can still count them
 * without a rainbow implying a per-exon category.
 */
export function msaAlignedExonArchitectureFigureSpec({
  gene, exons = [], nColumns = 0, boundaryGroups = [], inventory = [],
  presetName = "full",
}) {
  const P = preset(presetName);
  const bySp = new Map();
  for (const e of exons) {
    if (e.msa_start_column == null || e.msa_end_column == null) continue;
    if (!bySp.has(e.species_id)) bySp.set(e.species_id, []);
    bySp.get(e.species_id).push(e);
  }
  for (const list of bySp.values()) {
    list.sort((a, b) => Number(a.msa_start_column) - Number(b.msa_start_column));
  }
  const invById = new Map(inventory.map((r) => [r.species_id, r]));
  const species = speciesOrder([...bySp.keys()].map((id) => {
    const first = bySp.get(id)[0];
    return {
      species_id: id,
      scientific_name: invById.get(id)?.scientific_name || first.scientific_name,
      protein_id: invById.get(id)?.protein_id || first.protein_id,
    };
  }));
  const cols = Math.max(1, nColumns
    || Math.max(0, ...exons.map((e) => Number(e.msa_end_column) || 0)));

  // An exon is "comparable" when every other species has an exon overlapping the
  // same aligned interval. That is an alignment observation, not an orthology
  // claim, which is why the legend says "aligned interval matched".
  const overlaps = (a, b) => Number(a.msa_start_column) <= Number(b.msa_end_column)
    && Number(b.msa_start_column) <= Number(a.msa_end_column);
  const isShared = (sid, exon) => [...bySp.keys()].every((other) => other === sid
    || (bySp.get(other) || []).some((o) => overlaps(exon, o)));

  const laneH = 34;
  const gutter = speciesGutter(P, species);
  const boundaryTrackH = boundaryGroups.length ? 16 : 0;
  const fig = createFigure({
    preset: presetName,
    height: species.length * laneH + boundaryTrackH + 240,
  });
  let y = headerBlock(fig, {
    title: gene,
    subtitle: "MSA-aligned exon architecture · one primary protein per species",
    question: "Where do the coding exons of the selected primary proteins fall on a "
      + "common alignment coordinate, and which exon intervals are matched across "
      + "species?",
  });
  y += 6;
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right - 34;
  const scale = (c) => x0 + ((Number(c) - 1) / cols) * (x1 - x0);

  // Supported comparable-boundary groups above the lanes: they are the positions
  // the Boundary Explorer treats as the same boundary in every species.
  const supported = boundaryGroups.filter((g) => isSupported(g.mapping_status)
    && g.msa_column != null);
  if (boundaryTrackH) {
    fig.text(x0 - 6, y + P.font.small, "comparable boundaries",
      { size: "small", anchor: "end", fill: PALETTE.muted });
    for (const g of supported) {
      const x = scale(g.msa_column);
      fig.line(x, y + 3, x, y + 11, { stroke: PALETTE.boundary, lw: P.lw.rule });
      fig.circle(x, y + 3, P.marker * 0.6, { fill: PALETTE.boundary, stroke: "none" });
    }
    y += boundaryTrackH;
  }

  const top = y;
  species.forEach((s, i) => {
    const cy = top + i * laneH + laneH * 0.42;
    drawSpeciesLabel(fig, x0 - 8, cy, s);
    // Alignment extent of this species: gap columns stay visible as bare rule.
    fig.rect(x0, cy - 1.5, x1 - x0, 3, { fill: PALETTE.grid, stroke: "none" });
    const mapped = bySp.get(s.species_id) || [];
    const blocks = mapped.map((e) => ({
      x0: scale(e.msa_start_column),
      x1: scale(e.msa_end_column),
      label: e.exon_label || "",
      shared: isShared(s.species_id, e),
    }));
    blocks.forEach((b) => {
      const style = b.shared ? FEATURE_STYLES.shared_exon : FEATURE_STYLES.alternative_exon;
      fig.rect(b.x0, cy - 7, Math.max(0.8, b.x1 - b.x0), 14, {
        fill: style.fill, stroke: style.stroke, lw: P.lw.thin,
      });
      // Exon-boundary markers: the ends of every coding exon on the shared axis.
      fig.line(b.x1, cy - 9.5, b.x1, cy + 9.5,
        { stroke: FEATURE_STYLES.exon_boundary_tick.stroke, lw: P.lw.thin });
    });
    // Labels only where they fit; a leader row keeps the rest readable.
    const placed = placeBlockLabels(blocks, { size: P.font.small, rows: 1 });
    placed.forEach((b) => {
      if (!b.label) return;
      if (b.mode === "inside") {
        fig.text(b.labelX, cy + P.font.small * 0.35, b.label,
          { size: "small", anchor: "middle", fill: PALETTE.ink });
      } else if (b.mode === "below") {
        fig.text(b.labelX, cy + 9 + P.font.small, b.label,
          { size: "small", anchor: "middle", fill: PALETTE.muted });
      }
    });
    fig.text(x1 + 4, cy + P.font.small * 0.35, `${mapped.length} exons`,
      { size: "small", fill: PALETTE.muted });
  });

  y = top + species.length * laneH + 6;
  const { major } = axisTicks(1, cols, 10);
  fig.line(x0, y, x1, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of major) {
    if (t < 1 || t > cols) continue;
    fig.line(scale(t), y, scale(t), y + 3.5, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(scale(t), y + 4 + P.font.tick, String(Math.round(t)),
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  }
  y += 4 + P.font.tick * 2 + 6;
  fig.text((x0 + x1) / 2, y, "MSA column of the cross-species primary-protein alignment",
    { size: "label", anchor: "middle" });
  y += P.font.label + 8;

  const legend = [
    [FEATURE_STYLES.shared_exon.fill, "Coding exon, aligned interval matched in every species"],
    [FEATURE_STYLES.alternative_exon.fill, "Coding exon without a matched interval in some species"],
    [PALETTE.grid, "Alignment extent (gap columns carry no exon)"],
  ];
  if (supported.length) legend.push([PALETTE.boundary, "Supported comparable boundary"]);
  y = fig.legend(P.margin.left, y, legend);
  fig.text(P.margin.left, y + P.font.small,
    "A shared column means the residues were aligned, not that they are functionally "
    + "equivalent. Vertical ticks mark coding-exon ends on the alignment axis.",
    { size: "small", fill: PALETTE.muted });
  return finalise(fig, y + P.font.small + 4);
}

// --------------------------------------------------------------------------- //
// 4B. Native-coordinate exon architecture (secondary / supplement)
// --------------------------------------------------------------------------- //
export function nativeExonArchitectureFigureSpec({
  gene, models = [], presetName = "full",
}) {
  const P = preset(presetName);
  const species = speciesOrder(models);
  const maxLen = Math.max(1, ...species.map((m) => Number(m.protein_length) || 1));
  const laneH = 34;
  const gutter = speciesGutter(P, species);
  const fig = createFigure({
    preset: presetName, height: species.length * laneH + 220,
  });
  let y = headerBlock(fig, {
    title: gene,
    subtitle: "Native-coordinate exon architecture · each primary protein on its own "
      + "amino-acid axis",
    question: "How long is each primary protein, and how are its coding exons and exon "
      + "boundaries distributed along its own sequence?",
  });
  y += 6;
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right - 40;
  const scale = (aa) => x0 + (Number(aa) / maxLen) * (x1 - x0);
  const top = y;
  species.forEach((m, i) => {
    const cy = top + i * laneH + laneH * 0.42;
    const len = Number(m.protein_length) || 1;
    drawSpeciesLabel(fig, x0 - 8, cy, m);
    fig.rect(x0, cy - 1.5, Math.max(1, scale(len) - x0), 3,
      { fill: PALETTE.grid, stroke: "none" });
    const blocks = (m.exons || []).map((e) => ({
      x0: scale(e.start), x1: scale(e.end), label: e.label || "",
    }));
    blocks.forEach((b) => {
      fig.rect(b.x0, cy - 7, Math.max(0.8, b.x1 - b.x0), 14, {
        fill: FEATURE_STYLES.coding_exon.fill,
        stroke: FEATURE_STYLES.coding_exon.stroke, lw: P.lw.thin,
      });
    });
    for (const b of m.exon_boundaries || []) {
      const pos = b.protein_position ?? b.boundary_position_aa;
      if (pos == null) continue;
      fig.line(scale(pos), cy - 10, scale(pos), cy + 10,
        { stroke: FEATURE_STYLES.exon_boundary_tick.stroke, lw: P.lw.thin });
    }
    placeBlockLabels(blocks, { size: P.font.small, rows: 1 }).forEach((b) => {
      if (!b.label || b.mode === "none") return;
      if (b.mode === "inside") {
        fig.text(b.labelX, cy + P.font.small * 0.35, b.label,
          { size: "small", anchor: "middle", fill: PALETTE.ink });
      } else {
        fig.text(b.labelX, cy + 10 + P.font.small, b.label,
          { size: "small", anchor: "middle", fill: PALETTE.muted });
      }
    });
    fig.text(scale(len) + 4, cy + P.font.small * 0.35, `${len} aa`,
      { size: "small", fill: PALETTE.muted });
  });
  y = top + species.length * laneH + 6;
  const { major } = axisTicks(0, maxLen, 10);
  fig.line(x0, y, x1, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of major) {
    if (t < 0 || t > maxLen) continue;
    fig.line(scale(t), y, scale(t), y + 3.5, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(scale(t), y + 4 + P.font.tick, String(t),
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  }
  y += 4 + P.font.tick * 2 + 6;
  fig.text((x0 + x1) / 2, y, "Amino-acid position on each species' own primary protein",
    { size: "label", anchor: "middle" });
  y += P.font.label + 8;
  y = fig.legend(P.margin.left, y, [
    [FEATURE_STYLES.coding_exon.fill, "Coding exon"],
    [FEATURE_STYLES.exon_boundary_tick.stroke, "Exon boundary"],
  ]);
  fig.text(P.margin.left, y + P.font.small,
    "Native coordinates are not homologous alignment coordinates: the same x position "
    + "in two lanes is not the same residue. Use the MSA-aligned exon architecture for "
    + "positional comparison.",
    { size: "small", fill: PALETTE.muted });
  return finalise(fig, y + P.font.small + 4);
}

// --------------------------------------------------------------------------- //
// 4C. Cross-species primary-protein MSA overview
// --------------------------------------------------------------------------- //
/**
 * Column-resolved overview of the cross-species primary-protein alignment.
 *
 * Panels, top to bottom:
 *   1. per-species residue / gap coverage on the alignment axis
 *   2. column state track: identical, mismatch, indel
 *   3. windowed conservation curve
 *   4. variable blocks (runs of consecutive non-identical columns) with labels
 *   5. exon boundaries and representative-domain intervals, where supplied
 *
 * With two species the state track is literally the pairwise comparison, so the
 * identity metric printed in the header is the same number the identity table
 * carries — displayed as a metric of this figure rather than as its own page.
 */
export function primaryMsaOverviewFigureSpec({
  gene, alignmentText = "", records: recordsIn = null, inventory = [],
  exons = [], domains = [], minVariableBlock = 5, presetName = "full",
}) {
  const P = preset(presetName);
  const records = recordsIn && recordsIn.length
    ? recordsIn : parseAlignedFasta(alignmentText);
  const byProtein = new Map(inventory.map((r) => [r.protein_id, r]));
  const rows = records.map((r) => {
    const inv = byProtein.get(r.id) || {};
    return {
      species_id: inv.species_id || r.id,
      scientific_name: inv.scientific_name || speciesLabel(inv.species_id || r.id),
      protein_id: r.id,
      seq: r.seq,
    };
  });
  const n = rows[0]?.seq?.length || 0;

  const fig = createFigure({
    preset: presetName, height: rows.length * 20 + 340,
  });
  if (!n || rows.length < 2) {
    const y0 = headerBlock(fig, {
      title: gene,
      subtitle: "Cross-species primary-protein MSA overview",
      question: "Where do the selected primary proteins agree, vary, or carry gaps?",
    });
    fig.text(P.margin.left, y0 + 10,
      "No cross-species primary-protein alignment with at least two sequences is "
      + "available for this run.", { size: "label", fill: PALETTE.muted });
    return finalise(fig, y0 + 26);
  }

  // Column states. A column is identical when every sequence carries the same
  // residue, an indel when at least one sequence is gapped, otherwise a mismatch.
  const state = new Array(n);
  let nIdentical = 0;
  let nMismatch = 0;
  let nIndel = 0;
  for (let c = 0; c < n; c += 1) {
    let gap = false;
    let same = true;
    const first = rows[0].seq[c];
    for (const r of rows) {
      const ch = r.seq[c];
      if (ch === "-" || ch === "." || ch === undefined) gap = true;
      else if (ch !== first) same = false;
    }
    state[c] = gap ? "indel" : same ? "identical" : "mismatch";
    if (state[c] === "identical") nIdentical += 1;
    else if (state[c] === "mismatch") nMismatch += 1;
    else nIndel += 1;
  }
  const ungapped = nIdentical + nMismatch;
  const percentIdentity = ungapped ? (nIdentical / ungapped) * 100 : 0;

  let y = headerBlock(fig, {
    title: gene,
    subtitle: `Cross-species primary-protein MSA overview · ${rows.length} species, `
      + `${n} alignment columns`,
    question: "Where do the selected primary proteins agree, where do residues differ, "
      + "and where does one species carry an insertion or deletion?",
  });
  fig.text(P.margin.left, y + P.font.label,
    `Pairwise identity ${percentIdentity.toFixed(1)}% over ${ungapped} ungapped columns `
    + `· ${nMismatch} mismatch columns · ${nIndel} indel columns`,
    { size: "label", weight: "bold" });
  y += P.font.label + 8;

  const gutter = speciesGutter(P, rows);
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right;
  const scale = (c) => x0 + (c / n) * (x1 - x0);
  // One mark per pixel column rather than per alignment column: at 822 columns
  // over ~600 pt a per-column rect would be sub-hairline and bloat the PDF.
  const bins = Math.min(n, Math.max(120, Math.round(x1 - x0)));
  const binOf = (c) => Math.floor((c / n) * bins);

  // --- panel 1: per-species residue / gap coverage --------------------------- //
  const laneH = 20;
  rows.forEach((r, i) => {
    const cy = y + i * laneH + laneH * 0.42;
    drawSpeciesLabel(fig, x0 - 8, cy, r);
    fig.rect(x0, cy - 5, x1 - x0, 10,
      { fill: PALETTE.grid, stroke: "none" });
    // Contiguous residue runs, so a gap is a visible hole rather than a colour.
    let runStart = null;
    for (let c = 0; c <= n; c += 1) {
      const ch = c < n ? r.seq[c] : "-";
      const isRes = ch !== "-" && ch !== "." && ch !== undefined;
      if (isRes && runStart == null) runStart = c;
      if (!isRes && runStart != null) {
        fig.rect(scale(runStart), cy - 5, Math.max(0.5, scale(c) - scale(runStart)), 10,
          { fill: FEATURE_STYLES.primary_sequence.fill, stroke: "none" });
        runStart = null;
      }
    }
  });
  y += rows.length * laneH + 6;

  // --- panel 2: column state track ------------------------------------------ //
  fig.text(x0 - 8, y + 7, "column state", { size: "small", anchor: "end", fill: PALETTE.muted });
  const binState = new Array(bins).fill("identical");
  for (let c = 0; c < n; c += 1) {
    const b = binOf(c);
    if (state[c] === "indel") binState[b] = "indel";
    else if (state[c] === "mismatch" && binState[b] !== "indel") binState[b] = "mismatch";
  }
  const stateFill = {
    identical: FEATURE_STYLES.conserved_region.fill,
    mismatch: FEATURE_STYLES.variable_region.fill,
    indel: PALETTE.unavailable_or_uncertain,
  };
  const binW = (x1 - x0) / bins;
  for (let b = 0; b < bins; b += 1) {
    fig.rect(x0 + b * binW, y, Math.max(0.5, binW), 11,
      { fill: stateFill[binState[b]], stroke: "none" });
  }
  y += 15;

  // --- panel 3: windowed conservation ---------------------------------------- //
  const win = Math.max(5, Math.round(n / 60));
  const consH = 26;
  fig.text(x0 - 8, y + consH / 2, "conservation", { size: "small", anchor: "end", fill: PALETTE.muted });
  fig.rect(x0, y, x1 - x0, consH, { fill: PALETTE.paper, stroke: PALETTE.grid, lw: P.lw.thin });
  let prev = null;
  for (let b = 0; b < bins; b += 1) {
    const centre = Math.round(((b + 0.5) / bins) * n);
    const lo = Math.max(0, centre - win);
    const hi = Math.min(n, centre + win);
    let ident = 0;
    for (let c = lo; c < hi; c += 1) if (state[c] === "identical") ident += 1;
    const frac = hi > lo ? ident / (hi - lo) : 0;
    const px = x0 + (b + 0.5) * binW;
    const py = y + consH - frac * consH;
    if (prev) fig.line(prev[0], prev[1], px, py, { stroke: PALETTE.identity, lw: P.lw.rule });
    prev = [px, py];
  }
  fig.text(x1 + 2, y + P.font.small, "1", { size: "small", fill: PALETTE.muted });
  fig.text(x1 + 2, y + consH, "0", { size: "small", fill: PALETTE.muted });
  y += consH + 8;

  // --- panel 4: variable blocks ---------------------------------------------- //
  const blocks = [];
  let bStart = null;
  for (let c = 0; c <= n; c += 1) {
    const variable = c < n && state[c] !== "identical";
    if (variable && bStart == null) bStart = c;
    if (!variable && bStart != null) {
      if (c - bStart >= minVariableBlock) blocks.push({ start: bStart, end: c - 1 });
      bStart = null;
    }
  }
  fig.text(x0 - 8, y + 8, "variable blocks", { size: "small", anchor: "end", fill: PALETTE.muted });
  const labelled = placeBlockLabels(blocks.map((b, i) => ({
    x0: scale(b.start), x1: scale(b.end + 1), label: `V${i + 1}`,
  })), { size: P.font.small, rows: 2 });
  blocks.forEach((b, i) => {
    fig.rect(scale(b.start), y, Math.max(0.8, scale(b.end + 1) - scale(b.start)), 9,
      { fill: FEATURE_STYLES.variable_region.fill,
        stroke: FEATURE_STYLES.variable_region.stroke, lw: P.lw.thin });
    const place = labelled[i];
    if (!place || place.mode === "none") return;
    if (place.mode === "inside") {
      fig.text(place.labelX, y + 7, place.label,
        { size: "small", anchor: "middle", fill: PALETTE.ink });
    } else {
      fig.text(place.labelX, y + 11 + (place.row + 1) * (P.font.small + 1), place.label,
        { size: "small", anchor: "middle", fill: PALETTE.muted });
    }
  });
  y += 9 + (blocks.length ? P.font.small * 2 + 6 : 4);

  // --- panel 5: exon boundaries and domains on the alignment axis ------------ //
  const mappedExons = exons.filter((e) => e.msa_end_column != null);
  if (mappedExons.length) {
    fig.text(x0 - 8, y + 6, "exon boundaries", { size: "small", anchor: "end", fill: PALETTE.muted });
    for (const e of mappedExons) {
      fig.line(scale(Number(e.msa_end_column)), y, scale(Number(e.msa_end_column)), y + 8,
        { stroke: FEATURE_STYLES.exon_boundary_tick.stroke, lw: P.lw.thin });
    }
    y += 12;
  }
  const mappedDomains = domains.filter((d) => d.msa_start_column != null);
  if (mappedDomains.length) {
    fig.text(x0 - 8, y + 7, "domains", { size: "small", anchor: "end", fill: PALETTE.muted });
    const seen = new Set();
    for (const d of mappedDomains) {
      const key = `${d.interpro_accession}:${d.instance_number}`;
      if (seen.has(key)) continue;
      seen.add(key);
      const a = scale(Number(d.msa_start_column));
      const b = scale(Number(d.msa_end_column));
      fig.rect(a, y, Math.max(0.8, b - a), 9, {
        fill: domainInstanceFill(d.order_along_protein), stroke: PALETTE.ink,
        lw: P.lw.thin, opacity: 0.85,
      });
      if (b - a > textWidth(d.label || "", P.font.small) + 4) {
        fig.text((a + b) / 2, y + 7, d.label,
          { size: "small", anchor: "middle", fill: PALETTE.paper });
      }
    }
    y += 13;
  }

  // --- axis and legend -------------------------------------------------------- //
  const { major } = axisTicks(1, n, 10);
  fig.line(x0, y, x1, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of major) {
    if (t < 1 || t > n) continue;
    fig.line(scale(t), y, scale(t), y + 3.5, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(scale(t), y + 4 + P.font.tick, String(t),
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  }
  y += 4 + P.font.tick * 2 + 6;
  fig.text((x0 + x1) / 2, y, "Alignment column", { size: "label", anchor: "middle" });
  y += P.font.label + 8;
  y = fig.legend(P.margin.left, y, [
    [FEATURE_STYLES.primary_sequence.fill, "Aligned residues of that species"],
    [FEATURE_STYLES.conserved_region.fill, "Identical column"],
    [FEATURE_STYLES.variable_region.fill, "Mismatch column / variable block"],
    [PALETTE.unavailable_or_uncertain, "Indel column (gap in at least one species)"],
  ]);
  fig.text(P.margin.left, y + P.font.small,
    `Variable blocks are runs of at least ${minVariableBlock} consecutive non-identical `
    + "columns. Identity is computed on ungapped columns only and describes sequence "
    + "similarity, not functional equivalence.",
    { size: "small", fill: PALETTE.muted });
  return finalise(fig, y + P.font.small + 4);
}

// --------------------------------------------------------------------------- //
// 4D. Comparative domain architecture
// --------------------------------------------------------------------------- //
/**
 * Representative domain instances of one primary protein per species.
 *
 * @param mode "native" (each species' own amino-acid axis) or "msa" (shared
 *             alignment columns). Repeated accessions stay separate instances.
 */
export function comparativeDomainArchitectureFigureSpec({
  gene, models = [], domains = [], mode = "native", nColumns = 0,
  exons = [], presetName = "full",
}) {
  const P = preset(presetName);
  const species = speciesOrder(models);
  const useMsa = mode === "msa";
  const maxX = useMsa
    ? Math.max(1, nColumns || Math.max(0, ...domains.map((d) => Number(d.msa_end_column) || 0)))
    : Math.max(1, ...species.map((m) => Number(m.protein_length) || 1));
  const laneH = 44;
  const gutter = speciesGutter(P, species);
  const fig = createFigure({
    preset: presetName, height: species.length * laneH + 260,
  });
  let y = headerBlock(fig, {
    title: gene,
    subtitle: useMsa
      ? "Comparative domain architecture · MSA-aligned coordinates"
      : "Comparative domain architecture · native amino-acid coordinates",
    question: "Which representative domain instances and membrane segments does each "
      + "primary protein carry, and do they occupy comparable positions?",
  });
  y += 6;
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right - 40;
  const scale = (v) => x0 + (Number(v) / maxX) * (x1 - x0);
  const top = y;

  const bySpDom = new Map();
  for (const d of domains) {
    if (!bySpDom.has(d.species_id)) bySpDom.set(d.species_id, []);
    bySpDom.get(d.species_id).push(d);
  }
  const bySpExon = new Map();
  for (const e of exons) {
    if (!bySpExon.has(e.species_id)) bySpExon.set(e.species_id, []);
    bySpExon.get(e.species_id).push(e);
  }

  // The legend carries the full label of every instance drawn, so a short label
  // inside a block never has to stand on its own.
  const legendEntries = new Map();

  species.forEach((m, i) => {
    const cy = top + i * laneH + laneH * 0.36;
    drawSpeciesLabel(fig, x0 - 8, cy, m);
    const len = useMsa ? maxX : (Number(m.protein_length) || 1);
    fig.rect(x0, cy - 2, Math.max(1, scale(len) - x0), 4,
      { fill: FEATURE_STYLES.protein_backbone.fill, stroke: "none", opacity: 0.55 });

    const pending = m.status && m.status !== "available";
    if (pending) {
      fig.text(x0 + 4, cy + P.font.small * 0.35,
        `domain annotation ${m.status}`, { size: "small", fill: PALETTE.muted, italic: true });
    }

    // Instances come from the canonical model in native mode and from the
    // comparative MSA projection in aligned mode; both keep repeated accessions
    // apart by instance number.
    const inst = useMsa
      ? (bySpDom.get(m.species_id) || [])
        .filter((d) => d.msa_start_column != null)
        .map((d) => ({
          x0: scale(d.msa_start_column), x1: scale(d.msa_end_column),
          label: d.label || d.interpro_accession,
          order: d.order_along_protein || 1,
          full: `${d.label || d.interpro_accession} · columns `
            + `${d.msa_start_column}–${d.msa_end_column}`,
        }))
      : domainInstances(m.representative_domains || []).map((d) => ({
        x0: scale(d.start), x1: scale(d.end),
        label: d.short_label, order: d.display_order, full: d.full_label,
      }));

    inst.forEach((d) => {
      fig.rect(d.x0, cy - 9, Math.max(0.8, d.x1 - d.x0), 18, {
        fill: domainInstanceFill(d.order), stroke: PALETTE.ink, lw: P.lw.thin,
      });
      legendEntries.set(d.label, domainInstanceFill(d.order));
    });
    placeBlockLabels(inst.map((d) => ({ x0: d.x0, x1: d.x1, label: d.label })),
      { size: P.font.small, rows: 2 }).forEach((b, bi) => {
      if (b.mode === "inside") {
        fig.text(b.labelX, cy + P.font.small * 0.35, b.label,
          { size: "small", anchor: "middle", fill: PALETTE.paper, weight: "bold" });
      } else if (b.mode === "below") {
        const ly = cy + 11 + (b.row + 1) * (P.font.small + 1.5);
        fig.line((inst[bi].x0 + inst[bi].x1) / 2, cy + 9,
          (inst[bi].x0 + inst[bi].x1) / 2, ly - P.font.small,
          { stroke: PALETTE.axis, lw: P.lw.thin, opacity: 0.8 });
        fig.text(b.labelX, ly, b.label,
          { size: "small", anchor: "middle", fill: PALETTE.ink });
      }
    });

    // Transmembrane helices sit below the domain track; in MSA mode they are only
    // drawn when the comparative layer supplied aligned coordinates for them.
    if (!useMsa) {
      for (const tm of m.tm_regions || []) {
        const a = scale(tm.start ?? tm.start_aa);
        const b = scale(tm.end ?? tm.end_aa);
        if (!Number.isFinite(a) || !Number.isFinite(b)) continue;
        fig.rect(a, cy + 10, Math.max(1.2, b - a), 4,
          { fill: FEATURE_STYLES.tm_helix.fill, stroke: "none" });
      }
    }
    // Optional exon-boundary markers, so domain edges can be read against them.
    const exonRows = bySpExon.get(m.species_id) || [];
    for (const e of exonRows) {
      const at = useMsa ? e.msa_end_column : e.native_end;
      if (at == null) continue;
      fig.line(scale(at), cy - 12, scale(at), cy - 9.5,
        { stroke: FEATURE_STYLES.exon_boundary_tick.stroke, lw: P.lw.thin });
    }
  });

  y = top + species.length * laneH + 4;
  const { major } = axisTicks(0, maxX, 10);
  fig.line(x0, y, x1, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of major) {
    if (t < 0 || t > maxX) continue;
    fig.line(scale(t), y, scale(t), y + 3.5, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(scale(t), y + 4 + P.font.tick, String(t),
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  }
  y += 4 + P.font.tick * 2 + 6;
  fig.text((x0 + x1) / 2, y,
    useMsa ? "MSA column of the cross-species primary-protein alignment"
      : "Amino-acid position on each species' own primary protein",
    { size: "label", anchor: "middle" });
  y += P.font.label + 8;

  const legend = [...legendEntries.entries()].map(([label, colour]) => [colour, label]);
  if (!useMsa) legend.push([FEATURE_STYLES.tm_helix.fill, "Transmembrane helix"]);
  legend.push([FEATURE_STYLES.exon_boundary_tick.stroke, "Exon boundary"]);
  y = fig.legend(P.margin.left, y, legend);
  fig.text(P.margin.left, y + P.font.small,
    useMsa
      ? "Repeated domains of one accession stay separate instances, identified by "
        + "instance number and aligned interval."
      : "Native coordinates are species-specific; the MSA-aligned panel shows the same "
        + "instances on a common alignment axis.",
    { size: "small", fill: PALETTE.muted });
  return finalise(fig, y + P.font.small + 4);
}

// --------------------------------------------------------------------------- //
// 4E. Combined exon + domain architecture, grouped per species
// --------------------------------------------------------------------------- //
/**
 * Domains and coding exons of one primary protein per species, in one figure.
 *
 * The exon architecture and the domain architecture answered halves of the same
 * question in separate figures, and reading them together meant holding one
 * panel in memory while looking at the other. Here each species owns a small
 * two-track group on a single axis — domains above, exons below — so the
 * relationship the whole analysis is about, where an exon boundary falls
 * relative to a domain edge, is visible in one glance and per species.
 *
 * @param mode "native" (each species' own amino-acid axis) or "msa" (shared
 *             alignment columns). Both are built by the same code from the same
 *             data, so no gene has a second implementation.
 */
export function comparativeExonDomainArchitectureFigureSpec({
  gene, models = [], domains = [], exons = [], mode = "native", nColumns = 0,
  presetName = "full",
}) {
  const P = preset(presetName);
  const useMsa = mode === "msa";
  const species = speciesOrder(models);
  if (!species.length) return null;

  const maxX = useMsa
    ? Math.max(1, nColumns || Math.max(0, ...exons.map((e) => Number(e.msa_end_column) || 0),
      ...domains.map((d) => Number(d.msa_end_column) || 0)))
    : Math.max(1, ...species.map((m) => Number(m.protein_length) || 1));

  const byId = (rows) => {
    const map = new Map();
    for (const r of rows) {
      if (!map.has(r.species_id)) map.set(r.species_id, []);
      map.get(r.species_id).push(r);
    }
    return map;
  };
  const domainsBySp = byId(domains);
  const exonsBySp = byId(exons);

  // Two tracks plus breathing room; the group height is what makes the pairing
  // readable, so it is generous rather than tight.
  const domainTrackH = 20;
  const exonTrackH = 14;
  const groupH = domainTrackH + exonTrackH + 30;
  // The gutter also carries each species' exon and domain counts, so no label
  // has to sit past the right edge of the plotting area.
  const counts = new Map(species.map((m) => {
    const nExons = useMsa
      ? (exonsBySp.get(m.species_id) || []).filter((e) => e.msa_start_column != null).length
      : (m.exons || []).length;
    const nDomains = useMsa
      ? (domainsBySp.get(m.species_id) || []).filter((d) => d.msa_start_column != null).length
      : domainInstances(m.representative_domains || []).length;
    return [m.species_id, `${nExons} exons · ${nDomains} domains`];
  }));
  const gutter = Math.max(
    speciesGutter(P, species),
    ...[...counts.values()].map((s) => textWidth(s, P.font.small) + 12),
  );
  const fig = createFigure({
    preset: presetName, height: species.length * groupH + 280,
  });
  let y = headerBlock(fig, {
    title: gene,
    subtitle: useMsa
      ? "Exon and domain architecture per species · MSA-aligned coordinates"
      : "Exon and domain architecture per species · native amino-acid coordinates",
    question: "For each species, where do the coding exon boundaries of the primary "
      + "protein fall relative to the edges of its annotated domains?",
  });
  y += 6;
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right - 42;
  const scale = (v) => x0 + (Number(v) / maxX) * (x1 - x0);
  const top = y;
  const legendEntries = new Map();
  const classesSeen = new Set();

  species.forEach((m, i) => {
    const gy = top + i * groupH;
    const domainY = gy + 12;
    const exonY = domainY + domainTrackH + 6;
    const len = useMsa ? maxX : (Number(m.protein_length) || 1);

    // The label sits between the two tracks it names, so the pairing is visible
    // even when a species has no domain annotation yet.
    drawSpeciesLabel(fig, x0 - 8, (domainY + exonY) / 2 - 3, m);
    fig.text(x0 - 8, exonY + P.font.small + 3, counts.get(m.species_id),
      { size: "small", anchor: "end", fill: PALETTE.muted });

    // Separator above every group but the first: the eye needs to see where one
    // species ends when two tracks belong together.
    if (i > 0) {
      fig.line(x0 - gutter, gy - 2, x1, gy - 2,
        { stroke: PALETTE.grid, lw: P.lw.thin, opacity: 0.7 });
    }

    fig.rect(x0, domainY - 1.5, Math.max(1, scale(len) - x0), 3,
      { fill: FEATURE_STYLES.protein_backbone.fill, stroke: "none", opacity: 0.55 });

    const pending = m.status && m.status !== "available";
    if (pending) {
      fig.text(x0 + 4, domainY + P.font.small * 0.35,
        `domain annotation ${m.status}`,
        { size: "small", fill: PALETTE.muted, italic: true });
    }

    const inst = useMsa
      ? (domainsBySp.get(m.species_id) || [])
        .filter((d) => d.msa_start_column != null)
        .map((d) => ({
          x0: scale(d.msa_start_column), x1: scale(d.msa_end_column),
          label: d.label || d.interpro_accession, order: d.order_along_protein || 1,
        }))
      : domainInstances(m.representative_domains || []).map((d) => ({
        x0: scale(d.start), x1: scale(d.end),
        label: d.short_label, order: d.display_order,
      }));

    inst.forEach((d) => {
      fig.rect(d.x0, domainY - 8, Math.max(0.8, d.x1 - d.x0), 16, {
        fill: domainInstanceFill(d.order), stroke: PALETTE.ink, lw: P.lw.thin,
      });
      legendEntries.set(d.label, domainInstanceFill(d.order));
    });
    placeBlockLabels(inst.map((d) => ({ x0: d.x0, x1: d.x1, label: d.label })),
      { size: P.font.small, rows: 1 }).forEach((b) => {
      if (b.mode === "inside") {
        fig.text(b.labelX, domainY + P.font.small * 0.35, b.label,
          { size: "small", anchor: "middle", fill: PALETTE.paper, weight: "bold" });
      }
    });
    if (!useMsa) {
      for (const tm of m.tm_regions || []) {
        const a = scale(tm.start ?? tm.start_aa);
        const b = scale(tm.end ?? tm.end_aa);
        if (!Number.isFinite(a) || !Number.isFinite(b)) continue;
        fig.rect(a, domainY - 12, Math.max(1.2, b - a), 3,
          { fill: FEATURE_STYLES.tm_helix.fill, stroke: "none" });
      }
    }

    // Exon track: the same species' coding exons on the same axis.
    const exonRows = useMsa
      ? (exonsBySp.get(m.species_id) || []).filter((e) => e.msa_start_column != null)
        .map((e) => ({ a: e.msa_start_column, b: e.msa_end_column, label: e.exon_label }))
      : (m.exons || []).map((e) => ({ a: e.start, b: e.end, label: e.label }));

    fig.rect(x0, exonY - 1, Math.max(1, scale(len) - x0), 2,
      { fill: PALETTE.grid, stroke: "none" });
    // A single exon colour with a visible stroke: adjacent exons stay countable
    // without a second fill that would read as a different kind of exon.
    exonRows.forEach((e) => {
      const a = scale(e.a);
      const b = scale(e.b);
      fig.rect(a, exonY - 5, Math.max(0.8, b - a), 10, {
        fill: FEATURE_STYLES.coding_exon.fill,
        stroke: FEATURE_STYLES.coding_exon.stroke, lw: P.lw.thin,
      });
    });

    // The connector is the point of the figure: each exon boundary is carried up
    // through the domain track so its position relative to a domain edge is read
    // directly instead of being estimated across two separate panels. The class
    // comes from the same boundary analysis the Boundary Explorer shows, so the
    // two views cannot give different accounts of one boundary.
    for (const b of m.exon_boundaries || []) {
      const at = useMsa ? b.msa_column : (b.protein_position ?? b.boundary_position_aa);
      if (at == null) continue;
      const cls = canonClass(b.boundary_class || b.class);
      const colour = cls ? boundaryClassColour(cls) : FEATURE_STYLES.exon_boundary_tick.stroke;
      if (cls) classesSeen.add(cls);
      const x = scale(at);
      fig.line(x, domainY - 9, x, exonY + 6,
        { stroke: colour, lw: P.lw.thin, opacity: cls ? 0.9 : 0.5 });
      fig.circle(x, exonY + 6, P.marker * 0.45, { fill: colour, stroke: "none" });
    }
  });

  y = top + species.length * groupH + 4;
  const { major } = axisTicks(0, maxX, 10);
  fig.line(x0, y, x1, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of major) {
    if (t < 0 || t > maxX) continue;
    fig.line(scale(t), y, scale(t), y + 3.5, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(scale(t), y + 4 + P.font.tick, String(Math.round(t)),
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  }
  y += 4 + P.font.tick * 2 + 6;
  fig.text((x0 + x1) / 2, y,
    useMsa ? "MSA column of the cross-species primary-protein alignment"
      : "Amino-acid position on each species' own primary protein",
    { size: "label", anchor: "middle" });
  y += P.font.label + 8;

  const legend = [...legendEntries.entries()].map(([label, colour]) => [colour, label]);
  legend.push([FEATURE_STYLES.coding_exon.fill, "Coding exon"]);
  if (!useMsa) legend.push([FEATURE_STYLES.tm_helix.fill, "Transmembrane helix"]);
  for (const cls of [...classesSeen].sort()) {
    legend.push([boundaryClassColour(cls), `Exon boundary · ${boundaryClassLabel(cls)}`]);
  }
  if (!classesSeen.size) {
    legend.push([FEATURE_STYLES.exon_boundary_tick.stroke, "Exon boundary"]);
  }
  y = fig.legend(P.margin.left, y, legend);
  fig.text(P.margin.left, y + P.font.small,
    useMsa
      ? "Both tracks of a species share one alignment axis, so a vertical connector "
        + "reads directly as the distance between an exon boundary and a domain edge. "
        + "An aligned column means residues were aligned, not that they are equivalent."
      : "Both tracks of a species share that species' own amino-acid axis. Native axes "
        + "are not comparable position by position between species; the MSA-aligned "
        + "panel places every species on one axis.",
    { size: "small", fill: PALETTE.muted });
  return finalise(fig, y + P.font.small + 4);
}

// --------------------------------------------------------------------------- //
// 4F. Domain annotation matrix (supplement)
// --------------------------------------------------------------------------- //
export function domainAnnotationMatrixFigureSpec({
  gene, matrix = [], groups = [], presetName = "double",
}) {
  const P = preset(presetName);
  const speciesIds = [...new Set(matrix.map((r) => r.species_id))];
  const species = speciesOrder(speciesIds.map((id) => {
    const hit = matrix.find((r) => r.species_id === id);
    return { species_id: id, scientific_name: hit?.scientific_name || id };
  }));
  const groupRows = groups.length
    ? groups
    : [...new Set(matrix.map((r) => r.comparable_domain_group_id))].map((id) => ({
      comparable_domain_group_id: id,
      label: matrix.find((r) => r.comparable_domain_group_id === id)?.label || id,
    }));
  const cols = Math.max(1, groupRows.length);

  const gutter = Math.max(80, ...species.map((s) =>
    textWidth(speciesLabel(s.scientific_name), P.font.label) + 10));
  const rowH = 20;
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right;
  const colW = (x1 - x0) / cols;
  // Two header lines: the group id, and as much of the label as fits the column.
  const headH = P.font.small * 2 + 8;
  const fig = createFigure({
    preset: presetName, height: species.length * rowH + headH + 230,
  });
  let y = headerBlock(fig, {
    title: gene,
    subtitle: "Domain annotation matrix · comparable representative-domain instances",
    question: "Which comparable domain instance is detected in which species, and where "
      + "is the annotation still pending?",
  });
  y += 8;

  const clip = (s, w) => {
    let out = String(s || "");
    while (out.length > 3 && textWidth(out, P.font.small) > w) out = out.slice(0, -1);
    return out.length < String(s || "").length ? `${out.slice(0, -1)}…` : out;
  };
  groupRows.forEach((g, ci) => {
    const cx = x0 + ci * colW + colW / 2;
    fig.text(cx, y + P.font.small, String(g.comparable_domain_group_id || ""),
      { size: "small", anchor: "middle", weight: "bold" });
    fig.text(cx, y + P.font.small * 2 + 2, clip(g.label, colW - 4),
      { size: "small", anchor: "middle", fill: PALETTE.muted });
  });
  y += headH;

  const cell = new Map(matrix.map((r) =>
    [`${r.species_id}|${r.comparable_domain_group_id}`, r]));
  const statesUsed = new Set();
  species.forEach((s, ri) => {
    const cy = y + ri * rowH + rowH / 2;
    fig.text(x0 - 8, cy + P.font.label * 0.35, speciesLabel(s.scientific_name),
      { size: "label", anchor: "end", italic: true });
    groupRows.forEach((g, ci) => {
      const r = cell.get(`${s.species_id}|${g.comparable_domain_group_id}`);
      const state = r?.state || "unavailable";
      statesUsed.add(state);
      fig.rect(x0 + ci * colW + 1.5, cy - 7, colW - 3, 14, {
        fill: STATE_FILL[state] || PALETTE.grid,
        stroke: PALETTE.axis, lw: P.lw.thin,
      });
      if (state === "detected" && r?.native_start != null) {
        const text = `${r.native_start}–${r.native_end}`;
        if (textWidth(text, P.font.small) < colW - 8) {
          fig.text(x0 + ci * colW + colW / 2, cy + P.font.small * 0.35, text,
            { size: "small", anchor: "middle", fill: PALETTE.paper });
        }
      }
    });
  });
  y += species.length * rowH + 10;

  y = fig.legend(P.margin.left, y,
    STATE_ORDER.filter((s) => statesUsed.has(s)).map((s) => [STATE_FILL[s], s]));
  // Full group labels, because a column header may have been clipped.
  const key = groupRows
    .map((g) => `${g.comparable_domain_group_id} = ${g.label}`).join("  ·  ");
  fig.text(P.margin.left, y + P.font.small, key, { size: "small", fill: PALETTE.muted });
  y += P.font.small + 3;
  fig.text(P.margin.left, y + P.font.small,
    "Cell text is the native amino-acid interval of the detected instance. "
    + "‘Not detected’ means no representative annotation in this run, not biological "
    + "absence.",
    { size: "small", fill: PALETTE.muted });
  return finalise(fig, y + P.font.small + 4);
}

// --------------------------------------------------------------------------- //
// 4H. Boundary-position consistency — compact multi-panel summary
// --------------------------------------------------------------------------- //
/**
 * One row per comparable boundary group, three aligned panels:
 *
 *   A  the raw signed distance of every species (the observations themselves)
 *   B  the cross-species absolute difference
 *   C  class agreement, mapping confidence and species coverage
 *
 * At n = 2 panel A is the pair of observations and panel B is their difference,
 * so nothing in the figure dresses a single subtraction up as a distribution.
 */
export function boundaryConsistencyPanelFigureSpec({
  gene, stats = [], groups = [], selectedGroupId = null,
  nearEdgeThreshold = 5, presetName = "full",
}) {
  const P = preset(presetName);
  const rows = groups
    .map((g) => {
      const s = stats.find((x) => x.comparable_boundary_group_id
        === g.comparable_boundary_group_id) || {};
      const obs = (g.per_species_native_positions || [])
        .filter((o) => o.signed_distance != null);
      const classes = [...new Set(obs.map((o) => canonClass(o.boundary_class)))];
      const dists = obs.map((o) => Number(o.signed_distance));
      return {
        id: g.comparable_boundary_group_id,
        label: String(g.comparable_boundary_group_id || "").replace(/^CBG/, "CBG "),
        obs,
        diff: dists.length > 1 ? Math.max(...dists) - Math.min(...dists) : null,
        classes,
        classAgreement: classes.length <= 1,
        confidence: s.mapping_confidence ?? g.confidence,
        coverage: s.species_coverage ?? g.species_coverage,
        nObserved: obs.length,
        nSpecies: s.n_species_available ?? g.n_species ?? obs.length,
        supported: isSupported(g.mapping_status),
      };
    })
    .filter((r) => r.obs.length);

  const speciesOrderIds = orderSpeciesIds(rows.flatMap((r) => r.obs.map((o) => o.species_id)));
  const speciesName = new Map();
  for (const r of rows) for (const o of r.obs) speciesName.set(o.species_id, o.scientific_name);

  const rowH = Math.max(15, P.font.label + 8);
  const gutter = Math.max(56, ...rows.map((r) => textWidth(r.label, P.font.label) + 10));
  const left = P.margin.left + gutter;
  const usable = P.widthPt - P.margin.right - left;
  // Panel widths: the raw observations get the most room, the status strip the least.
  const wA = usable * 0.46;
  const wB = usable * 0.28;
  const wC = usable * 0.20;
  const gap = usable * 0.03;
  const aX0 = left;
  const aX1 = aX0 + wA;
  const bX0 = aX1 + gap;
  const bX1 = bX0 + wB;
  const cX0 = bX1 + gap;
  const cX1 = cX0 + wC;

  const allDist = rows.flatMap((r) => r.obs.map((o) => Number(o.signed_distance)));
  const lim = Math.max(10, Math.ceil(Math.max(...allDist.map(Math.abs), 10) / 10) * 10);
  const sA = (v) => aX0 + ((v + lim) / (2 * lim)) * (aX1 - aX0);
  const maxDiff = Math.max(1, ...rows.map((r) => r.diff || 0));
  const limB = Math.ceil(maxDiff / 5) * 5 || 5;
  const sB = (v) => bX0 + (v / limB) * (bX1 - bX0);

  const fig = createFigure({ preset: presetName, height: rows.length * rowH + 290 });
  let y = headerBlock(fig, {
    title: gene,
    subtitle: "Boundary-position consistency · raw observations, cross-species "
      + "difference and mapping status",
    question: "For each comparable boundary, where does each species place it relative "
      + "to the nearest domain edge, how far apart are those placements, and how well "
      + "supported is the comparison?",
  });
  y += 6;

  // Panel captions.
  fig.text((aX0 + aX1) / 2, y + P.font.small, "A · raw signed distance per species (aa)",
    { size: "small", anchor: "middle", weight: "bold" });
  fig.text((bX0 + bX1) / 2, y + P.font.small, "B · cross-species difference (aa)",
    { size: "small", anchor: "middle", weight: "bold" });
  fig.text((cX0 + cX1) / 2, y + P.font.small, "C · class · confidence · coverage",
    { size: "small", anchor: "middle", weight: "bold" });
  y += P.font.small + 6;

  const top = y;
  const bodyH = rows.length * rowH;
  // Panel A background: the near-edge band and the domain-edge zero line.
  fig.rect(sA(-nearEdgeThreshold), top - 2,
    sA(nearEdgeThreshold) - sA(-nearEdgeThreshold), bodyH + 4,
    { fill: PALETTE.grid, stroke: "none", opacity: 0.85 });
  fig.line(sA(0), top - 2, sA(0), top + bodyH + 2, { stroke: PALETTE.ink, lw: P.lw.rule });

  rows.forEach((r, ri) => {
    const cy = top + ri * rowH + rowH / 2;
    const sel = selectedGroupId && r.id === selectedGroupId;
    if (sel) {
      fig.rect(P.margin.left, cy - rowH / 2, cX1 - P.margin.left, rowH,
        { fill: PALETTE.grid, stroke: "none", opacity: 0.5 });
    }
    fig.text(left - 6, cy + P.font.label * 0.35, r.label,
      { size: "label", anchor: "end", fill: sel ? PALETTE.ink : PALETTE.muted,
        weight: sel ? "bold" : "normal" });

    // --- A: the observations themselves ------------------------------------- //
    const laneGap = Math.min(6, rowH / Math.max(1, speciesOrderIds.length));
    const pts = r.obs.map((o) => {
      const li = Math.max(0, speciesOrderIds.indexOf(o.species_id));
      const dy = speciesOrderIds.length > 1
        ? (li - (speciesOrderIds.length - 1) / 2) * laneGap : 0;
      return { x: sA(Number(o.signed_distance)), y: cy + dy, o };
    });
    if (pts.length > 1) {
      const xs = pts.map((p) => p.x);
      fig.line(Math.min(...xs), cy, Math.max(...xs), cy, {
        stroke: r.supported ? PALETTE.boundary : PALETTE.muted,
        lw: r.supported ? P.lw.rule : P.lw.thin,
        dash: r.supported ? undefined : "1,2", opacity: 0.85,
      });
    }
    for (const p of pts) {
      const colour = boundaryClassColour(canonClass(p.o.boundary_class));
      if (p.o.nearest_edge === "start") {
        fig.circle(p.x, p.y, P.marker * 0.9,
          { fill: PALETTE.paper, stroke: colour, lw: P.lw.outline });
      } else {
        fig.circle(p.x, p.y, P.marker * 0.9, { fill: colour, stroke: "none" });
      }
    }

    // --- B: the difference --------------------------------------------------- //
    if (r.diff != null) {
      const w = Math.max(0.6, sB(r.diff) - bX0);
      fig.rect(bX0, cy - rowH * 0.24, w, rowH * 0.48, {
        fill: r.classAgreement ? PALETTE.boundary : PALETTE.outside_annotated_domains,
        stroke: "none", opacity: 0.85,
      });
      fig.text(bX0 + w + 3, cy + P.font.small * 0.35, `${r.diff}`,
        { size: "small", fill: PALETTE.muted });
    } else {
      fig.text(bX0, cy + P.font.small * 0.35, "single observation",
        { size: "small", fill: PALETTE.muted, italic: true });
    }

    // --- C: agreement, confidence, coverage ---------------------------------- //
    const glyphX = cX0 + 5;
    if (r.classAgreement) {
      fig.rect(glyphX - 3.5, cy - 3.5, 7, 7,
        { fill: PALETTE.identity, stroke: "none" });
    } else {
      fig.rect(glyphX - 3.5, cy - 3.5, 7, 7,
        { fill: PALETTE.paper, stroke: PALETTE.outside_annotated_domains, lw: P.lw.outline });
    }
    const confW = Math.max(0, Math.min(1, Number(r.confidence) || 0)) * (wC * 0.4);
    fig.rect(glyphX + 8, cy - 3, wC * 0.4, 6,
      { fill: PALETTE.paper, stroke: PALETTE.grid, lw: P.lw.thin });
    fig.rect(glyphX + 8, cy - 3, Math.max(0.5, confW), 6,
      { fill: PALETTE.domainAlt, stroke: "none" });
    fig.text(cX1, cy + P.font.small * 0.35, `${r.nObserved}/${r.nSpecies}`,
      { size: "small", anchor: "end", fill: PALETTE.muted });
  });

  y = top + bodyH + 8;
  // Axis for panel A.
  const { major: majA } = axisTicks(-lim, lim, 6);
  fig.line(aX0, y, aX1, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of majA) {
    if (t < -lim || t > lim) continue;
    fig.line(sA(t), y, sA(t), y + 3, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(sA(t), y + 4 + P.font.tick, t > 0 ? `+${t}` : String(t),
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  }
  // Axis for panel B.
  const { major: majB } = axisTicks(0, limB, 4);
  fig.line(bX0, y, bX1, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of majB) {
    if (t < 0 || t > limB) continue;
    fig.line(sB(t), y, sB(t), y + 3, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(sB(t), y + 4 + P.font.tick, String(t),
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  }
  y += 4 + P.font.tick * 2 + 6;
  fig.text((aX0 + aX1) / 2, y, "0 = nearest representative-domain edge",
    { size: "small", anchor: "middle", fill: PALETTE.muted });
  fig.text((bX0 + bX1) / 2, y, "absolute difference",
    { size: "small", anchor: "middle", fill: PALETTE.muted });
  y += P.font.small + 8;

  const classes = [...new Set(rows.flatMap((r) => r.classes))];
  y = fig.legend(P.margin.left, y,
    classes.map((c) => [boundaryClassColour(c), boundaryClassLabel(c)]));
  if (speciesOrderIds.length > 1) {
    const lanes = speciesOrderIds.map((sid, i) => `${i === 0 ? "upper" : "lower"} lane: `
      + speciesTag(sid, speciesName.get(sid))).join(" · ");
    fig.text(P.margin.left, y + P.font.small, `Panel A — ${lanes}`,
      { size: "small", fill: PALETTE.muted });
    y += P.font.small + 3;
  }
  fig.text(P.margin.left, y + P.font.small,
    "Panel A: open marker = distance to a domain start edge, filled = end edge; solid "
    + "connector = supported mapping, dotted = tentative. Panel C: filled square = all "
    + "species agree on the boundary class; bar = mapping confidence; text = species "
    + "observed / analysed.",
    { size: "small", fill: PALETTE.muted });
  y += P.font.small + 3;
  fig.text(P.margin.left, y + P.font.small,
    "A small difference is evidence of a consistent boundary position, not of "
    + "evolutionary conservation.",
    { size: "small", fill: PALETTE.muted });
  return finalise(fig, y + P.font.small + 4);
}

// --------------------------------------------------------------------------- //
// 4J. Comparative synteny
// --------------------------------------------------------------------------- //
/**
 * Local genomic neighbourhood of the target locus, one row per species.
 *
 * This is an adapter, not a renderer: the drawing lives in syntenyFigures.js and
 * is the same code the interactive viewer exports, so the comparative gallery
 * card and the on-screen view cannot drift apart. It accepts the canonical
 * per-species contract rows published in `synteny_neighbourhood`.
 */
export function comparativeSyntenyFigureSpec({
  gene, syntenyNeighbourhood = [], datasetSpecies = [], presetName = "full",
}) {
  return syntenyNeighbourhoodFigureSpec({
    gene,
    rows: comparativeSyntenyRows({ gene, syntenyNeighbourhood, datasetSpecies }),
    presetName,
  });
}

/**
 * Canonical row set for a comparative synteny figure.
 *
 * Every dataset species keeps a row. A species the synteny index does not cover
 * becomes an explicit unresolved row rather than disappearing, because a
 * comparative figure that silently contains fewer species than the dataset
 * misrepresents the comparison.
 */
export function comparativeSyntenyRows({
  gene, syntenyNeighbourhood = [], datasetSpecies = [],
}) {
  const rows = normaliseSyntenyIndex({
    species: syntenyNeighbourhood, gene_symbol: gene,
  });
  const covered = new Set(rows.map((r) => r.speciesId));
  const missing = (datasetSpecies || [])
    .filter((s) => s && !covered.has(s.species_id || s.speciesId))
    .map((s) => unresolvedSpeciesRow({
      speciesId: s.species_id || s.speciesId,
      displayName: s.scientific_name || s.display_species_name
        || s.species_id || s.speciesId,
      gene,
    }));
  return [...rows, ...missing];
}

/** Neighbour-conservation matrix companion, from the same canonical rows. */
export function comparativeSyntenyMatrixFigureSpec({
  gene, syntenyNeighbourhood = [], datasetSpecies = [], presetName = "full",
}) {
  return neighbourConservationMatrixFigureSpec({
    gene,
    rows: comparativeSyntenyRows({ gene, syntenyNeighbourhood, datasetSpecies }),
    presetName,
  });
}

// --------------------------------------------------------------------------- //
// 4K. Isoform-diversity summary
// --------------------------------------------------------------------------- //
/**
 * Three aligned panels per species:
 *
 *   A  protein-model count split into curated and predicted models
 *   B  the annotated protein-length range with the primary protein marked
 *   C  exploratory candidate / variable-block counts
 *
 * Every value is printed next to its mark, so the figure carries the numbers a
 * prose summary would have carried.
 */
export function isoformDiversityFigureSpec({
  gene, rows = [], presetName = "full",
}) {
  const P = preset(presetName);
  const species = speciesOrder(rows.map((r) => ({
    ...r, protein_id: r.primary_protein_id,
  })));
  const fig = createFigure({
    preset: presetName, height: Math.max(1, species.length) * 46 + 260,
  });
  let y = headerBlock(fig, {
    title: gene,
    subtitle: "Comparative isoform diversity · annotated protein models per species",
    question: "How many protein models does each species have, how are they curated, "
      + "how much do their lengths vary, and how many exploratory candidates were "
      + "derived from them?",
  });
  y += 6;
  if (!species.length) {
    fig.text(P.margin.left, y + 8, "No isoform summary rows are available for this run.",
      { size: "label", fill: PALETTE.muted });
    return finalise(fig, y + 24);
  }

  const gutter = speciesGutter(P, species);
  const left = P.margin.left + gutter;
  const usable = P.widthPt - P.margin.right - left;
  const wA = usable * 0.24;
  const wB = usable * 0.48;
  const wC = usable * 0.20;
  const gap = usable * 0.04;
  const aX0 = left;
  const bX0 = aX0 + wA + gap;
  const cX0 = bX0 + wB + gap;

  const maxModels = Math.max(1, ...species.map((r) => Number(r.n_protein_models) || 0));
  const lenLo = Math.min(...species.map((r) => Number(r.protein_length_min) || Infinity));
  const lenHi = Math.max(...species.map((r) => Number(r.protein_length_max) || 0));
  const padLo = Math.max(0, Math.floor((lenLo - (lenHi - lenLo) * 0.08) / 10) * 10);
  const padHi = Math.ceil((lenHi + (lenHi - lenLo) * 0.08) / 10) * 10;
  const sA = (v) => aX0 + (v / maxModels) * wA;
  const sB = (v) => bX0 + ((v - padLo) / Math.max(1, padHi - padLo)) * wB;
  const maxCand = Math.max(1, ...species.map((r) => Number(r.n_exploratory_candidates) || 0));
  const sC = (v) => cX0 + (v / maxCand) * wC;

  fig.text(aX0, y + P.font.small, "A · protein models",
    { size: "small", weight: "bold" });
  fig.text(bX0, y + P.font.small, "B · model length range (aa)",
    { size: "small", weight: "bold" });
  fig.text(cX0, y + P.font.small, "C · exploratory candidates",
    { size: "small", weight: "bold" });
  y += P.font.small + 8;

  // Every panel keeps its numeric read-out on a second text line under its own
  // mark, so a long count never runs into the neighbouring panel and the length
  // range never has to fight the model count for the same strip of paper.
  const laneH = 46;
  const top = y;
  species.forEach((r, i) => {
    const cy = top + i * laneH + 14;
    const readOut = cy + 9 + P.font.small;
    drawSpeciesLabel(fig, left - 8, cy, r);

    // --- A: curated / predicted composition ---------------------------------- //
    const curated = Number(r.n_curated_models) || 0;
    const predicted = Number(r.n_predicted_models) || 0;
    const total = Number(r.n_protein_models) || (curated + predicted);
    fig.rect(aX0, cy - 6, Math.max(0.5, sA(curated) - aX0), 12,
      { fill: FEATURE_STYLES.primary_sequence.fill, stroke: "none" });
    fig.rect(sA(curated), cy - 6, Math.max(0.5, sA(curated + predicted) - sA(curated)), 12,
      { fill: FEATURE_STYLES.alternative_sequence.fill,
        stroke: FEATURE_STYLES.alternative_sequence.stroke, lw: P.lw.thin });
    fig.text(aX0, readOut,
      `${total} model${total === 1 ? "" : "s"}`
      + ` · ${curated} curated${predicted ? ` · ${predicted} predicted` : ""}`,
      { size: "small", fill: PALETTE.muted });

    // --- B: length range with the primary marked ------------------------------ //
    const lo = Number(r.protein_length_min);
    const hi = Number(r.protein_length_max);
    const primary = Number(r.primary_protein_length);
    if (Number.isFinite(lo) && Number.isFinite(hi)) {
      fig.line(sB(lo), cy, sB(hi), cy, { stroke: PALETTE.axis, lw: P.lw.rule });
      fig.line(sB(lo), cy - 4, sB(lo), cy + 4, { stroke: PALETTE.axis, lw: P.lw.rule });
      fig.line(sB(hi), cy - 4, sB(hi), cy + 4, { stroke: PALETTE.axis, lw: P.lw.rule });
      fig.text(sB(lo), readOut, lo === hi ? `${lo} aa` : `${lo}–${hi} aa`,
        { size: "small", anchor: lo === hi ? "middle" : "start", fill: PALETTE.muted });
    }
    if (Number.isFinite(primary)) {
      fig.circle(sB(primary), cy, P.marker,
        { fill: PALETTE.domain, stroke: PALETTE.paper, lw: P.lw.thin });
      fig.text(sB(primary), cy - 8, `primary ${primary}`,
        { size: "small", anchor: "middle", fill: PALETTE.ink });
    }

    // --- C: candidates and variable blocks ------------------------------------ //
    const cand = Number(r.n_exploratory_candidates) || 0;
    fig.rect(cX0, cy - 5, Math.max(0.5, sC(cand) - cX0), 10,
      { fill: FEATURE_STYLES.candidate_region.fill,
        stroke: FEATURE_STYLES.candidate_region.stroke, lw: P.lw.thin });
    const blocks = r.n_variable_alignment_blocks;
    fig.text(cX0, readOut,
      `${cand} candidate${cand === 1 ? "" : "s"}`
      + `${blocks != null ? ` · ${blocks} variable block${blocks === 1 ? "" : "s"}` : ""}`,
      { size: "small", fill: PALETTE.muted });
  });

  y = top + species.length * laneH + 4;
  const { major } = axisTicks(padLo, padHi, 6);
  fig.line(bX0, y, bX0 + wB, y, { stroke: PALETTE.axis, lw: P.lw.rule });
  for (const t of major) {
    if (t < padLo || t > padHi) continue;
    fig.line(sB(t), y, sB(t), y + 3, { stroke: PALETTE.axis, lw: P.lw.rule });
    fig.text(sB(t), y + 4 + P.font.tick, String(t),
      { size: "tick", anchor: "middle", fill: PALETTE.muted });
  }
  y += 4 + P.font.tick * 2 + 6;
  fig.text(bX0 + wB / 2, y, "Annotated protein-model length (aa)",
    { size: "label", anchor: "middle" });
  y += P.font.label + 8;
  y = fig.legend(P.margin.left, y, [
    [FEATURE_STYLES.primary_sequence.fill, "Curated protein model"],
    [FEATURE_STYLES.alternative_sequence.fill, "Predicted protein model"],
    [PALETTE.domain, "Selected primary protein"],
    [FEATURE_STYLES.candidate_region.fill, "Exploratory candidate interval"],
  ]);
  fig.text(P.margin.left, y + P.font.small,
    "Model counts are annotation counts, not validated splice products. Exploratory "
    + "candidates are intervals proposed from within-species isoform comparison and are "
    + "not validated events.",
    { size: "small", fill: PALETTE.muted });
  return finalise(fig, y + P.font.small + 4);
}
