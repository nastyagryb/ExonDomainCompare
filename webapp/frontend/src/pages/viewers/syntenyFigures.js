// Publication figures for local synteny, built from the shared display model.
//
// The interactive viewer and these exports read the same normalised rows from
// syntenyModel.js, so a figure in the paper shows exactly the loci a reader saw
// on screen, in the same order, with the same styles and the same target slot.
// Nothing here re-derives orthology, ranks or the target position.

import { createFigure, preset, PALETTE, textWidth } from "./figureSpec.js";
import { headerBlock, finalise, speciesLabel } from "./mainFigures.js";
import {
  legendEntries, orthologyStyle, slotGrid, syntenyCoverage, syntenyRowsTsv,
} from "./syntenyModel.js";

export { syntenyRowsTsv };

function speciesGutter(P, rows) {
  return Math.max(96, ...rows.map(
    (r) => textWidth(speciesLabel(r.displayName), P.font.label) + 12));
}

/**
 * Draw a note across the full text column, wrapping instead of running past the
 * right margin. Returns the y below the last line.
 */
function note(fig, x, y, text, width) {
  const P = fig.preset;
  const pt = P.font.small;
  const words = String(text).split(/\s+/);
  let line = "";
  let cursor = y;
  for (const word of words) {
    const next = line ? `${line} ${word}` : word;
    if (line && textWidth(next, pt) > width) {
      fig.text(x, cursor + pt, line, { size: "small", fill: PALETTE.muted });
      cursor += pt + 2;
      line = word;
    } else {
      line = next;
    }
  }
  if (line) {
    fig.text(x, cursor + pt, line, { size: "small", fill: PALETTE.muted });
    cursor += pt + 2;
  }
  return cursor;
}

/** Shorten a label until it fits its slot; the full symbol stays in the table. */
function fitLabel(text, width, pt) {
  let label = String(text || "");
  if (textWidth(label, pt) <= width) return label;
  while (label.length > 2 && textWidth(`${label}…`, pt) > width) {
    label = label.slice(0, -1);
  }
  return `${label}…`;
}

/** One locus: a gene box plus an arrowhead on the annotated strand. */
function drawLocus(fig, { x, w, cy, locus, P, labelDy }) {
  const style = orthologyStyle(locus.orthology_class);
  const h = locus.is_target ? 17 : 13;
  fig.rect(x, cy - h / 2, w, h, {
    fill: style.fill, stroke: style.stroke,
    lw: locus.is_target ? P.lw.outline : P.lw.thin,
  });
  if (locus.strand) {
    // Arrowhead inside the box, so transcription direction survives export
    // without needing a polygon mark the PDF writer does not have.
    const dir = locus.strand === "-" ? -1 : 1;
    const tipX = dir > 0 ? x + w - 1.5 : x + 1.5;
    const back = tipX - dir * Math.min(4.5, w * 0.3);
    const colour = locus.is_target ? PALETTE.paper : style.stroke;
    fig.line(tipX, cy, back, cy - 3.2, { stroke: colour, lw: P.lw.thin });
    fig.line(tipX, cy, back, cy + 3.2, { stroke: colour, lw: P.lw.thin });
  }
  if (labelDy) {
    fig.line(x + w / 2, cy + h / 2, x + w / 2, cy + h / 2 + labelDy - P.font.small,
      { stroke: PALETTE.axis, lw: P.lw.thin });
  }
  fig.text(x + w / 2, cy + h / 2 + P.font.small + 2 + (labelDy || 0),
    fitLabel(locus.symbol, w + 4, P.font.small), {
      size: "small", anchor: "middle",
      fill: locus.is_target ? PALETTE.ink : style.text,
      weight: locus.is_target ? "bold" : "normal",
      italic: Boolean(style.italic),
    });
}

/**
 * Target-centred local neighbourhood, one row per species.
 *
 * Every row uses the same slot grid, so the target column lines up down the
 * figure and a species with fewer real neighbours simply leaves the outer slot
 * empty. No locus is invented and none is dropped.
 */
export function syntenyNeighbourhoodFigureSpec({
  gene, rows = [], presetName = "full", title, subtitle, question,
}) {
  const P = preset(presetName);
  const species = rows.filter(Boolean);
  const fig = createFigure({
    preset: presetName, height: Math.max(1, species.length) * 46 + 250,
  });
  let y = headerBlock(fig, {
    title: title || `${gene} · local genomic neighbourhood`,
    subtitle: subtitle || (species.length > 1
      ? "One row per species, target locus centred on a shared slot grid"
      : species[0]?.countsLabel || ""),
    question: question || "Which annotated loci flank the target gene, in which "
      + "transcription direction, and are the same neighbours present across species?",
  });
  y += 8;
  if (!species.length) {
    fig.text(P.margin.left, y + 8,
      "No local gene neighbourhood is available for this run.",
      { size: "label", fill: PALETTE.muted });
    return finalise(fig, y + 24);
  }

  const grid = slotGrid(species);
  const gutter = speciesGutter(P, species);
  const x0 = P.margin.left + gutter;
  const x1 = P.widthPt - P.margin.right;
  const colW = (x1 - x0) / grid.columns;
  // When a symbol cannot fit its slot even truncated, labels alternate between
  // two rows with a short leader, so neighbouring labels never run together.
  const stagger = species.some((r) => r.loci.some(
    (n) => textWidth(String(n.symbol), P.font.small) > colW - 3));
  const laneH = stagger ? 46 + P.font.small : 46;
  const top = y;
  const targetX = x0 + grid.targetColumn * colW;

  fig.rect(targetX, top - 6, colW, species.length * laneH + 4,
    { fill: PALETTE.domain, stroke: "none", opacity: 0.07 });

  species.forEach((row, i) => {
    const cy = top + i * laneH + 10;
    fig.text(x0 - 8, cy + 3, speciesLabel(row.displayName),
      { size: "label", anchor: "end", italic: true });
    fig.line(x0, cy, x1, cy, { stroke: PALETTE.grid, lw: P.lw.thin });
    for (const locus of row.loci) {
      const col = grid.columnOf(locus);
      drawLocus(fig, {
        x: x0 + col * colW + 1.5, w: colW - 3, cy, locus, P,
        labelDy: stagger && col % 2 === 1 ? P.font.small + 2 : 0,
      });
    }
  });

  y = top + species.length * laneH + 4;
  fig.text(x0, y + P.font.small, `upstream (rank ${grid.perSide} → 1)`,
    { size: "small", fill: PALETTE.muted });
  fig.text(targetX + colW / 2, y + P.font.small, "target locus",
    { size: "small", anchor: "middle", fill: PALETTE.muted });
  fig.text(x1, y + P.font.small, `downstream (rank 1 → ${grid.perSide})`,
    { size: "small", anchor: "end", fill: PALETTE.muted });
  y += P.font.small + 10;

  // Only the classes actually drawn appear in the legend.
  y = fig.legend(P.margin.left, y,
    legendEntries(species).map((e) => [e.fill, e.label]));
  const textWidthAvailable = P.widthPt - P.margin.left - P.margin.right;
  y = note(fig, P.margin.left, y,
    "Loci are placed by annotated neighbour rank, not by genomic coordinate; the "
    + "spacing carries no distance information. The arrowhead is the annotated "
    + "transcription direction. Local synteny is supporting evidence for locus "
    + "identity and genomic context; it does not by itself establish isoform or "
    + "event identity.", textWidthAvailable);
  if (species.length > 1) {
    const cov = syntenyCoverage(species);
    y = note(fig, P.margin.left, y,
      `Species shown: ${cov.shown} of ${cov.requested} in this dataset · `
      + `complete neighbourhood: ${cov.complete.length ? cov.complete.join(", ") : "none"} · `
      + `partial: ${cov.partial.length ? cov.partial.join(", ") : "none"} · `
      + `no neighbourhood available: `
      + `${cov.unavailable.length ? cov.unavailable.join(", ") : "none"}.`,
      textWidthAvailable);
  }
  const incomplete = species.filter((r) => r.truncationStatus === "fewer_available");
  if (incomplete.length) {
    const detail = incomplete.length === 1
      ? `${speciesLabel(incomplete[0].displayName)}: ${incomplete[0].countsLabel}.`
      : incomplete.map((r) => `${speciesLabel(r.displayName)}: `
        + `${r.counts.displayedUpstream} upstream, `
        + `${r.counts.displayedDownstream} downstream`).join("; ") + ".";
    y = note(fig, P.margin.left, y,
      `${detail} Those slots are left empty rather than filled; not every assembly `
      + `provides ${grid.perSide} annotated loci on both sides.`, textWidthAvailable);
  }
  return finalise(fig, y + 4);
}

/**
 * Neighbour-conservation matrix: which flanking symbols recur across species.
 *
 * A cell is filled only where that species really carries the symbol at some
 * rank; a species that never resolved the locus stays an explicit blank.
 */
export function neighbourConservationMatrixFigureSpec({
  gene, rows = [], presetName = "full",
}) {
  const P = preset(presetName);
  const species = rows.filter(Boolean);
  const norm = (s) => String(s || "").toLowerCase();

  const bySymbol = new Map();
  for (const row of species) {
    for (const locus of row.loci) {
      if (locus.is_target || locus.placeholder || !locus.symbol) continue;
      const key = norm(locus.symbol);
      if (!bySymbol.has(key)) bySymbol.set(key, { label: locus.symbol, hits: new Map() });
      bySymbol.get(key).hits.set(row.speciesId, locus);
    }
  }
  // Reading order: most widely shared first, then alphabetically.
  const groups = [...bySymbol.values()]
    .filter((g) => g.hits.size > 1 || species.length === 1)
    .sort((a, b) => b.hits.size - a.hits.size || a.label.localeCompare(b.label));

  const cell = 18;
  const gutter = speciesGutter(P, species);
  const fig = createFigure({
    preset: presetName, height: species.length * cell + 260,
  });
  let y = headerBlock(fig, {
    title: `${gene} · neighbour conservation`,
    subtitle: "Flanking symbols recovered per species, canonical species order",
    question: "Which flanking loci recur across the analysed species, and where is "
      + "the neighbourhood incomplete?",
  });
  y += 8;
  if (!groups.length || !species.length) {
    fig.text(P.margin.left, y + 8,
      "No flanking symbol is shared by more than one species in this dataset.",
      { size: "label", fill: PALETTE.muted });
    return finalise(fig, y + 24);
  }

  const x0 = P.margin.left + gutter;
  const colW = Math.max(14, Math.min(34,
    (P.widthPt - P.margin.right - x0) / groups.length));
  // Column labels alternate across two rows so long symbols stay readable
  // without rotated text, which the PDF writer does not carry.
  const top = y + 2 * (P.font.small + 3) + 6;
  groups.forEach((g, j) => {
    const cx = x0 + j * colW + colW / 2;
    const row = j % 2;
    const ly = y + P.font.small + row * (P.font.small + 3);
    fig.text(cx, ly, fitLabel(g.label, colW * 2 - 2, P.font.small),
      { size: "small", anchor: "middle", fill: PALETTE.ink });
    if (row === 1) {
      fig.line(cx, ly + 2, cx, top - 2, { stroke: PALETTE.grid, lw: P.lw.thin });
    }
  });
  species.forEach((row, i) => {
    const cy = top + i * cell;
    fig.text(x0 - 8, cy + cell * 0.68, speciesLabel(row.displayName),
      { size: "label", anchor: "end", italic: true });
    groups.forEach((g, j) => {
      const locus = g.hits.get(row.speciesId);
      const style = locus ? orthologyStyle(locus.orthology_class) : null;
      fig.rect(x0 + j * colW + 1, cy + 1, colW - 2, cell - 2, {
        fill: style ? style.fill : PALETTE.paper,
        stroke: style ? style.stroke : PALETTE.grid, lw: P.lw.thin,
      });
    });
  });
  y = top + species.length * cell + 10;
  // Only the classes that a cell actually carries; placeholder loci never enter
  // the matrix, so their swatch must not appear in its legend.
  const cellClasses = new Set(groups.flatMap(
    (g) => [...g.hits.values()].map((l) => l.orthology_class)));
  y = fig.legend(P.margin.left, y,
    [...legendEntries(species).filter((e) => cellClasses.has(e.cls))
      .map((e) => [e.fill, e.label]),
    [PALETTE.paper, "Symbol not recovered in this species"]]);
  y = note(fig, P.margin.left, y,
    "A filled cell means the annotation of that species carries the same gene "
    + "symbol somewhere in the flanking window. A shared symbol is a nomenclature "
    + "match supporting locus context, not an orthology assignment.",
    P.widthPt - P.margin.left - P.margin.right);
  return finalise(fig, y + 4);
}
