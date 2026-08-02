// Behavioural checks for the shared synteny display model and renderer.
//
// The layout decisions that broke the old view live in JavaScript — how many
// loci are shown, where the target sits, whether anything falls outside the
// drawing area — so they are exercised here in Node against a real index rather
// than re-described in Python. The pytest wrapper in
// test_shared_synteny_contract.py runs this file and fails on a non-zero exit.
//
//   node tests/check_synteny_renderer.mjs <synteny_locus_index.json>

import { readFileSync } from "node:fs";
import {
  normaliseSyntenyIndex, normaliseSpeciesRow, slotGrid, legendEntries,
  syntenyRowsTsv, isPlaceholderLocus,
} from "../webapp/frontend/src/pages/viewers/syntenyModel.js";
import {
  syntenyNeighbourhoodFigureSpec, neighbourConservationMatrixFigureSpec,
} from "../webapp/frontend/src/pages/viewers/syntenyFigures.js";

const [indexPath] = process.argv.slice(2);
if (!indexPath) {
  console.error("usage: node tests/check_synteny_renderer.mjs <synteny_locus_index.json>");
  process.exit(2);
}

const index = JSON.parse(readFileSync(indexPath, "utf8"));
const gene = index.gene_symbol || "gene";
const rows = normaliseSyntenyIndex(index);

let failures = 0;
const check = (name, ok, detail = "") => {
  if (ok) {
    console.log(`  ok    ${name}`);
  } else {
    failures += 1;
    console.log(`  FAIL  ${name}${detail ? ` — ${detail}` : ""}`);
  }
};

// --------------------------------------------------------------------------- //
// the display model
// --------------------------------------------------------------------------- //
check("the index normalises to at least one species row", rows.length > 0);

for (const row of rows) {
  const flanking = [...row.upstream, ...row.downstream];
  check(`${row.speciesId}: the target is not one of the neighbours`,
    flanking.every((n) => !n.is_target));
  check(`${row.speciesId}: exactly one target locus`,
    row.loci.filter((n) => n.is_target).length === 1);
  check(`${row.speciesId}: displayed counts match the rendered loci`,
    row.upstream.length === row.counts.displayedUpstream
    && row.downstream.length === row.counts.displayedDownstream);
  check(`${row.speciesId}: loci are in slot order`,
    row.loci.map((n) => n.slot_x).every((v, i, a) => i === 0 || a[i - 1] < v));
  check(`${row.speciesId}: the counts label states what is shown`,
    row.countsLabel.startsWith(`${row.upstream.length + row.downstream.length} `)
    || row.countsLabel === "No flanking loci available");
}

// --------------------------------------------------------------------------- //
// the slot grid: the target is centred and nothing is clipped
// --------------------------------------------------------------------------- //
const grid = slotGrid(rows);
check("the grid has an odd column count so a centre column exists",
  grid.columns % 2 === 1);
check("the target column is the exact centre",
  grid.targetColumn === (grid.columns - 1) / 2);

for (const row of rows) {
  const columns = row.loci.map((n) => grid.columnOf(n));
  check(`${row.speciesId}: every locus lands inside the grid`,
    columns.every((c) => c >= 0 && c < grid.columns),
    `columns ${columns.join(",")} of ${grid.columns}`);
  check(`${row.speciesId}: the target lands on the centre column`,
    grid.columnOf(row.target) === grid.targetColumn);
  check(`${row.speciesId}: no two loci share a column`,
    new Set(columns).size === columns.length);
}

// A species with fewer real neighbours must not shift the target off centre.
const short = normaliseSpeciesRow({
  species_id: "short_species",
  display_species_name: "Short species",
  gene_symbol: gene,
  loci: [
    ...[5, 4, 3, 2, 1].map((r) => ({
      slot_x: -r, side: "upstream", rank: r, symbol: `UP${r}`,
      orthology_class: "exact", strand: "+",
    })),
    { slot_x: 0, side: "target", rank: 0, symbol: gene, is_target: true, strand: "+" },
    ...[1, 2, 3, 4].map((r) => ({
      slot_x: r, side: "downstream", rank: r, symbol: `DN${r}`,
      orthology_class: "exact", strand: "-",
    })),
  ],
});
const mixedGrid = slotGrid([short, ...rows]);
check("an unequal neighbourhood keeps the target on the centre column",
  mixedGrid.columnOf(short.target) === mixedGrid.targetColumn);
check("an unequal neighbourhood reports nine loci, not ten",
  short.countsLabel === "9 flanking loci shown · 5 upstream · 4 downstream",
  short.countsLabel);
check("an unequal neighbourhood is flagged as incomplete rather than padded",
  short.truncationStatus === "fewer_available" && short.loci.length === 10);

// --------------------------------------------------------------------------- //
// legend and placeholders
// --------------------------------------------------------------------------- //
const present = new Set(rows.flatMap((r) => r.classesPresent));
const legend = legendEntries(rows);
check("the legend lists only classes actually drawn",
  legend.every((e) => present.has(e.cls)) && legend.length === present.size);
check("every legend entry carries a definition",
  legend.every((e) => e.definition && e.definition.length > 10));
check("LOC identifiers are recognised as placeholder labels",
  isPlaceholderLocus("LOC121107413") && !isPlaceholderLocus("TACC1"));

// --------------------------------------------------------------------------- //
// exports draw exactly the displayed set
// --------------------------------------------------------------------------- //
const fig = syntenyNeighbourhoodFigureSpec({ gene, rows });
check("the figure renders without a layout warning", fig.warnings.length === 0,
  fig.warnings.join("; "));

const svg = fig.toSvg();
check("the exported SVG is self-contained (no CSS variables)",
  !svg.includes("var(--"));
for (const row of rows) {
  for (const locus of row.loci) {
    // Labels may be truncated to fit a slot, so match on the stable prefix.
    const stem = String(locus.symbol).slice(0, 4);
    check(`${row.speciesId}: ${locus.symbol} appears in the exported figure`,
      svg.includes(stem), stem);
  }
}

const tsvText = syntenyRowsTsv(rows);
const tsvLines = tsvText.trim().split("\n").slice(1);
const renderedLoci = rows.reduce((acc, r) => acc + r.loci.length, 0);
check("the source TSV has one row per rendered locus",
  tsvLines.length === renderedLoci, `${tsvLines.length} vs ${renderedLoci}`);
check("the source TSV marks the target rows",
  tsvLines.filter((l) => l.split("\t")[4] === "true").length === rows.length);

if (rows.length > 1) {
  const matrix = neighbourConservationMatrixFigureSpec({ gene, rows });
  check("the conservation matrix renders without a layout warning",
    matrix.warnings.length === 0, matrix.warnings.join("; "));
}

// A single-species run renders the same way, with no comparative wording.
const oneSvg = syntenyNeighbourhoodFigureSpec({ gene, rows: [rows[0]] }).toSvg();
check("a single-species figure makes no cross-species claim",
  !oneSvg.includes("One row per species"));

console.log(failures ? `\n${failures} check(s) failed` : "\nall checks passed");
process.exit(failures ? 1 : 0);
