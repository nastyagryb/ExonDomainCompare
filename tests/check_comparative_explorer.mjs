// Behavioural checks for the Comparative Exon–Domain Boundary Explorer's logic.
//
// The filtering, the matrix cell values, the pair-connection rule and the export
// tables live in JavaScript, so they are exercised here in Node against the real
// canonical comparative index rather than re-described in Python. The pytest wrapper
// in test_comparative_boundary_explorer.py runs this file and fails on a non-zero exit.
//
//   node tests/check_comparative_explorer.mjs <coordinateModel.json>

import { readFileSync } from "node:fs";
import {
  EMPTY_FILTERS, activeFilterCount, filterComparativeDataset,
} from "../webapp/frontend/src/pages/viewers/comparativeFilters.js";
import {
  comparativeMatrixFigureSpec, pairedSignedDistanceFigureSpec,
  consistencySummaryFigureSpec, comparativeArchitectureFigureSpec,
  comparativeLongTsv, comparativeMatrixTsv, comparableMappingTsv,
  isSupported, matrixCellFill, speciesTag,
} from "../webapp/frontend/src/pages/viewers/comparativeFigures.js";

const [modelPath] = process.argv.slice(2);
if (!modelPath) {
  console.error("usage: node tests/check_comparative_explorer.mjs <model.json>");
  process.exit(2);
}

const index = JSON.parse(readFileSync(modelPath, "utf8"));
const models = index.models || [];
const multi = index.boundary_dashboard?.multi_species || {};
const gene = index.gene_symbol || "gene";

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
// the dataset under test
// --------------------------------------------------------------------------- //
check("the canonical comparative index is available", multi.available === true);
check("it carries comparable boundary groups",
  (multi.comparable_boundary_groups || []).length > 0);
check("it carries one matrix row per species",
  (multi.boundary_matrix || []).length === models.length,
  `${(multi.boundary_matrix || []).length} rows vs ${models.length} models`);

const unfiltered = filterComparativeDataset(multi, EMPTY_FILTERS);
const allGroups = unfiltered.comparable_boundary_groups;

// --------------------------------------------------------------------------- //
// matrix cells carry the species' own observation, never a fabricated value
// --------------------------------------------------------------------------- //
{
  let mismatched = 0;
  let fabricated = 0;
  for (const row of unfiltered.boundary_matrix) {
    for (const cell of row.cells) {
      const group = allGroups.find(
        (g) => g.comparable_boundary_group_id === cell.comparable_boundary_group_id);
      const obs = (group?.per_species_native_positions || []).find(
        (o) => o.species_id === row.species_id);
      if (obs) {
        if (cell.signed_distance !== obs.signed_distance
          || cell.native_position !== obs.native_position) mismatched += 1;
      } else if (cell.observed || cell.signed_distance != null) {
        // A value in a cell with no backing observation would be invented data.
        fabricated += 1;
      }
    }
  }
  check("every observed cell equals the species' own observation", mismatched === 0,
    `${mismatched} cells disagree with the group detail`);
  check("no cell carries a value without a backing observation", fabricated === 0,
    `${fabricated} fabricated cells`);
}

// A missing observation must not be rendered as a boundary class, and above all not
// as a distance of zero, which is a real class ("sits on the domain edge").
{
  const missingStates = unfiltered.boundary_matrix
    .flatMap((r) => r.cells).filter((c) => !c.observed);
  check("unobserved cells carry no distance",
    missingStates.every((c) => c.signed_distance == null),
    `${missingStates.filter((c) => c.signed_distance != null).length} unobserved cells have a value`);
  check("unobserved cells are painted in a non-class fill",
    missingStates.every((c) => ["#eef1f4", "#f6efe0"].includes(matrixCellFill(c.state))));
}

// --------------------------------------------------------------------------- //
// pair connections: only supported mappings may be drawn as confirmed pairs
// --------------------------------------------------------------------------- //
{
  const tentative = allGroups.filter((g) => !isSupported(g.mapping_status));
  check("the dataset contains at least one tentative group to test against",
    tentative.length > 0);
  check("no tentative group is marked connectable",
    tentative.every((g) => g.connectable === false),
    `${tentative.filter((g) => g.connectable).length} tentative groups are connectable`);
  const supportedPairs = allGroups.filter(
    (g) => isSupported(g.mapping_status) && g.per_species_native_positions.length >= 2);
  check("supported groups with two observations are connectable",
    supportedPairs.every((g) => g.connectable === true));
}

// --------------------------------------------------------------------------- //
// filters feed every view from one dataset
// --------------------------------------------------------------------------- //
{
  const cls = multi.filter_options?.boundary_classes?.[0];
  const filtered = filterComparativeDataset(multi, { boundaryClasses: [cls] });
  const ids = new Set(filtered.comparable_boundary_groups.map(
    (g) => g.comparable_boundary_group_id));

  check("a class filter reduces the group set",
    ids.size > 0 && ids.size < allGroups.length,
    `${ids.size} of ${allGroups.length} groups kept for class ${cls}`);
  check("the matrix shows exactly the filtered groups",
    filtered.boundary_matrix.every((r) => r.cells.every(
      (c) => ids.has(c.comparable_boundary_group_id))));
  check("the consistency summary shows exactly the filtered groups",
    filtered.distance_statistics.every((s) => ids.has(s.comparable_boundary_group_id))
    && filtered.distance_statistics.length === ids.size);
  check("the inspection cases are restricted to the filtered groups",
    filtered.inspection_cases.every((c) => ids.has(c.comparable_boundary_group_id)));
  check("the exported long table contains only filtered observations", (() => {
    const body = comparativeLongTsv(filtered.comparable_boundary_groups)
      .trim().split("\n").slice(1);
    return body.length === filtered.counts.visible_observations;
  })());
  check("the exported matrix table has one column per filtered group", (() => {
    const head = comparativeMatrixTsv(filtered.boundary_matrix,
      filtered.comparable_boundary_groups).split("\n")[0].split("\t");
    return head.length === ids.size + 2;
  })());
}

// A species filter must not silently drop the other species' rows without saying so.
{
  const sid = models[0].species_id;
  const one = filterComparativeDataset(multi, { species: [sid] });
  check("a species filter keeps only that species' matrix row",
    one.boundary_matrix.length === 1 && one.boundary_matrix[0].species_id === sid);
  check("filtering to one species leaves no connectable pairs",
    one.comparable_boundary_groups.every((g) => g.connectable === false),
    "a single observation must never be presented as a confirmed cross-species pair");
}

// Hiding unmapped cells must not turn a filtered-out observation into "absent".
{
  const cls = multi.filter_options?.boundary_classes?.[0];
  const kept = filterComparativeDataset(multi, { boundaryClasses: [cls] });
  const filteredOut = kept.boundary_matrix.flatMap((r) => r.cells)
    .filter((c) => c.state === "filtered_out");
  check("an observation removed by a filter is labelled filtered_out, not absent",
    filteredOut.length > 0,
    "expected at least one hidden-by-filter cell to be distinguishable");
}

check("the active-filter count reflects the chosen filters",
  activeFilterCount({ ...EMPTY_FILTERS, species: ["x"], exactNearOnly: true }) === 2);
check("no filters means no active filters",
  activeFilterCount(EMPTY_FILTERS) === 0);

// --------------------------------------------------------------------------- //
// selected-group detail
// --------------------------------------------------------------------------- //
{
  const g = allGroups[0];
  const stat = unfiltered.distance_statistics.find(
    (s) => s.comparable_boundary_group_id === g.comparable_boundary_group_id);
  const required = ["species_id", "protein_id", "exon_transition", "native_position",
    "msa_column", "nearest_edge", "signed_distance", "boundary_class",
    "mapping_method", "mapping_confidence"];
  check("every detail row carries the full observation contract",
    g.per_species_native_positions.every(
      (o) => required.every((k) => Object.prototype.hasOwnProperty.call(o, k))),
    `missing: ${required.filter((k) => !(k in g.per_species_native_positions[0]))}`);
  check("the two-species statistic reports the raw pair",
    stat.primary_statistic === "raw_pair"
    && (stat.raw_signed_distances || []).length === 2);
  check("the statistic reports the cross-species difference",
    stat.cross_species_difference != null);
}

// --------------------------------------------------------------------------- //
// inspection cases point at real groups and stay non-judgemental
// --------------------------------------------------------------------------- //
{
  const cases = multi.inspection_cases || [];
  const ids = new Set(allGroups.map((g) => g.comparable_boundary_group_id));
  check("there is at least one inspection case", cases.length > 0);
  check("every inspection case points at a real comparable group",
    cases.every((c) => ids.has(c.comparable_boundary_group_id)));
  check("no inspection case calls a discrepancy an error",
    cases.every((c) => !/\berror\b|\bbug\b|\bwrong\b|\bfail/i.test(
      `${c.label} ${c.detail}`)),
    cases.map((c) => c.label).join(", "));
}

// --------------------------------------------------------------------------- //
// exports really render
// --------------------------------------------------------------------------- //
{
  const stats = unfiltered.distance_statistics;
  const figures = {
    matrix: comparativeMatrixFigureSpec({
      gene, matrix: unfiltered.boundary_matrix, groups: allGroups }),
    paired: pairedSignedDistanceFigureSpec({ gene, groups: allGroups, stats }),
    consistency: consistencySummaryFigureSpec({ gene, stats, groups: allGroups }),
    architecture: comparativeArchitectureFigureSpec({
      gene, group: allGroups[0], models }),
  };
  for (const [name, fig] of Object.entries(figures)) {
    const svg = fig.toSvg();
    const pdf = fig.toPdf();
    check(`${name} figure draws marks`, fig.marks.length > 0);
    check(`${name} figure produces standalone SVG`,
      svg.startsWith("<?xml") || svg.startsWith("<svg"));
    check(`${name} figure has no layout warnings`, fig.warnings.length === 0,
      fig.warnings.join("; "));
    check(`${name} figure produces a vector PDF`,
      pdf.length > 800 && new TextDecoder().decode(pdf.slice(0, 8)).startsWith("%PDF-"));
    // A CSS-dependent SVG renders black once it leaves the browser.
    check(`${name} SVG carries explicit fills`,
      !svg.includes("class=\"") && /fill="#/.test(svg));
  }

  // Every matrix mode must be exportable, since the toggle changes what is shown.
  for (const mode of ["signed", "absolute", "class"]) {
    const fig = comparativeMatrixFigureSpec({
      gene, matrix: unfiltered.boundary_matrix, groups: allGroups, mode });
    check(`matrix exports in ${mode} mode`, fig.marks.length > 0 && fig.warnings.length === 0);
  }

  const mapping = comparableMappingTsv(allGroups, stats);
  check("the mapping table has one row per comparable group",
    mapping.trim().split("\n").length === allGroups.length + 1);
  check("the mapping table records how each group was matched",
    mapping.split("\n")[0].includes("mapping_method")
    && mapping.includes("msa_aligned_position"));
  check("the mapping table never claims exon rank as evidence",
    !/exon_rank|exon_number_only/.test(mapping));
}

check("species tags abbreviate the genus", speciesTag("gallus_gallus", "Gallus gallus")
  === "G. gallus");

console.log(failures === 0
  ? "OK — comparative explorer behaviour verified"
  : `${failures} check(s) failed`);
process.exit(failures === 0 ? 0 : 1);
