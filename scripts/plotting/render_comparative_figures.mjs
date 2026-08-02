// Render the comparative multi-species boundary figures from the coordinate model.
//
// Like render_main_figures.mjs, this imports the very same builders the browser
// exports through (comparativeFigures.js) and feeds them the very same canonical
// comparative index the React explorer reads (boundary_dashboard.multi_species). A
// figure produced here and the same figure exported from the Explorer are therefore
// one implementation, not two that have to be kept in agreement by hand.
//
// Output per figure: standalone SVG plus true vector PDF. There is no HTML, no CSS and
// no screenshot step anywhere in this path.
//
// Usage:
//   node scripts/plotting/render_comparative_figures.mjs <coordinateModel.json> <outDir>
//        [--selected-group=CBG3] [--mode=signed|absolute|class]

import { mkdirSync, readFileSync, writeFileSync } from "node:fs";
import { join } from "node:path";
import {
  comparativeMatrixFigureSpec, pairedSignedDistanceFigureSpec,
  consistencySummaryFigureSpec, comparativeArchitectureFigureSpec,
  comparativeLongTsv, comparativeMatrixTsv, comparableMappingTsv,
} from "../../webapp/frontend/src/pages/viewers/comparativeFigures.js";

const args = process.argv.slice(2);
const [modelPath, outDir] = args.filter((a) => !a.startsWith("--"));
const flag = (name, fallback = null) => {
  const hit = args.find((a) => a.startsWith(`--${name}=`));
  return hit ? hit.split("=").slice(1).join("=") : fallback;
};
if (!modelPath || !outDir) {
  console.error("usage: node scripts/plotting/render_comparative_figures.mjs "
    + "<model.json> <outDir> [--selected-group=CBG1] [--mode=signed|absolute|class]");
  process.exit(2);
}

const index = JSON.parse(readFileSync(modelPath, "utf8"));
const models = index.models || [];
const dash = index.boundary_dashboard || {};
const multi = dash.multi_species || {};
const gene = index.gene_symbol || models[0]?.gene_symbol || "gene";
const groups = multi.comparable_boundary_groups || [];
const stats = multi.distance_statistics || [];
const matrix = multi.boundary_matrix || [];
const threshold = dash.near_edge_threshold_aa ?? 5;
const mode = flag("mode", "signed");

if (!multi.available || !groups.length) {
  // Refusing to render is the correct outcome: an empty comparative figure would be
  // read as "no boundaries are comparable", which is a scientific claim this run has
  // no evidence for either way.
  console.error("no comparable boundary groups in this model — nothing to render");
  process.exit(1);
}
mkdirSync(outDir, { recursive: true });

const written = [];
function emit(stem, fig, { tsv } = {}) {
  writeFileSync(join(outDir, `${stem}.svg`), fig.toSvg(), "utf8");
  writeFileSync(join(outDir, `${stem}.pdf`), Buffer.from(fig.toPdf()));
  if (tsv) writeFileSync(join(outDir, `${stem}.tsv`), tsv, "utf8");
  written.push({
    stem, width: fig.width, height: fig.height,
    marks: fig.marks.length, warnings: fig.warnings, has_tsv: Boolean(tsv),
  });
}

// A selected group is an on-demand export driven by a user selection, so it defaults to
// the first group only to keep the renderer usable from the command line.
const selectedGroupId = flag("selected-group", groups[0].comparable_boundary_group_id);
const selectedGroup = groups.find(
  (g) => g.comparable_boundary_group_id === selectedGroupId) || groups[0];

emit("comparative_boundary_matrix",
  comparativeMatrixFigureSpec({ gene, matrix, groups, mode, nearEdgeThreshold: threshold }),
  { tsv: comparativeMatrixTsv(matrix, groups, mode) });

emit("comparative_paired_signed_distance",
  pairedSignedDistanceFigureSpec({ gene, groups, stats, nearEdgeThreshold: threshold }),
  { tsv: comparativeLongTsv(groups) });

emit("comparative_boundary_consistency",
  consistencySummaryFigureSpec({ gene, stats, groups }),
  { tsv: comparableMappingTsv(groups, stats) });

emit(`ondemand_comparative_architecture_${selectedGroup.comparable_boundary_group_id}`,
  comparativeArchitectureFigureSpec({ gene, group: selectedGroup, models }));

writeFileSync(join(outDir, "comparative_render_summary.json"),
  JSON.stringify(written, null, 2));
for (const w of written) {
  console.log(`  ${w.stem}: ${w.width}x${w.height}pt, ${w.marks} marks`
    + (w.warnings.length ? ` (${w.warnings.length} layout warnings)` : ""));
}
console.log(`OK — rendered ${written.length} comparative figure(s) to ${outDir}`);
