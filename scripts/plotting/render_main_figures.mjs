// Render the single-species main figures from the validated coordinate model.
//
// This is the production renderer for the Figure Gallery's main figures. It
// imports the very same adapter (`figureData.js`) and figure builders
// (`mainFigures.js`) that the Gene Explorer's interactive views export through,
// and it is handed the same model objects the React components receive. A Gallery
// figure and the corresponding Gene Explorer figure are therefore the same figure,
// not two implementations that have to be kept in agreement by hand.
//
// Output per figure: standalone SVG, true vector PDF, and a source TSV. PNG is
// rasterised by the calling Python stage, which owns the 300 dpi conversion.
//
// Usage:
//   node scripts/plotting/render_main_figures.mjs <coordinateModel.json> <outDir>

import { mkdirSync, readFileSync, writeFileSync } from "node:fs";
import { join } from "node:path";
import {
  boundaryClassSummaryFigure, boundaryFigure, boundaryTsv, domainArchitectureFigure,
  domainArchitectureTsv, exonMapFigure, exonMapTsv, signedDistanceFigure,
} from "../../webapp/frontend/src/pages/viewers/figureData.js";

const args = process.argv.slice(2);
const [modelPath, outDir] = args.filter((a) => !a.startsWith("--"));
// The selected-boundary figure is an on-demand export, never a permanent Gallery
// card: a Gallery figure has no user selection to depict.
const selectFlag = args.find((a) => a.startsWith("--selected-boundary="));
const selectedBoundary = selectFlag ? selectFlag.split("=").slice(1).join("=") : null;
if (!modelPath || !outDir) {
  console.error("usage: node scripts/plotting/render_main_figures.mjs <model.json> <outDir> "
    + "[--selected-boundary=<transition|boundary id>]");
  process.exit(2);
}

const model = JSON.parse(readFileSync(modelPath, "utf8"));
const models = model.models || [];
if (!models.length) {
  console.error("no models in the coordinate model");
  process.exit(1);
}
mkdirSync(outDir, { recursive: true });

/** Only classified boundaries carry a distance; pre-cluster runs have none. */
const classifiedBoundaries = (m) =>
  (m.exon_boundaries || []).filter((b) => b.signed_distance != null);

const written = [];

function emit(stem, fig, { tsv, model, figureType } = {}) {
  writeFileSync(join(outDir, `${stem}.svg`), fig.toSvg(), "utf8");
  writeFileSync(join(outDir, `${stem}.pdf`), Buffer.from(fig.toPdf()));
  if (tsv) writeFileSync(join(outDir, `${stem}.tsv`), tsv, "utf8");
  written.push({
    stem, width: fig.width, height: fig.height,
    marks: fig.marks.length, warnings: fig.warnings,
    has_tsv: Boolean(tsv),
    // Recorded so nothing downstream has to recover the identity by parsing the
    // file name. These files are derived assets; the catalogue that decides which
    // of them a reader ever sees reads this summary, not the directory listing.
    figure_type: figureType || "",
    model_id: model?.model_id || "",
    model_role: model?.model_role || "",
    species_id: model?.species_id || "",
    isoform: model?.final_isoform_label || model?.isoform || "",
    is_primary_reference: Boolean(model?.is_primary_reference),
  });
}

// Every model must state its own identity. A renderer that works out which protein
// it is drawing by counting entries or by parsing a file name draws the right
// picture only until the input order changes, and then fails silently.
for (const m of models) {
  if (!m.model_id || !m.model_role) {
    console.error(`model for ${m.species_id || "?"} / ${m.protein_id || "?"} is `
      + "missing model_id or model_role; identity must be explicit");
    process.exit(1);
  }
}

// The file name is *derived from* that identity, never the source of it. A species'
// primary reference keeps the plain species stem, which is what every single-primary
// gene has always produced; a further model of the same species carries its role, so
// FGFR2's IIIc figures cannot overwrite its IIIb figures.
function figureKey(m) {
  const sp = m.species_id || "sp";
  if (m.is_primary_reference) return sp;
  const suffix = m.isoform || m.final_isoform_label || m.model_role;
  return `${sp}_${String(suffix).replace(/[^\w.-]+/g, "_")}`;
}

for (const m of models) {
  const sp = figureKey(m);
  const boundaries = classifiedBoundaries(m);
  const hasDomains = (m.representative_domains || []).length > 0;
  const available = m.status === "available";

  // The exon-to-protein projection needs no domain layer, so it is also the one
  // main figure a pre-cluster run can legitimately produce.
  emit(`main_${sp}_primary_exon_projection`, exonMapFigure(m),
    { tsv: exonMapTsv(m), model: m, figureType: "primary_exon_projection" });

  if (!available || !hasDomains) {
    console.log(`${sp}: domain layer pending — only the exon projection was rendered`);
    continue;
  }

  emit(`main_${sp}_integrated_domain_architecture`, domainArchitectureFigure(m),
    { tsv: domainArchitectureTsv(m), model: m,
      figureType: "integrated_domain_architecture" });

  if (!boundaries.length) {
    console.log(`${sp}: no classified exon boundaries — boundary figures skipped`);
    continue;
  }

  // All three boundary figures are views of the same classified boundaries, so
  // each ships that table: a reader who downloads one figure gets the numbers
  // behind it without having to find a different card first.
  const bndTsv = boundaryTsv(boundaries);
  emit(`main_${sp}_boundary_on_architecture`, boundaryFigure(m, boundaries),
    { tsv: bndTsv, model: m, figureType: "boundary_on_architecture" });
  emit(`main_${sp}_signed_boundary_distances`, signedDistanceFigure(m, boundaries),
    { tsv: bndTsv, model: m, figureType: "signed_boundary_distances" });
  emit(`main_${sp}_boundary_class_summary`, boundaryClassSummaryFigure(m, boundaries),
    { tsv: bndTsv, model: m, figureType: "boundary_class_summary" });

  if (selectedBoundary) {
    const hit = boundaries.find((b) => b.label === selectedBoundary
      || b.exon_boundary_id === selectedBoundary || b.boundary_id === selectedBoundary);
    if (!hit) {
      console.error(`${sp}: no boundary matches "${selectedBoundary}"`);
    } else {
      const id = hit.exon_boundary_id || hit.boundary_id;
      emit(`ondemand_${sp}_selected_boundary_detail`,
        boundaryFigure(m, boundaries, { selectedBoundaryId: id }),
        { model: m, figureType: "selected_boundary_detail" });
      emit(`ondemand_${sp}_selected_signed_distance`,
        signedDistanceFigure(m, boundaries, { selectedBoundaryId: id }),
        { model: m, figureType: "selected_signed_distance" });
    }
  }
}

writeFileSync(join(outDir, "render_summary.json"), JSON.stringify(written, null, 2));
for (const w of written) {
  console.log(`  ${w.stem}: ${w.width}x${w.height}pt, ${w.marks} marks`
    + (w.warnings.length ? ` (${w.warnings.length} layout warnings)` : ""));
}
console.log(`OK — rendered ${written.length} main figure(s) to ${outDir}`);
