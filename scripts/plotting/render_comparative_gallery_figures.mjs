// Headless renderer for Comparative Figure Gallery cards (Parts 3–5).
//
// Imports the same builders the Gallery exports through and feeds them the
// comparative_dataset_index + coordinate model. Boundary figures reuse
// comparativeFigures.js so Explorer exports and Gallery cards cannot diverge.
//
// Only figures that carry a scientific visualisation are emitted; the numbers
// that used to be their own text page (pairwise identity, boundary-alignment
// coverage) are marks inside the MSA overview and rows in the source tables.
//
// Usage:
//   node scripts/plotting/render_comparative_gallery_figures.mjs \
//        <coordinateModel.json> <comparativeDataset.json> <outDir>

import { mkdirSync, readFileSync, writeFileSync, existsSync } from "node:fs";
import { join, resolve } from "node:path";
import {
  comparativeMatrixFigureSpec, pairedSignedDistanceFigureSpec,
  comparativeArchitectureFigureSpec,
  comparativeLongTsv, comparativeMatrixTsv, comparableMappingTsv,
} from "../../webapp/frontend/src/pages/viewers/comparativeFigures.js";
import {
  msaAlignedExonArchitectureFigureSpec, nativeExonArchitectureFigureSpec,
  comparativeDomainArchitectureFigureSpec, domainAnnotationMatrixFigureSpec,
  comparativeExonDomainArchitectureFigureSpec,
  isoformDiversityFigureSpec, comparativeSyntenyFigureSpec,
  comparativeSyntenyMatrixFigureSpec, comparativeSyntenyRows,
  primaryMsaOverviewFigureSpec, boundaryConsistencyPanelFigureSpec,
} from "../../webapp/frontend/src/pages/viewers/comparativeGalleryFigures.js";
import { tsv } from "../../webapp/frontend/src/pages/viewers/mainFigures.js";

const [modelPath, comparativePath, outDir] = process.argv.slice(2);
if (!modelPath || !comparativePath || !outDir) {
  console.error("usage: node scripts/plotting/render_comparative_gallery_figures.mjs "
    + "<model.json> <comparative.json> <outDir>");
  process.exit(2);
}

const index = JSON.parse(readFileSync(modelPath, "utf8"));
const comparative = JSON.parse(readFileSync(comparativePath, "utf8"));
const models = index.models || [];
const gene = index.gene_symbol || comparative.gene_symbol || "gene";
const multi = (index.boundary_dashboard || {}).multi_species || {};
const threshold = (index.boundary_dashboard || {}).near_edge_threshold_aa ?? 5;

if ((comparative.n_species || models.length) < 2) {
  console.error("comparative gallery requires at least two species");
  process.exit(1);
}

mkdirSync(outDir, { recursive: true });
const written = [];

function emit(stem, fig, { tsv: table } = {}) {
  if (!fig) return;
  writeFileSync(join(outDir, `${stem}.svg`), fig.toSvg(), "utf8");
  writeFileSync(join(outDir, `${stem}.pdf`), Buffer.from(fig.toPdf()));
  if (table) writeFileSync(join(outDir, `${stem}.tsv`), table, "utf8");
  written.push({
    stem, width: fig.width, height: fig.height,
    marks: fig.marks.length, warnings: fig.warnings || [], has_tsv: Boolean(table),
  });
}

const nColumns = comparative.msa?.n_columns || 0;
const inventory = comparative.species_inventory || [];
const alignedExons = comparative.msa_aligned_exons || [];
const alignedDomains = comparative.msa_aligned_domains || [];
const boundaryGroups = multi.comparable_boundary_groups
  || comparative.comparable_boundary_groups || [];

// --- Comparative exon structure --------------------------------------------- //
emit("cmp_msa_aligned_exon_architecture",
  msaAlignedExonArchitectureFigureSpec({
    gene, exons: alignedExons, nColumns, boundaryGroups, inventory,
  }),
  { tsv: tsv(alignedExons, ["species_id", "scientific_name", "protein_id",
    "exon_id", "exon_label", "native_start", "native_end",
    "msa_start_column", "msa_end_column", "msa_mapping_status"]) });

emit("cmp_native_exon_architecture",
  nativeExonArchitectureFigureSpec({ gene, models }),
  { tsv: tsv(models.flatMap((m) => (m.exons || []).map((e) => ({
    species_id: m.species_id, scientific_name: m.scientific_name,
    protein_id: m.protein_id, protein_length: m.protein_length,
    exon_label: e.label, native_start: e.start, native_end: e.end,
  }))), ["species_id", "scientific_name", "protein_id", "protein_length",
    "exon_label", "native_start", "native_end"]) });

// --- Comparative sequence analysis ------------------------------------------ //
// Pairwise identity is a metric of this figure, not a card of its own, so the
// identity table travels with it as the source table.
const alnPath = comparative.msa?.alignment_file;
if (alnPath) {
  const abs = existsSync(alnPath) ? alnPath : resolve(process.cwd(), alnPath);
  if (existsSync(abs)) {
    emit("cmp_primary_msa_overview",
      primaryMsaOverviewFigureSpec({
        gene, alignmentText: readFileSync(abs, "utf8"), inventory,
        exons: alignedExons, domains: alignedDomains,
      }),
      { tsv: tsv(comparative.pairwise_identity || [],
        ["species_a", "protein_a", "species_b", "protein_b",
          "n_compared_columns", "n_identical", "percent_identity"]) });
  }
}

// --- Comparative domain architecture ---------------------------------------- //
const domainTsv = tsv(alignedDomains, ["species_id", "scientific_name", "protein_id",
  "domain_instance_id", "interpro_accession", "label", "instance_number",
  "native_start", "native_end", "msa_start_column", "msa_end_column",
  "msa_mapping_status", "order_along_protein"]);

emit("cmp_domain_architecture_native",
  comparativeDomainArchitectureFigureSpec({
    gene, models, domains: alignedDomains, mode: "native", exons: alignedExons,
  }),
  { tsv: domainTsv });

if (nColumns) {
  emit("cmp_domain_architecture_msa",
    comparativeDomainArchitectureFigureSpec({
      gene, models, domains: alignedDomains, mode: "msa", nColumns,
      exons: alignedExons,
    }),
    { tsv: domainTsv });
}

// The exon and domain tracks of one species side by side: the relationship the
// boundary analysis is about, shown per species rather than across two figures.
emit("cmp_exon_domain_architecture_native",
  comparativeExonDomainArchitectureFigureSpec({
    gene, models, domains: alignedDomains, exons: alignedExons, mode: "native",
  }),
  { tsv: domainTsv });

if (nColumns) {
  emit("cmp_exon_domain_architecture_msa",
    comparativeExonDomainArchitectureFigureSpec({
      gene, models, domains: alignedDomains, exons: alignedExons, mode: "msa",
      nColumns,
    }),
    { tsv: domainTsv });
}

if ((comparative.domain_annotation_matrix || []).length) {
  emit("cmp_domain_annotation_matrix",
    domainAnnotationMatrixFigureSpec({
      gene,
      matrix: comparative.domain_annotation_matrix || [],
      groups: comparative.comparable_domain_groups || [],
    }),
    { tsv: tsv(comparative.domain_annotation_matrix || [],
      ["species_id", "scientific_name", "comparable_domain_group_id",
        "interpro_accession", "label", "state", "domain_instance_id",
        "native_start", "native_end", "msa_start_column", "msa_end_column"]) });
}

// --- Comparative isoform diversity ------------------------------------------ //
if ((comparative.isoform_diversity || []).length) {
  emit("cmp_isoform_diversity",
    isoformDiversityFigureSpec({ gene, rows: comparative.isoform_diversity }),
    { tsv: tsv(comparative.isoform_diversity, ["species_id", "scientific_name",
      "primary_protein_id", "n_protein_models", "n_curated_models",
      "n_predicted_models", "primary_protein_length", "protein_length_min",
      "protein_length_max", "n_variable_alignment_blocks",
      "n_exploratory_candidates", "max_difference_from_primary_aa"]) });
}

// --- Comparative genomic context -------------------------------------------- //
const syntenyNeighbourhood = comparative.synteny_neighbourhood || [];
const syntenyTsvColumns = ["species_id", "scientific_name", "target_gene", "slot_x",
  "is_target", "neighbour_symbol", "side", "order", "orientation", "gene_id",
  "placeholder", "classification", "status", "orthology_confidence",
  "distance_to_target"];
// Every dataset species keeps a row: a comparative figure that silently drops a
// species would understate the comparison, so unresolved species are disclosed.
const syntenyRows = comparativeSyntenyRows({
  gene, syntenyNeighbourhood, datasetSpecies: inventory,
});
if (syntenyRows.length > 1) {
  emit("cmp_comparative_synteny",
    comparativeSyntenyFigureSpec({ gene, syntenyNeighbourhood,
      datasetSpecies: inventory }),
    { tsv: tsv(comparative.synteny || [], syntenyTsvColumns) });
  if (syntenyNeighbourhood.length > 1) {
    emit("cmp_synteny_neighbour_conservation",
      comparativeSyntenyMatrixFigureSpec({ gene, syntenyNeighbourhood,
        datasetSpecies: inventory }),
      { tsv: tsv(comparative.synteny || [], syntenyTsvColumns) });
  }
}

// --- Comparative exon–domain boundaries ------------------------------------- //
// Same builders as the accepted Comparative Boundary Explorer exports, so the
// Gallery cards carry the Explorer's values by construction.
if (multi.available && boundaryGroups.length) {
  const stats = multi.distance_statistics || [];
  const matrix = multi.boundary_matrix || [];
  emit("cmp_boundary_matrix",
    comparativeMatrixFigureSpec({
      gene, matrix, groups: boundaryGroups, mode: "signed",
      nearEdgeThreshold: threshold, presetName: "full",
    }),
    { tsv: comparativeMatrixTsv(matrix, boundaryGroups, "signed") });
  emit("cmp_paired_signed_distance",
    pairedSignedDistanceFigureSpec({
      gene, groups: boundaryGroups, stats, nearEdgeThreshold: threshold,
      presetName: "full",
    }),
    { tsv: comparativeLongTsv(boundaryGroups) });
  emit("cmp_boundary_position_consistency",
    boundaryConsistencyPanelFigureSpec({
      gene, stats, groups: boundaryGroups, nearEdgeThreshold: threshold,
    }),
    { tsv: comparableMappingTsv(boundaryGroups, stats) });
  const selected = boundaryGroups[0];
  if (selected) {
    emit("cmp_local_boundary_architecture",
      comparativeArchitectureFigureSpec({
        gene, group: selected, models, presetName: "full",
      }));
  }
}

writeFileSync(join(outDir, "render_summary.json"),
  JSON.stringify(written, null, 2), "utf8");
console.log(`rendered ${written.length} comparative gallery figure(s) → ${outDir}`);
