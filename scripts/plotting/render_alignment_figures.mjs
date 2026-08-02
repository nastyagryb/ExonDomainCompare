// Render the isoform-alignment figures from a real run, outside any browser.
//
// This is the production renderer for the Gallery's isoform-analysis figures and
// the one the test suite validates, so a Gallery figure and the figure the Gene
// Explorer exports come from a single implementation. It imports the same builders
// (`alignmentFigure.js`) the interactive viewer exports through, and running
// headlessly also proves the figures carry no browser-only state.
//
// Usage:
//   node scripts/plotting/render_alignment_figures.mjs <alignment_index.json> <out_dir>
//                                                     [--gallery-stems]
//
// With --gallery-stems the files are named `main_<species_id>_<kind>`, which is what
// the Gallery cards point at; without it the plain figure names are used.

import { mkdirSync, readFileSync, writeFileSync, existsSync } from "node:fs";
import { join, dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const VIEWERS = resolve(here, "../../webapp/frontend/src/pages/viewers");
const mod = await import(join(VIEWERS, "alignmentFigure.js"));
const { renderPdfPages } = await import(join(VIEWERS, "figureSpec.js"));

const argv = process.argv.slice(2);
const galleryStems = argv.includes("--gallery-stems");
const wanted = (argv.find((a) => a.startsWith("--species=")) || "").slice(10);
const [indexPath, outDir] = argv.filter((a) => !a.startsWith("--"));
if (!indexPath || !outDir) {
  console.error("usage: render_alignment_figures.mjs <index.json> <out_dir> "
    + "[--gallery-stems] [--species=<species_id>]");
  process.exit(2);
}
mkdirSync(outDir, { recursive: true });

const index = JSON.parse(readFileSync(indexPath, "utf8"));
// A multi-species run holds one within-species alignment per species. Which one is
// rendered is the caller's decision; without --species the first available one is
// used, which is the single-species behaviour.
const available = (index.species || []).filter((s) => s.status === "available");
const entry = wanted
  ? available.find((s) => (s.species_id || s.species) === wanted)
  : available[0];
if (!entry) {
  console.error(wanted
    ? `no available alignment for species ${wanted}`
    : "no available species alignment in index");
  process.exit(3);
}

// Row metadata (curation status, protein length, coding exons) comes from the
// validated coordinate model next to the alignment index, exactly as the
// application resolves it — never guessed from the accession prefix.
const coordPath = resolve(dirname(indexPath), "generic", "protein_coordinate_model.json");
const coord = existsSync(coordPath) ? JSON.parse(readFileSync(coordPath, "utf8")) : null;
const speciesId = entry.species_id || entry.species;
const coordModel = (coord?.models || [])
  .find((m) => (m.species_id || m.species) === speciesId) || (coord?.models || [])[0] || null;
const byProtein = new Map();
for (const t of coordModel?.transcript_models || []) {
  if (t.protein_id) byProtein.set(t.protein_id, t);
}

const rows = mod.orderRows((entry.sequences || []).map((s) => {
  const aligned = s.aligned_sequence || s.seq || "";
  const tm = byProtein.get(s.protein_id);
  return {
    protein_id: s.protein_id,
    transcript_id: s.transcript_id || tm?.transcript_id || null,
    is_primary: Boolean(s.is_primary),
    seq: aligned,
    protein_length: tm?.protein_length ?? aligned.replace(/-/g, "").length,
    curation_status: tm?.curation_status
      || (/^NP_/.test(s.protein_id || "") ? "curated" : "predicted"),
    display_species_name: entry.display_species_name,
  };
}));

const nCols = entry.alignment_length || rows[0].seq.length;
const primary = rows.find((r) => r.is_primary) || rows[0];
const gene = index.gene_symbol || "GENE";
const species = entry.display_species_name || speciesId;
const tool = index.tool || "MAFFT";
const exons = (coordModel?.exons || []).map((e) => ({
  label: e.label, start: e.start, end: e.end,
}));

// Candidate identity comes from the validated coordinate model, which numbers the
// candidates by position and is what every other figure and table labels them by.
// The alignment index lists them in scan order and carries no rank, so numbering
// them here by list position would give the same region two different names in two
// figures of the same Gallery.
const canonicalCandidates = new Map(
  (coordModel?.candidate_regions || []).map((c) => [`${c.start}-${c.end}`, c]));

// Candidate intervals, mapped from primary-protein residues onto alignment
// columns with the same helper the interactive viewer uses.
// A candidate the model does not list is a difference between two *alternative*
// isoforms: its coordinates are in that pair's frame, not the primary's, so it has
// no position on this alignment's residue axis. Numbering it positionally invented a
// canonical name for a region no other figure or table knows, so it is skipped and
// counted instead.
const offAxisCandidates = [];
const candidates = (entry.candidates || [])
  .slice()
  .sort((a, b) => (a.aa_start - b.aa_start) || (a.aa_end - b.aa_end))
  .map((c) => {
    const canon = canonicalCandidates.get(`${c.aa_start}-${c.aa_end}`);
    if (!canon) {
      offAxisCandidates.push(`${c.aa_start}-${c.aa_end}`);
      return null;
    }
    return {
      label: canon.id,
      candidate_id: c.candidate_id,
      aa_start: c.aa_start,
      aa_end: c.aa_end,
      col_start: mod.aaToColumn(primary.seq, c.aa_start),
      col_end: mod.aaToColumn(primary.seq, c.aa_end),
    };
  })
  .filter((c) => c && c.col_start != null && c.col_end != null);
if (offAxisCandidates.length) {
  console.error(`NOTE: ${offAxisCandidates.length} candidate interval(s) have no `
    + `primary-protein coordinate and are not drawn: ${offAxisCandidates.join(", ")}`);
}
const candidate = candidates[0] || null;

/** Gallery files are named per species; the bare names are used for inspection. */
const stem = (name) => (galleryStems ? `main_${speciesId}_${name}` : name);

const write = (name, fig) => {
  writeFileSync(join(outDir, `${stem(name)}.svg`), fig.toSvg(), "utf8");
  writeFileSync(join(outDir, `${stem(name)}.pdf`), Buffer.from(fig.toPdf()));
  return { name, stem: stem(name), width: fig.width, height: fig.height,
    marks: fig.marks.length, warnings: fig.warnings };
};

const figureStats = [];

// --- Figure 1: the complete alignment at column resolution -------------------
const overview = mod.alignmentOverviewFigureSpec({
  rows, nCols, gene, species, primaryId: primary.protein_id,
  transcriptId: primary.transcript_id, candidates, tool,
});
figureStats.push(write("full_isoform_alignment", overview));

// --- Figure 2: the same alignment at residue resolution, wrapped -------------
const wrapped = mod.wrappedAlignmentFigureSpecs({
  rows, nCols, gene, species, primaryId: primary.protein_id,
  transcriptId: primary.transcript_id, candidates, tool,
});
const layout = mod.wrappedAlignmentLayout({
  nRows: rows.length, nCols, hasCandidates: candidates.length > 0,
});
wrapped.forEach((fig, i) => {
  const name = `wrapped_alignment_p${i + 1}`;
  writeFileSync(join(outDir, `${stem(name)}.svg`), fig.toSvg(), "utf8");
  figureStats.push({ name, stem: stem(name), width: fig.width,
    height: fig.height, marks: fig.marks.length, warnings: fig.warnings });
});
// The wrapped alignment is one document with several pages, not several files.
writeFileSync(join(outDir, `${stem("wrapped_alignment")}.pdf`),
  Buffer.from(renderPdfPages(wrapped)));
// The Gallery preview needs a single image, so page one doubles as the card's SVG.
if (galleryStems && wrapped.length) {
  writeFileSync(join(outDir, `${stem("wrapped_alignment")}.svg`),
    wrapped[0].toSvg(), "utf8");
}

// --- Figure 3: the candidate interval at residue resolution ------------------
const summary = {
  gene, species, primary: primary.protein_id, transcript: primary.transcript_id,
  n_rows: rows.length, n_columns: nCols, tool,
  candidate: candidate ? candidate.label : null,
  n_candidates: candidates.length,
  protein_ids: rows.map((r) => r.protein_id),
  identities: Object.fromEntries(rows.map((r) => [r.protein_id,
    r.is_primary ? 100 : mod.identityPct(r.seq, primary.seq)])),
  curation: Object.fromEntries(rows.map((r) => [r.protein_id, r.curation_status])),
  wrapped: { ...layout, n_pages_written: wrapped.length },
  n_exons: exons.length,
  figures: ["full_isoform_alignment", "wrapped_alignment", ...(candidate ? ["candidate_alignment_detail"] : [])],
};

if (candidate) {
  const affected = mod.affectedProteins(rows, primary.seq, candidate.col_start, candidate.col_end);
  const detail = mod.candidateAlignmentFigureSpec({
    rows, nCols, gene, species, primaryId: primary.protein_id,
    candidate, affected, exons,
  });
  figureStats.push(write("candidate_alignment_detail", detail));
  summary.candidate_columns = [candidate.col_start + 1, candidate.col_end + 1];
  summary.candidate_aa = [candidate.aa_start, candidate.aa_end];
  summary.affected = [...affected].sort();
  summary.unaffected = rows.filter((r) => !r.is_primary && !affected.has(r.protein_id))
    .map((r) => r.protein_id).sort();
}

// --- source data that ships with the figures ---------------------------------
// Every alignment figure is a view of one alignment, so each ships the same table
// and FASTA: a reader who downloads one figure gets the numbers behind it.
const alnTsv = mod.tsv(mod.alignmentSummaryRows({ rows, nCols, candidates }),
  mod.ALIGNMENT_SUMMARY_COLUMNS);
writeFileSync(join(outDir, `${stem("alignment")}.fasta`), mod.alignmentFasta(rows), "utf8");
if (galleryStems) {
  // One table per card, not per rendered page.
  for (const name of ["full_isoform_alignment", "wrapped_alignment",
    ...(candidate ? ["candidate_alignment_detail"] : [])]) {
    writeFileSync(join(outDir, `${stem(name)}.tsv`), alnTsv, "utf8");
  }
} else {
  writeFileSync(join(outDir, "alignment_summary.tsv"), alnTsv, "utf8");
}

summary.render = figureStats;
summary.species_id = speciesId;
summary.gallery_stems = galleryStems;
writeFileSync(join(outDir, "alignment_render_summary.json"),
  JSON.stringify(summary, null, 2));
// Kept under the historical name too, so existing consumers do not break.
writeFileSync(join(outDir, "summary.json"), JSON.stringify(summary, null, 2));
for (const f of figureStats) {
  console.log(`${f.name}: ${f.width}x${f.height}pt, ${f.marks} marks`
    + (f.warnings.length ? ` (${f.warnings.length} layout warnings)` : ""));
}
console.log(JSON.stringify({ gene: summary.gene, species: summary.species,
  n_rows: summary.n_rows, n_columns: summary.n_columns, primary: summary.primary,
  candidate: summary.candidate, candidate_columns: summary.candidate_columns,
  wrapped: summary.wrapped }));
