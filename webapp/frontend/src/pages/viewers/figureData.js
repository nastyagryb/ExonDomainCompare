// Adapter from the validated protein-coordinate model to the publication figure
// builders.
//
// The interactive viewers and the exported figures read the same model records
// through this one adapter, so a figure downloaded from the Gene Explorer and the
// corresponding Figure Gallery figure cannot disagree about coordinates, feature
// order, labels or selection.

import {
  boundaryClassSummarySpec, boundaryFigureSpec, domainArchitectureFigureSpec,
  exonMapFigureSpec, signedDistanceFigureSpec, tsv,
} from "./mainFigures.js";

/** Common identity fields every single-species figure needs. */
function identity(speciesModel) {
  return {
    gene: speciesModel?.gene_symbol || "",
    species: speciesModel?.scientific_name || speciesModel?.species_id || "",
    proteinId: speciesModel?.protein_id || "",
    transcriptId: speciesModel?.transcript_id || "",
    proteinLength: speciesModel?.protein_length || 1,
    // Part of the identity wherever the analysed gene has more than one primary
    // protein per species. Empty for the usual gene, so nothing changes there.
    isoform: speciesModel?.final_isoform_label || speciesModel?.isoform || "",
  };
}

/**
 * Established events of this protein — not exploratory candidates.
 *
 * The model records a validated event where one exists (FGFR2's IIIb/IIIc cassette
 * exon), with the claim status the freeze assigned. An exploratory candidate region
 * is a different thing and stays in its own track with its own style.
 */
const eventsOf = (m) => (m?.cassette_regions || m?.validated_events || []).map((e) => ({
  ...e,
  event_label: m?.validated_event?.event || e.label || "",
  validated: m?.validated_event?.validated ?? true,
}));

/** Coordinate-model exons carry `start` / `end`; the builders expect the same. */
const exonsOf = (m) => (m?.exons || []).map((e) => ({
  ...e,
  feature_type: "coding_exon",
  label: e.label || `E${e.tooltip?.exon_number ?? ""}`,
  exon_number: e.tooltip?.exon_number,
  genomic_start: e.tooltip?.genomic_start,
  genomic_end: e.tooltip?.genomic_end,
  transcript_id: e.tooltip?.transcript_id,
}));

// The display-cluster identity, rank and lane come from the coordinate model, so the
// exported figure packs candidates exactly like the interactive track does.
const candidatesOf = (m) => (m?.candidate_regions || []).map((c) => ({
  ...c,
  candidate_id: c.id,
  rank_label: c.id,
  aa_start: c.start,
  aa_end: c.end,
  length: c.end - c.start + 1,
  status: c.status,
}));

export function exonMapFigure(speciesModel, { selectedCandidateId = null } = {}) {
  return exonMapFigureSpec({
    ...identity(speciesModel),
    exons: exonsOf(speciesModel),
    candidates: candidatesOf(speciesModel),
    events: eventsOf(speciesModel),
    selectedCandidateId,
    presetName: "full",
  });
}

export function exonMapTsv(speciesModel) {
  return tsv(exonsOf(speciesModel).map((e) => ({
    exon_label: e.label,
    exon_number: e.exon_number,
    protein_start_aa: e.start,
    protein_end_aa: e.end,
    length_aa: e.end - e.start + 1,
    genomic_start: e.genomic_start,
    genomic_end: e.genomic_end,
    transcript_id: e.transcript_id,
  })), ["exon_label", "exon_number", "protein_start_aa", "protein_end_aa",
    "length_aa", "genomic_start", "genomic_end", "transcript_id"]);
}

export function domainArchitectureFigure(speciesModel,
  { selectedDomainInstanceId = null, showAllCandidates = false } = {}) {
  return domainArchitectureFigureSpec({
    ...identity(speciesModel),
    domains: speciesModel?.representative_domains || [],
    families: speciesModel?.families_superfamilies || [],
    tm: speciesModel?.tm_regions || [],
    exons: exonsOf(speciesModel),
    candidates: candidatesOf(speciesModel),
    events: eventsOf(speciesModel),
    selectedDomainInstanceId,
    // The export shows the candidate lanes the reader has open on screen.
    showAllCandidates,
    presetName: "full",
  });
}

export function domainArchitectureTsv(speciesModel) {
  const rows = [];
  for (const d of speciesModel?.representative_domains || []) {
    rows.push({
      track: "representative_domain",
      domain_instance_id: d.domain_instance_id || "",
      label: d.full_label || d.label,
      interpro_accession: d.interpro_accession || "",
      instance_number: d.instance_number ?? "",
      start_aa: d.start, end_aa: d.end, source: d.source || "",
    });
  }
  for (const f of speciesModel?.families_superfamilies || []) {
    rows.push({
      track: "family_superfamily", domain_instance_id: "",
      label: f.label, interpro_accession: f.interpro_accession || "",
      instance_number: "", start_aa: f.start, end_aa: f.end, source: f.source || "",
    });
  }
  for (const t of speciesModel?.tm_regions || []) {
    rows.push({
      track: "tm_region", domain_instance_id: "", label: t.label,
      interpro_accession: "", instance_number: "",
      start_aa: t.start, end_aa: t.end, source: t.source || "",
    });
  }
  return tsv(rows, ["track", "domain_instance_id", "label", "interpro_accession",
    "instance_number", "start_aa", "end_aa", "source"]);
}

/**
 * Boundary figures always take the *filtered* boundary set, so the figure shows
 * exactly the rows the table and the summary counts show.
 */
export function boundaryFigure(speciesModel, boundaries,
  { selectedBoundaryId = null, showCandidates = true,
    selectedCandidateId = null, showAllCandidates = false } = {}) {
  return boundaryFigureSpec({
    ...identity(speciesModel),
    domains: speciesModel?.representative_domains || [],
    exons: exonsOf(speciesModel),
    boundaries,
    candidates: candidatesOf(speciesModel),
    events: eventsOf(speciesModel),
    selectedBoundaryId,
    nearEdgeThreshold: speciesModel?.near_edge_threshold_aa ?? 5,
    showCandidates,
    // The export shows the candidate lanes the reader has open on screen.
    selectedCandidateId,
    showAllCandidates,
    presetName: "full",
  });
}

export function signedDistanceFigure(speciesModel, boundaries,
  { selectedBoundaryId = null, groupByDomain = false } = {}) {
  return signedDistanceFigureSpec({
    ...identity(speciesModel),
    boundaries,
    selectedBoundaryId,
    groupByDomain,
    nearEdgeThreshold: speciesModel?.near_edge_threshold_aa ?? 5,
    presetName: "double",
  });
}

export function boundaryClassSummaryFigure(speciesModel, boundaries) {
  return boundaryClassSummarySpec({
    ...identity(speciesModel),
    boundaries,
    nearEdgeThreshold: speciesModel?.near_edge_threshold_aa ?? 5,
    presetName: "compact",
  });
}

export function boundaryTsv(boundaries) {
  return tsv((boundaries || []).map((b) => ({
    transition: b.label || b.transition_label || "",
    boundary_position_aa: b.boundary_position_aa ?? b.protein_position,
    left_exon: b.left_exon_label || "",
    right_exon: b.right_exon_label || "",
    nearest_domain_instance_id: b.nearest_domain_instance_id || "",
    nearest_domain_accession: b.nearest_domain_accession || "",
    nearest_domain_label: b.nearest_domain_full_label || b.nearest_domain_short_label || "",
    nearest_domain_start: b.nearest_domain_start ?? "",
    nearest_domain_end: b.nearest_domain_end ?? "",
    nearest_edge: b.nearest_edge || b.domain_edge_type || "",
    nearest_edge_position: b.nearest_edge_position ?? "",
    signed_distance_aa: b.signed_distance_aa ?? b.signed_distance ?? "",
    absolute_distance_aa: b.absolute_distance_aa ?? "",
    boundary_class: b.classification || b.category || "",
  })), ["transition", "boundary_position_aa", "left_exon", "right_exon",
    "nearest_domain_instance_id", "nearest_domain_accession", "nearest_domain_label",
    "nearest_domain_start", "nearest_domain_end", "nearest_edge",
    "nearest_edge_position", "signed_distance_aa", "absolute_distance_aa",
    "boundary_class"]);
}
