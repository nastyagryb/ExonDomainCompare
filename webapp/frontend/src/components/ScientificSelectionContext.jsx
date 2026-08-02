/* eslint-disable react-refresh/only-export-components */
// One authoritative selection state for the shared exploratory Gene Explorer.
//
// Every linked tab (Candidate ranking, Transcript & Exon Structure, Protein
// Architecture, Isoform Alignment, Evidence) reads from THIS object and never
// keeps a competing local selection. Relationships between entities
// (transcript ↔ protein, exon → projected aa range, candidate → affected
// proteins/transcripts/exons/coordinates) are derived here once from the
// CanonicalDatasetModel so the whole workspace stays consistent.
import { createContext, useContext, useMemo, useState, useCallback, useEffect } from "react";

const ScientificSelectionContext = createContext(null);

// Session-scoped persistence so a linked selection (candidate / exon / domain /
// boundary) survives a full page refresh, keyed per run + species. Tab changes
// already persist via the provider; this adds refresh-survival (Part 7).
function persistKey(model, speciesId) {
  const run = model?.run_id || model?.analysis?.run_id || "run";
  return `edc.selection.${run}.${speciesId || "sp"}`;
}
function loadPersisted(model, speciesId) {
  try {
    const raw = sessionStorage.getItem(persistKey(model, speciesId));
    return raw ? JSON.parse(raw) : null;
  } catch { return null; }
}

function proteinsOf(species) {
  return species?.proteins || species?.protein_models
    || species?.gene_explorer?.isoforms || species?.isoforms || [];
}

// Rank candidates by score (highest first); assign stable C1…Cn labels.
export function rankCandidates(candidates) {
  return [...(candidates || [])]
    .sort((a, b) => (b.overall_score || 0) - (a.overall_score || 0))
    .map((c, i) => ({ ...c, rank: i + 1, rank_label: `C${i + 1}` }));
}

function transcriptIdForProtein(transcripts, proteins, proteinId) {
  const tx = (transcripts || []).find((t) => t.protein_id === proteinId);
  if (tx) return tx.transcript_id;
  const p = (proteins || []).find((x) => x.protein_id === proteinId);
  return p?.transcript_id || "";
}

// Pick the section for `speciesId` from a species-aware index (top-level
// `species: [...]`), falling back to the reference/top-level fields for
// single-species or legacy indices. Fully generic (keyed only on species_id).
function sliceBySpecies(index, speciesId, key) {
  const arr = index?.species;
  if (Array.isArray(arr) && arr.length) {
    const hit = arr.find((s) => (s.species_id || s.species) === speciesId) || arr[0];
    return hit?.[key] || [];
  }
  return index?.[key] || [];
}

export function ScientificSelectionProvider({ children, species, model }) {
  const speciesId = species?.species || species?.species_id || "";
  const proteins = useMemo(() => proteinsOf(species), [species]);
  const transcripts = useMemo(
    () => sliceBySpecies(model?.transcript_exon_structure, speciesId, "transcripts"),
    [model, speciesId]);
  const rankedCandidates = useMemo(
    () => rankCandidates(sliceBySpecies(model?.candidate_evidence, speciesId, "candidates")),
    [model, speciesId]);

  const primaryProteinId = useMemo(() =>
    species?.selected_primary_protein
    || model?.selected_primary_protein
    || proteins.find((p) => p.is_primary || p.primary_status === "primary")?.protein_id
    || proteins[0]?.protein_id
    || "", [species, model, proteins]);
  const primaryTranscriptId = useMemo(
    () => transcriptIdForProtein(transcripts, proteins, primaryProteinId),
    [transcripts, proteins, primaryProteinId]);

  // Defaults: the top-ranked candidate (C1) is selected, which pins its
  // reference protein + projected aa range across the whole workspace.
  // A persisted session selection (page refresh) overrides the defaults.
  const persisted = useMemo(() => loadPersisted(model, speciesId), [model, speciesId]);
  const topCandidate = rankedCandidates[0] || null;
  const [selectedCandidateId, setSelectedCandidateId] = useState(
    persisted?.selectedCandidateId || topCandidate?.candidate_id || null);
  const [selectedProteinId, setSelectedProteinId] = useState(
    persisted?.selectedProteinId || topCandidate?.reference_protein || primaryProteinId);
  const [selectedTranscriptId, setSelectedTranscriptId] = useState(
    persisted?.selectedTranscriptId || transcriptIdForProtein(transcripts, proteins,
      topCandidate?.reference_protein || primaryProteinId));
  const [selectedExonId, setSelectedExonId] = useState(persisted?.selectedExonId || null);
  const [alignmentStart, setAlignmentStart] = useState(topCandidate?.aa_start ?? null);
  const [alignmentEnd, setAlignmentEnd] = useState(topCandidate?.aa_end ?? null);
  // Extended linked-selection state (Exon Map / Domain Architecture / Boundary).
  const [selectedDomainId, setSelectedDomainId] = useState(persisted?.selectedDomainId || null);
  const [selectedFeatureId, setSelectedFeatureId] = useState(null);
  const [selectedSignatureId, setSelectedSignatureId] = useState(null);
  const [selectedBoundaryId, setSelectedBoundaryId] = useState(persisted?.selectedBoundaryId || null);
  // Comparative selection (multi-species Boundary Explorer). A comparable-boundary
  // group is a cross-species object, so it cannot be expressed by the per-species
  // boundary id alone; keeping it here rather than in the component means the matrix,
  // the paired plot, the detail panel and the inspection cases are driven by one
  // selection instead of four that can drift apart.
  const [selectedComparableGroupId, setSelectedComparableGroupId] = useState(null);
  const [comparativeSpeciesId, setComparativeSpeciesId] = useState(null);
  const [coordinateStart, setCoordinateStart] = useState(null);
  const [coordinateEnd, setCoordinateEnd] = useState(null);
  const [coordinateMode, setCoordinateMode] = useState("native");
  const [visibleTracks, setVisibleTracks] = useState(null);
  // Part 6 linked interactions shared across Gallery, Downloads and Explorers.
  const [selectedFigureScope, setSelectedFigureScope] = useState(
    persisted?.selectedFigureScope || "comparative");
  const [selectedMSAColumn, setSelectedMSAColumn] = useState(
    persisted?.selectedMSAColumn || null);

  const candidateById = useCallback(
    (id) => rankedCandidates.find((c) => c.candidate_id === id) || null,
    [rankedCandidates]);

  // ---- relationship helpers (candidate → affected entities) --------------- //
  const affectedProteinsFor = useCallback((id) => {
    const c = candidateById(id);
    return new Set(c?.protein_isoform_evidence?.affected_proteins || []);
  }, [candidateById]);
  const affectedTranscriptsFor = useCallback((id) => {
    const c = candidateById(id);
    return new Set(c?.protein_isoform_evidence?.affected_transcripts || []);
  }, [candidateById]);
  const supportingExonIdsFor = useCallback((id) => {
    const c = candidateById(id);
    return new Set(c?.exon_evidence?.exon_ids || []);
  }, [candidateById]);

  // ---- selection actions (maintain relationships) ------------------------- //
  const selectCandidate = useCallback((candidateOrId) => {
    const id = typeof candidateOrId === "string"
      ? candidateOrId : candidateOrId?.candidate_id;
    const c = candidateById(id) || (typeof candidateOrId === "object" ? candidateOrId : null);
    setSelectedCandidateId(id || null);
    if (c) {
      const ref = c.reference_protein;
      if (ref) {
        setSelectedProteinId(ref);
        setSelectedTranscriptId(transcriptIdForProtein(transcripts, proteins, ref));
      }
      if (c.aa_start != null && c.aa_end != null) {
        setAlignmentStart(c.aa_start);
        setAlignmentEnd(c.aa_end);
      }
      setSelectedExonId(null);
    }
  }, [candidateById, transcripts, proteins]);

  const selectProtein = useCallback((proteinOrId) => {
    const id = typeof proteinOrId === "string" ? proteinOrId : proteinOrId?.protein_id;
    if (!id) return;
    setSelectedProteinId(id);
    setSelectedTranscriptId(transcriptIdForProtein(transcripts, proteins, id));
    setSelectedExonId(null);
  }, [transcripts, proteins]);

  const selectTranscript = useCallback((txOrId) => {
    const id = typeof txOrId === "string" ? txOrId : txOrId?.transcript_id;
    const tx = transcripts.find((t) => t.transcript_id === id)
      || (typeof txOrId === "object" ? txOrId : null);
    if (!id) return;
    setSelectedTranscriptId(id);
    if (tx?.protein_id) setSelectedProteinId(tx.protein_id);
    setSelectedExonId(null);
  }, [transcripts]);

  const selectExon = useCallback((exon) => {
    if (!exon) return;
    setSelectedExonId(exon.exon_id || null);
    if (exon.transcript_id) setSelectedTranscriptId(exon.transcript_id);
    if (exon.protein_id) setSelectedProteinId(exon.protein_id);
    if (exon.protein_start_aa != null && exon.protein_end_aa != null) {
      setAlignmentStart(exon.protein_start_aa);
      setAlignmentEnd(exon.protein_end_aa);
    }
  }, []);

  const selectDomain = useCallback((domainOrId) => {
    const id = typeof domainOrId === "string" ? domainOrId : domainOrId?.id;
    setSelectedDomainId(id || null);
    setSelectedFeatureId(id || null);
    if (typeof domainOrId === "object" && domainOrId?.start != null && domainOrId?.end != null) {
      setAlignmentStart(domainOrId.start);
      setAlignmentEnd(domainOrId.end);
    }
  }, []);
  // Generic non-domain feature (family / signature / site / disorder / TM).
  const selectFeature = useCallback((featureOrId) => {
    const id = typeof featureOrId === "string" ? featureOrId : featureOrId?.id;
    setSelectedFeatureId(id || null);
    if (typeof featureOrId === "object" && featureOrId?.feature_type === "member_signature") {
      setSelectedSignatureId(id || null);
    }
    if (typeof featureOrId === "object" && featureOrId?.start != null && featureOrId?.end != null) {
      setAlignmentStart(featureOrId.start);
      setAlignmentEnd(featureOrId.end);
    }
  }, []);
  const selectSignature = useCallback((sigOrId) => {
    const id = typeof sigOrId === "string" ? sigOrId : sigOrId?.id;
    setSelectedSignatureId(id || null);
    setSelectedFeatureId(id || null);
  }, []);
  const selectBoundary = useCallback((boundaryOrId) => {
    const id = typeof boundaryOrId === "string"
      ? boundaryOrId : (boundaryOrId?.id || boundaryOrId?.boundary_id);
    setSelectedBoundaryId(id || null);
    // Cross-tab links: selecting a boundary also pins its nearest representative
    // domain (Domain Architecture highlight) and centres the shared aa range
    // (Exon Map). No fabrication — only echoes fields already on the boundary.
    if (typeof boundaryOrId === "object" && boundaryOrId) {
      if (boundaryOrId.nearest_domain_id) setSelectedDomainId(boundaryOrId.nearest_domain_id);
      const pos = boundaryOrId.protein_position ?? boundaryOrId.start;
      if (pos != null) {
        const pad = 40;
        setAlignmentStart(Math.max(1, pos - pad));
        setAlignmentEnd(pos + pad);
      }
    }
  }, []);
  /**
   * Select a comparable-boundary group, optionally through one species' observation.
   *
   * Passing the observation lets a single matrix-cell click set the comparative group,
   * the species whose row was clicked and that species' own boundary at once, which is
   * what keeps the detail panel and the comparative architecture in step with the cell
   * the user actually pointed at.
   */
  const selectComparableGroup = useCallback((groupOrId, observation = null) => {
    const id = typeof groupOrId === "string"
      ? groupOrId : groupOrId?.comparable_boundary_group_id;
    setSelectedComparableGroupId(id || null);
    if (observation) {
      if (observation.species_id) setComparativeSpeciesId(observation.species_id);
      if (observation.boundary_id) setSelectedBoundaryId(observation.boundary_id);
      if (observation.msa_column != null) setSelectedMSAColumn(observation.msa_column);
      if (observation.native_position != null) {
        const pad = 40;
        setAlignmentStart(Math.max(1, observation.native_position - pad));
        setAlignmentEnd(observation.native_position + pad);
      }
    }
  }, []);

  /** Open a species-specific scientific page while preserving linked selection. */
  const openSpeciesView = useCallback((targetSpeciesId, viewHint = null) => {
    if (targetSpeciesId) {
      setComparativeSpeciesId(targetSpeciesId);
      setSelectedFigureScope(targetSpeciesId);
    }
    if (viewHint === "msa") setCoordinateMode("msa");
    if (viewHint === "native") setCoordinateMode("native");
  }, []);

  const setCoordinateRange = useCallback((start, end) => {
    setCoordinateStart(start);
    setCoordinateEnd(end);
  }, []);

  const selectAlignmentRegion = useCallback((start, end) => {
    setAlignmentStart(start);
    setAlignmentEnd(end);
    // Select an overlapping candidate if one exists (alignment → candidate).
    const overlap = rankedCandidates.find(
      (c) => c.aa_start <= end && c.aa_end >= start);
    if (overlap) setSelectedCandidateId(overlap.candidate_id);
  }, [rankedCandidates]);

  // Resolved objects (kept as convenience for consumers).
  const selectedCandidate = useMemo(
    () => candidateById(selectedCandidateId), [candidateById, selectedCandidateId]);
  const selectedExon = useMemo(() => {
    if (!selectedExonId) return null;
    for (const tx of transcripts) {
      const ex = (tx.exons || []).find((e) => e.exon_id === selectedExonId);
      if (ex) return ex;
    }
    return null;
  }, [selectedExonId, transcripts]);

  // Persist the linked selection for refresh-survival (session-scoped).
  useEffect(() => {
    try {
      sessionStorage.setItem(persistKey(model, speciesId), JSON.stringify({
        selectedCandidateId, selectedProteinId, selectedTranscriptId,
        selectedExonId, selectedDomainId, selectedBoundaryId,
        selectedComparableGroupId, selectedFigureScope, selectedMSAColumn,
        coordinateMode, comparativeSpeciesId,
      }));
    } catch { /* storage unavailable — selection still persists across tabs */ }
  }, [model, speciesId, selectedCandidateId, selectedProteinId, selectedTranscriptId,
    selectedExonId, selectedDomainId, selectedBoundaryId, selectedComparableGroupId,
    selectedFigureScope, selectedMSAColumn, coordinateMode, comparativeSpeciesId]);

  const value = useMemo(() => ({
    // authoritative ids
    selectedSpeciesId: speciesId,
    selectedTranscriptId,
    selectedProteinId,
    selectedExonId,
    selectedCandidateId,
    selectedAlignmentStart: alignmentStart,
    selectedAlignmentEnd: alignmentEnd,
    // back-compat aliases (existing consumers)
    selectedSpecies: speciesId,
    selectedProtein: selectedProteinId,
    selectedTranscript: selectedTranscriptId,
    selectedExon,
    selectedCandidate,
    selectedAlignmentRegion:
      alignmentStart != null && alignmentEnd != null ? [alignmentStart, alignmentEnd] : null,
    // extended Exon Map / Domain / Boundary linked state (persists across tabs)
    selectedDomainId, selectedFeatureId, selectedSignatureId, selectedBoundaryId,
    coordinateStart, coordinateEnd, coordinateMode, visibleTracks,
    // Part 17 canonical alias names (global Exon–Domain-Boundaries page)
    selectedCoordinateStart: coordinateStart, selectedCoordinateEnd: coordinateEnd,
    selectedVisibleTracks: visibleTracks,
    // data + derived
    model, species, proteins, transcripts, rankedCandidates,
    primaryProteinId, primaryTranscriptId,
    candidateById, affectedProteinsFor, affectedTranscriptsFor, supportingExonIdsFor,
    // actions
    selectCandidate, selectProtein, selectTranscript, selectExon,
    selectDomain, selectFeature, selectSignature, selectBoundary,
    selectedComparableGroupId,
    selectedComparableBoundaryGroupId: selectedComparableGroupId,
    selectedComparativeSpeciesId: comparativeSpeciesId,
    selectComparableGroup, selectComparativeSpecies: setComparativeSpeciesId,
    selectedFigureScope, setSelectedFigureScope,
    selectedMSAColumn, setSelectedMSAColumn,
    selectedCoordinateMode: coordinateMode,
    openSpeciesView,
    setCoordinateRange, setCoordinateMode, setVisibleTracks,
    selectAlignmentRegion, setSelectedAlignmentRegion: (r) => {
      if (Array.isArray(r)) selectAlignmentRegion(r[0], r[1]);
    },
  }), [speciesId, selectedTranscriptId, selectedProteinId, selectedExonId,
    selectedCandidateId, alignmentStart, alignmentEnd, selectedExon, selectedCandidate,
    selectedDomainId, selectedFeatureId, selectedSignatureId, selectedBoundaryId,
    coordinateStart, coordinateEnd, coordinateMode,
    visibleTracks, model, species, proteins, transcripts, rankedCandidates, primaryProteinId,
    primaryTranscriptId, candidateById, affectedProteinsFor, affectedTranscriptsFor,
    supportingExonIdsFor, selectCandidate, selectProtein, selectTranscript, selectExon,
    selectDomain, selectFeature, selectSignature, selectBoundary, setCoordinateRange,
    setCoordinateMode, setVisibleTracks, selectAlignmentRegion,
    selectedComparableGroupId, comparativeSpeciesId, selectComparableGroup,
    selectedFigureScope, selectedMSAColumn, openSpeciesView]);

  return (
    <ScientificSelectionContext.Provider value={value}>
      {children}
    </ScientificSelectionContext.Provider>
  );
}

export function useScientificSelection() {
  return useContext(ScientificSelectionContext);
}
