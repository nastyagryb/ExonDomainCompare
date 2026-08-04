import { useEffect, useMemo, useState } from "react";
import { fileUrl } from "../api";
import { Badge, Field, IsoBadge, Menu, Spinner, Empty } from "../ui";
import { readinessLabel, mainDisplayStatus } from "../uiStatus";
import { downloadBlob } from "./viewers/plotExport";
import { CHROME, featureProps, featureStyle, textProps } from "./viewers/semanticStyles";
import {
  downloadFigurePdf, downloadFigurePdfPages, downloadFigurePng, downloadFigureSvg,
  downloadFigureTsv,
} from "./viewers/figureExport";
import {
  aaToColumn, alignmentFasta, ALIGNMENT_SUMMARY_COLUMNS, alignmentOverviewFigureSpec,
  alignmentSummaryRows, candidateAlignmentFigureSpec, orderRows, wrappedAlignmentFigureSpecs,
} from "./viewers/alignmentFigure";
import {
  GeneExplorerShell, SpeciesPanel, SpeciesCard, WorkspaceHeader, ExplorerTabs, SummaryCard,
  EvidenceSummary, ProteinIsoformTable, PendingAnalysisCard,
} from "../components/shared";
import { ScientificSelectionProvider, useScientificSelection } from "../components/ScientificSelectionContext";
import {
  DATASET_STATUS_META, datasetStatusFromModel, datasetStatusLabel, missingPostClusterItems,
} from "../datasetStatus";
import CassetteExplorer from "./viewers/CassetteExplorer";
import CoordinateTrack from "./viewers/CoordinateTrack";
import ExonMap from "./viewers/ExonMap";
import ProteinArchitecture from "./viewers/ProteinArchitecture";
import BoundaryExplorer from "./viewers/BoundaryExplorer";
import MsaExplorer from "./viewers/MsaExplorer";
import SyntenyViewer from "./viewers/SyntenyViewer";
import DomainArchitecture from "./viewers/DomainArchitecture";
import DataDownloads from "./DataDownloads";
import BoundaryDetailTrack from "./viewers/BoundaryDetailTrack";

const TAXON_FILTERS = [
  ["all", "All"],
  ["Primates", "Primates"],
  ["Other mammals", "Mammals"],
  ["Birds", "Birds"],
  ["Reptiles", "Reptiles"],
  ["Amphibians", "Amphibians"],
  ["Teleost fish", "Fish"],
];

export default function GeneExplorer({ model, target, labels, setPage }) {
  const species = useMemo(() => (model?.species || []).map((sp) => ({
    ...sp,
    species: sp.species || sp.species_id,
    isoforms: sp.isoforms || sp.gene_explorer?.isoforms || [],
    proteins: sp.proteins || sp.protein_architecture?.proteins
      || sp.gene_explorer?.isoforms || [],
    overall_status: sp.overall_status || "accepted",
    taxon_group: sp.taxon_group || "Analysed species",
  })), [model]);
  const downloads = model?.downloads || [];
  const eventType = model?.event_layer?.type || "none";
  const [selected, setSelected] = useState(() => species[0] || null);
  const [search, setSearch] = useState("");
  const [readiness, setReadiness] = useState("all"); // all|primary|review
  const [taxon, setTaxon] = useState("all");
  const [tab, setTab] = useState("summary");
  const eventIndices = model?.validated_event_indices || {};
  const sharedEventIndices = model?.shared_event_indices || {};
  const isExploratory = (model?.event_layer?.type || "none") !== "validated";
  const idx = isExploratory ? {
    coordinateModel: model?.protein_coordinate_model,
    coordinates: sharedEventIndices.coordinates || model?.shared_indices?.coordinate_track_index,
    msa: sharedEventIndices.msa || model?.shared_indices?.msa_index,
    synteny: sharedEventIndices.synteny || model?.shared_indices?.synteny_locus_index,
    domainArch: sharedEventIndices.domainArch || model?.domain_architecture,
    boundaryMatrix: sharedEventIndices.boundaryMatrix || model?.boundary,
  } : {
    cassette: eventIndices.cassette || eventIndices.cassette_residue_index,
    coordinates: eventIndices.coordinates || eventIndices.coordinate_track_index,
    msa: eventIndices.msa || eventIndices.msa_index,
    synteny: eventIndices.synteny || model?.synteny,
    domainArch: eventIndices.domainArch || model?.legacy_fgfr2_indices?.species_domain_architecture,
    boundaryMatrix: eventIndices.boundaryMatrix || eventIndices.boundary_consistency_matrix,
  };

  // deep-link: select the requested species / tab when navigated from elsewhere
  useEffect(() => {
    if (!target?.species || !species?.length) return;
    const match = species.find((s) => s.species === target.species);
    if (match) {
      // Deep-link navigation intentionally synchronises local selection.
      // eslint-disable-next-line react-hooks/set-state-in-effect
      setSelected(match);
      if (target.tab) setTab(target.tab);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [target?._t, species]);

  const filtered = useMemo(() => {
    if (!species) return [];
    return species.filter((s) => {
      if (taxon !== "all" && s.taxon_group !== taxon) return false;
      if (readiness === "primary" && s.overall_readiness === "review") return false;
      if (readiness === "review" && s.overall_readiness !== "review") return false;
      if (search && !`${s.display_species_name} ${s.species}`.toLowerCase().includes(search.toLowerCase())) return false;
      return true;
    });
  }, [species, taxon, readiness, search]);

  if (!model) return <section className="page"><Spinner /></section>;

  const speciesPanel = (
    <SpeciesPanel
      filters={<>
        <input className="search" placeholder="Search species…" value={search} onChange={(e) => setSearch(e.target.value)} />
        <div className="chip-row">
          {[["all", "All"], ["primary", "Primary"], ["review", "Review"]].map(([id, l]) => (
            <button key={id} className={readiness === id ? "chip sel" : "chip"} onClick={() => setReadiness(id)}>{l}</button>
          ))}
        </div>
        <div className="chip-row">
          {TAXON_FILTERS.map(([id, l]) => (
            <button key={id} className={taxon === id ? "chip sel" : "chip"} onClick={() => setTaxon(id)}>{l}</button>
          ))}
        </div>
      </>}
    >
      {filtered.map((s) => (
        <SpeciesCard
          key={s.species}
          selected={selected?.species === s.species}
          onClick={() => { setSelected(s); setTab("summary"); }}
          title={s.display_species_name}
          badge={<Badge cls={s.overall_readiness || s.overall_status || "accepted"}>
            {readinessLabel(s.overall_readiness || s.overall_status || "accepted")}
          </Badge>}
          sub={s.taxon_group || s.scientific_name || "genome annotation"}
          extra={eventType === "validated"
            ? <span className="iso-mini">
                {(s.isoforms || []).slice(0, 4).map((i) => (
                  <span key={i.isoform || i.protein_id}
                    className={`iso iso-${String(i.isoform || "neutral").toLowerCase()} ${i.readiness_class === "review" ? "warn" : ""}`}>
                    {i.isoform || (i.role === "primary" ? "primary" : "alt")}
                  </span>
                ))}
              </span>
            : <span className="iso-mini">
                {s.common_name && <span className="iso iso-neutral">{s.common_name}</span>}
                {s.clade && <span className="iso iso-neutral">{s.taxonomic_group || s.clade}</span>}
              </span>}
        />
      ))}
      {filtered.length === 0 && <p className="muted pad">No species match.</p>}
    </SpeciesPanel>
  );

  return (
    <>
      <GeneExplorerShell sidebar={speciesPanel}>
        {selected && (
          <div className="workspace">
            <ScientificSelectionProvider key={selected.species} species={selected} model={model}>
              <Workspace
                sp={selected}
                tab={tab}
                setTab={setTab}
                downloads={downloads}
                idx={idx}
                labels={labels}
                model={model}
                eventType={eventType}
                setPage={setPage}
              />
            </ScientificSelectionProvider>
          </div>
        )}
      </GeneExplorerShell>
    </>
  );
}

const PHASE2_TABS = [
  ["exon", "Exon map"],
  ["cassette", "Cassette"],
  ["architecture", "Domain architecture"],
  ["boundary", "Boundary"],
  ["msa", "MSA"],
  ["synteny", "Synteny"],
];

function Workspace({ sp, tab, setTab, downloads, idx, labels, model, eventType, setPage }) {
  const isValidated = eventType === "validated";
  const header = (sp.isoforms || []).reduce((acc, i) => { acc[i.isoform] = i; return acc; }, {});
  // Config-driven tab labels with FGFR2 wording as fallback (event_region -> "Cassette").
  const phase2Tabs = PHASE2_TABS.map(([id, l]) =>
    id === "cassette" ? [id, labels?.eventRegion || l] : [id, l]);
  const validatedTabs = [
    ["summary", "Summary"], ["isoforms", "Isoforms"], ["evidence", "Evidence"],
    ["files", "Data & Downloads"],
    ...phase2Tabs.map(([id, l], i) => (i === 0 ? { id, label: l, separatorBefore: true } : { id, label: l })),
  ];
  // Non-FGFR2 genes reuse the EXACT FGFR2 base module order / ids / components
  // (Summary · Isoforms · Evidence · Files | Exon map · Domain architecture ·
  // Boundary · Alignment · Synteny). The single Alignment tab resolves its mode
  // and visible label from the dataset (isoform_alignment vs cross_species_msa).
  // The only biological difference is that the validated FGFR2 event layer
  // (Cassette and IIIb/IIIc) is replaced by one optional exploratory tab.
  // Transcript/Exon comparison lives inline inside Exon map; Figure Gallery is
  // the global navigation page only (never duplicated here).
  // Canonical status: the domain/boundary layers are "pending" iff the validated
  // coordinate model (single source of truth, also read by the pages themselves)
  // is not yet available. This guarantees the tab badge and the page can never
  // disagree (Part B).
  const coordPending = ((model?.protein_coordinate_model?.models?.[0]?.status) ?? "pending") !== "available";
  // Cross-species MSA is a genuine multi-species product; it is shown as its own
  // top-level tab only when ≥2 species are present. Within-species isoform
  // alignment is NOT a top-level tab any more — it lives inside Candidate
  // Evidence (Part 6). The generic Evidence tab is removed (Part 5).
  const nSpecies = (model?.species || []).length;
  const multiSpecies = nSpecies > 1;
  const sharedTabs = [
    ["summary", "Summary"], ["isoforms", "Isoforms"], ["files", "Data & Downloads"],
    { id: "exon", label: "Exon map", separatorBefore: true },
    { id: "architecture", label: "Domain architecture", pending: coordPending },
    { id: "boundary", label: "Boundary", pending: coordPending },
    ...(multiSpecies ? [{ id: "msa", label: "MSA" }] : []),
    ["synteny", "Synteny"],
    // Optional exploratory extension (only for event_layer.type = exploratory).
    { id: "candidates", label: "Exploratory Candidate Evidence", separatorBefore: true },
  ];
  // Redirect removed generic tabs to a safe destination (Part 5): old Evidence /
  // Isoform Alignment URLs land on Summary / MSA / Candidate Evidence.
  const sharedIds = new Set(["summary", "isoforms", "files", "exon", "architecture",
    "boundary", "msa", "synteny", "candidates"]);
  const effTab = isValidated ? tab
    : (tab === "evidence" ? "summary"
      : tab === "alignment" ? (multiSpecies ? "msa" : "candidates")
        : sharedIds.has(tab) ? tab : "summary");
  const primary = pickPrimary(sp);
  return (
    <>
      <WorkspaceHeader
        title={sp.display_species_name}
        sub={isValidated ? sp.taxon_group : `${model?.analysis?.gene_symbol || "gene"} · gene-level analysis`}
        badges={isValidated ? <>
          <Badge cls={sp.overall_readiness}>{readinessLabel(sp.overall_readiness)}</Badge>
          {["IIIb", "IIIc"].map((iso) => header[iso] && (
            <span key={iso} className="ws-iso">
              <span className={`iso iso-${iso.toLowerCase()}`}>{iso}</span>
              <Badge cls={header[iso].readiness_class} soft>{header[iso].readiness_class === "review" ? "review" : "primary"}</Badge>
            </span>
          ))}
        </> : <>
          <Badge cls={sp.overall_status || "accepted"}>
            {readinessLabel(sp.overall_status || "accepted")}
          </Badge>
          {primary && <Badge cls="accepted" soft>primary selected</Badge>}
        </>}
      />

      <ExplorerTabs tabs={isValidated ? validatedTabs : sharedTabs} active={effTab} onSelect={setTab} />

      <div className="tab-body">
        {effTab === "summary" && <SummaryTab sp={sp} model={model} eventType={eventType} />}
        {effTab === "isoforms" && <IsoformsTab sp={sp} model={model} eventType={eventType} />}
        {effTab === "evidence" && isValidated && <EvidenceTab sp={sp} />}
        {effTab === "files" && <DataDownloads downloads={downloads} eventType={eventType}
                                             availability={model?.analysis_availability} />}
        {isValidated ? <>
          {effTab === "exon" && <CoordinateTrack preloaded={idx?.coordinates || {}} species={sp.species} embedded />}
          {effTab === "cassette" && <CassetteExplorer preloaded={idx?.cassette || {}} species={sp.species} embedded />}
          {effTab === "architecture" && <DomainArchitecture preloaded={idx?.domainArch || {}} species={sp.species} embedded />}
          {effTab === "boundary" && <BoundarySpeciesPanel sp={sp} idx={idx} />}
          {effTab === "msa" && <MsaExplorer preloaded={idx?.msa || {}} species={sp.species} embedded />}
          {effTab === "synteny" && <SyntenyViewer preloaded={idx?.synteny || {}} species={sp.species} embedded />}
        </> : <SharedTabBody sp={sp} tab={effTab} model={model} idx={idx} setPage={setPage} multiSpecies={multiSpecies} />}
      </div>
    </>
  );
}

const BOUNDARY_CLASS_META = [
  ["exact_edge", "Exact domain edge", "accepted"],
  ["near_edge", "Near a domain edge", "minor"],
  ["inside_domain", "Inside a domain", "info"],
  ["outside_domain", "Outside annotated domains", "neutral"],
  ["unknown", "Unavailable / uncertain", "neutral"],
];

// Compact, species-specific exon–domain boundary summary for the selected primary
// protein. Values recompute whenever the selected species (sp) changes.
function GenericBoundarySummary({ sp, model, setPage }) {
  const boundary = model?.boundary || {};
  const speciesProteins = (boundary.proteins || []).filter(
    (p) => (p.species_id || p.species) === (sp.species || sp.species_id));
  const boundaries = speciesProteins.flatMap((p) =>
    (p.boundaries || []).map((b) => ({ ...b, protein_id: p.protein_id })));
  const counts = boundaries.reduce((acc, b) => {
    acc[b.category] = (acc[b.category] || 0) + 1;
    return acc;
  }, {});
  const threshold = boundary.near_edge_threshold_aa ?? 5;
  const proteinLabel = speciesProteins.map((p) => p.protein_id).join(", ") || "—";

  if (!boundaries.length) {
    return <Empty title="No boundary summary for this species"
      hint="No internal coding-exon boundary was computed for the selected primary protein." />;
  }
  return (
    <div className="viewer">
      <div className="arch-note info">
        Exon–domain boundary summary for <b>{sp.display_species_name}</b> · primary protein
        <code> {proteinLabel}</code>. Distances are to the nearest representative InterPro
        domain edge (structural domains only; near-edge threshold {threshold} aa).
        Internal coding-exon junctions only.
      </div>
      <div className="kpi-grid">
        <div className="kpi"><div className="kpi-value">{boundaries.length}</div>
          <div className="kpi-label">Internal boundaries</div></div>
        {BOUNDARY_CLASS_META.map(([key, label, cls]) => (
          <div className="kpi" key={key}>
            <div className="kpi-value"><Badge cls={cls} soft>{counts[key] || 0}</Badge></div>
            <div className="kpi-label">{label}</div>
          </div>
        ))}
      </div>
      <div className="table-scroll"><table className="mini-tbl">
        <thead><tr><th>Boundary</th><th>Position (aa)</th><th>Nearest domain</th>
          <th>Edge</th><th>Distance (aa)</th><th>Class</th></tr></thead>
        <tbody>{boundaries.map((b) => (
          <tr key={`${b.protein_id}:${b.exon_boundary_id}`}>
            <td>{b.exon_boundary_id}</td>
            <td>{b.boundary_position_aa}</td>
            <td>{b.nearest_domain_name || b.nearest_domain_accession || b.nearest_domain_id || "—"}</td>
            <td>{b.nearest_edge || b.domain_edge_type || "—"}</td>
            <td>{b.absolute_distance_aa ?? b.distance_aa ?? "—"}</td>
            <td><Badge cls={(BOUNDARY_CLASS_META.find((m) => m[0] === b.category) || [])[2] || "neutral"} soft>
              {b.category}</Badge></td>
          </tr>
        ))}</tbody>
      </table></div>
      {setPage && (
        <p className="muted sm">
          <button className="btn ghost sm" onClick={() => setPage("boundary")}>
            Open the full Exon–Domain Boundaries page →</button>
        </p>
      )}
    </div>
  );
}

function SharedTabBody({ sp, tab, model, idx, setPage, multiSpecies }) {
  const pendingDomain = !model?.domain_architecture?.available;
  const pendingBoundary = !model?.boundary?.available;
  return <>
        {/* --- Shared FGFR2 base modules (same components + data contracts) --- */}
        {tab === "exon" && (
          idx?.coordinateModel?.models?.length
            ? <ExonMap model={idx.coordinateModel} species={sp.species} />
            : idx?.coordinates?.available
              ? <CoordinateTrack preloaded={idx.coordinates} species={sp.species} embedded
                  compareModels />
              : <Empty title="Exon map not available"
                  hint="No validated protein-coordinate model was built for this run." />
        )}
        {tab === "architecture" && (
          idx?.coordinateModel?.models?.length
            ? <ProteinArchitecture model={idx.coordinateModel} species={sp.species} />
            : pendingDomain
              ? <PendingAnalysisCard title="Domain architecture pending"
                  badge={<Badge cls="neutral" soft>pending cluster</Badge>}
                  description="Real InterProScan / pyTMHMM annotation has not completed for this gene. No
                    domain result is shown; run the cluster round-trip to populate this same module."
                  command={model?.analysis_stage?.cluster_command} />
              : <DomainArchitecture preloaded={model?.domain_architecture?.species_index
                  || model?.domain_architecture} species={sp.species} embedded />
        )}
        {tab === "boundary" && (
          idx?.coordinateModel?.models?.length
            ? <BoundaryExplorer model={idx.coordinateModel} species={sp.species} />
            : pendingBoundary
              ? <PendingAnalysisCard title="Boundary analysis pending"
                  badge={<Badge cls="neutral" soft>pending cluster</Badge>}
                  description="Exon–domain boundary distances are computed only after real domain
                    annotation. This same Boundary module is enriched once the cluster round-trip runs."
                  command={model?.analysis_stage?.cluster_command} />
              : <GenericBoundarySummary sp={sp} model={model} setPage={setPage} />
        )}
        {/* Cross-species MSA — top-level tab only for multi-species runs (Part 6C).
            One selected primary protein per species; within-species isoform
            alignment lives inside Candidate Evidence instead. */}
        {tab === "msa" && (
          multiSpecies && idx?.msa?.available
            ? <MsaExplorer preloaded={idx.msa} species={sp.species} embedded isoformMode />
            : <Empty title="Cross-species MSA not available"
                hint="A cross-species alignment needs ≥2 species with one primary protein each." />
        )}
        {tab === "synteny" && (
          idx?.synteny?.available
            ? <SyntenyViewer preloaded={idx.synteny} species={sp.species} embedded />
            : <Empty title="Synteny not available" />
        )}
        {/* --- Optional exploratory extension (replaces the FGFR2 event layer) --- */}
        {tab === "candidates" && <CandidateEvidencePanel model={model} idx={idx}
          multiSpecies={multiSpecies} sp={sp} />}
    </>;
}

// Real post-cluster domain context for a candidate, derived from the validated
// protein-coordinate model (Part 7). While the domain layer is not yet available
// this returns a "pending" marker so no stale or fabricated context is shown.
function candidateDomainContext(coordModel, speciesId, cand) {
  const models = coordModel?.models || [];
  const m = models.find((x) => (x.species_id || x.species) === speciesId) || models[0];
  if (!m) return { status: "unavailable" };
  if (m.status !== "available") return { status: "pending" };
  const domains = m.representative_domains || [];
  const s = cand.aa_start, e = cand.aa_end;
  const dstart = (d) => d.start ?? d.start_aa;
  const dend = (d) => d.end ?? d.end_aa;
  const dlabel = (d) => d.label || d.name || d.interpro_accession || "domain";
  const overlaps = domains
    .map((d) => ({ d, ds: dstart(d), de: dend(d),
      overlap: Math.min(e, dend(d)) - Math.max(s, dstart(d)) + 1 }))
    .filter((x) => x.overlap > 0)
    .sort((a, b) => b.overlap - a.overlap);
  const boundaries = (m.exon_boundaries || [])
    .filter((b) => b.protein_position >= s && b.protein_position <= e)
    .map((b) => b.exon_boundary_id || b.boundary_id || `aa ${b.protein_position}`);
  if (overlaps.length) {
    const top = overlaps[0];
    return { status: "available", overlapping: true,
      label: dlabel(top.d), accession: top.d.interpro_accession || top.d.accession || null,
      domainStart: top.ds, domainEnd: top.de, overlapLen: top.overlap, boundaries };
  }
  let nearest = null, best = Infinity;
  for (const d of domains) {
    const ds = dstart(d), de = dend(d);
    const dist = s > de ? s - de : ds > e ? ds - e : 0;
    if (dist < best) { best = dist; nearest = { d, ds, de, dist }; }
  }
  return { status: "available", overlapping: false,
    label: nearest ? dlabel(nearest.d) : null,
    accession: nearest?.d?.interpro_accession || null,
    nearestDist: nearest?.dist, boundaries };
}

function CandidateEvidencePanel({ model, idx, multiSpecies, sp }) {
  const selection = useScientificSelection();
  const coordModel = model?.protein_coordinate_model;
  const sidebarSpecies = sp?.species || selection?.selectedSpeciesId;
  const ceSpecies = model?.candidate_evidence?.species;
  const speciesOptions = (multiSpecies && Array.isArray(ceSpecies))
    ? ceSpecies.map((s) => ({ id: s.species_id || s.species,
        name: s.display_species_name || s.species_id || s.species,
        candidates: s.candidates || [] }))
    : [];
  const [activeSpecies, setActiveSpecies] = useState(sidebarSpecies);
  const [localSelId, setLocalSelId] = useState(null);
  const activeIsSidebar = !multiSpecies || activeSpecies === sidebarSpecies;

  const candidates = activeIsSidebar
    ? (selection?.rankedCandidates || [])
    : (speciesOptions.find((o) => o.id === activeSpecies)?.candidates || []);
  const isoformRows = model?.shared_indices?.gene_explorer?.isoforms
    || model?.protein_models || [];
  const activeIsoforms = isoformRows
    .filter((row) => (row.species_id || row.species) === activeSpecies);
  const activeProteinCount = new Set(activeIsoforms
    .map((row) => row.protein_id)
    .filter(Boolean)).size;
  const activePrimaryProtein = (activeIsoforms.find((row) =>
    row.primary_status === "primary" || row.is_primary)?.protein_id
    || activeIsoforms[0]?.protein_id || "");
  if (!candidates.length && !multiSpecies) {
    return <Empty title="Exploratory Candidate Evidence unavailable" />;
  }
  const selectedId = activeIsSidebar ? selection?.selectedCandidateId : localSelId;
  const selected = candidates.find((c) => c.candidate_id === selectedId) || candidates[0];
  const onSelect = (c) => activeIsSidebar ? selection?.selectCandidate(c) : setLocalSelId(c.candidate_id);
  const refProtein = selected?.reference_protein
    || (activeIsSidebar ? selection?.primaryProteinId : activePrimaryProtein);
  const coord = selection?.model?.shared_indices?.coordinate_track_index
    || selection?.model?.shared_event_indices?.coordinates;
  const coordSpecies = (coord?.species || []).find((s) => (s.species === activeSpecies || s.species_id === activeSpecies))
    || coord?.species?.[0];
  const primaryPanel = coordSpecies?.panels?.primary;
  const proteinLength = primaryPanel?.protein_length || coord?.protein_length
    || Math.max(...candidates.map((c) => c.aa_end || 0), 1);
  const exonBoundaries = (primaryPanel?.blocks || []).flatMap((b) => [b.start, b.end])
    .filter((v) => v != null);
  const domainCtx = selected ? candidateDomainContext(coordModel, activeSpecies, selected) : null;

  return (
    <div className="candidate-master-detail">
      <div className="arch-note info">
        Exploratory evidence only — no candidate is a validated splicing event. Reference protein:
        <code> {refProtein || "n/a"}</code> (canonical primary). Candidates are ranked by a transparent
        evidence score; the top-ranked candidate is selected by default and drives every linked view.
      </div>
      {multiSpecies && speciesOptions.length > 0 && (
        <div className="msa-iso-filter">
          <span className="muted small">Species:</span>
          {speciesOptions.map((o) => (
            <button key={o.id} className={`chip sm${activeSpecies === o.id ? " sel" : ""}`}
              onClick={() => { setActiveSpecies(o.id); setLocalSelId(null); }}>{o.name}</button>
          ))}
        </div>
      )}
      {!candidates.length ? (
        // An empty list can mean one available model or several identical proteins.
        <Empty title="No protein-isoform difference block was detected."
          hint={activeProteinCount === 1
            ? "Only one protein isoform was available for this species, so no within-species protein-isoform comparison could be made."
            : "The compared isoforms of this species encode identical protein sequences, so there is no difference block to explore. Nothing was removed by a length threshold."} />
      ) : (
      <div className="cand-page">
        {/* A — compact horizontal candidate strip. The ranking no longer reserves
            half the page width, so the alignment can use the full content width. */}
        <CandidateStrip candidates={candidates} selected={selected} onSelect={onSelect} />

        {/* B — full-width analysis workspace for the selected candidate */}
        {selected && <CandidateWorkspace c={selected} proteinLength={proteinLength}
          exonBoundaries={exonBoundaries} candidates={candidates} domainCtx={domainCtx}
          model={model} idx={idx} activeSpecies={activeSpecies} />}
      </div>
      )}

    </div>
  );
}

// Compact, horizontally scrollable candidate selector (Part 1A). Shows only the
// ranking essentials; every further detail appears after selection.
function CandidateStrip({ candidates, selected, onSelect }) {
  return (
    <div className="cand-strip-wrap">
      <div className="cand-strip-head">
        <b className="cand-strip-title">Exploratory candidates</b>
        <span className="muted small">{candidates.length} ranked by evidence score · select one to analyse</span>
      </div>
      <div className="cand-strip" role="tablist" aria-label="Exploratory candidate ranking">
        {candidates.map((c) => {
          const isSel = selected?.candidate_id === c.candidate_id;
          const affected = c.protein_isoform_evidence?.affected_proteins?.length ?? 0;
          return (
            <button key={c.candidate_id} role="tab" aria-selected={isSel}
              className={`cand-chipcard${isSel ? " sel" : ""}`} onClick={() => onSelect(c)}>
              <span className="ccc-rank">{c.rank_label}</span>
              <span className="ccc-region">aa {c.aa_start}–{c.aa_end}</span>
              <span className="ccc-meta">{c.length} aa · {affected} isoform{affected === 1 ? "" : "s"}</span>
              <span className="ccc-score">{c.overall_score}<i>/100</i></span>
            </button>
          );
        })}
      </div>
    </div>
  );
}

// Adapt one species entry of the shared within-species isoform-alignment index
// into the payload shape MsaExplorer's isoform view already consumes, so every
// species gets its own alignment without a second alignment renderer (Part 14).
// Returns null when the species has no alignment product, or
// { singleModel: true } when only one valid protein model exists.
function withinSpeciesAlignment(model, speciesId) {
  const index = model?.isoform_alignment;
  const entries = index?.species || [];
  if (!entries.length) return null;
  const entry = entries.find((s) => (s.species_id || s.species) === speciesId)
    || (entries.length === 1 ? entries[0] : null);
  if (!entry || entry.status !== "available") return null;
  const seqs = entry.sequences || [];
  if (seqs.length < 2) return { singleModel: true };

  // Row metadata (curation status, protein length) comes from the validated
  // coordinate model rather than being guessed from the accession prefix.
  const coord = model?.protein_coordinate_model?.models || [];
  const coordSpecies = coord.find((m) => (m.species_id || m.species) === speciesId) || coord[0];
  const byProtein = new Map();
  for (const t of coordSpecies?.transcript_models || []) {
    if (t.protein_id) byProtein.set(t.protein_id, t);
  }

  const speciesKey = entry.species_id || entry.species || speciesId;
  const rows = seqs.map((s) => {
    const aligned = s.aligned_sequence || s.seq || "";
    const tm = byProtein.get(s.protein_id);
    return {
      protein_id: s.protein_id,
      transcript_id: s.transcript_id || tm?.transcript_id || null,
      is_primary: Boolean(s.is_primary),
      seq: aligned,
      species: speciesKey,
      display_species_name: entry.display_species_name || speciesKey,
      curation_status: tm?.curation_status || null,
      protein_length: tm?.protein_length ?? (aligned.replace(/-/g, "").length || null),
      exon_count: (tm?.blocks || []).length || null,
    };
  });
  return {
    available: true,
    mode: "isoform_alignment",
    disclaimer: index.disclaimer
      || "This is a protein-isoform alignment within one species, not a cross-species conservation analysis.",
    tabs: [{ key: "isoform", label: "Within-species isoform alignment" }],
    alignments: {
      isoform: {
        available: true,
        label: "Within-species isoform alignment",
        rows,
        n_columns: entry.alignment_length || rows[0].seq.length,
        sequence_count: entry.sequence_count || rows.length,
        reference_sequence: entry.reference_sequence || null,
        tool: index.tool || "MAFFT",
        file: entry.alignment_file || null,
      },
    },
  };
}

// Within-species isoform alignment embedded in Candidate Evidence (Parts 13/14).
// Every species with at least two valid protein models gets its own alignment,
// including inside multi-species datasets, because a cross-species MSA does not
// replace within-species isoform analysis.
function IntegratedIsoformAlignment({ model, idx, activeSpecies, focusCandidate }) {
  const perSpecies = withinSpeciesAlignment(model, activeSpecies);
  // Fall back to the run-level MSA index only when it really is a
  // within-species isoform alignment (single-species runs).
  const msa = idx?.msa;
  const fallback = (!perSpecies && msa?.available
    && (msa.mode === "isoform_alignment" || msa?.alignments?.isoform?.available))
    ? msa : null;
  const payload = (perSpecies && !perSpecies.singleModel) ? perSpecies : fallback;

  return (
    <div className={`cand-alignment${focusCandidate ? " focused" : " full-width"}`}>
      <div className="cand-aln-head">
        <h4 className="ev-section">{focusCandidate
          ? "Candidate-focused alignment"
          : "Full Isoform Alignment · within species"}</h4>
        {payload && <AlignmentExportMenu payload={payload} focusCandidate={focusCandidate} />}
      </div>
      {payload ? (
        <MsaExplorer key={`${activeSpecies}-${focusCandidate ? "focus" : "full"}`}
          preloaded={payload} species={activeSpecies}
          embedded isoformMode focusCandidate={focusCandidate} />
      ) : (
        <p className="muted sm">Isoform alignment unavailable: only one protein model is available.</p>
      )}
    </div>
  );
}

function prettyClass(cls) {
  return String(cls || "").replace(/_/g, " ") || "isoform difference";
}

function evidenceStrength(score) {
  if (score >= 60) return "strong";
  if (score >= 35) return "moderate";
  return "weak";
}

// Mini full-protein location plot for the selected candidate (PART 7B).
function CandidateLocationPlot({ c, proteinLength, exonBoundaries, candidates }) {
  const W = 900, PAD = 8, len = Math.max(1, proteinLength);
  const x = (aa) => PAD + (Math.max(0, aa) / len) * (W - 2 * PAD);
  return (
    <svg className="cand-mini-plot" viewBox={`0 0 ${W} 60`} preserveAspectRatio="xMidYMid meet"
      role="img" aria-label={`Location of ${c.rank_label} on the primary protein`}>
      <line x1={PAD} y1="34" x2={W - PAD} y2="34" stroke={CHROME.rule} />
      {/* nearest exon boundaries */}
      {exonBoundaries.map((pos, i) => (
        <line key={`b${i}`} x1={x(pos)} y1="28" x2={x(pos)} y2="40"
          stroke={featureStyle("coding_exon").stroke} strokeWidth="0.6" />
      ))}
      {/* other candidates faint */}
      {candidates.filter((o) => o.candidate_id !== c.candidate_id).map((o) => (
        <rect key={o.candidate_id} x={x(o.aa_start)} y="22" width={Math.max(2, x(o.aa_end) - x(o.aa_start))}
          height="24" rx="2" {...featureProps("candidate_region", { faint: true })}>
          <title>{o.rank_label} · aa {o.aa_start}–{o.aa_end}</title></rect>
      ))}
      {/* selected candidate strong */}
      <rect x={x(c.aa_start)} y="18" width={Math.max(2, x(c.aa_end) - x(c.aa_start))} height="32" rx="3"
        {...featureProps("candidate_region")} />
      <text x={x(c.aa_start)} y="14" {...textProps("candidateLabel")}>{c.rank_label} · {c.aa_start}–{c.aa_end}</text>
      <text x={PAD} y="58" {...textProps("axis")}>1</text>
      <text x={W - PAD} y="58" textAnchor="end" {...textProps("axis")}>{len} aa</text>
    </svg>
  );
}

// Compact scientific summary cards for a candidate. Generic (reads only the
// candidate evidence contract) and reused by both the Evidence tab and the
// Exploratory Candidate Evidence detail so the two views stay connected.
function EvCard({ label, value, hint, tone }) {
  return (
    <div className={`ev-card${tone ? ` ev-${tone}` : ""}`}>
      <div className="ev-card-label">{label}</div>
      <div className="ev-card-value">{value}</div>
      {hint && <div className="ev-card-hint">{hint}</div>}
    </div>
  );
}

function CandidateEvidenceCards({ c, domainCtx }) {
  const iso = c.protein_isoform_evidence || {};
  const ex = c.exon_evidence || {};
  const al = c.alignment_evidence || {};
  const affProteins = iso.affected_proteins?.length || 0;
  const affTranscripts = iso.affected_transcripts?.length || 0;
  const pairs = iso.supporting_isoform_pairs?.length || 0;
  const exonLabels = (ex.transcript_exon_numbers || []).map((n) => `E${n}`).join(", ");
  const alignedIso = al.supporting_aligned_isoforms || 0;
  // Real domain context (Part 7): overlap/nearest come from the post-cluster
  // coordinate model; pending only when the domain layer is genuinely not ready.
  const dc = domainCtx || { status: "pending" };
  const domValue = dc.status === "available"
    ? (dc.overlapping ? `${dc.label} (${dc.overlapLen} aa)` : (dc.label ? `nearest ${dc.label}` : "no domain"))
    : dc.status === "pending" ? "pending" : "unavailable";
  const domHint = dc.status === "available"
    ? (dc.overlapping ? `overlaps aa ${dc.domainStart}–${dc.domainEnd}${dc.accession ? ` · ${dc.accession}` : ""}`
      : (dc.nearestDist != null ? `${dc.nearestDist} aa away` : "domain context"))
    : dc.status === "pending" ? "after InterProScan/pyTMHMM" : "no coordinate model";
  return (
    <div className="ev-cards">
      <EvCard label="Affected isoforms" value={affProteins || "—"}
        hint={`${affTranscripts || 0} transcript(s)`} tone={affProteins ? "ok" : null} />
      <EvCard label="Affected exons" value={exonLabels || "—"}
        hint={ex.exon_aligned ? "exon-aligned" : "not exon-aligned"} tone={ex.exon_aligned ? "ok" : null} />
      <EvCard label="Protein coordinates" value={`aa ${c.aa_start}–${c.aa_end}`}
        hint={`${c.length} aa · ${c.reference_protein}`} />
      <EvCard label="Alignment support" value={alignedIso || "—"}
        hint={alignedIso ? `cols ${al.alignment_start ?? c.aa_start}–${al.alignment_end ?? c.aa_end}` : "no aligned isoform"}
        tone={alignedIso ? "ok" : null} />
      <EvCard label="Protein-model support" value={pairs || "—"}
        hint="pairwise comparison(s)" tone={pairs ? "ok" : null} />
      <EvCard label="Biological validation" value="not validated"
        hint="exploratory candidate" tone="neutral" />
      <EvCard label="Domain context" value={domValue} hint={domHint}
        tone={dc.status === "available" && dc.overlapping ? "ok" : dc.status === "pending" ? "pending" : "neutral"} />
    </div>
  );
}

// One compact Export menu for the alignment (Part 5). SVG, PDF and PNG are three
// download formats of the same figure source, never separate figures.
function AlignmentExportMenu({ payload, focusCandidate }) {
  const selection = useScientificSelection();
  const aln = payload?.alignments?.isoform || {};
  // Same row order as the interactive view and the figures: primary first.
  const rows = useMemo(() => orderRows(aln.rows || []), [aln.rows]);
  const nCols = aln.n_columns || rows[0]?.seq?.length || 0;
  const primary = rows.find((r) => r.is_primary) || rows[0];
  const gene = selection?.model?.dataset_info?.gene_symbol
    || selection?.model?.gene_symbol || "gene";
  const species = rows[0]?.display_species_name || rows[0]?.species || "";
  const cand = selection?.selectedCandidate;
  const candBand = cand && primary
    ? [aaToColumn(primary.seq, cand.aa_start), aaToColumn(primary.seq, cand.aa_end)]
    : null;
  const candidateForFig = (candBand && candBand[0] != null && candBand[1] != null)
    ? { label: cand.rank_label, aa_start: cand.aa_start, aa_end: cand.aa_end,
        col_start: Math.min(...candBand), col_end: Math.max(...candBand) }
    : null;

  const common = { rows, nCols, gene, species, primaryId: primary?.protein_id,
    transcriptId: primary?.transcript_id };
  const candidates = candidateForFig ? [candidateForFig] : [];
  // Coding exons of the primary protein, in protein coordinates, so the candidate
  // figure can state which exon each affected residue belongs to. Read from the
  // validated coordinate model — never derived from the alignment.
  const exons = useMemo(() => {
    const models = selection?.model?.protein_coordinate_model?.models || [];
    const m = models.find((x) => x.protein_id === primary?.protein_id) || models[0];
    return (m?.exons || []).map((e) => ({ label: e.label, start: e.start, end: e.end }));
  }, [selection?.model, primary?.protein_id]);

  // Every format is built from a figure specification, so the SVG, the PDF and the
  // PNG are three views of one figure instead of three separate pipelines.
  const buildFigure = () => (focusCandidate && candidateForFig
    ? candidateAlignmentFigureSpec({ ...common, candidate: candidateForFig, exons,
        affected: selection?.affectedProteinsFor(cand.candidate_id) || null })
    : alignmentOverviewFigureSpec({ ...common, candidates, tool: aln.tool || "MAFFT" }));

  const stem = focusCandidate
    ? `isoform_alignment_${primary?.protein_id || "primary"}_${cand?.rank_label || "candidate"}`
    : `isoform_alignment_${primary?.protein_id || "primary"}_overview`;
  const label = focusCandidate ? "Candidate alignment" : "Full Isoform Alignment";

  const exportSvg = () => downloadFigureSvg(buildFigure(), stem);
  const exportPdf = () => downloadFigurePdf(buildFigure(), stem);
  const exportPng = () => downloadFigurePng(buildFigure(), stem);

  // The residue-level alignment is far too long for one page, so it is wrapped
  // into blocks and spread over as many pages as the layout needs. The PDF is one
  // document; SVG and PNG have no multi-page form, so those write one file per page.
  const wrappedStem = `isoform_alignment_${primary?.protein_id || "primary"}_residues`;
  const buildWrapped = () => wrappedAlignmentFigureSpecs({
    ...common, candidates, tool: aln.tool || "MAFFT" });
  const exportWrappedPdf = () => downloadFigurePdfPages(buildWrapped(), wrappedStem);
  const exportWrappedSvg = () => buildWrapped()
    .forEach((f, i) => downloadFigureSvg(f, `${wrappedStem}_p${i + 1}`));
  const exportWrappedPng = () => buildWrapped()
    .forEach((f, i) => downloadFigurePng(f, `${wrappedStem}_p${i + 1}`));

  const exportTsv = () => downloadFigureTsv(
    [ALIGNMENT_SUMMARY_COLUMNS.join("\t"),
      ...alignmentSummaryRows({ rows, nCols, candidates }).map((r) => r.join("\t"))].join("\n"),
    stem);
  const exportFasta = () => downloadBlob(
    new Blob([alignmentFasta(rows)], { type: "text/plain;charset=utf-8" }),
    `${stem}.fasta`);

  return (
    <Menu label="Export" title="Export alignment figure and data" align="right">
      <button className="menu-item" onClick={exportSvg}>{label} — SVG (vector)</button>
      <button className="menu-item" onClick={exportPdf}>{label} — PDF (vector)</button>
      <button className="menu-item" onClick={exportPng}>{label} — PNG (300 dpi)</button>
      <div className="menu-sep" />
      <button className="menu-item" onClick={exportWrappedPdf}>
        Residue-level alignment — multi-page PDF (vector)</button>
      <button className="menu-item" onClick={exportWrappedSvg}>
        Residue-level alignment — SVG per page</button>
      <button className="menu-item" onClick={exportWrappedPng}>
        Residue-level alignment — PNG per page (300 dpi)</button>
      <div className="menu-sep" />
      {aln.file
        ? <a className="menu-item" href={fileUrl(aln.file)}>Alignment FASTA</a>
        : <button className="menu-item" onClick={exportFasta}>Alignment FASTA</button>}
      <button className="menu-item" onClick={exportTsv}>
        {focusCandidate ? "Candidate-region TSV" : "Alignment summary TSV"}</button>
    </Menu>
  );
}

// Raw pairwise comparisons are supporting bookkeeping, not the conclusion of the
// page, so they stay collapsed behind a compact summary (Part 4).
function PairwiseComparisonDetails({ pairs, c }) {
  if (!pairs.length) return null;
  const supporting = pairs.filter((p) => {
    const s = p.region_start_aa, e = p.region_end_aa;
    return s != null && e != null && s <= c.aa_end && e >= c.aa_start;
  }).length;
  const types = [...new Set(pairs.map((p) => p.event_type_candidate).filter(Boolean))];
  const exonAligned = pairs.filter((p) => p.exon_aligned === true
    || String(p.exon_aligned).toLowerCase() === "true").length;
  return (
    <details className="tech-details pairwise-details">
      <summary>Pairwise comparison details</summary>
      <p className="muted sm pairwise-summary">
        <b>{pairs.length}</b> pairwise comparison{pairs.length === 1 ? "" : "s"} ·
        {" "}<b>{supporting}</b> overlapping {c.rank_label} ·
        {" "}type{types.length === 1 ? "" : "s"} {types.map(prettyClass).join(", ") || "—"} ·
        {" "}<b>{exonAligned}</b> exon-aligned
      </p>
      <table className="cand-score-table">
        <thead><tr><th>Protein A</th><th>Protein B</th><th>Region (aa)</th><th>Type</th><th>Exon-aligned</th></tr></thead>
        <tbody>
          {pairs.map((p, i) => (
            <tr key={i}>
              <td>{p.protein_a}</td><td>{p.protein_b}</td>
              <td>{p.region_start_aa}–{p.region_end_aa}</td>
              <td>{prettyClass(p.event_type_candidate)}</td>
              <td>{String(p.exon_aligned)}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </details>
  );
}

// The Full Isoform Alignment is a first-class analysis view, not a subsection of
// the candidate's isoform evidence, because it describes every protein isoform of
// the species independently of any single candidate.
const CAND_SECTIONS = [
  ["overview", "Candidate overview"],
  ["alignment", "Full Isoform Alignment"],
  ["focused", "Candidate-focused alignment"],
  ["context", "Exon & domain context"],
];

// Candidate workspace: a stable header (id, interval, score, strength, validation)
// plus four full-width analysis views.
function CandidateWorkspace({ c, proteinLength, exonBoundaries, candidates, domainCtx,
  model, idx, activeSpecies }) {
  const [sect, setSect] = useState("overview");
  const sc = c.score_components || {};
  // Score components shown to the user. Generic external (UniProt) annotation is
  // intentionally excluded from the exploratory UI (Part 2); it remains in the
  // backend score and downloadable tables for reproducibility.
  const scoreRows = [
    ["Isoform support", sc.isoform_support_score, "Distinct isoform pairs showing the region"],
    ["Exon boundary", sc.exon_boundary_score, "Region boundaries match exon boundaries"],
    ["Alignment", sc.alignment_score, "Region is supported in the isoform alignment"],
    ["Domain context", sc.domain_context_score, "Overlap with a representative InterPro domain"],
    ["Penalty", sc.penalty != null ? -Math.abs(sc.penalty) : null, "Deductions for weak/conflicting evidence"],
  ];
  const section = (title, tone, children) => (
    <div className="cand-source">
      <div className="cand-source-head">
        <b>{title}</b>
        <Badge cls={tone === "supported" ? "accepted" : tone === "pending" ? "neutral" : "minor"} soft>
          {tone === "supported" ? "supported" : tone === "pending" ? "pending cluster" : "neutral"}
        </Badge>
      </div>
      <div className="muted small">{children}</div>
    </div>
  );
  const ex = c.exon_evidence || {};
  const al = c.alignment_evidence || {};
  const iso = c.protein_isoform_evidence || {};
  const pairs = iso.supporting_isoform_pairs || [];
  const strength = evidenceStrength(c.overall_score);

  return (
    <div className="cand-workspace">
      {/* Selected-candidate header (Part 1B) — clearly exploratory, never framed
          as a confirmed splice event. */}
      <div className="cand-sel-head">
        <div className="csh-main">
          <span className="csh-rank">{c.rank_label}</span>
          <div className="csh-titles">
            <h3>aa {c.aa_start}–{c.aa_end} on {c.reference_protein}</h3>
            <span className="muted small">{c.length} aa · {prettyClass(c.candidate_class)}
              {" "}· primary-protein coordinates</span>
          </div>
        </div>
        <div className="csh-badges">
          <Badge cls="accepted" soft>Evidence score: {c.overall_score}/100</Badge>
          <Badge cls={strength === "strong" ? "accepted" : strength === "moderate" ? "minor" : "neutral"} soft>
            Evidence strength: {strength}</Badge>
          <Badge cls="minor" soft>Biological validation: not validated</Badge>
        </div>
      </div>

      <div className="seg cand-sect-seg">
        {CAND_SECTIONS.map(([id, label]) => (
          <button key={id} className={`seg-btn${sect === id ? " on" : ""}`}
            onClick={() => setSect(id)}>{label}</button>
        ))}
      </div>

      {sect === "overview" && (
        <>
          <CandidateLocationPlot c={c} proteinLength={proteinLength}
            exonBoundaries={exonBoundaries} candidates={candidates} />

          <CandidateEvidenceCards c={c} domainCtx={domainCtx} />

          <div className="cand-sources">
            {section("Observation", "supported",
              `An isoform-difference region at aa ${c.aa_start}–${c.aa_end} (${c.length} aa) on the canonical `
              + `primary protein ${c.reference_protein}. ${c.interpretation || ""}`)}
            {section("Protein-model evidence", pairs.length > 0 ? "supported" : "neutral",
              `${pairs.length} pairwise comparison(s); ${iso.affected_proteins?.length || 0} proteins / `
              + `${iso.affected_transcripts?.length || 0} transcripts affected. `
              + `Classification: ${prettyClass(iso.classification)}.`)}
          </div>

          <details className="tech-details">
            <summary>Transparent score calculation (overall {c.overall_score}/100)</summary>
            <table className="cand-score-table">
              <tbody>
                {scoreRows.map(([label, val, desc]) => (
                  <tr key={label}>
                    <td>{label}</td>
                    <td className="cand-score-val">{val == null ? "—" : val}</td>
                    <td className="muted small">{desc}</td>
                  </tr>
                ))}
                <tr className="cand-score-total">
                  <td><b>Evidence score</b></td><td className="cand-score-val"><b>{c.overall_score}/100</b></td>
                  <td className="muted small">Evidence strength: {strength}</td>
                </tr>
              </tbody>
            </table>
          </details>

          <p className="cand-interp"><b>Interpretation:</b> {c.confidence_reason}
            {" "}This is exploratory evidence only; the region is <b>not a validated event</b>.</p>
        </>
      )}

      {/* Full Isoform Alignment — an independent analysis of every protein isoform
          of the species, useful even with no candidate selected (Part 2). */}
      {sect === "alignment" && (
        <div className="cand-aln-view">
          <p className="muted sm">Complete within-species alignment of all valid protein models.
            This view is independent of the selected candidate; {c.rank_label} is shown only as an
            optional overlay.</p>
          <IntegratedIsoformAlignment model={model} idx={idx} activeSpecies={activeSpecies} />
        </div>
      )}

      {/* Candidate-focused alignment — the complementary, candidate-specific view. */}
      {sect === "focused" && (
        <div className="cand-aln-view">
          <p className="muted sm">Which isoforms support {c.rank_label}, and how do their aligned
            residues differ? Restricted to the candidate interval; the full alignment stays
            available in <b>Full Isoform Alignment</b>.</p>
          <div className="cand-sources">
            {section("Alignment evidence", (al.supporting_aligned_isoforms || 0) > 0 ? "supported" : "neutral",
              `${al.supporting_aligned_isoforms || 0} aligned isoform(s) support columns `
              + `${al.alignment_start ?? c.aa_start}–${al.alignment_end ?? c.aa_end}; pattern `
              + `${prettyClass(al.gap_divergence_pattern)}.`)}
          </div>
          <IntegratedIsoformAlignment model={model} idx={idx} activeSpecies={activeSpecies}
            focusCandidate />
          <PairwiseComparisonDetails pairs={pairs} c={c} />
        </div>
      )}

      {sect === "context" && (
        <div className="cand-sources">
          {section("Exon evidence", ex.exon_aligned === true ? "supported" : "neutral",
            `Supporting exons ${(ex.transcript_exon_numbers || []).map((n) => `E${n}`).join(", ") || "—"} `
            + `(${(ex.exon_ids || []).join(", ") || "no exon id"}); start boundary `
            + `${ex.start_boundary_matches_exon_boundary ? "matches" : "does not match"}, end boundary `
            + `${ex.end_boundary_matches_exon_boundary ? "matches" : "does not match"} an exon boundary.`)}
          {section("Domain context",
            domainCtx?.status === "pending" ? "pending"
              : (domainCtx?.overlapping ? "supported" : "neutral"),
            domainCtx?.status === "available"
              ? (domainCtx.overlapping
                ? `Overlaps the representative domain ${domainCtx.label} (aa ${domainCtx.domainStart}–`
                  + `${domainCtx.domainEnd}${domainCtx.accession ? `, ${domainCtx.accession}` : ""}); `
                  + `overlap ${domainCtx.overlapLen} aa`
                  + (domainCtx.boundaries?.length ? `; internal exon boundaries: ${domainCtx.boundaries.join(", ")}.` : ".")
                : (domainCtx.label
                  ? `No overlap with a representative domain; nearest is ${domainCtx.label}`
                    + (domainCtx.nearestDist != null ? ` (${domainCtx.nearestDist} aa away).` : ".")
                  : "No representative domain annotated on this protein."))
              : domainCtx?.status === "pending"
                ? "Domain overlap is computed only after real InterProScan / pyTMHMM results arrive."
                : "No coordinate model available for this candidate.")}
          {section("Boundary relationship",
            domainCtx?.boundaries?.length ? "supported" : "neutral",
            domainCtx?.boundaries?.length
              ? `Internal coding-exon boundaries inside the candidate interval: ${domainCtx.boundaries.join(", ")}.`
              : "No internal coding-exon boundary falls inside this candidate interval.")}
        </div>
      )}
    </div>
  );
}

function normalizeProtein(p) {
  const isPrimary = p.is_primary === true
    || p.primary_status === "primary" || p.role === "primary";
  return {
    ...p,
    protein_length: p.protein_length ?? p.length_aa,
    is_primary: isPrimary,
    primary_status: isPrimary ? "primary" : "alternative",
    role: isPrimary ? "primary" : "alternative",
  };
}

// Single canonical primary protein for a species. The canonical model always
// exposes selected_primary_protein (from primary_selection_report.json); fall
// back to per-protein flags only if it is absent.
function pickPrimary(sp) {
  const proteins = sp?.proteins || sp?.protein_models || sp?.isoforms || [];
  const selectedId = sp?.selected_primary_protein;
  if (selectedId) {
    const match = proteins.find((p) => p.protein_id === selectedId);
    if (match) return match;
  }
  return proteins.find((p) => p.is_primary === true
    || p.role === "primary" || p.primary_status === "primary") || proteins[0];
}

function sharedStatusClass(status) {
  const value = String(status || "").toLowerCase();
  if (["accepted", "available", "complete", "found"].includes(value)) return "accepted";
  if (["exploratory", "pending", "pending_cluster", "medium"].includes(value)) return "minor";
  if (["failed", "excluded", "not_found"].includes(value)) return "excluded";
  return "neutral";
}

function BoundarySpeciesPanel({ sp, idx }) {
  const matrix = idx?.boundaryMatrix;
  const arch = idx?.domainArch;
  if (!matrix?.available) {
    return <Empty title="Boundary consistency not available"
      hint="The Module 1 boundary-consistency analysis (step 16) was not found for this run." />;
  }
  const rows = (matrix.rows || []).filter((r) => r.species === sp.species);
  if (!rows.length) {
    return <Empty title="No cassette boundary data" hint={`No cassette boundary rows for ${sp.display_species_name}.`} />;
  }
  const panelFor = (iso) => arch?.species?.find((s) => s.species === sp.species)?.panels?.[iso] || null;
  return (
    <div className="bc-species-panel">
      <p className="muted sm">Cassette start/end boundaries vs. the nearest protein-domain boundary. Full explorer: Boundary Consistency tab.</p>
      {rows.map((r) => (
        <BoundaryDetailTrack key={r.isoform} pd={panelFor(r.isoform)} row={r} isoform={r.isoform} embedded />
      ))}
    </div>
  );
}

const LAYER_ORDER = { neutral: 0, accepted: 1, minor: 2, review: 3, excluded: 4 };

function mergeLayer(sp, key) {
  const cells = sp.isoforms.map((i) => i.layers?.[key]).filter(Boolean);
  if (!cells.length) return null;
  return cells.slice().sort((a, b) => (LAYER_ORDER[b.class] || 0) - (LAYER_ORDER[a.class] || 0))[0];
}

function SummaryTab({ sp, model, eventType }) {
  if (eventType !== "validated") {
    const proteins = sp.proteins || sp.isoforms || [];
    const primary = pickPrimary(sp);
    const cap = model?.analysis_stage || model?.dataset_summary?.capability || {};
    const multiSpecies = (model?.species || []).length > 1;
    // Canonical run stage from the validated coordinate model (single source of
    // truth), resolved for the *selected* species so a multi-species Summary can
    // never report another species' state.
    const coordModel = model?.protein_coordinate_model;
    const coordModels = coordModel?.models || [];
    const coordSpecies = coordModels.find((m) => (m.species_id || m.species) === sp.species)
      || coordModels[0];
    const stage = datasetStatusFromModel(coordModel, model?.dataset_info);
    const domainsReady = coordSpecies?.status === "available";
    const domVal = domainsReady ? "available" : (cap.domain_architecture || "pending_cluster");
    const bndVal = domainsReady ? "available" : (cap.exon_domain_boundaries || "pending_cluster");
    const missing = domainsReady ? [] : missingPostClusterItems(coordModel)
      .filter((e) => e.species === sp.species || coordModels.length === 1);
    const alignLabel = multiSpecies ? "Cross-species MSA" : "Isoform alignment";
    const items = [
      { id: "models", label: "Protein models", value: `${proteins.length}`, cls: "accepted", soft: true },
      { id: "exon", label: "Exon mapping", value: sp.exon_map_status || cap.exon_map || "available", cls: sharedStatusClass(sp.exon_map_status || "available"), soft: true },
      { id: "msa", label: alignLabel, value: sp.msa_status || (model?.shared_indices?.msa_index?.available ? "available" : "available"), cls: sharedStatusClass("available"), soft: true },
      { id: "syn", label: "Synteny", value: sp.synteny_status || cap.synteny || "available", cls: sharedStatusClass("available"), soft: true },
      { id: "dom", label: "InterProScan", value: domVal, cls: sharedStatusClass(domVal), soft: true },
      { id: "tm", label: "pyTMHMM", value: domVal, cls: sharedStatusClass(domVal), soft: true },
      { id: "bnd", label: "Boundary", value: bndVal, cls: sharedStatusClass(bndVal), soft: true },
    ];
    // Only render metadata cards that carry a real value (hide "Assembly: —",
    // empty scaffold, etc. — Part 3).
    const cards = [
      ["Gene", "accepted", model?.analysis?.gene_symbol],
      ["Species", "accepted", sp.display_species_name || sp.species],
      ["Common name", "neutral", sp.common_name],
      ["Taxonomic group", "neutral", sp.taxon_group || sp.taxonomic_group],
      ["Assembly", "neutral", sp.assembly_accession],
      ["Primary protein", "accepted", primary?.protein_id || model?.selected_primary_protein],
      ["Primary transcript", "accepted", primary?.transcript_id || model?.selected_primary_transcript],
      ["Protein length", "accepted",
        (primary?.length_aa || primary?.protein_length) ? `${primary.length_aa || primary.protein_length} aa` : null],
    ].filter(([, , v]) => v != null && v !== "" && v !== "—");
    return (
      <div className="summary-tab">
        <div className="summary-cards">
          <SummaryCard label="Dataset status" cls={DATASET_STATUS_META[stage]?.[0] || "neutral"}
            value={datasetStatusLabel(stage)} />
          {cards.map(([label, cls, value]) => (
            <SummaryCard key={label} label={label} cls={cls} value={value} />
          ))}
        </div>
        <EvidenceSummary items={items} />
        {missing.length > 0 && (
          <p className="muted sm">Still missing for a complete post-cluster result:
            {" "}{missing.map((e) => e.missing.join(", ")).join("; ")}.
            {" "}The cluster round-trip is started in <b>My Runs</b>.</p>
        )}
      </div>
    );
  }
  const iiib = sp.isoforms.find((x) => x.isoform === "IIIb");
  const iiic = sp.isoforms.find((x) => x.isoform === "IIIc");
  const coord = mergeLayer(sp, "coordinates");
  const fmsa = mergeLayer(sp, "full_msa");
  const cass = mergeLayer(sp, "cassette_msa");
  const syn = mergeLayer(sp, "synteny");
  const inclusion = sp.isoforms.map((i) => `${i.isoform} ${i.interpro_included === "supplement" ? "supplement" : "primary"}`).join(" · ");
  const anyPrimary = sp.isoforms.some((i) => i.interpro_included === "primary");
  const reviewNotes = sp.isoforms.filter((i) => i.review_explanation);

  const overall = sp.overall_readiness;
  const overallNote = overall === "review"
    ? "At least one isoform is kept as an optional supplement/review row."
    : overall === "minor"
      ? "Accepted into the primary InterPro set (with minor annotation flags)."
      : "Accepted into the primary InterPro set.";

  const isoCard = (iso, i) => i
    ? <SummaryCard label={`${iso} status`} cls={i.readiness_class} value={readinessLabel(i.readiness_class)} note={i.layers?.readiness?.note} />
    : <SummaryCard key={iso} label={`${iso} status`} cls="neutral" value="no row" note={`No ${iso} row for this species.`} />;

  return (
    <div className="summary-tab">
      <div className="summary-cards">
        <SummaryCard label="Final dataset status" cls={overall} value={readinessLabel(overall)} note={overallNote} />
        {isoCard("IIIb", iiib)}
        {isoCard("IIIc", iiic)}
        {coord && <SummaryCard label="Coordinate validation" cls={coord.class} value={coord.value} note={coord.note} />}
        {fmsa && <SummaryCard label="MSA support" cls={fmsa.class} value={fmsa.value} note={fmsa.note} />}
        {cass && <SummaryCard label="Cassette / residue agreement" cls={cass.class} value={cass.value} note={cass.note} />}
        {syn && <SummaryCard label="Synteny support" cls={syn.class} value={syn.value} note={syn.note} />}
        <SummaryCard label="InterPro input" cls={anyPrimary ? "accepted" : "review"} value={inclusion}
          note={anyPrimary ? "Included in the primary InterPro FASTA." : "Kept in the review-included (optional) FASTA only."} />
      </div>

      {reviewNotes.length > 0 && (
        <div className="review-notes">
          {reviewNotes.map((i) => (
            <div key={i.isoform} className="review-note warn">
              <b>{i.isoform}:</b> {i.review_explanation}
            </div>
          ))}
        </div>
      )}
    </div>
  );
}

function IsoformsTab({ sp, model, eventType }) {
  if (eventType !== "validated") {
    const isoforms = (sp.proteins || sp.isoforms || []).map(normalizeProtein);
    const candidates = model?.candidate_evidence?.candidates || [];
    // Real per-model exon counts come from the coordinate-track models (coding exons).
    const coordModels = model?.shared_indices?.coordinate_track_index?.models
      || model?.shared_event_indices?.coordinates?.models || [];
    const exonCountFor = (pid) => {
      const m = coordModels.find((x) => x.protein_id === pid);
      return m?.blocks?.length ?? null;
    };
    const primaryLen = isoforms.find((x) => x.is_primary)?.protein_length;
    const enriched = isoforms.map((iso) => {
      const affected = candidates.filter((c) =>
        (c.protein_isoform_evidence?.affected_proteins || []).includes(iso.protein_id));
      const len = iso.protein_length ?? iso.length_aa;
      const lenDiff = (!iso.is_primary && primaryLen && len)
        ? (len === primaryLen ? "same length as primary"
          : `${len > primaryLen ? "+" : ""}${len - primaryLen} aa vs primary`)
        : "";
      return {
        ...iso,
        exon_count: iso.exon_count ?? exonCountFor(iso.protein_id),
        selection_reason: iso.is_primary
          ? (iso.selection_reason || model?.primary_selection_evidence?.explanation || "")
          : (iso.selection_reason || ""),
        completeness: iso.completeness || (iso.curation_status === "curated"
          ? "curated RefSeq (complete model)" : "predicted model"),
        diff_from_primary: iso.is_primary ? "— (this is the primary model)"
          : [lenDiff, affected.length ? `${affected.length} candidate region(s)` : ""].filter(Boolean).join(" · ") || "—",
        affected_candidates: affected.map((c) => c.rank_label || c.candidate_id).join(", ") || "—",
      };
    });
    return (
      <div>
        <ProteinIsoformTable isoforms={enriched}
          selectionMethod={model?.primary_selection_evidence?.selection_rule_label
            || model?.analysis?.selection_method} expandable />
      </div>
    );
  }
  return (
    <div>
      <div className="iso-grid">
        {["IIIb", "IIIc"].map((iso) => {
          const i = sp.isoforms.find((x) => x.isoform === iso);
          if (!i) return <div key={iso} className="iso-col muted">No {iso} row</div>;
          const inc = i.interpro_included;
          return (
            <div key={iso} className={`iso-col iso-frame-${iso.toLowerCase()}`}>
              <div className="iso-col-head">
                <IsoBadge iso={iso} />
                <Badge cls={i.readiness_class}>{readinessLabel(i.readiness_class)}</Badge>
              </div>
              <Field label="Final isoform label">{i.final_isoform_label || "—"}</Field>
              <Field label="Status">
                {(() => {
                  const st = mainDisplayStatus(i.layers?.readiness?.value || i.final_claim_status_after_rescue);
                  return <Badge cls={st.cls} soft>{st.label}</Badge>;
                })()}
              </Field>
              <Field label="InterPro input included?">
                <Badge cls={inc === "primary" ? "accepted" : inc === "supplement" ? "review" : "neutral"} soft>
                  {inc === "primary" ? "yes · primary (main)" : inc === "supplement" ? "supplement (optional)" : "—"}
                </Badge>
              </Field>
              <Field label="Transcript ID"><code>{i.transcript_id || "—"}</code></Field>
              <Field label="Protein ID"><code>{i.protein_id || "—"}</code></Field>
              <Field label="Protein length">{i.protein_length ? `${i.protein_length} aa` : "—"}</Field>
              <Field label="Coordinate status">
                <Badge cls={i.layers?.coordinates?.class} soft>{i.layers?.coordinates?.value || "—"}</Badge>
              </Field>
              <Field label="MSA status"><Badge cls={i.layers?.full_msa?.class} soft>{i.layers?.full_msa?.value || "—"}</Badge></Field>
              {i.review_explanation && i.readiness_class === "review" && (
                <div className="review-note warn">{i.review_explanation}</div>
              )}
            </div>
          );
        })}
      </div>
    </div>
  );
}

// Calm main-view label per status colour family. Used only when a layer carries
// a provenance/offset tone, so that reconciled/rescued/native-offset layers read
// as a plain status instead of surfacing the technical qualifier in the main grid.
const CLASS_MAIN_LABEL = {
  accepted: "Accepted",
  minor: "Supported with minor note",
  inspection: "Inspection note",
  review: "Supplement / review",
  excluded: "Failed",
  neutral: "—",
};

// Validated (FGFR2) evidence layers only — the generic Evidence tab was removed.
function EvidenceTab({ sp }) {
  const layers = [
    ["Final isoform assignment", "label"],
    ["Rescue / provenance", "rescue"],
    ["Orthology / paralog", "orthology"],
    ["Coordinate mapping", "coordinates"],
    ["Full-length MSA", "full_msa"],
    ["Cassette MSA / residue agreement", "cassette_msa"],
    ["Synteny / locus", "synteny"],
    ["InterPro input", "readiness"],
  ];
  return (
    <div className="evidence-list">
      {sp.isoforms.map((i) => (
        <div key={i.isoform} className="evidence-iso">
          <h4><span className={`iso iso-${i.isoform.toLowerCase()}`}>{i.isoform}</span></h4>
          <div className="evidence-rows">
            {layers.map(([label, key]) => {
              const c = i.layers?.[key];
              // If a layer only carries a provenance/offset tone, show the calm
              // status label rather than the technical value ("native offset",
              // "corrected", …). No orange tone markers in the main grid.
              const shown = c?.tone ? (CLASS_MAIN_LABEL[c.class] || c.value) : (c?.value || "—");
              return (
                <div key={key} className="evidence-row">
                  <span className="er-label">{label}</span>
                  <span className="er-status">
                    <Badge cls={c?.class} soft>{shown}</Badge>
                  </span>
                </div>
              );
            })}
          </div>
        </div>
      ))}
    </div>
  );
}
