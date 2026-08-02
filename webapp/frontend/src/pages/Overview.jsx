import { useMemo, useState } from "react";
import { fileUrl } from "../api";
import { Badge, Drawer, Field, Spinner, Modal } from "../ui";
import { DatasetPageHeader, KpiGrid } from "../components/shared";
import CassetteExplorer from "./viewers/CassetteExplorer";
import CoordinateTrack from "./viewers/CoordinateTrack";
import MsaExplorer from "./viewers/MsaExplorer";
import SyntenyViewer from "./viewers/SyntenyViewer";
import SpeciesStory from "./viewers/SpeciesStory";
import HumanReferenceBadge from "./HumanReferenceBadge";

const VIEWER_TITLES = {
  cassette: "Cassette Sequence Explorer",
  coordinates: "Exon → Protein Coordinate Track",
  msa: "MSA Explorer",
  synteny: "Synteny Locus Viewer",
  story: "Evidence Story",
};

// which interactive viewer best matches each evidence layer
const LAYER_VIEWER = {
  label: "story", rescue: "story", readiness: "story",
  coordinates: "coordinates", full_msa: "msa",
  cassette_msa: "cassette", synteny: "synteny", orthology: "synteny",
};

export default function Overview({ model, datasetInfo, isExample, setPage, openGene }) {
  return <UnifiedOverview model={model} datasetInfo={datasetInfo} isExample={isExample}
    setPage={setPage} openGene={openGene} />;
}

// Shown only when the assembly calls the gene something else than the user did. The
// requested symbol stays the headline; this explains which locus it reached and how, so
// the result is checkable against the source without the LOC id replacing the gene name.
const RESOLUTION_ROUTE = {
  annotation_exact_symbol: "matched the annotation's own gene symbol",
  annotation_case_normalized_symbol: "matched the annotation's gene symbol",
  annotation_provided_alias: "matched an alias in the annotation",
  ncbi_gene_official_symbol: "matched the official NCBI Gene symbol",
  ncbi_alias_to_geneid_to_annotation:
    "resolved through NCBI Gene as an alias, then mapped into the assembly by GeneID",
  ncbi_alias_to_official_symbol:
    "resolved through NCBI Gene as an alias, then matched by its official symbol",
};

// Concise reasons for the two analyses that are commonly not applicable, so Summary can
// state them in a few words instead of the full sentence the analysis page shows.
const SCOPE_REASON = {
  single_coding_exon: "No internal coding-exon boundaries",
  no_internal_coding_exon_boundaries: "No internal coding-exon boundaries",
  single_unique_protein_sequence: "One distinct protein sequence",
  single_species_dataset: "Single-species dataset",
};

/**
 * Which analyses do not apply to this gene, and why, in one line each.
 *
 * Only the settled not-applicable states are listed. Analyses that are available need no
 * explanation, and anything still pending or missing is a run-status matter reported by the
 * status stepper, not an "analysis scope" question.
 */
function AnalysisScopeCard({ availability }) {
  const notApplicable = (availability?.analyses || [])
    .filter((a) => a.status === "not_applicable");
  if (!notApplicable.length) return null;
  return (
    <div className="card analysis-scope">
      <h3>Analysis scope</h3>
      <p className="muted small">
        These analyses do not apply to the recovered gene model. The run is complete
        without them.
      </p>
      <dl className="scope-list">
        {notApplicable.map((a) => (
          <div key={a.analysis_name} className="scope-row">
            <dt>{a.label}</dt>
            <dd>
              <span className="scope-value">Not applicable</span>
              <span className="muted small">
                {SCOPE_REASON[a.reason_code]
                  || `${a.prerequisite_name}: ${a.prerequisite_count ?? 0}`}
              </span>
            </dd>
          </div>
        ))}
      </dl>
    </div>
  );
}

function GeneResolutionNote({ identity }) {
  if (!identity?.any_symbol_differs_from_source) return null;
  const bySpecies = identity.by_species || {};
  const entries = Object.entries(identity.source_symbols_by_species || {});
  if (!entries.length) return null;
  return (
    <div className="card gene-resolution">
      <h3>Gene identity</h3>
      <p className="muted small">
        Requested gene <strong>{identity.requested_gene_symbol}</strong>. The assembly
        annotates it under a different symbol, shown per species below.
      </p>
      <ul className="gene-resolution-list">
        {entries.map(([species, symbol]) => {
          const rec = bySpecies[species] || {};
          const route = RESOLUTION_ROUTE[rec.resolution_method] || rec.resolution_method;
          return (
            <li key={species}>
              <span className="mono">{species.replace(/_/g, " ")}</span>
              {" — NCBI annotation symbol "}
              <strong className="mono">{symbol}</strong>
              {rec.resolved_gene_id && <> (GeneID {rec.resolved_gene_id})</>}
              {rec.source_description && <>: {rec.source_description}</>}
              {route && <div className="muted small">{route}</div>}
            </li>
          );
        })}
      </ul>
    </div>
  );
}

function UnifiedOverview({ model, datasetInfo, isExample, setPage, openGene }) {
  const validated = model?.event_layer?.type === "validated";
  const vi = model?.validated_event_indices || {};
  const summary = model?.overview?.summary || model?.overview?.run_index || model?.overview || {};
  const stack = normalizeEvidenceStack(model);
  const arch = model?.domain_architecture?.summary
    || model?.legacy_fgfr2_indices?.domain_architecture_summary || null;
  const bc = model?.boundary?.validated?.summary || model?.boundary?.summary
    || model?.validated_event_indices?.boundary_consistency_summary || null;
  const [filterGroup, setFilterGroup] = useState("all");
  const [filterClaim, setFilterClaim] = useState("all");
  const [search, setSearch] = useState("");
  const [cell, setCell] = useState(null);
  const [viewer, setViewer] = useState(null);   // {kind, species, isoform}
  const k = summary.kpi || summary.kpis || {};
  const validatedKpis = [
    ["Species", k.species, "ortholog panel", ""],
    ["Isoform rows", k.isoform_rows, "IIIb + IIIc", ""],
    ["Primary-ready", k.primary_ready, "InterPro-ready", "accepted"],
    ["Review-only", k.review_only, "supplement", k.review_only ? "review" : "accepted"],
    ["Rescued & validated", k.rescued_validated, "accepted as primary", "accepted"],
    ["Freeze", summary.freeze_ready ? "ready" : "—", `${k.primary_fasta_sequences} sequences`, summary.freeze_ready ? "accepted" : "neutral"],
  ];
  const sharedKpis = [
    ["Species", k.species_analysed ?? model?.species?.length ?? "—", "ortholog panel", ""],
    ["Isoform rows", k.protein_isoforms ?? "—", "protein isoforms", ""],
    ["Primary-ready", k.primary_proteins ?? "—", "InterPro-ready", "accepted"],
    ["Review-only", k.review_only ?? 0, "supplement", "accepted"],
    ["Synteny evidence", k.synteny_neighbours ?? "—", "neighbours resolved", k.synteny_neighbours ? "accepted" : "neutral"],
    ["Event layer", model?.event_layer?.type === "exploratory" ? "exploratory" : "none",
      `${k.exploratory_event_candidates ?? 0} candidates`,
      model?.event_layer?.type === "exploratory" ? "minor" : "neutral"],
  ];
  const kpis = validated ? validatedKpis : sharedKpis;

  const groups = useMemo(() => {
    if (!stack) return [];
    return Array.from(new Set(stack.rows.map((r) => r.taxon_group))).sort();
  }, [stack]);

  const rows = useMemo(() => {
    if (!stack) return [];
    return stack.rows.filter((r) => {
      if (filterGroup !== "all" && r.taxon_group !== filterGroup) return false;
      if (filterClaim === "primary" && r.row_class === "review") return false;
      if (filterClaim === "review" && r.row_class !== "review") return false;
      if (search && !`${r.display_species_name} ${r.species}`.toLowerCase().includes(search.toLowerCase())) return false;
      return true;
    });
  }, [stack, filterGroup, filterClaim, search]);

  return (
    <section className="page">
      <DatasetPageHeader
        eyebrow="Mission control"
        title="Overview"
        badges={validated ? <>
            {!isExample && <HumanReferenceBadge humanReference={datasetInfo?.human_reference} />}
            <Badge cls={summary.gate_status === "pass" ? "accepted" : "excluded"} soft>
              Gate {summary.gate_status} · {summary.gate_checks_total - summary.gate_checks_failed}/{summary.gate_checks_total}
            </Badge>
            <Badge cls={summary.full_clean_run_completed ? "accepted" : "info"} soft>
              {summary.full_clean_run_completed
                ? "Full clean run"
                : (isExample ? "Freeze-based" : "Fresh run-local analysis")}
            </Badge>
          </> : <>
            <Badge cls="accepted" soft>Gene models ready</Badge>
            <Badge cls={model?.domain_architecture?.available ? "accepted" : "neutral"} soft>
              {model?.domain_architecture?.available ? "Domain annotation complete" : "Domain annotation pending"}
            </Badge>
          </>}
      />

      <KpiGrid items={kpis} />

      <GeneResolutionNote identity={summary.gene_identity
        || model?.gene_identity || datasetInfo?.gene_identity} />

      <AnalysisScopeCard availability={model?.analysis_availability} />


      {setPage && (
        <div className="card run-cta">
          <div>
            <h3>Create or open run</h3>
            <p className="muted small">
              {validated
                ? "Start a new local FGFR2 run or open an existing one, and copy the local pipeline commands."
                : "Start a new local run or open an existing one, and copy the local pipeline commands."}
            </p>
          </div>
          <button className="btn primary sm" onClick={() => setPage("runs")}>Run Workflow →</button>
        </div>
      )}

      {arch?.available && <ArchitectureCard arch={arch} setPage={setPage} />}

      {validated && bc?.available && <BoundaryCard bc={bc} setPage={setPage} />}

      <div className="card">
        <div className="card-head">
          <h3>Evidence stack</h3>
          <div className="filters">
            <input className="search" placeholder="Search species…" value={search} onChange={(e) => setSearch(e.target.value)} />
            <select value={filterGroup} onChange={(e) => setFilterGroup(e.target.value)}>
              <option value="all">All taxa</option>
              {groups.map((g) => <option key={g} value={g}>{g}</option>)}
            </select>
            <select value={filterClaim} onChange={(e) => setFilterClaim(e.target.value)}>
              <option value="all">All rows</option>
              <option value="primary">Primary only</option>
              <option value="review">Review only</option>
            </select>
          </div>
        </div>

        {!stack ? <Spinner /> : (
          <Heatmap stack={stack} rows={rows} onCell={setCell} />
        )}
        <Legend />
      </div>

      <Drawer
        open={Boolean(cell)}
        onClose={() => setCell(null)}
        title={cell ? `${cell.row.display_species_name} · ${cell.row.isoform}` : ""}
        subtitle={cell ? cell.col.label : ""}
      >
        {cell && (
          <>
            <div className="drawer-badges">
              <Badge cls={cell.data.class}>{cell.data.value || "—"}</Badge>
            </div>
            {cell.data.note && <p className="drawer-lead">{cell.data.note}</p>}

            <Field label="Final isoform assignment">{cell.row.final_isoform_label || "—"}</Field>
            <Field label="Validated exon type">{cell.row.validated_exon_type || "—"}</Field>
            <Field label="InterPro input">
              <Badge cls={cell.row.interpro_included === "primary" ? "accepted" : cell.row.interpro_included === "supplement" ? "review" : "neutral"} soft>
                {cell.row.interpro_included === "primary" ? "Primary (main)" : cell.row.interpro_included === "supplement" ? "Review-included (optional)" : "—"}
              </Badge>
            </Field>
            <Field label="Taxon group">{cell.row.taxon_group}</Field>
            {cell.data.source_table && (
              <Field label="Source table" wide>
                <a className="link" href={fileUrl(cell.data.source_table)}>{cell.data.source_table.split("/").pop()}</a>
              </Field>
            )}

            {cell.col.key === "label" && (
              <p className="drawer-lead sm">
                {cell.row.upstream_corrected
                  ? "Upstream annotation was reconciled with sequence evidence; the final label above is the accepted result."
                  : cell.row.upstream_consistent
                    ? "Upstream annotation was already consistent with the sequence-calibrated label."
                    : "Final label was assigned from sequence evidence."}
              </p>
            )}

            {validated && <div className="drawer-actions">
              <span className="field-label">Open interactive evidence</span>
              <div className="action-chips">
                {[["cassette", "Cassette"], ["coordinates", "Coordinates"], ["msa", "MSA"], ["synteny", "Synteny"], ["story", "Story"]].map(([k, l]) => (
                  <button key={k} className="btn ghost sm" onClick={() => setViewer({ kind: k, species: cell.row.species, isoform: cell.row.isoform })}>{l}</button>
                ))}
              </div>
            </div>}
            <button className="btn ghost" onClick={() => { const kind = LAYER_VIEWER[cell.col.key] || "story"; setCell(null); openGene({ species: cell.row.species, isoform: cell.row.isoform, tab: kind === "story" ? "evidence" : kind }); }}>
              Open in Gene Explorer →
            </button>
          </>
        )}
      </Drawer>

      <Modal open={Boolean(viewer)} onClose={() => setViewer(null)}
        title={viewer ? VIEWER_TITLES[viewer.kind] : ""}
        subtitle={viewer ? `${viewer.species.replaceAll("_", " ")} · ${viewer.isoform}` : ""}>
        {viewer?.kind === "cassette" && <CassetteExplorer preloaded={vi.cassette || vi.cassette_residue_index || {}} species={viewer.species} />}
        {viewer?.kind === "coordinates" && <CoordinateTrack preloaded={vi.coordinates || vi.coordinate_track_index || {}} species={viewer.species} />}
        {viewer?.kind === "msa" && <MsaExplorer preloaded={vi.msa || vi.msa_index || {}} species={viewer.species} />}
        {viewer?.kind === "synteny" && <SyntenyViewer preloaded={vi.synteny || model?.synteny || {}} species={viewer.species} />}
        {viewer?.kind === "story" && <SpeciesStory preloaded={vi.story || vi.species_story_index || {}} species={viewer.species} isoform={viewer.isoform} />}
      </Modal>
    </section>
  );
}

function normalizeEvidenceStack(model) {
  const existing = model?.evidence?.stack || model?.overview?.evidence_stack;
  if (existing?.columns && existing?.rows) return existing;
  if (model?.evidence?.columns && model?.evidence?.rows) return model.evidence;
  const items = model?.evidence?.items || existing?.items || model?.overview?.evidence_summary || [];
  if (!items.length) return null;
  const columns = items.map((it) => ({ key: it.id, label: it.title || it.label || it.id }));
  const species = model?.species || [];
  const rows = species.map((sp) => {
    const local = sp.evidence || items;
    const byId = Object.fromEntries(local.map((it) => [it.id, it]));
    return {
      species: sp.species || sp.species_id,
      display_species_name: sp.display_species_name || sp.display_name || sp.species,
      taxon_group: sp.taxon_group || "Analysed species",
      isoform: "gene",
      row_class: "primary",
      final_isoform_label: "",
      validated_exon_type: "",
      interpro_included: "primary",
      provenance: {},
      cells: Object.fromEntries(items.map((it) => {
        const value = byId[it.id] || it;
        return [it.id, {
          class: evidenceCellClass(value.status),
          value: value.status || "—",
          raw: value.status || "",
          note: value.explanation || value.note || "",
          source_table: value.source_file || "",
        }];
      })),
    };
  });
  return { columns, rows };
}

function evidenceCellClass(status) {
  const value = String(status || "").toLowerCase();
  if (["accepted", "available", "complete", "found"].includes(value)) return "accepted";
  if (["exploratory", "pending", "pending_cluster", "medium"].includes(value)) return "minor";
  if (["failed", "excluded", "not_found"].includes(value)) return "excluded";
  return "neutral";
}

function ArchitectureCard({ arch, setPage }) {
  const stats = [
    ["Proteins annotated", arch.proteins_annotated, ""],
    ["Kinase detected", arch.with_kinase, "accepted"],
    ["TM predicted (pyTMHMM)", arch.with_tm_pytmhmm, "accepted"],
    ["Supported architectures", arch.supported, "accepted"],
    ["Minor flags", arch.minor_flags, arch.minor_flags ? "minor" : "accepted"],
    ["Review warnings", arch.review_warnings, arch.review_warnings ? "review" : "accepted"],
  ];
  return (
    <div className="card">
      <div className="card-head">
        <h3>Post-InterPro architecture</h3>
        <div className="head-badges">
          <Badge cls="info" soft>TM layer: pyTMHMM</Badge>
          {setPage && <button className="btn ghost sm" onClick={() => setPage("figures")}>Figure gallery →</button>}
        </div>
      </div>
      <div className="arch-stat-grid">
        {stats.map(([l, v, c]) => (
          <div key={l} className="arch-stat">
            <span className="arch-stat-value"><Badge cls={c || "neutral"} soft>{v ?? "—"}</Badge></span>
            <span className="arch-stat-label">{l}</span>
          </div>
        ))}
      </div>
      {arch.review_warnings > 0 && (
        <div className={`arch-note ${arch.review_all_coordinate_artifacts ? "info" : "warn"}`}>
          <b>{arch.review_warnings} domain-order warning{arch.review_warnings > 1 ? "s" : ""}</b>
          {arch.review_all_coordinate_artifacts
            ? " — all audited as low-confidence exon-block coordinate artifacts (native codon-phase offset), not biological anomalies. IIIb/IIIc labels and membership unchanged."
            : " — see the QC audit for details."}
          {arch.review_cases?.length > 0 && (
            <ul className="arch-review-list">
              {arch.review_cases.map((c) => (
                <li key={`${c.species}-${c.isoform}`}>
                  <b>{c.display_species_name} {c.isoform}</b>: {c.likely_issue.replaceAll("_", " ")}
                </li>
              ))}
            </ul>
          )}
        </div>
      )}
      {arch.report && (
        <div className="arch-links">
          <a className="btn ghost sm" href={fileUrl(arch.report)}>Architecture summary report</a>
        </div>
      )}
    </div>
  );
}

function BoundaryCard({ bc, setPage }) {
  const c = bc.boundary_class_counts || {};
  const total = bc.total_primary_proteins ?? 0;
  const stats = [
    ["Primary proteins", total, ""],
    ["Missing data", (bc.proteins_with_cassette_data === total && bc.proteins_with_interpro_domain_data === total && bc.proteins_with_tm_data === total) ? 0 : "—", "accepted"],
    ["Between-domain", c.between_domains ?? 0, (c.between_domains ? "minor" : "accepted")],
    ["Aligned", c.aligned_to_domain_boundary ?? 0, "accepted"],
    ["Near boundary", c.near_domain_boundary ?? 0, "accepted"],
    ["Within domain", c.within_domain ?? 0, "neutral"],
  ];
  return (
    <div className="card">
      <div className="card-head">
        <h3>Boundary consistency</h3>
        <div className="head-badges">
          <Badge cls="info" soft>Cassette ↔ domain boundary</Badge>
          {setPage && <button className="btn primary sm" onClick={() => setPage("boundary")}>Open explorer →</button>}
        </div>
      </div>
      <div className="arch-stat-grid">
        {stats.map(([l, v, cls]) => (
          <div key={l} className="arch-stat">
            <span className="arch-stat-value"><Badge cls={cls || "neutral"} soft>{v ?? "—"}</Badge></span>
            <span className="arch-stat-label">{l}</span>
          </div>
        ))}
      </div>
      <div className="arch-note info">
        Cassette-end boundaries sit near the Ig-like domain boundary; cassette-start boundaries
        usually lie within the Ig-like region — consistent with robust exon–domain boundary
        identification across the primary vertebrate FGFR2 dataset.
      </div>
    </div>
  );
}

function Heatmap({ stack, rows, onCell }) {
  return (
    <div className="heatmap-wrap">
      <table className="heatmap">
        <thead>
          <tr>
            <th className="sticky-col">Species</th>
            <th>Iso</th>
            {stack.columns.map((c) => <th key={c.key} title={c.label}>{c.label}</th>)}
          </tr>
        </thead>
        <tbody>
          {rows.map((r) => (
            <tr key={`${r.species}-${r.isoform}`} className={r.row_class === "review" ? "row-review" : ""}>
              <td className="sticky-col species-name">
                {r.display_species_name}
                <small>{r.taxon_group}</small>
              </td>
              <td><span className={`iso iso-${r.isoform.toLowerCase()}`}>{r.isoform}</span></td>
              {stack.columns.map((c) => {
                const d = r.cells[c.key];
                return (
                  <td key={c.key} className="hm-cell">
                    <button
                      className={`cell st-${d.class}${d.tone ? " tone-" + d.tone : ""}`}
                      title={`${c.label}: ${d.value}${d.note ? " — " + d.note : ""}`}
                      onClick={() => onCell({ row: r, col: c, data: d })}
                    >
                      <span className="cell-dot" />
                    </button>
                  </td>
                );
              })}
            </tr>
          ))}
        </tbody>
      </table>
      {rows.length === 0 && <p className="muted pad">No rows match the current filters.</p>}
    </div>
  );
}

function Legend() {
  const items = [
    ["accepted", "Accepted / supported"],
    ["minor", "Supported · minor note"],
    ["review", "Supplement / review"],
    ["excluded", "Failed"],
    ["neutral", "Not applicable"],
  ];
  return (
    <div className="legend legend-stack">
      <div className="legend">
        {items.map(([c, l]) => (
          <span key={c} className="legend-item"><span className={`cell-dot st-${c}`} />{l}</span>
        ))}
      </div>
    </div>
  );
}
