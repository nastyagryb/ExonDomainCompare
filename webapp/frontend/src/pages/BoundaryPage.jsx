import { useMemo, useState } from "react";
import { fileUrl } from "../api";
import { Badge, Empty, Modal } from "../ui";
import { DatasetPageHeader, KpiGrid, PendingAnalysisCard } from "../components/shared";
import BoundarySummaryCards from "./viewers/BoundarySummaryCards";
import BoundaryHeatmap from "./viewers/BoundaryHeatmap";
import BoundaryDistancePlot from "./viewers/BoundaryDistancePlot";
import BoundaryOutlierTable from "./viewers/BoundaryOutlierTable";
import BoundaryDetailTrack from "./viewers/BoundaryDetailTrack";
import SignedDistancePlot from "./viewers/SignedDistancePlot";
import GlobalBoundaryDashboard from "./GlobalBoundaryDashboard";
import { canonClass, CANON_CLASS_COLOR, CANON_CLASS_LABEL } from "./viewers/boundaryClasses";

// The single boundary controller for every dataset. The model decides whether
// the shared page receives validated event-boundary results (FGFR2 Boundary
// Consistency — unchanged), the redesigned generic global dashboard (driven by
// the validated coordinate model + boundary_dashboard contract), or a pending
// stage gate.
export default function BoundaryPage({ model, setPage, openGene, labels }) {
  const eventType = model?.event_layer?.type || "none";

  // Generic (non-event) datasets with a validated coordinate model render the
  // new global Exon–Domain-Boundaries dashboard (single-species / pending /
  // multi-species). Page mode is resolved from the canonical dataset — no
  // gene-symbol branching. FGFR2's validated event page below is untouched.
  // This is evaluated before any hooks so hook order stays stable per mount
  // (the key on <BoundaryPage> remounts it when the active dataset changes).
  if (eventType !== "validated" && model?.protein_coordinate_model?.models?.length) {
    return <GlobalBoundaryDashboard model={model} setPage={setPage} openGene={openGene} />;
  }
  return <ValidatedOrPendingBoundaryPage model={model} setPage={setPage} labels={labels} />;
}

function ValidatedOrPendingBoundaryPage({ model, setPage, labels }) {
  const [sel, setSel] = useState(null);
  const eventType = model?.event_layer?.type || "none";
  const boundary = model?.boundary || {};
  const validated = boundary.validated || model?.validated_event_indices?.boundary || {};
  const summary = validated.summary || boundary.summary
    || model?.validated_event_indices?.boundary_consistency_summary || null;
  const matrix = validated.matrix || boundary.matrix
    || model?.validated_event_indices?.boundary_consistency_matrix || null;
  const outliers = validated.outliers || boundary.outliers
    || model?.validated_event_indices?.boundary_consistency_outliers || null;
  const index = validated.index || boundary.index
    || model?.validated_event_indices?.boundary_consistency_index || boundary || null;
  const arch = validated.architecture || model?.domain_architecture?.species_index
    || model?.validated_event_indices?.domain_architecture_species
    || model?.legacy_fgfr2_indices?.species_domain_architecture || null;

  const selPanel = useMemo(() => {
    if (!sel || !arch?.species) return null;
    const sp = arch.species.find((s) => s.species === sel.species);
    return sp?.panels?.[sel.isoform] || null;
  }, [sel, arch]);
  const selRow = useMemo(() => {
    if (!sel || !matrix?.rows) return null;
    return matrix.rows.find((r) => r.species === sel.species && r.isoform === sel.isoform) || null;
  }, [sel, matrix]);

  if (eventType !== "validated") {
    const generic = boundary.generic || boundary;
    if (generic.available) {
      return <GenericBoundaryAnalysis generic={generic}
        domains={model?.domain_architecture || model?.domain_features}
        gene={generic.gene_symbol || model?.analysis?.gene_symbol || "gene"} />;
    }
    const command = model?.stages?.cluster_command || model?.dataset?.cluster_command;
    return (
      <section className="page">
        <DatasetPageHeader eyebrow={`Mission control · ${model?.analysis?.gene_symbol || "gene"}`}
          title="Exon–Domain Analysis"
          badges={<Badge cls="neutral" soft>pending domain annotation</Badge>} />
        <div className="arch-note info">
          <b>All-exon boundary mode.</b> Domain boundary distances are not computed yet.
          They become available after InterProScan / pyTMHMM; no result is inferred before then.
        </div>
        <PendingAnalysisCard title="Analysis inputs"
          badge={<Badge cls="neutral" soft>pending cluster</Badge>}
          description="Run the cluster round-trip from a local terminal. This same page will render the resulting all-exon boundary analysis."
          command={command} />
      </section>
    );
  }

  if (!summary?.available) {
    return (
      <section className="page">
        <Empty title="Boundary consistency not available"
          hint="The validated event-boundary analysis was not found for this dataset." />
      </section>
    );
  }

  const figs = index?.source_figures || summary.figure_links || {};
  return (
    <section className="page bc-page">
      <div className="page-head">
        <div>
          <p className="eyebrow">Final thesis analysis · Module 1</p>
          <h2>Boundary Consistency Explorer</h2>
        </div>
        <div className="head-badges">
          <Badge cls="info" soft>Cassette ↔ domain boundary</Badge>
          {setPage && <button className="btn ghost sm" onClick={() => setPage("figures")}>Figure gallery →</button>}
        </div>
      </div>
      <BoundarySummaryCards summary={summary} />
      <div className="card">
        <div className="card-head">
          <h3>{labels?.domainRelationDescription || "Cassette-to-domain boundary consistency"}</h3>
          <span className="muted sm">Rows ordered taxonomically · click a cell for the protein-level detail track</span>
        </div>
        <BoundaryHeatmap matrix={matrix} onSelect={(r) => setSel({ species: r.species, isoform: r.isoform })}
          selected={sel} />
      </div>
      <div className="card">
        <div className="card-head">
          <h3>Distance to nearest protein-domain boundary</h3>
          <span className="muted sm">Interactive distribution · hover a point for species detail</span>
        </div>
        <BoundaryDistancePlot matrix={matrix} />
      </div>
      <div className="card">
        <BoundaryOutlierTable outliers={outliers} onSelect={(t) => t && setSel({ species: t.species, isoform: t.isoform })} />
      </div>
      <div className="card bc-static-figs">
        <div className="card-head">
          <h3>Static thesis figures</h3>
          <Badge cls="neutral" soft>downloadable previews</Badge>
        </div>
        <p className="muted sm">The interactive explorer above is the primary view. These are the publication-ready static figures.</p>
        <div className="arch-links">
          {["heatmap", "distance_distribution"].map((k) => {
            const f = figs[k];
            if (!f) return null;
            const label = k === "heatmap" ? "Figure 11 · heatmap" : "Figure 12 · distance distribution";
            return (
              <span key={k} className="bc-fig-links">
                <span className="bc-fig-name">{label}:</span>
                {f.png && <a className="btn ghost sm" href={fileUrl(f.png, true)} target="_blank" rel="noreferrer">PNG</a>}
                {f.svg && <a className="btn ghost sm" href={fileUrl(f.svg)}>SVG</a>}
                {f.pdf && <a className="btn ghost sm" href={fileUrl(f.pdf)}>PDF</a>}
              </span>
            );
          })}
          {index?.report && <a className="btn ghost sm" href={fileUrl(index.report)}>Analysis report (MD)</a>}
        </div>
      </div>
      <Modal open={Boolean(sel)} onClose={() => setSel(null)}
        title={selRow ? selRow.display_label : (sel ? `${sel.species.replaceAll("_", " ")} · ${sel.isoform}` : "")}
        subtitle="Protein-level boundary detail · cassette markers vs. nearest domain boundary">
        {sel && <BoundaryDetailTrack pd={selPanel} row={selRow} isoform={sel.isoform} />}
      </Modal>
    </section>
  );
}

function GenericBoundaryAnalysis({ generic, domains, gene }) {
  const [selBoundary, setSelBoundary] = useState(null);
  const bySpecies = useMemo(() => {
    const map = new Map();
    for (const p of (generic.proteins || [])) {
      const sid = p.species_id || p.species || "—";
      if (!map.has(sid)) map.set(sid, []);
      for (const b of (p.boundaries || [])) {
        map.get(sid).push({ ...b, protein_id: p.protein_id, species_id: sid });
      }
    }
    return map;
  }, [generic]);
  const rows = [...bySpecies.values()].flat();
  // Canonical (5-class) counts derived from real rows; falls back to index counts.
  const counts = useMemo(() => {
    const c = { exact_domain_edge: 0, near_domain_edge: 0, inside_domain: 0,
      outside_annotated_domains: 0, unavailable_or_uncertain: 0 };
    for (const r of rows) c[canonClass(r.category || r.class)] += 1;
    return c;
  }, [rows]);
  const threshold = generic.near_edge_threshold_aa ?? 5;
  const speciesScope = generic.species_scope || [...bySpecies.keys()];
  const displaySpecies = (s) => (s || "").replaceAll("_", " ");
  const selId = selBoundary && (selBoundary.exon_boundary_id || selBoundary.id);

  return (
    <section className="page">
      <DatasetPageHeader eyebrow={`Gene analysis · ${gene}`} title="Exon–Domain Boundaries"
        badges={<Badge cls="accepted" soft>real InterProScan coordinates</Badge>} />
      {/* Explicit scope: never present one aggregate number without saying what it covers. */}
      <div className="arch-note info">
        <b>Scope.</b> Gene <b>{gene}</b> · {speciesScope.length} species
        (<i>{speciesScope.map(displaySpecies).join(", ") || "—"}</i>) ·
        {" "}protein scope <b>{generic.protein_scope || "primary_only"}</b>
        {" "}({generic.n_proteins ?? new Set(rows.map((r) => r.protein_id)).size} proteins) ·
        {" "}internal coding-exon boundaries · domain layer <b>{generic.domain_layer || "representative_domain"}</b>
        {" "}(representative InterPro structural domains only; families, sites and disorder are excluded) ·
        {" "}near-edge threshold <b>{threshold} aa</b>.
      </div>
      <KpiGrid items={[
        ["Boundaries", generic.n_boundaries ?? rows.length, `across ${generic.n_proteins ?? "—"} primary proteins`, ""],
        ["Exact domain edge", counts.exact_domain_edge, "0 aa", "accepted"],
        ["Near domain edge", counts.near_domain_edge, `≤ ${threshold} aa`, "minor"],
        ["Inside domain", counts.inside_domain, "within a domain", "info"],
        ["Outside annotated", counts.outside_annotated_domains, "beyond annotated domains", "neutral"],
        ["Unavailable", counts.unavailable_or_uncertain, "no domain / uncertain", "neutral"],
      ]} />
      {[...bySpecies.entries()].map(([sid, srows]) => (
        <div className="card" key={sid}>
          <div className="card-head"><h3><i>{displaySpecies(sid)}</i></h3>
            <span className="muted sm">{srows.length} internal boundaries · signed distance to nearest representative-domain edge (0 = edge)</span></div>
          <SignedDistancePlot rows={srows} threshold={threshold}
            selectedId={selId} onSelect={(r) => setSelBoundary(r)} />
        </div>
      ))}
      <div className="card">
        <div className="card-head"><h3>Boundary-class distribution</h3>
          <span className="muted sm">representative-domain layer · mutually exclusive classes</span></div>
        <div className="class-dist">
          {Object.keys(CANON_CLASS_LABEL).map((c) => {
            const total = rows.length || 1;
            const pct = Math.round((counts[c] / total) * 100);
            return (
              <div className="class-dist-row" key={c}>
                <span className="class-dist-label">{CANON_CLASS_LABEL[c]}</span>
                <span className="class-dist-track">
                  <span className="class-dist-fill" style={{ width: `${pct}%`, background: CANON_CLASS_COLOR[c] }} />
                </span>
                <span className="class-dist-num">{counts[c]} · {pct}%</span>
              </div>
            );
          })}
        </div>
      </div>
      <div className="card">
        <div className="card-head"><h3>Per-boundary evidence and provenance</h3>
          <span className="muted sm">Click a point above or a row here — the selection is shared. Species · protein · nearest domain edge · signed distance</span></div>
        <div className="table-scroll"><table className="mini-tbl">
          <thead><tr><th>Species</th><th>Protein</th><th>Boundary</th><th>Position</th><th>Nearest domain</th>
            <th>Edge</th><th>Signed</th><th>|dist|</th><th>Class</th></tr></thead>
          <tbody>{rows.map((row) => {
            const cls = canonClass(row.category || row.class);
            const rid = `${row.protein_id}:${row.exon_boundary_id}`;
            const isSel = selId && (row.exon_boundary_id === selId || row.id === selId);
            return (
            <tr key={`${rid}:table`} className={isSel ? "row-selected" : ""}
              onClick={() => setSelBoundary(row)} style={{ cursor: "pointer" }}>
              <td><i>{displaySpecies(row.species_id)}</i></td>
              <td>{row.protein_id}</td><td>{row.exon_boundary_id}</td>
              <td>{row.boundary_position_aa}</td>
              <td>{row.nearest_domain_name || row.nearest_domain_accession || row.nearest_domain_id || "—"}
                {row.nearest_domain_type ? <span className="muted sm"> · {row.nearest_domain_type}</span> : null}</td>
              <td>{row.nearest_edge || row.domain_edge_type || row.nearest_domain_boundary_type || "—"}</td>
              <td>{row.signed_distance_aa ?? "—"}</td>
              <td>{row.absolute_distance_aa ?? row.distance_aa ?? "—"}</td>
              <td><span className="class-chip" style={{ background: CANON_CLASS_COLOR[cls] }} />
                <span className="muted sm">{CANON_CLASS_LABEL[cls]}</span></td>
            </tr>);
          })}</tbody>
        </table></div>
        <details className="tech-details"><summary>Domain feature provenance</summary>
          <pre>{JSON.stringify(domains, null, 2)}</pre>
        </details>
      </div>
    </section>
  );
}
