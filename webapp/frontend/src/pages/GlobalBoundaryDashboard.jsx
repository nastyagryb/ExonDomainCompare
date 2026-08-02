import { useMemo, useRef, useState } from "react";
import { Badge, Menu } from "../ui";
import { DATASET_STATUS_META, datasetStatusFromModel } from "../datasetStatus";
import { DatasetPageHeader } from "../components/shared";
import {
  ScientificSelectionProvider, useScientificSelection,
} from "../components/ScientificSelectionContext";
import SignedDistancePlot from "./viewers/SignedDistancePlot";
import ComparativeBoundaryDashboard from "./viewers/ComparativeBoundaryDashboard";
import { CANON_CLASS_COLOR, CANON_CLASS_LABEL, canonClass } from "./viewers/boundaryClasses";
import { textProps } from "./viewers/semanticStyles";
import { downloadBlob } from "./viewers/plotExport";

// The global "Exon–Domain Boundaries" page is an OVERVIEW and comparison surface,
// not a second copy of the Gene Explorer Boundary tab. For a single-species,
// results-ready dataset it shows a compact heading, six summary metrics, one
// primary signed-distance visualization, a compact class distribution, a short
// inspection summary and a button that opens the full interactive Boundary tab in
// Gene Explorer (carrying the current selection through the shared context). The
// full filter toolbar, evidence table, editable caption and duplicated
// architecture live in Gene Explorer / a Methods drawer / the Figure Gallery.
// Page mode is resolved from the coordinate model — never from a gene symbol.

const CLASS_ORDER = ["exact_domain_edge", "near_domain_edge", "inside_domain",
  "outside_annotated_domains", "unavailable_or_uncertain"];
const METRIC_ORDER = [
  ["total", "Internal boundaries"],
  ["exact_domain_edge", "Exact edge"],
  ["near_domain_edge", "Near edge"],
  ["inside_domain", "Inside domain"],
  ["outside_annotated_domains", "Outside domains"],
  ["unavailable_or_uncertain", "Uncertain"],
];

// Axis ink and end-label typography come from the shared scientific specification and
// are written onto the marks as explicit SVG attributes, so the overview stays legible
// without the page stylesheet.
const AXIS = textProps("axis");
const AXIS_END = textProps("axisEmphasis");

const posOf = (b) => b.protein_position ?? b.boundary_position_aa ?? b.start;
const clsOf = (b) => canonClass(b.boundary_class || b.category || b.class);

export default function GlobalBoundaryDashboard({ model, setPage, openGene }) {
  const coord = useMemo(() => model?.protein_coordinate_model || {}, [model]);
  const dash = coord.boundary_dashboard || {};
  const models = useMemo(() => coord.models || [], [coord]);
  const mode = dash.page_mode;
  const gene = dash.gene_symbol || coord.gene_symbol || model?.analysis?.gene_symbol || "gene";
  const threshold = dash.near_edge_threshold_aa ?? 5;

  const speciesIds = useMemo(() => models.map((m) => m.species_id), [models]);
  const [activeSpecies, setActiveSpecies] = useState(
    (dash.species_available && dash.species_available[0]) || speciesIds[0] || null);

  if (mode === "generic_multi_species_results_ready") {
    // The comparative explorer runs inside the shared selection provider so that a
    // matrix-cell click carries the species, its own boundary and the aa window into
    // the rest of the workspace instead of staying local to this page.
    return (
      <section className="page boundary-dashboard">
        <DatasetPageHeader eyebrow={`Exon–domain boundaries · ${gene}`}
          title="Comparative exon–domain boundaries"
          badges={setPage && <button className="btn ghost sm" onClick={() => setPage("figures")}>Figure gallery →</button>} />
        <ScientificSelectionProvider species={{ species: speciesIds[0] }} model={model}>
          <ComparativeBoundaryDashboard multi={dash.multi_species} gene={gene}
            threshold={threshold} models={models} setPage={setPage} openGene={openGene} />
        </ScientificSelectionProvider>
      </section>
    );
  }

  const speciesObj = { species: activeSpecies };
  return (
    <section className="page boundary-dashboard boundary-overview">
      <ScientificSelectionProvider key={activeSpecies} species={speciesObj} model={model}>
        <SingleSpeciesOverview coord={coord} dash={dash} gene={gene}
          activeSpecies={activeSpecies} speciesIds={speciesIds}
          onSpecies={setActiveSpecies} setPage={setPage} openGene={openGene} />
      </ScientificSelectionProvider>
    </section>
  );
}

function SingleSpeciesOverview({ coord, dash, gene, activeSpecies, speciesIds, onSpecies, setPage, openGene }) {
  const selection = useScientificSelection();
  const speciesModel = useMemo(
    () => (coord.models || []).find((m) => m.species_id === activeSpecies) || (coord.models || [])[0],
    [coord, activeSpecies]);
  const ss = useMemo(() => dash.single_species || {}, [dash]);
  const header = ss.header || {};
  const summary = ss.summary || {};
  const pending = speciesModel?.status !== "available";
  const status = datasetStatusFromModel(coord);
  const [badgeCls, badgeLabel] = DATASET_STATUS_META[status] || DATASET_STATUS_META.unavailable;

  const boundaries = useMemo(() => speciesModel?.exon_boundaries || [], [speciesModel]);
  const boundaryById = useMemo(() => {
    const m = new Map();
    boundaries.forEach((b) => m.set(b.id || b.boundary_id, b));
    return m;
  }, [boundaries]);

  const [classFilter, setClassFilter] = useState(new Set());
  const toggleClass = (c) => {
    const next = new Set(classFilter);
    if (next.has(c)) next.delete(c); else next.add(c);
    setClassFilter(next);
  };

  const total = summary.total || boundaries.length || 0;
  const threshold = header.near_edge_threshold_aa ?? dash.near_edge_threshold_aa ?? 5;
  const length = speciesModel?.protein_length || 1;

  const plotRows = useMemo(() => (
    classFilter.size ? boundaries.filter((b) => classFilter.has(clsOf(b))) : boundaries
  ), [boundaries, classFilter]);

  const selBoundaryId = selection?.selectedBoundaryId || null;
  const clickBoundary = (b) => selection?.selectBoundary?.(b);

  // Inspection cases arrive already priority-ordered from the backend
  // (large distance → mapping → unavailable → incomplete → candidate). We show
  // the 4 highest-priority by default and let the user expand to all real cases.
  const inspection = ss.inspection_cases || [];
  const [showAllCases, setShowAllCases] = useState(false);
  const visibleCases = showAllCases ? inspection : inspection.slice(0, 4);
  const detailRef = useRef(null);

  const selBoundary = selBoundaryId ? boundaryById.get(selBoundaryId) : null;
  const selCase = selBoundaryId
    ? inspection.find((c) => c.boundary_id === selBoundaryId) || null : null;

  const selectCase = (c) => {
    const b = boundaryById.get(c.boundary_id);
    if (b) clickBoundary(b);
    // reveal the selected-boundary detail (Part 14: scroll to / reveal detail)
    window.requestAnimationFrame(() => detailRef.current?.scrollIntoView(
      { behavior: "smooth", block: "center" }));
  };

  const goExplore = () => (openGene
    ? openGene({ species: activeSpecies, tab: "boundary" })
    : setPage && setPage("gene"));
  const goExonMap = () => (openGene
    ? openGene({ species: activeSpecies, tab: "exon" })
    : setPage && setPage("gene"));

  // Compact export for the overview: the class-summary data only.
  const exportSummaryTsv = () => {
    const head = ["class", "count", "near_edge_threshold_aa"].join("\t");
    const rows = METRIC_ORDER.filter(([k]) => k !== "total")
      .map(([k]) => [k, summary[k] || 0, threshold].join("\t"));
    downloadBlob(new Blob([[head, ...rows].join("\n")], { type: "text/tab-separated-values" }),
      `boundary_class_summary_${header.protein_id || "protein"}.tsv`);
  };

  // ---- pre-cluster: compact metadata + positional overview only ---- //
  if (pending) {
    const positions = boundaries.map(posOf).filter((v) => v != null);
    return (
      <>
        <DatasetPageHeader eyebrow={`Exon–domain boundaries · ${gene}`}
          title={<span><i>{header.scientific_name || activeSpecies}</i> — exon–domain boundaries</span>}
          badges={<Badge cls={badgeCls} soft>{badgeLabel}</Badge>} />

        <div className="card bnd-overview-pending">
          <div className="bnd-dash-facts compact">
            <Fact k="Gene" v={header.gene || gene} />
            <Fact k="Species" v={<i>{header.scientific_name || activeSpecies}</i>} />
            <Fact k="Protein" v={<code>{header.protein_id}</code>} />
            <Fact k="Protein length" v={`${header.protein_length ?? length} aa`} />
            <Fact k="Coding exons" v={header.n_coding_exons} />
            <Fact k="Available boundary positions" v={positions.length} />
          </div>
          <p className="muted sm">{positions.length} internal coding-exon boundary positions are available.
            {" "}Domain-edge distances and classifications will be added after the real InterProScan results
            are imported. No signed-distance result or class is shown before the cluster round-trip.</p>

          <PositionOverview positions={positions} length={length} />

          <div className="bnd-overview-actions">
            <button className="btn primary sm" onClick={goExplore}>Open in Gene Explorer →</button>
            {setPage && <button className="btn ghost sm" onClick={() => setPage("runs")}>
              Run the cluster step in My Runs →</button>}
          </div>
        </div>
      </>
    );
  }

  // ---- results-ready single-species overview ---- //
  return (
    <>
      <DatasetPageHeader eyebrow={`Exon–domain boundaries · ${gene}`}
        title={<span><i>{header.scientific_name || activeSpecies}</i> — exon–domain boundaries</span>}
        badges={<>
          <Badge cls={badgeCls} soft>{badgeLabel}</Badge>
          {speciesIds.length > 1 && (
            <label className="bnd-species-sel">Species
              <select value={activeSpecies} onChange={(e) => onSpecies(e.target.value)}>
                {speciesIds.map((s) => <option key={s} value={s}>{s}</option>)}
              </select>
            </label>
          )}
          <Menu label="Export" title="Export overview" align="right">
            <button className="menu-item" onClick={exportSummaryTsv}>Class summary (TSV)</button>
            {setPage && <button className="menu-item" onClick={() => setPage("figures")}>Open Figure Gallery →</button>}
          </Menu>
        </>} />

      <p className="muted sm bnd-overview-sub">
        <code>{header.protein_id}</code>{header.transcript_id ? <> · <code>{header.transcript_id}</code></> : null}
        {" "}· {header.protein_length ?? length} aa · {header.n_coding_exons ?? "?"} coding exons ·
        {" "}{header.representative_domain_source || "representative InterPro domains"} · near-edge ±{threshold} aa
      </p>

      {/* Six compact summary metrics (click a class to focus the plot) */}
      <div className="bnd-metric-row">
        {METRIC_ORDER.map(([k, label]) => {
          const n = k === "total" ? total : (summary[k] || 0);
          const isClass = k !== "total";
          const on = !isClass || classFilter.size === 0 || classFilter.has(k);
          return (
            <button key={k} className={`bnd-metric${on ? "" : " off"}${isClass ? " clickable" : ""}`}
              onClick={isClass ? () => toggleClass(k) : undefined} disabled={!isClass}>
              {isClass && <span className="bnd-metric-sw" style={{ background: CANON_CLASS_COLOR[k] }} />}
              <span className="bnd-metric-n">{n}</span>
              <span className="bnd-metric-l">{label}</span>
            </button>
          );
        })}
      </div>

      {/* One primary visualization: zero-centred signed-distance plot */}
      <div className="card bnd-plot-card">
        <div className="card-head"><h4>Signed distance to nearest representative-domain edge</h4>
          <span className="muted sm">0 = domain edge · shaded band = ±{threshold} aa · click a point to select a boundary</span></div>
        <SignedDistancePlot rows={plotRows} threshold={threshold} sort="position"
          selectedId={selBoundaryId} onSelect={clickBoundary} />
      </div>

      {/* Compact boundary-class distribution */}
      <div className="card bnd-dist-card">
        <div className="card-head"><h4>Boundary-class distribution</h4>
          <span className="muted sm">click a class to focus the plot above</span></div>
        <div className="bnd-dist-bar">
          {CLASS_ORDER.map((c) => {
            const n = summary[c] || 0;
            if (!n) return null;
            const pct = (n / (total || 1)) * 100;
            const on = classFilter.size === 0 || classFilter.has(c);
            return (
              <button key={c} className={`bnd-dist-seg${on ? "" : " off"}`}
                style={{ width: `${pct}%`, background: CANON_CLASS_COLOR[c] }}
                onClick={() => toggleClass(c)}
                title={`${CANON_CLASS_LABEL[c]} · ${n} (${Math.round(pct)}%)`}>
                {pct > 8 ? `${n}` : ""}</button>
            );
          })}
        </div>
      </div>

      {/* Inspection cases — click selects the boundary and reveals its detail below */}
      <div className="card bnd-inspect-card">
        <div className="card-head"><h4>Inspection cases</h4>
          <span className="muted sm">
            {inspection.length
              ? `${visibleCases.length} of ${inspection.length} inspection case${inspection.length > 1 ? "s" : ""} shown — observations worth a closer look, not errors`
              : "No inspection cases for this protein"}</span></div>
        {inspection.length ? (
          <>
            <ul className="bnd-inspect-list compact">
              {visibleCases.map((c) => (
                <li key={c.case_id}>
                  <button
                    className={`bnd-inspect-item${c.boundary_id === selBoundaryId ? " sel" : ""}`}
                    onClick={() => selectCase(c)}
                    title="Select this boundary and reveal its detail (carried into Gene Explorer)">
                    <span className="bnd-inspect-kind">{c.label}</span>
                    <b>{c.boundary_label}</b> <span className="muted sm">· aa {c.protein_position}</span>
                    {c.detail && <span className="bnd-inspect-why muted sm">{c.detail}</span>}
                  </button>
                </li>
              ))}
            </ul>
            {inspection.length > 4 && (
              <button className="btn ghost sm" onClick={() => setShowAllCases((v) => !v)}>
                {showAllCases ? "Show fewer" : `Show all ${inspection.length}`}</button>
            )}
          </>
        ) : null}
      </div>

      {/* Selected-boundary detail (revealed on inspection-case / plot click) */}
      <div className="card bnd-detail-card" ref={detailRef}>
        <div className="card-head"><h4>Selected boundary</h4></div>
        {selBoundary ? (
          <div className="bnd-detail-grid">
            {selCase && (
              <div className="wide bnd-detail-reason">
                <b>Inspection reason:</b> {selCase.label} — {selCase.detail}</div>
            )}
            <div><span className="fld">Exon transition</span>{selBoundary.label
              || `${selBoundary.left_exon_label || "?"} → ${selBoundary.right_exon_label || "?"}`}</div>
            <div><span className="fld">Protein position</span>aa {posOf(selBoundary) ?? "—"}</div>
            <div><span className="fld">Nearest representative domain</span>
              {selBoundary.nearest_domain_short_label
                || selBoundary.nearest_domain_label || "none"}
              {selBoundary.nearest_domain_start != null
                ? ` (aa ${selBoundary.nearest_domain_start}–${selBoundary.nearest_domain_end})` : ""}</div>
            <div><span className="fld">Nearest edge</span>
              {selBoundary.nearest_edge_type || "—"}
              {selBoundary.nearest_edge_position != null ? ` @ aa ${selBoundary.nearest_edge_position}` : ""}</div>
            <div><span className="fld">Signed distance</span>
              {selBoundary.signed_distance ?? "—"}{selBoundary.signed_distance != null ? " aa" : ""}</div>
            <div><span className="fld">Absolute distance</span>
              {selBoundary.absolute_distance ?? (selBoundary.signed_distance != null
                ? Math.abs(selBoundary.signed_distance) : "—")}
              {selBoundary.absolute_distance != null || selBoundary.signed_distance != null ? " aa" : ""}</div>
            <div><span className="fld">Boundary class</span>
              {CANON_CLASS_LABEL[clsOf(selBoundary)] || clsOf(selBoundary) || "—"}</div>
            <div><span className="fld">Mapping status</span>{selBoundary.mapping_status || "—"}</div>
            <div className="wide bnd-detail-actions">
              <button className="btn primary sm" onClick={goExplore}>Open in Gene Explorer Boundary →</button>
              <button className="btn ghost sm" onClick={goExonMap}>Open in Exon Map →</button>
            </div>
          </div>
        ) : (
          <p className="muted sm">Select an inspection case above or a point in the signed-distance plot
            to see the full boundary detail here.</p>
        )}
      </div>

      <div className="bnd-overview-actions">
        <button className="btn primary sm" onClick={goExplore}>Explore individual boundaries →</button>
        <span className="muted sm">Opens the full interactive Boundary tab in Gene Explorer with the current selection.</span>
      </div>
    </>
  );
}

// Minimal positional overview for pre-cluster: protein axis + boundary ticks only.
function PositionOverview({ positions, length }) {
  const W = 960, PAD = 30, H = 60;
  const x = (aa) => PAD + (aa / Math.max(1, length)) * (W - 2 * PAD);
  return (
    <svg className="bnd-pos-overview" viewBox={`0 0 ${W} ${H}`} preserveAspectRatio="xMidYMid meet"
      role="img" aria-label="exon-boundary positions on the protein">
      <rect x="0" y="0" width={W} height={H} fill="#ffffff" />
      <line x1={PAD} x2={W - PAD} y1="34" y2="34" stroke={AXIS.fill} strokeWidth="1" />
      {positions.map((p, i) => (
        <line key={i} x1={x(p)} x2={x(p)} y1="20" y2="48" stroke="#8b98a8" strokeWidth="1.2" strokeDasharray="3 2">
          <title>{`internal coding-exon boundary · aa ${p}`}</title>
        </line>
      ))}
      <text x={PAD} y="14" className="em-axis-lbl"
        fill={AXIS_END.fill} fontSize={AXIS_END.fontSize} fontWeight={AXIS_END.fontWeight}>1</text>
      <text x={W - PAD} y="14" textAnchor="end" className="em-axis-lbl"
        fill={AXIS_END.fill} fontSize={AXIS_END.fontSize}
        fontWeight={AXIS_END.fontWeight}>{length} aa</text>
    </svg>
  );
}

function Fact({ k, v }) {
  return (
    <div className="bnd-fact">
      <span className="bnd-fact-k">{k}</span>
      <span className="bnd-fact-v">{v ?? "—"}</span>
    </div>
  );
}
