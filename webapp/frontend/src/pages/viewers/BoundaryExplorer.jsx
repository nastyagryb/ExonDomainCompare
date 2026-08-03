import { useMemo, useRef, useState } from "react";
import { Badge, Empty, Menu } from "../../ui";
import { useScientificSelection } from "../../components/ScientificSelectionContext";
import { CANON_CLASS_COLOR, CANON_CLASS_LABEL, canonClass } from "./boundaryClasses";
import {
  boundaryProps, domainInstanceFill, featureProps, textProps,
} from "./semanticStyles";
import SignedDistancePlot from "./SignedDistancePlot";
import {
  niceStep, overlaps, downloadBlob,
} from "./plotExport";
import {
  domainInstances, domainFilterOptions, filterBoundaries, boundaryInstanceId,
} from "./common";
import {
  boundaryClassSummaryFigure, boundaryFigure, boundaryTsv, signedDistanceFigure,
} from "./figureData";
import {
  downloadFigurePdf, downloadFigurePng, downloadFigureSvg, downloadFigureTsv,
} from "./figureExport";
import {
  candidateDisplayLayout, candidateLabelFits, candidateTooltip, laneDensity,
} from "./candidateDisplay";

// Integrated exon–domain-boundary explorer for the Gene Explorer Boundary tab.
// It is driven ENTIRELY by the validated protein-coordinate model
// (src/exondomaincompare/shared_gene_analysis/protein_coordinate_model.py) — the same single
// source of truth used by the Exon Map and Domain Architecture — and reuses the
// shared boundary classifier (boundary_classification.py) values, the shared
// ScientificSelectionContext and the shared track toolbar. No coordinate and no
// classification is recomputed in React; pre-cluster runs render an honest
// pending state (positions only, no signed distance / class).

const W = 1000;
// Left gutter reserved for track labels so a label can never overlap a block.
const LBLW = 86;
const PADR = 46;
const CLASS_ORDER = ["exact_domain_edge", "near_domain_edge", "inside_domain",
  "outside_annotated_domains", "unavailable_or_uncertain"];
// Text roles and mark paint come from the shared scientific specification and are
// written onto every mark as explicit SVG attributes, so the figure stays legible
// without the component stylesheet. The stylesheet keeps layout, cursor and hover.
const AXIS = textProps("axis");
const AXIS_END = textProps("axisEmphasis");
const LANE = textProps("trackLabel");
const EMPTY = textProps("empty");
const FEAT = textProps("featureLabel");
const EXON = featureProps("coding_exon");

const posOf = (b) => b.protein_position ?? b.boundary_position_aa ?? b.start;
const absOf = (b) => b.absolute_distance ?? b.absolute_distance_aa
  ?? (b.signed_distance != null ? Math.abs(b.signed_distance) : null);
const clsOf = (b) => canonClass(b.boundary_class || b.category || b.class);

export default function BoundaryExplorer({ model, species, classFilter: classFilterProp, onClassFilterChange }) {
  const selection = useScientificSelection();
  const svgRef = useRef(null);

  const speciesModel = useMemo(() => {
    const models = model?.models || [];
    if (!models.length) return null;
    return models.find((m) => (m.species_id || m.species) === (species || null)) || models[0];
  }, [model, species]);

  const length = speciesModel?.protein_length || 1;
  const exons = useMemo(() => speciesModel?.exons || [], [speciesModel]);
  const domains = useMemo(
    () => domainInstances(speciesModel?.representative_domains || []), [speciesModel]);
  // Domain filter options derived from the real domain instances of this model:
  // "All domains", one "All <entry>s" group per repeated InterPro entry (e.g. the
  // three FGFR1 Ig-like instances) and one option per individual instance. Every
  // option carries the explicit domain_instance_id set it selects.
  const domainOptions = useMemo(
    () => domainFilterOptions(speciesModel?.representative_domains || []), [speciesModel]);
  const domainGroups = useMemo(
    () => domainOptions.filter((o) => o.kind === "group"), [domainOptions]);
  const candidates = useMemo(() => speciesModel?.candidate_regions || [], [speciesModel]);
  const boundaries = useMemo(() => speciesModel?.exon_boundaries || [], [speciesModel]);
  const threshold = speciesModel?.near_edge_threshold_aa ?? 5;
  const pending = speciesModel?.status !== "available";

  // ---- view window (shared toolbar semantics) ---- //
  const [view, setView] = useState([1, length]);
  const [lo, hi] = view;
  const span = Math.max(1, hi - lo);
  const x = (aa) => LBLW + ((Number(aa) - lo) / span) * (W - LBLW - PADR);
  const clampView = ([a, b]) => {
    let s = Math.max(1, Math.round(a));
    let e = Math.min(length, Math.round(b));
    if (e - s < 8) { e = Math.min(length, s + 8); s = Math.max(1, e - 8); }
    return [s, e];
  };
  const zoomBy = (f) => { const c = (lo + hi) / 2, half = (span / f) / 2; setView(clampView([c - half, c + half])); };
  const pan = (frac) => {
    const d = span * frac;
    if (lo + d < 1) setView(clampView([1, 1 + span]));
    else if (hi + d > length) setView(clampView([length - span, length]));
    else setView(clampView([lo + d, hi + d]));
  };

  const [tracks, setTracks] = useState({ domains: true, exons: true, boundaries: true, candidates: true });
  const [showAllCandidates, setShowAllCandidates] = useState(false);
  const toggle = (k) => setTracks((t) => {
    const next = { ...t, [k]: !t[k] };
    selection?.setVisibleTracks?.(next);
    return next;
  });

  // ---- filters (Part 5) ---- //
  // classFilter can be optionally controlled by a parent (global dashboard) so
  // its summary strip + class-distribution bar share ONE filter with this view.
  const [classFilterState, setClassFilterState] = useState(new Set()); // empty == all
  const classFilter = classFilterProp ?? classFilterState;
  const setClassFilter = onClassFilterChange ?? setClassFilterState;
  // domainFilter holds "all", "inst:<domain_instance_id>" for one feature instance,
  // or "grp:<accession>" for every instance of one InterPro entry. Both non-"all"
  // forms resolve through the option's explicit domain_instance_id list.
  const [domainFilter, setDomainFilter] = useState("all");
  const [mappingFilter, setMappingFilter] = useState("all");        // all|mapped|unmapped
  const [exonFilter, setExonFilter] = useState("all");              // boundary id (E1 → E2)
  const [candOnly, setCandOnly] = useState(false);
  const [distMin, setDistMin] = useState("");
  const [distMax, setDistMax] = useState("");
  const [sort, setSort] = useState("position");                     // position|distance
  const resetFilters = () => {
    setClassFilter(new Set()); setDomainFilter("all"); setMappingFilter("all");
    setExonFilter("all"); setCandOnly(false); setDistMin(""); setDistMax("");
    setSort("position");
  };
  const activeFilterCount = (classFilter.size ? 1 : 0) + (domainFilter !== "all" ? 1 : 0)
    + (mappingFilter !== "all" ? 1 : 0) + (exonFilter !== "all" ? 1 : 0)
    + (candOnly ? 1 : 0) + (distMin !== "" ? 1 : 0) + (distMax !== "" ? 1 : 0);
  const toggleClass = (c) => {
    const next = new Set(classFilter);
    if (next.has(c)) next.delete(c); else next.add(c);
    setClassFilter(next);
  };

  // ---- selection (linked, persists across tabs) ---- //
  const selBoundaryId = selection?.selectedBoundaryId || null;
  const selExonId = selection?.selectedExonId || null;
  const selDomainId = selection?.selectedDomainId || null;
  const alignS = selection?.selectedAlignmentStart;
  const alignE = selection?.selectedAlignmentEnd;
  const selCandidate = useMemo(() => {
    if (alignS == null || alignE == null) return null;
    const reg = { start: alignS, end: alignE };
    return candidates.find((c) => overlaps(c, reg)) || null;
  }, [candidates, alignS, alignE]);

  const selBoundary = useMemo(
    () => boundaries.find((b) => (b.id || b.boundary_id) === selBoundaryId) || null,
    [boundaries, selBoundaryId]);

  const candForOverlap = selCandidate || candidates[0] || null;
  // ONE central filtered dataset (Part 11). Every downstream consumer — summary
  // counts, architecture plot, signed-distance plot, evidence table and the
  // visible TSV export — reads exactly this array, so no view can disagree.
  const filteredBoundaries = useMemo(() => filterBoundaries(boundaries, {
    domainFilter, domainOptions, mappingFilter, exonFilter, distMin, distMax,
    candidate: candForOverlap, candOnly, classFilter, sort, classOf: clsOf,
  }), [boundaries, domainFilter, domainOptions, mappingFilter, exonFilter, distMin, distMax,
    candOnly, candForOverlap, classFilter, sort]);
  const filtered = filteredBoundaries;

  const counts = useMemo(() => {
    const c = { total: filteredBoundaries.length };
    for (const k of CLASS_ORDER) c[k] = 0;
    for (const b of filteredBoundaries) c[clsOf(b)] += 1;
    return c;
  }, [filteredBoundaries]);

  // boundaries highlighted from a cross-tab selection (exon / domain / candidate)
  const linkedIds = useMemo(() => {
    const ids = new Set();
    if (selExonId) boundaries.filter((b) => b.left_exon_id === selExonId || b.right_exon_id === selExonId)
      .forEach((b) => ids.add(b.id));
    if (selDomainId) {
      const dom = domains.find((d) => d.id === selDomainId) || null;
      boundaries.filter((b) => b.nearest_domain_id === selDomainId
        || (dom && boundaryInstanceId(b) === dom.instanceId)).forEach((b) => ids.add(b.id));
    }
    return ids;
  }, [boundaries, domains, selExonId, selDomainId]);

  const clickBoundary = (b) => selection?.selectBoundary?.(b);

  if (!speciesModel) {
    return <Empty title="Boundary analysis not available"
      hint="No validated protein-coordinate model was built for this run." />;
  }

  // Repeated domain instances stay separable through the shared instance ramp, so a
  // domain keeps the same colour here and in the exported figure.
  const domColorMap = {};
  domains.forEach((d, i) => { domColorMap[d.id] = domainInstanceFill(i + 1); });

  // ---- architecture-on-boundary SVG geometry ---- //
  const RULER_H = 30, DOM_Y = 40, EX_Y = 70, LANE_H = 20, CONN_Y = 98, CAND_Y = 112;
  // Candidate clusters, lanes, labels and tooltips come from the one shared display
  // module the Domain Architecture track and the exported figures use, so a candidate
  // has the same cluster id and the same lane wherever it is drawn. Boundary is a
  // narrower track, so it grows for the lanes it shows rather than stacking boxes.
  const candLayout = candidateDisplayLayout(candidates, {
    selectedId: candForOverlap?.id || null, showAll: showAllCandidates,
  });
  const CAND_LANE_H = 10, CAND_LANE_GAP = 3;
  const candLaneY = (lane) => CAND_Y + lane * (CAND_LANE_H + CAND_LANE_GAP);
  const candExtra = tracks.candidates
    ? Math.max(0, candLayout.laneCount * (CAND_LANE_H + CAND_LANE_GAP) - CAND_LANE_GAP - CAND_LANE_H)
    : 0;
  const H = 150 + candExtra;
  const major = niceStep(span), minor = major / 5;
  const majors = [], minors = [];
  for (let t = Math.ceil(lo / minor) * minor; t <= hi; t += minor) {
    if (Math.abs(t / major - Math.round(t / major)) < 1e-6) majors.push(Math.round(t));
    else minors.push(Math.round(t));
  }
  const barW = (a, b) => Math.max(2, x(b) - x(a));
  const inView = (v) => v >= lo && v <= hi;

  // ---- exports (figure + visible boundary table) ---- //
  const stem = `boundary_${speciesModel.protein_id}`;
  // Publication output comes from the shared figure specification and is fed the
  // *filtered* boundary set, so a downloaded figure shows exactly the boundaries
  // the table and the summary counts show. Nothing here screenshots the live SVG.
  const buildBoundaryFig = () => boundaryFigure(speciesModel, filteredBoundaries,
    { selectedBoundaryId: selBoundaryId, selectedCandidateId: candForOverlap?.id || null,
      showAllCandidates });
  const buildSignedFig = () => signedDistanceFigure(speciesModel, filteredBoundaries,
    { selectedBoundaryId: selBoundaryId });
  const buildClassFig = () => boundaryClassSummaryFigure(speciesModel, filteredBoundaries);
  const exportSvg = () => downloadFigureSvg(buildBoundaryFig(), stem);
  const exportPdf = () => downloadFigurePdf(buildBoundaryFig(), stem);
  const exportPng = () => downloadFigurePng(buildBoundaryFig(), stem);
  const exportSignedSvg = () => downloadFigureSvg(buildSignedFig(), `${stem}_signed_distance`);
  const exportSignedPdf = () => downloadFigurePdf(buildSignedFig(), `${stem}_signed_distance`);
  const exportClassPdf = () => downloadFigurePdf(buildClassFig(), `${stem}_class_summary`);
  const exportVisibleTsv = () => downloadFigureTsv(boundaryTsv(filteredBoundaries), stem);
  const exportTsv = () => {
    const head = ["boundary_id", "label", "protein_position", "left_exon", "right_exon",
      "nearest_domain", "nearest_domain_instance_id", "nearest_domain_instance_number",
      "nearest_domain_accession", "nearest_domain_start", "nearest_domain_end",
      "nearest_edge_type", "nearest_edge_position", "signed_distance", "absolute_distance",
      "boundary_class", "near_threshold", "mapping_status", "source", "source_file"].join("\t");
    const rows = filtered.map((b) => [
      b.id, b.label, posOf(b), b.left_exon_label, b.right_exon_label,
      b.nearest_domain_full_label ?? b.nearest_domain_label ?? "",
      boundaryInstanceId(b) ?? "", b.nearest_domain_instance_number ?? "",
      b.nearest_domain_accession ?? "", b.nearest_domain_start ?? "",
      b.nearest_domain_end ?? "", b.nearest_edge_type ?? "", b.nearest_edge_position ?? "",
      b.signed_distance ?? "", absOf(b) ?? "", clsOf(b), b.near_threshold ?? threshold,
      b.mapping_status ?? "", b.source ?? "", b.source_file ?? "",
    ].join("\t"));
    downloadBlob(new Blob([[head, ...rows].join("\n")], { type: "text/tab-separated-values" }),
      `boundaries_${speciesModel.protein_id}.tsv`);
  };
  const nClassified = boundaries.filter((b) => b.signed_distance != null).length;

  return (
    <div className="viewer coord-viewer exon-map boundary-explorer">
      <div className="viewer-head">
        <div>
          <b>Exon–domain boundaries</b>
          <p className="muted sm">{speciesModel.scientific_name} ·
            {" "}<code>{speciesModel.protein_id}</code>
            {speciesModel.transcript_id ? <> · <code>{speciesModel.transcript_id}</code></> : null}
            {" "}· {length} aa · {boundaries.length} internal coding-exon boundaries</p>
        </div>
        {/* The badge reads the same filtered dataset as every plot and table, so
            a total can never be mixed with filtered class counts (Part 11). */}
        <Badge cls={pending ? "neutral" : "accepted"} soft>
          {pending ? "domain-edge distances pending cluster"
            : `${filteredBoundaries.length} of ${boundaries.length} shown · `
              + `${counts.inside_domain} inside · ${counts.near_domain_edge} near · `
              + `${counts.outside_annotated_domains} outside`}</Badge>
      </div>

      {/* Compact toolbar: Zoom | Fit protein | Tracks▾ | Export▾ */}
      <div className="em-toolbar compact-toolbar">
        <div className="seg">
          <button className="seg-btn" onClick={() => zoomBy(1.6)} title="Zoom in">＋</button>
          <button className="seg-btn" onClick={() => zoomBy(1 / 1.6)} title="Zoom out">－</button>
          <button className="seg-btn" onClick={() => pan(-0.25)} title="Pan left">◀</button>
          <button className="seg-btn" onClick={() => pan(0.25)} title="Pan right">▶</button>
        </div>
        <button className="seg-btn" onClick={() => setView([1, length])} title="Fit whole protein">Fit protein</button>
        {selBoundary && <button className="seg-btn"
          onClick={() => setView(clampView([posOf(selBoundary) - 60, posOf(selBoundary) + 60]))}
          title="Zoom to selected boundary">Zoom to {selBoundary.label}</button>}
        {tracks.candidates && candLayout.hiddenCount > 0 && (
          <button className="seg-btn" onClick={() => setShowAllCandidates(true)}
            title="Show every exploratory candidate cluster, including lower-ranked lanes">
            Show all candidates ({candLayout.hiddenCount} more)</button>
        )}
        {tracks.candidates && showAllCandidates && candLayout.total > candLayout.laneCount && (
          <button className="seg-btn on" onClick={() => setShowAllCandidates(false)}
            title="Show only the top-ranked candidate lanes">Top-ranked candidates only</button>
        )}
        <Menu label="Tracks" title="Show / hide tracks">
          {[["domains", "Representative domains"], ["exons", "Coding exons"],
            ["boundaries", "Exon boundaries"], ["candidates", "Candidate regions"]].map(([k, l]) => (
            <label key={k} className="menu-check">
              <input type="checkbox" checked={tracks[k]} onChange={() => toggle(k)} /><span>{l}</span>
            </label>
          ))}
        </Menu>
        <Menu label="Export" title="Export figure and data" align="right">
          <button className="menu-item" onClick={exportSvg}>Boundary figure — SVG (vector)</button>
          <button className="menu-item" onClick={exportPdf}>Boundary figure — PDF (vector)</button>
          <button className="menu-item" onClick={exportPng}>Boundary figure — PNG (300 dpi)</button>
          <div className="menu-sep" />
          <button className="menu-item" onClick={exportSignedSvg}>Signed distances — SVG (vector)</button>
          <button className="menu-item" onClick={exportSignedPdf}>Signed distances — PDF (vector)</button>
          <button className="menu-item" onClick={exportClassPdf}>Class summary — PDF (vector)</button>
          <div className="menu-sep" />
          <button className="menu-item" onClick={exportVisibleTsv}>Visible boundaries (TSV)</button>
          <button className="menu-item" onClick={exportTsv}>Full boundary table (TSV)</button>
        </Menu>
        <span className="spacer" />
        <span className="muted small">visible aa {lo}–{hi}</span>
      </div>

      {selBoundary && (
        <div className="em-selsum muted small">
          Selected: <b>{selBoundary.label}</b> · aa {posOf(selBoundary)}
          {!pending && <> · {CANON_CLASS_LABEL[clsOf(selBoundary)]} · signed {selBoundary.signed_distance} aa
            {" "}to the {selBoundary.nearest_edge_type} edge of
            {" "}{selBoundary.nearest_domain_full_label || selBoundary.nearest_domain_label}</>}
        </div>
      )}

      {pending ? (
        /* Honest pre-cluster state: short explanation + positional context only.
           No disabled class chips, no filters, no signed-distance legend, no
           classification, no large evidence table (Part D / Part 11). */
        <div className="pending-note compact">
          <Badge cls="neutral" soft>pending cluster</Badge>
          <span>{boundaries.length ? `${boundaries.length} internal coding-exon boundary positions are available.`
            : "Coding-exon boundary positions are available."}{" "}
            Domain-edge distances and classifications will be added after the real InterProScan results
            are imported. The cluster command is available in <b>My Runs</b>.</span>
        </div>
      ) : (
        /* Compact summary strip — class chips double as the shared class filter */
        <div className="bnd-summary">
          <button className={`bnd-chip${classFilter.size === 0 ? " on" : ""}`} onClick={() => setClassFilter(new Set())}>
            <span className="bnd-chip-n">{counts.total}</span><span className="bnd-chip-l">internal boundaries</span></button>
          {CLASS_ORDER.map((c) => (
            <button key={c} className={`bnd-chip${classFilter.has(c) ? " on" : ""}`}
              onClick={() => toggleClass(c)}>
              <span className="bnd-chip-sw" style={{ background: CANON_CLASS_COLOR[c] }} />
              <span className="bnd-chip-n">{counts[c]}</span>
              <span className="bnd-chip-l">{CANON_CLASS_LABEL[c]}</span>
            </button>
          ))}
        </div>
      )}

      {/* Integrated boundary-on-architecture (Part 2) */}
      <div className="em-canvas">
        <svg ref={svgRef} className="em-svg" viewBox={`0 0 ${W} ${H}`} role="img"
          aria-label={`${speciesModel.protein_id} exon–domain boundary architecture`}
          preserveAspectRatio="xMidYMid meet">
          <rect x="0" y="0" width={W} height={H} fill="#ffffff" />

          {/* aa axis ruler */}
          <line x1={LBLW} x2={W - PADR} y1={RULER_H - 8} y2={RULER_H - 8} stroke={AXIS.fill} strokeWidth="1" />
          {minors.filter(inView).map((t) => (
            <line key={`mn${t}`} x1={x(t)} x2={x(t)} y1={RULER_H - 11} y2={RULER_H - 5} stroke="#8b98a8" strokeWidth="0.6" />
          ))}
          {majors.filter(inView).map((t) => (
            <g key={`mj${t}`}>
              <line x1={x(t)} x2={x(t)} y1={RULER_H - 14} y2={RULER_H - 2} stroke={AXIS.fill} strokeWidth="1" />
              <text x={x(t)} y={RULER_H - 18} textAnchor="middle" className="em-axis-lbl"
                fill={AXIS.fill} fontSize={AXIS.fontSize}>{t}</text>
            </g>
          ))}

          <text x={LBLW - 10} y={DOM_Y + LANE_H / 2 + 3} textAnchor="end" className="da-lane-lbl"
            fill={LANE.fill} fontSize={LANE.fontSize}>Domains</text>
          <text x={LBLW - 10} y={EX_Y + LANE_H / 2 + 3} textAnchor="end" className="da-lane-lbl"
            fill={LANE.fill} fontSize={LANE.fontSize}>Exons</text>

          {/* representative domains */}
          {tracks.domains && domains.map((d) => {
            const sel = d.id === selDomainId
              || (selBoundary && boundaryInstanceId(selBoundary) === d.instanceId);
            return (
              <rect key={d.id} x={x(d.start)} y={DOM_Y} width={barW(d.start, d.end)} height={LANE_H} rx="3"
                fill={domColorMap[d.id]} opacity={sel ? 1 : 0.82}
                stroke={sel ? "#12151a" : "#2b2f36"} strokeWidth={sel ? 1.4 : 0.4}
                onClick={() => selection?.selectDomain?.(d)} style={{ cursor: "pointer" }}>
                <title>{`${d.instanceLabel} · ${d.interpro_accession || ""}\n`
                  + `instance ${d.instanceId}`}</title>
              </rect>
            );
          })}
          {tracks.domains && !domains.length && (
            <text x={W / 2} y={DOM_Y + LANE_H / 2 + 3} textAnchor="middle" className="da-empty"
              fill={EMPTY.fill} fontSize={EMPTY.fontSize} fontStyle={EMPTY.fontStyle}>
              {pending ? "representative domains pending post-cluster InterProScan"
                : "no representative InterPro domain"}</text>)}

          {/* coding exons */}
          {tracks.exons && exons.map((ex) => {
            const bx = x(ex.start), bw = barW(ex.start, ex.end);
            const adj = selBoundary && (ex.id === selBoundary.left_exon_id || ex.id === selBoundary.right_exon_id);
            const sel = ex.id === selExonId;
            // An exon keeps its scientific fill when it is selected or adjacent to the
            // selected boundary; the emphasis is the outline only.
            const paint = featureProps("coding_exon", { selected: sel || adj });
            return (
              <g key={ex.id}>
                <rect x={bx} y={EX_Y} width={bw} height={LANE_H} rx="2"
                  className={`em-exon${sel || adj ? " sel" : ""}`}
                  fill={paint.fill} fillOpacity={paint.fillOpacity}
                  stroke={paint.stroke} strokeWidth={paint.strokeWidth}
                  onClick={() => selection?.selectExon?.({
                    exon_id: ex.id, transcript_id: ex.tooltip?.transcript_id,
                    protein_id: speciesModel.protein_id, protein_start_aa: ex.start, protein_end_aa: ex.end })}
                  style={{ cursor: "pointer" }}>
                  <title>{`${ex.label} · protein aa ${ex.start}–${ex.end}`}</title>
                </rect>
                {bw >= 16 && <text x={bx + bw / 2} y={EX_Y + LANE_H / 2 + 3} textAnchor="middle"
                  className="em-exon-num" fill={FEAT.fill} fontSize={FEAT.fontSize}
                  fontWeight={FEAT.fontWeight} style={{ pointerEvents: "none" }}>{ex.label}</text>}
              </g>
            );
          })}

          {/* connector from selected boundary to its nearest domain edge */}
          {!pending && selBoundary && selBoundary.nearest_edge_position != null && inView(posOf(selBoundary)) && (
            <g>
              <line x1={x(posOf(selBoundary))} x2={x(selBoundary.nearest_edge_position)}
                y1={CONN_Y} y2={CONN_Y} stroke={boundaryProps(clsOf(selBoundary)).stroke} strokeWidth="1.6" />
              <circle cx={x(selBoundary.nearest_edge_position)} cy={CONN_Y} r="3"
                fill={boundaryProps(clsOf(selBoundary)).fill} />
              {/* The signed distance, nearest domain and nearest edge are reported in the
                  selected-boundary summary above the plot, the tooltip and the detail
                  panel — never drawn on top of domains or exons (Part 10). */}
              <title>{`${selBoundary.label} · signed ${selBoundary.signed_distance} aa to the `
                + `${selBoundary.nearest_edge_type} edge of `
                + `${selBoundary.nearest_domain_full_label || selBoundary.nearest_domain_label}`}</title>
            </g>
          )}

          {/* boundary markers — all active filters narrow the architecture view */}
          {tracks.boundaries && filtered.map((b) => {
            const p = posOf(b);
            if (!inView(p)) return null;
            const cls = clsOf(b);
            const sel = (b.id || b.boundary_id) === selBoundaryId;
            const linked = linkedIds.has(b.id);
            // Pre-cluster the class is not measured yet, so the marker stays neutral
            // grey and dashed instead of borrowing a class colour.
            const col = pending ? "#8b98a8" : boundaryProps(cls).fill;
            return (
              <g key={b.id} onClick={() => clickBoundary(b)} style={{ cursor: "pointer" }}>
                <line x1={x(p)} x2={x(p)} y1={DOM_Y - 6} y2={EX_Y + LANE_H + 6}
                  stroke={col} strokeWidth={sel ? 3 : linked ? 2 : 1.2}
                  opacity={sel ? 1 : linked ? 0.95 : 0.8} strokeDasharray={pending ? "3 3" : undefined} />
                <circle cx={x(p)} cy={DOM_Y - 8} r={sel ? 4.5 : 3} fill={col}
                  stroke={sel ? "#12151a" : "none"} strokeWidth="1" />
                <title>{`${b.label} · aa ${p}\n`
                  + (pending ? "domain-edge distance pending cluster"
                    : `nearest ${b.nearest_domain_full_label || b.nearest_domain_label} `
                      + `(aa ${b.nearest_domain_start}–${b.nearest_domain_end}) `
                      + `${b.nearest_edge_type} edge\nsigned ${b.signed_distance} aa · |${absOf(b)}| · `
                      + `${CANON_CLASS_LABEL[cls]}`)
                  + `\nsource ${b.source || "—"}`}</title>
              </g>
            );
          })}

          {/* optional candidate overlay — same display clusters and lanes as the
              Domain Architecture track and the exported figures */}
          {tracks.candidates && candLayout.visible.map((c) => {
            const sel = Boolean(candForOverlap && c.id === candForOverlap.id);
            const paint = featureProps("candidate_region", { selected: sel, opacity: sel ? 0.5 : 0.25 });
            const w = barW(c.start, c.end);
            const y = candLaneY(candLayout.laneOf(c));
            return (
              <g key={c.id}>
                <rect x={x(c.start)} y={y} width={w} height={CAND_LANE_H} rx="2"
                  fill={paint.fill} fillOpacity={paint.fillOpacity}
                  stroke={paint.stroke} strokeWidth={paint.strokeWidth}
                  style={{ cursor: "pointer" }}
                  onClick={() => { selection?.selectAlignmentRegion?.(c.start, c.end);
                    selection?.setCoordinateRange?.(c.start, c.end); }}>
                  <title>{candidateTooltip(c)}</title>
                </rect>
                {candidateLabelFits(w) && (
                  <text x={(x(c.start) + x(c.end)) / 2} y={y + CAND_LANE_H - 2}
                    textAnchor="middle" className="em-exon-num" fill={FEAT.fill}
                    fontSize={FEAT.fontSize} fontWeight={FEAT.fontWeight}
                    style={{ pointerEvents: "none" }}>{c.id}</text>
                )}
              </g>
            );
          })}
          {/* clusters too narrow to read at this scale are counted, not smeared */}
          {tracks.candidates && candLayout.byLane.map((laneItems, i) => {
            const d = laneDensity(laneItems, (c) => barW(c.start, c.end));
            if (d.narrow < 2) return null;
            return (
              <text key={`cd${i}`} x={x(d.end) + 6} y={candLaneY(i) + CAND_LANE_H - 2}
                className="em-axis-lbl" fill={FEAT.fill} fontSize={FEAT.fontSize}
                style={{ pointerEvents: "none" }}>{d.narrow} narrow clusters (aa {d.start}–{d.end})</text>
            );
          })}

          <text x={LBLW} y={H - 4} textAnchor="start" className="em-axis-end"
            fill={AXIS_END.fill} fontSize={AXIS_END.fontSize} fontWeight={AXIS_END.fontWeight}>1</text>
          <text x={W - PADR} y={H - 4} textAnchor="end" className="em-axis-end"
            fill={AXIS_END.fill} fontSize={AXIS_END.fontSize}
            fontWeight={AXIS_END.fontWeight}>{length} aa</text>
        </svg>
      </div>

      {!pending && (
        <div className="legend res-legend">
          {CLASS_ORDER.map((c) => (
            <span key={c} className="legend-item"><span className="pa-swatch" style={{ background: CANON_CLASS_COLOR[c] }} />
              {CANON_CLASS_LABEL[c]}</span>
          ))}
          <span className="muted small">0 = domain edge · ▶ start / ◀ end edge on the signed-distance plot</span>
        </div>
      )}

      {/* Zero-centred signed-distance plot (Part 3) */}
      {!pending && nClassified > 0 && (
        <div className="card bnd-plot-card">
          <div className="card-head"><h4>Signed distance to nearest representative-domain edge</h4>
            <span className="muted sm">0 = domain edge · shaded band = ±{threshold} aa · click a point to link everything</span></div>
          <SignedDistancePlot rows={filtered} threshold={threshold} sort={sort}
            selectedId={selBoundaryId} onSelect={(b) => clickBoundary(b)} />
        </div>
      )}

      {/* Boundary detail panel (Part 4) — default visible for a selection */}
      {!pending && selBoundary && (
        <BoundaryDetail b={selBoundary} model={speciesModel} exons={exons} domains={domains}
          candidates={candidates} threshold={threshold} pending={pending} domColorMap={domColorMap} />
      )}

      {/* Advanced analysis — filters, full evidence table and provenance are
          progressively disclosed so the default view stays focused (Part D). */}
      {!pending && (
        <details className="card bnd-advanced">
          <summary><b>Advanced analysis</b><span className="muted sm"> — filters and the full evidence table ·
            {" "}{filtered.length} of {boundaries.length} boundaries
            {activeFilterCount > 0 ? ` · ${activeFilterCount} active filter${activeFilterCount > 1 ? "s" : ""}` : ""}</span></summary>

          <div className="bnd-filters">
            {/* One option per real domain feature INSTANCE plus an "all instances of
                this InterPro entry" option for repeated entries. Selection always
                resolves through nearest_domain_instance_id, never an accession. */}
            <label>Domain
              <select value={domainFilter} onChange={(e) => setDomainFilter(e.target.value)} disabled={!domains.length}>
                <option value="all">All domains</option>
                {domainGroups.map((g) => (
                  <option key={g.value} value={g.value}>{g.label}</option>
                ))}
                {domains.map((d) => (
                  <option key={d.instanceId} value={`inst:${d.instanceId}`}>{d.instanceLabel}</option>
                ))}
              </select>
            </label>
            <label>Mapping
              <select value={mappingFilter} onChange={(e) => setMappingFilter(e.target.value)}>
                <option value="all">all</option><option value="mapped">mapped</option>
                <option value="unmapped">unmapped / pending</option>
              </select>
            </label>
            <label>Boundary
              <select value={exonFilter} onChange={(e) => setExonFilter(e.target.value)}>
                <option value="all">all</option>
                {boundaries.map((b) => (
                  <option key={b.id || b.boundary_id} value={b.id || b.boundary_id}>{b.label}</option>
                ))}
              </select>
            </label>
            <label>|dist| ≥ <input type="number" className="bnd-num" value={distMin}
              onChange={(e) => setDistMin(e.target.value)} /></label>
            <label>|dist| ≤ <input type="number" className="bnd-num" value={distMax}
              onChange={(e) => setDistMax(e.target.value)} /></label>
            {candForOverlap && <label className="check inline">
              <input type="checkbox" checked={candOnly} onChange={() => setCandOnly((v) => !v)} />
              <span>overlaps {candForOverlap.id}</span></label>}
            <label>Sort
              <select value={sort} onChange={(e) => setSort(e.target.value)}>
                <option value="position">protein position</option>
                <option value="distance">absolute distance</option>
              </select>
            </label>
            <button className="btn ghost sm" onClick={resetFilters} disabled={!activeFilterCount}>
              Reset all filters{activeFilterCount ? ` (${activeFilterCount})` : ""}</button>
            <span className="muted small"><b>{filtered.length}</b> of {boundaries.length} boundaries shown</span>
          </div>

          <div className="table-scroll"><table className="mini-tbl bnd-table">
            <thead><tr>
              <th>Boundary</th><th>Position</th><th>Nearest domain</th><th>Edge</th>
              <th>Signed</th><th>|dist|</th><th>Class</th><th>Mapping</th><th>Source</th>
            </tr></thead>
            <tbody>{filtered.map((b) => {
              const sel = (b.id || b.boundary_id) === selBoundaryId;
              const cls = clsOf(b);
              return (
                <tr key={b.id} className={sel ? "row-selected" : ""} onClick={() => clickBoundary(b)}
                  style={{ cursor: "pointer" }}>
                  <td>{b.label}</td>
                  <td>{posOf(b)}</td>
                  <td>{b.nearest_domain_short_label || b.nearest_domain_label || "—"}
                    {b.nearest_domain_start != null && <span className="muted sm"> (aa {b.nearest_domain_start}–{b.nearest_domain_end})</span>}</td>
                  <td>{b.nearest_edge_type || "—"}{b.nearest_edge_position != null ? ` @${b.nearest_edge_position}` : ""}</td>
                  <td>{b.signed_distance ?? "—"}</td>
                  <td>{absOf(b) ?? "—"}</td>
                  <td><span className="class-chip" style={{ background: CANON_CLASS_COLOR[cls] }} />
                    <span className="muted sm">{CANON_CLASS_LABEL[cls]}</span></td>
                  <td>{b.mapping_status || "—"}</td>
                  <td className="muted sm">{b.source || "—"}</td>
                </tr>);
            })}</tbody>
          </table></div>
          {!filtered.length && <p className="muted sm pad">No boundaries match the current filters.</p>}
        </details>
      )}
    </div>
  );
}

// ---- boundary detail panel with a local architecture zoom (Part 4) ---- //
function BoundaryDetail({ b, model, exons, domains, candidates, threshold, pending, domColorMap }) {
  const pos = posOf(b);
  const cls = clsOf(b);
  const dom = domains.find((d) => d.instanceId === boundaryInstanceId(b))
    || domains.find((d) => d.id === b.nearest_domain_id) || null;
  const left = exons.find((e) => e.id === b.left_exon_id) || null;
  const right = exons.find((e) => e.id === b.right_exon_id) || null;
  const candOverlap = candidates.find((c) => pos >= c.start && pos <= c.end) || null;

  // local zoom window around the boundary
  const halfWin = Math.max(40, (b.absolute_distance ?? b.absolute_distance_aa ?? 30) + 25);
  const zlo = Math.max(1, pos - halfWin), zhi = Math.min(model.protein_length, pos + halfWin);
  const ZW = 620, ZPAD = 30, ZH = 96;
  const zx = (aa) => ZPAD + ((aa - zlo) / Math.max(1, zhi - zlo)) * (ZW - 2 * ZPAD);

  // Name the domain feature INSTANCE, so a repeated InterPro entry can never make
  // the interpretation read as if a different instance had been measured.
  const domName = b.nearest_domain_full_label || b.nearest_domain_short_label
    || b.nearest_domain_label;
  const interp = pending
    ? `Exon boundary ${b.label} lies at protein position ${pos}. Its distance to any representative `
      + `domain edge is pending post-cluster InterProScan and is therefore not classified yet.`
    : cls === "inside_domain"
      ? `This coding-exon boundary lies ${b.absolute_distance ?? Math.abs(b.signed_distance)} amino acids inside `
        + `${domName} and is therefore classified as inside_domain.`
      : cls === "exact_domain_edge"
        ? `This coding-exon boundary coincides exactly with the ${b.nearest_edge_type} edge of `
          + `${domName} (0 aa) and is classified as exact_domain_edge.`
        : cls === "near_domain_edge"
          ? `This coding-exon boundary lies ${b.absolute_distance} amino acids from the ${b.nearest_edge_type} edge `
            + `of ${domName} (≤ ${threshold} aa) and is classified as near_domain_edge.`
          : cls === "outside_annotated_domains"
            ? `This coding-exon boundary lies ${b.absolute_distance} amino acids `
              + `${b.signed_distance < 0 ? "N-terminal of" : "C-terminal of"} the nearest edge of `
              + `${domName}, outside all annotated representative domains `
              + `(outside_annotated_domains).`
            : `No representative domain is available for this boundary (unavailable_or_uncertain).`;

  return (
    <div className="card bnd-detail">
      <div className="card-head">
        <h4>{b.label} · coding-exon boundary</h4>
        {!pending && <Badge cls="neutral" soft style={{ background: CANON_CLASS_COLOR[cls] }}>{CANON_CLASS_LABEL[cls]}</Badge>}
      </div>
      <div className="bnd-detail-grid">
        <div><span className="fld">Protein position</span>aa {pos}</div>
        <div><span className="fld">Protein / transcript</span><code>{model.protein_id}</code> · <code>{model.transcript_id || "—"}</code></div>
        <div><span className="fld">Adjacent exons</span>{left?.label || "?"} (aa {left?.start}–{left?.end}) → {right?.label || "?"} (aa {right?.start}–{right?.end})</div>
        <div><span className="fld">Nearest representative domain</span>
          {b.nearest_domain_short_label || b.nearest_domain_label || "—"}
          {b.nearest_domain_accession ? <> (<code>{b.nearest_domain_accession}</code>)</> : null}</div>
        <div><span className="fld">Domain instance</span>{b.nearest_domain_start != null
          ? `aa ${b.nearest_domain_start}–${b.nearest_domain_end}`
          + (b.nearest_domain_instance_number != null ? ` · instance ${b.nearest_domain_instance_number}` : "")
          : "—"}</div>
        <div><span className="fld">Nearest edge</span>{b.nearest_edge_type || "—"}{b.nearest_edge_position != null ? ` @ aa ${b.nearest_edge_position}` : ""}</div>
        <div><span className="fld">Signed distance</span>{b.signed_distance ?? "—"} aa</div>
        <div><span className="fld">Absolute distance</span>{absOf(b) ?? "—"} aa</div>
        <div><span className="fld">Classification</span>{pending ? "pending" : CANON_CLASS_LABEL[cls]}</div>
        <div><span className="fld">Near-edge threshold</span>{threshold} aa</div>
        <div><span className="fld">Mapping status</span>{b.mapping_status || "—"}</div>
        <div><span className="fld">Candidate overlap</span>{candOverlap ? `${candOverlap.id} (exploratory)` : "none"}</div>
      </div>

      {/* local architecture zoom around the boundary */}
      <svg className="bnd-zoom" viewBox={`0 0 ${ZW} ${ZH}`} preserveAspectRatio="xMidYMid meet"
        role="img" aria-label={`local architecture around ${b.label}`}>
        <rect x="0" y="0" width={ZW} height={ZH} fill="#ffffff" />
        <line x1={ZPAD} x2={ZW - ZPAD} y1="20" y2="20" stroke={AXIS.fill} strokeWidth="0.8" />
        <text x={ZPAD} y="12" className="em-axis-lbl"
          fill={AXIS.fill} fontSize={AXIS.fontSize}>{zlo}</text>
        <text x={ZW - ZPAD} y="12" textAnchor="end" className="em-axis-lbl"
          fill={AXIS.fill} fontSize={AXIS.fontSize}>{zhi}</text>
        {/* domain in window */}
        {dom && (
          <rect x={zx(Math.max(zlo, dom.start))} y="30" width={Math.max(2, zx(Math.min(zhi, dom.end)) - zx(Math.max(zlo, dom.start)))}
            height="18" rx="3" fill={domColorMap[dom.id] || domainInstanceFill(1)} opacity="0.85">
            <title>{dom.instanceLabel || `${dom.label} aa ${dom.start}–${dom.end}`}</title></rect>
        )}
        {/* adjacent exons */}
        {[left, right].filter(Boolean).map((ex) => (
          <rect key={ex.id} x={zx(Math.max(zlo, ex.start))} y="54"
            width={Math.max(2, zx(Math.min(zhi, ex.end)) - zx(Math.max(zlo, ex.start)))} height="16" rx="2"
            className="em-exon" fill={EXON.fill} fillOpacity={EXON.fillOpacity}
            stroke={EXON.stroke} strokeWidth={EXON.strokeWidth}>
            <title>{`${ex.label} aa ${ex.start}–${ex.end}`}</title></rect>
        ))}
        {[left, right].filter(Boolean).map((ex) => (
          <text key={`l${ex.id}`} x={(zx(Math.max(zlo, ex.start)) + zx(Math.min(zhi, ex.end))) / 2} y="66"
            textAnchor="middle" className="em-exon-num" fill={FEAT.fill} fontSize={FEAT.fontSize}
            fontWeight={FEAT.fontWeight}>{ex.label}</text>
        ))}
        {/* nearest edge + boundary + connector */}
        {!pending && b.nearest_edge_position != null && (
          <line x1={zx(b.nearest_edge_position)} x2={zx(b.nearest_edge_position)} y1="26" y2="52"
            stroke="#12151a" strokeWidth="1" strokeDasharray="2 2" />
        )}
        <line x1={zx(pos)} x2={zx(pos)} y1="26" y2="74" stroke={pending ? "#8b98a8" : boundaryProps(cls).stroke}
          strokeWidth="2.4" />
        {!pending && b.nearest_edge_position != null && (
          <line x1={zx(pos)} x2={zx(b.nearest_edge_position)} y1="80" y2="80"
            stroke={boundaryProps(cls).stroke} strokeWidth="1.4" />
        )}
        <text x={zx(pos)} y="88" textAnchor="middle" className="bnd-conn-lbl"
          fill={FEAT.fill} fontSize={FEAT.fontSize}
          fontWeight={FEAT.fontWeight}>{b.label} @ {pos}</text>
      </svg>

      <p className="cand-interp"><b>Interpretation:</b> {interp} This is a coordinate-level classification;
        no functional consequence is inferred.</p>
    </div>
  );
}
