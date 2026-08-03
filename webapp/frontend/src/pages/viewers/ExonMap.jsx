import { useMemo, useRef, useState } from "react";
import { fileUrl } from "../../api";
import { Badge, Empty, Menu } from "../../ui";
import { useScientificSelection } from "../../components/ScientificSelectionContext";
import { exonMapFigure, exonMapTsv } from "./figureData";
import {
  downloadFigurePdf, downloadFigurePng, downloadFigureSvg, downloadFigureTsv,
} from "./figureExport";
import { featureProps, featureStyle, textProps } from "./semanticStyles";

// Interactive Exon Map driven entirely by the validated protein-coordinate model
// (src/exondomaincompare/shared_gene_analysis/protein_coordinate_model.py). This is the single
// source of truth for the primary exon-to-protein projection, exploratory
// candidate overlays, optional real domain context, and inline transcript-model
// comparison. Fully gene-agnostic — no hard-coded gene / transcript IDs.

const W = 1000;      // logical px width of the plot area
const PADX = 46;     // left/right padding for the aa axis

// Paint for every scientific mark of this view, written onto the elements as
// explicit SVG attributes. The component stylesheet keeps layout, cursor,
// transition and hover; colour and stroke come from the shared specification, so
// a mark stays legible in a standalone SVG, where that stylesheet is absent.
const AXIS_TEXT = textProps("axis");
const AXIS_END_TEXT = textProps("axisEmphasis");
const BLOCK_TEXT = textProps("featureLabel");
const CAND_TEXT = textProps("candidateLabel");
const CAND_FAINT = featureProps("candidate_region", { faint: true });
const CAND_SELECTED = featureProps("candidate_region", { selected: true });
const MISSING_SEGMENT = featureProps("gap");
const ALT_TERMINUS = featureProps("alternative_exon");

function niceStep(span) {
  const raw = span / 9;
  const pow = Math.pow(10, Math.floor(Math.log10(Math.max(1, raw))));
  const n = raw / pow;
  const m = n >= 5 ? 5 : n >= 2 ? 2 : 1;
  return Math.max(5, m * pow);
}

export default function ExonMap({ model, species }) {
  const selection = useScientificSelection();
  const svgRef = useRef(null);
  const [cmpMode, setCmpMode] = useState("all");

  const speciesModel = useMemo(() => {
    const models = model?.models || [];
    if (!models.length) return null;
    return models.find((m) => (m.species_id || m.species) === (species || null)) || models[0];
  }, [model, species]);

  const length = speciesModel?.protein_length || 1;
  const exons = speciesModel?.exons || [];
  const candidates = speciesModel?.candidate_regions || [];
  const transcriptModels = speciesModel?.transcript_models || [];

  // ---- view window (zoom / pan / fit) ---- //
  const [view, setView] = useState([1, length]);
  const [lo, hi] = view;
  const span = Math.max(1, hi - lo);
  const x = (aa) => PADX + ((Number(aa) - lo) / span) * (W - 2 * PADX);
  const clampView = ([a, b]) => {
    let s = Math.max(1, Math.round(a));
    let e = Math.min(length, Math.round(b));
    if (e - s < 8) { e = Math.min(length, s + 8); s = Math.max(1, e - 8); }
    return [s, e];
  };
  const zoomBy = (factor) => {
    const c = (lo + hi) / 2;
    const half = (span / factor) / 2;
    setView(clampView([c - half, c + half]));
  };
  const pan = (frac) => {
    const d = span * frac;
    if (lo + d < 1) setView(clampView([1, 1 + span]));
    else if (hi + d > length) setView(clampView([length - span, length]));
    else setView(clampView([lo + d, hi + d]));
  };
  const fitProtein = () => setView([1, length]);
  const fitSelection = () => {
    const c = candidates.find((cc) => cc.id === selCandId);
    if (c) setView(clampView([c.start - span * 0.15, c.end + span * 0.15]));
  };

  // ---- track visibility (also captured in export state) ---- //
  const [tracks, setTracks] = useState({ exons: true, boundaries: true, candidates: true });
  const toggle = (k) => setTracks((t) => {
    const next = { ...t, [k]: !t[k] };
    selection?.setVisibleTracks?.(next);
    return next;
  });

  // ---- candidate selection (linked via alignment region) ---- //
  const [selCandId, setSelCandId] = useState(candidates[0]?.id || null);
  const selectCand = (c) => {
    setSelCandId(c.id);
    selection?.selectAlignmentRegion?.(c.start, c.end);
    selection?.setCoordinateRange?.(c.start, c.end);
  };

  const [tip, setTip] = useState(null); // {x, y, exon}
  const showTip = (ex, evt) => {
    const host = svgRef.current?.parentElement;
    if (!host) return;
    const r = host.getBoundingClientRect();
    setTip({ x: evt.clientX - r.left + 12, y: evt.clientY - r.top + 12, exon: ex });
  };

  const selExonId = selection?.selectedExonId;
  const clickExon = (ex) => {
    selection?.selectExon?.({
      exon_id: ex.id, transcript_id: ex.tooltip?.transcript_id,
      protein_id: speciesModel.protein_id,
      protein_start_aa: ex.start, protein_end_aa: ex.end,
    });
  };

  if (!speciesModel) {
    return <Empty title="Exon map not available"
      hint="No validated protein-coordinate model was built for this run." />;
  }

  // axis ticks (major + minor)
  const major = niceStep(span);
  const minor = major / 5;
  const majors = [], minors = [];
  for (let t = Math.ceil(lo / minor) * minor; t <= hi; t += minor) {
    if (Math.abs(t / major - Math.round(t / major)) < 1e-6) majors.push(Math.round(t));
    else minors.push(Math.round(t));
  }

  // Publication output is built from the coordinate model, not screenshotted from
  // the interactive SVG above: the same specification yields the standalone SVG,
  // a true vector PDF and a 300 dpi PNG. Built on click so the export reflects the
  // current candidate selection.
  const stem = `exon_map_${speciesModel.protein_id}`;
  const buildFigure = () => exonMapFigure(speciesModel, { selectedCandidateId: selCandId });
  const exportSvg = () => downloadFigureSvg(buildFigure(), stem);
  const exportPdf = () => downloadFigurePdf(buildFigure(), stem);
  const exportPng = () => downloadFigurePng(buildFigure(), stem);
  const exportTsv = () => downloadFigureTsv(exonMapTsv(speciesModel), stem);

  const tsvHref = speciesModel && model?.models
    ? fileUrl(speciesModel.exons?.[0]?.source_file
      || model.provenance?.source_files?.exon_protein_map) : null;

  return (
    <div className="viewer coord-viewer exon-map">
      <div className="viewer-head">
        <div>
          <b>Exon map</b>
          <p className="muted sm">Validated protein-coordinate model · {speciesModel.protein_id}
            {speciesModel.transcript_id ? ` · ${speciesModel.transcript_id}` : ""} · {length} aa ·
            {" "}{exons.length} coding exons</p>
        </div>
      </div>

      {/* zoom / navigation toolbar */}
      <div className="em-toolbar">
        <div className="seg">
          <button className="seg-btn" onClick={() => zoomBy(1.6)} title="Zoom in">＋</button>
          <button className="seg-btn" onClick={() => zoomBy(1 / 1.6)} title="Zoom out">－</button>
          <button className="seg-btn" onClick={fitProtein} title="Fit whole protein">Fit</button>
          <button className="seg-btn" onClick={() => setView([1, length])} title="Reset">Reset</button>
        </div>
        <div className="seg">
          <button className="seg-btn" onClick={() => pan(-0.25)} title="Pan left">◀</button>
          <button className="seg-btn" onClick={() => pan(0.25)} title="Pan right">▶</button>
          {selCandId && <button className="seg-btn" onClick={fitSelection}
            title="Zoom to selected candidate">Zoom to {selCandId}</button>}
        </div>
        <span className="muted small">visible aa {lo}–{hi}</span>
        <span className="spacer" />
        <div className="em-tracks">
          {[["exons", "exons"], ["boundaries", "boundaries"], ["candidates", "candidates"]].map(([k, l]) => (
            <label key={k} className="check inline">
              <input type="checkbox" checked={tracks[k]} onChange={() => toggle(k)} /><span>{l}</span>
            </label>
          ))}
        </div>
      </div>

      <div className="em-canvas">
      <svg ref={svgRef} className="em-svg" viewBox={`0 0 ${W} 200`} role="img"
        aria-label={`${speciesModel.protein_id} exon-to-protein map`}
        preserveAspectRatio="xMidYMid meet">
        <rect x="0" y="0" width={W} height="200" fill="#ffffff" />
        {/* Domain / TM layers deliberately live in Domain Architecture only (Part 7):
            the Exon Map answers which coding exons produce which protein regions. */}

        {/* candidate overlays (faint others + strong selected) */}
        {tracks.candidates && candidates.filter((c) => c.id !== selCandId).map((c) => (
          <rect key={c.id} x={x(c.start)} y="58" width={Math.max(2, x(c.end) - x(c.start))} height="58"
            fill={CAND_FAINT.fill} fillOpacity={CAND_FAINT.fillOpacity} stroke="none"
            style={{ cursor: "pointer" }} onClick={() => selectCand(c)}>
            <title>{`${c.label} · ${c.candidate_type || "candidate"} · ${c.confidence || ""} (exploratory)`}</title>
          </rect>
        ))}
        {tracks.candidates && candidates.filter((c) => c.id === selCandId).map((c) => (
          <g key={c.id}>
            <rect x={x(c.start)} y="54" width={Math.max(2, x(c.end) - x(c.start))} height="66" rx="3"
              fill={CAND_SELECTED.fill} fillOpacity={CAND_SELECTED.fillOpacity}
              stroke={CAND_SELECTED.stroke} strokeWidth={CAND_SELECTED.strokeWidth}
              style={{ cursor: "pointer" }} onClick={() => selectCand(c)}>
              <title>{`${c.label} · ${c.candidate_type || "candidate"} · ${c.confidence || ""} (selected, exploratory)`}</title>
            </rect>
            <text x={x(c.start)} y="50" className="em-cand-lbl" fill={CAND_TEXT.fill}
              fontSize={CAND_TEXT.fontSize} fontWeight={CAND_TEXT.fontWeight}>{c.id} · {c.start}–{c.end}
              {c.candidate_type ? ` · ${c.candidate_type}` : ""}</text>
          </g>
        ))}

        {/* protein backbone */}
        <line x1={x(Math.max(1, lo))} x2={x(Math.min(length, hi))} y1="88" y2="88"
          stroke={featureStyle("protein_backbone").fill} strokeWidth="2" />

        {/* exon boundary markers */}
        {tracks.boundaries && exons.map((ex) => (
          ex.start > 1 ? <line key={`b${ex.id}`} x1={x(ex.start)} x2={x(ex.start)} y1="66" y2="110"
            stroke="#8b98a8" strokeWidth="0.7" /> : null
        ))}

        {/* numbered coding exon blocks */}
        {tracks.exons && exons.map((ex) => {
          const bx = x(ex.start), bw = Math.max(1.5, x(ex.end) - x(ex.start));
          const sel = ex.id === selExonId;
          const paint = featureProps("coding_exon", { selected: sel });
          const tt = ex.tooltip || {};
          return (
            <g key={ex.id}>
              <rect x={bx} y="72" width={bw} height="32" rx="2.5"
                className={`em-exon${sel ? " sel" : ""}`}
                fill={paint.fill} fillOpacity={paint.fillOpacity}
                stroke={paint.stroke} strokeWidth={paint.strokeWidth}
                onClick={(e) => { clickExon(ex); showTip(ex, e); }} style={{ cursor: "pointer" }}
                onMouseMove={(e) => showTip(ex, e)}
                onMouseLeave={() => setTip(null)}>
                <title>{`${ex.label} · exon ${tt.exon_number ?? "?"}\n`
                  + `exon/CDS id ${ex.id}\n`
                  + `genomic ${tt.genomic_start ?? "—"}–${tt.genomic_end ?? "—"} · strand ${tt.strand || "?"}\n`
                  + `CDS ${tt.cds_start ?? "—"}–${tt.cds_end ?? "—"} · phase ${tt.phase ?? "n/a"}\n`
                  + `protein aa ${ex.start}–${ex.end}\n`
                  + `source ${ex.source || "—"} · shared group ${tt.shared_exon_group || "—"}`}</title>
              </rect>
              {bw >= 20 && <text x={bx + bw / 2} y="92" textAnchor="middle" className="em-exon-num"
                fill={BLOCK_TEXT.fill} fontSize={BLOCK_TEXT.fontSize} fontWeight={BLOCK_TEXT.fontWeight}
                style={{ pointerEvents: "none" }}>{ex.label}</text>}
            </g>
          );
        })}

        {/* amino-acid axis: minor + major ticks, start/end coordinates */}
        <line x1={PADX} x2={W - PADX} y1="150" y2="150" stroke={AXIS_TEXT.fill} strokeWidth="1" />
        {minors.filter((t) => t >= lo && t <= hi).map((t) => (
          <line key={`mn${t}`} x1={x(t)} x2={x(t)} y1="147" y2="153" stroke="#8b98a8" strokeWidth="0.6" />
        ))}
        {majors.filter((t) => t >= lo && t <= hi).map((t) => (
          <g key={`mj${t}`}>
            <line x1={x(t)} x2={x(t)} y1="144" y2="156" stroke={AXIS_TEXT.fill} strokeWidth="1" />
            <text x={x(t)} y="170" textAnchor="middle" className="em-axis-lbl"
              fill={AXIS_TEXT.fill} fontSize={AXIS_TEXT.fontSize}>{t}</text>
          </g>
        ))}
        <text x={PADX} y="188" textAnchor="start" className="em-axis-end" fill={AXIS_END_TEXT.fill}
          fontSize={AXIS_END_TEXT.fontSize} fontWeight={AXIS_END_TEXT.fontWeight}>1</text>
        <text x={W - PADX} y="188" textAnchor="end" className="em-axis-end" fill={AXIS_END_TEXT.fill}
          fontSize={AXIS_END_TEXT.fontSize} fontWeight={AXIS_END_TEXT.fontWeight}>{length} aa</text>
      </svg>
        {tip && (() => {
          const tt = tip.exon.tooltip || {};
          return (
            <div className="em-tip" style={{ left: tip.x, top: tip.y }}>
              <div className="em-tip-h">{tip.exon.label} · exon {tt.exon_number ?? "?"}</div>
              <div><span>exon/CDS id</span><code>{tip.exon.id}</code></div>
              <div><span>genomic</span>{tt.genomic_start ?? "—"}–{tt.genomic_end ?? "—"}</div>
              <div><span>CDS</span>{tt.cds_start ?? "—"}–{tt.cds_end ?? "—"}</div>
              <div><span>phase / strand</span>{tt.phase ?? "n/a"} / {tt.strand || "?"}</div>
              <div><span>protein aa</span>{tip.exon.start}–{tip.exon.end}</div>
              <div><span>source</span>{tip.exon.source || "—"}</div>
              <div><span>shared group</span><code>{tt.shared_exon_group || "—"}</code></div>
            </div>
          );
        })()}
      </div>

      <div className="legend res-legend">
        <span className="legend-item"><span className="pa-swatch exon" />coding exon (E1…En)</span>
        <span className="legend-item"><span className="pa-swatch boundary" />exon boundary</span>
        <span className="legend-item"><span className="pa-swatch cand" />selected candidate</span>
        <span className="legend-item"><span className="pa-swatch candfaint" />other candidates</span>
      </div>

      {/* export bar — one compact menu (Part 9) */}
      <div className="em-export">
        <Menu label="Export" title="Export figure and source table" align="right">
          <button className="menu-item" onClick={exportSvg}>Main figure — SVG (vector)</button>
          <button className="menu-item" onClick={exportPdf}>Main figure — PDF (vector)</button>
          <button className="menu-item" onClick={exportPng}>Main figure — PNG (300 dpi)</button>
          <div className="menu-sep" />
          <button className="menu-item" onClick={exportTsv}>Figure source table (TSV)</button>
          {tsvHref && <a className="menu-item" href={tsvHref}>Full exon map (TSV)</a>}
        </Menu>
      </div>

      {transcriptModels.length > 1 && (
        <CompareTranscriptModels models={transcriptModels} selection={selection}
          selectedCandidate={candidates.find((c) => c.id === selCandId)}
          mode={cmpMode} setMode={setCmpMode} />
      )}
    </div>
  );
}

// ---- inline "Compare transcript models" (coordinate-model driven) ---- //
const CW = 1000, CPAD = 8;

function groupKey(b) { return b.shared_exon_group_id || b.id || `${b.start}-${b.end}`; }
function overlaps(b, c) {
  if (!c || b.start == null || b.end == null) return false;
  return b.start <= c.end && b.end >= c.start;
}

// Exon identity is decided on genomic CDS evidence, never on protein coordinates:
// an insertion or deletion in an upstream exon shifts every downstream protein
// position, which must not mark the unchanged downstream exons as different.
// Returns "shared" | "shift" | "alt".
function classifyBlock(b, primaryByGroup) {
  const p = primaryByGroup.get(groupKey(b));
  if (!p) return "alt";                                  // exon absent from primary
  const sameSpan = p.genomic_start === b.genomic_start && p.genomic_end === b.genomic_end;
  const sameStrand = (p.strand || "") === (b.strand || "");
  const samePhase = String(p.phase ?? "") === String(b.phase ?? "");
  const sameLength = (b.end - b.start) === (p.end - p.start);
  return (sameSpan && sameStrand && samePhase && sameLength) ? "shared" : "shift";
}

const BLOCK_CLS = { shared: "cmp-blk-shared", shift: "cmp-blk-shift", alt: "cmp-blk-alt" };
// The difference category is a scientific statement, so its paint comes from the
// shared specification rather than from the row's stylesheet class.
const BLOCK_KEY = { shared: "shared_exon", shift: "shifted_boundary", alt: "alternative_exon" };
const BLOCK_WHY = {
  shared: "shared genomic exon",
  shift: "shared exon with shifted boundary",
  alt: "alternative exon / inserted CDS",
};

function CompareTranscriptModels({ models, selection, selectedCandidate, mode, setMode }) {
  const [hoverGroup, setHoverGroup] = useState(null);

  const ordered = useMemo(
    () => [...models].sort((a, b) => (b.is_primary ? 1 : 0) - (a.is_primary ? 1 : 0)), [models]);
  const primary = ordered.find((m) => m.is_primary) || ordered[0];
  const maxLen = Math.max(...ordered.map((m) => m.protein_length || 0), 1);
  const total = ordered.length;

  const groupCounts = useMemo(() => {
    const counts = new Map();
    for (const m of ordered) {
      const seen = new Set();
      for (const b of m.blocks || []) {
        const k = groupKey(b);
        if (!seen.has(k)) { seen.add(k); counts.set(k, (counts.get(k) || 0) + 1); }
      }
    }
    return counts;
  }, [ordered]);

  // Primary exon structure indexed by shared exon group (genomic identity).
  const primaryByGroup = useMemo(() => {
    const map = new Map();
    for (const b of primary?.blocks || []) map.set(groupKey(b), b);
    return map;
  }, [primary]);

  // Protein intervals of the primary that this model does not represent at all.
  // Drawn at primary coordinates as a dashed empty interval.
  const missingOf = (m) => {
    if (m.is_primary) return [];
    const have = new Set((m.blocks || []).map(groupKey));
    return (primary?.blocks || []).filter((b) => !have.has(groupKey(b)));
  };
  // Alternative N-/C-terminus: the model starts or ends on a different exon.
  const terminiOf = (m) => {
    const bl = m.blocks || [], pb = primary?.blocks || [];
    if (m.is_primary || !bl.length || !pb.length) return { altN: false, altC: false };
    return {
      altN: groupKey(bl[0]) !== groupKey(pb[0]),
      altC: groupKey(bl[bl.length - 1]) !== groupKey(pb[pb.length - 1]),
    };
  };
  const modelDiffers = (m) => {
    if (m.is_primary) return false;
    if (m.protein_length !== primary?.protein_length) return true;
    if (missingOf(m).length) return true;
    return (m.blocks || []).some((b) => classifyBlock(b, primaryByGroup) !== "shared");
  };

  const identicalCount = ordered.filter((m) => !m.is_primary && !modelDiffers(m)).length;

  const rows = ordered.filter((m) => {
    if (mode === "diff") return m.is_primary || modelDiffers(m);
    if (mode === "curated") return m.curation_status === "curated";
    return true;
  });

  const x = (aa) => CPAD + (Math.max(0, aa) / maxLen) * (CW - 2 * CPAD);
  const ticks = [];
  const step = niceStep(maxLen);
  for (let t = 0; t <= maxLen; t += step) ticks.push(t);
  const minorTicks = [];
  for (let t = step / 2; t <= maxLen; t += step) minorTicks.push(t);

  const MODES = [["all", "All models"], ["diff", "Differences only"], ["curated", "Curated only"]];

  return (
    <details className="tech-prov compare-models" open>
      <summary>Compare transcript models ({models.length})</summary>
      <div className="cmp-controls">
        <div className="seg">
          {MODES.map(([id, l]) => (
            <button key={id} className={`seg-btn${mode === id ? " on" : ""}`}
              onClick={() => setMode(id)}>{l}</button>
          ))}
        </div>
        <span className="muted small">Common aa scale · 1–{maxLen} aa · primary first
          {mode === "diff"
            ? ` · ${rows.length} of ${total} models shown${identicalCount ? ` · ${identicalCount} model${identicalCount > 1 ? "s" : ""} identical to primary` : ""}`
            : mode !== "all" ? ` · ${rows.length}/${total} models` : ""}</span>
      </div>

      <div className="legend cmp-legend">
        <span className="legend-item"><span className="cmp-sw shared" />shared genomic exon</span>
        <span className="legend-item"><span className="cmp-sw alt" />alternative exon / inserted CDS</span>
        <span className="legend-item"><span className="cmp-sw shift" />shifted exon boundary</span>
        <span className="legend-item"><span className="cmp-sw missing" />missing protein segment</span>
        <span className="legend-item"><span className="cmp-sw term" />alternative N-/C-terminus</span>
        <span className="legend-item"><span className="pa-swatch cand" />selected candidate</span>
      </div>

      <div className="cmp-row cmp-ruler-row">
        <div className="cmp-label cmp-ruler-label"><span className="muted small">aa</span></div>
        <svg className="cmp-ruler" viewBox={`0 0 ${CW} 18`} preserveAspectRatio="none">
          <line x1={x(0)} y1="4" x2={x(maxLen)} y2="4" stroke="#d2d9e6" strokeWidth="0.6" />
          {minorTicks.map((t) => (
            <line key={`m${t}`} x1={x(t)} y1="2" x2={x(t)} y2="5.5"
              stroke="#d2d9e6" strokeWidth="0.4" opacity="0.6" />
          ))}
          {ticks.map((t) => (
            <g key={t}>
              <line x1={x(t)} y1="1" x2={x(t)} y2="8" stroke="#d2d9e6" strokeWidth="0.7" />
              <text x={x(t)} y="16" textAnchor={t === 0 ? "start" : "middle"} className="cmp-tick"
                fill={AXIS_TEXT.fill} fontSize={AXIS_TEXT.fontSize}>
                {t === 0 ? 1 : t}</text>
            </g>
          ))}
          <line x1={x(maxLen)} y1="1" x2={x(maxLen)} y2="8" stroke="#d2d9e6" strokeWidth="0.7" />
          <text x={x(maxLen)} y="16" textAnchor="end" className="cmp-tick"
            fill={AXIS_TEXT.fill} fontSize={AXIS_TEXT.fontSize}>{maxLen} aa</text>
        </svg>
      </div>

      <div className="cmp-rows">
        {rows.map((m) => {
          const lenDiff = primary && m.protein_length != null && primary.protein_length != null
            && m.protein_length !== primary.protein_length
            ? `${m.protein_length > primary.protein_length ? "+" : ""}${m.protein_length - primary.protein_length} aa` : null;
          return (
            <div key={m.protein_id} className={`cmp-row${m.is_primary ? " is-primary" : ""}`}
              onClick={() => selection?.selectProtein?.(m.protein_id)} style={{ cursor: "pointer" }}>
              <div className="cmp-label">
                <b className="cmp-pid">{m.protein_id}{m.is_primary ? " ★" : ""}</b>
                <span className="cmp-tx">{m.transcript_id || "—"}</span>
                <span className="cmp-dim">{m.exon_count} exons • {m.protein_length ?? "—"} aa
                  {lenDiff ? ` (Δ ${lenDiff})` : ""}</span>
                <span className="cmp-tags">
                  <Badge cls={m.curation_status === "curated" ? "accepted" : "neutral"} soft>
                    {m.curation_status || "predicted"}</Badge>
                  <Badge cls={m.is_primary ? "accepted" : "neutral"} soft>
                    {m.is_primary ? "primary" : "alternative"}</Badge>
                </span>
              </div>
              <svg className="cmp-track" viewBox={`0 0 ${CW} 26`} preserveAspectRatio="none">
                <line x1={x(0)} y1="12" x2={x(m.protein_length || maxLen)} y2="12"
                  stroke="#d2d9e6" strokeWidth="0.6" />
                {selectedCandidate && (() => {
                  const rowOv = (m.blocks || []).some((b) => overlaps(b, selectedCandidate));
                  // The band marks the same candidate on every row; rows it does not
                  // touch keep it faint instead of dropping the context.
                  const band = featureProps("candidate_region", { faint: !rowOv });
                  return (
                    <rect x={x(selectedCandidate.start)} y="2"
                      width={Math.max(2, x(selectedCandidate.end) - x(selectedCandidate.start))} height="22"
                      rx="2" className={`cmp-cand${rowOv ? " overlaps" : ""}`}
                      fill={band.fill} fillOpacity={band.fillOpacity}
                      stroke={band.stroke} strokeWidth={band.strokeWidth}>
                      <title>{`${selectedCandidate.id} · aa ${selectedCandidate.start}–${selectedCandidate.end}`
                        + `${rowOv ? " · overlaps this model" : ""}`}</title>
                    </rect>
                  );
                })()}
                {/* protein regions of the primary that this model does not represent */}
                {missingOf(m).map((b) => (
                  <rect key={`miss-${groupKey(b)}`} x={x(b.start)} y="5"
                    width={Math.max(1.5, x(b.end) - x(b.start))} height="14" rx="1.5"
                    className="cmp-blk-missing" fill={MISSING_SEGMENT.fill}
                    stroke={MISSING_SEGMENT.stroke} strokeWidth={MISSING_SEGMENT.strokeWidth}>
                    <title>{`missing protein segment · ${b.label || b.id}\n`
                      + `absent from ${m.protein_id}; drawn at primary aa ${b.start}–${b.end}`}</title>
                  </rect>
                ))}
                {(m.blocks || []).map((b, i) => {
                  if (b.start == null || b.end == null) return null;
                  const bx = x(b.start), bw = Math.max(1.2, x(b.end) - x(b.start));
                  const gk = groupKey(b);
                  const shared = groupCounts.get(gk) || 1;
                  const isHl = hoverGroup && gk === hoverGroup;
                  const kind = m.is_primary ? "shared" : classifyBlock(b, primaryByGroup);
                  const candOv = overlaps(b, selectedCandidate);
                  // In differences-only mode unchanged exons stay readable but recede.
                  const muted = mode === "diff" && kind === "shared" && !m.is_primary && !candOv;
                  const exonNo = b.transcript_exon_number ?? b.exon_number;
                  // Exon numbers appear inside sufficiently wide blocks; narrow blocks
                  // keep the number in the hover tooltip only.
                  const showLabel = bw >= 16 && exonNo != null;
                  const paint = featureProps(BLOCK_KEY[kind]);
                  return (
                    <g key={i}
                      onMouseEnter={() => setHoverGroup(gk)}
                      onMouseLeave={() => setHoverGroup((c) => (c === gk ? null : c))}
                      onClick={(e) => { e.stopPropagation();
                        selection?.selectExon?.({ exon_id: b.id, transcript_id: b.transcript_id,
                          protein_id: m.protein_id, protein_start_aa: b.start, protein_end_aa: b.end }); }}>
                      <rect x={bx} y="4" width={bw} height="16" rx="1.5"
                        className={`cmp-blk ${BLOCK_CLS[kind]}${isHl ? " hl" : ""}`}
                        fill={paint.fill} fillOpacity={paint.fillOpacity}
                        stroke={paint.stroke} strokeWidth={paint.strokeWidth}
                        opacity={muted ? 0.45 : 1}>
                        <title>{`${b.label || b.id} · exon ${exonNo ?? "?"} · ${BLOCK_WHY[kind]}\n`
                          + `protein aa ${b.start}–${b.end}\n`
                          + `genomic ${b.genomic_start ?? "—"}–${b.genomic_end ?? "—"}`
                          + ` · strand ${b.strand || "n/a"} · phase ${b.phase ?? "n/a"}\n`
                          + `CDS ${b.cds_start ?? "—"}–${b.cds_end ?? "—"}\n`
                          + `shared exon group present in ${shared}/${total} models`}</title>
                      </rect>
                      {showLabel && (
                        <text x={bx + bw / 2} y="15" textAnchor="middle" className="cmp-exon-lbl"
                          fill={BLOCK_TEXT.fill} fontSize={BLOCK_TEXT.fontSize}
                          fontWeight={BLOCK_TEXT.fontWeight}
                          opacity={muted ? 0.5 : 1}>E{exonNo}</text>
                      )}
                    </g>
                  );
                })}
                {/* alternative N-/C-terminus markers */}
                {(() => {
                  const { altN, altC } = terminiOf(m);
                  const bl = m.blocks || [];
                  if (!bl.length) return null;
                  return (
                    <>
                      {altN && (
                        <rect x={Math.max(0, x(bl[0].start) - 3)} y="2" width="3" height="20"
                          className="cmp-term" fill={ALT_TERMINUS.fill}
                          stroke={ALT_TERMINUS.stroke} strokeWidth={ALT_TERMINUS.strokeWidth}>
                          <title>{`alternative N-terminus · starts on ${bl[0].label || bl[0].id}`}</title>
                        </rect>
                      )}
                      {altC && (
                        <rect x={x(bl[bl.length - 1].end)} y="2" width="3" height="20"
                          className="cmp-term" fill={ALT_TERMINUS.fill}
                          stroke={ALT_TERMINUS.stroke} strokeWidth={ALT_TERMINUS.strokeWidth}>
                          <title>{`alternative C-terminus · ends on ${bl[bl.length - 1].label || bl[bl.length - 1].id}`}</title>
                        </rect>
                      )}
                    </>
                  );
                })()}
              </svg>
              <span className="cmp-len">{m.protein_length ?? "—"} aa</span>
            </div>
          );
        })}
      </div>
      <p className="muted sm">Exon identity comes from genomic CDS coordinates, so an upstream
        insertion or deletion does not mark the unchanged downstream exons as different.
        Protein-sequence differences are shown in the isoform alignment under Exploratory
        Candidate Evidence.</p>
    </details>
  );
}
