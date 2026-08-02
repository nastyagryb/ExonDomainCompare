import { useMemo, useRef, useState } from "react";
import { Badge, Empty, Drawer, Field, Menu } from "../../ui";
import { useScientificSelection } from "../../components/ScientificSelectionContext";
import {
  candidateDisplayLayout, candidateLabelFits, candidateTooltip, laneDensity,
} from "./candidateDisplay";
import { domainInstances } from "./common";
import { domainArchitectureFigure, domainArchitectureTsv } from "./figureData";
import {
  downloadFigurePdf, downloadFigurePng, downloadFigureSvg, downloadFigureTsv,
} from "./figureExport";
import {
  domainInstanceFill, featureProps, featureStyle, textProps,
} from "./semanticStyles";

// Interactive Domain Architecture driven ENTIRELY by the validated
// protein-coordinate model (scripts/shared_gene_analysis/protein_coordinate_model.py),
// the exact same single source of truth used by the Exon Map. No coordinate is
// reconstructed in React and no biological feature is fabricated: pre-cluster
// runs render the same axis + exon/candidate context and report the domain /
// family / site / TM layers as explicitly *pending*.

const W = 1000;
// The left gutter is reserved exclusively for track labels so a label can never
// overlap a feature block; the plot area starts after it.
const LBLW = 136;
const PADR = 46;
// Domain instances, families and TM helices are named by the shared scientific
// specification, so their paint is read from it instead of being restated here.
// Member signatures and disorder regions have no key in either module yet and keep
// their local fill until one names them.
const FAMILY_FILL = featureStyle("family_superfamily").fill;
const SIG_FILL = featureStyle("member_signature").fill;
const SITE_FILL = featureStyle("functional_site").fill;
const DISORDER_FILL = featureStyle("disorder_region").fill;
const TM_FILL = featureStyle("tm_helix").fill;

// Paint written onto the marks as explicit SVG attributes. Component CSS keeps
// layout, cursor, transition and hover; colour and stroke come from the shared
// specification, so a mark stays legible without that stylesheet.
const AXIS_TEXT = textProps("axis");
const AXIS_END_TEXT = textProps("axisEmphasis");
const LANE_TEXT = textProps("trackLabel");
const ON_BLOCK_TEXT = textProps("onFeatureLabel");
const BLOCK_TEXT = textProps("featureLabel");
const CAND_TEXT = textProps("candidateLabel");
const EMPTY_TEXT = textProps("empty");
const DOMAIN_PAINT = featureStyle("representative_domain");
const TM_PAINT = featureProps("tm_helix");
const SELECTION = featureStyle("selected_feature");

function niceStep(span) {
  const raw = span / 9;
  const pow = Math.pow(10, Math.floor(Math.log10(Math.max(1, raw))));
  const n = raw / pow;
  const m = n >= 5 ? 5 : n >= 2 ? 2 : 1;
  return Math.max(5, m * pow);
}
function overlaps(a, b) {
  return a && b && a.start != null && b.start != null && a.start <= b.end && a.end >= b.start;
}
function downloadBlob(blob, name) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export default function ProteinArchitecture({ model, species }) {
  const selection = useScientificSelection();
  const svgRef = useRef(null);
  const [detail, setDetail] = useState(null);

  const speciesModel = useMemo(() => {
    const models = model?.models || [];
    if (!models.length) return null;
    return models.find((m) => (m.species_id || m.species) === (species || null)) || models[0];
  }, [model, species]);

  const length = speciesModel?.protein_length || 1;
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

  // Default view: representative domains + families + coding exons + boundaries
  // + TM + candidates visible; member signatures collapsed; sites/disorder off.
  const [vis, setVis] = useState({
    families: true, signatures: false, sites: false, disorder: false,
    tm: true, exons: true, boundaries: true, candidates: true,
  });
  // Overlapping candidate clusters are packed into lanes; by default only the
  // top-ranked lanes are drawn so the track stays readable, and this opens the rest.
  const [showAllCandidates, setShowAllCandidates] = useState(false);
  const setTrack = (k, v) => setVis((t) => {
    const next = { ...t, [k]: v ?? !t[k] };
    selection?.setVisibleTracks?.(next);
    return next;
  });
  const representativeOnly = () => {
    const next = { ...vis, families: false, signatures: false, sites: false, disorder: false };
    setVis(next); selection?.setVisibleTracks?.(next);
  };

  const domains = useMemo(
    () => domainInstances(speciesModel?.representative_domains || []), [speciesModel]);
  const families = speciesModel?.families_superfamilies || [];
  const signatures = speciesModel?.member_signatures || [];
  const sites = speciesModel?.functional_sites || [];
  const disorder = speciesModel?.disorder_regions || [];
  const tms = speciesModel?.tm_regions || [];
  const exons = useMemo(() => speciesModel?.exons || [], [speciesModel]);
  const candidates = useMemo(() => speciesModel?.candidate_regions || [], [speciesModel]);
  const tmInfo = speciesModel?.tm_analysis || {};
  const pending = speciesModel?.status !== "available";
  const pendingInfo = speciesModel?.pending_info || null;
  // Pre-cluster: never render multiple empty domain-layer rows with repeated
  // "pending" labels. Only exon / boundary / candidate context is real, so the
  // domain-family-site-disorder-TM lanes collapse to a single pending card above.
  const eff = pending
    ? { ...vis, families: false, signatures: false, sites: false, disorder: false, tm: false }
    : vis;

  // Linked selection (persisted across tabs via ScientificSelectionContext).
  const selDomainId = selection?.selectedDomainId || null;
  // Repeated domains share an InterPro accession, so the exported figure is told
  // which *instance* is selected rather than which accession.
  const selDomainInstanceId = useMemo(() => {
    const hit = (speciesModel?.representative_domains || [])
      .find((d) => d.id === selDomainId);
    return hit?.domain_instance_id || null;
  }, [speciesModel, selDomainId]);
  // Optional "all instances of one InterPro entry" selection (e.g. all Ig-like domains).
  const [domGroup, setDomGroup] = useState(null);
  const domainGroups = useMemo(() => {
    const by = new Map();
    for (const d of domains) {
      if (!by.has(d.groupKey)) by.set(d.groupKey, { key: d.groupKey, label: d.groupLabel, n: 0 });
      by.get(d.groupKey).n += 1;
    }
    return [...by.values()].filter((g) => g.n > 1);
  }, [domains]);
  const selExonId = selection?.selectedExonId || null;
  const alignS = selection?.selectedAlignmentStart;
  const alignE = selection?.selectedAlignmentEnd;
  const selCandidate = useMemo(() => {
    if (alignS == null || alignE == null) return candidates[0] || null;
    const reg = { start: alignS, end: alignE };
    return candidates.find((c) => overlaps(c, reg)) || null;
  }, [candidates, alignS, alignE]);

  const boundaryAAs = useMemo(
    () => exons.map((e) => e.start).filter((v) => v != null && v > 1), [exons]);

  const openDomainDetail = (d) => {
    selection?.selectDomain?.(d);
    const contribSigs = signatures.filter(
      (s) => (d.interpro_accession && s.interpro_accession === d.interpro_accession) || overlaps(s, d));
    const ovExons = exons.filter((e) => overlaps(e, d));
    const internalB = boundaryAAs.filter((b) => b > d.start && b < d.end);
    const nearestCand = candidates.map((c) => ({
      c, dist: overlaps(c, d) ? 0 : Math.min(Math.abs(c.start - d.end), Math.abs(d.start - c.end)),
    })).sort((a, b) => a.dist - b.dist)[0];
    setDetail({ kind: "domain", feature: d, contribSigs, ovExons, internalB, nearestCand });
  };
  const openFeatureDetail = (f, layer) => {
    selection?.selectFeature?.(f);
    const ovExons = exons.filter((e) => overlaps(e, f));
    setDetail({ kind: "feature", layer, feature: f, ovExons });
  };
  const clickExon = (ex) => {
    selection?.selectExon?.({
      exon_id: ex.id, transcript_id: ex.tooltip?.transcript_id,
      protein_id: speciesModel.protein_id, protein_start_aa: ex.start, protein_end_aa: ex.end,
    });
    const ovDomains = domains.filter((d) => overlaps(d, ex));
    const ovFeatures = [...families, ...signatures, ...sites, ...disorder, ...tms].filter((f) => overlaps(f, ex));
    const adjB = boundaryAAs.filter((b) => Math.abs(b - ex.start) < 2 || Math.abs(b - ex.end) < 2);
    setDetail({ kind: "exon", feature: ex, ovDomains, ovFeatures, adjB });
  };
  const clickCandidate = (c) => {
    selection?.selectAlignmentRegion?.(c.start, c.end);
    selection?.setCoordinateRange?.(c.start, c.end);
    const ovDomains = domains.filter((d) => overlaps(d, c));
    const nearestDomain = ovDomains[0] || domains.map((d) => ({
      d, dist: Math.min(Math.abs(d.start - c.end), Math.abs(c.start - d.end)),
    })).sort((a, b) => a.dist - b.dist)[0]?.d || null;
    setDetail({ kind: "candidate", feature: c, ovDomains, nearestDomain });
  };

  if (!speciesModel) {
    return <Empty title="Domain architecture not available"
      hint="No validated protein-coordinate model was built for this run." />;
  }

  // ---- dynamic stacked-track layout (only visible tracks take vertical space) ---- //
  const LANE_H = 20, GAP = 14, RULER_H = 38;
  const lanes = [];
  const addLane = (id, label, items) => lanes.push({ id, label, items });
  addLane("domains", "Representative domains", domains);        // always
  if (eff.families) addLane("families", "Families / superfamilies", families);
  if (eff.signatures) addLane("signatures", "Member signatures", signatures);
  if (eff.sites) addLane("sites", "Functional sites", sites);
  if (eff.disorder) addLane("disorder", "Disorder / other", disorder);
  if (eff.tm) addLane("tm", "pyTMHMM topology", tms);
  if (eff.exons) addLane("exons", "Coding exons", exons);
  if (eff.boundaries) addLane("boundaries", "Exon boundaries", boundaryAAs);
  if (eff.candidates) addLane("candidates", "Candidate regions", candidates);

  const laneY = {};
  lanes.forEach((ln, i) => { laneY[ln.id] = RULER_H + i * (LANE_H + GAP); });
  // Candidate clusters are packed into sub-lanes (shared with the exported figure via
  // the coordinate model's `display_lane`), so this track claims the height its lanes
  // need instead of stacking every box on one row.
  const candLayout = candidateDisplayLayout(candidates, {
    selectedId: selCandidate?.id || null, showAll: showAllCandidates,
  });
  const CAND_LANE_H = 13, CAND_LANE_GAP = 3;
  const candExtra = eff.candidates && candLayout.laneCount > 0
    ? Math.max(0, candLayout.laneCount * (CAND_LANE_H + CAND_LANE_GAP) - CAND_LANE_GAP - LANE_H)
    : 0;
  const candLaneY = (lane) => laneY.candidates + lane * (CAND_LANE_H + CAND_LANE_GAP);
  const H = RULER_H + lanes.length * (LANE_H + GAP) + candExtra + 22;

  // axis ticks
  const major = niceStep(span);
  const minor = major / 5;
  const majors = [], minors = [];
  for (let t = Math.ceil(lo / minor) * minor; t <= hi; t += minor) {
    if (Math.abs(t / major - Math.round(t / major)) < 1e-6) majors.push(Math.round(t));
    else minors.push(Math.round(t));
  }

  const laneColor = { families: FAMILY_FILL, signatures: SIG_FILL, sites: SITE_FILL, disorder: DISORDER_FILL, tm: TM_FILL };
  const domColorMap = {};
  domains.forEach((d, i) => { domColorMap[d.id] = domainInstanceFill(i + 1); });
  // The interactive view above is optimised for exploration; publication output is
  // rendered from the coordinate model through the shared figure specification, so
  // the exported SVG, PDF and PNG are one vector figure rather than a screenshot.
  const archStem = `domain_architecture_${speciesModel.protein_id}`;
  const buildArchFigure = () => domainArchitectureFigure(speciesModel,
    { selectedDomainInstanceId: selDomainInstanceId, showAllCandidates });
  const exportSvg = () => downloadFigureSvg(buildArchFigure(), archStem);
  const exportPdf = () => downloadFigurePdf(buildArchFigure(), archStem);
  const exportPng = () => downloadFigurePng(buildArchFigure(), archStem);
  const exportFigureTsv = () => downloadFigureTsv(
    domainArchitectureTsv(speciesModel), archStem);
  const exportFeatureTsv = () => {
    const head = ["feature_id", "label", "layer", "start_aa", "end_aa", "source",
      "interpro_accession", "status", "source_file"].join("\t");
    const rows = [];
    const push = (arr, layer) => (arr || []).forEach((f) => rows.push([
      f.id, f.label, layer, f.start, f.end, f.source || "",
      f.interpro_accession || "", f.status || "", f.source_file || ""].join("\t")));
    push(domains, "representative_domain"); push(families, "family_superfamily");
    push(signatures, "member_signature"); push(sites, "functional_site");
    push(disorder, "disorder_region"); push(tms, "tm_region");
    (exons || []).forEach((e) => rows.push([e.id, e.label, "coding_exon", e.start, e.end,
      e.source || "", "", e.status || "", e.source_file || ""].join("\t")));
    (candidates || []).forEach((c) => rows.push([c.id, c.label, "candidate_region", c.start, c.end,
      "exploratory", "", c.confidence || "", ""].join("\t")));
    downloadBlob(new Blob([[head, ...rows].join("\n")], { type: "text/tab-separated-values" }),
      `domain_features_${speciesModel.protein_id}.tsv`);
  };
  const barW = (a, b) => Math.max(2, x(b) - x(a));
  const nDom = domains.length, nFam = families.length, nSig = signatures.length;

  return (
    <div className="viewer coord-viewer exon-map domain-arch">
      <div className="viewer-head">
        <div>
          <b>Domain architecture</b>
          <p className="muted sm">{speciesModel.scientific_name} ·
            {" "}<code>{speciesModel.protein_id}</code>
            {speciesModel.transcript_id ? <> · <code>{speciesModel.transcript_id}</code></> : null}
            {" "}· {length} aa · {exons.length} coding exons</p>
        </div>
        <Badge cls={pending ? "neutral" : "accepted"} soft>
          {pending ? "domain layers pending cluster" : `${nDom} domains · ${nFam} families · ${nSig} signatures`}</Badge>
      </div>

      {/* Compact toolbar: Zoom | Fit protein | Tracks▾ | Export▾ — advanced
          track and export controls are progressively disclosed inside menus. */}
      <div className="em-toolbar compact-toolbar">
        <div className="seg">
          <button className="seg-btn" onClick={() => zoomBy(1.6)} title="Zoom in">＋</button>
          <button className="seg-btn" onClick={() => zoomBy(1 / 1.6)} title="Zoom out">－</button>
          <button className="seg-btn" onClick={() => pan(-0.25)} title="Pan left">◀</button>
          <button className="seg-btn" onClick={() => pan(0.25)} title="Pan right">▶</button>
        </div>
        <button className="seg-btn" onClick={() => setView([1, length])} title="Fit whole protein">Fit protein</button>
        {domains.length > 0 && (
          <select className="da-dom-select" title="Select a representative domain instance"
            value={domGroup ? `grp:${domGroup}` : (selDomainId || "")}
            onChange={(e) => {
              const v = e.target.value;
              if (!v) { setDomGroup(null); selection?.selectDomain?.(null); return; }
              if (v.startsWith("grp:")) {
                const g = v.slice(4);
                setDomGroup(g);
                selection?.selectDomain?.(null);
                const inst = domains.filter((d) => d.groupKey === g);
                if (inst.length) setView(clampView([
                  Math.min(...inst.map((d) => d.start)) - 10,
                  Math.max(...inst.map((d) => d.end)) + 10]));
                return;
              }
              setDomGroup(null);
              const d = domains.find((z) => z.id === v);
              if (d) { openDomainDetail(d); setView(clampView([d.start - 15, d.end + 15])); }
            }}>
            <option value="">Select domain instance…</option>
            {domains.map((d) => (
              <option key={d.id} value={d.id}>{d.instanceLabel}</option>
            ))}
            {domainGroups.map((g) => (
              <option key={`grp:${g.key}`} value={`grp:${g.key}`}>
                {g.label} ({g.n} instances)</option>
            ))}
          </select>
        )}
        {selCandidate && <button className="seg-btn"
          onClick={() => setView(clampView([selCandidate.start - span * 0.15, selCandidate.end + span * 0.15]))}
          title="Zoom to selected candidate">Zoom to {selCandidate.id}</button>}
        {eff.candidates && candLayout.hiddenCount > 0 && (
          <button className="seg-btn" onClick={() => setShowAllCandidates(true)}
            title="Show every exploratory candidate cluster, including lower-ranked lanes">
            Show all candidates ({candLayout.hiddenCount} more)</button>
        )}
        {eff.candidates && showAllCandidates && candLayout.total > candLayout.laneCount && (
          <button className="seg-btn on" onClick={() => setShowAllCandidates(false)}
            title="Show only the top-ranked candidate lanes">Top-ranked candidates only</button>
        )}
        <Menu label="Tracks" title="Show / hide feature tracks">
          <button className="menu-item" onClick={representativeOnly}>Representative domains only</button>
          <div className="menu-sep" />
          {[["families", "Families / superfamilies"], ["signatures", "Member signatures"],
            ["sites", "Functional sites"], ["disorder", "Disorder / other"], ["tm", "TM (pyTMHMM)"],
            ["exons", "Coding exons"], ["boundaries", "Exon boundaries"], ["candidates", "Candidate regions"]]
            .map(([k, l]) => (
            <label key={k} className="menu-check">
              <input type="checkbox" checked={vis[k]} onChange={() => setTrack(k)}
                disabled={pending && ["families", "signatures", "sites", "disorder", "tm"].includes(k)} />
              <span>{l}</span>
            </label>
          ))}
        </Menu>
        <Menu label="Export" title="Export figure and data" align="right">
          <button className="menu-item" onClick={exportSvg}>Main figure — SVG (vector)</button>
          <button className="menu-item" onClick={exportPdf}>Main figure — PDF (vector)</button>
          <button className="menu-item" onClick={exportPng}>Main figure — PNG (300 dpi)</button>
          <button className="menu-item" onClick={exportFigureTsv}>Figure source table (TSV)</button>
          <div className="menu-sep" />
          <button className="menu-item" onClick={exportFeatureTsv}>Feature table (TSV)</button>
        </Menu>
        <span className="spacer" />
        <span className="muted small">visible aa {lo}–{hi}</span>
      </div>

      {pending && (
        <div className="pending-note compact">
          <Badge cls="neutral" soft>pending cluster</Badge>
          <span>{pendingInfo?.message
            || "Representative domain, family, site and transmembrane layers are computed after the cluster InterProScan / pyTMHMM round-trip."}
            {" "}Coding exons, exon boundaries and exploratory candidate regions below are already
            projected on the real amino-acid axis. The cluster command is available in <b>My Runs</b>.</span>
        </div>
      )}

      <div className="em-canvas">
        <svg ref={svgRef} className="em-svg" viewBox={`0 0 ${W} ${H}`} role="img"
          aria-label={`${speciesModel.protein_id} domain architecture`} preserveAspectRatio="xMidYMid meet">
          <rect x="0" y="0" width={W} height={H} fill="#ffffff" />

          {/* 1 — shared amino-acid coordinate axis (top ruler) */}
          <line x1={LBLW} x2={W - PADR} y1={RULER_H - 8} y2={RULER_H - 8} stroke={AXIS_TEXT.fill} strokeWidth="1" />
          {minors.filter((t) => t >= lo && t <= hi).map((t) => (
            <line key={`mn${t}`} x1={x(t)} x2={x(t)} y1={RULER_H - 11} y2={RULER_H - 5} stroke={featureStyle("exon_boundary_tick").stroke} strokeWidth="0.6" />
          ))}
          {majors.filter((t) => t >= lo && t <= hi).map((t) => (
            <g key={`mj${t}`}>
              <line x1={x(t)} x2={x(t)} y1={RULER_H - 14} y2={RULER_H - 2} stroke={AXIS_TEXT.fill} strokeWidth="1" />
              <text x={x(t)} y={RULER_H - 18} textAnchor="middle" className="em-axis-lbl"
                fill={AXIS_TEXT.fill} fontSize={AXIS_TEXT.fontSize}>{t}</text>
            </g>
          ))}

          {/* Track lane labels in the reserved left gutter + guide baseline.
              Each track owns one row, so labels never sit on a feature block. */}
          {lanes.map((ln) => (
            <g key={ln.id}>
              <text x={LBLW - 10} y={laneY[ln.id] + LANE_H / 2 + 3} textAnchor="end"
                className="da-lane-lbl" fill={LANE_TEXT.fill}
                fontSize={LANE_TEXT.fontSize}>{ln.label}</text>
              <line x1={LBLW} x2={W - PADR} y1={laneY[ln.id] + LANE_H / 2}
                y2={laneY[ln.id] + LANE_H / 2} stroke="#e3e8ee" strokeWidth="0.7" />
            </g>
          ))}

          {/* selected-domain (instance or whole InterPro entry) highlight band */}
          {domains.filter((d) => (domGroup ? d.groupKey === domGroup : d.id === selDomainId)).map((d) => (
            <rect key={`hl${d.id}`} x={x(d.start)} y={RULER_H - 2} width={barW(d.start, d.end)}
              height={H - RULER_H - 20} fill={DOMAIN_PAINT.fill} opacity="0.08" />
          ))}

          {/* 2 — representative domains */}
          {domains.map((d) => {
            const sel = d.id === selDomainId;
            return (
              <g key={d.id}>
                <rect x={x(d.start)} y={laneY.domains} width={barW(d.start, d.end)} height={LANE_H} rx="3"
                  fill={domColorMap[d.id]} stroke={sel ? SELECTION.stroke : DOMAIN_PAINT.stroke}
                  strokeWidth={sel ? 1.4 : 0.4}
                  style={{ cursor: "pointer" }} onClick={() => openDomainDetail(d)}>
                  <title>{`${d.instanceLabel}\nInterPro short name: ${d.label}\n`
                    + `InterPro ${d.interpro_accession || "—"} · ${d.tooltip?.interpro_type || "DOMAIN"}\n`
                    + `aa ${d.start}–${d.end} (${d.end - d.start + 1} aa)\ninstance id: ${d.id}`}</title>
                </rect>
                {barW(d.start, d.end) > 44 && (
                  <text x={(x(d.start) + x(d.end)) / 2} y={laneY.domains + LANE_H / 2 + 3}
                    textAnchor="middle" className="da-blk-lbl" fill={ON_BLOCK_TEXT.fill}
                    fontSize={ON_BLOCK_TEXT.fontSize}
                    fontWeight={ON_BLOCK_TEXT.fontWeight}>{d.shortLabel}</text>)}
              </g>
            );
          })}
          {!domains.length && (
            <text x={W / 2} y={laneY.domains + LANE_H / 2 + 3} textAnchor="middle" className="da-empty"
              fill={EMPTY_TEXT.fill} fontSize={EMPTY_TEXT.fontSize} fontStyle={EMPTY_TEXT.fontStyle}>
              {pending ? "Representative domains pending post-cluster InterProScan"
                : "No representative InterPro domain for this protein"}</text>)}

          {/* 3–6 — families / signatures / sites / disorder lanes */}
          {["families", "signatures", "sites", "disorder"].filter((k) => eff[k]).map((k) => {
            const items = { families, signatures, sites, disorder }[k];
            return (
              <g key={k}>
                {items.map((f) => (
                  <rect key={f.id} x={x(f.start)} y={laneY[k] + (k === "families" ? 5 : 6)}
                    width={barW(f.start, f.end)} height={k === "families" ? LANE_H - 10 : LANE_H - 12} rx="2"
                    fill={laneColor[k]} opacity={selection?.selectedFeatureId === f.id ? 1 : 0.78}
                    stroke={selection?.selectedFeatureId === f.id ? SELECTION.stroke : "none"}
                    strokeWidth="0.8"
                    style={{ cursor: "pointer" }} onClick={() => openFeatureDetail(f, k)}>
                    <title>{`${f.label} · ${f.source || ""}\n${f.tooltip?.interpro_type || f.feature_type}`
                      + `${f.interpro_accession ? ` · ${f.interpro_accession}` : ""}\naa ${f.start}–${f.end}`}</title>
                  </rect>
                ))}
                {!items.length && (
                  <text x={W / 2} y={laneY[k] + LANE_H / 2 + 3} textAnchor="middle" className="da-empty"
                    fill={EMPTY_TEXT.fill} fontSize={EMPTY_TEXT.fontSize} fontStyle={EMPTY_TEXT.fontStyle}>
                    {pending ? "pending post-cluster InterProScan" : "none detected for this protein"}</text>)}
              </g>
            );
          })}

          {/* 7 — pyTMHMM topology (explicit zero / pending messaging, never blank) */}
          {eff.tm && (
            <g>
              {tms.map((t, i) => (
                <rect key={`tm${i}`} x={x(t.start)} y={laneY.tm + 4} width={barW(t.start, t.end)} height={LANE_H - 8}
                  rx="2" fill={TM_PAINT.fill} fillOpacity={TM_PAINT.fillOpacity}
                  stroke={TM_PAINT.stroke} strokeWidth="0.5"
                  style={{ cursor: "pointer" }} onClick={() => openFeatureDetail(t, "tm")}>
                  <title>{`pyTMHMM transmembrane helix\naa ${t.start}–${t.end}`}</title>
                </rect>
              ))}
              {!tms.length && (
                <text x={W / 2} y={laneY.tm + LANE_H / 2 + 3} textAnchor="middle" className="da-empty"
                  fill={EMPTY_TEXT.fill} fontSize={EMPTY_TEXT.fontSize} fontStyle={EMPTY_TEXT.fontStyle}>
                  {tmInfo.pending ? "pyTMHMM pending post-cluster InterProScan"
                    : "No transmembrane region predicted by pyTMHMM"}</text>)}
            </g>
          )}

          {/* 8 — coding exons (projected) */}
          {eff.exons && exons.map((ex) => {
            const sel = ex.id === selExonId;
            const paint = featureProps("coding_exon", { selected: sel });
            const bx = x(ex.start), bw = barW(ex.start, ex.end);
            const tt = ex.tooltip || {};
            return (
              <g key={ex.id}>
                <rect x={bx} y={laneY.exons + 3} width={bw} height={LANE_H - 6} rx="2"
                  className={`em-exon${sel ? " sel" : ""}`}
                  fill={paint.fill} fillOpacity={paint.fillOpacity}
                  stroke={paint.stroke} strokeWidth={paint.strokeWidth} style={{ cursor: "pointer" }}
                  onClick={() => clickExon(ex)}>
                  <title>{`${ex.label} · exon ${tt.exon_number ?? "?"}\nprotein aa ${ex.start}–${ex.end}`
                    + `\ngenomic ${tt.genomic_start ?? "—"}–${tt.genomic_end ?? "—"} · strand ${tt.strand || "?"}`
                    + `\nCDS ${tt.cds_start ?? "—"}–${tt.cds_end ?? "—"} · phase ${tt.phase ?? "n/a"}`}</title>
                </rect>
                {bw >= 18 && <text x={bx + bw / 2} y={laneY.exons + LANE_H / 2 + 3} textAnchor="middle"
                  className="em-exon-num" fill={BLOCK_TEXT.fill} fontSize={BLOCK_TEXT.fontSize}
                  fontWeight={BLOCK_TEXT.fontWeight}
                  style={{ pointerEvents: "none" }}>{ex.label}</text>}
              </g>
            );
          })}

          {/* 9 — exon boundaries */}
          {eff.boundaries && (
            <g>
              {boundaryAAs.map((b, i) => (
                <line key={`b${i}`} x1={x(b)} x2={x(b)} y1={laneY.boundaries + 2} y2={laneY.boundaries + LANE_H - 2}
                  stroke={featureStyle("exon_boundary_tick").stroke} strokeWidth="1" strokeDasharray="2 2" />
              ))}
              {!boundaryAAs.length && (
                <text x={W / 2} y={laneY.boundaries + LANE_H / 2 + 3} textAnchor="middle" className="da-empty"
                  fill={EMPTY_TEXT.fill} fontSize={EMPTY_TEXT.fontSize}
                  fontStyle={EMPTY_TEXT.fontStyle}>no internal coding-exon boundaries</text>)}
            </g>
          )}

          {/* 10 — exploratory candidate clusters, one packed lane per row */}
          {eff.candidates && candLayout.visible.map((c) => {
            const sel = Boolean(selCandidate && c.id === selCandidate.id);
            const paint = featureProps("candidate_region", { selected: sel, faint: !sel });
            const w = barW(c.start, c.end);
            const y = candLaneY(candLayout.laneOf(c));
            return (
              <g key={c.id}>
                <rect x={x(c.start)} y={y} width={w} height={CAND_LANE_H} rx="2"
                  fill={paint.fill} fillOpacity={paint.fillOpacity}
                  stroke={paint.stroke} strokeWidth={paint.strokeWidth}
                  style={{ cursor: "pointer" }} onClick={() => clickCandidate(c)}>
                  <title>{candidateTooltip(c)}</title>
                </rect>
                {candidateLabelFits(w) && (
                  <text x={(x(c.start) + x(c.end)) / 2} y={y + CAND_LANE_H - 3}
                    textAnchor="middle" className="em-exon-num" fill={CAND_TEXT.fill}
                    fontSize={CAND_TEXT.fontSize} fontWeight={CAND_TEXT.fontWeight}
                    style={{ pointerEvents: "none" }}>{c.id}</text>
                )}
              </g>
            );
          })}
          {/* density indication: clusters too narrow to read at this scale are
              counted rather than drawn as an unreadable smear of boxes */}
          {eff.candidates && candLayout.byLane.map((laneItems, i) => {
            const d = laneDensity(laneItems, (c) => barW(c.start, c.end));
            if (d.narrow < 2) return null;
            return (
              <text key={`cd${i}`} x={x(d.end) + 6} y={candLaneY(i) + CAND_LANE_H - 3}
                className="em-axis-lbl" fill={CAND_TEXT.fill} fontSize={CAND_TEXT.fontSize}
                style={{ pointerEvents: "none" }}>{d.narrow} narrow clusters (aa {d.start}–{d.end})</text>
            );
          })}

          <text x={LBLW} y={H - 6} textAnchor="start" className="em-axis-end" fill={AXIS_END_TEXT.fill}
            fontSize={AXIS_END_TEXT.fontSize} fontWeight={AXIS_END_TEXT.fontWeight}>1</text>
          <text x={W - PADR} y={H - 6} textAnchor="end" className="em-axis-end" fill={AXIS_END_TEXT.fill}
            fontSize={AXIS_END_TEXT.fontSize} fontWeight={AXIS_END_TEXT.fontWeight}>{length} aa</text>
        </svg>
      </div>

      <div className="legend res-legend">
        {domains.map((d) => (
          <span key={d.id} className="legend-item"><span className="pa-swatch" style={{ background: domColorMap[d.id] }} />
            {d.label}{d.interpro_accession ? ` (${d.interpro_accession})` : ""}</span>
        ))}
        {eff.families && families.length > 0 && <span className="legend-item"><span className="pa-swatch" style={{ background: FAMILY_FILL }} />family / superfamily</span>}
        {eff.signatures && signatures.length > 0 && <span className="legend-item"><span className="pa-swatch" style={{ background: SIG_FILL }} />member signature</span>}
        {eff.sites && sites.length > 0 && <span className="legend-item"><span className="pa-swatch" style={{ background: SITE_FILL }} />functional site</span>}
        {eff.disorder && disorder.length > 0 && <span className="legend-item"><span className="pa-swatch" style={{ background: DISORDER_FILL }} />disorder / other</span>}
        {eff.tm && tms.length > 0 && <span className="legend-item"><span className="pa-swatch" style={{ background: TM_FILL }} />TM helix (pyTMHMM)</span>}
        <span className="legend-item"><span className="pa-swatch exon" />coding exon</span>
        <span className="legend-item"><span className="pa-swatch cand" />candidate (exploratory)</span>
      </div>

      <DomainDetailDrawer detail={detail} onClose={() => setDetail(null)}
        speciesModel={speciesModel} />
    </div>
  );
}

function DomainDetailDrawer({ detail, onClose, speciesModel }) {
  if (!detail) return null;
  const f = detail.feature;
  const sub = `${speciesModel.scientific_name} · ${speciesModel.protein_id}`;
  if (detail.kind === "domain") {
    const sigByIpr = new Map();
    for (const s of detail.contribSigs) {
      const key = s.interpro_accession || "(unintegrated)";
      if (!sigByIpr.has(key)) sigByIpr.set(key, []);
      sigByIpr.get(key).push(s);
    }
    return (
      <Drawer open title={f.label} subtitle={sub} onClose={onClose}>
        <Field label="Feature type">Representative InterPro domain ({f.tooltip?.interpro_type || "DOMAIN"})</Field>
        <Field label="InterPro accession">{f.interpro_accession
          ? <a href={`https://www.ebi.ac.uk/interpro/entry/InterPro/${f.interpro_accession}/`} target="_blank" rel="noreferrer"><code>{f.interpro_accession}</code></a>
          : "—"}</Field>
        <Field label="Region (aa)">{f.start}–{f.end} · {f.end - f.start + 1} aa</Field>
        <Field label="Source">{f.source || "InterProScan"}</Field>
        <Field label="Contributing member signatures" wide>
          {sigByIpr.size ? [...sigByIpr.entries()].map(([ipr, sigs]) => (
            <div key={ipr} className="sig-group">
              <b>{ipr === "(unintegrated)" ? "unintegrated" : ipr}</b>
              <ul>{sigs.map((s) => (
                <li key={s.id}>{s.source} · <code>{s.signature_accession || "—"}</code>
                  {" "}{s.tooltip?.signature_name || s.label} · aa {s.start}–{s.end}</li>))}</ul>
            </div>
          )) : <span className="muted">none grouped under this entry</span>}
        </Field>
        <Field label="Overlapping coding exons" wide>
          {detail.ovExons.length ? detail.ovExons.map((e) => `${e.label} (aa ${e.start}–${e.end})`).join(", ") : "none"}</Field>
        <Field label="Internal exon boundaries" wide>
          {detail.internalB.length ? detail.internalB.map((b) => `aa ${b}`).join(", ") : "none inside this domain"}</Field>
        <Field label="Nearest candidate region">
          {detail.nearestCand ? `${detail.nearestCand.c.id} (aa ${detail.nearestCand.c.start}–${detail.nearestCand.c.end}`
            + `${detail.nearestCand.dist === 0 ? " · overlaps" : ` · ${detail.nearestCand.dist} aa away`}) — exploratory` : "none"}</Field>
      </Drawer>
    );
  }
  if (detail.kind === "exon") {
    const tt = f.tooltip || {};
    return (
      <Drawer open title={`${f.label} · coding exon`} subtitle={sub} onClose={onClose}>
        <Field label="Projected protein interval">aa {f.start}–{f.end} · {f.end - f.start + 1} aa</Field>
        <Field label="Genomic">{tt.genomic_start ?? "—"}–{tt.genomic_end ?? "—"} · strand {tt.strand || "?"}</Field>
        <Field label="CDS / phase">{tt.cds_start ?? "—"}–{tt.cds_end ?? "—"} · phase {tt.phase ?? "n/a"}</Field>
        <Field label="Overlapping domains / features" wide>
          {detail.ovDomains.length || detail.ovFeatures.length
            ? [...detail.ovDomains, ...detail.ovFeatures].map((d) => `${d.label} (aa ${d.start}–${d.end})`).join(", ")
            : "none"}</Field>
        <Field label="Adjacent exon boundaries">{detail.adjB.length ? detail.adjB.map((b) => `aa ${b}`).join(", ") : "terminal / none"}</Field>
      </Drawer>
    );
  }
  if (detail.kind === "candidate") {
    return (
      <Drawer open title={`${f.label} · exploratory candidate`} subtitle={sub} onClose={onClose}>
        <Field label="Status"><Badge cls="neutral" soft>exploratory candidate — not a confirmed event</Badge></Field>
        <Field label="Region (aa)">{f.start}–{f.end} · {f.end - f.start + 1} aa</Field>
        <Field label="Type / confidence">{f.candidate_type || "candidate"} · {f.confidence || "—"}</Field>
        <Field label="Overlapping domains" wide>
          {detail.ovDomains.length ? detail.ovDomains.map((d) => `${d.label} (aa ${d.start}–${d.end})`).join(", ")
            : detail.nearestDomain ? `none overlapping · nearest: ${detail.nearestDomain.label}` : "no domain context yet (pending cluster)"}</Field>
        <Field label="Affected proteins" wide>{(f.affected_proteins || []).join(", ") || "—"}</Field>
      </Drawer>
    );
  }
  // generic feature (family / signature / site / disorder / tm)
  return (
    <Drawer open title={f.label} subtitle={sub} onClose={onClose}>
      <Field label="Layer">{detail.layer}</Field>
      <Field label="Feature type">{f.tooltip?.interpro_type || f.feature_type || detail.layer}</Field>
      <Field label="Region (aa)">{f.start}–{f.end} · {f.end - f.start + 1} aa</Field>
      <Field label="Source database">{f.source || "—"}</Field>
      <Field label="InterPro accession">{f.interpro_accession
        ? <a href={`https://www.ebi.ac.uk/interpro/entry/InterPro/${f.interpro_accession}/`} target="_blank" rel="noreferrer"><code>{f.interpro_accession}</code></a> : "—"}</Field>
      <Field label="Signature">{f.signature_accession ? <code>{f.signature_accession}</code> : "—"}</Field>
      <Field label="Overlapping coding exons" wide>
        {detail.ovExons?.length ? detail.ovExons.map((e) => `${e.label} (aa ${e.start}–${e.end})`).join(", ") : "none"}</Field>
    </Drawer>
  );
}
