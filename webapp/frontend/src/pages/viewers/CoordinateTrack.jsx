import { useEffect, useMemo, useState } from "react";
import { fileUrl } from "../../api";
import { Badge, Spinner, Empty, Field, Drawer, AvailabilityState } from "../../ui";
import { useIndex, unavailableState } from "./common";
import { useScientificSelection } from "../../components/ScientificSelectionContext";
import { domainInstanceFill, featureProps, featureStyle, textProps } from "./semanticStyles";
import { FGFR2_CASSETTE_FILL } from "./fgfr2Styles";

// Paint for the scientific marks of this view, written onto the elements as
// explicit SVG attributes. Component CSS keeps layout, cursor, transition and
// hover; colour and stroke come from the shared specification, so a mark stays
// legible in a standalone SVG, where that stylesheet is absent.
const AXIS_TEXT = textProps("axis");
const AXIS_END_TEXT = textProps("axisEmphasis");
const BLOCK_TEXT = textProps("featureLabel");
const CAND_TEXT = textProps("candidateLabel");
const CAND_FAINT = featureProps("candidate_region", { faint: true });
const CAND_SELECTED = featureProps("candidate_region", { selected: true });
const EXON_PAINT = featureProps("coding_exon");
const SELECTION = featureStyle("selected_feature");

// The IIIb / IIIc cassette belongs to the frozen FGFR2 vocabulary, not to the
// gene-agnostic specification, so its two colours come from that module — which
// mirrors the --iiib / --iiic design tokens the stylesheet uses for the same
// distinction. A panel outside that vocabulary (a generic primary panel) has no
// cassette colour of its own; there the band marks an alternatively spliced
// variable region, which the shared specification does name.
function cassetteColour(panel) {
  return FGFR2_CASSETTE_FILL[String(panel || "").toLowerCase()]
    || featureStyle("variable_region").fill;
}

export default function CoordinateTrack({ preloaded, species, embedded, genericProteins, scientificIndex, compareModels }) {
  const { data, loading } = useIndex((client) => client.coordinates(),
    (genericProteins || scientificIndex) ? {} : preloaded);
  const scientificSelection = useScientificSelection();
  const [hZoom, setHZoom] = useState(1);
  const panelKeys = useMemo(() => {
    const sp = (data?.species || [])[0];
    return Object.keys(sp?.panels || {});
  }, [data]);
  const isFgfr2Iso = panelKeys.includes("IIIb") || panelKeys.includes("IIIc");
  const defaultPanel = isFgfr2Iso ? "IIIb" : (panelKeys[0] || "primary");
  const [panel, setPanel] = useState(defaultPanel);
  const [sel, setSel] = useState(species || null);
  const [compare, setCompare] = useState(false);
  const [zoom, setZoom] = useState(false);
  const [block, setBlock] = useState(null);
  useEffect(() => {
    // eslint-disable-next-line react-hooks/set-state-in-effect
    if (species) setSel(species);
  }, [species]);
  useEffect(() => {
    // eslint-disable-next-line react-hooks/set-state-in-effect
    setPanel(defaultPanel);
  }, [defaultPanel]);

  const speciesList = useMemo(() => (data?.species || []), [data]);
  const current = useMemo(
    () => speciesList.find((s) => s.species === (sel || species)) || speciesList[0],
    [speciesList, sel, species]
  );

  if (scientificIndex?.proteins?.length) {
    const proteins = scientificIndex.proteins;
    const selectedProtein = proteins.find((p) => p.protein_id === scientificSelection?.selectedProteinId)
      || proteins.find((p) => p.protein_id === scientificSelection?.selectedProtein)
      || proteins.find((p) => p.is_primary) || proteins[0];
    const length = selectedProtein.length_aa || 1;
    const x = (pos) => 4 + (Number(pos || 0) / length) * 92;
    const hasDomains = scientificIndex.domain_annotation_status === "available"
      && (selectedProtein.domains?.length || selectedProtein.tm_regions?.length);
    // Candidate overlays are drawn from candidate evidence (reference-protein
    // coordinates), not faked into the index. Only overlay candidates whose
    // reference protein is the one currently on the axis.
    const candidates = (scientificSelection?.rankedCandidates || [])
      .filter((c) => c.reference_protein === selectedProtein.protein_id);
    const selId = scientificSelection?.selectedCandidateId;
    const selExon = scientificSelection?.selectedExon;
    const showExon = selExon && selExon.protein_id === selectedProtein.protein_id
      && selExon.protein_start_aa != null;
    return (
      <div className="viewer coord-viewer">
        <div className="viewer-head">
          <div><b>Exon map</b><p className="muted sm">
            Canonical primary protein backbone with exon-boundary markers and exploratory candidate overlays.
          </p></div>
          <Badge cls={hasDomains ? "accepted" : "neutral"} soft>
            {hasDomains ? "real domain/TM annotation" : "Domain and TM annotation pending cluster"}
          </Badge>
        </div>
        <svg viewBox="0 0 100 26" className="protein-arch-svg" role="img"
          aria-label={`${selectedProtein.protein_id} protein architecture`}>
          <line x1="4" x2="96" y1="12" y2="12" stroke={featureStyle("protein_backbone").fill} strokeWidth="1.3" />
          {/* faint overlays for non-selected candidates */}
          {candidates.filter((c) => c.candidate_id !== selId).map((c) =>
            <rect key={c.candidate_id} x={x(c.aa_start)} y="8"
              width={Math.max(0.7, x(c.aa_end) - x(c.aa_start))} height="8"
              fill={featureStyle("candidate_region").fill} opacity=".16">
              <title>{c.rank_label} · aa {c.aa_start}–{c.aa_end} (exploratory)</title>
            </rect>)}
          {/* strong overlay + label for the selected candidate */}
          {candidates.filter((c) => c.candidate_id === selId).map((c) => (
            <g key={c.candidate_id}>
              <rect x={x(c.aa_start)} y="6" width={Math.max(0.8, x(c.aa_end) - x(c.aa_start))}
                height="12" fill={featureStyle("candidate_region").fill} opacity=".38"
                stroke={featureStyle("candidate_region").stroke} strokeWidth=".4">
                <title>{c.rank_label} · aa {c.aa_start}–{c.aa_end} (selected, exploratory)</title>
              </rect>
              <text x={x(c.aa_start)} y="4.2" fontSize="2.4" fill={textProps("candidateLabel").fill}>
                {c.rank_label} · {c.aa_start}–{c.aa_end}</text>
            </g>
          ))}
          {/* selected exon projection */}
          {showExon &&
            <rect x={x(selExon.protein_start_aa)} y="9.5"
              width={Math.max(0.6, x(selExon.protein_end_aa) - x(selExon.protein_start_aa))}
              height="5" fill={featureStyle("coding_exon").fill} opacity=".55"
              stroke={featureStyle("selected_feature").stroke} strokeWidth=".2">
              <title>Selected exon {selExon.exon_label} · aa {selExon.protein_start_aa}–{selExon.protein_end_aa}</title>
            </rect>}
          {(selectedProtein.domains || []).map((d, i) =>
            <rect key={`${d.signature_accession}:${i}`} x={x(d.start_aa)} y="8"
              width={Math.max(0.8, x(d.end_aa) - x(d.start_aa))} height="8"
              rx="1" fill={domainInstanceFill(i + 1)}>
              <title>{d.source} {d.signature_accession}: {d.description}</title>
            </rect>)}
          {(selectedProtein.tm_regions || []).map((tm, i) =>
            <rect key={`tm:${i}`} x={x(tm.start_aa)} y="6"
              width={Math.max(0.7, x(tm.end_aa) - x(tm.start_aa))} height="12"
              fill={featureStyle("tm_helix").fill}>
              <title>pyTMHMM {tm.start_aa}–{tm.end_aa}</title></rect>)}
          {(selectedProtein.exon_boundaries || []).map((pos) =>
            <line key={pos} x1={x(pos)} x2={x(pos)} y1="5" y2="19"
              stroke={featureStyle("exon_boundary_tick").stroke} strokeWidth=".35" />)}
          <text x="4" y="24" fontSize="2.3" fill={AXIS_TEXT.fill}>1</text>
          <text x="96" y="24" textAnchor="end" fontSize="2.3" fill={AXIS_TEXT.fill}>{length} aa</text>
        </svg>
        <div className="isoform-switch">
          {proteins.map((p) => <button key={p.protein_id}
            className={p.protein_id === selectedProtein.protein_id ? "chip sel" : "chip"}
            onClick={() => scientificSelection?.selectProtein(p)}>
            {p.protein_id}{p.is_primary ? " ★" : ""}
          </button>)}
        </div>
        <div className="legend res-legend">
          <span className="legend-item"><span className="pa-swatch backbone" />protein backbone</span>
          <span className="legend-item"><span className="pa-swatch boundary" />exon boundary</span>
          <span className="legend-item"><span className="pa-swatch cand" />selected candidate</span>
          <span className="legend-item"><span className="pa-swatch candfaint" />other candidates</span>
          <span className="legend-item"><span className="pa-swatch exon" />selected exon</span>
          {hasDomains && <span className="legend-item"><span className="pa-swatch dom" />InterPro domain</span>}
          {hasDomains && <span className="legend-item"><span className="pa-swatch tm" />TM helix</span>}
        </div>
        {!hasDomains && <p className="muted sm">
          No domain or TM feature is displayed before real InterProScan/pyTMHMM parsing. Run the
          cluster round-trip to add domain and TM tracks here.
        </p>}
      </div>
    );
  }
  if (genericProteins) {
    return <Empty title="Coordinate data not available"
      hint="Use coordinate_track_index (FGFR2-compatible) instead of the legacy generic protein track." />;
  }
  if (loading) return <Spinner label="Loading coordinate map…" />;
  if (!data?.available) {
    const why = unavailableState(data, "Coordinate data");
    return <AvailabilityState why={why} />;
  }

  const pd0 = current?.panels?.[panel] || Object.values(current?.panels || {})[0];
  const panels = compare && isFgfr2Iso ? ["IIIb", "IIIc"] : [panel];

  return (
    <div className="viewer coord-viewer">
      <div className="viewer-controls">
        {isFgfr2Iso && !compare && (
          <div className="seg">
            {["IIIb", "IIIc"].map((p) => (
              <button key={p} className={`seg-btn iso-tint-${p.toLowerCase()}${panel === p ? " on" : ""}`} onClick={() => setPanel(p)}>{p}</button>
            ))}
          </div>
        )}
        {!isFgfr2Iso && pd0 && (
          <Badge cls="accepted" soft>{pd0.protein_id || "Primary"}
            {pd0.transcript_id ? ` · ${pd0.transcript_id}` : ""} · {pd0.protein_length || "—"} aa</Badge>
        )}
        {!embedded && (
          <select value={current?.species || ""} onChange={(e) => setSel(e.target.value)}>
            {speciesList.map((s) => <option key={s.species} value={s.species}>{s.display_species_name}</option>)}
          </select>
        )}
        {isFgfr2Iso && <>
          <label className="check inline"><input type="checkbox" checked={compare} onChange={(e) => setCompare(e.target.checked)} /><span>Compare IIIb / IIIc</span></label>
          <label className="check inline"><input type="checkbox" checked={zoom} onChange={(e) => setZoom(e.target.checked)} /><span>Zoom to cassette</span></label>
        </>}
        {!isFgfr2Iso && <>
          <span className="muted small">Expand</span>
          <button className="btn ghost sm" onClick={() => setHZoom((z) => Math.max(1, z - 0.5))} disabled={hZoom <= 1}>−</button>
          <span className="muted small">{hZoom.toFixed(1)}×</span>
          <button className="btn ghost sm" onClick={() => setHZoom((z) => Math.min(5, z + 0.5))} disabled={hZoom >= 5}>+</button>
        </>}
        <span className="spacer" />
        {Object.values(current?.panels || {})[0]?.source_table && (
          <a className="btn ghost sm" href={fileUrl(Object.values(current.panels)[0].source_table)}>Source TSV</a>
        )}
      </div>

      <div className={isFgfr2Iso ? undefined : "coord-scroll"}>
        <div style={isFgfr2Iso ? undefined : { minWidth: `${hZoom * 100}%` }}>
          {panels.map((p) => {
            const pd = current?.panels?.[p];
            return <ProteinTrack key={p} panel={p} pd={pd} zoom={zoom} onBlock={setBlock} species={current}
              selection={scientificSelection} isPrimaryMode={!isFgfr2Iso} />;
          })}
        </div>
      </div>

      {!isFgfr2Iso && (
        <div className="legend res-legend">
          <span className="legend-item"><span className="pa-swatch exon" />coding exon (E1…En)</span>
          <span className="legend-item"><span className="pa-swatch boundary" />exon boundary</span>
          <span className="legend-item"><span className="pa-swatch cand" />selected candidate</span>
          <span className="legend-item"><span className="pa-swatch candfaint" />other candidates</span>
        </div>
      )}

      <div className="domain-placeholder">
        <span className="dp-dot" /> Domain layer pending InterProScan — functional domains are added in a later phase.
      </div>

      {compareModels && !isFgfr2Iso && (current?.models?.length > 1) && (
        <CompareTranscriptModels models={current.models} selection={scientificSelection} />
      )}

      <Drawer open={Boolean(block)} onClose={() => setBlock(null)}
        title={block ? `${block.feature_type} · ${block.label || block.id}` : ""}
        subtitle={block ? `${block.species} · ${block.panel}` : ""}>
        {block && (
          <>
            <div className="drawer-badges">
              {block.in_cassette && <Badge cls="info" soft>in cassette</Badge>}
              {(block.is_iiib_cassette || block.is_iiic_cassette) && <Badge cls="accepted" soft>cassette exon</Badge>}
            </div>
            <Field label="Transcript-relative exon number">{block.transcript_exon_number ?? "—"}</Field>
            <Field label="Exon ID"><code>{block.id || "—"}</code></Field>
            <Field label="Shared exon group ID"><code>{block.shared_exon_group_id || "—"}</code></Field>
            <Field label="Protein start / end (AA)">{block.start ?? "—"}–{block.end ?? "—"}</Field>
            <Field label="Genomic start / end">{block.genomic_start ?? "—"}–{block.genomic_end ?? "—"}</Field>
            <Field label="CDS start / end">{block.cds_start ?? "—"}–{block.cds_end ?? "—"}</Field>
            <Field label="Phase">{block.phase || "—"}</Field>
            <Field label="Strand">{block.strand || "—"}</Field>
            {block.cassette_available !== false && (block.cassette_start_aa != null) && (
              <Field label="Cassette span (AA)">{block.cassette_start_aa}–{block.cassette_end_aa}</Field>
            )}
            <Field label="Transcript"><code>{block.transcript || "—"}</code></Field>
            <Field label="Source">{block.source || block.plot_status || "—"}</Field>
            <Field label="Source table" wide><code>{block.source_table}</code></Field>
          </>
        )}
      </Drawer>
    </div>
  );
}

const W = 960, PAD = 8;

function ProteinTrack({ panel, pd, zoom, onBlock, species, selection, isPrimaryMode }) {
  if (!pd) return <div className="track-card muted">No {panel} isoform for {species?.display_species_name}.</div>;
  const length = pd.protein_length || (pd.blocks.reduce((m, b) => Math.max(m, b.end || 0), 0)) || 1;
  const cassOk = pd.cassette_available;
  const lo = zoom && cassOk ? Math.max(0, pd.cassette_start_aa - 12) : 0;
  const hi = zoom && cassOk ? Math.min(length, pd.cassette_end_aa + 12) : length;
  const span = Math.max(1, hi - lo);
  const x = (aa) => PAD + ((aa - lo) / span) * (W - 2 * PAD);

  const indexCandidates = pd.candidate_regions || [];
  const ctxCandidates = (selection?.rankedCandidates || []).filter(
    (c) => !pd.protein_id || c.reference_protein === pd.protein_id);
  const candidates = ctxCandidates.length ? ctxCandidates : indexCandidates;
  const selId = selection?.selectedCandidateId
    || (selection?.selectedCandidate?.candidate_id);
  const panelLabel = isPrimaryMode ? (pd.final_isoform_label || "Primary") : panel;

  const ticks = [];
  const step = niceStep(span);
  for (let t = Math.ceil(lo / step) * step; t <= hi; t += step) ticks.push(t);

  return (
    <div className={`track-card iso-frame-${isPrimaryMode ? "primary" : panel.toLowerCase()}`}>
      <div className="track-head">
        <span className={`iso iso-${isPrimaryMode ? "neutral" : panel.toLowerCase()}`}>{panelLabel}</span>
        <span className="track-meta">{length} aa
          {pd.transcript_id ? ` · ${pd.transcript_id}` : ""}
          {cassOk ? ` · cassette ${pd.cassette_start_aa}–${pd.cassette_end_aa}` : ""}</span>
        {pd.is_review ? <Badge cls="review" soft>review</Badge>
          : <Badge cls={pd.claim_class || "accepted"} soft>{pd.final_plot_status || "primary"}</Badge>}
      </div>
      <svg className="track-svg" viewBox={`0 0 ${W} 92`} preserveAspectRatio="xMidYMid meet">
        {/* axis */}
        <line x1={PAD} y1="70" x2={W - PAD} y2="70" stroke="#d2d9e6" />
        {ticks.map((t) => (
          <g key={t}>
            <line x1={x(t)} y1="66" x2={x(t)} y2="74" stroke="#d2d9e6" />
            <text x={x(t)} y="86" textAnchor="middle" className="axis-label"
              fill={AXIS_TEXT.fill} fontSize={AXIS_TEXT.fontSize}>{t}</text>
          </g>
        ))}
        {/* exploratory candidate overlays (faint + selected) */}
        {candidates.filter((c) => (c.candidate_id || c.rank_label) !== selId
          && c.candidate_id !== selId).map((c) => {
          const start = c.aa_start, end = c.aa_end;
          if (start == null || end == null) return null;
          return (
            <rect key={c.candidate_id || c.rank_label} x={x(start)} y="26"
              width={Math.max(2, x(end) - x(start))} height="30" rx="2"
              className="cand-overlay faint" fill={CAND_FAINT.fill}
              fillOpacity={CAND_FAINT.fillOpacity} stroke="none">
              <title>{`${c.rank_label || c.candidate_id} · aa ${start}–${end} (exploratory)`}</title>
            </rect>
          );
        })}
        {candidates.filter((c) => c.candidate_id === selId
          || (c.rank_label && selection?.selectedCandidate?.rank_label === c.rank_label)).map((c) => {
          const start = c.aa_start, end = c.aa_end;
          if (start == null || end == null) return null;
          return (
            <g key={`sel-${c.candidate_id || c.rank_label}`}>
              <rect x={x(start)} y="22" width={Math.max(2, x(end) - x(start))} height="38" rx="3"
                className="cand-overlay selected" fill={CAND_SELECTED.fill}
                fillOpacity={CAND_SELECTED.fillOpacity} stroke={CAND_SELECTED.stroke}
                strokeWidth={CAND_SELECTED.strokeWidth} />
              <text x={x(start)} y="18" className="cass-label" fill={CAND_TEXT.fill}
                fontSize={CAND_TEXT.fontSize} fontWeight={CAND_TEXT.fontWeight}>
                {c.rank_label || "C1"} · {start}–{end}</text>
            </g>
          );
        })}
        {/* exon boundary markers (thin ticks at every projected boundary) */}
        {isPrimaryMode && pd.blocks.map((b, idx) => (
          b.start != null ? <line key={`bl${idx}`} x1={x(b.start)} y1="28" x2={x(b.start)} y2="54"
            stroke={featureStyle("exon_boundary_tick").stroke} strokeWidth="0.4" /> : null
        ))}
        {/* exon/CDS blocks (with E1..En labels when width allows) */}
        {pd.blocks.map((b, idx) => {
          if (b.start == null || b.end == null) return null;
          const bx = x(b.start), bw = Math.max(1.5, x(b.end) - x(b.start));
          const isSelExon = isPrimaryMode && selection?.selectedExonId
            && (b.id === selection.selectedExonId);
          const showLabel = isPrimaryMode && bw >= 22 && b.label;
          return (
            <g key={idx}>
              <rect x={bx} y="30" width={bw} height="22" rx="2"
                className={`exon-block${b.in_cassette ? " in-cass" : ""}${isSelExon ? " sel-exon" : ""}`}
                fill={EXON_PAINT.fill} fillOpacity={EXON_PAINT.fillOpacity}
                stroke={isSelExon ? SELECTION.stroke : EXON_PAINT.stroke}
                strokeWidth={isSelExon ? SELECTION.strokeWidth : EXON_PAINT.strokeWidth}
                // .exon-block carries a stroke of its own, and a class rule outranks a
                // presentation attribute, so the selection outline has to be inline to
                // remain visible on screen.
                style={isSelExon ? { stroke: SELECTION.stroke,
                  strokeWidth: SELECTION.strokeWidth } : undefined}
                onClick={() => {
                  onBlock({ ...b, panel, species: species.display_species_name,
                    cassette_start_aa: pd.cassette_start_aa, cassette_end_aa: pd.cassette_end_aa,
                    plot_status: pd.final_plot_status, transcript: b.transcript_id || pd.transcript_id || "",
                    source_table: pd.source_table });
                  if (isPrimaryMode && selection?.selectExon) {
                    selection.selectExon({ exon_id: b.id, transcript_id: b.transcript_id || pd.transcript_id,
                      protein_id: pd.protein_id, protein_start_aa: b.start, protein_end_aa: b.end });
                  }
                }}>
                <title>{`${b.label || b.id} · exon ${b.transcript_exon_number ?? "?"}\n`
                  + `AA ${b.start}–${b.end} · phase ${b.phase || "n/a"} · ${b.strand || ""}\n`
                  + `genomic ${b.genomic_start ?? "—"}–${b.genomic_end ?? "—"}\n`
                  + `CDS ${b.cds_start ?? "—"}–${b.cds_end ?? "—"}\n`
                  + `exon id ${b.id || "—"} · shared group ${b.shared_exon_group_id || "—"}`}</title>
              </rect>
              {showLabel && <text x={bx + bw / 2} y="45" textAnchor="middle" className="exon-lbl"
                fontSize="9" fill={BLOCK_TEXT.fill} fontWeight={BLOCK_TEXT.fontWeight}
                style={{ pointerEvents: "none" }}>{b.label}</text>}
            </g>
          );
        })}
        {/* cassette overlay + boundary markers */}
        {cassOk && (
          <>
            <rect x={x(pd.cassette_start_aa)} y="24" width={Math.max(2, x(pd.cassette_end_aa) - x(pd.cassette_start_aa))} height="34" rx="3"
              className={`cassette-band band-${panel.toLowerCase()}`}
              fill={cassetteColour(panel)} fillOpacity="0.28"
              stroke={cassetteColour(panel)} strokeWidth="1">
              <title>{`${panel} cassette · AA ${pd.cassette_start_aa}–${pd.cassette_end_aa}`}</title>
            </rect>
            <line x1={x(pd.cassette_start_aa)} y1="18" x2={x(pd.cassette_start_aa)} y2="62" className={`boundary b-${panel.toLowerCase()}`}
              stroke={cassetteColour(panel)} strokeWidth="1.5" />
            <line x1={x(pd.cassette_end_aa)} y1="18" x2={x(pd.cassette_end_aa)} y2="62" className={`boundary b-${panel.toLowerCase()}`}
              stroke={cassetteColour(panel)} strokeWidth="1.5" />
            <text x={x(pd.cassette_start_aa)} y="14" textAnchor="middle" className="cass-label"
              fill={AXIS_END_TEXT.fill} fontSize={AXIS_END_TEXT.fontSize}
              fontWeight={AXIS_END_TEXT.fontWeight}>{pd.cassette_start_aa}</text>
            <text x={x(pd.cassette_end_aa)} y="14" textAnchor="middle" className="cass-label"
              fill={AXIS_END_TEXT.fill} fontSize={AXIS_END_TEXT.fontSize}
              fontWeight={AXIS_END_TEXT.fontWeight}>{pd.cassette_end_aa}</text>
          </>
        )}
      </svg>
      {isPrimaryMode && candidates.length > 0 && (
        <div className="legend res-legend">
          <span className="legend-item"><span className="pa-swatch cand" />selected candidate (exploratory)</span>
          <span className="legend-item"><span className="pa-swatch candfaint" />other candidates</span>
        </div>
      )}
    </div>
  );
}

function niceStep(span) {
  const raw = span / 8;
  const pow = Math.pow(10, Math.floor(Math.log10(raw)));
  const n = raw / pow;
  const m = n >= 5 ? 5 : n >= 2 ? 2 : 1;
  return Math.max(10, m * pow);
}

// Inline "Compare transcript models" panel. Projects each protein model's coding
// exons onto a COMMON absolute amino-acid scale so differing model ends are
// visible, with a synchronized aa ruler, per-row candidate band, data-driven
// filters (differences-only / affected-by-candidate), cross-row hover linking and
// a rich exon tooltip. All logic is generic (no gene- or transcript-specific IDs).
const CMP_W = 960, CMP_PAD = 8;

function exonGroupKey(b) {
  return b.shared_exon_group_id || b.id || `${b.start}-${b.end}`;
}

function overlapsCandidate(block, cand) {
  if (!cand || block.start == null || block.end == null) return false;
  return block.start <= cand.aa_end && block.end >= cand.aa_start;
}

function CompareTranscriptModels({ models, selection }) {
  const [open, setOpen] = useState(true);
  const [mode, setMode] = useState("all"); // all | diff | affected
  const [hoverGroup, setHoverGroup] = useState(null);
  const [hoverModel, setHoverModel] = useState(null);
  const cand = selection?.selectedCandidate;
  const affected = cand ? selection.affectedProteinsFor(cand.candidate_id) : null;

  const ordered = useMemo(() => [...models].sort(
    (a, b) => (b.is_primary ? 1 : 0) - (a.is_primary ? 1 : 0)), [models]);
  const primary = ordered.find((m) => m.is_primary) || ordered[0];
  const maxLen = Math.max(...ordered.map((m) => m.protein_length || 0), 1);
  const totalModels = ordered.length;

  // How many models share each exon group — drives the "shared by x/y" tooltip
  // and the data-driven "differences only" filter. Fully generic.
  const groupCounts = useMemo(() => {
    const counts = new Map();
    for (const m of ordered) {
      const seen = new Set();
      for (const b of m.blocks || []) {
        const k = exonGroupKey(b);
        if (seen.has(k)) continue;
        seen.add(k);
        counts.set(k, (counts.get(k) || 0) + 1);
      }
    }
    return counts;
  }, [ordered]);

  const primaryGroups = useMemo(
    () => new Set((primary?.blocks || []).map(exonGroupKey)), [primary]);
  const differsFromPrimary = (m) => {
    if (m.is_primary) return false;
    if ((m.protein_length ?? null) !== (primary?.protein_length ?? null)) return true;
    const g = (m.blocks || []).map(exonGroupKey);
    if (g.length !== primaryGroups.size) return true;
    return g.some((k) => !primaryGroups.has(k));
  };
  const modelAffected = (m) => {
    if (m.is_primary) return true;
    if (affected?.has(m.protein_id)) return true;
    return (m.blocks || []).some((b) => overlapsCandidate(b, cand));
  };

  const rows = ordered.filter((m) => {
    if (mode === "affected") return modelAffected(m);
    if (mode === "diff") return m.is_primary || differsFromPrimary(m);
    return true;
  });

  const x = (aa) => CMP_PAD + (Math.max(0, aa) / maxLen) * (CMP_W - 2 * CMP_PAD);
  const ticks = [];
  const step = niceStep(maxLen);
  for (let t = 0; t <= maxLen; t += step) ticks.push(t);

  return (
    <details className="tech-prov compare-models" open={open}
      onToggle={(e) => setOpen(e.currentTarget.open)}>
      <summary>Compare transcript models ({models.length})</summary>
      <div className="cmp-controls">
        <div className="seg">
          {[["all", "All models"], ["diff", "Differences only"], ["affected", "Affected by candidate"]].map(([id, l]) => (
            <button key={id} className={`seg-btn${mode === id ? " on" : ""}`}
              onClick={() => setMode(id)} disabled={id === "affected" && !cand}>{l}</button>
          ))}
        </div>
        <span className="muted small">Common absolute aa scale · 1–{maxLen} aa · primary first
          {mode !== "all" ? ` · ${rows.length}/${totalModels} models` : ""}</span>
      </div>

      {/* Synchronized amino-acid ruler (shares the x scale with every model row). */}
      <div className="cmp-row cmp-ruler-row">
        <div className="cmp-label cmp-ruler-label"><span className="muted small">aa</span></div>
        <svg className="cmp-ruler" viewBox={`0 0 ${CMP_W} 18`} preserveAspectRatio="none">
          <line x1={x(0)} y1="4" x2={x(maxLen)} y2="4" stroke="#d2d9e6" strokeWidth="0.6" />
          {ticks.map((t) => (
            <g key={t}>
              <line x1={x(t)} y1="2" x2={x(t)} y2="7" stroke="#d2d9e6" strokeWidth="0.6" />
              <text x={x(t)} y="15" textAnchor={t === 0 ? "start" : "middle"} className="cmp-tick"
                fill={AXIS_TEXT.fill} fontSize={AXIS_TEXT.fontSize}>{t}</text>
            </g>
          ))}
          <text x={x(maxLen)} y="15" textAnchor="end" className="cmp-tick"
            fill={AXIS_TEXT.fill} fontSize={AXIS_TEXT.fontSize}>{maxLen} aa</text>
        </svg>
      </div>

      <div className="cmp-rows">
        {rows.map((m) => {
          const exonCount = (m.blocks || []).length;
          const lenDiff = primary && m.protein_length != null && primary.protein_length != null
            && m.protein_length !== primary.protein_length
            ? `${m.protein_length > primary.protein_length ? "+" : ""}${m.protein_length - primary.protein_length} aa`
            : null;
          const isHoverModel = hoverModel === m.protein_id;
          return (
            <div key={m.protein_id}
              className={`cmp-row${m.is_primary ? " is-primary" : ""}${isHoverModel ? " row-hover" : ""}`}
              onMouseEnter={() => setHoverModel(m.protein_id)}
              onMouseLeave={() => setHoverModel((cur) => (cur === m.protein_id ? null : cur))}>
              <div className="cmp-label">
                <b className="cmp-pid">{m.protein_id}{m.is_primary ? " ★" : ""}</b>
                <span className="cmp-tx">{m.transcript_id || "—"}</span>
                <span className="cmp-dim">{exonCount} exons • {m.protein_length ?? "—"} aa
                  {lenDiff ? ` (Δ ${lenDiff})` : ""}</span>
                <span className="cmp-tags">
                  <Badge cls={m.curation_status === "curated" ? "accepted" : "neutral"} soft>
                    {m.curation_status || "predicted"}</Badge>
                  <Badge cls={m.is_primary ? "accepted" : "neutral"} soft>
                    {m.is_primary ? "primary" : "alternative"}</Badge>
                </span>
              </div>
              <svg className="cmp-track" viewBox={`0 0 ${CMP_W} 26`} preserveAspectRatio="none">
                <line x1={x(0)} y1="12" x2={x(m.protein_length || maxLen)} y2="12"
                  stroke="#d2d9e6" strokeWidth="0.6" />
                {/* per-row candidate band (identical region on every model row) */}
                {cand && (() => {
                  const rowOverlaps = (m.blocks || []).some((b) => overlapsCandidate(b, cand));
                  // The band marks the same candidate on every row; rows it does not
                  // touch keep it faint instead of dropping the context.
                  const band = featureProps("candidate_region", { faint: !rowOverlaps });
                  return (
                    <rect x={x(cand.aa_start)} y="2"
                      width={Math.max(2, x(cand.aa_end) - x(cand.aa_start))} height="22" rx="2"
                      className={`cmp-cand${rowOverlaps ? " overlaps" : ""}`}
                      fill={band.fill} fillOpacity={band.fillOpacity}
                      stroke={band.stroke} strokeWidth={band.strokeWidth}>
                      <title>{`${cand.rank_label} · aa ${cand.aa_start}–${cand.aa_end}`
                        + `${rowOverlaps ? " · overlaps this model" : ""}`}</title>
                    </rect>
                  );
                })()}
                {(m.blocks || []).map((b, i) => {
                  if (b.start == null || b.end == null) return null;
                  const bx = x(b.start), bw = Math.max(1.2, x(b.end) - x(b.start));
                  const gk = exonGroupKey(b);
                  const shared = groupCounts.get(gk) || 1;
                  const isHl = hoverGroup && gk === hoverGroup;
                  const len = (b.end - b.start) + 1;
                  return (
                    <rect key={i} x={bx} y="4" width={bw} height="16" rx="1.5"
                      className={`exon-block${isHl ? " hl" : ""}`}
                      fill={EXON_PAINT.fill} fillOpacity={EXON_PAINT.fillOpacity}
                      stroke={EXON_PAINT.stroke} strokeWidth={EXON_PAINT.strokeWidth}
                      onMouseEnter={() => setHoverGroup(gk)}
                      onMouseLeave={() => setHoverGroup((cur) => (cur === gk ? null : cur))}
                      onClick={() => selection?.selectExon({ exon_id: b.id, transcript_id: b.transcript_id,
                        protein_id: m.protein_id, protein_start_aa: b.start, protein_end_aa: b.end })}>
                      <title>{`${b.label || b.id} · exon ${b.transcript_exon_number ?? "?"}\n`
                        + `protein aa ${b.start}–${b.end} · length ${len} aa\n`
                        + `genomic ${b.genomic_start ?? "—"}–${b.genomic_end ?? "—"}\n`
                        + `CDS ${b.cds_start ?? "—"}–${b.cds_end ?? "—"}\n`
                        + `phase ${b.phase ?? "n/a"} · strand ${b.strand || "?"}\n`
                        + `shared by ${shared}/${totalModels} transcript models`}</title>
                    </rect>
                  );
                })}
              </svg>
            </div>
          );
        })}
      </div>
      <p className="muted sm">Each row is one protein model on the same absolute amino-acid axis (top
        ruler). The orange band marks the selected exploratory candidate on every model; a solid band
        means the region overlaps that model. Hover an exon to highlight the same exon group across all
        models; the tooltip reports coordinates, phase, strand and how many models share the exon.</p>
    </details>
  );
}
