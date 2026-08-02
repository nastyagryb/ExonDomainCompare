import { useMemo, useState } from "react";
import { fileUrl } from "../../api";
import { Badge, Spinner, Empty, Field, Drawer } from "../../ui";
import { architectureStatusLabel } from "../../uiStatus";
import { useIndex, unavailableState } from "./common";
import {
  CHROME, domainInstanceFill, featureProps, featureStyle, textProps,
} from "./semanticStyles";
import {
  FGFR2_BLOCK_LABEL_INK, FGFR2_BLOCK_OUTLINE, FGFR2_CASSETTE_FILL,
  FGFR2_DOMAIN_FALLBACK, FGFR2_DOMAIN_FILL, FGFR2_DOMAIN_LABEL, FGFR2_TM_FILL,
} from "./fgfr2Styles";
import { useScientificSelection } from "../../components/ScientificSelectionContext";

// Post-InterPro / pyTMHMM domain architecture, interactive per species/isoform.
//
// Hard rules honoured in the UI:
//   * IIIb/IIIc labels come from the final truth table (final_isoform_label);
//     InterProScan / pyTMHMM never relabel a cassette.
//   * pyTMHMM is the transmembrane layer (InterProScan annotates no TM helix).
//   * QC warnings are shown but not exaggerated: the audited review cases are
//     surfaced as a minor flag with an explicit coordinate note.

// The FGFR2 architecture, topology and cassette vocabularies are frozen and live in
// fgfr2Styles.js; the generic feature kinds come from semanticStyles.js.
const DOMAIN_FILL = FGFR2_DOMAIN_FILL;
const DOMAIN_LABEL = FGFR2_DOMAIN_LABEL;
const TM_FILL = FGFR2_TM_FILL;
const cassetteFill = (iso) => FGFR2_CASSETTE_FILL[iso] || textProps("muted").fill;

// Text roles and mark paint come from the shared scientific specification and are
// written onto every mark as explicit SVG attributes, so a track stays legible
// without the component stylesheet. The stylesheet keeps layout, cursor and hover.
const AXIS = textProps("axis");
const AXIS_END = textProps("axisEmphasis");
const FEAT = textProps("featureLabel");
const CAND = textProps("candidateLabel");
const EXON = featureProps("coding_exon");

// A dataset is FGFR2-specialised only when the index carries the validated
// cassette panels (panels.{IIIb,IIIc}). Every other gene uses the generic mode
// driven purely by real InterProScan + pyTMHMM annotations.
function isFgfr2CassetteIndex(data) {
  if (!data) return false;
  if (data.mode === "generic") return false;
  return (data.species || []).some((s) => s?.panels && (s.panels.IIIb || s.panels.IIIc));
}

export default function DomainArchitecture(props) {
  const { preloaded } = props;
  const { data, loading } = useIndex((client) => client.domainArchitectureSpecies(), preloaded);
  if (loading) return <Spinner label="Loading domain architecture…" />;
  if (!data?.available) {
    // The run's own verdict, so a run awaiting its cluster round-trip reads as pending work
    // rather than as annotation that could not be found.
    const why = unavailableState(data, "Domain architecture");
    return <Empty title={why.title} hint={why.hint} />;
  }
  if (!isFgfr2CassetteIndex(data)) {
    return <GenericDomainArchitecture data={data} {...props} />;
  }
  return <Fgfr2DomainArchitecture data={data} {...props} />;
}

function Fgfr2DomainArchitecture({ data, species, embedded }) {
  const [panel, setPanel] = useState("IIIb");
  const [sel, setSel] = useState(null);
  const [detail, setDetail] = useState(null);

  const speciesList = useMemo(() => (data?.species || []), [data]);
  const current = useMemo(
    () => speciesList.find((s) => s.species === (sel || species)) || speciesList[0],
    [speciesList, sel, species]
  );

  const available = ["IIIb", "IIIc"].filter((p) => current?.panels?.[p]);
  const activePanel = current?.panels?.[panel] ? panel : (available[0] || "IIIb");
  const pd = current?.panels?.[activePanel];

  return (
    <div className="viewer arch-viewer">
      <div className="viewer-controls">
        <div className="seg">
          {["IIIb", "IIIc"].map((p) => (
            <button key={p} disabled={!current?.panels?.[p]}
              className={`seg-btn iso-tint-${p.toLowerCase()}${activePanel === p ? " on" : ""}`}
              onClick={() => setPanel(p)}>{p}</button>
          ))}
        </div>
        {!embedded && (
          <select value={current?.species || ""} onChange={(e) => setSel(e.target.value)}>
            {speciesList.map((s) => <option key={s.species} value={s.species}>{s.display_species_name}</option>)}
          </select>
        )}
        <span className="spacer" />
        <span className="muted sm">TM layer: pyTMHMM · labels: final truth table</span>
      </div>

      {pd ? (
        <ArchTrack pd={pd} isoform={activePanel} onFeature={setDetail} species={current} />
      ) : (
        <div className="track-card muted">No {activePanel} isoform for {current?.display_species_name}.</div>
      )}

      <ArchLegend />

      <Drawer open={Boolean(detail)} onClose={() => setDetail(null)}
        title={detail?.title || ""} subtitle={detail?.subtitle || ""}>
        {detail?.rows?.map((r) => <Field key={r.k} label={r.k} wide={r.wide}>{r.v}</Field>)}
      </Drawer>
    </div>
  );
}

// --------------------------------------------------------------------------- //
// Generic domain architecture (non-FGFR2): real InterProScan + pyTMHMM only.
// No IIIb/IIIc, no cassette terminology. Species- and protein-selectable.
// --------------------------------------------------------------------------- //
// Colours are assigned to distinct InterPro domain entries in order of appearance —
// no domain-name assumptions. The ramp is the shared instance ramp, so a domain has
// the same colour here and in the exported figure.
const FAMILY_FILL = featureStyle("family_superfamily").fill;  // separate lane
const FEATURE_FILL = featureStyle("functional_site").fill;
const GENERIC_TM_FILL = featureStyle("tm_helix").fill;        // pyTMHMM topology layer

// Build a stable accession -> colour map for one protein's representative domains.
function domainColorMap(domains) {
  const map = {};
  let i = 0;
  for (const d of domains || []) {
    const key = d.interpro_accession || d.domain_id || d.interpro_name;
    if (key && !(key in map)) {
      i += 1;
      map[key] = domainInstanceFill(i);
    }
  }
  return map;
}

function GenericDomainArchitecture({ data, species, embedded }) {
  const selection = useScientificSelection();
  const [detail, setDetail] = useState(null);
  const [showFeatures, setShowFeatures] = useState(false);
  const [showRaw, setShowRaw] = useState(false);
  // User overrides (null = follow the incoming species / default primary). The whole
  // subtree remounts when the selected species changes (provider key), so derivation
  // with a fallback keeps species/protein in sync without setState-in-effect.
  const [selSpecies, setSelSpecies] = useState(null);
  const [selPid, setSelPid] = useState(null);
  const speciesList = useMemo(() => (data?.species || []), [data]);

  const current = useMemo(
    () => speciesList.find((s) => (s.species || s.species_id) === (selSpecies || species))
      || speciesList[0],
    [speciesList, selSpecies, species]
  );
  const proteinList = useMemo(() => (current?.proteins || []), [current]);
  // Only proteins with real InterProScan / pyTMHMM results are selectable — the
  // workflow submits primary proteins only, so unsubmitted isoforms must never
  // appear as though they carry domain annotations.
  const annotatedProteins = useMemo(
    () => proteinList.filter((p) => p.annotated), [proteinList]);
  const primaryPid = useMemo(
    () => (annotatedProteins.find((p) => p.role === "primary")
      || annotatedProteins[0])?.protein_id || null,
    [annotatedProteins]
  );
  const activePid = annotatedProteins.some((p) => p.protein_id === selPid) ? selPid : primaryPid;
  const protein = useMemo(
    () => annotatedProteins.find((p) => p.protein_id === activePid) || annotatedProteins[0],
    [annotatedProteins, activePid]
  );

  if (!current || !protein) {
    return <Empty title="Domain architecture not available"
      hint="No annotated protein was found for the selected species. Only primary proteins are submitted to InterProScan / pyTMHMM." />;
  }

  // Candidate overlay: only on the primary protein of the currently selected
  // species, and only when the linked selection actually points at this species.
  const cand = selection?.selectedCandidate;
  const candRegion = (protein.role === "primary"
    && cand && cand.aa_start != null && cand.aa_end != null
    && (selection?.selectedSpeciesId
        ? selection.selectedSpeciesId === (current.species || current.species_id)
        : true))
    ? { start: cand.aa_start, end: cand.aa_end, label: cand.rank_label || cand.candidate_id }
    : null;

  return (
    <div className="viewer arch-viewer">
      <div className="viewer-controls">
        {!embedded && (
          <select value={current?.species || ""} onChange={(e) => { setSelSpecies(e.target.value); setSelPid(null); }}>
            {speciesList.map((s) => <option key={s.species} value={s.species}>
              {s.display_species_name || s.species}</option>)}
          </select>
        )}
        {/* Show a protein selector only when 2+ annotated proteins exist; otherwise a
            static label (never imply unsubmitted isoforms carry annotations). */}
        {annotatedProteins.length > 1
          ? <select value={protein.protein_id} onChange={(e) => setSelPid(e.target.value)}>
              {annotatedProteins.map((p) => <option key={p.protein_id} value={p.protein_id}>
                {p.protein_id}{p.role === "primary" ? " · primary" : ""}
              </option>)}
            </select>
          : <span className="static-select" title="Only the primary protein was submitted to InterProScan / pyTMHMM">
              <code>{protein.protein_id}</code>{protein.role === "primary" ? " · primary" : ""}
            </span>}
        <span className="spacer" />
        <label className="chk sm"><input type="checkbox" checked={showFeatures}
          onChange={(e) => setShowFeatures(e.target.checked)} /> Sites &amp; disorder</label>
        <label className="chk sm"><input type="checkbox" checked={showRaw}
          onChange={(e) => setShowRaw(e.target.checked)} /> Raw signatures</label>
      </div>
      <p className="muted sm">Representative InterPro domains · families/superfamilies and
        pyTMHMM topology shown separately · real coordinates.</p>

      <GenericArchTrack protein={protein} species={current} candRegion={candRegion}
        showFeatures={showFeatures} onFeature={setDetail} />

      <GenericArchLegend protein={protein} showFeatures={showFeatures} />

      {showRaw && <RawSignatureTable protein={protein} />}

      <Drawer open={Boolean(detail)} onClose={() => setDetail(null)}
        title={detail?.title || ""} subtitle={detail?.subtitle || ""}>
        {detail?.rows?.map((r) => <Field key={r.k} label={r.k} wide={r.wide}>{r.v}</Field>)}
      </Drawer>
    </div>
  );
}

// Raw-signature provenance layer: every member-database hit, for traceability.
function RawSignatureTable({ protein }) {
  const rows = protein.raw_signatures || [];
  if (!rows.length) return null;
  return (
    <details className="raw-sig" open>
      <summary>{rows.length} raw member-database signatures (InterProScan)</summary>
      <div className="table-scroll"><table className="mini-tbl">
        <thead><tr><th>Member DB</th><th>Signature</th><th>Name</th><th>InterPro</th>
          <th>Type</th><th>Region (aa)</th></tr></thead>
        <tbody>{rows.map((r, i) => (
          <tr key={`${r.signature_accession}-${r.start_aa}-${i}`}>
            <td>{r.member_database}</td><td><code>{r.signature_accession}</code></td>
            <td>{r.signature_name || "—"}</td>
            <td>{r.interpro_accession ? <code>{r.interpro_accession}</code> : "—"}</td>
            <td>{r.interpro_type || (r.is_integrated ? "—" : "unintegrated")}</td>
            <td>{r.start_aa}–{r.end_aa}</td>
          </tr>
        ))}</tbody>
      </table></div>
    </details>
  );
}

function GenericArchTrack({ protein, species, candRegion, showFeatures, onFeature }) {
  const speciesName = species.display_species_name || species.species;
  const exons = (protein.exons || []).slice().sort(
    (a, b) => (a.protein_start_aa || 0) - (b.protein_start_aa || 0));
  const domains = protein.domains || [];
  const families = protein.families || [];
  const features = protein.features || [];
  const tms = protein.tm || [];
  const colorMap = domainColorMap(domains);
  const length = protein.length_aa
    || Math.max(1, ...exons.map((e) => e.protein_end_aa || 0),
                ...domains.map((d) => d.end_aa || 0),
                ...families.map((d) => d.end_aa || 0));
  const x = (aa) => PAD + (aa / length) * (W - 2 * PAD);
  const ticks = [];
  const step = niceStep(length);
  for (let t = 0; t <= length; t += step) ticks.push(t);
  if (ticks[ticks.length - 1] !== length) ticks.push(length);
  // internal coding-exon boundaries (exclude the C-terminal end of the last exon)
  const internalBoundaries = exons.slice(0, -1)
    .map((e) => e.protein_end_aa).filter((v) => v != null);

  const openFeature = (title, subtitle, rows) => onFeature({ title, subtitle, rows });
  const domColor = (d) => colorMap[d.interpro_accession || d.domain_id || d.interpro_name] || FGFR2_DOMAIN_FALLBACK;
  const domName = (d) => d.interpro_name || d.domain_name || d.domain_id || "domain";

  return (
    <div className="track-card">
      <div className="track-head">
        <span className="iso">{protein.role === "primary" ? "primary" : "alternative"}</span>
        <span className="track-meta">{length} aa · <code>{protein.protein_id}</code>
          {protein.transcript_id ? <> · <code>{protein.transcript_id}</code></> : null}</span>
        <Badge cls="accepted" soft>{domains.length} domains · {families.length} families · {tms.length} TM</Badge>
      </div>

      <svg className="track-svg" viewBox={`0 0 ${W} 150`} preserveAspectRatio="xMidYMid meet">
        <text x={PAD} y="14" className="axis-label" textAnchor="start"
          fill={AXIS_END.fill} fontSize={AXIS_END.fontSize} fontWeight={AXIS_END.fontWeight}>N</text>
        <text x={W - PAD} y="14" className="axis-label" textAnchor="end"
          fill={AXIS_END.fill} fontSize={AXIS_END.fontSize} fontWeight={AXIS_END.fontWeight}>C</text>
        <line x1={PAD} y1={Y_DOM + H_DOM / 2} x2={W - PAD} y2={Y_DOM + H_DOM / 2} stroke={CHROME.grid} strokeWidth="1" />

        {/* Family / superfamily lane (separate; never a structural domain) */}
        {families.map((d, i) => (
          <rect key={`fam${i}`} x={x(d.start_aa)} y={Y_FAM}
            width={Math.max(2, x(d.end_aa) - x(d.start_aa))} height={H_FAM} rx="2"
            fill={FAMILY_FILL} opacity="0.7"
            onClick={() => openFeature(domName(d), `${speciesName} · ${protein.protein_id}`,
              [{ k: "Layer", v: `Family (${d.interpro_type || "family"})` },
               { k: "Region (aa)", v: `${d.start_aa}–${d.end_aa}` },
               { k: "InterPro", v: d.interpro_accession || "—" },
               { k: "Member DBs", v: d.member_databases || "—" }])}>
            <title>{`${domName(d)} — family/superfamily\nInterPro ${d.interpro_accession}\naa ${d.start_aa}–${d.end_aa}`}</title>
          </rect>
        ))}

        {/* Representative InterPro domains (structural) */}
        {domains.map((d, i) => (
          <g key={`d${i}`}>
            <rect x={x(d.start_aa)} y={Y_DOM} width={Math.max(2, x(d.end_aa) - x(d.start_aa))} height={H_DOM}
              rx="3" fill={domColor(d)} stroke={FGFR2_BLOCK_OUTLINE} strokeWidth="0.4"
              onClick={() => openFeature(domName(d), `${speciesName} · ${protein.protein_id}`,
                [{ k: "Layer", v: `Representative domain (${d.interpro_type || "DOMAIN"})` },
                 { k: "Region (aa)", v: `${d.start_aa}–${d.end_aa}` },
                 { k: "InterPro", v: d.interpro_accession || "—" },
                 { k: "Member DBs", v: d.member_databases || d.domain_source || "—" },
                 { k: "Supporting entries", v: d.supporting_interpro || "—", wide: true }])}>
              <title>{`${domName(d)}\nInterPro ${d.interpro_accession} · ${d.interpro_type || "DOMAIN"}\naa ${d.start_aa}–${d.end_aa}`}</title>
            </rect>
            {x(d.end_aa) - x(d.start_aa) > 30 && (
              <text x={(x(d.start_aa) + x(d.end_aa)) / 2} y={Y_DOM + H_DOM / 2 + 3}
                textAnchor="middle" className="blk-label" fill={FGFR2_BLOCK_LABEL_INK}>{domName(d)}</text>
            )}
          </g>
        ))}

        {/* pyTMHMM transmembrane helices (topology layer) */}
        {tms.map((t, i) => (
          <rect key={`t${i}`} x={x(t.start_aa)} y={Y_DOM - 3}
            width={Math.max(2, x(t.end_aa) - x(t.start_aa))} height={H_DOM + 6} rx="2"
            fill={GENERIC_TM_FILL} stroke={FGFR2_BLOCK_OUTLINE} strokeWidth="0.5"
            onClick={() => openFeature("Transmembrane helix (pyTMHMM)", `${speciesName} · ${protein.protein_id}`,
              [{ k: "Layer", v: "Topology (pyTMHMM)" }, { k: "TM (aa)", v: `${t.start_aa}–${t.end_aa}` }])}>
            <title>{`pyTMHMM transmembrane\naa ${t.start_aa}–${t.end_aa}`}</title>
          </rect>
        ))}

        {/* Sites / disorder feature markers (optional layer) */}
        {showFeatures && features.map((f, i) => (
          <rect key={`f${i}`} x={x(f.start_aa)} y={Y_FEAT}
            width={Math.max(2, x(f.end_aa) - x(f.start_aa))} height={H_FEAT} rx="1"
            fill={FEATURE_FILL} opacity="0.75"
            onClick={() => openFeature(f.interpro_name || "Feature", `${speciesName} · ${protein.protein_id}`,
              [{ k: "Layer", v: `Feature (${f.interpro_type || "site"})` },
               { k: "Region (aa)", v: `${f.start_aa}–${f.end_aa}` },
               { k: "InterPro", v: f.interpro_accession || "—" }])}>
            <title>{`${f.interpro_name} — ${f.interpro_type || "feature"}\naa ${f.start_aa}–${f.end_aa}`}</title>
          </rect>
        ))}

        {/* internal coding-exon boundary ticks */}
        {internalBoundaries.map((b, i) => (
          <line key={`b${i}`} x1={x(b)} y1={Y_EXON - 4} x2={x(b)} y2={Y_EXON + H_EXON + 4}
            stroke={CHROME.rule} strokeWidth="0.6" strokeDasharray="2 2" />
        ))}
        {/* numbered coding exon blocks */}
        {exons.map((e, i) => (
          <g key={`e${i}`}>
            <rect x={x(e.protein_start_aa)} y={Y_EXON}
              width={Math.max(1.5, x(e.protein_end_aa) - x(e.protein_start_aa))} height={H_EXON}
              rx="2" className="exon-block" fill={EXON.fill} fillOpacity={EXON.fillOpacity}
              stroke={EXON.stroke} strokeWidth={EXON.strokeWidth}
              onClick={() => openFeature(`Coding exon ${e.exon_number ?? ""}`, `${speciesName} · ${protein.protein_id}`,
                [{ k: "Exon", v: e.exon_id || String(e.exon_number ?? "") },
                 { k: "Region (aa)", v: `${e.protein_start_aa}–${e.protein_end_aa}` }])}>
              <title>{`exon ${e.exon_number ?? ""}\naa ${e.protein_start_aa}–${e.protein_end_aa}`}</title>
            </rect>
            {x(e.protein_end_aa) - x(e.protein_start_aa) > 12 && e.exon_number != null && (
              <text x={(x(e.protein_start_aa) + x(e.protein_end_aa)) / 2} y={Y_EXON + H_EXON / 2 + 3}
                textAnchor="middle" className="exon-num" fill={FEAT.fill} fontSize={FEAT.fontSize}
                fontWeight={FEAT.fontWeight}>{e.exon_number}</text>
            )}
          </g>
        ))}

        {/* selected exploratory candidate region overlay (primary protein) */}
        {candRegion && (
          <>
            <rect x={x(candRegion.start)} y={Y_DOM - 6}
              width={Math.max(2, x(candRegion.end) - x(candRegion.start))}
              height={(Y_EXON + H_EXON) - (Y_DOM - 6)} rx="3"
              fill={featureStyle("candidate_region").fill} fillOpacity="0.16"
              stroke={CHROME.rule} strokeDasharray="3 2">
              <title>{`Selected candidate ${candRegion.label}\naa ${candRegion.start}–${candRegion.end}`}</title>
            </rect>
            <text x={(x(candRegion.start) + x(candRegion.end)) / 2} y={Y_EXON + H_EXON + 12}
              textAnchor="middle" className="axis-label" fill={CAND.fill} fontSize={CAND.fontSize}
              fontWeight={CAND.fontWeight}>candidate {candRegion.label}</text>
          </>
        )}

        <line x1={PAD} y1={Y_AXIS} x2={W - PAD} y2={Y_AXIS} stroke={CHROME.rule} />
        {ticks.map((t) => (
          <g key={t}>
            <line x1={x(t)} y1={Y_AXIS - 4} x2={x(t)} y2={Y_AXIS + 4} stroke={CHROME.rule} />
            <text x={x(t)} y={Y_AXIS + 16} textAnchor="middle" className="axis-label"
              fill={AXIS.fill} fontSize={AXIS.fontSize}>{t}</text>
          </g>
        ))}
      </svg>
    </div>
  );
}

// Legend generated dynamically from the protein's actual annotations (no
// hard-coded domain names): one swatch per distinct representative domain, plus
// whichever of family / TM / feature layers are present.
function GenericArchLegend({ protein, showFeatures }) {
  const domains = protein.domains || [];
  const colorMap = domainColorMap(domains);
  const seen = new Set();
  const domItems = [];
  for (const d of domains) {
    const key = d.interpro_accession || d.domain_id || d.interpro_name;
    if (key && !seen.has(key)) {
      seen.add(key);
      domItems.push([colorMap[key], `${d.interpro_name || d.domain_name} (${key})`]);
    }
  }
  return (
    <div className="arch-legend">
      {domItems.map(([c, l]) => (
        <span key={l} className="leg-item"><span className="leg-swatch" style={{ background: c }} />{l}</span>
      ))}
      {(protein.families || []).length > 0 && (
        <span className="leg-item"><span className="leg-swatch" style={{ background: FAMILY_FILL }} />family / superfamily</span>)}
      {(protein.tm || []).length > 0 && (
        <span className="leg-item"><span className="leg-swatch" style={{ background: GENERIC_TM_FILL }} />TM helix (pyTMHMM)</span>)}
      {showFeatures && (protein.features || []).length > 0 && (
        <span className="leg-item"><span className="leg-swatch" style={{ background: FEATURE_FILL }} />site / disorder</span>)}
      <span className="leg-item"><span className="leg-swatch exon-block" />coding exon (numbered)</span>
    </div>
  );
}

const W = 980, PAD = 12;
const Y_FAM = 24, H_FAM = 6;      // family / superfamily lane (above domains)
const Y_DOM = 36, H_DOM = 26;     // representative domain + TM band
const Y_FEAT = 66, H_FEAT = 6;    // sites / disorder lane
const Y_EXON = 92, H_EXON = 20;   // numbered exon band
const Y_AXIS = 130;

function ArchTrack({ pd, isoform, onFeature, species }) {
  const length = pd.protein_length || pd.axis?.end || 1;
  const x = (aa) => PAD + (aa / length) * (W - 2 * PAD);
  const qc = pd.qc || {};
  const status = architectureStatusLabel(qc.display_qc_status || qc.final_qc_status);

  const ticks = [];
  const step = niceStep(length);
  for (let t = 0; t <= length; t += step) ticks.push(t);
  if (ticks[ticks.length - 1] !== length) ticks.push(length);

  const openFeature = (title, subtitle, rows) => onFeature({ title, subtitle, rows });

  return (
    <div className={`track-card iso-frame-${isoform.toLowerCase()}`}>
      <div className="track-head">
        <span className={`iso iso-${isoform.toLowerCase()}`}>{pd.final_isoform_label || isoform}</span>
        <span className="track-meta">
          {length} aa · <code>{pd.protein_id || "—"}</code>
          {qc.kinase_found ? " · kinase ✓" : ""}{qc.pytmhmm_tm_found ? " · TM ✓" : ""}
        </span>
        <Badge cls={status.cls} soft title={qc.display_note || undefined}>
          {status.label}
        </Badge>
      </div>

      <svg className="track-svg" viewBox={`0 0 ${W} 132`} preserveAspectRatio="xMidYMid meet">
        {/* N / C terminus labels */}
        <text x={PAD} y="18" className="axis-label" textAnchor="start"
          fill={AXIS_END.fill} fontSize={AXIS_END.fontSize} fontWeight={AXIS_END.fontWeight}>N</text>
        <text x={W - PAD} y="18" className="axis-label" textAnchor="end"
          fill={AXIS_END.fill} fontSize={AXIS_END.fontSize} fontWeight={AXIS_END.fontWeight}>C</text>

        {/* backbone */}
        <line x1={PAD} y1={Y_DOM + H_DOM / 2} x2={W - PAD} y2={Y_DOM + H_DOM / 2} stroke={CHROME.grid} strokeWidth="1" />

        {/* faint full-length family bar (other_domain / FGFR fingerprint) */}
        {(pd.domains || []).filter((d) => d.class === "other_domain").map((d, i) => (
          <rect key={`o${i}`} x={x(d.start)} y={Y_DOM + 4} width={Math.max(2, x(d.end) - x(d.start))} height={H_DOM - 8}
            rx="3" fill={DOMAIN_FILL.other_domain} opacity="0.5"
            onClick={() => openFeature(`${DOMAIN_LABEL.other_domain} · ${d.label}`, `${species.display_species_name} · ${isoform}`,
              [{ k: "Region (AA)", v: `${d.start}–${d.end}` }, { k: "Source", v: d.source }])}>
            <title>{`${DOMAIN_LABEL.other_domain}: ${d.label}\nAA ${d.start}–${d.end}`}</title>
          </rect>
        ))}

        {/* Ig / kinase / signal-peptide domain boxes */}
        {(pd.domains || []).filter((d) => d.class !== "other_domain").map((d, i) => (
          <g key={`d${i}`}>
            <rect x={x(d.start)} y={Y_DOM} width={Math.max(2, x(d.end) - x(d.start))} height={H_DOM}
              rx="3" fill={DOMAIN_FILL[d.class] || FGFR2_DOMAIN_FALLBACK} stroke={FGFR2_BLOCK_OUTLINE} strokeWidth="0.4"
              onClick={() => openFeature(`${DOMAIN_LABEL[d.class] || d.class} · ${d.label}`, `${species.display_species_name} · ${isoform}`,
                [{ k: "Region (AA)", v: `${d.start}–${d.end}` }, { k: "Source", v: d.source }])}>
              <title>{`${DOMAIN_LABEL[d.class] || d.class}: ${d.label}\nAA ${d.start}–${d.end}`}</title>
            </rect>
            {x(d.end) - x(d.start) > 26 && (
              <text x={(x(d.start) + x(d.end)) / 2} y={Y_DOM + H_DOM / 2 + 3} textAnchor="middle"
                className="blk-label" fill={FGFR2_BLOCK_LABEL_INK}>{d.label}</text>
            )}
          </g>
        ))}

        {/* pyTMHMM transmembrane helix + N-terminal anchors */}
        {(pd.tm || []).map((t, i) => (
          <rect key={`t${i}`} x={x(t.start)} y={Y_DOM - 3} width={Math.max(2, x(t.end) - x(t.start))} height={H_DOM + 6}
            rx="2" fill={TM_FILL[t.status] || TM_FILL.receptor_tm} stroke={FGFR2_BLOCK_OUTLINE} strokeWidth="0.5"
            opacity={t.status === "n_terminal_signal_anchor" ? 0.7 : 1}
            onClick={() => openFeature(`Transmembrane (pyTMHMM) · ${t.status}`, `${species.display_species_name} · ${isoform}`,
              [{ k: "TM (AA)", v: `${t.start}–${t.end}` }, { k: "Status", v: t.status }, { k: "Source", v: "pyTMHMM" }])}>
            <title>{`pyTMHMM ${t.status}\nAA ${t.start}–${t.end}`}</title>
          </rect>
        ))}

        {/* numbered coding exon blocks */}
        {(pd.exons || []).filter((e) => !e.is_cassette).map((e, i) => (
          <g key={`e${i}`}>
            <rect x={x(e.start)} y={Y_EXON} width={Math.max(1.5, x(e.end) - x(e.start))} height={H_EXON}
              rx="2" className="exon-block" fill={EXON.fill} fillOpacity={EXON.fillOpacity}
              stroke={EXON.stroke} strokeWidth={EXON.strokeWidth}
              onClick={() => openFeature(`Coding exon ${e.number ?? ""}`, `${species.display_species_name} · ${isoform}`,
                [{ k: "Exon", v: e.label }, { k: "Region (AA)", v: `${e.start}–${e.end}` }, { k: "Source", v: "figure3C exon→protein map" }])}>
              <title>{`${e.label}\nAA ${e.start}–${e.end}`}</title>
            </rect>
            {x(e.end) - x(e.start) > 12 && e.number != null && (
              <text x={(x(e.start) + x(e.end)) / 2} y={Y_EXON + H_EXON / 2 + 3} textAnchor="middle"
                className="exon-num" fill={FEAT.fill} fontSize={FEAT.fontSize}
                fontWeight={FEAT.fontWeight}>{e.number}</text>
            )}
          </g>
        ))}

        {/* IIIb/IIIc cassette slot (highlighted, spanning both bands) */}
        {pd.cassette && pd.cassette.start != null && (
          <>
            <rect x={x(pd.cassette.start)} y={Y_DOM - 6} width={Math.max(2, x(pd.cassette.end) - x(pd.cassette.start))}
              height={(Y_EXON + H_EXON) - (Y_DOM - 6)} rx="3"
              className={`cassette-band band-${isoform.toLowerCase()}`}
              fill={cassetteFill(isoform.toLowerCase())} fillOpacity="0.16"
              stroke={CHROME.rule} strokeDasharray="3 2"
              onClick={() => onFeature({
                title: `${pd.cassette.slot_type} · ${pd.cassette.label}`,
                subtitle: `${species.display_species_name} · ${isoform}`,
                rows: [
                  { k: "Cassette (AA)", v: `${pd.cassette.start}–${pd.cassette.end}` },
                  { k: "Label source", v: "final truth table (InterPro/pyTMHMM never relabel)" },
                  ...(pd.audit ? [{ k: "QC audit", v: pd.audit.final_interpretation, wide: true }] : []),
                ],
              })}>
              <title>{`${pd.cassette.label}\nAA ${pd.cassette.start}–${pd.cassette.end}`}</title>
            </rect>
            <text x={(x(pd.cassette.start) + x(pd.cassette.end)) / 2} y={Y_EXON + H_EXON + 12} textAnchor="middle"
              className={`cass-label iso-${isoform.toLowerCase()}`}
              fill={cassetteFill(isoform.toLowerCase())}>{isoform} cassette{pd.cassette.number ? ` · exon ${pd.cassette.number}` : ""}</text>
          </>
        )}

        {/* AA axis */}
        <line x1={PAD} y1={Y_AXIS} x2={W - PAD} y2={Y_AXIS} stroke={CHROME.rule} />
        {ticks.map((t) => (
          <g key={t}>
            <line x1={x(t)} y1={Y_AXIS - 4} x2={x(t)} y2={Y_AXIS + 4} stroke={CHROME.rule} />
            <text x={x(t)} y={Y_AXIS + 16} textAnchor="middle" className="axis-label"
              fill={AXIS.fill} fontSize={AXIS.fontSize}>{t}</text>
          </g>
        ))}
      </svg>

      {/* Annotation flags stay visible as scientific content; the collapsed
          bookkeeping accordion that used to wrap them is gone (Part 1). */}
      {(qc.display_note || (qc.warnings || []).length > 0) && (
        <div className="arch-qc-note">
          {qc.display_note && <p className="muted sm">{qc.display_note}</p>}
          {(qc.warnings || []).length > 0 && (
            <>
              <span className="field-label">Annotation flags</span>
              <ul className="tp-flags">{qc.warnings.map((w, i) => <li key={i}>{w}</li>)}</ul>
            </>
          )}
        </div>
      )}

      {/* Generated figure variants live in the Figure Gallery only (Part 17). */}
      <div className="arch-links">
        {pd.source_table && <a className="btn ghost sm" href={fileUrl(pd.source_table)}>Feature table (TSV)</a>}
      </div>
    </div>
  );
}

function ArchLegend() {
  const items = [
    [DOMAIN_FILL.ig_like_domain, "Ig-like (InterPro)"],
    [DOMAIN_FILL.kinase_domain, "Kinase (InterPro)"],
    [TM_FILL.receptor_tm, "TM helix (pyTMHMM)"],
    [TM_FILL.n_terminal_signal_anchor, "N-term anchor (pyTMHMM)"],
    [DOMAIN_FILL.other_domain, "FGFR family (InterPro)"],
  ];
  return (
    <div className="arch-legend">
      {items.map(([c, l]) => (
        <span key={l} className="leg-item"><span className="leg-swatch" style={{ background: c }} />{l}</span>
      ))}
      <span className="leg-item"><span className="leg-swatch exon-block" />coding exon (numbered)</span>
      <span className="leg-item"><span className="leg-swatch cass-swatch" />IIIb/IIIc cassette slot</span>
    </div>
  );
}

function niceStep(span) {
  const raw = span / 8;
  const pow = Math.pow(10, Math.floor(Math.log10(raw)));
  const n = raw / pow;
  const m = n >= 5 ? 5 : n >= 2 ? 2 : 1;
  return Math.max(50, m * pow);
}
