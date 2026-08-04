import { Badge } from "../../ui";
import { architectureStatusLabel } from "../../uiStatus";
import { classInfo, boundaryTypeLabel } from "./boundary";
import { CHROME, featureProps, textProps } from "./semanticStyles";
import {
  FGFR2_BLOCK_LABEL_INK, FGFR2_BLOCK_OUTLINE, FGFR2_CASSETTE_FILL,
  FGFR2_DOMAIN_FALLBACK, FGFR2_DOMAIN_FILL, FGFR2_DOMAIN_LABEL, FGFR2_TM_FILL,
} from "./fgfr2Styles";

// Protein-level detail track for one (species, isoform). Reuses the Domain
// Architecture panel data (domains / TM / kinase / numbered coding exons /
// cassette slot) and overlays the cassette boundary markers with their distance
// to the nearest protein-domain boundary.
const DOMAIN_FILL = FGFR2_DOMAIN_FILL;
const DOMAIN_LABEL = FGFR2_DOMAIN_LABEL;
const TM_FILL = FGFR2_TM_FILL;
const cassetteFill = (iso) => FGFR2_CASSETTE_FILL[iso] || textProps("muted").fill;

// Text roles and mark paint come from the shared scientific specification and are
// written onto every mark as explicit SVG attributes, so the track stays legible
// without the component stylesheet.
const AXIS = textProps("axis");
const AXIS_END = textProps("axisEmphasis");
const FEAT = textProps("featureLabel");
const EXON = featureProps("coding_exon");

const W = 980, PAD = 14;
const Y_DOM = 42, H_DOM = 26;
const Y_EXON = 86, H_EXON = 18;
const Y_AXIS = 124;

function niceStep(span) {
  const raw = span / 8;
  const pow = Math.pow(10, Math.floor(Math.log10(raw)));
  const n = raw / pow;
  const m = n >= 5 ? 5 : n >= 2 ? 2 : 1;
  return Math.max(50, m * pow);
}

export default function BoundaryDetailTrack({ pd, row, isoform, embedded }) {
  const cells = row?.cells || [];
  const length = pd?.protein_length || row?.protein_length || pd?.axis?.end || 1;
  const x = (aa) => PAD + (aa / length) * (W - 2 * PAD);
  const isoLc = String(isoform || row?.isoform || "").toLowerCase();

  const ticks = [];
  const step = niceStep(length);
  for (let t = 0; t <= length; t += step) ticks.push(t);
  if (ticks[ticks.length - 1] !== length) ticks.push(length);

  const domains = pd?.domains || [];
  const tms = pd?.tm || [];
  const exons = (pd?.exons || []).filter((e) => !e.is_cassette);
  const cassette = pd?.cassette;
  const qc = pd?.qc || {};

  return (
    <div className={`viewer bc-detail${embedded ? " embedded" : ""}`}>
      <div className="track-card">
        <div className="track-head">
          <span className={`iso iso-${isoLc}`}>{pd?.final_isoform_label || row?.display_label || isoform}</span>
          <span className="track-meta">{length} aa · <code>{pd?.protein_id || row?.protein_id || "—"}</code></span>
          {(row?.architecture_qc_status || qc.display_qc_status) && (() => {
            const st = architectureStatusLabel(qc.display_qc_status || row.architecture_qc_status);
            return <Badge cls={st.cls} soft>{st.label}</Badge>;
          })()}
        </div>

        <svg className="track-svg" viewBox={`0 0 ${W} 150`} preserveAspectRatio="xMidYMid meet">
          <text x={PAD} y="18" className="axis-label" textAnchor="start"
            fill={AXIS_END.fill} fontSize={AXIS_END.fontSize} fontWeight={AXIS_END.fontWeight}>N</text>
          <text x={W - PAD} y="18" className="axis-label" textAnchor="end"
            fill={AXIS_END.fill} fontSize={AXIS_END.fontSize} fontWeight={AXIS_END.fontWeight}>C</text>
          <line x1={PAD} y1={Y_DOM + H_DOM / 2} x2={W - PAD} y2={Y_DOM + H_DOM / 2} stroke={CHROME.grid} strokeWidth="1" />

          {/* family bar */}
          {domains.filter((d) => d.class === "other_domain").map((d, i) => (
            <rect key={`o${i}`} x={x(d.start)} y={Y_DOM + 4} width={Math.max(2, x(d.end) - x(d.start))}
              height={H_DOM - 8} rx="3" fill={DOMAIN_FILL.other_domain} opacity="0.5" />
          ))}

          {/* Ig / kinase / signal domains */}
          {domains.filter((d) => d.class !== "other_domain").map((d, i) => (
            <g key={`d${i}`}>
              <rect x={x(d.start)} y={Y_DOM} width={Math.max(2, x(d.end) - x(d.start))} height={H_DOM}
                rx="3" fill={DOMAIN_FILL[d.class] || FGFR2_DOMAIN_FALLBACK} stroke={FGFR2_BLOCK_OUTLINE} strokeWidth="0.4">
                <title>{`${DOMAIN_LABEL[d.class] || d.class}: ${d.label}\nAA ${d.start}–${d.end}`}</title>
              </rect>
              {x(d.end) - x(d.start) > 26 && (
                <text x={(x(d.start) + x(d.end)) / 2} y={Y_DOM + H_DOM / 2 + 3} textAnchor="middle"
                  className="blk-label" fill={FGFR2_BLOCK_LABEL_INK}>{d.label}</text>
              )}
            </g>
          ))}

          {/* pyTMHMM TM */}
          {tms.map((t, i) => (
            <rect key={`t${i}`} x={x(t.start)} y={Y_DOM - 3} width={Math.max(2, x(t.end) - x(t.start))}
              height={H_DOM + 6} rx="2" fill={TM_FILL[t.status] || TM_FILL.receptor_tm}
              stroke={FGFR2_BLOCK_OUTLINE} strokeWidth="0.5" opacity={t.status === "n_terminal_signal_anchor" ? 0.7 : 1}>
              <title>{`pyTMHMM ${t.status}\nAA ${t.start}–${t.end}`}</title>
            </rect>
          ))}

          {/* numbered coding exons */}
          {exons.map((e, i) => (
            <g key={`e${i}`}>
              <rect x={x(e.start)} y={Y_EXON} width={Math.max(1.5, x(e.end) - x(e.start))} height={H_EXON}
                rx="2" className="exon-block" fill={EXON.fill} fillOpacity={EXON.fillOpacity}
                stroke={EXON.stroke} strokeWidth={EXON.strokeWidth}>
                <title>{`${e.label}\nAA ${e.start}–${e.end}`}</title>
              </rect>
              {x(e.end) - x(e.start) > 12 && e.number != null && (
                <text x={(x(e.start) + x(e.end)) / 2} y={Y_EXON + H_EXON / 2 + 3} textAnchor="middle"
                  className="exon-num" fill={FEAT.fill} fontSize={FEAT.fontSize}
                  fontWeight={FEAT.fontWeight}>{e.number}</text>
              )}
            </g>
          ))}

          {/* cassette slot band */}
          {cassette && cassette.start != null && (
            <rect x={x(cassette.start)} y={Y_DOM - 6} width={Math.max(2, x(cassette.end) - x(cassette.start))}
              height={(Y_EXON + H_EXON) - (Y_DOM - 6)} rx="3"
              className={`cassette-band band-${isoLc}`}
              fill={cassetteFill(isoLc)} fillOpacity="0.16"
              stroke={CHROME.rule} strokeDasharray="3 2">
              <title>{`${cassette.label}\nAA ${cassette.start}–${cassette.end}`}</title>
            </rect>
          )}

          {/* cassette boundary markers + nearest-domain-boundary connectors */}
          {cells.map((c, i) => {
            const bx = x(c.boundary_aa);
            const info = classInfo(c.boundary_class);
            // nearest domain edge closest to this boundary
            const edges = [c.nearest_domain_start_aa, c.nearest_domain_end_aa].filter((v) => v != null);
            const nearEdge = edges.length
              ? edges.reduce((a, b) => (Math.abs(b - c.boundary_aa) < Math.abs(a - c.boundary_aa) ? b : a))
              : null;
            return (
              <g key={`m${i}`}>
                {nearEdge != null && (
                  <line x1={bx} y1={Y_DOM - 12} x2={x(nearEdge)} y2={Y_DOM - 12}
                    stroke={info.color} strokeWidth="2" opacity="0.8" />
                )}
                <line x1={bx} y1={Y_DOM - 14} x2={bx} y2={Y_EXON + H_EXON + 4}
                  stroke={info.color} strokeWidth="1.6" strokeDasharray="2 2" />
                <circle cx={bx} cy={Y_DOM - 14} r="3.4" fill={info.color} stroke="#fff" strokeWidth="0.8">
                  <title>{`${boundaryTypeLabel(c.boundary_type)} · AA ${c.boundary_aa}\n${info.label} · ${c.distance_to_nearest_domain_boundary} aa from ${c.nearest_domain_label}`}</title>
                </circle>
                <text x={bx} y={Y_EXON + H_EXON + 16} textAnchor="middle" className="bc-marker-label" fill={info.color}>
                  {c.boundary_type === "cassette_start" ? "start" : "end"} · {c.distance_to_nearest_domain_boundary}aa
                </text>
              </g>
            );
          })}

          {/* axis */}
          <line x1={PAD} y1={Y_AXIS} x2={W - PAD} y2={Y_AXIS} stroke={CHROME.rule} />
          {ticks.map((t) => (
            <g key={t}>
              <line x1={x(t)} y1={Y_AXIS - 4} x2={x(t)} y2={Y_AXIS + 4} stroke={CHROME.rule} />
              <text x={x(t)} y={Y_AXIS + 16} textAnchor="middle" className="axis-label"
                fill={AXIS.fill} fontSize={AXIS.fontSize}>{t}</text>
            </g>
          ))}
        </svg>

        {/* per-boundary detail rows */}
        <div className="bc-detail-rows">
          {cells.map((c, i) => {
            const info = classInfo(c.boundary_class);
            return (
              <div key={i} className="bc-detail-row">
                <span className="bc-detail-type">{boundaryTypeLabel(c.boundary_type)}</span>
                <span className="bc-detail-badge">
                  <span className="bc-swatch" style={{ background: info.color }} /> {info.label}
                </span>
                <span className="bc-detail-meta">
                  AA {c.boundary_aa} · {c.distance_to_nearest_domain_boundary} aa from {c.nearest_domain_label}
                  {c.nearest_domain_class ? ` (${c.nearest_domain_class})` : ""}
                </span>
              </div>
            );
          })}
        </div>

        {!pd && (
          <div className="arch-note info">
            Cassette boundaries shown from the validated reference coordinates; full
            per-protein domain track is available in the Domain Architecture tab.
          </div>
        )}
      </div>
    </div>
  );
}
