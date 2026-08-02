import { useMemo, useState } from "react";
import { CANON_CLASS_COLOR, CANON_CLASS_LABEL, canonClass } from "./boundaryClasses";
import { boundaryProps, featureStyle, textProps } from "./semanticStyles";

// Zero-centred signed-distance plot for the GENERIC exon-domain boundary page.
// Each point is an exon boundary positioned by its signed distance (aa) to the
// nearest representative-domain edge (0 = domain edge). Colour = canonical class,
// marker direction = nearest edge (start ▶ / end ◀). This replaces the previous
// uninformative absolute-distance bars. It never touches the validated FGFR2
// Boundary Consistency explorer.

// Tick ink, row-label ink and the near-edge band come from the shared scientific
// specification and are written onto the marks as explicit SVG attributes, so the plot
// stays legible without the component stylesheet.
const TICK = textProps("axis");
const ROW_LABEL = textProps("muted");
const NEAR_BAND = featureStyle("boundary_near").fill;

const signedOf = (r) =>
  r.signed_distance_aa != null ? Number(r.signed_distance_aa)
    : (r.signed_distance != null ? Number(r.signed_distance) : 0);

export default function SignedDistancePlot({ rows, threshold = 5, selectedId, onSelect, sort = "position" }) {
  const [hover, setHover] = useState(null);
  const posOf = (r) => r.protein_position ?? r.boundary_position_aa ?? 0;
  const absOf = (r) => (r.absolute_distance ?? r.absolute_distance_aa ?? Math.abs(signedOf(r)));
  const data = useMemo(
    () => (rows || [])
      .filter((r) => signedOf(r) != null && !Number.isNaN(signedOf(r)))
      .map((r) => ({ ...r, _signed: signedOf(r), _class: canonClass(r.boundary_class || r.category || r.class) }))
      .sort((a, b) => (sort === "distance"
        ? absOf(a) - absOf(b)
        : posOf(a) - posOf(b))),
    [rows, sort],
  );
  if (!data.length) return <p className="muted sm">No signed-distance data for this dataset.</p>;

  const W = 900;
  const rowH = 30;
  const padTop = 22;
  const padBottom = 34;
  const padL = 20;
  const padR = 20;
  const H = padTop + padBottom + data.length * rowH;
  const maxAbs = Math.max(threshold + 2, ...data.map((r) => Math.abs(r._signed)));
  const dom = Math.ceil(maxAbs * 1.12);
  const x = (v) => padL + ((v + dom) / (2 * dom)) * (W - padL - padR);
  const yFor = (i) => padTop + i * rowH + rowH / 2;

  // gridline ticks at nice intervals
  const step = dom <= 20 ? 5 : dom <= 60 ? 20 : 50;
  const ticks = [];
  for (let t = -Math.floor(dom / step) * step; t <= dom; t += step) ticks.push(t);

  return (
    <div className="signed-dist">
      <svg viewBox={`0 0 ${W} ${H}`} className="signed-dist-svg" preserveAspectRatio="xMidYMid meet"
        role="img" aria-label="Signed distance of exon boundaries to nearest domain edge">
        {/* near-edge band — the ±threshold zone is the near_domain_edge class itself */}
        <rect x={x(-threshold)} y={padTop - 6} width={x(threshold) - x(-threshold)}
          height={H - padTop - padBottom + 12} fill={NEAR_BAND} opacity="0.22" />
        {/* ticks */}
        {ticks.map((t) => (
          <g key={t}>
            <line x1={x(t)} x2={x(t)} y1={padTop - 6} y2={H - padBottom + 4}
              stroke={t === 0 ? "#1F2933" : "#E4E7EB"} strokeWidth={t === 0 ? 1.4 : 1} />
            <text x={x(t)} y={H - padBottom + 18} textAnchor="middle" className="sd-tick"
              fill={TICK.fill} fontSize={TICK.fontSize}>{t}</text>
          </g>
        ))}
        {/* rows */}
        {data.map((r, i) => {
          const y = yFor(i);
          const id = r.id || r.exon_boundary_id;
          const sel = selectedId && id === selectedId;
          const col = boundaryProps(r._class).fill;
          const edge = (r.nearest_edge || r.nearest_edge_type || "").toLowerCase();
          const dir = edge === "start" ? "▶" : edge === "end" ? "◀" : "●";
          const lbl = r.label || (id || "").split(":").pop();
          return (
            <g key={id} className={`sd-row${sel ? " sel" : ""}`}
              onMouseEnter={() => setHover(r)} onMouseLeave={() => setHover(null)}
              onClick={() => onSelect && onSelect(r)} style={{ cursor: onSelect ? "pointer" : "default" }}>
              <line x1={x(0)} x2={x(r._signed)} y1={y} y2={y} stroke={col} strokeWidth={sel ? 2.4 : 1.3} opacity="0.7" />
              <text x={x(r._signed)} y={y + 4} textAnchor="middle" fontSize={sel ? 16 : 13} fill={col}
                stroke="#1F2933" strokeWidth="0.4">{dir}</text>
              <text x={padL} y={y - 8} className={`sd-label${sel ? " sel" : ""}`}
                fill={ROW_LABEL.fill} fontSize={ROW_LABEL.fontSize}>{lbl}</text>
            </g>
          );
        })}
      </svg>
      {hover && (
        <div className="sd-tip">
          <b>{hover.label || (hover.id || hover.exon_boundary_id || "").split(":").pop()}</b>
          {" · "}aa {hover.protein_position ?? hover.boundary_position_aa}
          {" · "}signed {hover._signed} aa (|{Math.abs(hover._signed)}|)
          {" · "}nearest {hover.nearest_domain_label || hover.nearest_domain_name || hover.nearest_domain_id || "—"}
          {" "}{hover.nearest_edge || hover.nearest_edge_type || ""}
          {" · "}<span style={{ color: CANON_CLASS_COLOR[hover._class] }}>{CANON_CLASS_LABEL[hover._class]}</span>
        </div>
      )}
      <div className="sd-legend">
        {Object.keys(CANON_CLASS_LABEL).map((c) => (
          <span key={c} className="sd-legend-item">
            <span className="sd-swatch" style={{ background: CANON_CLASS_COLOR[c] }} />
            {CANON_CLASS_LABEL[c]}
          </span>
        ))}
        <span className="muted sm">0 = domain edge · shaded band = ±{threshold} aa near-edge · ▶ start edge · ◀ end edge</span>
      </div>
    </div>
  );
}
