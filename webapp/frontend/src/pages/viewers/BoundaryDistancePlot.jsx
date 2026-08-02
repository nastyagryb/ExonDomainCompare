import { useMemo, useState } from "react";
import { classInfo, boundaryTypeLabel, isLowConfidence } from "./boundary";
import { textProps } from "./semanticStyles";

// Interactive distance distribution: distance to nearest protein-domain boundary,
// grouped by isoform or taxon group, filterable by isoform / boundary type / taxon /
// inspection-only. Hover shows species, isoform, distance, nearest domain.
const TAXA = ["Primates", "Other mammals", "Birds", "Reptiles", "Amphibians", "Teleost fish"];

// Text ink and size come from the shared scientific specification and are written onto
// every label as explicit SVG attributes, so the plot stays legible without the
// component stylesheet.
const AXIS = textProps("axis");
const GROUP = textProps("axisEmphasis");
const SUB = textProps("muted");

export default function BoundaryDistancePlot({ matrix }) {
  const [isoFilter, setIsoFilter] = useState("both");     // both | IIIb | IIIc
  const [btFilter, setBtFilter] = useState("both");       // both | cassette_start | cassette_end
  const [groupBy, setGroupBy] = useState("isoform");      // isoform | taxon
  const [inspectOnly, setInspectOnly] = useState(false);

  const points = useMemo(() => {
    const out = [];
    for (const r of matrix?.rows || []) {
      if (isoFilter !== "both" && r.isoform !== isoFilter) continue;
      if (inspectOnly && !isLowConfidence(r.exon_block_display_status)) continue;
      for (const c of r.cells || []) {
        if (btFilter !== "both" && c.boundary_type !== btFilter) continue;
        if (c.distance_to_nearest_domain_boundary == null) continue;
        out.push({
          species: r.display_species_name,
          isoform: r.isoform,
          taxon: r.taxon_group,
          boundary_type: c.boundary_type,
          distance: c.distance_to_nearest_domain_boundary,
          nearest: c.nearest_domain_label,
          cls: c.boundary_class,
        });
      }
    }
    return out;
  }, [matrix, isoFilter, btFilter, inspectOnly]);

  const groups = groupBy === "isoform"
    ? ["IIIb", "IIIc"].filter((g) => points.some((p) => p.isoform === g))
    : TAXA.filter((g) => points.some((p) => p.taxon === g));
  const groupOf = (p) => (groupBy === "isoform" ? p.isoform : p.taxon);

  const maxD = Math.max(20, ...points.map((p) => p.distance));
  const W = 720, H = 300, padL = 46, padB = 62, padT = 14, padR = 12;
  const plotH = H - padB - padT, plotW = W - padL - padR;
  const bandW = plotW / Math.max(1, groups.length);
  const y = (d) => padT + plotH - (d / maxD) * plotH;

  const yTicks = [];
  const yStep = maxD <= 40 ? 10 : maxD <= 120 ? 30 : 60;
  for (let t = 0; t <= maxD; t += yStep) yTicks.push(t);

  const stats = useMemo(() => {
    const m = {};
    for (const g of groups) {
      const ds = points.filter((p) => groupOf(p) === g).map((p) => p.distance).sort((a, b) => a - b);
      if (!ds.length) continue;
      const med = ds[Math.floor(ds.length / 2)];
      m[g] = { median: med, n: ds.length };
    }
    return m;
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [points, groups, groupBy]);

  return (
    <div className="bc-distplot">
      <div className="bc-filters">
        <Seg label="Isoform" value={isoFilter} onChange={setIsoFilter}
          opts={[["both", "IIIb + IIIc"], ["IIIb", "IIIb"], ["IIIc", "IIIc"]]} />
        <Seg label="Boundary" value={btFilter} onChange={setBtFilter}
          opts={[["both", "Start + End"], ["cassette_start", "Start"], ["cassette_end", "End"]]} />
        <Seg label="Group by" value={groupBy} onChange={setGroupBy}
          opts={[["isoform", "Isoform"], ["taxon", "Taxon"]]} />
        <label className="check inline">
          <input type="checkbox" checked={inspectOnly} onChange={(e) => setInspectOnly(e.target.checked)} />
          <span>Inspection cases only</span>
        </label>
      </div>

      <svg className="bc-distplot-svg" viewBox={`0 0 ${W} ${H}`} preserveAspectRatio="xMidYMid meet">
        {/* y grid + ticks */}
        {yTicks.map((t) => (
          <g key={t}>
            <line x1={padL} y1={y(t)} x2={W - padR} y2={y(t)} stroke="#e3e8f0" strokeWidth="1" />
            <text x={padL - 8} y={y(t) + 3} textAnchor="end" className="axis-label"
              fill={AXIS.fill} fontSize={AXIS.fontSize}>{t}</text>
          </g>
        ))}
        <text x={12} y={padT + plotH / 2} className="axis-label"
          fill={AXIS.fill} fontSize={AXIS.fontSize}
          transform={`rotate(-90 12 ${padT + plotH / 2})`} textAnchor="middle">
          distance to nearest domain boundary (aa)
        </text>

        {groups.map((g, gi) => {
          const cx = padL + bandW * gi + bandW / 2;
          const gpts = points.filter((p) => groupOf(p) === g);
          return (
            <g key={g}>
              <text x={cx} y={H - padB + 20} textAnchor="middle" className="axis-label bc-group-x"
                fill={GROUP.fill} fontSize={GROUP.fontSize} fontWeight={GROUP.fontWeight}>{g}</text>
              <text x={cx} y={H - padB + 36} textAnchor="middle" className="bc-group-x-sub"
                fill={SUB.fill} fontSize={SUB.fontSize}>
                n={stats[g]?.n ?? 0} · med {stats[g]?.median ?? "—"}
              </text>
              {/* median line */}
              {stats[g] && (
                <line x1={cx - bandW * 0.32} y1={y(stats[g].median)} x2={cx + bandW * 0.32} y2={y(stats[g].median)}
                  stroke="#1c2433" strokeWidth="1.6" />
              )}
              {/* jittered points */}
              {gpts.map((p, i) => {
                const jitter = (((i * 73) % 100) / 100 - 0.5) * bandW * 0.5;
                const info = classInfo(p.cls);
                return (
                  <circle key={i} cx={cx + jitter} cy={y(p.distance)} r="4"
                    fill={info.color} stroke="#fff" strokeWidth="0.7" opacity="0.9">
                    <title>{`${p.species} ${p.isoform} · ${boundaryTypeLabel(p.boundary_type)}\n${p.distance} aa from ${p.nearest} · ${info.label}`}</title>
                  </circle>
                );
              })}
            </g>
          );
        })}
        <line x1={padL} y1={padT + plotH} x2={W - padR} y2={padT + plotH} stroke="#d2d9e6" />
      </svg>
      {points.length === 0 && <p className="muted pad">No boundaries match the current filters.</p>}
    </div>
  );
}

function Seg({ label, value, onChange, opts }) {
  return (
    <div className="bc-seg-group">
      <span className="bc-seg-label">{label}</span>
      <div className="seg">
        {opts.map(([v, l]) => (
          <button key={v} className={`seg-btn${value === v ? " on" : ""}`} onClick={() => onChange(v)}>{l}</button>
        ))}
      </div>
    </div>
  );
}
