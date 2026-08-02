import { useMemo } from "react";
import { CANON_CLASS_COLOR, CANON_CLASS_LABEL, canonClass } from "./boundaryClasses";
import { isSupported, speciesTag } from "./comparativeFigures";
// Structural colours come from the shared figure palette, never from literals here:
// the exported figure reads the same tokens, and a second copy of a muted grey in a
// component is how the screen and the paper drifted apart before.
import { PALETTE } from "./figureSpec";
import { orderSpeciesIds, speciesCompare } from "./speciesOrder";

// Paired dot (dumbbell) plot of signed boundary-to-domain-edge distances.
//
// One row per comparable-boundary group, one dot per species observation. With the two
// species of the current datasets the two raw observations and the gap between them are
// the entire result, so this is a paired plot rather than a boxplot: a box drawn from
// n = 2 would show quartiles that do not exist.
//
// A connector between two dots asserts "this is the same boundary in both species".
// It is therefore drawn solid only when the backend called the mapping supported;
// tentative groups get a dotted hint, which states that the positions are close without
// claiming equivalence. The same rule is implemented once, in comparativeFilters.js,
// which publishes `connectable` per group so this view and the exported figure cannot
// disagree about which pairs are confirmed.

const PLOT_W = 860;
const ROW_H = 22;
const PAD = { top: 26, right: 34, bottom: 54, left: 78 };

const signedLabel = (v) => (v > 0 ? `+${v}` : String(v));

export default function ComparativePairedPlot({
  groups = [], stats = [], threshold = 5, selectedGroupId = null,
  onSelectGroup, height,
}) {
  const rows = useMemo(() => {
    const statById = new Map(stats.map((s) => [s.comparable_boundary_group_id, s]));
    return groups.map((g) => ({
      id: g.comparable_boundary_group_id,
      mappingStatus: g.mapping_status,
      supported: isSupported(g.mapping_status),
      connectable: Boolean(g.connectable),
      msaColumn: g.msa_column,
      stat: statById.get(g.comparable_boundary_group_id) || {},
      obs: (g.per_species_native_positions || [])
        .filter((o) => o.signed_distance != null)
        .slice()
        .sort((a, b) => speciesCompare(a.species_id, b.species_id)),
    })).filter((r) => r.obs.length);
  }, [groups, stats]);

  if (!rows.length) {
    return (
      <div className="arch-note info">
        No comparable boundary with a measurable distance matches the active filters.
      </div>
    );
  }

  // Species are identified by a fixed lane inside each row, exactly as in the exported
  // figure. Labelling every dot instead puts two labels on top of the markers they
  // explain whenever the species agree, which is the common case in real data.
  const speciesOrder = orderSpeciesIds(rows.flatMap((r) => r.obs.map((o) => o.species_id)));
  const speciesName = new Map();
  for (const r of rows) {
    for (const o of r.obs) speciesName.set(o.species_id, o.scientific_name);
  }

  const nSpeciesMax = Math.max(...rows.map((r) => r.obs.length));
  // The summary annotation grows with the species count without changing the plot type:
  // raw observations always stay visible.
  const summaryMode = nSpeciesMax <= 2 ? "pair" : nSpeciesMax <= 4 ? "range" : "box";

  const maxAbs = Math.max(10, ...rows.flatMap((r) => r.obs.map(
    (o) => Math.abs(Number(o.signed_distance)))));
  const lim = Math.ceil(maxAbs / 10) * 10;
  const innerW = PLOT_W - PAD.left - PAD.right;
  const bodyH = rows.length * ROW_H;
  const H = height || bodyH + PAD.top + PAD.bottom;
  const x = (v) => PAD.left + ((v + lim) / (2 * lim)) * innerW;

  const ticks = [];
  const step = lim <= 20 ? 5 : lim <= 60 ? 10 : 20;
  for (let t = -lim; t <= lim; t += step) ticks.push(t);

  const classesPresent = [...new Set(rows.flatMap((r) => r.obs.map(
    (o) => canonClass(o.boundary_class))))];

  return (
    <div className="cbe-plot-wrap">
      <svg className="cbe-paired-plot" viewBox={`0 0 ${PLOT_W} ${H}`}
        role="img"
        aria-label="Paired signed distances from comparable exon boundaries to the nearest representative-domain edge">
        {/* near-edge band and the domain-edge zero line */}
        <rect x={x(-threshold)} y={PAD.top - 4} width={x(threshold) - x(-threshold)}
          height={bodyH + 8} fill={PALETTE.grid} />
        <text x={x(0)} y={PAD.top - 8} textAnchor="middle" fontSize="9" fill={PALETTE.muted}>
          ±{threshold} aa
        </text>
        <line x1={x(0)} y1={PAD.top - 4} x2={x(0)} y2={PAD.top + bodyH + 4}
          stroke={PALETTE.ink} strokeWidth="1" />

        {rows.map((r, ri) => {
          const cy = PAD.top + ri * ROW_H + ROW_H / 2;
          const sel = r.id === selectedGroupId;
          const xs = r.obs.map((o) => x(Number(o.signed_distance)));
          const lo = Math.min(...xs);
          const hi = Math.max(...xs);
          const med = r.stat.median_signed_distance;
          return (
            <g key={r.id} className={`cbe-plot-row${sel ? " sel" : ""}`}
              onClick={() => onSelectGroup?.(r.id, r.obs[0])}
              style={{ cursor: "pointer" }}>
              {sel && (
                <rect x={PAD.left - 74} y={cy - ROW_H / 2} width={PLOT_W - PAD.left + 60}
                  height={ROW_H} fill={PALETTE.grid} opacity="0.6" />
              )}
              <text x={PAD.left - 8} y={cy + 3.5} textAnchor="end" fontSize="10"
                fill={sel ? PALETTE.ink : PALETTE.muted}
                fontWeight={sel ? "600" : "400"}>
                {String(r.id).replace(/^CBG/, "CBG ")}
              </text>

              {summaryMode === "pair" && r.obs.length === 2 && (
                r.connectable ? (
                  <line x1={lo} y1={cy} x2={hi} y2={cy} stroke={PALETTE.boundary}
                    strokeWidth="1.6" opacity="0.85" />
                ) : (
                  <line x1={lo} y1={cy} x2={hi} y2={cy} stroke={PALETTE.muted}
                    strokeWidth="1" strokeDasharray="2,2" opacity="0.9" />
                )
              )}
              {summaryMode !== "pair" && (
                <>
                  <line x1={lo} y1={cy} x2={hi} y2={cy} stroke={PALETTE.axis} strokeWidth="1" />
                  {med != null && (
                    <line x1={x(med)} y1={cy - 6} x2={x(med)} y2={cy + 6}
                      stroke={PALETTE.ink} strokeWidth="1.4" />
                  )}
                </>
              )}

              {r.obs.map((o, oi) => {
                const cx = x(Number(o.signed_distance));
                const dy = r.obs.length > 1
                  ? (oi - (r.obs.length - 1) / 2) * Math.min(5, ROW_H / (r.obs.length + 1))
                  : 0;
                const colour = CANON_CLASS_COLOR[canonClass(o.boundary_class)];
                const rad = sel ? 4.6 : 3.6;
                const title = `${o.scientific_name || o.species_id} · ${o.exon_transition}`
                  + ` · native aa ${o.native_position}`
                  + (o.msa_column != null ? ` · column ${o.msa_column}` : "")
                  + ` · ${signedLabel(Number(o.signed_distance))} aa to the `
                  + `${o.nearest_edge} edge of `
                  + `${o.nearest_domain_label || "no annotated domain"}`
                  + ` · ${CANON_CLASS_LABEL[canonClass(o.boundary_class)]}`
                  + ` · ${o.mapping_status === "tentative" ? "tentative" : "supported"} `
                  + `mapping, present in `
                  + `${Math.round((o.mapping_confidence ?? 0) * 100)}% of analysed species`;
                return (
                  <g key={`${r.id}:${o.species_id}`}>
                    <title>{title}</title>
                    {/* Open marker: measured to a domain start edge. Filled: end edge. */}
                    <circle cx={cx} cy={cy + dy} r={rad}
                      fill={o.nearest_edge === "start" ? PALETTE.paper : colour}
                      stroke={o.nearest_edge === "start" ? colour
                        : (sel ? PALETTE.ink : "none")}
                      strokeWidth={o.nearest_edge === "start" ? 1.4 : 0.8} />
                    {r.obs.length <= 2 && (
                      <text x={cx + rad + 3} y={cy + dy + 3} fontSize="8.5" fill={PALETTE.muted}>
                        {speciesTag(o.species_id, o.scientific_name).split(" ")[0]}
                      </text>
                    )}
                  </g>
                );
              })}
            </g>
          );
        })}

        {/* axis */}
        <line x1={PAD.left} y1={PAD.top + bodyH + 10} x2={PLOT_W - PAD.right}
          y2={PAD.top + bodyH + 10} stroke={PALETTE.axis} strokeWidth="1" />
        {ticks.map((t) => (
          <g key={t}>
            <line x1={x(t)} y1={PAD.top + bodyH + 10} x2={x(t)} y2={PAD.top + bodyH + 14}
              stroke={PALETTE.axis} strokeWidth="1" />
            <text x={x(t)} y={PAD.top + bodyH + 25} textAnchor="middle" fontSize="9"
              fill={PALETTE.muted}>{signedLabel(t)}</text>
          </g>
        ))}
        <text x={PAD.left + innerW / 2} y={PAD.top + bodyH + 42} textAnchor="middle"
          fontSize="10" fill={PALETTE.ink}>
          Signed distance to nearest representative-domain edge (aa) · 0 = domain edge
        </text>
      </svg>

      <div className="legend res-legend">
        {classesPresent.map((c) => (
          <span key={c} className="legend-item">
            <span className="pa-swatch" style={{ background: CANON_CLASS_COLOR[c] }} />
            {CANON_CLASS_LABEL[c]}
          </span>
        ))}
      </div>
      {speciesOrder.length > 1 && (
        <div className="muted sm">
          Marker lane per species —{" "}
          {speciesOrder.map((sid, i) => (
            <span key={sid}>
              {i > 0 && " · "}
              {i === 0 ? "upper" : i === speciesOrder.length - 1
                ? "lower" : `lane ${i + 1}`}:{" "}
              <i>{speciesName.get(sid) || sid}</i>
            </span>
          ))}
        </div>
      )}
      <div className="muted sm">
        {summaryMode === "pair"
          ? "Solid connector: supported cross-species mapping. Dotted connector: "
            + "tentative mapping — the positions are close but equivalence is not "
            + "established, so the two dots are not a confirmed pair."
          : "Bar: observed range across species. Vertical tick: median. Raw "
            + "observations are always shown."}
        {" "}Open marker: distance measured to a domain start edge; filled marker: to a
        domain end edge. Negative values lie upstream of that edge.
      </div>
    </div>
  );
}
