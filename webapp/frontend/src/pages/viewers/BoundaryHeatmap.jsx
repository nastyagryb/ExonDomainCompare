import { Fragment } from "react";
import { classInfo, boundaryTypeLabel, isLowConfidence, BOUNDARY_ORDER, BOUNDARY_CLASS } from "./boundary";

// Interactive boundary-consistency heatmap.
//   * rows ordered taxonomically, grouped with subtle separators + group labels
//   * columns: cassette_start / cassette_end (+ any extra relation columns present)
//   * cells coloured by boundary_class; hover tooltip; click opens detail track
export default function BoundaryHeatmap({ matrix, onSelect, selected }) {
  if (!matrix?.available) return null;
  const cols = matrix.columns || [];
  const rows = matrix.rows || [];

  // group boundaries by taxon_group (rows are already taxonomically ordered), so a
  // group label belongs on the first row of every run of equal groups
  const showsGroupLabel = rows.map((r, i) => r.taxon_group !== (i === 0 ? null : rows[i - 1].taxon_group));

  return (
    <div className="bc-heatmap-wrap">
      <table className="heatmap bc-heatmap">
        <thead>
          <tr>
            <th className="sticky-col">Species</th>
            <th>Iso</th>
            {cols.map((c) => <th key={c.key} title={c.label}>{c.label}</th>)}
          </tr>
        </thead>
        <tbody>
          {rows.map((r, rowIndex) => {
            const cellByType = Object.fromEntries((r.cells || []).map((c) => [c.boundary_type, c]));
            const lowConf = isLowConfidence(r.exon_block_display_status);
            const showGroup = showsGroupLabel[rowIndex];
            const key = `${r.species}-${r.isoform}`;
            const isSel = selected && selected.species === r.species && selected.isoform === r.isoform;
            return (
              <Fragment key={key}>
                {showGroup && (
                  <tr className="bc-group-row">
                    <td className="sticky-col bc-group-label" colSpan={2 + cols.length}>{r.taxon_group}</td>
                  </tr>
                )}
                <tr className={isSel ? "bc-row-sel" : ""}>
                  <td className="sticky-col species-name">
                    {r.display_species_name}
                    {lowConf && <span className="bc-star" title="Reconstructed / cassette-only display — inspection case">*</span>}
                  </td>
                  <td><span className={`iso iso-${r.isoform.toLowerCase()}`}>{r.isoform}</span></td>
                  {cols.map((col) => {
                    const cell = cellByType[col.key];
                    if (!cell) return <td key={col.key} className="hm-cell"><span className="bc-cell-empty" /></td>;
                    const info = classInfo(cell.boundary_class);
                    const tip = `${r.display_species_name} ${r.isoform} · ${boundaryTypeLabel(cell.boundary_type)}\n`
                      + `AA ${cell.boundary_aa} · ${info.label}\n`
                      + `nearest: ${cell.nearest_domain_label} (${cell.distance_to_nearest_domain_boundary} aa)`;
                    return (
                      <td key={col.key} className="hm-cell">
                        <button
                          className={`bc-cell${lowConf ? " low-conf" : ""}`}
                          style={{ background: info.color }}
                          title={tip}
                          onClick={() => onSelect?.(r, cell)}
                        >
                          <span className="bc-cell-dist">{cell.distance_to_nearest_domain_boundary}</span>
                        </button>
                      </td>
                    );
                  })}
                </tr>
              </Fragment>
            );
          })}
        </tbody>
      </table>
      <BoundaryLegend />
    </div>
  );
}

function BoundaryLegend() {
  return (
    <div className="arch-legend bc-legend">
      {BOUNDARY_ORDER.map((k) => (
        <span key={k} className="leg-item" title={BOUNDARY_CLASS[k].tip}>
          <span className="leg-swatch" style={{ background: BOUNDARY_CLASS[k].color }} />
          {BOUNDARY_CLASS[k].label}
        </span>
      ))}
      <span className="leg-item"><span className="bc-star">*</span>reconstructed / cassette-only display (inspect)</span>
      <span className="leg-item bc-legend-hint">cell value = distance to nearest domain boundary (aa)</span>
    </div>
  );
}
