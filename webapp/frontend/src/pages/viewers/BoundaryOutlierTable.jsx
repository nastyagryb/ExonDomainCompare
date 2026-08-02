import { Badge } from "../../ui";
import { classInfo, boundaryTypeLabel } from "./boundary";

// "Inspection cases" — display-coordinate confidence flags, never biological
// failures. Click a row to open its protein-level detail track.
export default function BoundaryOutlierTable({ outliers, onSelect }) {
  if (!outliers?.available) return null;
  const rows = outliers.outliers || [];
  return (
    <div className="bc-outliers">
      <div className="card-head">
        <h3>Inspection cases</h3>
        <Badge cls="minor" soft>{rows.length} display-confidence flags</Badge>
      </div>
      <p className="muted sm">{outliers.note}</p>
      <div className="heatmap-wrap">
        <table className="heatmap bc-outlier-table">
          <thead>
            <tr>
              <th className="sticky-col" style={{ textAlign: "left" }}>Species</th>
              <th>Iso</th>
              <th>Boundary</th>
              <th>Distance</th>
              <th>Nearest domain</th>
              <th>Display status</th>
              <th style={{ textAlign: "left" }}>Interpretation</th>
              <th></th>
            </tr>
          </thead>
          <tbody>
            {rows.map((o, i) => {
              const info = classInfo(o.boundary_class);
              return (
                <tr key={i}>
                  <td className="sticky-col" style={{ textAlign: "left" }}>{o.display_species_name}</td>
                  <td><span className={`iso iso-${o.isoform.toLowerCase()}`}>{o.isoform}</span></td>
                  <td>{boundaryTypeLabel(o.boundary_type)}</td>
                  <td>
                    <span className="bc-swatch" style={{ background: info.color }} /> {o.distance} aa
                  </td>
                  <td>{o.nearest_domain_label}</td>
                  <td><Badge cls="minor" soft>{(o.exon_block_display_status || "").replaceAll("_", " ")}</Badge></td>
                  <td style={{ textAlign: "left", maxWidth: 320 }}><small>{o.interpretation}</small></td>
                  <td><button className="btn ghost sm" onClick={() => onSelect?.(o.link_target)}>Detail →</button></td>
                </tr>
              );
            })}
            {rows.length === 0 && <tr><td colSpan={8} className="muted pad">No inspection cases.</td></tr>}
          </tbody>
        </table>
      </div>
    </div>
  );
}
