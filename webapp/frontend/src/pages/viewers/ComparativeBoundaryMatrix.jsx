import { useState } from "react";
import { CANON_CLASS_COLOR, CANON_CLASS_LABEL, canonClass } from "./boundaryClasses";
import { MATRIX_STATE_FILL, isSupported, speciesTag } from "./comparativeFigures";
import { PALETTE } from "./figureSpec";

// Species × comparable-boundary-group matrix.
//
// The cells are the backend's observations, not a re-derivation: each cell carries the
// species' own boundary record, so the hover, the detail panel and the exported table
// are the same numbers. Cells without an observation are painted in an unsaturated
// grey/sand and carry no value — a zero would read as "sits exactly on a domain edge",
// which is a real class and the opposite of missing data.

const STATE_LABEL = {
  boundary_absent_or_unmapped: "No comparable boundary mapped in this species",
  result_pending: "Species analysis still pending",
  filtered_out: "Observation exists but is hidden by the active filters",
};

const cellFill = (state) => MATRIX_STATE_FILL[state]
  || CANON_CLASS_COLOR[canonClass(state)] || MATRIX_STATE_FILL.boundary_absent_or_unmapped;

const DARK_CLASSES = ["exact_domain_edge", "near_domain_edge", "inside_domain"];

function cellText(cell, mode) {
  if (!cell.observed || mode === "class") return "";
  const v = mode === "absolute" ? cell.absolute_distance : cell.signed_distance;
  if (v == null) return "";
  return mode === "absolute" ? String(v) : (v > 0 ? `+${v}` : String(v));
}

/** Hover card with every field the scientific reading of a cell needs. */
function CellHoverCard({ hover }) {
  if (!hover) return null;
  const { row, cell } = hover;
  const o = cell.observation;
  return (
    <div className="cbe-hover-card" role="status">
      <div className="cbe-hover-head">
        <i>{row.scientific_name || row.species_id}</i>
        <span className="muted sm">{cell.comparable_boundary_group_id}</span>
      </div>
      {!cell.observed ? (
        <div className="muted sm">{STATE_LABEL[cell.state] || cell.state}</div>
      ) : (
        <dl className="cbe-hover-dl">
          <dt>Protein</dt><dd><code>{o?.protein_id || row.protein_id || "—"}</code></dd>
          <dt>Exon transition</dt><dd>{o?.exon_transition || "—"}</dd>
          <dt>Native position</dt><dd>{o?.native_position ?? "—"} aa</dd>
          <dt>Alignment column</dt><dd>{o?.msa_column ?? "not mapped"}</dd>
          <dt>Nearest domain</dt>
          <dd>{o?.nearest_domain_full_label || o?.nearest_domain_label
            || "no annotated domain nearby"}</dd>
          <dt>Nearest edge</dt>
          <dd>{o?.nearest_edge ? `${o.nearest_edge} edge` : "—"}
            {o?.nearest_edge_position != null ? ` at aa ${o.nearest_edge_position}` : ""}</dd>
          <dt>Signed distance</dt>
          <dd>{o?.signed_distance == null ? "—"
            : `${o.signed_distance > 0 ? "+" : ""}${o.signed_distance} aa`}</dd>
          <dt>Boundary class</dt>
          <dd>{CANON_CLASS_LABEL[canonClass(o?.boundary_class)]}</dd>
          <dt>Mapping method</dt><dd>{o?.mapping_method || cell.mapping_method || "—"}</dd>
          {/* Mapping confidence is the qualitative evidence status, not the coverage
              fraction: a group can be mapped in every species (coverage 100 %) and
              still be tentative, so reporting the number as "confidence" would read
              as certainty the evidence does not support. */}
          <dt>Mapping confidence</dt>
          <dd>{o?.mapping_status === "tentative"
            ? "tentative — columns close, equivalence not established"
            : "supported"}</dd>
          <dt>Species coverage</dt>
          <dd>{o?.mapping_confidence != null
            ? `${Math.round(o.mapping_confidence * 100)}% of analysed species` : "—"}</dd>
        </dl>
      )}
    </div>
  );
}

export default function ComparativeBoundaryMatrix({
  matrix = [], groups = [], mode = "signed", selectedGroupId = null,
  selectedSpeciesId = null, onSelectCell, threshold = 5,
}) {
  const [hover, setHover] = useState(null);

  if (!matrix.length || !groups.length) {
    return (
      <div className="arch-note info">
        No comparable boundaries match the active filters. Reset the filters to see the
        full comparative matrix.
      </div>
    );
  }

  return (
    <div className="cbe-matrix-wrap">
      <div className="table-scroll">
        <table className="bnd-heatmap cbe-matrix">
          <thead>
            <tr>
              <th className="cbe-corner">Species</th>
              {groups.map((g) => {
                const sel = g.comparable_boundary_group_id === selectedGroupId;
                return (
                  <th key={g.comparable_boundary_group_id}
                    className={sel ? "sel" : ""}
                    title={`${g.comparable_boundary_group_id} · ${g.mapping_method}`
                      + ` · ${g.mapping_status}`
                      + (g.msa_column != null ? ` · alignment column ${g.msa_column}` : "")}
                    onClick={() => onSelectCell?.(null, {
                      comparable_boundary_group_id: g.comparable_boundary_group_id })}>
                    {String(g.comparable_boundary_group_id).replace(/^CBG/, "")}
                    {/* A tentative column is marked in the header so the caveat is
                        visible without hovering every cell. */}
                    {!isSupported(g.mapping_status) && <sup title="tentative mapping">t</sup>}
                  </th>
                );
              })}
            </tr>
          </thead>
          <tbody>
            {matrix.map((row) => (
              <tr key={row.species_id}
                className={row.species_id === selectedSpeciesId ? "row-selected" : ""}>
                <td className="cbe-row-head">
                  <i>{speciesTag(row.species_id, row.scientific_name)}</i>
                  {row.taxonomic_group && (
                    <span className="muted sm"> · {row.taxonomic_group}</span>)}
                </td>
                {(row.cells || []).map((cell) => {
                  const cls = canonClass(cell.state);
                  const dark = cell.observed && DARK_CLASSES.includes(cls);
                  const sel = cell.comparable_boundary_group_id === selectedGroupId;
                  const tentative = cell.observed && !isSupported(cell.mapping_status);
                  return (
                    <td key={cell.comparable_boundary_group_id}
                      className={`bnd-heat-cell${sel ? " sel" : ""}`
                        + `${tentative ? " cbe-tentative" : ""}`}
                      style={{ background: cellFill(cell.state),
                        color: dark ? PALETTE.paper : PALETTE.ink }}
                      tabIndex={0}
                      onMouseEnter={() => setHover({ row, cell })}
                      onFocus={() => setHover({ row, cell })}
                      onMouseLeave={() => setHover(null)}
                      onBlur={() => setHover(null)}
                      onClick={() => onSelectCell?.(cell.observation, cell)}
                      onKeyDown={(e) => {
                        if (e.key === "Enter" || e.key === " ") {
                          e.preventDefault();
                          onSelectCell?.(cell.observation, cell);
                        }
                      }}>
                      {cellText(cell, mode)}
                    </td>
                  );
                })}
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <CellHoverCard hover={hover} />

      <div className="legend res-legend">
        {[...new Set(matrix.flatMap((r) => (r.cells || [])
          .filter((c) => c.observed).map((c) => canonClass(c.state))))].map((c) => (
          <span key={c} className="legend-item">
            <span className="pa-swatch" style={{ background: CANON_CLASS_COLOR[c] }} />
            {CANON_CLASS_LABEL[c]}
          </span>
        ))}
        <span className="legend-item">
          <span className="pa-swatch"
            style={{ background: MATRIX_STATE_FILL.boundary_absent_or_unmapped }} />
          No comparable boundary mapped
        </span>
        <span className="legend-item">
          <span className="pa-swatch"
            style={{ background: MATRIX_STATE_FILL.result_pending }} />
          Species analysis pending
        </span>
        <span className="legend-item cbe-legend-tentative">
          <span className="pa-swatch cbe-swatch-tentative" />
          Hatched border: tentative mapping
        </span>
      </div>
      <div className="muted sm">
        Cell value: {mode === "class" ? "class only"
          : mode === "absolute" ? "absolute distance"
            : "signed distance"} to the nearest representative-domain edge ·
        near-edge threshold ±{threshold} aa · columns are comparable-boundary groups,
        never exon rank.
      </div>
    </div>
  );
}
