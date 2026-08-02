import { Kpi } from "../../ui";
import { BOUNDARY_CLASS } from "./boundary";

// Compact scientific summary of the boundary-consistency analysis.
export default function BoundarySummaryCards({ summary }) {
  if (!summary?.available) return null;
  const c = summary.boundary_class_counts || {};
  const total = summary.total_primary_proteins ?? 0;

  const coverage = [
    ["Primary proteins", total, "vertebrate FGFR2 panel", ""],
    ["Cassette data", `${summary.proteins_with_cassette_data ?? 0}/${total}`, "cassette coordinates", "accepted"],
    ["InterPro domains", `${summary.proteins_with_interpro_domain_data ?? 0}/${total}`, "domain annotation", "accepted"],
    ["pyTMHMM TM", `${summary.proteins_with_tm_data ?? 0}/${total}`, "transmembrane layer", "accepted"],
    ["Missing / review", c.review_or_missing ?? 0, "no coordinate", (c.review_or_missing ? "review" : "accepted")],
    ["Between domains", c.between_domains ?? 0, "cassette in linker", (c.between_domains ? "minor" : "accepted")],
  ];

  const classes = [
    ["aligned_to_domain_boundary", c.aligned_to_domain_boundary ?? 0],
    ["near_domain_boundary", c.near_domain_boundary ?? 0],
    ["within_domain", c.within_domain ?? 0],
  ];

  return (
    <div className="bc-summary">
      <div className="kpi-grid">
        {coverage.map(([l, v, s, cls]) => <Kpi key={l} label={l} value={v} sub={s} cls={cls} />)}
      </div>

      <div className="bc-class-strip">
        {classes.map(([key, n]) => {
          const info = BOUNDARY_CLASS[key];
          return (
            <div key={key} className="bc-class-chip" title={info.tip}>
              <span className="bc-swatch" style={{ background: info.color }} />
              <strong>{n}</strong>
              <span className="bc-class-label">{info.label}</span>
            </div>
          );
        })}
        <div className="bc-class-note">
          n = {summary.n_cassette_boundaries ?? 0} cassette boundaries · median{" "}
          {summary.median_distance ?? "—"} aa to nearest domain boundary
        </div>
      </div>

      {summary.key_interpretation && (
        <p className="bc-interpretation">{summary.key_interpretation}</p>
      )}
    </div>
  );
}
