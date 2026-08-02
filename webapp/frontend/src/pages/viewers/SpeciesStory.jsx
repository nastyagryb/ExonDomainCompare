import { useMemo, useState } from "react";
import { fileUrl } from "../../api";
import { Badge, Spinner, Empty } from "../../ui";
import { readinessLabel } from "../../uiStatus";
import { useIndex } from "./common";

export default function SpeciesStory({ preloaded, species, isoform, compact }) {
  const { data, loading } = useIndex((client) => client.story(), preloaded);
  const [sel, setSel] = useState(species || null);
  const [panel, setPanel] = useState(isoform || "IIIb");
  // A species handed down from the page (deep link, cross-view jump) wins over the
  // local dropdown choice. React's documented way to adjust state on a changed prop
  // is during render, not from an effect.
  const [syncedSpecies, setSyncedSpecies] = useState(species);
  if (species && species !== syncedSpecies) {
    setSyncedSpecies(species);
    setSel(species);
  }

  const speciesList = useMemo(() => (data?.species || []), [data]);
  const current = useMemo(
    () => speciesList.find((s) => s.species === (sel || species)) || speciesList[0],
    [speciesList, sel, species]
  );

  if (loading) return <Spinner label="Loading evidence story…" />;
  if (!data?.available || !current) return <Empty title="Evidence story not available" />;

  const panels = Object.keys(current.panels || {});
  const pd = current.panels?.[panel] || current.panels?.[panels[0]];

  return (
    <div className="viewer story-viewer">
      <div className="viewer-controls">
        {!compact && (
          <select value={current.species} onChange={(e) => setSel(e.target.value)}>
            {speciesList.map((s) => <option key={s.species} value={s.species}>{s.display_species_name}</option>)}
          </select>
        )}
        <div className="seg">
          {panels.map((p) => (
            <button key={p} className={`seg-btn iso-tint-${p.toLowerCase()}${(pd?.isoform === p) ? " on" : ""}`} onClick={() => setPanel(p)}>{p}</button>
          ))}
        </div>
        <span className="spacer" />
        {pd && (pd.is_review
          ? <Badge cls="review">Review / supplement</Badge>
          : <Badge cls={pd.overall}>{pd.rescued ? "Accepted · rescued" : readinessLabel(pd.overall)}</Badge>)}
      </div>

      {pd && (
        <>
          {pd.rescued && !pd.is_review && (
            <div className="story-note accepted">Rescued &amp; validated — accepted as a primary row. Provenance retained (see Rescue step).</div>
          )}
          {pd.is_review && pd.review_explanation && (
            <div className="story-note review">{pd.review_explanation}</div>
          )}
          <ol className="story-timeline">
            {pd.steps.map((s) => (
              <li key={s.key} className={`story-step st-${s.class}`}>
                <span className="story-dot" />
                <div className="story-card">
                  <div className="story-top">
                    <b>{s.title}</b>
                    <Badge cls={s.class} soft>{labelFor(s.class)}</Badge>
                  </div>
                  <p className="story-text">{s.text}</p>
                  <div className="story-meta">
                    {Object.entries(s.ids || {}).filter(([, v]) => v !== "" && v != null).slice(0, 4).map(([k, v]) => (
                      <span key={k} className="id-chip"><em>{k.replace(/_/g, " ")}</em> <code>{String(v)}</code></span>
                    ))}
                  </div>
                  <div className="story-links">
                    {s.source_table && <a className="btn ghost sm" href={fileUrl(s.source_table)}>Source table</a>}
                    {s.figure && <span className="fig-ref">{s.figure}</span>}
                  </div>
                </div>
              </li>
            ))}
          </ol>
        </>
      )}
    </div>
  );
}

function labelFor(cls) {
  return { accepted: "pass", minor: "minor", review: "review", excluded: "fail", neutral: "info" }[cls] || cls;
}
