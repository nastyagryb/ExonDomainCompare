import { Badge, Spinner } from "../ui";

// PART 1 — a calm homepage with exactly three actions:
//   1. Explore Example Dataset   2. Create New Run   3. Open My Runs
export default function StartPage({
  health, datasets, loading,
  onExploreExample, onCreateRun, onOpenRuns, onExploreDataset,
}) {
  const apiOk = health?.status === "ok";
  const runs = (datasets || []).filter((d) => d.kind === "run" && !d.bundled_example);
  const exampleAvailable = (datasets || []).some((d) => d.kind === "example") || health?.example_available;

  return (
    <section className="page start-page">
      <div className="start-hero">
        <h1>ExonDomainCompare</h1>
        <p>Annotation-aware comparative exon–protein analysis · FGFR2 IIIb/IIIc</p>
        {!apiOk && <Badge cls="excluded" soft>Backend offline — start the API on :8000</Badge>}
      </div>

      <div className="start-cards">
        <ActionCard
          title="Explore Example Dataset"
          desc="The validated 30-species FGFR2 thesis dataset."
          badge="Validated · read-only"
          accent="info"
          disabled={!apiOk || !exampleAvailable}
          note={!exampleAvailable ? "Example dataset not found." : ""}
          cta="Explore example"
          onClick={onExploreExample}
        />
        <ActionCard
          title="Create New Run"
          desc="Enter your own species and start the analysis with one click."
          badge="Custom species"
          accent="accepted"
          disabled={!apiOk}
          cta="Create new run"
          onClick={onCreateRun}
        />
        <ActionCard
          title="Open My Runs"
          desc="See your local runs, their status and next action."
          badge={runs.length ? `${runs.length} run${runs.length === 1 ? "" : "s"}` : "No runs yet"}
          accent="neutral"
          disabled={!apiOk}
          cta="Open my runs"
          onClick={onOpenRuns}
        />
      </div>

      {loading && <Spinner label="Loading dataset…" />}

      {runs.length > 0 && (
        <div className="start-runs">
          <h3>My runs</h3>
          <div className="start-run-list">
            {runs.slice(0, 6).map((r) => (
              <button
                key={r.id}
                className="start-run-row"
                onClick={() => (r.explorable ? onExploreDataset(r.id) : onOpenRuns())}
                title={r.explorable ? "Explore this run" : "Open in My Runs"}
              >
                <span className="start-run-name">{r.label}</span>
                <span className="start-run-meta muted small">
                  {r.species_count ? `${r.species_count} species · ` : ""}{r.status_label || r.status}
                </span>
                <Badge cls={r.status === "results_ready" ? "accepted" : r.status === "failed" ? "excluded" : "neutral"} soft>
                  {r.explorable ? "Explore →" : "Details"}
                </Badge>
              </button>
            ))}
          </div>
        </div>
      )}
    </section>
  );
}

function ActionCard({ title, desc, badge, accent, onClick, cta, disabled, note }) {
  return (
    <article className={`action-card${disabled ? " disabled" : ""}`}>
      <div className="action-top">
        <Badge cls={accent} soft>{badge}</Badge>
      </div>
      <h3>{title}</h3>
      <p>{desc}</p>
      {note && <small className="muted">{note}</small>}
      <button className="btn primary" onClick={onClick} disabled={disabled}>{cta}</button>
    </article>
  );
}
