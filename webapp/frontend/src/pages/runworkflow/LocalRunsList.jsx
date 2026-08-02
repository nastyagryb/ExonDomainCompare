import { Badge } from "../../ui";

function fmtDate(s) {
  if (!s) return "—";
  const d = new Date(s);
  return Number.isNaN(d.getTime()) ? s : d.toLocaleString();
}

const STATUS_CLS = {
  created: "neutral", complete: "accepted", running: "minor",
  failed: "excluded", queued: "neutral", unknown: "neutral",
};

export default function LocalRunsList({
  runs = [], selectedRunId, onSelect, onLoadExample, exampleActive, loading,
}) {
  return (
    <div className="runs-list">
      <div className="runs-list-head">
        <h3>Runs</h3>
        {loading && <span className="muted small">refreshing…</span>}
      </div>

      {/* Example / freeze dataset — always available, read-only */}
      <div className={`run-row example${exampleActive ? " selected" : ""}`}>
        <div className="run-row-main">
          <div className="run-row-title">
            <b>Example: FGFR2 30-species freeze</b>
            <Badge cls="info" soft>read-only example</Badge>
          </div>
          <code className="run-row-path">results/final_30_until_interpro_prepare/</code>
        </div>
        <button className="btn ghost" onClick={onLoadExample}>Open example</button>
      </div>

      {runs.length === 0 ? (
        <div className="runs-empty muted">
          No local runs yet. Create one with the panel on the left.
        </div>
      ) : (
        runs.map((r) => (
          <div
            key={r.run_id}
            className={`run-row${selectedRunId === r.run_id ? " selected" : ""}`}
            onClick={() => onSelect?.(r.run_id)}
          >
            <div className="run-row-main">
              <div className="run-row-title">
                <b>{r.run_name || r.run_id}</b>
                <Badge cls={STATUS_CLS[r.status] || "neutral"} soft>{r.status || "unknown"}</Badge>
              </div>
              <div className="run-row-meta">
                <span>{r.run_id}</span>
                <span>· {r.species_count ?? "—"} species</span>
                <span>· created {fmtDate(r.created_at)}</span>
              </div>
            </div>
            <button
              className={`btn ${selectedRunId === r.run_id ? "primary" : "ghost"}`}
              onClick={(e) => { e.stopPropagation(); onSelect?.(r.run_id); }}
            >
              {selectedRunId === r.run_id ? "Selected" : "Open"}
            </button>
          </div>
        ))
      )}
    </div>
  );
}
