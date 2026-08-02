import { api } from "../../api";
import { Badge, Spinner } from "../../ui";
import { RUNNING_STATUSES, statusBadgeCls } from "./runStatus";

// What a card offers depends on where the run stands, not on how much internal
// detail exists. A finished run needs a way in; a broken one needs a way
// forward; neither needs a stage checklist or a terminal log.

function formatDate(iso) {
  if (!iso) return "";
  const d = new Date(iso);
  if (Number.isNaN(d.getTime())) return "";
  return d.toLocaleDateString(undefined, { day: "numeric", month: "long", year: "numeric" });
}

function Meta({ run }) {
  const bits = [
    run.gene_symbol,
    run.species_count ? `${run.species_count} species` : "",
    formatDate(run.created_at),
  ].filter(Boolean);
  return (
    <div className="run-card-meta muted small">
      <span>{bits.join(" · ")}</span>
      {run.species_summary && <span className="run-card-species">{run.species_summary}</span>}
      <span className="run-card-id">Run ID: <code>{run.run_id}</code></span>
    </div>
  );
}

export default function RunCard({ run, selected, busyAction, onSelect, onExplore,
                                  onStop, onDelete, onRetry }) {
  const status = run.status;
  const group = run.group || "completed";
  const done = status === "results_ready";
  const attention = group === "attention";

  return (
    <div className={`run-card${selected ? " selected" : ""} run-card-${group}`}>
      <div className="run-card-main" onClick={() => onSelect?.(run.run_id)}>
        <div className="run-card-title">
          <b>{run.display_name || run.run_id}</b>
          <Badge cls={statusBadgeCls(status)} soft>{run.status_label || status}</Badge>
        </div>
        <Meta run={run} />
        {done && run.completion_summary && (
          <p className="run-card-summary muted small">{run.completion_summary}</p>
        )}
        {attention && run.failure_summary && (
          <p className="run-card-summary warn small">{run.failure_summary}</p>
        )}
        {group === "active" && (
          <div className="run-card-progress">
            <Spinner label={run.current_step || run.status_label || "Working…"} />
          </div>
        )}
      </div>

      <div className="run-card-actions">
        {run.explorable && (
          <button className="btn primary small" onClick={() => onExplore?.(`run:${run.run_id}`)}>
            Open results
          </button>
        )}
        {!run.explorable && group === "active" && (
          <button className="btn ghost small" onClick={() => onSelect?.(run.run_id)}>
            View status
          </button>
        )}
        {attention && onRetry && run.core_only && (
          <button className="btn ghost small" onClick={() => onRetry(run.run_id)}>Retry</button>
        )}
        {attention && (
          <a className="btn ghost small" href={api.runDiagnosticsUrl(run.run_id)}
             download>Download diagnostics</a>
        )}
        {RUNNING_STATUSES.has(status) && (
          <button className="btn ghost small" disabled={busyAction === "stop"}
                  onClick={() => onStop?.(run.run_id)}>Cancel</button>
        )}
        <button className="btn ghost small danger" disabled={busyAction === "delete"}
                onClick={() => onDelete?.(run.run_id)}>Delete</button>
      </div>
    </div>
  );
}
