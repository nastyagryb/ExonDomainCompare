import { useCallback, useEffect, useRef, useState } from "react";
import { api } from "../api";
import { Badge, Drawer, Empty, Field } from "../ui";

const STEP_CLASS = {
  pass: "accepted", running: "info", pending: "neutral", queued: "neutral",
  warning: "minor", failed: "excluded", fail: "excluded", skipped: "neutral",
};

const NO_LOGS = [];

export default function RunMonitor({ activeRunId, setActiveRunId, summary, onRunFinished }) {
  const [runs, setRuns] = useState([]);
  const [polledStatus, setPolledStatus] = useState(null);
  const [polledLogs, setPolledLogs] = useState([]);
  // Both values belong to the poller below, so with no run being monitored the page
  // must not keep showing the last poll. Deriving that here rather than clearing it
  // from the effect keeps the poll effect free of synchronous state writes.
  const status = activeRunId ? polledStatus : null;
  const logs = activeRunId ? polledLogs : NO_LOGS;
  const [stepDrawer, setStepDrawer] = useState(null);
  const [showLogs, setShowLogs] = useState(false);
  const finishedRef = useRef(false);

  const refreshRuns = useCallback(() => { api.runs().then(setRuns).catch(() => setRuns([])); }, []);
  useEffect(() => { refreshRuns(); }, [refreshRuns]);

  // poll active run
  useEffect(() => {
    if (!activeRunId) return;
    finishedRef.current = false;
    let stop = false;
    async function poll() {
      try {
        const st = await api.runStatus(activeRunId);
        const lg = await api.runLogs(activeRunId, 600);
        if (stop) return;
        setPolledStatus(st);
        setPolledLogs(lg.lines || []);
        if ((st.status === "finished") && !finishedRef.current) {
          finishedRef.current = true;
          refreshRuns();
          onRunFinished?.();
        }
        if (st.status === "finished" || st.status === "failed") return; // stop polling
      } catch { /* keep trying */ }
      if (!stop) setTimeout(poll, 2000);
    }
    poll();
    return () => { stop = true; };
  }, [activeRunId, onRunFinished, refreshRuns]);

  // when nothing actively monitored, show the loaded run's steps
  const steps = status?.steps || summary?.steps || [];
  const runMode = status?.run_mode || summary?.run_mode;
  const isCached = runMode && (runMode.used_cached_v3_outputs || runMode.used_cached_msa_outputs);
  const headerStatus = status?.status || (summary ? "finished" : null);

  return (
    <section className="page">
      <div className="page-head">
        <div>
          <p className="eyebrow">Pipeline execution monitor</p>
          <h2>Run Monitor</h2>
        </div>
        <div className="head-badges">
          {headerStatus && <Badge cls={STEP_CLASS[headerStatus] || "neutral"} soft>{headerStatus}</Badge>}
          {isCached && <Badge cls="info" soft>cached / debug — not full clean</Badge>}
          {runMode?.full_clean_run_completed && <Badge cls="accepted" soft>full clean run</Badge>}
        </div>
      </div>

      <div className="monitor-layout">
        <aside className="run-list">
          <h4>Runs</h4>
          {runs.map((r) => (
            <button
              key={r.run_id}
              className={`run-item${(activeRunId === r.run_id) ? " sel" : ""}`}
              onClick={() => setActiveRunId(r.kind === "example" ? null : r.run_id)}
              disabled={r.kind === "example"}
            >
              <Badge cls={STEP_CLASS[r.status] || (r.kind === "example" ? "info" : "neutral")} soft>{r.status || r.kind}</Badge>
              <span>{r.label || r.run_id}</span>
              {r.species_count && <small>{r.species_count} species</small>}
            </button>
          ))}
          {runs.length === 0 && <p className="muted">No runs yet.</p>}
        </aside>

        <div className="timeline-wrap">
          {steps.length === 0 ? (
            <Empty title="No active run" hint="Start a run from the Start page, or load a run to see its step timeline." />
          ) : (
            <>
              <ol className="timeline">
                {steps.map((s) => {
                  const cls = STEP_CLASS[s.status] || "neutral";
                  return (
                    <li key={s.step_id} className={`tl-step st-${cls}`}>
                      <span className="tl-dot" />
                      <button className="tl-body" onClick={() => setStepDrawer(s)}>
                        <div className="tl-top">
                          <b>{s.label}</b>
                          <Badge cls={cls} soft>{s.status || "pending"}</Badge>
                        </div>
                        <div className="tl-sub">
                          <small>{s.step_id}</small>
                          {s.runtime_seconds != null && <small>{fmtRuntime(s.runtime_seconds)}</small>}
                          {s.output_files && <small className="mono">{s.output_files.split(/[;\s]/)[0]}</small>}
                          {s.warning_summary && <small className="warn">{s.warning_summary.slice(0, 50)}</small>}
                        </div>
                      </button>
                    </li>
                  );
                })}
              </ol>

              <div className="log-toggle">
                <button className="btn ghost sm" onClick={() => setShowLogs(!showLogs)}>
                  {showLogs ? "Hide" : "Show"} raw log
                </button>
              </div>
              {showLogs && (
                <pre className="log-pre">{logs.length ? logs.join("\n") : "No log lines yet."}</pre>
              )}
              {status?.status === "failed" && (
                <div className="run-error">
                  <Badge cls="excluded">failed</Badge>
                  <span>{status.error || "Pipeline failed. See raw log."}</span>
                </div>
              )}
            </>
          )}
        </div>
      </div>

      <Drawer
        open={Boolean(stepDrawer)}
        onClose={() => setStepDrawer(null)}
        title={stepDrawer?.label}
        subtitle={stepDrawer ? `${stepDrawer.step_id} · ${stepDrawer.step_name}` : ""}
      >
        {stepDrawer && (
          <>
            <div className="drawer-badges">
              <Badge cls={STEP_CLASS[stepDrawer.status] || "neutral"}>{stepDrawer.status || "pending"}</Badge>
              {stepDrawer.return_code != null && <Badge cls={stepDrawer.return_code === 0 ? "accepted" : "excluded"} soft>rc={stepDrawer.return_code}</Badge>}
            </div>
            <Field label="Runtime">{stepDrawer.runtime_seconds != null ? fmtRuntime(stepDrawer.runtime_seconds) : "—"}</Field>
            <Field label="Command" wide><code className="block">{stepDrawer.command || "—"}</code></Field>
            <Field label="Output files" wide><code className="block">{stepDrawer.output_files || "—"}</code></Field>
            {stepDrawer.warning_summary && <Field label="Warnings" wide><span className="warn">{stepDrawer.warning_summary}</span></Field>}
          </>
        )}
      </Drawer>
    </section>
  );
}

function fmtRuntime(s) {
  if (s == null) return "—";
  if (s < 60) return `${s}s`;
  const m = Math.floor(s / 60), r = s % 60;
  return `${m}m ${r}s`;
}
