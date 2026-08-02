import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { api } from "../api";
import { Badge, Spinner, Empty } from "../ui";
import CreateRunPanel from "./runworkflow/CreateRunPanel";
import RunStatusStepper from "./runworkflow/RunStatusStepper";
import CopyButton from "./runworkflow/CopyButton";
import RunCard from "./runworkflow/RunCard";
import { statusBadgeCls } from "./runworkflow/runStatus";
import { isActiveRunState, sortNewestFirst } from "../runStates";
import HumanReferenceBadge from "./HumanReferenceBadge";

/** How often the selected run's detail is re-read while it is still progressing. */
const RUN_POLL_MS = 5000;
const RUN_DISCOVERY_POLL_MS = 1000;
const RUN_DISCOVERY_ATTEMPTS = 30;
const RUNNING = new Set(["pre_interpro_running", "post_interpro_running", "running"]);
// Terminal pre-cluster failure states (shared/core runner + FGFR2 pipeline).
const FAILED = new Set(["failed", "core_model_collection_failed", "incomplete"]);

//: Runs are grouped by what they need from the reader, newest first inside each
//: group. The backend already returns them in order, so grouping preserves it.
const GROUPS = [
  ["active", "Active"],
  ["attention", "Attention required"],
  ["completed", "Completed"],
];

// Build a calm fallback detail from the run-list entry when status.json cannot be
// read yet (run is starting / status file mid-write). Never show a scary error.
function fallbackDetail(listed, err) {
  return {
    _runId: listed.run_id,
    _fallback: true,
    _error: err?.message || "",
    config: { run_name: listed.run_name, species_count: listed.species_count },
    status: {},
    model: {
      run_id: listed.run_id,
      status: listed.status || "created",
      status_label: listed.status_label || "Run is starting…",
      current_step: listed.current_step || "",
      next_action: listed.next_action || "",
      available_views: listed.available_views || {},
      human_reference: listed.human_reference || {},
      primary_fasta_count: listed.primary_fasta_count || 0,
    },
    files: {},
  };
}

export default function RunWorkflowPage({ activeId, onExploreDataset, onDatasetsChanged,
                                          onRefreshAll, onActiveDeleted }) {
  const [runs, setRuns] = useState([]);
  const [loadingRuns, setLoadingRuns] = useState(false);
  const [starting, setStarting] = useState(false);
  const [selectedId, setSelectedId] = useState(null);
  const [detail, setDetail] = useState(null);   // {config,status,model,files,_runId,_fallback}
  const [commands, setCommands] = useState(null);
  const [loadingDetail, setLoadingDetail] = useState(false);
  const [busy, setBusy] = useState(false);
  const [action, setAction] = useState("");     // "stop" | "delete" | "refresh" | ""
  const pollRef = useRef(null);
  const runsRef = useRef([]);

  const loadRuns = useCallback(async ({ preserveIds = [] } = {}) => {
    setLoadingRuns(true);
    try {
      const list = await api.localRuns();
      const persisted = Array.isArray(list) ? list : [];
      const persistedIds = new Set(persisted.map((run) => run.run_id));
      const pending = runsRef.current.filter(
        (run) => preserveIds.includes(run.run_id) && !persistedIds.has(run.run_id),
      );
      const arr = sortNewestFirst([...pending, ...persisted]);
      setRuns(arr);
      runsRef.current = arr;
      return { visible: arr, persisted };
    } catch {
      if (preserveIds.length === 0) {
        setRuns([]);
        runsRef.current = [];
      }
      return { visible: runsRef.current, persisted: [] };
    } finally {
      setLoadingRuns(false);
    }
  }, []);

  /**
   * Put a just-created run into the list straight away, newest first.
   *
   * The subsequent registry refetch replaces this entry with the persisted one, so the
   * optimistic row is a bridge across the create/list gap rather than a second source of
   * truth. It is keyed by run_id, so the refetch cannot produce a duplicate.
   */
  const insertRun = useCallback((created) => {
    const runId = created?.run_id;
    if (!runId) return;
    setRuns((prev) => {
      const rest = prev.filter((r) => r.run_id !== runId);
      const existing = prev.find((r) => r.run_id === runId) || {};
      const next = sortNewestFirst([{
        ...existing,
        run_id: runId,
        run_name: created.run_name || existing.run_name || "",
        gene_symbol: created.gene_symbol || existing.gene_symbol || "",
        species_count: created.species_count ?? existing.species_count,
        status: created.status || existing.status || "created",
        group: created.group || existing.group || "active",
      }, ...rest]);
      runsRef.current = next;
      return next;
    });
  }, []);

  const loadDetail = useCallback(async (runId, { full = false } = {}) => {
    if (!runId) return;
    setLoadingDetail(true);
    try {
      const [ref, cmds] = await Promise.all([
        full ? api.localRunRefreshAll(runId) : api.localRunRefresh(runId),
        api.localRunCommands(runId).catch(() => null),
      ]);
      setDetail({
        _runId: runId,
        config: ref.config, status: ref.status,
        model: ref.status_model, files: ref.files,
      });
      if (ref.summary) {
        setRuns((prev) => {
          const next = sortNewestFirst(prev.map((run) => (
            run.run_id === runId ? { ...run, ...ref.summary } : run
          )));
          runsRef.current = next;
          return next;
        });
      }
      if (cmds) setCommands(cmds);
      onDatasetsChanged?.();
      return ref;
    } catch (err) {
      // Never wipe a good detail on a transient poll error, and never show a
      // scary "unavailable" panel while the run exists: fall back to the
      // run-list summary so a starting/running run always has a useful panel.
      setDetail((prev) => {
        if (prev && prev._runId === runId && !prev._fallback) return prev;
        const listed = runsRef.current.find((r) => r.run_id === runId);
        return listed ? fallbackDetail(listed, err) : prev;
      });
      return null;
    } finally {
      setLoadingDetail(false);
    }
  }, [onDatasetsChanged]);

  useEffect(() => { loadRuns(); }, [loadRuns]);
  useEffect(() => { if (selectedId) loadDetail(selectedId); }, [selectedId, loadDetail]);

  // Auto-refresh until the run reaches a stable state. Polling only the three local
  // "running" states meant a run waiting for the external cluster round-trip was never
  // re-read, so its completion appeared only after a browser refresh. Every non-terminal
  // state is polled now, including cluster_required, which is exactly the state a run sits
  // in while the round-trip runs outside the app.
  useEffect(() => {
    const status = detail?.model?.status;
    if (selectedId && isActiveRunState(status)) {
      pollRef.current = setInterval(() => loadDetail(selectedId), RUN_POLL_MS);
      return () => clearInterval(pollRef.current);
    }
    if (pollRef.current) clearInterval(pollRef.current);
    return undefined;
  }, [selectedId, detail?.model?.status, loadDetail]);

  async function handleStart(payload) {
    setStarting(true);
    try {
      const knownIds = new Set(runsRef.current.map((run) => run.run_id));
      const res = await api.localRunStart(payload);
      // Show the run before asking the registry for it. The backend answers as soon as
      // the run is created, which can be a moment before it is fully listable, and
      // waiting for the list is what made a new run seem to need a browser refresh.
      if (res?.run_id) {
        insertRun(res);
        setSelectedId(res.run_id);
      }
      const preserveIds = res?.run_id ? [res.run_id] : [];
      for (let attempt = 0; attempt < RUN_DISCOVERY_ATTEMPTS; attempt += 1) {
        const { persisted } = await loadRuns({ preserveIds });
        const created = persisted.find((run) => (
          res?.run_id ? run.run_id === res.run_id : !knownIds.has(run.run_id)
        ));
        if (created) {
          setSelectedId(created.run_id);
          break;
        }
        await new Promise((resolve) => window.setTimeout(resolve, RUN_DISCOVERY_POLL_MS));
      }
      onDatasetsChanged?.();
    } finally {
      setStarting(false);
    }
  }

  async function runPostInterpro() {
    if (!selectedId) return;
    setBusy(true);
    try {
      await api.startPostInterpro(selectedId);
      await loadDetail(selectedId);
    } catch (e) {
      alert(e?.message || "Could not start Post-InterPro.");
    } finally {
      setBusy(false);
    }
  }

  // Full refresh: rebuild run-local indices, reload this run's detail, the run
  // list, and every dataset-scoped page (so Boundary Consistency etc. appear
  // without a browser reload).
  async function refreshEverything() {
    setAction("refresh");
    try {
      await loadRuns();
      if (selectedId) await loadDetail(selectedId, { full: true });
      await onRefreshAll?.(activeId);
    } finally {
      setAction("");
    }
  }

  async function handleStop(runId) {
    if (!runId) return;
    setAction("stop");
    try {
      const res = await api.stopLocalRun(runId);
      await loadDetail(runId);
      await loadRuns();
      onDatasetsChanged?.();
      if (res?.stopped === false) {
        alert(res?.cluster_note || "No local process was running to stop.");
      }
    } catch (e) {
      alert(e?.message || "Could not stop the run.");
    } finally {
      setAction("");
    }
  }

  async function handleRetry(runId) {
    if (!runId) return;
    setBusy(true);
    try {
      const res = await api.retryPreCluster(runId);
      // The retry rebuilds this run, so selection stays put; older non-in-place
      // launches still hand back a different id to follow.
      const nextId = res?.launch?.run_id;
      await loadRuns();
      onDatasetsChanged?.();
      if (nextId) setSelectedId(nextId);
      else alert(res?.launch?.note || "Retry started; the run is provisioning.");
    } catch (e) {
      alert(e?.message || "Could not retry this run.");
    } finally {
      setBusy(false);
    }
  }

  // Rebuilds the failed local stage of this run. The run id, name and species stay
  // as they were, so one request never turns into two runs.
  async function handleRetryLocal(runId) {
    if (!runId) return;
    setBusy(true);
    try {
      const res = await api.retryLocalPreparation(runId);
      await loadRuns();
      onDatasetsChanged?.();
      if (res?.note) alert(res.note);
    } catch (e) {
      alert(e?.message || "Could not restart local preparation for this run.");
    } finally {
      setBusy(false);
    }
  }

  async function handleDelete(runId) {
    if (!runId) return;
    const ok = window.confirm(
      "Delete this run and all run-local outputs? This cannot be undone.");
    if (!ok) return;
    setAction("delete");
    try {
      await api.deleteLocalRun(runId);
      if (selectedId === runId) { setSelectedId(null); setDetail(null); setCommands(null); }
      await loadRuns();
      onDatasetsChanged?.();
      if (activeId === `run:${runId}`) onActiveDeleted?.();
    } catch (e) {
      alert(e?.message || "Could not delete the run.");
    } finally {
      setAction("");
    }
  }

  const model = detail?.model;
  const files = detail?.files;
  const status = model?.status;
  const humanRef = model?.human_reference || detail?.config?.human_reference;
  const selectedRun = useMemo(
    () => runs.find((r) => r.run_id === selectedId) || null, [runs, selectedId]);
  const grouped = useMemo(() => GROUPS
    .map(([key, label]) => [key, label, runs.filter((r) => (r.group || "completed") === key)])
    .filter(([, , items]) => items.length > 0), [runs]);

  return (
    <section className="page run-workflow">
      <div className="page-head">
        <div>
          <h1>My Runs</h1>
          <p className="muted">
            Create a run with your own species and follow it here. The example dataset is always
            available and read-only; custom runs write only under <code>runs/&lt;run_id&gt;/</code>.
          </p>
        </div>
      </div>

      <div className="rw-grid">
        <CreateRunPanel onStart={handleStart} starting={starting} />

        <div className="local-runs card">
          <div className="lr-head">
            <h3>My runs</h3>
            <button className="btn ghost small" onClick={loadRuns} disabled={loadingRuns}>Refresh</button>
          </div>
          {loadingRuns && runs.length === 0 ? (
            <Spinner label="Loading runs…" />
          ) : runs.length === 0 ? (
            <p className="muted small">No custom runs yet. Create one on the left.</p>
          ) : (
            <div className="lr-list">
              {grouped.map(([key, label, items]) => (
                <div key={key} className="lr-group">
                  {grouped.length > 1 && <h4 className="lr-group-head">{label}</h4>}
                  {items.map((r) => (
                    <RunCard key={r.run_id} run={r} selected={selectedId === r.run_id}
                             busyAction={action}
                             onSelect={setSelectedId}
                             onExplore={onExploreDataset}
                             onStop={handleStop}
                             onDelete={handleDelete}
                             onRetry={handleRetry} />
                  ))}
                </div>
              ))}
            </div>
          )}
        </div>
      </div>

      {selectedId && (
        <div className="rw-detail card">
          {loadingDetail && !detail ? (
            <Spinner label="Loading run…" />
          ) : detail && model ? (
            <>
              <div className="rw-detail-head">
                <div>
                  <h2>{selectedRun?.display_name || detail.config?.run_name || selectedId}</h2>
                  <div className="rw-detail-meta muted">
                    <span>
                      {[selectedRun?.gene_symbol,
                        `${selectedRun?.species_count ?? detail.config?.species_count ?? "—"} species`,
                        selectedRun?.analysis_mode].filter(Boolean).join(" · ")}
                    </span>
                    <span>Run ID: <code>{selectedId}</code></span>
                  </div>
                </div>
                <div className="rw-detail-actions">
                  <Badge cls={statusBadgeCls(status)} soft>{model.status_label}</Badge>
                  <HumanReferenceBadge humanReference={humanRef} />
                  <button className="btn ghost small" onClick={refreshEverything}
                          disabled={loadingDetail || action === "refresh"}>
                    {action === "refresh" ? "Refreshing…" : "Refresh"}
                  </button>
                  {RUNNING.has(status) && (
                    <button className="btn ghost small" onClick={() => handleStop(selectedId)}
                            disabled={action === "stop"}>
                      {action === "stop" ? "Stopping…" : "Stop run"}
                    </button>
                  )}
                  <button className="btn ghost small danger" onClick={() => handleDelete(selectedId)}
                          disabled={action === "delete"}>
                    {action === "delete" ? "Deleting…" : "Delete"}
                  </button>
                </div>
              </div>

              {detail._fallback && (
                <p className="muted small rw-fallback-note">
                  Run is starting; status file not available yet. Showing the latest known summary —
                  this will update automatically.
                </p>
              )}

              {/* A finished run has nothing left to track; the checklist would
                  only restate that every stage is done. */}
              {status !== "results_ready" && (
                <RunStatusStepper status={detail.status} files={files} model={model} />
              )}

              <NextAction
                model={model}
                commands={commands}
                busy={busy}
                runId={selectedId}
                onExplore={() => onExploreDataset(`run:${selectedId}`)}
                onRunPost={runPostInterpro}
                onRetry={() => handleRetry(selectedId)}
                onRetryLocal={() => handleRetryLocal(selectedId)}
                onDelete={() => handleDelete(selectedId)}
                onRefresh={() => loadDetail(selectedId)}
                refreshing={loadingDetail}
              />
            </>
          ) : (
            <div className="rw-detail-recover">
              <Empty
                title="Run is starting; status file not available yet."
                hint="The run folder exists but its status could not be read this moment. This usually resolves on its own." />
              <div className="rw-detail-recover-actions">
                <button className="btn ghost small" onClick={() => loadDetail(selectedId)}
                        disabled={loadingDetail}>
                  {loadingDetail ? "Retrying…" : "Retry"}
                </button>
                <code className="muted small">runs/{selectedId}/status.json</code>
              </div>
            </div>
          )}
        </div>
      )}
    </section>
  );
}

// What went wrong and what to do about it. The full log is a download rather
// than a panel: a terminal dump is the wrong first thing to show someone whose
// run failed, and it is still there for whoever needs to investigate.
function RunFailurePanel({ model, runId, busy, onRetry, onRetryLocal, onDelete }) {
  const stage = model.failed_stage || "";
  const reason = model.last_error || "";
  // A local data-acquisition failure is repaired in this run, not by asking for a
  // second one. The older retry, which scaffolds a fresh run, stays for the cases it
  // was written for.
  const canRetryLocal = model.can_retry_local_preparation && Boolean(onRetryLocal);
  const canRetry = !canRetryLocal && model.core_only && Boolean(onRetry);

  const STAGE_LABEL = {
    pre_cluster_data_acquisition: "Pre-cluster data acquisition failed",
    core_model_collection: "Gene and protein model collection failed",
    gene_locus_resolution: "Gene locus resolution failed",
    core_primary_fasta: "Primary protein FASTA could not be generated",
  };
  const who = model.failed_species
    ? model.failed_species.replace(/_/g, " ").replace(/^./, (c) => c.toUpperCase())
    : "";
  const what = STAGE_LABEL[stage] || (stage ? stage.replace(/_/g, " ") : "");
  const headline = what
    ? (who ? `${who}: ${what.charAt(0).toLowerCase()}${what.slice(1)}` : what)
    : "The run stopped before producing results.";

  return (
    <div className="next-action warn">
      <h3>{headline}</h3>
      {/* One concise cause from the recorded backend status, never a traceback. */}
      {reason && <p className="muted small rw-fail-reason">{reason.split("\n")[0]}</p>}
      <div className="rw-fail-actions">
        {canRetryLocal && (
          <button className="btn primary small" onClick={onRetryLocal} disabled={busy}>
            {busy ? "Restarting…" : "Retry local preparation"}
          </button>
        )}
        {canRetry && (
          <button className="btn primary small" onClick={onRetry} disabled={busy}>
            {busy ? "Starting retry…" : "Retry with the same input"}
          </button>
        )}
        <a className="btn ghost small" href={api.runDiagnosticsUrl(runId)} download>
          Download diagnostics
        </a>
        {onDelete && (
          <button className="btn ghost small" onClick={onDelete} disabled={busy}>
            Delete run
          </button>
        )}
      </div>
    </div>
  );
}

/**
 * Which post-cluster analyses this run really has, named from the backend's own verdict.
 *
 * This sentence used to assert that domain architecture and boundary consistency were
 * available whenever the card was shown, which is how a run awaiting its cluster round-trip
 * came to advertise two analyses it had not performed.
 */
function analysesSentence(model) {
  const views = model?.view_availability || {};
  const named = [["domain_architecture", "Domain architecture"],
                 ["boundary_consistency", "boundary consistency"]];
  const have = named.filter(([k]) => views[k]?.available).map(([, label]) => label);
  if (have.length === named.length) return "Domain architecture and boundary consistency are available.";
  if (!have.length) return "Every required view has current data.";
  return `${have.join(" and ")} ${have.length > 1 ? "are" : "is"} available.`;
}

// The single most relevant next step, in calm language (PART 3).
function NextAction({ model, commands, busy, runId, onExplore, onRunPost, onRetry,
                     onRetryLocal, onDelete, onRefresh, refreshing }) {
  const status = model.status;
  const cmd = commands?.cluster_roundtrip?.portable_command
    || model.cluster_command
    || commands?.cluster_roundtrip?.command;

  if (status === "running") {
    return (
      <div className="next-action running">
        <Spinner label={model.current_step
          || "Running the local pre-InterPro pipeline…"} />
      </div>
    );
  }

  if (status === "results_ready") {
    return (
      <div className="next-action ok">
        <div>
          <h3>Results ready</h3>
          <p className="muted small">
            {model.primary_fasta_count} primary protein(s) analysed. {analysesSentence(model)}
          </p>
        </div>
        <button className="btn primary" onClick={onExplore}>Open results →</button>
      </div>
    );
  }

  if (status === "pre_interpro_running") {
    return (
      <div className="next-action running">
        <Spinner label="Preparing the primary protein FASTA — this can take a while." />
      </div>
    );
  }

  if (status === "post_interpro_running") {
    return (
      <div className="next-action running">
        <Spinner label="Post-InterPro analysis running…" />
      </div>
    );
  }

  if (status === "cluster_fetch_complete") {
    return (
      <div className="next-action">
        <h3>Cluster results fetched — run Post-InterPro analysis</h3>
        <p className="muted small">
          The one-command round-trip normally does this automatically after fetch. If it did not,
          you can run it locally now.
        </p>
        <button className="btn primary" onClick={onRunPost} disabled={busy}>
          {busy ? "Starting…" : "Run Post-InterPro locally"}
        </button>
      </div>
    );
  }

  if (status === "stopped") {
    return (
      <div className="next-action warn">
        <h3>Stopped by user</h3>
        <p className="muted small">
          This run was stopped. Partial outputs are kept. You can restart the local pre-InterPro
          step from the commands below, or delete this run from the list on the left.
        </p>
      </div>
    );
  }

  if (FAILED.has(status)) {
    return <RunFailurePanel model={model} runId={runId} busy={busy} onRetry={onRetry}
                            onRetryLocal={onRetryLocal} onDelete={onDelete} />;
  }

  // A freshly created run needs the local pre-InterPro step first, not the cluster.
  if (status === "created" && !model.cluster_required) {
    const preCmd = `python scripts/run_pre_interpro_for_run.py --run-id ${model.run_id}`;
    return (
      <div className="next-action">
        <h3>Prepare FASTA (local pre-InterPro)</h3>
        <p className="muted small">
          Run the local pre-InterPro pipeline to produce this run's primary FASTA. No cluster
          login is needed yet.
        </p>
        <div className="cmd-line"><code>{preCmd}</code><CopyButton text={preCmd} /></div>
      </div>
    );
  }

  // cluster_required / cluster_running → exactly ONE recommended terminal command.
  const running = status === "cluster_running";
  const gene = model.gene_symbol || "gene";
  return (
    <div className="next-action">
      <h3>{running ? "Cluster annotation running" : "Cluster input ready"}</h3>
      <p className="muted small">
        {running
          ? "The cluster round-trip is running in your terminal."
          : `The local ${gene} preparation is complete. Run the cluster roundtrip to `
            + "generate InterProScan, pyTMHMM, Domain Architecture and Boundary results."}
      </p>
      {cmd && (
        <div className="cmd-line">
          <code>{cmd}</code>
          <CopyButton text={cmd} />
        </div>
      )}
      <p className="muted small">
        LRZ password / MFA happens in your terminal; no credentials are stored in the webapp.
        The script reuses the SSH session to avoid repeated prompts, and Post-InterPro runs
        automatically once the annotation is back.
      </p>
      {onRefresh && (
        <button className="btn ghost small" onClick={onRefresh} disabled={refreshing}>
          {refreshing ? "Refreshing…" : "Refresh status"}
        </button>
      )}
    </div>
  );
}
