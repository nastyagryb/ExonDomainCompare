// Derive a stepper state (done | active | failed | pending) from status.json
// fields plus the file-existence checks returned by /refresh.
function tokenState(raw) {
  const v = String(raw || "").toLowerCase();
  if (v === "complete" || v === "available" || v === "finished") return "done";
  if (v === "running") return "active";
  if (v === "failed") return "failed";
  if (v === "submitted") return "active";
  return "pending";
}

// A run is in its FINAL state when the authoritative model / status.json say so. In that
// case every earlier pipeline step is complete by definition — no intermediate cluster step
// (InterProScan / pyTMHMM submitted) may remain blue/current just because a stale per-job
// field was never cleared. We therefore derive step states from the final status, not from
// partially stale job fields.
function isFinalStatus(st, model) {
  const s = (v) => String(v || "").toLowerCase();
  const m = model || {};
  return (
    s(m.status) === "results_ready" ||
    s(st.status) === "results_ready" ||
    s(m.post_interpro_status) === "complete" ||
    s(st.post_interpro_status) === "complete" ||
    s(m.next_action) === "open_results" ||
    s(st.next_action) === "open_results"
  );
}

function buildSteps(status, files, model) {
  const f = files || {};
  const st = status || {};
  const primary = f.primary_fasta || String(st.primary_fasta_status).toLowerCase() === "available";
  const fetched = Boolean(f.interproscan_output && f.pytmhmm_output);
  const indices = f.website_indices || String(st.website_indices_status).toLowerCase() === "complete";
  const final = isFinalStatus(st, model);

  const steps = [
    { key: "created", label: "Run created", state: "done" },
    { key: "pre", label: "Pre-InterPro pipeline", state: tokenState(st.pre_interpro_status) },
    { key: "primary", label: "Primary FASTA ready", state: primary ? "done" : "pending" },
    { key: "ips", label: "InterProScan submitted", state: tokenState(st.interproscan_status) },
    { key: "tm", label: "pyTMHMM submitted", state: tokenState(st.pytmhmm_status) },
    { key: "fetch", label: "Cluster outputs fetched", state: fetched ? "done" : "pending" },
    { key: "post", label: "Post-InterPro analysis", state: tokenState(st.post_interpro_status) },
    { key: "indices", label: "Website indices ready", state: indices ? "done" : "pending" },
    { key: "explore", label: "Explore results", state: f.explorable ? "done" : "pending" },
  ];

  if (final) {
    // Everything up to and including results is complete; nothing stays active/pending.
    return steps.map((s) => ({ ...s, state: "done" }));
  }

  // A failed step is TERMINAL for the pre-cluster phase: the pipeline stopped, so
  // no later step may be shown as "active" (that previously made a failed run look
  // like it was still working, e.g. "Primary FASTA ready" glowing after a failure).
  const failedIdx = steps.findIndex((s) => s.state === "failed");
  if (failedIdx >= 0) {
    return steps.map((s, i) => (i > failedIdx && s.state !== "done"
      ? { ...s, state: "blocked" } : s));
  }

  // mark the first non-done, non-failed step as "active" for a clear focus point
  const hasActive = steps.some((s) => s.state === "active");
  if (!hasActive) {
    const idx = steps.findIndex((s) => s.state === "pending");
    if (idx >= 0) steps[idx] = { ...steps[idx], state: "active" };
  }
  return steps;
}

const GLYPH = { done: "✓", active: "●", failed: "!", blocked: "·", pending: "" };

export default function RunStatusStepper({ status, files, model }) {
  const steps = buildSteps(status, files, model);
  return (
    <ol className="run-stepper">
      {steps.map((s) => (
        <li key={s.key} className={`step step-${s.state}`}>
          <span className="step-dot">{GLYPH[s.state]}</span>
          <span className="step-label">{s.label}</span>
        </li>
      ))}
    </ol>
  );
}
