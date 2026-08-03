export const API_BASE = import.meta.env?.VITE_API_BASE || "http://localhost:8000";

// --- global active dataset -------------------------------------------------
// The whole app explores exactly one dataset at a time:
//   "example"        -> validated 30-species freeze (read-only)
//   "run:<run_id>"   -> a local custom run under runs/<run_id>/
// Every explorable data endpoint carries this as ?dataset=…, so all pages
// (Overview, Gene Explorer, Domain Architecture, Boundary Consistency, …) show
// the selected dataset without a backend-global switch.
let activeDataset = "example";
const listeners = new Set();

export function getActiveDataset() {
  return activeDataset;
}
export function setActiveDataset(ds) {
  activeDataset = ds || "example";
  listeners.forEach((fn) => {
    try { fn(activeDataset); } catch { /* ignore */ }
  });
}
export function onActiveDatasetChange(fn) {
  listeners.add(fn);
  return () => listeners.delete(fn);
}

function withDataset(path, dataset) {
  const sep = path.includes("?") ? "&" : "?";
  const ds = dataset || activeDataset;
  return `${path}${sep}dataset=${encodeURIComponent(ds)}`;
}

async function jget(path, opts) {
  const res = await fetch(`${API_BASE}${path}`, { signal: opts?.signal });
  if (!res.ok) {
    let detail = res.statusText;
    try { detail = (await res.json()).detail || detail; } catch { /* ignore */ }
    const err = new Error(detail);
    err.status = res.status;
    throw err;
  }
  return res.json();
}

async function jpost(path, body, opts) {
  const res = await fetch(`${API_BASE}${path}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: body ? JSON.stringify(body) : undefined,
    signal: opts?.signal,
  });
  const data = await res.json().catch(() => ({}));
  if (!res.ok) {
    const err = new Error(data?.detail || res.statusText);
    err.status = res.status;
    throw err;
  }
  return data;
}

async function jdelete(path) {
  const res = await fetch(`${API_BASE}${path}`, { method: "DELETE" });
  const data = await res.json().catch(() => ({}));
  if (!res.ok) {
    const err = new Error(data?.detail || res.statusText);
    err.status = res.status;
    throw err;
  }
  return data;
}

// --- per-call dataset scope ------------------------------------------------
// `activeDataset` alone cannot express *which* dataset a request was made for: it is a
// mutable global read at request time, so a response arriving after the user switched
// datasets was indistinguishable from a fresh one and got rendered under the new heading.
// A scope pins the dataset (and an abort signal) for the duration of one call. It is safe
// to keep in a module variable because every api method builds its URL synchronously,
// before the first await, so two overlapping calls can never observe each other's scope.
let callScope = null;

function scoped(dataset, signal, fn) {
  const previous = callScope;
  callScope = { dataset, signal };
  try {
    return fn();
  } finally {
    callScope = previous;
  }
}

// Dataset-scoped GET/POST: appends the explicit ?dataset= of the current scope, falling
// back to the active dataset for unscoped legacy callers.
const dget = (path) => jget(withDataset(path, callScope?.dataset), callScope);
const dpost = (path, body) => jpost(withDataset(path, callScope?.dataset), body, callScope);

export const api = {
  health: () => jget("/api/health"),
  presets: () => jget("/api/presets"),
  parseSpecies: (text) => jpost("/api/species/parse", { text }),

  // --- datasets (switcher) -------------------------------------------------
  datasets: () => jget("/api/datasets"),
  datasetStatus: (id) => jget(`/api/datasets/${encodeURIComponent(id)}/status`),
  humanReference: () => jget("/api/datasets/human-reference"),

  // --- local run workflow (runs/<run_id>/) — safe: no SSH/SLURM ------------
  localRuns: () => jget("/api/local-runs"),
  localRun: (runId) => jget(`/api/local-runs/${encodeURIComponent(runId)}`),
  localRunCreate: (payload) => jpost("/api/local-runs/create", payload),
  localRunStart: (payload) => jpost("/api/local-runs/start", payload),
  analysisRouter: (gene, mode = "auto") =>
    jget(`/api/analysis-router?gene=${encodeURIComponent(gene || "")}&mode=${encodeURIComponent(mode)}`),
  resolveRunInputs: (payload) => jpost("/api/runs/resolve-inputs", payload),
  localRunCommands: (runId) => jget(`/api/local-runs/${encodeURIComponent(runId)}/commands`),
  localRunStatus: (runId) => jget(`/api/local-runs/${encodeURIComponent(runId)}/status`),
  localRunRefresh: (runId) => jpost(`/api/local-runs/${encodeURIComponent(runId)}/refresh`),
  localRunRefreshAll: (runId) => jpost(`/api/local-runs/${encodeURIComponent(runId)}/refresh-all`),
  stopLocalRun: (runId) => jpost(`/api/local-runs/${encodeURIComponent(runId)}/stop`),
  deleteLocalRun: (runId) => jdelete(`/api/local-runs/${encodeURIComponent(runId)}`),
  legacyWebRuns: () => jget("/api/web-runs/legacy"),
  startPreInterpro: (runId, body) =>
    jpost(`/api/local-runs/${encodeURIComponent(runId)}/start-preinterpro`, body || {}),
  stopPreInterpro: (runId) =>
    jpost(`/api/local-runs/${encodeURIComponent(runId)}/stop-preinterpro`),
  startPostInterpro: (runId) =>
    jpost(`/api/local-runs/${encodeURIComponent(runId)}/start-post-interpro`),
  preInterproLogs: (runId, tail = 200) =>
    jget(`/api/local-runs/${encodeURIComponent(runId)}/logs/preinterpro?tail=${tail}`),
  coreRunLogs: (runId, tail = 400) =>
    jget(`/api/local-runs/${encodeURIComponent(runId)}/logs/core?tail=${tail}`),
  retryPreCluster: (runId) =>
    jpost(`/api/local-runs/${encodeURIComponent(runId)}/retry-precluster`),
  // Rebuilds the failed local stage of the SAME run, so one request stays one run.
  retryLocalPreparation: (runId) =>
    jpost(`/api/local-runs/${encodeURIComponent(runId)}/retry-local-preparation`),
  // Logs stay available for diagnostics without living inside a run card.
  runDiagnosticsUrl: (runId) =>
    `${API_BASE}/api/local-runs/${encodeURIComponent(runId)}/diagnostics`,

  loadExample: () => jpost("/api/runs/example/load"),
  openRun: (path) => jpost("/api/runs/open", { path }),
  runStatus: (runId) => jget(`/api/runs/${runId}/status`),
  runLogs: (runId, tail = 600) => jget(`/api/runs/${runId}/logs?tail=${tail}`),

  // --- dataset-scoped explorable data (auto ?dataset=) ---------------------
  // Canonical application model. Overview, Gene Explorer, Figure Gallery and
  // Boundary all consume this single endpoint for every gene/event layer.
  datasetModel: () => dget("/api/runs/current/dataset-model"),
  summary: (rebuild = false) => dget(`/api/runs/current/summary${rebuild ? "?rebuild=true" : ""}`),
  evidenceStack: () => dget("/api/runs/current/evidence-stack"),
  species: () => dget("/api/runs/current/species"),
  speciesOne: (sp) => dget(`/api/runs/current/species/${encodeURIComponent(sp)}`),
  figures: () => dget("/api/runs/current/figures"),
  downloads: () => dget("/api/runs/current/downloads"),
  // The canonical availability contract for Data & Downloads; the page renders
  // exactly what this returns and never infers availability from a filename.
  packageCapabilities: (scope) => dget(
    `/api/runs/current/package-capabilities${scope ? `?scope=${encodeURIComponent(scope)}` : ""}`),
  packageCatalogue: () => dget("/api/runs/current/package-catalogue"),
  createPackage: (selection) => dpost("/api/runs/current/packages", selection),
  packageStatus: (jobId) => dget(`/api/runs/current/packages/${encodeURIComponent(jobId)}`),
  freeze: () => dget("/api/runs/current/freeze"),
  cassette: () => dget("/api/runs/current/cassette"),
  coordinates: () => dget("/api/runs/current/coordinates"),
  msa: () => dget("/api/runs/current/msa"),
  synteny: () => dget("/api/runs/current/synteny"),
  story: () => dget("/api/runs/current/story"),
  domainArchitecture: () => dget("/api/runs/current/domain-architecture"),
  domainArchitectureSummary: () => dget("/api/runs/current/domain-architecture/summary"),
  domainArchitectureSpecies: () => dget("/api/runs/current/domain-architecture/species"),
  domainArchitectureQc: () => dget("/api/runs/current/domain-architecture/qc"),
  boundaryConsistency: () => dget("/api/runs/current/boundary-consistency"),
  boundaryConsistencySummary: () => dget("/api/runs/current/boundary-consistency/summary"),
  boundaryConsistencyMatrix: () => dget("/api/runs/current/boundary-consistency/matrix"),
  boundaryConsistencyOutliers: () => dget("/api/runs/current/boundary-consistency/outliers"),

  // --- shared-pipeline indices (one contract for every gene, PART 7) -------
  // These read runs/<id>/website_indices/ (rich canonical indices) so the same
  // shared renderer works for FGFR1 and any future shared_gene_pipeline dataset.
  sharedOverview: () => dget("/api/runs/current/shared/overview"),
  sharedEvidenceStack: () => dget("/api/runs/current/shared/evidence-stack"),
  sharedGeneExplorer: () => dget("/api/runs/current/shared/gene-explorer"),
  sharedProteinArchitecture: () => dget("/api/runs/current/shared/protein-architecture"),
  sharedSynteny: () => dget("/api/runs/current/shared/synteny"),
  sharedEventEvidence: () => dget("/api/runs/current/shared/event-evidence"),
  sharedDomainArchitecture: () => dget("/api/runs/current/shared/domain-architecture"),
  sharedExonDomainBoundaries: () => dget("/api/runs/current/shared/exon-domain-boundaries"),
  sharedFigures: () => dget("/api/runs/current/shared/figures"),
  sharedAvailableViews: () => dget("/api/runs/current/shared/available-views"),

  // --- LEGACY core-only index endpoints (superseded by shared/* above) -----
  // Kept only for backward compatibility / tooling. The shared renderer no
  // longer uses these index methods; new shared_gene_pipeline runs go through
  // the shared/* contract. coreCapability is still used for cluster status.
  coreSummary: () => dget("/api/runs/current/core/summary"),
  coreEvidenceStack: () => dget("/api/runs/current/core/evidence-stack"),
  corePrimarySelection: () => dget("/api/runs/current/core/primary-selection"),
  coreGeneAnalysis: () => dget("/api/runs/current/core/gene-analysis"),
  coreExonProteinMap: () => dget("/api/runs/current/core/exon-protein-map"),
  coreProteinArchitecture: () => dget("/api/runs/current/core/protein-architecture"),
  coreEventEvidenceIndex: () => dget("/api/runs/current/core/event-evidence-index"),
  coreFigures: () => dget("/api/runs/current/core/figures"),
  coreDomainArchitecture: () => dget("/api/runs/current/core/domain-architecture"),
  coreSynteny: () => dget("/api/runs/current/core/synteny"),
  coreExonDomainBoundaries: () => dget("/api/runs/current/core/exon-domain-boundaries"),
  coreEventCandidates: () => dget("/api/runs/current/core/event-candidates"),
  coreEventEvidence: () => dget("/api/runs/current/core/event-evidence"),
  coreCapability: () => dget("/api/runs/current/core/capability"),
  coreValidate: (runId) => jget(`/api/local-runs/${encodeURIComponent(runId)}/core-validate`),
};

// --- explicitly bound client ----------------------------------------------
/**
 * The same api surface, bound to one dataset and one abort signal.
 *
 * Callers use this instead of `api` so a request can never silently pick up whichever
 * dataset happens to be selected by the time it runs. Every method keeps its own name and
 * arguments; only the dataset and the cancellation are supplied here.
 */
export function forDataset(datasetId, signal) {
  return new Proxy(api, {
    get(target, prop) {
      const value = target[prop];
      if (typeof value !== "function") return value;
      return (...args) => scoped(datasetId || "example", signal,
        () => value.apply(target, args));
    },
  });
}

/** The dataset a payload belongs to, or "" when the payload carries no identity. */
export function payloadDataset(payload) {
  if (!payload || typeof payload !== "object") return "";
  return payload.dataset_id || payload.dataset?.id || payload.id || "";
}

/** The run a payload belongs to, or "" when the payload carries no identity. */
export function payloadRunId(payload) {
  if (!payload || typeof payload !== "object") return "";
  return payload.run_id || payload.dataset?.run_id || "";
}

/** The run id inside a dataset id, so `run:<id>` and `<id>` compare equal. */
export function runIdOf(datasetId) {
  const id = String(datasetId || "");
  return id.startsWith("run:") ? id.slice(4) : "";
}

/**
 * The run a dataset id refers to. The validated freeze is its own run, named `example`.
 */
function expectedRunId(datasetId) {
  const id = String(datasetId || "");
  return runIdOf(id) || id;
}

/**
 * Whether a payload may be applied while `datasetId` is selected.
 *
 * A payload without identity is accepted: several small endpoints return bare arrays, and
 * for those the abort signal and the caller's own epoch check already do the work. A
 * payload that *does* name a dataset or a run must name the selected one — that is the
 * check whose absence let one gene's protein appear under another gene's heading.
 */
export function payloadMatchesDataset(payload, datasetId) {
  const wantDataset = String(datasetId || "");
  const wantRun = expectedRunId(wantDataset);
  const gotDataset = String(payloadDataset(payload) || "");
  const gotRun = String(payloadRunId(payload) || "");
  if (gotDataset && gotDataset !== wantDataset && gotDataset !== wantRun) return false;
  if (gotRun && wantRun && gotRun !== wantRun) return false;
  return true;
}

/** Whether a scientific model was built from the indices named by its status payload. */
export function payloadMatchesIndexVersion(payload, status) {
  const expected = String(status?.index_version || "");
  if (!expected) return true;
  return String(payload?.index_version || "") === expected;
}

export const fileUrl = (path, inline = false) =>
  `${API_BASE}/api/download?path=${encodeURIComponent(path)}${inline ? "&inline=true" : ""}`;

export const runFileUrl = (runId, path, inline = false) =>
  `${API_BASE}/api/runs/${encodeURIComponent(runId)}/files?path=${encodeURIComponent(path)}`
  + `${inline ? "&inline=true" : ""}`;

// Absolute URL for a server-provided (already API-rooted) url, e.g. figure urls.
export const absUrl = (u) => (!u ? u : (/^https?:\/\//.test(u) ? u : `${API_BASE}${u}`));

// Capped text preview of a small project file (provenance hub).
export const filePreview = (path, maxBytes = 20000) =>
  jget(`/api/file-preview?path=${encodeURIComponent(path)}&max_bytes=${maxBytes}`);
