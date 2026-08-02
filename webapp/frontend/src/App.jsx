import { useCallback, useEffect, useRef, useState } from "react";
import "./App.css";
import {
  api, forDataset, payloadMatchesDataset, setActiveDataset as setApiDataset,
} from "./api";
import { getDatasetLabels } from "./labels";
// A run that leaves a stable state (e.g. a retry) starts being polled again, because the
// polling effect depends on the status it last observed.
import { POLL_INTERVAL_MS, TERMINAL_RUN_STATES } from "./runStates";
import { isResolvedState } from "./pages/viewers/common";
import { Badge } from "./ui";
import StartPage from "./pages/StartPage";
import Overview from "./pages/Overview";
import GeneExplorer from "./pages/GeneExplorer";
import FigureGallery from "./pages/FigureGallery";
import FreezeViewer from "./pages/FreezeViewer";
import BoundaryPage from "./pages/BoundaryPage";
import RunWorkflowPage from "./pages/RunWorkflowPage";
import DatasetSwitcher from "./pages/DatasetSwitcher";
import { ScientificSelectionProvider } from "./components/ScientificSelectionContext";

function FigureGalleryWithSelection({ model, openBoundary }) {
  const models = model?.coordinate_models || model?.models
    || model?.protein_coordinate_model?.models || [];
  const first = models[0] || { species_id: "", species: "" };
  const species = {
    species: first.species_id || first.species || "",
    species_id: first.species_id || first.species || "",
    scientific_name: first.scientific_name,
    selected_primary_protein: first.protein_id,
  };
  return (
    <ScientificSelectionProvider species={species} model={model}>
      <FigureGallery model={model} openBoundary={openBoundary} />
    </ScientificSelectionProvider>
  );
}

const ACTIVE_DS_KEY = "edc.activeDataset";
const DEV = Boolean(import.meta.env?.DEV);

/**
 * Whether the model's availability manifest has settled this analysis.
 *
 * Used to keep a "pending" badge off a navigation entry whose analysis is finished. A
 * missing manifest returns false, so the previous, more cautious labelling stands.
 */
function isResolvedAnalysis(model, analysisName) {
  const analyses = model?.analysis_availability?.analyses || [];
  const found = analyses.find((a) => a.analysis_name === analysisName);
  return Boolean(found) && isResolvedState(found.status);
}

// Primary navigation — deliberately small (PART 1 / PART 11). Each explore page
// is gated by the active dataset's available_views; unavailable pages show a
// calm "pending" note instead of an error.
const EXPLORE_PAGES = [
  ["overview", "Overview", "overview"],
  ["gene", "Gene Explorer", "gene_explorer"],
  ["boundary", "Boundary Consistency", "boundary_consistency"],
  ["figures", "Figure Gallery", "figure_gallery"],
];

const VIEW_PENDING_MSG = {
  gene_explorer: "Gene Explorer becomes available once pre-InterPro has prepared the FASTA.",
  domain_architecture: "Domain annotation pending. Run cluster annotation to enable this view.",
  boundary_consistency: "Post-InterPro domain annotation pending. Run the cluster round-trip, then Post-InterPro.",
  exon_domain_boundaries: "Exon–domain boundary analysis becomes available after domain annotation and the exon map exist.",
  figure_gallery: "Figures appear after the analysis produces them.",
  overview: "Overview becomes available once this run has been prepared.",
};

export default function App() {
  const [page, setPage] = useState("start");
  const [health, setHealth] = useState(null);
  const [datasets, setDatasets] = useState([]);
  const [activeId, setActiveId] = useState("example");
  const [datasetInfo, setDatasetInfo] = useState(null); // status + available_views + human_reference
  const [datasetModel, setDatasetModel] = useState(null);
  const [loading, setLoading] = useState(false);
  const [geneTarget, setGeneTarget] = useState(null);
  // Bumped by a full refresh to force every dataset-scoped page to remount and
  // refetch (so e.g. Boundary Consistency appears after a run finishes, without
  // a browser reload).
  const [refreshNonce, setRefreshNonce] = useState(0);
  // Which load is current, and how to cancel the one before it. Refs, not state:
  // these must be readable by an in-flight promise without re-rendering.
  const epochRef = useRef(0);
  const abortRef = useRef(null);
  const pendingRef = useRef("example");

  const openGene = useCallback((target) => {
    setGeneTarget(target ? { ...target, _t: Date.now() } : null);
    setPage("gene");
  }, []);

  const refreshHealth = useCallback(async () => {
    try {
      const h = await api.health();
      setHealth(h);
      return h;
    } catch {
      setHealth({ status: "offline" });
      return null;
    }
  }, []);

  // --- load epochs ---------------------------------------------------------
  // One counter and one abort controller decide which reply is allowed to reach the
  // screen. `epochRef` holds the newest load; anything older is discarded even if it
  // resolves later, which is what stops a dataset the user has already left from
  // overwriting the one they are looking at.
  const beginLoad = useCallback((target) => {
    const epoch = epochRef.current + 1;
    epochRef.current = epoch;
    abortRef.current?.abort();
    abortRef.current = new AbortController();
    pendingRef.current = target;
    setLoading(true);
    return epoch;
  }, []);

  const endLoad = useCallback((epoch) => {
    if (epoch === epochRef.current) setLoading(false);
  }, []);

  /** Fetch a dataset's status and canonical model, applying them only if still current. */
  const loadInto = useCallback(async (target, epoch) => {
    const signal = abortRef.current?.signal;
    const client = forDataset(target, signal);
    const [info, model] = await Promise.all([
      client.datasetStatus(target).catch(() => null),
      client.datasetModel().catch(() => null),
    ]);
    // Two independent checks: the reply must belong to this load, and it must name the
    // dataset we asked about. Either alone would let a mismatched payload through.
    if (epoch !== epochRef.current) return null;
    const okInfo = payloadMatchesDataset(info, target);
    const okModel = payloadMatchesDataset(model, target);
    if (okInfo) setDatasetInfo(info);
    if (okModel) setDatasetModel(model);
    return okInfo ? info : null;
  }, []);

  const loadDatasets = useCallback(async () => {
    try {
      const res = await api.datasets();
      const list = res?.datasets || [];
      setDatasets(list);
      return list;
    } catch {
      setDatasets([]);
      return [];
    }
  }, []);

  // Point the API client + all pages at a dataset, then load its data.
  // This is the SINGLE source of truth for the active dataset: it updates the
  // API client, activeId, datasetInfo and persistence together so the dropdown
  // label can never disagree with the displayed content. Both the Dataset
  // dropdown and My Runs → Explore call this exact function.
  //
  // The previous dataset's status and model are cleared *before* the new request
  // starts. Keeping them while the new dataset loaded is what let one dataset's
  // transcript and protein values appear briefly under another dataset's name.
  // Each switch also takes the next epoch and aborts the one before it, so a slow
  // reply that lost the race is dropped instead of overwriting the newer data.
  const selectDataset = useCallback(async (id, goto) => {
    const target = id || "example";
    const epoch = beginLoad(target);
    setApiDataset(target);
    setActiveId(target);
    setDatasetInfo(null);
    setDatasetModel(null);
    setGeneTarget(null);
    try { localStorage.setItem(ACTIVE_DS_KEY, target); } catch { /* ignore */ }
    if (DEV) console.debug("[activeDataset] →", target, goto ? `(goto: ${goto})` : "");
    // Destination first: the page frame may switch immediately, it is only the
    // scientific content that has to wait for the identity check.
    if (goto) setPage(goto);
    try {
      // Make sure the dataset list contains the target so the dropdown option
      // exists (e.g. a run created during this session, opened via My Runs).
      if (target !== "example" && !datasets.some((d) => d.id === target)) {
        await loadDatasets();
      }
      const info = await loadInto(target, epoch);
      if (DEV && info) {
        console.debug("[activeDataset] loaded", {
          dataset_id: info.id, gene_symbol: info.gene_symbol,
          pipeline_type: info.pipeline_type, support_level: info.support_level,
          status: info.status,
        });
      }
      return info;
    } finally {
      endLoad(epoch);
    }
  }, [datasets, loadDatasets, beginLoad, endLoad, loadInto]);

  // Full refresh: reload the dataset list + the active dataset's status/summary,
  // then bump the remount nonce so every explore page refetches its own data
  // (Overview, Gene Explorer, Boundary Consistency, Figure Gallery, …).
  const refreshAll = useCallback(async (id) => {
    const target = id || activeId;
    const epoch = beginLoad(target);
    try {
      await loadDatasets();
      const info = await loadInto(target, epoch);
      if (epoch === epochRef.current) setRefreshNonce((n) => n + 1);
      return info;
    } finally {
      endLoad(epoch);
    }
  }, [activeId, loadDatasets, beginLoad, endLoad, loadInto]);

  useEffect(() => {
    (async () => {
      await refreshHealth();
      const list = await loadDatasets();
      // Restore the previously selected dataset if it still exists, otherwise
      // fall back consistently to the example dataset (no mixed state).
      let restore = "example";
      try {
        const saved = localStorage.getItem(ACTIVE_DS_KEY);
        if (saved && (saved === "example" || (list || []).some((d) => d.id === saved))) {
          restore = saved;
        }
      } catch { /* ignore */ }
      await selectDataset(restore);
    })();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Stage-aware runs refresh themselves after the external cluster roundtrip.
  // The polling endpoint is read-only; scientific outputs are rebuilt by the
  // roundtrip before it publishes results_ready/results_partial.
  //
  // Polling stops at a stable state and resumes if the run leaves one. On every
  // meaningful change the dataset list, the canonical model and the remount nonce are
  // refreshed together, so My Runs, the selector, Summary, availability, the Figure
  // Gallery and Data & Downloads all move at once — a finished roundtrip shows up
  // without the browser being reloaded.
  useEffect(() => {
    if (!activeId.startsWith("run:")) return undefined;
    if (TERMINAL_RUN_STATES.has(datasetInfo?.status)) return undefined;
    const timer = window.setInterval(async () => {
      const client = forDataset(activeId);
      const info = await client.datasetStatus(activeId).catch(() => null);
      // The interval outlives a dataset switch by up to one tick, so a reply is only
      // used while its dataset is still the selected one.
      if (!info || activeId !== pendingRef.current
          || !payloadMatchesDataset(info, activeId)) return;
      const changed = info.status !== datasetInfo?.status
        || info.current_step !== datasetInfo?.current_step
        || info.index_version !== datasetInfo?.index_version;
      setDatasetInfo(info);
      if (changed || TERMINAL_RUN_STATES.has(info.status)) {
        const model = await client.datasetModel().catch(() => null);
        if (activeId !== pendingRef.current) return;
        if (model && payloadMatchesDataset(model, activeId)) setDatasetModel(model);
        await loadDatasets();
        setRefreshNonce((n) => n + 1);
      }
    }, POLL_INTERVAL_MS);
    return () => window.clearInterval(timer);
  }, [activeId, datasetInfo?.status, datasetInfo?.current_step,
      datasetInfo?.index_version, loadDatasets]);

  const views = datasetModel?.available_views || datasetInfo?.available_views || {};
  const isExample = activeId === "example";
  const eventLayerType = datasetModel?.event_layer?.type || "none";
  const hasEvent = eventLayerType === "validated";
  // Config-driven display labels (FGFR2 wording is the default fallback, so the
  // UI is unchanged for FGFR2). Maps generic view ids -> visible nav labels.
  const labels = getDatasetLabels(datasetInfo);
  const navLabelFor = (id, fallback) =>
    id === "gene" ? labels.geneExplorer
      : id === "boundary" ? labels.boundaryRelation
        : fallback;

  // The boundary nav item is event-specific when an event region is configured
  // (FGFR2 -> "Boundary Consistency"), otherwise it becomes the generic
  // all-exon "Exon–Domain Boundaries" view. This is how the UI stops assuming a
  // cassette exists for every gene.
  const boundaryView = hasEvent ? "boundary_consistency" : "exon_domain_boundaries";
  // Core-only navigation (PART 9): Home, Overview, Gene Explorer,
  // Exon–Domain Boundaries, My Runs. Figure Gallery is hidden (not a disabled
  // tab), and the Exon–Domain Boundaries tab is always openable — it shows a
  // calm "pending cluster" page before annotation instead of being disabled.
  const navItems = EXPLORE_PAGES
    .map(([id, label, view]) => {
      if (id === "boundary") {
        // Core-only pre-cluster: keep the tab openable but flag it as pending so
        // it never reads as an existing analysis result (PART 6).
        //
        // "Pending" means still waiting. A boundary analysis that resolved as not
        // applicable — a protein encoded by one coding exon has no internal boundaries —
        // is finished, so the badge is dropped and the page explains the state instead.
        const boundaryPending = !hasEvent && !views[boundaryView]
          && !isResolvedAnalysis(datasetModel, "boundary_analysis");
        const boundaryLabel = hasEvent ? labels.boundaryRelation : (
          boundaryPending
            ? <span className="nav-with-badge">Exon–Domain Boundaries<span className="nav-badge">pending</span></span>
            : "Exon–Domain Boundaries");
        return [id, boundaryLabel, boundaryView, !hasEvent ? true : undefined];
      }
      // For a core-only dataset every generic nav item is openable (Overview,
      // Gene Explorer, Exon–Domain Boundaries, Figure Gallery). None are disabled
      // "Not available" tabs — each renders a first-class or calm-pending page.
      const forceEnabled = Boolean(datasetModel) && (id === "overview" || id === "gene" || id === "figures");
      return [id, navLabelFor(id, label), view, forceEnabled ? true : undefined];
    });

  return (
    <div className="shell">
      <header className="topbar">
        <div className="brand" onClick={() => setPage("start")}>
          <span className="brand-mark">EDC</span>
          <div className="brand-text">
            <b>ExonDomainCompare</b>
            <small>Annotation-aware comparative exon–protein analysis</small>
          </div>
        </div>
        <nav className="nav">
          <button className={page === "start" ? "active" : ""} onClick={() => setPage("start")}>Home</button>
          {navItems.map(([id, label, view, forceEnabled]) => {
            const enabled = forceEnabled || Boolean(views[view]);
            return (
              <button
                key={id}
                className={page === id ? "active" : ""}
                onClick={() => setPage(id)}
                disabled={!enabled}
                title={!enabled ? (VIEW_PENDING_MSG[view] || "Not available yet") : ""}
              >
                {label}
              </button>
            );
          })}
          <button className={page === "runs" ? "active" : ""} onClick={() => setPage("runs")}>My Runs</button>
        </nav>
        <div className="run-context">
          {datasets.length > 0 ? (
            <DatasetSwitcher
              datasets={datasets}
              activeId={activeId}
              onChange={(id) => selectDataset(id, page === "start" ? "overview" : undefined)}
            />
          ) : (
            <Badge cls={health?.status === "ok" ? "neutral" : "excluded"} soft>
              {health?.status === "ok" ? "Loading datasets…" : "API offline"}
            </Badge>
          )}
        </div>
      </header>

      <main className="content">
        {page === "start" && (
          <StartPage
            health={health}
            datasets={datasets}
            loading={loading}
            onExploreExample={() => selectDataset("example", "overview")}
            onCreateRun={() => setPage("runs")}
            onOpenRuns={() => setPage("runs")}
            onExploreDataset={(id) => selectDataset(id, "overview")}
          />
        )}
        {page === "runs" && (
          <RunWorkflowPage
            activeId={activeId}
            onExploreDataset={(id) => selectDataset(id, "overview")}
            onDatasetsChanged={loadDatasets}
            onRefreshAll={refreshAll}
            onActiveDeleted={() => selectDataset("example")}
          />
        )}
        {page === "overview" && (
          <Overview key={`${activeId}:${refreshNonce}`} model={datasetModel} datasetInfo={datasetInfo} isExample={isExample}
                    setPage={setPage} openGene={openGene} />
        )}
        {page === "gene" && (
          <GeneExplorer key={`${activeId}:${refreshNonce}`} model={datasetModel} target={geneTarget} labels={labels} setPage={setPage} />
        )}
        {page === "boundary" && (
          <BoundaryPage key={`${activeId}:${refreshNonce}`} model={datasetModel} setPage={setPage}
            openGene={openGene} labels={labels} />
        )}
        {page === "figures" && (
          <FigureGalleryWithSelection
            key={`${activeId}:${refreshNonce}`}
            model={datasetModel}
            openBoundary={() => setPage("boundary")}
          />
        )}
        {page === "freeze" && <FreezeViewer key={`${activeId}:${refreshNonce}`} />}
      </main>
    </div>
  );
}
