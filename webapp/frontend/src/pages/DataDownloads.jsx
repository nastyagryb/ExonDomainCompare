import { useEffect, useMemo, useState } from "react";
import { api, fileUrl } from "../api";
import { Badge } from "../ui";
import { useScientificSelection } from "../components/ScientificSelectionContext";

/**
 * Data & Downloads.
 *
 * The page renders exactly one thing: the backend capability object for the
 * selected Scope. Availability, reasons, dependencies, presets and estimated
 * sizes all come from there, so the page can never offer a file the run cannot
 * deliver, and a single-species run never sees a comparative option.
 *
 * My Runs stays responsible for pipeline control; this page only packages
 * scientific tables, sequences, alignments, figures and workbooks.
 */
export default function DataDownloads({ downloads, eventType, availability }) {
  // FGFR2 validated example keeps its immutable Files view.
  if (eventType === "validated") {
    return (
      <>
        <LegacyFgfr2Files downloads={downloads} />
        <AvailabilityManifest availability={availability} />
      </>
    );
  }
  return (
    <>
      <PackageBuilder downloads={downloads} />
      <AvailabilityManifest availability={availability} />
    </>
  );
}

/** Human-readable names for the canonical availability states. */
const STATUS_LABEL = {
  available: "Available",
  not_applicable: "Not applicable",
  scientifically_unavailable: "Not supported by the annotation",
  pending: "Pending",
  technically_missing: "Missing",
  stale: "Out of date",
  failed: "Failed",
};

/**
 * What each analysis produced, and for anything absent, the prerequisite that decided it.
 *
 * A status record rather than a data product: an analysis that does not apply gets a row
 * here instead of an empty scientific table written solely so the run looks complete.
 */
function AvailabilityManifest({ availability }) {
  const analyses = availability?.analyses || [];
  if (!analyses.length) return null;
  return (
    <div className="card availability-manifest">
      <div className="card-head">
        <h3>Analysis availability</h3>
        <Badge cls={availability.ready ? "accepted" : "info"} soft>
          {availability.ready ? "All applicable analyses resolved" : "In progress"}
        </Badge>
      </div>
      <table className="table small">
        <thead>
          <tr>
            <th>Analysis</th><th>Status</th><th>Prerequisite</th>
            <th className="num">Count</th><th>Reason</th>
          </tr>
        </thead>
        <tbody>
          {analyses.map((a) => (
            <tr key={a.analysis_name}>
              <td><code>{a.analysis_name}</code></td>
              <td>
                <span className={a.status === "not_applicable" ? "muted" : ""}>
                  {STATUS_LABEL[a.status] || a.status}
                </span>
              </td>
              <td><code>{a.prerequisite_name || "—"}</code></td>
              <td className="num">{a.prerequisite_count ?? "—"}</td>
              <td className="muted small">{a.user_message || a.reason_code || "—"}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function LegacyFgfr2Files({ downloads }) {
  const list = Array.isArray(downloads) ? downloads : [];
  const truth = list.find((d) => d.label?.includes("truth"));
  const manifest = list.find((d) => d.label?.includes("manifest"));
  const primaryFasta = list.find((d) => d.label?.includes("Primary protein"));
  const legacyItems = [primaryFasta, manifest, truth].filter(Boolean);
  return (
    <div className="files-tab">
      <div className="file-list">
        {legacyItems.map((d) => (
          <a key={d.path} className="file-row" href={fileUrl(d.path)}>
            <span>{d.label}</span>
            <small>{d.size_human}</small>
          </a>
        ))}
      </div>
    </div>
  );
}

function humanSize(n) {
  if (!n) return "—";
  if (n < 1024) return `${n} B`;
  if (n < 1024 * 1024) return `${(n / 1024).toFixed(1)} KB`;
  return `${(n / (1024 * 1024)).toFixed(1)} MB`;
}

const BUILD_STATES = {
  idle: "Ready to build",
  validating: "Validating selection",
  building: "Building package",
  ready: "Package ready",
  failed: "Build failed",
  expired: "Package expired",
};

function PackageBuilder({ downloads }) {
  const selection = useScientificSelection();
  const [caps, setCaps] = useState(null);
  // Linked selection: the Scope chosen in the Figure Gallery opens this page on
  // the matching species. The backend falls back to its own default if the run
  // does not have that Scope, so an unknown value is harmless.
  const [scope, setScope] = useState(() => {
    const figScope = selection?.selectedFigureScope;
    return figScope && figScope !== "comparative" ? figScope : "";
  });
  const [preset, setPreset] = useState("recommended");
  const [customItems, setCustomItems] = useState([]);
  const [selectedSpecies, setSelectedSpecies] = useState([]);
  const [showReview, setShowReview] = useState(false);
  const [showDeps, setShowDeps] = useState(false);
  const [job, setJob] = useState(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");

  // Scope drives everything: a new Scope means a new capability object.
  useEffect(() => {
    let cancelled = false;
    api.packageCapabilities(scope || undefined)
      .then((c) => {
        if (cancelled) return;
        setCaps(c);
        setSelectedSpecies(c.selected_species || []);
        setError("");
      })
      .catch((e) => { if (!cancelled) setError(e.message || String(e)); });
    return () => { cancelled = true; };
  }, [scope]);

  const activePreset = caps?.presets?.[preset];
  // Both feed the dependency resolution below, so they need a stable identity across
  // renders — a fresh `{}` or `[]` each render would re-resolve on every keystroke.
  const items = useMemo(() => caps?.items || {}, [caps]);
  const requested = useMemo(
    () => (preset === "custom" ? customItems : (activePreset?.items || [])),
    [preset, customItems, activePreset]
  );

  // Dependencies are resolved against the capability object, so an unavailable
  // dependency can never sneak into the selection.
  const resolved = useMemo(() => {
    const seen = new Set();
    const out = [];
    const add = (id) => {
      const it = items[id];
      if (!it || seen.has(id) || !it.available) return;
      seen.add(id);
      for (const dep of it.dependencies || []) add(dep);
      if (!out.includes(id)) out.push(id);
    };
    for (const id of requested) add(id);
    return out;
  }, [requested, items]);

  const dependencyOnly = resolved.filter((id) => !requested.includes(id));
  const estimate = resolved.reduce((s, id) => s + (items[id]?.estimated_bytes || 0), 0);
  const state = busy ? (job?.status || "building") : (job?.status || "idle");

  const toggleItem = (id) => {
    if (!items[id]?.available) return;
    setCustomItems((prev) => {
      const base = prev.length || preset === "custom" ? prev : (activePreset?.items || []);
      return base.includes(id) ? base.filter((x) => x !== id) : [...base, id];
    });
    setPreset("custom");
  };

  const build = async () => {
    setBusy(true);
    setError("");
    setJob({ status: "validating", message: "Validating selection" });
    try {
      const res = await api.createPackage({
        preset, scope, items: resolved, species: selectedSpecies,
      });
      setJob(res);
      let n = 0;
      while (res.job_id && res.status === "building" && n < 40) {
        await new Promise((r) => setTimeout(r, 250));
        const st = await api.packageStatus(res.job_id);
        setJob(st);
        if (st.status === "ready" || st.status === "failed") break;
        n += 1;
      }
    } catch (e) {
      setError(e.message || String(e));
      setJob({ status: "failed", error: e.message || String(e) });
    } finally {
      setBusy(false);
    }
  };

  if (error && !caps) {
    return (
      <div className="files-tab">
        <p className="error">{error}</p>
        <ProvenanceFallback downloads={downloads} />
      </div>
    );
  }
  if (!caps) {
    return <div className="files-tab"><p className="muted">Loading downloads…</p></div>;
  }

  const groups = caps.groups || [];
  const scopeLabel = (caps.scopes || []).find((s) => s.id === caps.scope)?.label
    || caps.scope;

  return (
    <div className="files-tab data-downloads">
      <div className="dd-head">
        <div>
          <h3 className="dd-title">Data &amp; Downloads</h3>
          <p className="muted sm">
            {caps.multi_species
              ? `Scientific downloads for ${caps.gene_symbol}. Choose a Scope: the `
                + "comparative products, all species, or one species."
              : `Scientific downloads for ${caps.gene_symbol} · ${scopeLabel}. `
                + "This run has one species, so it has no comparative products."}
          </p>
        </div>
        {caps.multi_species && (
          <label className="dd-scope">
            <span className="muted sm">Scope</span>
            <select value={caps.scope} onChange={(e) => { setScope(e.target.value); setJob(null); setPreset("recommended"); setCustomItems([]); }}
              aria-label="Package scope">
              {caps.scopes.map((s) => (
                <option key={s.id} value={s.id}>{s.label}</option>
              ))}
            </select>
          </label>
        )}
      </div>

      <DirectDownloads caps={caps} />

      <div className="dd-builder">
        <h4 className="dd-sub">Package builder</h4>

        <div className="dd-presets">
          {Object.values(caps.presets || {}).map((p) => (
            <button key={p.id}
              className={`btn sm ${preset === p.id ? "primary" : "ghost"}`}
              title={p.description}
              onClick={() => { setPreset(p.id); setCustomItems(p.items || []); setJob(null); }}>
              {p.label}
              {p.id !== "custom" && (
                <span className="muted sm"> · {p.items.length}</span>
              )}
            </button>
          ))}
        </div>
        {activePreset?.description && preset !== "custom" && (
          <p className="muted sm dd-preset-note">
            {activePreset.description}
            {activePreset.unavailable_items?.length > 0 && (
              <> {activePreset.unavailable_items.length} item(s) in this preset are
                not available for this run and were left out.</>
            )}
          </p>
        )}

        {caps.multi_species && caps.scope === "all" && (
          <div className="dd-species">
            <span className="muted sm">Species in this package</span>
            <div className="dd-species-list">
              {caps.species.map((r) => {
                const checked = selectedSpecies.includes(r.species_id);
                return (
                  <label key={r.species_id} className="dd-species-row">
                    <input type="checkbox" checked={checked}
                      onChange={() => setSelectedSpecies((prev) => (
                        checked ? prev.filter((x) => x !== r.species_id)
                          : [...prev, r.species_id]))} />
                    <span>{r.scientific_name || r.species_id}</span>
                    <Badge cls={r.analysis_status === "available" ? "accepted" : "neutral"} soft>
                      {r.analysis_status || "unknown"}
                    </Badge>
                  </label>
                );
              })}
            </div>
          </div>
        )}

        <div className="dd-summary">
          {groups.map((g) => {
            const chosen = g.items.filter((i) => resolved.includes(i));
            return (
              <div key={g.id} className="dd-summary-chip">
                <b>{g.label}</b>
                <span className="muted sm">
                  {chosen.length} / {g.items.filter((i) => items[i]?.available).length}
                </span>
              </div>
            );
          })}
          <div className="dd-summary-chip">
            <b>Estimated size</b>
            <span className="muted sm">≈ {humanSize(estimate)}</span>
          </div>
        </div>

        {preset === "custom" && (
          <div className="dd-items">
            {groups.map((g) => (
              <div key={g.id} className="dd-group">
                <h5 className="file-group-head">{g.label}</h5>
                {g.items.map((id) => {
                  const it = items[id];
                  const on = resolved.includes(id);
                  const depOnly = on && !customItems.includes(id);
                  return (
                    <label key={id}
                      className={`dd-item${it.available ? "" : " unavailable"}`}
                      title={it.available ? (it.description || it.label) : it.reason}>
                      <input type="checkbox" checked={on && it.available}
                        disabled={!it.available}
                        onChange={() => toggleItem(id)} />
                      <span>
                        {it.label}
                        {depOnly && <em className="muted sm"> · required by another item</em>}
                      </span>
                      {it.available
                        ? <span className="muted sm">{humanSize(it.estimated_bytes)}</span>
                        : <Badge cls="neutral" soft>unavailable</Badge>}
                      {!it.available && (
                        <span className="dd-reason muted sm">{it.reason}</span>
                      )}
                    </label>
                  );
                })}
              </div>
            ))}
          </div>
        )}

        <details className="dd-deps" open={showDeps}
          onToggle={(e) => setShowDeps(e.currentTarget.open)}>
          <summary>
            Dependency resolution · {resolved.length} file group(s)
            {dependencyOnly.length > 0 && ` · ${dependencyOnly.length} added automatically`}
          </summary>
          <ul className="dd-dep-list">
            {resolved.map((id) => (
              <li key={id} className="muted sm">
                {items[id]?.label || id}
                {dependencyOnly.includes(id) && " · added as a dependency"}
              </li>
            ))}
          </ul>
        </details>

        <div className="dd-actions">
          <button className="btn ghost" onClick={() => setShowReview((v) => !v)}>
            {showReview ? "Hide contents" : "Review contents"}
          </button>
          <button className="btn primary" disabled={busy || !resolved.length}
            onClick={build}>
            {busy ? "Building package…" : "Build package"}
          </button>
          <Badge cls={state === "ready" ? "accepted" : state === "failed" ? "rejected" : "neutral"} soft>
            {BUILD_STATES[state] || state}
          </Badge>
        </div>

        {busy && (
          <div className="dd-progress">
            <div className="dd-progress-bar"
              style={{ width: `${Math.round((job?.progress || 0.1) * 100)}%` }} />
            <span className="muted sm">{job?.message || "Working…"}</span>
          </div>
        )}

        {showReview && (
          <div className="dd-review">
            <h5>Package contents</h5>
            <ul>
              {resolved.map((id) => (
                <li key={id}>
                  {items[id]?.label}
                  <span className="muted sm"> · {items[id]?.group_label}
                    {" · "}{humanSize(items[id]?.estimated_bytes)}</span>
                </li>
              ))}
            </ul>
            {Object.values(items).some((i) => !i.available) && (
              <>
                <h5>Not available for this run</h5>
                <ul>
                  {Object.values(items).filter((i) => !i.available).map((i) => (
                    <li key={i.id} className="muted sm">{i.label}: {i.reason}</li>
                  ))}
                </ul>
              </>
            )}
          </div>
        )}

        {job?.status === "ready" && (
          <div className="dd-ready">
            <h5>{job.package_name}</h5>
            <p className="muted sm">
              {job.n_files} file(s) · {humanSize(job.estimated_bytes)}
              {job.updated_at
                && ` · generated ${new Date(job.updated_at).toLocaleString()}`}
            </p>
            {job.warnings?.length > 0 && job.warnings.map((w, i) => (
              <p key={i} className="warn sm">{w}</p>
            ))}
            <div className="dd-actions">
              <a className="btn primary" href={fileUrl(job.zip_path)}>Download ZIP</a>
              <button className="btn ghost" onClick={() => setJob(null)}>
                Build another package
              </button>
            </div>
            {job.manifest?.omitted?.length > 0 && (
              <details className="dd-omitted">
                <summary>Not included ({job.manifest.omitted.length})</summary>
                <ul>
                  {job.manifest.omitted.map((o, i) => (
                    <li key={i} className="muted sm">
                      {o.label || o.item}: {o.reason}
                    </li>
                  ))}
                </ul>
              </details>
            )}
          </div>
        )}
        {job?.status === "failed" && (
          <p className="error sm">{job.error || "Package build failed."}</p>
        )}
        {job?.status === "expired" && (
          <p className="warn sm">This package has expired. Build it again.</p>
        )}
        {error && <p className="error sm">{error}</p>}
      </div>

      <details className="dd-provenance">
        <summary>Individual provenance files</summary>
        <ProvenanceFallback downloads={downloads} />
      </details>
    </div>
  );
}

/** Direct, one-click downloads for the available items of the current Scope. */
function DirectDownloads({ caps }) {
  const groups = caps.groups || [];
  const items = caps.items || {};
  return (
    <div className="dd-direct">
      <h4 className="dd-sub">Direct downloads</h4>
      {groups.map((g) => {
        const rows = g.items.map((id) => items[id]).filter((i) => i && i.available);
        if (!rows.length) return null;
        return (
          <div key={g.id} className="dd-group">
            <h5 className="file-group-head">{g.label}</h5>
            {rows.map((it) => (
              it.path ? (
                <a key={it.id} className="file-row" href={fileUrl(it.path)}
                  title={it.description || it.label}>
                  <span>{it.label}</span>
                  <small>{humanSize(it.estimated_bytes)}</small>
                </a>
              ) : (
                <div key={it.id} className="file-row muted"
                  title={it.description || it.label}>
                  <span>{it.label}</span>
                  <small>
                    {it.n_files > 1 ? `${it.n_files} files · package only`
                      : "package only"}
                  </small>
                </div>
              )
            ))}
          </div>
        );
      })}
    </div>
  );
}

function ProvenanceFallback({ downloads }) {
  const all = Array.isArray(downloads)
    ? downloads
    : (Array.isArray(downloads?.items) ? downloads.items : []);
  if (!all.length) {
    return <p className="muted">No individual artefacts published yet.</p>;
  }
  return (
    <div className="prov-list">
      {all.filter((d) => d.path).slice(0, 40).map((d) => (
        <a key={d.path} className="file-row" href={fileUrl(d.path)}>
          <span>{d.label || d.name}</span>
          <small>{d.size_human}</small>
        </a>
      ))}
    </div>
  );
}
