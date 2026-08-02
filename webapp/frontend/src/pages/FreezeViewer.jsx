import { useEffect, useState } from "react";
import { api, fileUrl } from "../api";
import { Badge, Spinner } from "../ui";

export default function FreezeViewer() {
  const [freeze, setFreeze] = useState(null);
  const [downloads, setDownloads] = useState([]);

  useEffect(() => {
    api.freeze().then(setFreeze).catch(() => setFreeze(null));
    api.downloads().then(setDownloads).catch(() => setDownloads([]));
  }, []);

  if (!freeze) return <section className="page"><Spinner /></section>;

  const reports = downloads.filter((d) => d.group === "Reports");
  const archive = downloads.filter((d) => d.group === "Archive");

  return (
    <section className="page">
      <div className="page-head">
        <div>
          <p className="eyebrow">Reproducibility & InterPro-ready export</p>
          <h2>Freeze Viewer</h2>
        </div>
        <div className="head-badges">
          <Badge cls={freeze.run_mode_label === "full clean run" ? "accepted" : "info"} soft>
            {freeze.run_mode_label}
          </Badge>
          <Badge cls={freeze.checksum_gate === "pass" ? "accepted" : "neutral"} soft>
            checksums {freeze.checksum_gate} · {freeze.checksum_count}
          </Badge>
        </div>
      </div>

      {freeze.interpro_policy && (
        <div className="card interpro-policy">
          <div className="ip-head">
            <h3>InterPro input policy</h3>
            <div className="ip-badges">
              <Badge cls="accepted" soft>Primary · {freeze.interpro_policy.primary_count} · main</Badge>
              <Badge cls="minor" soft>Review-included · {freeze.interpro_policy.review_included_count} · optional</Badge>
            </div>
          </div>
          <ul className="ip-notes">
            {freeze.interpro_policy.notes.map((n, i) => <li key={i}>{n}</li>)}
          </ul>
        </div>
      )}

      <div className="freeze-grid">
        {freeze.cards.map((c) => (
          <article key={c.key} className={`freeze-card${c.available ? "" : " na"}${c.key === "primary_fasta" ? " fc-primary" : ""}`}>
            <div className="fc-head">
              <h3>{c.title}</h3>
              <Badge cls={c.key === "primary_fasta" ? "accepted" : c.key === "review_fasta" ? "minor" : c.available ? "accepted" : "neutral"} soft>
                {c.key === "primary_fasta" ? "main input" : c.key === "review_fasta" ? "optional" : c.available ? "available" : "missing"}
              </Badge>
            </div>
            <p className="fc-role">{c.role}</p>
            <div className="fc-meta">
              {c.sequences != null && <span><b>{c.sequences}</b> sequences</span>}
              <span>{c.size_human}</span>
            </div>
            <code className="fc-file">{c.name}</code>
            {c.available && (
              <a className="btn primary sm" href={fileUrl(c.path)}>Download</a>
            )}
          </article>
        ))}
      </div>

      <div className="freeze-cols">
        <div className="card">
          <h3>Final reports</h3>
          <div className="file-list">
            {reports.map((r) => (
              <a key={r.path} className="file-row" href={fileUrl(r.path)}>
                <span>{r.label}</span><small>{r.size_human}</small>
              </a>
            ))}
            {reports.length === 0 && <p className="muted">No reports.</p>}
          </div>
        </div>
        <div className="card">
          <h3>Final archive</h3>
          <div className="file-list">
            {archive.map((a) => (
              <a key={a.path} className="file-row" href={fileUrl(a.path)}>
                <span>{a.label}</span><small>{a.size_human}</small>
              </a>
            ))}
            {archive.length === 0 && <p className="muted">No archive.</p>}
          </div>
          <div className="run-mode-flags">
            <Flag ok={freeze.run_mode.full_clean_run_completed} label="full_clean_run_completed" />
            <Flag ok={!freeze.run_mode.used_cached_v3_outputs} label="fresh v3 outputs" />
            <Flag ok={!freeze.run_mode.used_cached_msa_outputs} label="fresh MSA outputs" />
          </div>
        </div>
      </div>
    </section>
  );
}

function Flag({ ok, label }) {
  return (
    <span className="flag">
      <span className={`cell-dot st-${ok ? "accepted" : "minor"}`} />{label}
    </span>
  );
}
