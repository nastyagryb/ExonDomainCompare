/* eslint-disable react-refresh/only-export-components */
// Shared ExonDomainCompare UI components.
//
// These are the ACTUAL building blocks of the FGFR2 pages, extracted verbatim so
// that FGFR2 and any other shared-pipeline gene (e.g. FGFR1) render through the
// *same* React components — not look-alike copies. FGFR2 keeps its validated
// event-specific sections by passing them as children/slots; other genes leave
// those slots empty and fill the exploratory-evidence slot instead.
//
// Rule: no gene-specific vocabulary lives here. Callers pass data + labels.
import { Fragment, useState } from "react";
import { Badge, Kpi, Field } from "../../ui";
import { fileUrl } from "../../api";
import { textProps } from "../../pages/viewers/semanticStyles";

// Label ink comes from the shared scientific specification and is written onto the
// marks as explicit attributes; the stylesheet keeps size and layout. A label that
// is painted only by CSS turns black as soon as the SVG leaves the page.
const AXIS_TEXT = textProps("axis");
const CAND_TEXT = textProps("candidateLabel");

// --------------------------------------------------------------------------- //
// Page + layout primitives (extracted from Overview.jsx / GeneExplorer.jsx)
// --------------------------------------------------------------------------- //
export function DatasetPageHeader({ eyebrow, title, badges }) {
  return (
    <div className="page-head">
      <div>
        {eyebrow && <p className="eyebrow">{eyebrow}</p>}
        <h2>{title}</h2>
      </div>
      {badges && <div className="head-badges">{badges}</div>}
    </div>
  );
}

export function KpiGrid({ items }) {
  return (
    <div className="kpi-grid">
      {items.map((it) => {
        const [label, value, sub, cls] = Array.isArray(it)
          ? it : [it.label, it.value, it.sub, it.cls];
        return <Kpi key={label} label={label} value={value} sub={sub} cls={cls} />;
      })}
    </div>
  );
}

// The gene-page frame: <section.page.gene-page> <div.gene-layout> sidebar + body.
export function GeneExplorerShell({ sidebar, children }) {
  return (
    <section className="page gene-page">
      <div className="gene-layout">
        {sidebar}
        {children}
      </div>
    </section>
  );
}

export function SpeciesPanel({ filters, children }) {
  return (
    <aside className="species-panel">
      {filters}
      <div className="species-list">{children}</div>
    </aside>
  );
}

// Single species button in the left panel (identical DOM to the FGFR2 card).
export function SpeciesCard({ title, badge, sub, extra, selected, onClick }) {
  return (
    <button className={`species-card${selected ? " sel" : ""}`} onClick={onClick}>
      <div className="sc-top">
        <b>{title}</b>
        {badge}
      </div>
      <div className="sc-bottom">
        <small>{sub}</small>
        {extra}
      </div>
    </button>
  );
}

export function WorkspaceHeader({ title, sub, badges }) {
  return (
    <div className="ws-header">
      <div>
        <h2>{title}</h2>
        {sub && <span className="muted">{sub}</span>}
      </div>
      {badges && <div className="ws-badges">{badges}</div>}
    </div>
  );
}

// Tab bar. `tabs` items: [id, label] or {id, label, pending, separatorBefore}.
export function ExplorerTabs({ tabs, active, onSelect }) {
  return (
    <div className="tabs">
      {tabs.map((t) => {
        const id = Array.isArray(t) ? t[0] : t.id;
        const label = Array.isArray(t) ? t[1] : t.label;
        const pending = Array.isArray(t) ? false : t.pending;
        const sep = Array.isArray(t) ? false : t.separatorBefore;
        return (
          <Fragment key={id}>
            {sep && <span className="tab-sep" />}
            <button className={active === id ? "tab sel" : "tab"} onClick={() => onSelect(id)}>
              {label}{pending && <span className="tab-pending"> ·</span>}
            </button>
          </Fragment>
        );
      })}
    </div>
  );
}

// --------------------------------------------------------------------------- //
// Evidence summary (extracted from GeneExplorer.jsx SummaryCard/SummaryTab)
// --------------------------------------------------------------------------- //
export function SummaryCard({ label, cls, value, note, soft, children }) {
  return (
    <div className="summary-card">
      <span className="mini-label">{label}</span>
      <div className="sc-badge">
        <Badge cls={cls} soft={soft}>{value || "—"}</Badge>
        {children}
      </div>
      {note && <p className="sc-note">{note}</p>}
    </div>
  );
}

// The `.summary-cards` grid. `items`: {label,cls,value,note,soft} plus optional
// `extra` node appended inside the card (used for "View →" links).
export function EvidenceSummary({ items, footer, children }) {
  return (
    <>
      <div className="summary-cards">
        {items?.map((it) => (
          <SummaryCard key={it.id || it.label} label={it.label} cls={it.cls}
                       value={it.value} note={it.note} soft={it.soft}>
            {it.extra}
          </SummaryCard>
        ))}
        {children}
      </div>
      {footer}
    </>
  );
}

// Evidence rows list (extracted from GeneExplorer.jsx EvidenceTab markup).
export function EvidenceRow({ label, badge, note, onClick, selected = false }) {
  return (
    <div className={`evidence-row${selected ? " selected" : ""}`}
      role={onClick ? "button" : undefined} tabIndex={onClick ? 0 : undefined}
      onClick={onClick} onKeyDown={onClick ? (e) => e.key === "Enter" && onClick() : undefined}>
      <span className="er-label">{label}</span>
      <span className="er-status">{badge}</span>
      {note && <span className="er-note">{note}</span>}
    </div>
  );
}

// --------------------------------------------------------------------------- //
// Figure gallery (extracted from FigureGallery.jsx FigCard/fig-grid)
// --------------------------------------------------------------------------- //
export function FigureCardGrid({ children }) {
  return <div className="fig-grid">{children}</div>;
}

// Flexible figure card matching the FGFR2 fig-card DOM: thumb (button/link/node)
// + meta (title-row + badge, caption, actions row).
export function FigureCard({ thumb, thumbHref, thumbOnClick, thumbStyle,
  title, badge, caption, actions }) {
  let thumbEl;
  if (thumbHref) {
    thumbEl = <a className="fig-thumb" href={thumbHref} target="_blank" rel="noreferrer">{thumb}</a>;
  } else if (thumbOnClick) {
    thumbEl = <button className="fig-thumb" onClick={thumbOnClick}>{thumb}</button>;
  } else {
    thumbEl = <div className="fig-thumb" style={thumbStyle}>{thumb}</div>;
  }
  return (
    <article className="fig-card">
      {thumbEl}
      <div className="fig-meta">
        <div className="fig-title-row"><b>{title}</b>{badge}</div>
        {caption && <p className="fig-cap">{caption}</p>}
        {actions && <div className="fig-actions">{actions}</div>}
      </div>
    </article>
  );
}

// --------------------------------------------------------------------------- //
// Pending analysis gate (shared stage gate; used before the cluster step)
// --------------------------------------------------------------------------- //
export function PendingAnalysisCard({ title, badge, description, command, children }) {
  const [copied, setCopied] = useState(false);
  const copy = async () => {
    try { await navigator.clipboard.writeText(command); setCopied(true); setTimeout(() => setCopied(false), 1500); }
    catch { /* ignore */ }
  };
  return (
    <div className="card">
      <div className="card-head">
        <h3>{title}</h3>
        {badge}
      </div>
      {description && <p className="muted small">{description}</p>}
      {command && (
        <div className="run-cta">
          <code style={{ userSelect: "all" }}>{command}</code>
          <button className="btn ghost sm" onClick={copy}>{copied ? "Copied" : "Copy"}</button>
        </div>
      )}
      {children}
    </div>
  );
}

// --------------------------------------------------------------------------- //
// Protein isoform table (shared; used by the generic Gene Explorer)
// --------------------------------------------------------------------------- //
// A model row has meaningful details only when at least one real (non-placeholder)
// field beyond the base table columns is present. Never render an empty drawer.
function hasIsoformDetails(iso) {
  return [iso.exon_count, iso.selection_reason, iso.completeness,
    iso.diff_from_primary, iso.affected_candidates]
    .some((v) => v != null && v !== "" && v !== "—");
}

export function ProteinIsoformTable({ isoforms, selectionMethod, expandable }) {
  const [open, setOpen] = useState({});
  if (!isoforms?.length) return <p className="muted">No protein isoforms.</p>;
  const toggle = (id) => setOpen((o) => ({ ...o, [id]: !o[id] }));
  return (
    <table className="mini-table">
      <thead><tr>
        <th>Protein</th><th>Transcript</th><th>Length (aa)</th><th>Source</th><th>Role</th><th>Status</th><th>Primary</th>
        {expandable && <th />}
      </tr></thead>
      <tbody>
        {isoforms.map((iso, i) => {
          const primary = iso.is_primary === true
            || String(iso.primary_status || "").toLowerCase() === "primary";
          const curated = String(iso.curation_status || "").toLowerCase() === "curated";
          const predicted = String(iso.curation_status || "").toLowerCase() === "predicted";
          const pid = iso.protein_id || `row-${i}`;
          const isOpen = open[pid];
          // Only offer Details when there is real content to show.
          const showDetails = expandable && hasIsoformDetails(iso);
          return (
            <Fragment key={pid}>
            <tr className={primary ? "row-primary" : ""}>
              <td>{iso.protein_id}{primary ? " ★" : ""}</td>
              <td>{iso.transcript_id || "—"}</td>
              <td>{iso.protein_length ?? iso.length_aa ?? "—"}</td>
              <td>{iso.curation_status
                ? <Badge cls={curated ? "accepted" : predicted ? "neutral" : "minor"} soft>
                    {curated ? "curated RefSeq" : predicted ? "predicted" : iso.curation_status}
                  </Badge>
                : "—"}</td>
              <td><Badge cls={primary ? "accepted" : "neutral"} soft>{primary ? "primary" : "alternative"}</Badge></td>
              <td className="muted small">{iso.status || (primary ? "primary" : "alternative")}</td>
              <td className="muted small">{primary ? (iso.selection_reason || selectionMethod || "—") : ""}</td>
              {expandable && (
                <td>{showDetails
                  ? <button className="btn ghost sm" onClick={() => toggle(pid)}>{isOpen ? "Hide" : "Details"}</button>
                  : null}</td>
              )}
            </tr>
            {showDetails && isOpen && (
              <tr key={`${pid}-detail`} className="row-detail">
                <td colSpan={8}>
                  <div className="field-grid compact">
                    <Field label="Transcript ID"><code>{iso.transcript_id || "—"}</code></Field>
                    <Field label="Protein ID"><code>{iso.protein_id || "—"}</code></Field>
                    <Field label="Curated / predicted">{curated ? "curated" : predicted ? "predicted" : "—"}</Field>
                    <Field label="Role">{primary ? "primary" : "alternative"}</Field>
                    <Field label="Protein length">{iso.protein_length ?? iso.length_aa
                      ? `${iso.protein_length ?? iso.length_aa} aa` : "—"}</Field>
                    <Field label="Exon count">{iso.exon_count ?? "—"}</Field>
                    <Field label="Model completeness">{iso.completeness || "—"}</Field>
                    <Field label="Primary-selection reason" wide>{iso.selection_reason || "—"}</Field>
                    <Field label="Difference from primary">{iso.diff_from_primary || "—"}</Field>
                    <Field label="Affected candidate regions">{iso.affected_candidates || "—"}</Field>
                  </div>
                </td>
              </tr>
            )}
            </Fragment>
          );
        })}
      </tbody>
    </table>
  );
}

// --------------------------------------------------------------------------- //
// Exon/protein architecture track (SVG). Gene-agnostic; exploratory candidate
// regions are overlays only, never validated events.
// --------------------------------------------------------------------------- //
export function ProteinArchitectureTrack({ proteins, primaryOnly }) {
  const [sel, setSel] = useState(null);
  if (!proteins?.length) return null;
  const primary = proteins.find((p) => p.role === "primary") || proteins[0];
  const active = sel ? proteins.find((p) => p.protein_id === sel) || primary : primary;
  return (
    <>
      <ProteinTrack protein={active} highlight={active.role === "primary"} />
      {active.candidate_regions?.length > 0 && (
        <p className="muted small" style={{ marginTop: 6 }}>
          Orange overlays are <b>exploratory</b> isoform-difference candidate regions on the primary
          protein — not validated events.
        </p>
      )}
      {!primaryOnly && proteins.length > 1 && (
        <div className="chip-row" style={{ marginTop: 8 }}>
          <span className="muted small" style={{ marginRight: 6 }}>Isoform:</span>
          {proteins.map((p) => (
            <button key={p.protein_id}
                    className={`chip${(active.protein_id === p.protein_id) ? " sel" : ""}`}
                    onClick={() => setSel(p.protein_id)}>
              {p.protein_id}{p.role === "primary" ? " ★" : ""}
            </button>
          ))}
        </div>
      )}
    </>
  );
}

function ProteinTrack({ protein, highlight }) {
  const W = 900, padL = 6, padR = 6;
  const maxLen = Math.max(protein.length_aa || 0,
    ...(protein.exons || []).map((e) => e.protein_end_aa || 0), 1);
  const exonY = 30, exonH = 20, candY = 8, candH = 14;
  const domains = protein.domains || [];
  const tms = protein.tm_regions || [];
  const hasDom = domains.length > 0 || tms.length > 0;
  const domY = 54, domH = 12;
  const axisY = hasDom ? 74 : 58;
  const H = axisY + 16;
  const xs = (aa) => padL + (Math.max(0, (aa - 1)) / maxLen) * (W - padL - padR);
  const ticks = [];
  const step = maxLen > 800 ? 200 : maxLen > 300 ? 100 : 50;
  for (let t = 0; t <= maxLen; t += step) ticks.push(t);
  return (
    <div className="track-wrap">
      <svg viewBox={`0 0 ${W} ${H}`} className="exon-track" preserveAspectRatio="xMidYMid meet"
           role="img" aria-label={`Exon/protein track for ${protein.protein_id}`}>
        <rect x={xs(1)} y={exonY + exonH / 2 - 1} width={xs(maxLen) - xs(1)} height={2} fill="#cbd5e1" />
        {(protein.exons || []).map((e, i) => {
          const x = xs(e.protein_start_aa || 1);
          const w = Math.max(1.5, xs((e.protein_end_aa || e.protein_start_aa || 1) + 1) - x);
          return (
            <rect key={i} x={x} y={exonY} width={w} height={exonH} rx={2}
                  fill={i % 2 ? "#93a4bd" : "#c3cede"}
                  stroke={highlight ? "#475569" : "#94a3b8"} strokeWidth={highlight ? 1 : 0.6}>
              <title>{`Exon ${e.exon_number ?? i + 1}: aa ${e.protein_start_aa}–${e.protein_end_aa}${e.confidence ? ` (${e.confidence})` : ""}`}</title>
            </rect>
          );
        })}
        {(protein.candidate_regions || []).map((c, i) => {
          const x = xs(c.start_aa || 1);
          const w = Math.max(2, xs((c.end_aa || c.start_aa || 1) + 1) - x);
          return (
            <g key={i}>
              <rect x={x} y={candY} width={w} height={candH} rx={2}
                    fill="rgba(234,88,12,0.30)" stroke="#ea580c" strokeWidth={1} strokeDasharray="3 2">
                <title>{`Exploratory candidate ${c.start_aa}–${c.end_aa} aa · support ${c.support_count ?? "?"} · ${c.confidence || ""}${c.exon_aligned ? " · exon-aligned" : ""} (not a validated event)`}</title>
              </rect>
              <text x={x + 2} y={candY - 1} className="track-cand-label"
                    fill={CAND_TEXT.fill} fontSize={CAND_TEXT.fontSize}
                    fontWeight={CAND_TEXT.fontWeight}>{c.start_aa}–{c.end_aa}</text>
            </g>
          );
        })}
        {domains.map((d, i) => {
          const x = xs(d.start_aa || 1);
          const w = Math.max(2, xs((d.end_aa || d.start_aa || 1) + 1) - x);
          return (
            <rect key={`d${i}`} x={x} y={domY} width={w} height={domH} rx={2}
                  fill="rgba(13,148,136,0.35)" stroke="#0d9488" strokeWidth={0.8}>
              <title>{`${d.domain_name || d.domain_id} (${d.domain_source}) ${d.start_aa}–${d.end_aa} aa`}</title>
            </rect>
          );
        })}
        {tms.map((t, i) => {
          const x = xs(t.start_aa || 1);
          const w = Math.max(2, xs((t.end_aa || t.start_aa || 1) + 1) - x);
          return (
            <rect key={`t${i}`} x={x} y={domY} width={w} height={domH} rx={2}
                  fill="rgba(124,58,237,0.35)" stroke="#7c3aed" strokeWidth={0.8}>
              <title>{`TM helix ${t.start_aa}–${t.end_aa} aa`}</title>
            </rect>
          );
        })}
        {ticks.map((t) => (
          <g key={`ax${t}`}>
            <line x1={xs(t || 1)} y1={axisY} x2={xs(t || 1)} y2={axisY + 4} stroke="#94a3b8" strokeWidth={0.8} />
            <text x={xs(t || 1)} y={axisY + 13} className="track-axis-label"
                  fill={AXIS_TEXT.fill} fontSize={AXIS_TEXT.fontSize}>{t}</text>
          </g>
        ))}
      </svg>
      <div className="track-legend muted small">
        <span><i className="lg-exon" /> exon</span>
        {(protein.candidate_regions || []).length > 0 && <span><i className="lg-cand" /> exploratory candidate</span>}
        {domains.length > 0 && <span><i className="lg-dom" /> domain</span>}
        {tms.length > 0 && <span><i className="lg-tm" /> TM helix</span>}
        {!hasDom && <span className="muted">domain context pending cluster</span>}
        <span className="track-meta">{protein.protein_id} · {protein.length_aa || "?"} aa · {(protein.exons || []).length} exons</span>
      </div>
    </div>
  );
}

// The synteny neighbourhood strip that used to live here was a second drawing
// implementation with its own ordering, styling and target handling. It is gone:
// every synteny view now renders through pages/viewers/SyntenyNeighbourhood.jsx
// on the shared `shared_synteny_v1` contract.

export const fileUrlShared = fileUrl;
