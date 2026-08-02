import { useEffect, useRef, useState } from "react";
import {
  TONE_GLYPH, TONE_LABEL, statusClass,
} from "./uiStatus";

export function ToneMark({ tone }) {
  if (!tone || !TONE_GLYPH[tone]) return null;
  return <span className={`tone-chip tone-${tone}`} title={TONE_LABEL[tone]}>{TONE_GLYPH[tone]}</span>;
}

export function Badge({ children, cls, title, soft }) {
  const c = cls || statusClass(children);
  return (
    <span className={`badge st-${c}${soft ? " soft" : ""}`} title={title}>
      {children}
    </span>
  );
}

export function IsoBadge({ iso, active = true, children }) {
  const t = String(iso || "").toLowerCase();
  const c = t.includes("iiib") ? "iiib" : t.includes("iiic") ? "iiic" : "neutral";
  return <span className={`iso iso-${c}${active ? "" : " off"}`}>{children || iso}</span>;
}

export function Dot({ cls }) {
  return <span className={`dot st-${cls || "neutral"}`} />;
}

export function Kpi({ label, value, sub, cls }) {
  return (
    <div className={`kpi${cls ? " kpi-" + cls : ""}`}>
      <span className="kpi-label">{label}</span>
      <strong className="kpi-value">{value}</strong>
      {sub && <em className="kpi-sub">{sub}</em>}
    </div>
  );
}

export function Drawer({ open, onClose, title, subtitle, children, footer }) {
  useEffect(() => {
    function onKey(e) { if (e.key === "Escape") onClose?.(); }
    if (open) window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [open, onClose]);
  if (!open) return null;
  return (
    <div className="drawer-overlay" onClick={onClose}>
      <aside className="drawer" onClick={(e) => e.stopPropagation()}>
        <header className="drawer-head">
          <div>
            <h3>{title}</h3>
            {subtitle && <p>{subtitle}</p>}
          </div>
          <button className="icon-btn" onClick={onClose} aria-label="Close">×</button>
        </header>
        <div className="drawer-body">{children}</div>
        {footer && <footer className="drawer-foot">{footer}</footer>}
      </aside>
    </div>
  );
}

export function Modal({ open, onClose, title, subtitle, children }) {
  useEffect(() => {
    function onKey(e) { if (e.key === "Escape") onClose?.(); }
    if (open) window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [open, onClose]);
  if (!open) return null;
  return (
    <div className="modal-overlay" onClick={onClose}>
      <div className="modal" onClick={(e) => e.stopPropagation()}>
        <header className="modal-head">
          <div>
            <h3>{title}</h3>
            {subtitle && <p>{subtitle}</p>}
          </div>
          <button className="icon-btn" onClick={onClose} aria-label="Close">×</button>
        </header>
        <div className="modal-body">{children}</div>
      </div>
    </div>
  );
}

export function Field({ label, children, wide }) {
  return (
    <div className={`field${wide ? " wide" : ""}`}>
      <span className="field-label">{label}</span>
      <div className="field-value">{children ?? "—"}</div>
    </div>
  );
}

export function Empty({ title, hint, action }) {
  return (
    <div className="empty-state">
      <div className="empty-mark">◇</div>
      <h3>{title}</h3>
      {hint && <p>{hint}</p>}
      {action}
    </div>
  );
}

/**
 * An analysis that cannot be performed for this gene, stated plainly.
 *
 * Deliberately not an `Empty`: no diamond mark, no warning or pending colour and no retry
 * prompt. A chicken MC1R protein encoded by one coding exon has no internal coding-exon
 * boundaries, so there is nothing to fix and nothing to wait for — showing that as an
 * error told readers their run had failed when it had simply finished.
 */
export function NotApplicable({ title, reason, badge, prerequisite, count }) {
  return (
    <div className="not-applicable" role="note">
      <div className="na-head">
        <h3>{title}</h3>
        {badge && <span className="na-badge">{badge}</span>}
      </div>
      {reason && <p className="na-reason">{reason}</p>}
      {prerequisite && (
        <p className="na-prereq muted small">
          <code>{prerequisite}</code>: {count ?? 0}
        </p>
      )}
    </div>
  );
}

/** `NotApplicable` for a settled state, otherwise the neutral empty state. */
export function AvailabilityState({ why, action }) {
  if (why?.notApplicable) {
    return (
      <NotApplicable title={why.title} reason={why.hint} badge={why.badge}
                     prerequisite={why.prerequisiteName} count={why.prerequisiteCount} />
    );
  }
  return <Empty title={why?.title} hint={why?.hint} action={action} />;
}

export function Spinner({ label }) {
  return <div className="spinner"><span className="spin" />{label || "Loading…"}</div>;
}

// Compact dropdown menu used to collapse dense toolbars (Tracks / Export) into a
// single button, so scientific pages show one focal plot and a restrained set of
// controls. Closes on outside click or Escape.
export function Menu({ label, title, align = "left", children, disabled }) {
  const [open, setOpen] = useState(false);
  const ref = useRef(null);
  useEffect(() => {
    function onDoc(e) { if (ref.current && !ref.current.contains(e.target)) setOpen(false); }
    function onKey(e) { if (e.key === "Escape") setOpen(false); }
    if (open) { document.addEventListener("mousedown", onDoc); document.addEventListener("keydown", onKey); }
    return () => { document.removeEventListener("mousedown", onDoc); document.removeEventListener("keydown", onKey); };
  }, [open]);
  return (
    <div className={`ui-menu${open ? " open" : ""}`} ref={ref}>
      <button type="button" className="seg-btn menu-btn" title={title} disabled={disabled}
        onClick={() => setOpen((v) => !v)} aria-expanded={open}>{label} <span className="menu-caret">▾</span></button>
      {open && <div className={`ui-menu-panel align-${align}`}
        onClick={(e) => e.stopPropagation()}>{children}</div>}
    </div>
  );
}

