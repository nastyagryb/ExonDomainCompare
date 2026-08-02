import { useState } from "react";

// Small calm badge shown for custom runs that reuse the curated human FGFR2
// IIIb/IIIc reference from the validated example dataset. Human is NOT counted
// as an analysed species unless the user explicitly selected homo_sapiens.
export default function HumanReferenceBadge({ humanReference }) {
  const [open, setOpen] = useState(false);
  const hr = humanReference || {};
  if (!hr.enabled && !hr.human_role) return null;
  const inPanel = hr.homo_sapiens_in_panel;
  const label = inPanel
    ? "Human: analysed + reference control"
    : "Human reference control enabled";
  return (
    <span className="human-ref-badge" onMouseLeave={() => setOpen(false)}>
      <button
        type="button"
        className="human-ref-chip"
        onClick={() => setOpen((v) => !v)}
        onMouseEnter={() => setOpen(true)}
      >
        <span className="hr-dot" aria-hidden>◆</span>
        {label}
      </button>
      {open && (
        <span className="human-ref-tip">
          This run uses the curated human FGFR2 IIIb/IIIc sequences from the validated
          example dataset as a fixed reference for marker comparison and IIIb/IIIc
          orientation.{" "}
          {inPanel
            ? "homo_sapiens is also part of the analysed panel."
            : "Human is not added to the analysed species panel."}
        </span>
      )}
    </span>
  );
}
