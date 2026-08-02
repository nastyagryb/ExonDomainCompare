/**
 * The display vocabulary the UI speaks: how a raw pipeline or QC token becomes a
 * calm, colour-blind-safe status, and how a species name is written out.
 *
 * These are pure functions and constants, not components. They live beside `ui.jsx`
 * rather than inside it because a module that exports both components and plain
 * values cannot be hot-reloaded as a unit, and because the status vocabulary is
 * consumed by viewers that do not render a badge themselves.
 */

// Map any pipeline status token to a calm, color-blind-safe class.
export function statusClass(value) {
  const v = String(value || "").toLowerCase();
  if (!v) return "unknown";
  if (v.includes("excluded") || v === "fail" || v.includes("_fail")) return "excluded";
  if (v.includes("supplement") || v.includes("review") || v.includes("unresolved") || v.includes("manual"))
    return "review";
  if (v.includes("with_minor_flags") || v === "minor" || v.includes("moderate") || v.includes("warn"))
    return "minor";
  if (v === "pass" || v === "strong" || v === "robust" || v.includes("rescued_ok") ||
      v.includes("ready") || v.includes("primary") || v.includes("supported") ||
      v.includes("high_confidence") || v.includes("confirmed") || v === "true")
    return "accepted";
  return "neutral";
}

const READINESS_LABEL = {
  accepted: "Primary",
  minor: "Primary · minor flags",
  review: "Review",
  excluded: "Excluded",
  neutral: "—",
  unknown: "—",
};

export function readinessLabel(cls) {
  return READINESS_LABEL[cls] || cls;
}

// --------------------------------------------------------------------------- //
// Main display status normalization.
//
// Scientific principle: the UI shows the calm biological/architecture status.
// Pipeline bookkeeping (upstream correction, rescue validity, native offset,
// fallback / sanitized coordinates, display-coordinate confidence) is never a
// display status — it stays in the downloadable tables and internal JSON.
// This helper maps any raw pipeline / QC token to one of:
//   Accepted · Supported · Supported with minor note · Inspection note ·
//   Supplement / review · Failed · Not available
// with a calm, color-blind-safe badge class.
const MAIN_STATUS_RULES = [
  // true biological / coordinate gate failure (red) — narrow on purpose
  [/(^fail$|_fail\b|failed|gate_failed|excluded)/, { label: "Failed", cls: "excluded" }],
  // supplement / review membership (amber) — a dataset decision, not a failure
  [/(supplement|review_included|review_only|review_supplement|^review$|is_review|unresolved|manual_review)/,
    { label: "Supplement / review", cls: "review" }],
  // display-coordinate confidence only (soft amber outline) — inspection, not error
  [/(cassette_only|hidden_untrusted|untrusted|unusual_domain_order|low.?confidence|inspection)/,
    { label: "Inspection note", cls: "inspection" }],
  // accepted but with a minor note (calm blue)
  [/(with_minor_flags|minor_length_clamped|native_exon_blocks_reconstructed|reconstructed|artifact_resolved|display.*resolved|clamp|minor)/,
    { label: "Supported with minor note", cls: "minor" }],
  // final accepted / supported (green) — incl. provenance-only corrections/rescues
  [/(architecture_supported|validated|accepted|primary_ready|primary|supported|rescued|reconciled|corrected|high_confidence|confirmed|robust|strong|ready|^pass$|^true$)/,
    { label: "Accepted", cls: "accepted" }],
];

export function mainDisplayStatus(raw) {
  const v = String(raw || "").toLowerCase().trim();
  if (!v || v === "neutral" || v === "unknown" || v === "not_applicable" ||
      v === "na" || v === "none" || v === "—") {
    return { label: "Not available", cls: "neutral" };
  }
  for (const [re, out] of MAIN_STATUS_RULES) if (re.test(v)) return out;
  return { label: "Not available", cls: "neutral" };
}

// architecture_supported keeps the calmer "Supported" wording (green), while
// generic "accepted/validated" rows read "Accepted".
export function architectureStatusLabel(raw) {
  const v = String(raw || "").toLowerCase();
  const base = mainDisplayStatus(raw);
  if (base.cls === "accepted" && v.includes("architecture_supported")) {
    return { label: "Supported", cls: "accepted" };
  }
  return base;
}

// Provenance / note refinements (a native offset or an upstream correction is a
// framework success, not biological uncertainty).
export const TONE_GLYPH = { corrected: "✎", rescued: "✚", offset: "≈" };
export const TONE_LABEL = {
  corrected: "corrected by sequence evidence",
  rescued: "rescued & validated",
  offset: "native offset (note)",
};

export function prettySpecies(name) {
  if (!name) return "—";
  const s = String(name).replaceAll("_", " ").trim();
  return s.charAt(0).toUpperCase() + s.slice(1);
}
