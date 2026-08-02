// Shared UI mapping for the exon-domain boundary consistency explorer (Module 1).
// Colours are color-blind-safe and match the static thesis figures.
// Semantics are "inspection", never "error / bad".

export const BOUNDARY_CLASS = {
  aligned_to_domain_boundary: {
    label: "Aligned",
    color: "#1B7837",
    tip: "Aligned to domain boundary (0–3 aa) — coincides with a protein-domain boundary.",
  },
  near_domain_boundary: {
    label: "Near boundary",
    color: "#A6DBA0",
    tip: "Near domain boundary (4–15 aa) — close to a protein-domain boundary.",
  },
  within_domain: {
    label: "Within domain",
    color: "#FDB863",
    tip: "Within domain (>15 aa from its edges) — inside a protein domain.",
  },
  between_domains: {
    label: "Between domains",
    color: "#B2ABD2",
    tip: "Between domains — in a linker region between two domains.",
  },
  review_or_missing: {
    label: "Missing / review",
    color: "#D9D9D9",
    tip: "No cassette / domain coordinate available for this protein.",
  },
};

export const BOUNDARY_ORDER = [
  "aligned_to_domain_boundary",
  "near_domain_boundary",
  "within_domain",
  "between_domains",
  "review_or_missing",
];

// exon-block display statuses that mark a low-confidence (inspection) display case
export const LOW_CONF_DISPLAY = new Set([
  "cassette_only_high_confidence",
  "native_exon_blocks_reconstructed",
]);

export function classInfo(cls) {
  return BOUNDARY_CLASS[cls] || BOUNDARY_CLASS.review_or_missing;
}

export function boundaryTypeLabel(bt) {
  if (bt === "cassette_start") return "Cassette start";
  if (bt === "cassette_end") return "Cassette end";
  return bt || "—";
}

export function isLowConfidence(exonBlockStatus) {
  return LOW_CONF_DISPLAY.has(exonBlockStatus || "");
}
