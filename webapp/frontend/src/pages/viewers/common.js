import { useEffect, useState } from "react";

import { forDataset, getActiveDataset, onActiveDatasetChange, payloadMatchesDataset }
  from "../../api.js";
import { prettyDomainName } from "./domainNames.js";

export { prettyDomainName };

/** The selected dataset, re-rendering the caller whenever the selection changes. */
export function useActiveDataset() {
  const [dataset, setDataset] = useState(getActiveDataset());
  useEffect(() => onActiveDatasetChange(setDataset), []);
  return dataset;
}

/**
 * Fetch one index for the selected dataset; expose `{data, error, loading, dataset}`.
 *
 * `loader` receives a client already bound to the dataset and the request's abort signal,
 * so it cannot accidentally read a dataset the user has since left:
 *
 *     useIndex((client) => client.msa())
 *
 * Three things make the result trustworthy after a dataset switch:
 *
 * 1. the effect re-runs when the dataset changes — it used to depend only on `preloaded`,
 *    so switching datasets left the previous dataset's data sitting in state;
 * 2. `data` is cleared before the new request starts, so no value from the previous
 *    dataset is ever shown under the new one's heading;
 * 3. a response is applied only if it names the dataset that is still selected, which
 *    discards a slow reply that lost the race against a faster switch.
 *
 * `preloaded` still short-circuits the fetch, letting the Gene Explorer hydrate viewers.
 */
export function useIndex(loader, preloaded) {
  const dataset = useActiveDataset();
  const [state, setState] = useState(() => ({
    data: preloaded || null, error: null, loading: !preloaded, dataset,
  }));

  useEffect(() => {
    if (preloaded) {
      // eslint-disable-next-line react-hooks/set-state-in-effect
      setState({ data: preloaded, error: null, loading: false, dataset });
      return undefined;
    }
    const controller = new AbortController();
    let alive = true;
    // Clear first: a stale value must never be visible next to a new dataset's identity.
    setState({ data: null, error: null, loading: true, dataset });
    const client = forDataset(dataset, controller.signal);
    Promise.resolve()
      .then(() => loader(client, controller.signal))
      .then((payload) => {
        if (!alive) return;
        if (!payloadMatchesDataset(payload, dataset)) return;
        setState({ data: payload, error: null, loading: false, dataset });
      })
      .catch((err) => {
        if (!alive || err?.name === "AbortError") return;
        setState({ data: null, error: err, loading: false, dataset });
      });
    return () => { alive = false; controller.abort(); };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [preloaded, dataset]);

  return state;
}

// The canonical availability states, written by scripts/shared_gene_analysis/
// run_availability.py and carried in each index's `availability` block.
export const AVAILABILITY = {
  available: "available",
  notApplicable: "not_applicable",
  scientificallyUnavailable: "scientifically_unavailable",
  technicallyMissing: "technically_missing",
  failed: "failed",
  pending: "pending",
  stale: "stale",
};

// Headline per state. Only `scientifically_unavailable` may speak about the biology;
// every other state is about the run, and saying "no cassette was found" for a file
// that was never written told users the species lacked a cassette it in fact has.
const HEADLINE = {
  [AVAILABILITY.notApplicable]: "Not applicable to this run",
  [AVAILABILITY.scientificallyUnavailable]: "No supported result for this run",
  [AVAILABILITY.technicallyMissing]: "Expected outputs were not generated",
  [AVAILABILITY.failed]: "This analysis step failed",
  [AVAILABILITY.pending]: "Analysis still running",
  [AVAILABILITY.stale]: "Results are out of date",
};

const FALLBACK_HINT = {
  [AVAILABILITY.technicallyMissing]: "Retry local analysis to rebuild them.",
  [AVAILABILITY.stale]: "Rebuild the run's indices to see the current analysis.",
  [AVAILABILITY.pending]: "This view appears as soon as the step finishes.",
};

/**
 * Title and hint for an unavailable view, taken from the index's own availability
 * block. `label` names the view for states that have no recorded reason.
 *
 * `not_applicable` is titled after the analysis itself — "Exon–domain boundary analysis
 * not applicable" — rather than with the generic "Not applicable to this run". The reader
 * needs to know *which* analysis does not apply and why, and a single-exon gene is a
 * property of the gene, not a shortcoming of the run.
 */
export function unavailableState(data, label) {
  const block = data?.availability || {};
  const state = block.state || AVAILABILITY.technicallyMissing;
  const notApplicable = state === AVAILABILITY.notApplicable;
  const name = block.label || label;
  return {
    state,
    notApplicable,
    // A recorded label names the state exactly ("Pending cluster annotation"); the generic
    // headline is only the fallback for an index that carries no availability block.
    title: notApplicable
      ? `${name} not applicable`
      : (block.label || HEADLINE[state] || `${label} not available`),
    hint: block.reason || FALLBACK_HINT[state]
      || `No ${label.toLowerCase()} data was found for this run.`,
    badge: block.badge || "",
    reasonCode: block.reason_code || "",
    prerequisiteName: block.prerequisite_name || "",
    prerequisiteCount: block.prerequisite_count ?? null,
    missingInputs: block.missing_inputs || [],
  };
}

/**
 * Whether an availability state is settled, i.e. the question was asked and answered.
 *
 * Used to decide that a navigation entry no longer needs a PENDING badge: an analysis that
 * resolved as not applicable is finished, not waiting.
 */
export function isResolvedState(state) {
  return state === AVAILABILITY.available
    || state === AVAILABILITY.notApplicable
    || state === AVAILABILITY.scientificallyUnavailable;
}

// Repeated InterPro entries are real separate domain instances (e.g. the three
// Ig-like domains of FGFR1) and must never be merged. Instance identity is
// established in the coordinate model (domain_instance_id = accession:start-end);
// this only mirrors it for older models that predate that field, and adds the
// group key used for "all instances of this entry" selection.
export function instanceIdOf(d) {
  return d?.domain_instance_id
    || `${d?.interpro_accession || "NA"}:${d?.start}-${d?.end}`;
}

export function domainInstances(domains) {
  const counts = new Map();
  for (const d of domains || []) {
    const k = d.interpro_accession || d.label;
    counts.set(k, (counts.get(k) || 0) + 1);
  }
  const seen = new Map();
  return (domains || []).map((d) => {
    const groupKey = d.interpro_accession || d.label;
    const many = (counts.get(groupKey) || 1) > 1;
    const n = (seen.get(groupKey) || 0) + 1;
    seen.set(groupKey, n);
    const base = prettyDomainName(d.label);
    const short = d.short_label || (many ? `${base} ${n}` : base);
    return {
      ...d,
      instanceId: instanceIdOf(d),
      instanceIndex: d.instance_number ?? (many ? n : null),
      groupKey,
      groupLabel: `All ${base}s`,
      shortLabel: short,
      instanceLabel: d.full_label || `${short} · aa ${d.start}–${d.end}`,
    };
  });
}

// Domain filter options built from the REAL domain instances of one model.
// "all" · one group option per repeated InterPro entry (the union of its
// instances) · one option per individual instance. Every option carries the
// explicit set of domain_instance_ids it selects, so no option can ever be
// resolved by InterPro accession alone.
export function domainFilterOptions(domains) {
  const insts = domainInstances(domains);
  const options = [{ value: "all", label: "All domains", kind: "all", instanceIds: [] }];
  const groups = new Map();
  for (const d of insts) {
    if (!groups.has(d.groupKey)) groups.set(d.groupKey, { label: d.groupLabel, ids: [] });
    groups.get(d.groupKey).ids.push(d.instanceId);
  }
  for (const [key, g] of groups) {
    if (g.ids.length < 2) continue;
    options.push({
      value: `grp:${key}`, label: `${g.label} (${g.ids.length})`,
      kind: "group", accession: key, instanceIds: [...g.ids],
    });
  }
  for (const d of insts) {
    options.push({
      value: `inst:${d.instanceId}`, label: d.instanceLabel, kind: "instance",
      accession: d.interpro_accession || null, instanceIds: [d.instanceId],
    });
  }
  return options;
}

export function boundaryInstanceId(b) {
  return b?.nearest_domain_instance_id ?? null;
}

// A boundary matches a domain filter only when the instance it was actually
// measured against is one of the instances the option selects.
export function matchesDomainFilter(b, value, options) {
  if (!value || value === "all") return true;
  const opt = (options || []).find((o) => o.value === value);
  if (!opt) return false;
  const iid = boundaryInstanceId(b);
  return iid != null && opt.instanceIds.includes(iid);
}

export const boundaryPosition = (b) => b?.protein_position
  ?? b?.boundary_position_aa ?? b?.start ?? null;

export const boundaryAbsDistance = (b) => b?.absolute_distance
  ?? b?.absolute_distance_aa
  ?? (b?.signed_distance != null ? Math.abs(b.signed_distance) : null);

// THE single boundary-filtering rule of the Boundary view. Every linked view —
// summary counts, architecture figure, signed-distance plot, selected-boundary
// list, evidence table and the exported visible TSV — reads the array this
// produces, so no two views can disagree.
export function filterBoundaries(boundaries, {
  domainFilter = "all", domainOptions = [], mappingFilter = "all", exonFilter = "all",
  distMin = "", distMax = "", candidate = null, candOnly = false, classFilter = null,
  sort = "position", classOf = (b) => b.boundary_class || b.category || b.class,
} = {}) {
  const kept = (boundaries || []).filter((b) => {
    if (!matchesDomainFilter(b, domainFilter, domainOptions)) return false;
    if (mappingFilter === "mapped" && b.mapping_status !== "mapped") return false;
    if (mappingFilter === "unmapped" && b.mapping_status === "mapped") return false;
    if (exonFilter !== "all" && (b.id || b.boundary_id) !== exonFilter) return false;
    const ad = boundaryAbsDistance(b);
    if (distMin !== "" && distMin != null && (ad == null || ad < Number(distMin))) return false;
    if (distMax !== "" && distMax != null && (ad == null || ad > Number(distMax))) return false;
    if (candOnly) {
      const p = boundaryPosition(b);
      if (!(candidate && p != null && p >= candidate.start && p <= candidate.end)) return false;
    }
    if (classFilter && classFilter.size && !classFilter.has(classOf(b))) return false;
    return true;
  });
  const key = sort === "distance"
    ? (b) => boundaryAbsDistance(b) ?? Number.POSITIVE_INFINITY
    : (b) => boundaryPosition(b) ?? Number.POSITIVE_INFINITY;
  return kept.sort((a, b) => key(a) - key(b) || (boundaryPosition(a) - boundaryPosition(b)));
}

// Neutral residue-agreement colour classes (no functional-effect claims).
export const RESIDUE_LABEL = {
  identical: "Identical to human",
  conservative: "Conservative substitution",
  nonconservative: "Non-conservative substitution",
  gap: "Gap / missing",
};

export const RESIDUE_LEGEND = [
  ["identical", "Identical"],
  ["conservative", "Conservative"],
  ["nonconservative", "Non-conservative"],
  ["gap", "Gap / missing"],
];

export function isoColorClass(iso) {
  const t = String(iso || "").toLowerCase();
  return t.includes("iiib") ? "iiib" : t.includes("iiic") ? "iiic" : "neutral";
}

export function pretty(name) {
  if (!name) return "—";
  const s = String(name).replaceAll("_", " ").trim();
  return s.charAt(0).toUpperCase() + s.slice(1);
}
