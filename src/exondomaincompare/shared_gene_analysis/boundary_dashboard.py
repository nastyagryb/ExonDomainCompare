"""boundary_dashboard.py — derived analysis dashboard over the shared coordinate model.

This module turns the validated protein-coordinate model (built by
``protein_coordinate_model.py``) into the authoritative data contract consumed by
the **global "Exon–Domain Boundaries" page**. It performs no coordinate maths and
never fabricates domains, boundaries or comparative results — it only *aggregates*
and *classifies views* over fields already present on the model:

  * page-mode resolution (single-species / multi-species / pending / unavailable)
  * single-species summary (class counts) + inspection cases + auto captions
  * the multi-species comparative data contract and comparable-boundary evidence priority

The comparative sections are intentionally *empty* unless real, mutually
comparable multi-species evidence exists (shared exon groups or an MSA-aligned
position). Boundaries are never compared only because both are called "E3→E4".

The frozen FGFR2 Boundary Consistency vocabulary is never touched here.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence
from .boundary_classification import CANONICAL_CLASSES, DEFAULT_NEAR_EDGE_THRESHOLD_AA, UNAVAILABLE, canonical_class

SCHEMA_VERSION = 1

# Page modes.
PAGE_VALIDATED_EVENT = "validated_event"
PAGE_SINGLE = "generic_single_species_results_ready"
PAGE_MULTI = "generic_multi_species_results_ready"
PAGE_PENDING = "pending_cluster"
PAGE_UNAVAILABLE = "unavailable"

# Comparable-boundary matching.
# Evidence priority. Exon rank is *secondary descriptive* evidence only and must
# never be the sole basis for calling two boundaries comparable.
MAPPING_METHOD_PRIORITY = [
    "shared_exon_group",
    "msa_aligned_position",
    "conserved_local_alignment",
    "domain_relative_context",
    "exon_rank",
]
COMPARABLE_STATES = (
    "high_confidence_comparable",
    "supported_comparable",
    "tentative",
    "unmapped",
    "unavailable",
)

# Descriptive inspection heuristic (NOT an error threshold): a coding-exon
# boundary this far from the nearest representative-domain edge is flagged for
# a closer look, not marked as wrong.
LARGE_DISTANCE_AA = 40

# Alignment columns this far apart are treated as the same junction shifted by a local
# indel, and the resulting group is marked tentative rather than supported. Two columns
# is deliberately tight: a wider window starts merging genuinely different junctions.
NEAR_COLUMN_TOLERANCE = 2


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _pos(b: Dict[str, Any]) -> Optional[int]:
    return b.get("protein_position") if b.get("protein_position") is not None \
        else b.get("boundary_position_aa")


def _abs(b: Dict[str, Any]) -> Optional[int]:
    if b.get("absolute_distance") is not None:
        return b.get("absolute_distance")
    if b.get("absolute_distance_aa") is not None:
        return b.get("absolute_distance_aa")
    if b.get("signed_distance") is not None:
        return abs(b["signed_distance"])
    return None


def _cls(b: Dict[str, Any]) -> str:
    return canonical_class(b.get("boundary_class") or b.get("category") or b.get("class"))


# --------------------------------------------------------------------------- #
# Page-mode resolution.
# --------------------------------------------------------------------------- #
def resolve_page_mode(index: Dict[str, Any], event_layer_type: Optional[str] = None) -> str:
    """Resolve the global-page mode purely from the coordinate model (no gene names)."""
    if (event_layer_type or "").lower() == "validated":
        return PAGE_VALIDATED_EVENT
    models = index.get("models") or []
    if not models:
        return PAGE_UNAVAILABLE
    available = [m for m in models if m.get("status") == "available"]
    if not available:
        return PAGE_PENDING
    if len(available) >= 2:
        return PAGE_MULTI
    return PAGE_SINGLE


# --------------------------------------------------------------------------- #
# single-species summary + inspection + captions (Parts 3, 4, 9, 20)
# --------------------------------------------------------------------------- #
def class_summary(model: Dict[str, Any]) -> Dict[str, int]:
    counts = {c: 0 for c in CANONICAL_CLASSES}
    bnds = model.get("exon_boundaries") or []
    for b in bnds:
        counts[_cls(b)] += 1
    return {"total": len(bnds), **counts}


def _status_badge(model: Dict[str, Any]) -> str:
    if model.get("status") != "available":
        return "pending_cluster"
    bnds = model.get("exon_boundaries") or []
    if not bnds:
        return "uncertain_mapping"
    unresolved = [b for b in bnds if (b.get("mapping_status") in ("unavailable", "unmapped")
                                      or _cls(b) == UNAVAILABLE)]
    if unresolved and len(unresolved) == len(bnds):
        return "uncertain_mapping"
    if unresolved:
        return "partial"
    return "results_ready"


def build_inspection_cases(model: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Return real applicable inspection cases without labelling them errors."""
    if model.get("status") != "available":
        return []
    bnds = model.get("exon_boundaries") or []
    candidates = model.get("candidate_regions") or []
    cases: List[Dict[str, Any]] = []

    def _add(kind: str, label: str, b: Dict[str, Any], detail: str):
        cases.append({
            "case_id": f"{kind}:{b.get('id')}",
            "kind": kind,
            "label": label,
            "boundary_id": b.get("id"),
            "boundary_label": b.get("label"),
            "protein_position": _pos(b),
            "detail": detail,
        })

    # 1. Large domain-edge distances (descriptive)
    classified = [b for b in bnds if b.get("signed_distance") is not None and _abs(b) is not None]
    for b in sorted(classified, key=lambda x: -(_abs(x) or 0)):
        if (_abs(b) or 0) >= LARGE_DISTANCE_AA:
            _add("large_domain_edge_distance", "Large domain-edge distance", b,
                 f"{_abs(b)} aa from the {b.get('nearest_domain_label') or 'nearest'} "
                 f"{b.get('nearest_edge_type') or 'domain'} edge.")

    # 2. Mapping requires inspection (mapped to a domain but flagged unmapped/uncertain)
    for b in bnds:
        if b.get("mapping_status") == "unmapped":
            _add("mapping_requires_inspection", "Mapping requires inspection", b,
                 "A representative domain exists but this boundary could not be mapped to one.")

    # 3. Representative annotation unavailable
    for b in bnds:
        if b.get("mapping_status") == "unavailable" or (_cls(b) == UNAVAILABLE
                                                        and b.get("mapping_status") != "unmapped"):
            _add("representative_annotation_unavailable", "Representative annotation unavailable", b,
                 "No representative InterPro domain is available to classify this boundary.")

    # 4. Incomplete coordinate evidence
    for b in bnds:
        if _pos(b) is None or (b.get("signed_distance") is not None and b.get("nearest_edge_position") is None):
            _add("incomplete_evidence", "Incomplete evidence", b,
                 "Coordinate evidence for this boundary is incomplete.")

    # 5. Candidate-overlapping boundaries
    for b in bnds:
        p = _pos(b)
        if p is None:
            continue
        hit = next((c for c in candidates if c.get("start") is not None
                    and c.get("start") <= p <= c.get("end")), None)
        if hit:
            _add("candidate_associated", "Candidate-associated boundary", b,
                 f"Falls within exploratory candidate {hit.get('id')} "
                 f"(aa {hit.get('start')}–{hit.get('end')}).")

    return cases


def generate_caption(model: Dict[str, Any], gene_symbol: Optional[str] = None,
                     *, species_scope: str = "single_species") -> Dict[str, Any]:
    """Return an editable, cautious automatic caption."""
    gene = (gene_symbol or model.get("gene_symbol") or "gene")
    sci = model.get("scientific_name") or model.get("species_id") or "the species"
    pid = model.get("protein_id") or "the primary protein"
    thr = model.get("near_edge_threshold_aa", DEFAULT_NEAR_EDGE_THRESHOLD_AA)
    coord = model.get("coordinate_system") or "protein_1_based_inclusive"
    available = model.get("status") == "available"
    stage = "results ready" if available else "pending cluster"
    if available:
        text = (
            f"Exon–domain boundary analysis of {sci} {gene}. Internal coding-exon "
            f"boundaries were projected onto amino-acid coordinates of {pid} and classified "
            f"relative to representative InterPro domains using a near-edge threshold of "
            f"{thr} amino acids ({coord}). Classifications describe positional overlap only "
            f"and do not by themselves imply a functional effect."
        )
    else:
        text = (
            f"Exon–domain boundary preparation for {sci} {gene}. Coding-exon boundary "
            f"positions on {pid} are available; exon–domain classifications relative to "
            f"representative InterPro domains ({coord}) will be calculated after the cluster "
            f"InterProScan results are returned."
        )
    return {
        "text": text,
        "fields": {
            "gene": gene,
            "species_scope": species_scope,
            "selected_primary_protein": pid,
            "domain_annotation_source": "representative InterPro domains",
            "near_edge_threshold_aa": thr,
            "coordinate_system": coord,
            "analysis_status": model.get("status"),
            "analysis_stage": stage,
        },
    }


def build_single_species_dashboard(model: Dict[str, Any],
                                   gene_symbol: Optional[str] = None) -> Dict[str, Any]:
    exons = model.get("exons") or []
    bnds = model.get("exon_boundaries") or []
    domains = model.get("representative_domains") or []
    return {
        "header": {
            "gene": gene_symbol or model.get("gene_symbol"),
            "scientific_name": model.get("scientific_name"),
            "species_id": model.get("species_id"),
            "protein_id": model.get("protein_id"),
            "transcript_id": model.get("transcript_id"),
            "protein_length": model.get("protein_length"),
            "n_coding_exons": len(exons),
            "n_internal_boundaries": len(bnds),
            "representative_domain_source": ("representative InterPro domains"
                                             if domains else "pending post-cluster InterProScan"),
            "near_edge_threshold_aa": model.get("near_edge_threshold_aa",
                                                DEFAULT_NEAR_EDGE_THRESHOLD_AA),
            "analysis_stage": ("results ready" if model.get("status") == "available"
                               else "pending cluster"),
            "status_badge": _status_badge(model),
            "run_id": (model.get("provenance") or {}).get("run_id"),
        },
        "summary": class_summary(model),
        "inspection_cases": build_inspection_cases(model),
        "caption": generate_caption(model, gene_symbol),
    }


# --------------------------------------------------------------------------- #
# multi-species contract + comparable-boundary matching (Parts 12, 13, 16)
# --------------------------------------------------------------------------- #
def _exon_shared_groups(model: Dict[str, Any]) -> Dict[str, Optional[str]]:
    out: Dict[str, Optional[str]] = {}
    for e in model.get("exons") or []:
        out[e.get("id")] = (e.get("tooltip") or {}).get("shared_exon_group")
    return out


def _observation(model: Dict[str, Any], boundary: Dict[str, Any], method: str,
                 mapping_status: str, coverage: float) -> Dict[str, Any]:
    """One species' real observation of a comparable boundary.

    Every field is copied from the species' own boundary record; nothing is derived
    across species here. This is the single record the matrix cell, the hover, the
    detail row and the export table all read, so a cell can never show a value the
    detail view contradicts.

    ``domain_annotation_available`` is false when the species has no representative
    domain near this boundary — a real annotation gap that must stay distinguishable
    from a distance of zero.
    """
    has_domain = boundary.get("nearest_domain_instance_id") is not None
    return {
        "species_id": model.get("species_id"),
        "scientific_name": model.get("scientific_name"),
        "taxonomic_group": (model.get("provenance") or {}).get("taxonomic_group"),
        "protein_id": model.get("protein_id"),
        "transcript_id": model.get("transcript_id"),
        "boundary_id": boundary.get("id"),
        "exon_transition": (boundary.get("label")
                            or (boundary.get("tooltip") or {}).get("adjacent_exon_transition")),
        # The source exon identifiers, kept beside the display labels: a label is a
        # rendering choice ("E2→E3") while the id is what the annotation called the exon,
        # and only the latter lets a reader check the group against the source.
        "left_exon_id": boundary.get("left_exon_id"),
        "right_exon_id": boundary.get("right_exon_id"),
        "left_exon_label": boundary.get("left_exon_label"),
        "right_exon_label": boundary.get("right_exon_label"),
        "native_position": _pos(boundary),
        "msa_column": boundary.get("msa_column"),
        "msa_mapping_status": boundary.get("msa_mapping_status"),
        "nearest_domain_instance_id": boundary.get("nearest_domain_instance_id"),
        "nearest_domain_accession": boundary.get("nearest_domain_accession"),
        "nearest_domain_label": boundary.get("nearest_domain_short_label"),
        "nearest_domain_full_label": boundary.get("nearest_domain_full_label"),
        "nearest_domain_start": boundary.get("nearest_domain_start"),
        "nearest_domain_end": boundary.get("nearest_domain_end"),
        "nearest_edge": boundary.get("nearest_edge") or boundary.get("nearest_edge_type"),
        "nearest_edge_position": boundary.get("nearest_edge_position"),
        "signed_distance": boundary.get("signed_distance"),
        "absolute_distance": _abs(boundary),
        "boundary_class": _cls(boundary),
        "domain_annotation_available": has_domain,
        "mapping_method": method,
        "mapping_status": mapping_status,
        # The qualitative evidence status is the mapping confidence. The fraction below
        # is species coverage and is kept under its old key for compatibility: a group
        # can be mapped in every species and still be tentative, so presenting the
        # fraction as "confidence" would claim certainty the evidence does not carry.
        "mapping_confidence": round(coverage, 3),
        "species_coverage": round(coverage, 3),
    }


def match_comparable_boundaries(models: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Group boundaries across species using the evidence priority.

    Returns comparable-boundary groups ONLY when ≥2 available species contribute
    to a group through strong evidence (a shared exon group, or an MSA-aligned
    protein position). Never groups by exon rank alone. Returns ``[]`` for a
    single available species or when no shared/aligned evidence exists.
    """
    available = [m for m in models if m.get("status") == "available"]
    if len(available) < 2:
        return []

    # A shared exon group is a hash of the *genomic* interval, so it identifies the same
    # exon across isoforms of one species and can never coincide between species sitting
    # on different assemblies. Using it as the first-choice cross-species key put each
    # species in its own bucket and silently discarded every group. It is therefore
    # accepted here only for the values that demonstrably occur in more than one
    # species; anything else falls through to the alignment column.
    group_species: Dict[str, set] = {}
    for m in available:
        for g in _exon_shared_groups(m).values():
            if g:
                group_species.setdefault(g, set()).add(m.get("species_id"))
    cross_species_groups = {g for g, sp in group_species.items() if len(sp) > 1}

    buckets: Dict[tuple, Dict[str, Any]] = {}
    for m in available:
        groups = _exon_shared_groups(m)
        for b in m.get("exon_boundaries") or []:
            lg = groups.get(b.get("left_exon_id"))
            rg = groups.get(b.get("right_exon_id"))
            msa_col = b.get("msa_column")
            if lg and rg and lg in cross_species_groups and rg in cross_species_groups:
                key = ("shared_exon_group", lg, rg)
                method = "shared_exon_group"
                evidence = {"shared_exon_group": [lg, rg]}
            elif msa_col is not None:
                key = ("msa_aligned_position", msa_col)
                method = "msa_aligned_position"
                evidence = {"msa_column": msa_col}
            else:
                # No strong comparative evidence — do NOT compare by exon rank alone.
                continue
            bucket = buckets.setdefault(key, {"method": method, "members": [], "evidence": evidence})
            bucket["members"].append((m, b))

    # Boundaries that agree to within a residue or two of alignment column describe the
    # same junction shifted by a local indel, and dropping them would hide real
    # observations. They are merged into one bucket but never promoted to
    # "supported": the group is reported as tentative so a reader can see that the
    # columns were close rather than identical.
    single = [k for k, v in buckets.items()
              if k[0] == "msa_aligned_position"
              and len({m.get("species_id") for m, _ in v["members"]}) < 2]
    for a in sorted(single, key=lambda k: k[1]):
        if a not in buckets:
            continue
        a_species = {m.get("species_id") for m, _ in buckets[a]["members"]}
        for b_key in sorted(single, key=lambda k: k[1]):
            if b_key == a or b_key not in buckets or a not in buckets:
                continue
            if abs(int(a[1]) - int(b_key[1])) > NEAR_COLUMN_TOLERANCE:
                continue
            if {m.get("species_id") for m, _ in buckets[b_key]["members"]} & a_species:
                continue
            buckets[a]["members"].extend(buckets[b_key]["members"])
            buckets[a]["method"] = "msa_aligned_position"
            buckets[a]["tentative"] = True
            buckets[a]["evidence"] = {
                "msa_column": int(a[1]),
                "msa_columns_merged": sorted({int(a[1]), int(b_key[1])}),
                "column_offset": abs(int(a[1]) - int(b_key[1])),
            }
            del buckets[b_key]
            a_species = {m.get("species_id") for m, _ in buckets[a]["members"]}

    result: List[Dict[str, Any]] = []
    for i, (key, bucket) in enumerate(buckets.items()):
        species_ids = {m.get("species_id") for m, _ in bucket["members"]}
        if len(species_ids) < 2:
            continue  # present in a single species → not comparable
        method = bucket["method"]
        coverage = len(species_ids) / len(available)
        if bucket.get("tentative"):
            state = "tentative"
        elif method == "shared_exon_group" and coverage >= 1.0:
            state = "high_confidence_comparable"
        elif method == "shared_exon_group":
            state = "supported_comparable"
        elif method == "msa_aligned_position":
            state = "supported_comparable"
        else:
            state = "tentative"
        per_species = [_observation(m, b, method, state, coverage)
                       for m, b in bucket["members"]]
        result.append({
            "comparable_boundary_group_id": f"CBG{i + 1}",
            "label": f"Comparable boundary group {i + 1}",
            "mapping_method": method,
            "mapping_status": state,
            "confidence": round(coverage, 3),
            # How many of the dataset's species this group was actually observed in, and
            # which ones. Exported tables read these; without them a group table states a
            # comparison without naming its participants.
            "n_species": len(species_ids),
            "species_coverage": ",".join(sorted(s for s in species_ids if s)),
            "msa_column": bucket["evidence"].get("msa_column"),
            "shared_exon_group": bucket["evidence"].get("shared_exon_group"),
            "domain_relative_context": None,
            "per_species_native_positions": per_species,
            "supporting_evidence": bucket["evidence"],
        })
    return result


def boundary_position_consistency(models: Sequence[Dict[str, Any]],
                                  comparable: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Calculate comparable-boundary-group consistency statistics.

    Uses cautious wording (boundary-position consistency / domain-edge proximity
    consistency); never asserts "conserved" without meeting explicit criteria.
    Empty unless real comparable groups exist.
    """
    stats: List[Dict[str, Any]] = []
    n_available = len([m for m in models if m.get("status") == "available"]) or 1
    for g in comparable:
        members = g.get("per_species_native_positions") or []
        signed = [m["signed_distance"] for m in members if m.get("signed_distance") is not None]
        classes = [m["boundary_class"] for m in members if m.get("boundary_class")]
        mapped = len(members)
        exact_near = sum(1 for c in classes if c in ("exact_domain_edge", "near_domain_edge"))
        median = None
        spread = None
        if signed:
            s = sorted(signed)
            mid = len(s) // 2
            median = s[mid] if len(s) % 2 else (s[mid - 1] + s[mid]) / 2
            spread = max(s) - min(s)
        # A set's iteration order is process-randomized. Resolve equal counts
        # deterministically so two scientifically identical rebuilds cannot emit
        # different displayed summaries. Lexical max preserves the accepted
        # projection for existing ties while leaving every count untouched.
        dominant = max(
            sorted(set(classes)), key=lambda value: (classes.count(value), value)
        ) if classes else None

        # With two species a median is just the mean of the only two observations and a
        # "distribution" is a straight line, so the honest summary is the raw pair and
        # the difference between them. Median and spread stay available for larger sets
        # but are marked as not the primary statistic here.
        raw = [{"species_id": m.get("species_id"),
                "signed_distance": m.get("signed_distance"),
                "boundary_class": m.get("boundary_class"),
                "domain_annotation_available": m.get("domain_annotation_available")}
               for m in members]
        pair_difference = None
        if len(signed) == 2:
            pair_difference = abs(signed[0] - signed[1])
        elif len(signed) > 2:
            pair_difference = max(signed) - min(signed)

        confidences = [m.get("mapping_confidence") for m in members
                       if m.get("mapping_confidence") is not None]
        stats.append({
            "comparable_boundary_group_id": g.get("comparable_boundary_group_id"),
            "species_coverage": round(mapped / n_available, 3),
            "species_with_mapped_boundary": mapped,
            "n_species_available": n_available,
            "mapping_coverage": round(mapped / n_available, 3),
            "exact_or_near_proportion": round(exact_near / mapped, 3) if mapped else 0,
            "raw_signed_distances": raw,
            "cross_species_difference": pair_difference,
            "distance_range": [min(signed), max(signed)] if signed else None,
            "median_signed_distance": median,
            "distance_spread": spread,
            "primary_statistic": "raw_pair" if len(signed) <= 2 else "distribution",
            "dominant_class": dominant,
            "classes_differ": len(set(classes)) > 1,
            "domain_annotation_available_in_all": all(
                m.get("domain_annotation_available") for m in members) if members else False,
            "mapping_confidence": g.get("confidence"),
            "mapping_confidence_distribution": (
                {"min": min(confidences), "max": max(confidences)} if confidences else None),
            "mapping_status": g.get("mapping_status"),
            "metric_label": "Boundary-position consistency",
        })
    return stats


def build_boundary_matrix(models: Sequence[Dict[str, Any]],
                          comparable: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Build species-by-comparable-boundary-group cell states.

    A cell carries the species' own observation verbatim, so the matrix, its hover and
    the detail panel are the same numbers rather than three re-derivations. Cells with
    no observation carry an explicit state and null values — never a zero distance,
    which would read as a boundary sitting exactly on a domain edge.
    """
    matrix: List[Dict[str, Any]] = []
    if not comparable:
        return matrix
    for m in models:
        row = {"species_id": m.get("species_id"),
               "scientific_name": m.get("scientific_name"),
               "taxonomic_group": (m.get("provenance") or {}).get("taxonomic_group"),
               "protein_id": m.get("protein_id"),
               "analysis_status": m.get("status"),
               "cells": []}
        for g in comparable:
            member = next((x for x in (g.get("per_species_native_positions") or [])
                           if x.get("species_id") == m.get("species_id")), None)
            if m.get("status") != "available":
                state = "result_pending"
            elif member is None:
                state = "boundary_absent_or_unmapped"
            else:
                state = member.get("boundary_class") or "unavailable_or_uncertain"
            cell = {
                "comparable_boundary_group_id": g.get("comparable_boundary_group_id"),
                "state": state,
                "observed": member is not None,
                "native_position": member.get("native_position") if member else None,
                "signed_distance": member.get("signed_distance") if member else None,
                "absolute_distance": member.get("absolute_distance") if member else None,
                "mapping_status": (member or {}).get("mapping_status") or g.get("mapping_status"),
                "mapping_method": (member or {}).get("mapping_method") or g.get("mapping_method"),
                "mapping_confidence": (member or {}).get("mapping_confidence"),
                "observation": member,
            }
            row["cells"].append(cell)
        matrix.append(row)
    return matrix


def build_comparative_inspection_cases(models: Sequence[Dict[str, Any]],
                                       comparable: Sequence[Dict[str, Any]],
                                       stats: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Groups worth a second look, each pointing at a real comparable-boundary group.

    Every case is phrased as an observation, not a verdict: a large cross-species
    distance difference or a missing domain annotation may be biology, an annotation
    gap, or an alignment artefact, and this layer cannot tell which. Wording therefore
    never calls a discrepancy an error.
    """
    cases: List[Dict[str, Any]] = []
    by_id = {s.get("comparable_boundary_group_id"): s for s in stats}
    available_species = {m.get("species_id") for m in models if m.get("status") == "available"}

    for g in comparable:
        gid = g.get("comparable_boundary_group_id")
        st = by_id.get(gid, {})
        members = g.get("per_species_native_positions") or []
        member_species = {m.get("species_id") for m in members}

        def add(case_type: str, label: str, detail: str, severity: str = "review") -> None:
            cases.append({
                "case_id": f"{case_type}:{gid}",
                "case_type": case_type,
                "comparable_boundary_group_id": gid,
                "label": label,
                "detail": detail,
                "severity": severity,
                "species_ids": sorted(member_species),
                "msa_column": g.get("msa_column"),
            })

        diff = st.get("cross_species_difference")
        if diff is not None and diff >= LARGE_DISTANCE_AA:
            add("large_distance_difference",
                "Large cross-species distance difference",
                f"The signed distances to the nearest representative-domain edge differ "
                f"by {diff} aa between species. This may reflect a real difference in "
                f"exon structure or domain annotation extent.")

        if st.get("classes_differ"):
            classes = sorted({m.get("boundary_class") for m in members if m.get("boundary_class")})
            add("different_boundary_classes",
                "Different boundary classes",
                f"The same aligned position is classified as {' vs '.join(classes)} in "
                f"different species.")

        missing_domain = [m.get("species_id") for m in members
                          if not m.get("domain_annotation_available")]
        if missing_domain:
            add("missing_domain_annotation",
                "Domain annotation missing in one species",
                f"No representative domain was annotated near this boundary in "
                f"{', '.join(missing_domain)}, so no distance can be measured there. "
                f"Absence of annotation is not evidence of absence of the domain.")

        if g.get("mapping_status") == "tentative":
            ev = g.get("supporting_evidence") or {}
            offset = ev.get("column_offset")
            add("tentative_mapping",
                "Tentative alignment mapping",
                f"Boundaries were merged across an alignment-column offset of "
                f"{offset} column(s); the positions are close but not identical, so "
                f"equivalence is not established.",
                severity="caution")

        if member_species and member_species != available_species:
            only = sorted(available_species - member_species)
            add("single_species_observation",
                "Observation missing in at least one species",
                f"No comparable boundary was mapped at this position in "
                f"{', '.join(only)}.")

        instances = {m.get("nearest_domain_accession") for m in members
                     if m.get("nearest_domain_accession")}
        if len(instances) > 1:
            add("inconsistent_domain_instance",
                "Boundary assigned to different domain families",
                f"The nearest representative domain belongs to different InterPro "
                f"entries across species ({', '.join(sorted(instances))}).")

    # Candidate-associated boundaries are a species-level property, so they are matched
    # onto the comparative groups by boundary id rather than recomputed.
    candidate_boundaries = set()
    for m in models:
        for region in m.get("candidate_regions") or []:
            for bid in region.get("associated_boundary_ids") or []:
                candidate_boundaries.add(bid)
    if candidate_boundaries:
        for g in comparable:
            hit = [m for m in (g.get("per_species_native_positions") or [])
                   if m.get("boundary_id") in candidate_boundaries]
            if hit:
                gid = g.get("comparable_boundary_group_id")
                cases.append({
                    "case_id": f"candidate_associated:{gid}",
                    "case_type": "candidate_associated",
                    "comparable_boundary_group_id": gid,
                    "label": "Candidate-associated boundary",
                    "detail": ("At least one species' boundary in this group lies within a "
                               "candidate region flagged by the exploratory event layer."),
                    "severity": "review",
                    "species_ids": sorted({m.get("species_id") for m in hit}),
                    "msa_column": g.get("msa_column"),
                })
    return cases


def build_multi_species_contract(index: Dict[str, Any]) -> Dict[str, Any]:
    """Build the comparative data structure. Comparative arrays fill only
    when real comparable evidence exists; otherwise they stay empty (honest)."""
    models = index.get("models") or []
    species_rows = [{
        "species_id": m.get("species_id"),
        "scientific_name": m.get("scientific_name"),
        "taxonomic_group": (m.get("provenance") or {}).get("taxonomic_group"),
        "primary_protein": m.get("protein_id"),
        "transcript": m.get("transcript_id"),
        "protein_length": m.get("protein_length"),
        "analysis_status": m.get("status"),
    } for m in models]
    comparable = match_comparable_boundaries(models)
    stats = boundary_position_consistency(models, comparable)
    cases = build_comparative_inspection_cases(models, comparable, stats)

    # Filter vocabularies are published with the data so the frontend offers exactly the
    # values that occur, instead of a hardcoded list that can drift out of step.
    domain_groups: Dict[str, Dict[str, Any]] = {}
    for g in comparable:
        for obs in g.get("per_species_native_positions") or []:
            acc = obs.get("nearest_domain_accession")
            if acc:
                domain_groups.setdefault(acc, {
                    "interpro_accession": acc,
                    "label": obs.get("nearest_domain_label") or acc,
                })
    return {
        "available": bool(comparable),
        "n_species": len([m for m in models if m.get("status") == "available"]),
        "species_rows": species_rows,
        "comparable_boundary_groups": comparable,
        "mapping_methods": MAPPING_METHOD_PRIORITY,
        "boundary_matrix": build_boundary_matrix(models, comparable),
        "distance_statistics": stats,
        "inspection_cases": cases,
        "filter_options": {
            "species": [{"species_id": r["species_id"],
                         "scientific_name": r["scientific_name"],
                         "taxonomic_group": r["taxonomic_group"]} for r in species_rows],
            "taxonomic_groups": sorted({r["taxonomic_group"] for r in species_rows
                                        if r.get("taxonomic_group")}),
            "boundary_classes": sorted({obs.get("boundary_class")
                                        for g in comparable
                                        for obs in g.get("per_species_native_positions") or []
                                        if obs.get("boundary_class")}),
            "representative_domain_groups": sorted(domain_groups.values(),
                                                   key=lambda d: d["label"]),
            "mapping_statuses": sorted({g.get("mapping_status") for g in comparable
                                        if g.get("mapping_status")}),
            "edges": sorted({obs.get("nearest_edge")
                             for g in comparable
                             for obs in g.get("per_species_native_positions") or []
                             if obs.get("nearest_edge")}),
            "inspection_case_types": sorted({c["case_type"] for c in cases}),
        },
        "near_edge_band_aa": [-5, 5],
    }


# --------------------------------------------------------------------------- #
# top-level builder (injected into the coordinate-model index)
# --------------------------------------------------------------------------- #
def build_boundary_dashboard(index: Dict[str, Any],
                             event_layer_type: Optional[str] = None) -> Dict[str, Any]:
    mode = resolve_page_mode(index, event_layer_type)
    models = index.get("models") or []
    primary = next((m for m in models if m.get("status") == "available"),
                   models[0] if models else None)
    gene = index.get("gene_symbol")
    out: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "page_mode": mode,
        "gene_symbol": gene,
        "run_id": index.get("run_id"),
        "coordinate_system": index.get("coordinate_system"),
        "near_edge_threshold_aa": (primary or {}).get("near_edge_threshold_aa",
                                                      DEFAULT_NEAR_EDGE_THRESHOLD_AA),
        "species_available": [m.get("species_id") for m in models
                              if m.get("status") == "available"],
        "generated_by": "src/exondomaincompare/shared_gene_analysis/boundary_dashboard.py",
    }
    if mode in (PAGE_SINGLE, PAGE_MULTI, PAGE_PENDING) and primary is not None:
        out["single_species"] = build_single_species_dashboard(primary, gene)
    else:
        out["single_species"] = None
    out["multi_species"] = build_multi_species_contract(index)
    return out
