"""boundary_classification.py — canonical exon-boundary vs domain classifier.

Single source of truth for the **generic** (gene-agnostic) exon-boundary
classification used by the coordinate model, the Exon–Domain Boundaries page and
the signed-distance plots. It uses *representative / integrated InterPro domains*
only (never raw member signatures), and a mutually-exclusive priority.

This module does not touch the frozen FGFR2 Boundary Consistency vocabulary
(`aligned_to_domain_boundary` / `near_domain_boundary` / …), which stays in
`scripts/analyze_exon_domain_boundary_consistency.py`.

Canonical classes (mutually exclusive, in priority order):
    exact_domain_edge          absolute distance == 0
    near_domain_edge           0 < absolute distance <= threshold (default 5 aa)
    inside_domain              boundary strictly inside a representative domain
    outside_annotated_domains  representative domains exist, boundary outside all
    unavailable_or_uncertain   no representative domain / incomplete coordinates

Sign convention (single definition for the whole project):
    ``signed_distance = boundary_position - nearest_edge_position``
so a boundary N-terminal of a domain start is negative, a boundary C-terminal of a
domain end is positive, and a boundary inside a domain measured against that
domain's end is negative.

Repeated InterPro entries are distinct feature *instances*. Every domain passed in
is identified by ``<accession>:<start>-<end>`` and the classifier always reports the
instance it actually measured — an accession alone can never identify a feature.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Sequence

DEFAULT_NEAR_EDGE_THRESHOLD_AA = 5

EXACT = "exact_domain_edge"
NEAR = "near_domain_edge"
INSIDE = "inside_domain"
OUTSIDE = "outside_annotated_domains"
UNAVAILABLE = "unavailable_or_uncertain"

CANONICAL_CLASSES = (EXACT, NEAR, INSIDE, OUTSIDE, UNAVAILABLE)

# Map from the legacy Core-runner category names to the canonical vocabulary so
# existing indices keep working while the canonical class is added alongside.
LEGACY_TO_CANONICAL = {
    "exact_edge": EXACT,
    "near_edge": NEAR,
    "inside_domain": INSIDE,
    "outside_domain": OUTSIDE,
    "unknown": UNAVAILABLE,
    "": UNAVAILABLE,
}



def canonical_class(legacy_or_canonical: Optional[str]) -> str:
    """Normalize any incoming class string to the canonical vocabulary."""
    v = (legacy_or_canonical or "").strip()
    if v in CANONICAL_CLASSES:
        return v
    return LEGACY_TO_CANONICAL.get(v, UNAVAILABLE)


def domain_instance_id(accession: Optional[str], start: Any, end: Any) -> str:
    """Stable identity of ONE domain feature instance: ``<accession>:<start>-<end>``.

    Two instances of the same InterPro entry (e.g. FGFR1 ``IPR007110`` at 33–118 and
    at 145–244) get different ids, so no lookup can ever collapse them.
    """
    return f"{accession or 'NA'}:{start}-{end}"


def instance_id_of(domain: Dict[str, Any]) -> str:
    return domain.get("domain_instance_id") or domain_instance_id(
        domain.get("interpro_accession"), domain.get("start"), domain.get("end"))


def _nearest_edge(pos: int, domains: Sequence[Dict[str, Any]]):
    """Return (signed_distance, abs_distance, domain, edge_type) for nearest edge.

    ``signed_distance = pos - edge``: negative means the boundary is N-terminal
    (left) of that edge. Domains are scanned in start order and ties keep the
    earlier-starting instance, so the result is deterministic per *instance* — the
    winning entry is a concrete domain object, never an accession.
    """
    best = None
    ordered = sorted(
        (d for d in domains if d.get("start") is not None and d.get("end") is not None),
        key=lambda d: (d["start"], d["end"], str(d.get("interpro_accession") or "")))
    for d in ordered:
        for edge_type, edge in (("start", d["start"]), ("end", d["end"])):
            signed = pos - edge
            dist = abs(signed)
            if best is None or dist < best[1]:
                best = (signed, dist, d, edge_type)
    return best


def classify_boundary(
    boundary_position: int,
    domains: Sequence[Dict[str, Any]],
    *,
    threshold: int = DEFAULT_NEAR_EDGE_THRESHOLD_AA,
) -> Dict[str, Any]:
    """Classify one exon boundary against representative domains.

    ``domains``: sequence of dicts with ``start``/``end`` (1-based inclusive) and
    optional ``id``/``label``. Returns the stored boundary record fields.
    """
    rec: Dict[str, Any] = {
        "boundary_position": boundary_position,
        "nearest_domain_id": None,
        "nearest_domain_instance_id": None,
        "nearest_domain_instance_number": None,
        "nearest_domain_accession": None,
        "nearest_domain_label": None,
        "nearest_domain_start": None,
        "nearest_domain_end": None,
        "nearest_edge_type": None,
        "nearest_edge_position": None,
        "signed_distance": None,
        "absolute_distance": None,
        "class": UNAVAILABLE,
        "threshold": threshold,
        "source": "representative_domain",
    }
    valid = [d for d in domains if d.get("start") is not None and d.get("end") is not None]
    if boundary_position is None or not valid:
        return rec

    near = _nearest_edge(boundary_position, valid)
    if near is None:
        return rec
    signed, dist, dom, edge_type = near
    edge_pos = dom.get("start") if edge_type == "start" else dom.get("end")
    rec.update({
        "nearest_domain_id": dom.get("id") or instance_id_of(dom),
        "nearest_domain_instance_id": instance_id_of(dom),
        "nearest_domain_instance_number": dom.get("instance_number"),
        "nearest_domain_accession": dom.get("interpro_accession"),
        "nearest_domain_label": dom.get("label") or dom.get("interpro_name"),
        "nearest_domain_start": dom.get("start"),
        "nearest_domain_end": dom.get("end"),
        "nearest_edge_type": edge_type,
        "nearest_edge_position": edge_pos,
        "signed_distance": signed,
        "absolute_distance": dist,
    })

    inside = any(d["start"] <= boundary_position <= d["end"] for d in valid)
    if dist == 0:
        rec["class"] = EXACT
    elif dist <= threshold:
        rec["class"] = NEAR
    elif inside:
        rec["class"] = INSIDE
    else:
        rec["class"] = OUTSIDE
    if rec["class"] == INSIDE and not (dom["start"] <= boundary_position <= dom["end"]):
        # The boundary is inside SOME domain but the nearest edge belongs to another
        # instance. Report the containing instance so the persisted coordinates always
        # describe the instance the class refers to.
        host = min((d for d in valid if d["start"] <= boundary_position <= d["end"]),
                   key=lambda d: min(abs(boundary_position - d["start"]),
                                     abs(boundary_position - d["end"])))
        edge_type = ("start" if abs(boundary_position - host["start"])
                     <= abs(boundary_position - host["end"]) else "end")
        edge_pos = host["start"] if edge_type == "start" else host["end"]
        rec.update({
            "nearest_domain_id": host.get("id") or instance_id_of(host),
            "nearest_domain_instance_id": instance_id_of(host),
            "nearest_domain_instance_number": host.get("instance_number"),
            "nearest_domain_accession": host.get("interpro_accession"),
            "nearest_domain_label": host.get("label") or host.get("interpro_name"),
            "nearest_domain_start": host.get("start"),
            "nearest_domain_end": host.get("end"),
            "nearest_edge_type": edge_type,
            "nearest_edge_position": edge_pos,
            "signed_distance": boundary_position - edge_pos,
            "absolute_distance": abs(boundary_position - edge_pos),
        })
    return rec
