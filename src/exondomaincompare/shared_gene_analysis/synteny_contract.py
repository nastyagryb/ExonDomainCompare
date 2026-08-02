"""One canonical local-synteny contract for every gene and every dataset.

Every synteny consumer — the generic single-species view, the generic
multi-species comparative view, the FGFR2 view and every publication export —
reads the shape built here. Gene-specific behaviour is supplied as *data*
(target symbol, orthology classes, review state), never as a second renderer.

Three rules make the contract honest:

* The target locus is a first-class object with its own slot ``0``. It is never
  one of the flanking neighbours and is never counted as one, so a species with
  five upstream and five downstream neighbours reports ten flanking loci, not
  eleven.
* Availability is stated, not implied. ``upstream_count_available`` and
  ``downstream_count_available`` are what the annotation actually yielded;
  ``displayed_*`` is what the view shows. When they differ the row carries a
  ``truncation_status`` and an ``omission_reason``, so a genome with only four
  downstream genes reads as "4 downstream", never as a rendering bug and never
  as a fabricated fifth gene.
* Internal state strings stay internal. Every status and orthology class also
  carries a readable label and an exact definition for a tooltip.
"""
from __future__ import annotations

import re
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

#: How many flanking loci a view asks for on each side. The annotation decides
#: how many actually exist; this is a request, not a promise.
REQUESTED_NEIGHBOUR_COUNT = 5

TARGET_SLOT = 0

_PLACEHOLDER_SYMBOL = re.compile(r"^(?:LOC\d+|GENE\d+|ENS\w*G\d+)$", re.IGNORECASE)

#: Internal orthology class -> (readable label, exact definition).
ORTHOLOGY_DISPLAY: Dict[str, Tuple[str, str]] = {
    "target": ("Target gene",
               "The gene this analysis is about, shown in its own central slot."),
    "exact": ("Resolved ortholog",
              "The neighbouring locus carries a curated gene symbol that matches "
              "the reference symbol exactly."),
    "curated": ("Curated ortholog",
                "The neighbouring locus was assigned by a curated orthology "
                "resource rather than by symbol matching alone."),
    "rbh": ("Best reciprocal hit",
            "The neighbouring locus is the reciprocal best protein hit against "
            "the reference proteome."),
    "weak": ("Weak hit",
             "A protein hit exists but falls below the identity or coverage "
             "threshold used for a confident assignment."),
    "ambiguous": ("Ambiguous paralog",
                  "Several loci hit the same reference gene, so the assignment "
                  "cannot be resolved to one ortholog."),
    "placeholder": ("Placeholder locus",
                    "Placeholder locus label; curated gene symbol unavailable. "
                    "The genomic position is known."),
    "unresolved": ("Unresolved locus",
                   "The locus is annotated in the assembly but no orthology "
                   "assignment was attempted or succeeded."),
}

#: Internal synteny status -> (readable label, exact definition).
STATUS_DISPLAY: Dict[str, Tuple[str, str]] = {
    "synteny_strong": (
        "Synteny strong",
        "All five flanking loci on both sides match the reference neighbourhood "
        "in order and orientation. Supporting evidence for locus identity, not "
        "for event or isoform identity."),
    "synteny_supported": (
        "Synteny supported",
        "The flanking loci match the reference neighbourhood in order and "
        "orientation. Supporting evidence for locus identity, not for event "
        "or isoform identity."),
    "synteny_supported_with_minor_rearrangement": (
        "Supported, minor rearrangement",
        "The same neighbouring loci are present but at least one is reordered "
        "or inverted relative to the reference."),
    "synteny_partial": (
        "Partially supported",
        "Only part of the reference neighbourhood is recovered in this "
        "assembly."),
    "synteny_partial_blast_supported": (
        "Partially supported (protein hits)",
        "Part of the neighbourhood is recovered, and the remaining loci are "
        "supported only by protein similarity, not by a curated symbol."),
    "local_neighbourhood": (
        "Local neighbourhood",
        "Loci flanking the target in this one assembly. No cross-species "
        "conservation is claimed."),
    "unavailable": (
        "Not available",
        "No local gene neighbourhood could be extracted for this species."),
    "review": (
        "Needs review",
        "The neighbourhood was extracted but flagged for manual review."),
}


def is_placeholder_locus(symbol: Any) -> bool:
    """True when the annotation gave an identifier instead of a gene name.

    NCBI assigns ``LOC…`` identifiers to loci that have no approved symbol yet.
    That is a labelling limitation, not a failed orthology assignment, so such a
    locus keeps its own honest style rather than being called unresolved.
    """
    text = str(symbol or "").strip()
    return not text or bool(_PLACEHOLDER_SYMBOL.match(text))


def display_binomial(species_id: str) -> str:
    """``gallus_gallus`` -> ``Gallus gallus`` (genus capitalised, epithet not)."""
    text = str(species_id or "").replace("_", " ").strip()
    return text[:1].upper() + text[1:]


def status_display(status: str) -> Tuple[str, str]:
    """Readable label and exact definition for an internal status string."""
    key = str(status or "").strip()
    if key in STATUS_DISPLAY:
        return STATUS_DISPLAY[key]
    return (key.replace("_", " ").capitalize() or "Not available",
            "No definition is recorded for this status.")


def orthology_display(cls: str) -> Tuple[str, str]:
    """Readable label and exact definition for an internal orthology class."""
    key = str(cls or "").strip() or "unresolved"
    return ORTHOLOGY_DISPLAY.get(key, ORTHOLOGY_DISPLAY["unresolved"])


def _int(value: Any) -> Optional[int]:
    if value in (None, ""):
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _strand(value: Any) -> str:
    text = str(value or "").strip()
    if text.startswith("-"):
        return "-"
    if text.startswith("+"):
        return "+"
    return ""


def neighbour_locus(*, side: str, rank: int, source_symbol: str = "",
                    resolved_symbol: str = "", gene_id: str = "",
                    protein_id: str = "", strand: str = "",
                    orthology_class: str = "", mapping_confidence: str = "",
                    identity_status: str = "", percent_identity: Any = None,
                    coverage: Any = None, distance: Any = None,
                    seqid: str = "", start: Any = None, end: Any = None,
                    method: str = "") -> Dict[str, Any]:
    """One flanking locus. ``rank`` is 1-based distance from the target."""
    rank = abs(int(rank or 0))
    slot = -rank if side == "upstream" else rank
    display = str(resolved_symbol or source_symbol or "").strip()
    placeholder = is_placeholder_locus(display)
    cls = str(orthology_class or "").strip()
    if not cls:
        cls = "placeholder" if placeholder else "exact" if display else "unresolved"
    label, definition = orthology_display(cls)
    return {
        "slot_x": slot,
        "side": side,
        "rank": rank,
        "is_target": False,
        "symbol": display or (gene_id or ""),
        "source_symbol": str(source_symbol or "").strip(),
        "resolved_symbol": str(resolved_symbol or "").strip(),
        "gene_id": gene_id or "",
        "protein_id": protein_id or "",
        "placeholder": placeholder,
        "orthology_class": cls,
        "orthology_label": label,
        "orthology_definition": definition,
        "mapping_confidence": mapping_confidence or _confidence_for(cls),
        "identity_status": identity_status or "",
        "method": method or "",
        "percent_identity": percent_identity,
        "coverage": coverage,
        "strand": _strand(strand),
        "distance": _int(distance),
        "seqid": seqid or "",
        "genomic_start": _int(start),
        "genomic_end": _int(end),
    }


def target_locus(*, gene_symbol: str, gene_id: str = "", strand: str = "",
                 seqid: str = "", start: Any = None, end: Any = None,
                 protein_id: str = "",
                 coordinate_source: str = "annotation") -> Dict[str, Any]:
    """The target gene as its own slot-0 object, explicit and never a neighbour."""
    label, definition = orthology_display("target")
    symbol = str(gene_symbol or "").strip() or "target gene"
    return {
        "slot_x": TARGET_SLOT,
        "side": "target",
        "rank": 0,
        "is_target": True,
        "symbol": symbol,
        "source_symbol": symbol,
        "resolved_symbol": symbol,
        "gene_id": gene_id or "",
        "protein_id": protein_id or "",
        "placeholder": False,
        "orthology_class": "target",
        "orthology_label": label,
        "orthology_definition": definition,
        "mapping_confidence": "target",
        "identity_status": "target",
        "method": "target",
        "percent_identity": None,
        "coverage": None,
        "strand": _strand(strand),
        "distance": 0,
        "seqid": seqid or "",
        "genomic_start": _int(start),
        "genomic_end": _int(end),
        # "annotation" when the run recorded the locus itself; a legacy run that
        # predates the target table reconstructs the span from the neighbour
        # offsets, which is approximate and is labelled as such.
        "coordinate_source": coordinate_source,
    }


def _confidence_for(cls: str) -> str:
    return {
        "target": "target",
        "exact": "high",
        "curated": "high",
        "rbh": "medium",
        "weak": "low",
        "ambiguous": "low",
        "placeholder": "position_only",
        "unresolved": "none",
    }.get(cls, "none")


def counts_label(up: int, down: int) -> str:
    """The honest one-line availability statement shown above a track."""
    total = up + down
    if not total:
        return "No flanking loci available"
    return (f"{total} flanking {'locus' if total == 1 else 'loci'} shown · "
            f"{up} upstream · {down} downstream")


def species_row(species_id: str, *, gene_symbol: str,
                target: Dict[str, Any],
                neighbours: Sequence[Dict[str, Any]],
                display_name: str = "",
                requested: int = REQUESTED_NEIGHBOUR_COUNT,
                synteny_status: str = "local_neighbourhood",
                taxon_group: str = "",
                clade: str = "",
                comparison_available: bool = False,
                is_review: bool = False,
                is_human_reference_control: bool = False,
                extra: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Assemble one canonical species row from already-built locus objects.

    ``neighbours`` may arrive in any order and may contain more loci than
    ``requested``; the row keeps the real order by rank, truncates to the
    request and records what was left out and why.
    """
    up_all = sorted((n for n in neighbours if n.get("side") == "upstream"),
                    key=lambda n: n.get("rank") or 0)
    down_all = sorted((n for n in neighbours if n.get("side") == "downstream"),
                      key=lambda n: n.get("rank") or 0)
    up_avail, down_avail = len(up_all), len(down_all)
    limit = max(0, int(requested or 0)) or max(up_avail, down_avail)
    up = up_all[:limit]
    down = down_all[:limit]

    # Upstream reads outward-to-inward so the display order is left-to-right:
    # farthest upstream, …, nearest upstream, target, nearest downstream, ….
    loci = list(reversed(up)) + [target] + list(down)

    truncated = (len(up) < up_avail) or (len(down) < down_avail)
    short = (len(up) < limit) or (len(down) < limit)
    if truncated:
        status = "truncated_to_request"
        reason = (f"The assembly provides {up_avail} upstream and {down_avail} "
                  f"downstream loci; the view requests {limit} per side.")
    elif short:
        status = "fewer_available"
        reason = (f"The annotation provides only {len(up)} upstream and "
                  f"{len(down)} downstream protein-coding loci around "
                  f"{gene_symbol or 'the target'} on this scaffold.")
    else:
        status = "complete"
        reason = ""

    label, definition = status_display(synteny_status)
    row: Dict[str, Any] = {
        "species_id": species_id,
        "species": species_id,  # legacy key kept for existing consumers
        "display_species_name": display_name or display_binomial(species_id),
        "taxon_group": taxon_group,
        "clade": clade,
        "gene_symbol": gene_symbol,

        "target_gene_id": target.get("gene_id", ""),
        "target_symbol": target.get("symbol", ""),
        "target_slot": TARGET_SLOT,
        "target_strand": target.get("strand", ""),
        "target_position": _position_text(target),
        "target_coordinate_source": target.get("coordinate_source", "annotation"),
        "chromosome_or_scaffold": target.get("seqid", ""),
        "target": target,

        "upstream": list(reversed(up)),
        "downstream": list(down),
        "loci": loci,

        "upstream_count_available": up_avail,
        "downstream_count_available": down_avail,
        "requested_neighbour_count": limit,
        "displayed_upstream_count": len(up),
        "displayed_downstream_count": len(down),
        "displayed_flanking_count": len(up) + len(down),
        "truncation_status": status,
        "omission_reason": reason,
        "counts_label": counts_label(len(up), len(down)),

        "synteny_status": synteny_status,
        "synteny_status_label": label,
        "synteny_status_definition": definition,
        "comparison_available": bool(comparison_available),
        "is_review": bool(is_review),
        "is_human_reference_control": bool(is_human_reference_control),
        "orthology_classes_present": sorted(
            {n.get("orthology_class") for n in loci if n.get("orthology_class")}),
    }
    if extra:
        row.update(extra)
    return row


def _position_text(target: Dict[str, Any]) -> str:
    seqid = target.get("seqid") or ""
    start, end = target.get("genomic_start"), target.get("genomic_end")
    if seqid and start is not None and end is not None:
        return f"{seqid}:{start:,}–{end:,}"
    return seqid


def legacy_nodes(row: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The older flat neighbour shape, so existing readers keep working.

    The canonical ``loci`` array stays authoritative; this projection only adds
    the historical key names that the first synteny viewer was written against.
    """
    out: List[Dict[str, Any]] = []
    for n in row.get("loci") or []:
        out.append({
            **n,
            "raw_symbol": n.get("source_symbol") or n.get("symbol") or "",
            "method_class": _legacy_method_class(n),
            "is_anchor": bool(n.get("is_target")),
            "resolved": n.get("orthology_class") in ("target", "exact", "curated", "rbh"),
        })
    return out


def _legacy_method_class(node: Dict[str, Any]) -> str:
    cls = node.get("orthology_class") or "unresolved"
    return "anchor" if cls == "target" else cls


def summarise(rows: Iterable[Dict[str, Any]]) -> Dict[str, Any]:
    """Dataset-level totals derived from the per-species rows."""
    rows = list(rows)
    flanking = sum(r.get("displayed_flanking_count", 0) for r in rows)
    resolved = sum(
        1 for r in rows for n in (r.get("loci") or [])
        if not n.get("is_target")
        and n.get("orthology_class") in ("exact", "curated", "rbh")
    )
    return {
        "n_species": len(rows),
        "n_flanking_loci": flanking,
        "n_resolved_neighbours": resolved,
        "classes_present": sorted({c for r in rows
                                   for c in r.get("orthology_classes_present", [])}),
        "any_truncated": any(r.get("truncation_status") == "truncated_to_request"
                             for r in rows),
        "any_incomplete": any(r.get("truncation_status") == "fewer_available"
                              for r in rows),
    }
