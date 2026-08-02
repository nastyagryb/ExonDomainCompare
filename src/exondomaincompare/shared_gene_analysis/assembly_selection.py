"""assembly_selection.py — choose a genome assembly, and say why the others lost.

The Equus quagga run wrote two assembly tables with zero rows each. Zero rows was the
only thing they could say, because every failure path in the collector returned an
empty list: a rejected taxonomy query, a species with no assembly, a species whose
assemblies are all unannotated, a download that failed and an archive that would not
parse all produced the same empty table. The actual cause — NCBI refused the query
term ``equus_quagga`` — was invisible, and an annotated chromosome-level RefSeq
reference assembly (GCF_021613505.1) went unnoticed.

Those states are separated here, and every candidate keeps a row saying whether it was
selected or rejected and on what ground. An empty candidate list now means one thing
only: the source returned no assembly for this taxon.

The priority rule is deterministic and stated in one place:

1. an annotated assembly, always — an unannotated one has no genes to find;
2. RefSeq (``GCF_``) over GenBank (``GCA_``), because RefSeq is the annotated copy;
3. the registry's own assembly preference;
4. a reference or representative genome over an alternate one;
5. assembly level: complete genome > chromosome > scaffold > contig;
6. accession, so that ties break the same way on every run.

Nothing here is species-specific. The rule is the same for a zebra, a chicken and any
species added later.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

#: Outcomes. The point of naming this many is that each one implies a different next
#: action: fix the name, wait and retry, accept that the species has no annotation, or
#: investigate a parser.
NO_QUERY_TERM = "no_taxon_query_term"
TAXON_REJECTED = "taxon_not_recognised_by_source"
SERVICE_FAILED = "assembly_service_failed"
NO_ASSEMBLY = "no_assembly_returned"
NONE_ANNOTATED = "assemblies_returned_but_none_annotated"
SELECTED = "annotated_assembly_selected"
ALL_REJECTED = "annotated_assembly_returned_but_rejected"
DOWNLOAD_FAILED = "assembly_download_failed"
PARSE_FAILED = "assembly_archive_parse_failed"

#: Statuses that mean "there is nothing to download, and it is not our mistake".

_LEVEL_RANK = {
    "complete genome": 40,
    "chromosome": 35,
    "scaffold": 20,
    "contig": 10,
}

_SUPPRESSED = {"suppressed", "replaced", "withdrawn"}


@dataclass
class AssemblySelection:
    """The decision, the candidates and the reason — one object, one story."""

    status: str
    selected: Optional[Dict[str, Any]] = None
    candidates: List[Dict[str, Any]] = field(default_factory=list)
    detail: str = ""
    query_term: str = ""

    @property
    def n_annotated(self) -> int:
        return sum(1 for c in self.candidates if c.get("annotated") == "1")

    def summary(self) -> Dict[str, Any]:
        return {
            "status": self.status,
            "query_term": self.query_term,
            "n_candidates": len(self.candidates),
            "n_annotated": self.n_annotated,
            "selected_accession": (self.selected or {}).get("accession", ""),
            "detail": self.detail,
        }


def _norm(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip()).lower()


def normalise_candidate(report: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """One assembly report in the fields the selector needs.

    Tolerant about where a field sits, because the Datasets JSON has moved fields
    between releases before. Missing *optional* metadata is recorded as empty and does
    not disqualify a candidate; a missing accession does, since there is then nothing
    to download.
    """
    info = report.get("assembly_info") or report.get("assembly") or {}
    organism = report.get("organism") or {}
    annotation = report.get("annotation_info") or {}

    accession = (info.get("assembly_accession") or report.get("accession")
                 or report.get("assembly_accession") or "")
    if not accession:
        return None

    annotation_name = annotation.get("name") or annotation.get("release_name") or ""
    paired = (info.get("paired_assembly") or {}).get("accession", "")
    return {
        "accession": str(accession),
        "assembly_name": str(info.get("assembly_name") or report.get("assembly_name") or ""),
        "assembly_level": str(info.get("assembly_level") or report.get("assembly_level") or ""),
        "refseq_category": str(info.get("refseq_category") or report.get("refseq_category") or ""),
        "assembly_source": str(info.get("assembly_source") or report.get("assembly_source") or ""),
        "assembly_status": str(info.get("assembly_status") or report.get("assembly_status") or ""),
        "annotated": "1" if annotation_name else "0",
        "annotation_name": str(annotation_name),
        "annotation_release_date": str(annotation.get("release_date") or ""),
        "paired_accession": str(paired or ""),
        "organism_name": str(organism.get("organism_name") or ""),
        "tax_id": str(organism.get("tax_id") or ""),
    }


def parse_summary(reports: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for report in reports:
        row = normalise_candidate(report)
        if row is not None:
            out.append(row)
    return out


def _taxon_match(candidate: Dict[str, Any], requested_taxid: str,
                 requested_name: str) -> str:
    """How this assembly's organism relates to the requested taxon.

    A taxon query returns descendants, so a subspecies assembly can appear under a
    species request. That is legitimate and stays available, but it is labelled: a
    reader should be able to see that a chosen assembly is a subspecies rather than
    the nominate form.
    """
    taxid = str(candidate.get("tax_id") or "")
    if requested_taxid and taxid == str(requested_taxid):
        return "exact_taxon"
    name = _norm(candidate.get("organism_name"))
    request = _norm(requested_name)
    if request and name == request:
        return "exact_taxon"
    if request and name.startswith(request + " "):
        return "descendant_taxon"
    if requested_taxid and taxid and taxid != str(requested_taxid):
        return "related_taxon"
    return "unverified_taxon"


def score_candidate(candidate: Dict[str, Any], preference: str) -> Dict[str, Any]:
    accession = candidate.get("accession", "")
    is_refseq = accession.startswith("GCF_")
    is_genbank = accession.startswith("GCA_")
    annotated = candidate.get("annotated") == "1"
    category = _norm(candidate.get("refseq_category"))
    reference = "reference genome" in category or "representative genome" in category
    level = _LEVEL_RANK.get(_norm(candidate.get("assembly_level")), 0)
    pref = _norm(preference) or "refseq"
    pref_match = ((pref == "refseq" and is_refseq) or (pref == "genbank" and is_genbank)
                  or pref in {"any", "none", "best_available", ""})

    notes: List[str] = []
    if annotated:
        notes.append("annotated")
    else:
        notes.append("no_annotation")
    if is_refseq:
        notes.append("refseq")
    if pref_match:
        notes.append("matches_registry_assembly_preference")
    if reference:
        notes.append("reference_or_representative")
    notes.append(f"level={_norm(candidate.get('assembly_level')) or 'unknown'}")

    scored = dict(candidate)
    scored["assembly_score"] = str(
        annotated * 1000 + is_refseq * 500 + pref_match * 200 + reference * 100 + level)
    scored["assembly_decision_notes"] = ";".join(notes)
    return scored


def select(candidates: Sequence[Dict[str, Any]], *, preference: str = "RefSeq",
           requested_taxid: str = "", requested_name: str = "",
           query_term: str = "") -> AssemblySelection:
    """Choose one assembly, and record a decision for every candidate."""
    if not candidates:
        return AssemblySelection(
            status=NO_ASSEMBLY, query_term=query_term,
            detail=("the source returned no genome assembly for this taxon"))

    scored = [score_candidate(c, preference) for c in candidates]
    for row in scored:
        row["taxon_match"] = _taxon_match(row, requested_taxid, requested_name)

    def _rank(row: Dict[str, Any]):
        # Score and taxon match descending, then accession ascending so that a tie
        # resolves to the same assembly on every run rather than to whichever
        # accession happens to sort last.
        return (-int(row["assembly_score"]),
                0 if row["taxon_match"] == "exact_taxon" else 1,
                row["accession"])

    ordered = sorted(scored, key=_rank)

    # Suppressed, replaced and withdrawn assemblies are disqualified before ranking is
    # consulted. Deciding this inside the acceptance loop made the label depend on sort
    # position: a suppressed assembly that happened to rank below an accepted one was
    # reported as merely lower priority, hiding the fact that NCBI has withdrawn it.
    for row in ordered:
        status = _norm(row.get("assembly_status"))
        if status in _SUPPRESSED:
            row["decision"] = "rejected"
            row["rejection_reason"] = f"assembly_status_{status}"

    annotated = [r for r in ordered
                 if r.get("annotated") == "1" and "decision" not in r]
    if not annotated:
        # Not a defect: plenty of genomes are sequenced and never annotated. Selecting
        # one anyway would download an archive with no gene records and produce the
        # empty tables this module exists to explain.
        n_unannotated = 0
        for row in ordered:
            if "decision" in row:
                continue
            row["decision"] = "rejected"
            row["rejection_reason"] = "assembly_has_no_annotation_release"
            n_unannotated += 1
        if n_unannotated == 0:
            return AssemblySelection(
                status=ALL_REJECTED, candidates=ordered, query_term=query_term,
                detail=("every assembly for this taxon has been withdrawn or replaced "
                        "by the source: "
                        + "; ".join(f"{r['accession']}={r['rejection_reason']}"
                                    for r in ordered)))
        return AssemblySelection(
            status=NONE_ANNOTATED, candidates=ordered, query_term=query_term,
            detail=(f"{len(ordered)} assembly/assemblies exist for this taxon but none "
                    "carries a usable annotation release, so no gene model can be "
                    "derived from them"))

    accepted: Optional[Dict[str, Any]] = annotated[0]

    for row in ordered:
        if row is accepted:
            row["decision"] = "selected"
            row["rejection_reason"] = ""
        elif "decision" not in row:
            row["decision"] = "rejected"
            row["rejection_reason"] = (
                "assembly_has_no_annotation_release" if row.get("annotated") != "1"
                else "lower_priority_than_selected_assembly")

    if accepted is None:
        return AssemblySelection(
            status=ALL_REJECTED, candidates=ordered, query_term=query_term,
            detail=("every annotated assembly for this taxon was rejected: "
                    + "; ".join(f"{r['accession']}={r['rejection_reason']}"
                                for r in annotated)))

    detail = (f"selected {accepted['accession']} "
              f"({accepted.get('assembly_level') or 'unknown level'}, "
              f"annotation {accepted.get('annotation_name') or 'unnamed'}, "
              f"{accepted['taxon_match']})")
    if accepted["taxon_match"] == "descendant_taxon":
        detail += (f" — organism is {accepted.get('organism_name')}, a subspecies of "
                   "the requested taxon")
    return AssemblySelection(status=SELECTED, selected=accepted, candidates=ordered,
                             query_term=query_term, detail=detail)


def selection_rows(selection: AssemblySelection, species_input: str,
                   species_canonical: str, taxid: str) -> List[Dict[str, str]]:
    """The candidate table. Written even when nothing was selected, because the
    rejected candidates are the evidence for why."""
    rows: List[Dict[str, str]] = []
    for row in selection.candidates:
        rows.append({
            "species_input": species_input,
            "species_canonical": species_canonical,
            "taxid": taxid or row.get("tax_id", ""),
            "selection_status": selection.status,
            **{k: str(v) for k, v in row.items()},
        })
    return rows


def failure_row(status: str, detail: str, species_input: str, species_canonical: str,
                taxid: str, query_term: str) -> Dict[str, str]:
    """A single row for a failure that produced no candidates at all.

    Without this the table is empty and a reader cannot tell a rejected query term
    from a species that genuinely has no assembly — the exact ambiguity that hid the
    Equus quagga root cause.
    """
    return {
        "species_input": species_input,
        "species_canonical": species_canonical,
        "taxid": taxid,
        "selection_status": status,
        "accession": "",
        "assembly_name": "",
        "assembly_level": "",
        "refseq_category": "",
        "annotated": "0",
        "annotation_name": "",
        "organism_name": "",
        "tax_id": taxid,
        "taxon_match": "not_evaluated",
        "decision": "none_available",
        "rejection_reason": status,
        "assembly_score": "",
        "assembly_decision_notes": detail,
        "query_term": query_term,
    }
