#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

# Accession classification (RefSeq / Ensembl). Curated NM_/NP_ carry manual
# curation; XM_/XP_ are model-predicted; ENS* are Ensembl.
_RE_REFSEQ_CURATED = re.compile(r"^(NM_|NP_|NR_)", re.I)
_RE_REFSEQ_PREDICTED = re.compile(r"^(XM_|XP_|XR_)", re.I)
_RE_ENSEMBL = re.compile(r"^ENS[A-Z]*[TP]\d", re.I)

# Selection rules, in priority order, with a stable display + confidence.
SELECTION_RULES = [
    ("mane_select", "MANE Select", "high"),
    ("appris_principal", "APPRIS principal", "high"),
    ("ensembl_canonical", "Ensembl canonical", "high"),
    ("refseq_curated_over_predicted", "RefSeq curated (NM/NP) over predicted (XM/XP)", "high"),
    ("uniprot_reviewed", "UniProt reviewed / canonical", "medium"),
    ("longest_protein_fallback", "Longest protein (fallback)", "medium"),
]
_RULE_META = {rid: (label, conf) for rid, label, conf in SELECTION_RULES}


def classify_accession(protein_id: str, transcript_id: str = "") -> Dict[str, str]:
    pid = (protein_id or "").strip()
    tid = (transcript_id or "").strip()
    probe = pid or tid
    if _RE_REFSEQ_CURATED.match(probe) or _RE_REFSEQ_CURATED.match(tid):
        return {"source": "refseq_curated", "source_label": "RefSeq curated (NM/NP)", "curated": "yes"}
    if _RE_REFSEQ_PREDICTED.match(probe) or _RE_REFSEQ_PREDICTED.match(tid):
        return {"source": "refseq_predicted", "source_label": "RefSeq predicted (XM/XP)", "curated": "no"}
    if _RE_ENSEMBL.match(probe) or _RE_ENSEMBL.match(tid):
        return {"source": "ensembl", "source_label": "Ensembl", "curated": "unknown"}
    return {"source": "other", "source_label": "Other / unclassified", "curated": "unknown"}


def _to_int(v: Any) -> Optional[int]:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def build_primary_selection(
    isoform_rows: List[Dict[str, str]],
    collection_report: Optional[Dict[str, Any]] = None,
    tags_by_protein: Optional[Dict[str, Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    collection_report = collection_report or {}
    tags_by_protein = tags_by_protein or {}

    # Reference species = the run's first/reference species (single-species runs
    # have exactly one). The top-level primary reports THIS species so that all
    # downstream reference-protein consumers (coordinate overlays, candidate
    # evidence, canonical selected_primary) are mutually consistent.
    ref_species = str(collection_report.get("species_id", "") or "")

    proteins = []
    max_len = -1
    # per-species primary (first primary row seen for each species)
    species_primaries: Dict[str, Dict[str, Any]] = {}
    for r in isoform_rows:
        pid = r.get("protein_id", "")
        if not pid:
            continue
        length = _to_int(r.get("protein_length"))
        sid = str(r.get("species_id", "") or "")
        is_primary = str(r.get("primary_status", "")).lower() == "primary"
        cls = classify_accession(pid, r.get("transcript_id", ""))
        tags = tags_by_protein.get(pid, {})
        proteins.append({
            "protein_id": pid,
            "transcript_id": r.get("transcript_id", ""),
            "species_id": sid,
            "length_aa": length,
            "primary_status": "primary" if is_primary else "alternative",
            "source": cls["source"],
            "source_label": cls["source_label"],
            "curated": cls["curated"],
            "tags": tags,
        })
        if is_primary and sid not in species_primaries:
            species_primaries[sid] = {
                "species_id": sid, "primary_protein_id": pid,
                "primary_transcript_id": r.get("transcript_id", ""),
                "primary_length_aa": length,
            }
        if length is not None and length > max_len:
            max_len = length

    # Top-level primary = reference species' primary when known, else the first
    # primary row (single-species runs are unaffected).
    if ref_species and ref_species in species_primaries:
        primary_id = species_primaries[ref_species]["primary_protein_id"]
    else:
        primary_id = next((p["protein_id"] for p in proteins
                           if p["primary_status"] == "primary"), "")

    primary = next((p for p in proteins if p["protein_id"] == primary_id), None)

    # Determine which rule fired for the selected primary. Prefer explicit tags,
    # then accession curation, then longest fallback.
    def _rule_for(p: Dict[str, Any]) -> str:
        t = p.get("tags", {})
        if t.get("mane") or t.get("mane_select"):
            return "mane_select"
        if t.get("appris") == "principal" or t.get("appris_principal"):
            return "appris_principal"
        if t.get("ensembl_canonical") or t.get("canonical") is True:
            return "ensembl_canonical"
        if p.get("source") == "refseq_curated":
            return "refseq_curated_over_predicted"
        if t.get("uniprot_reviewed"):
            return "uniprot_reviewed"
        return "longest_protein_fallback"

    rule = _rule_for(primary) if primary else "longest_protein_fallback"
    rule_label, confidence = _RULE_META.get(rule, (rule, "medium"))

    # Availability of each hierarchy source (transparent recording).
    any_curated = any(p["source"] == "refseq_curated" for p in proteins)
    any_tags = any(p.get("tags") for p in proteins)
    sources_available = {
        "mane_select": bool(any(p.get("tags", {}).get("mane") for p in proteins)),
        "appris_principal": bool(any(p.get("tags", {}).get("appris") for p in proteins)),
        "ensembl_canonical": bool(any(p.get("tags", {}).get("ensembl_canonical") for p in proteins)),
        "refseq_curated": any_curated,
        "uniprot_reviewed": bool(any(p.get("tags", {}).get("uniprot_reviewed") for p in proteins)),
        "longest_protein": True,
    }

    # alternatives_considered with a reason each (esp. "longer but not curated").
    alternatives = []
    plen = primary["length_aa"] if primary and primary["length_aa"] is not None else None
    # In a multi-species run the alternatives compared against the reference
    # primary are the OTHER isoforms of the SAME (reference) species; other
    # species are reported via ``species_primaries`` and their own per-species view.
    for p in proteins:
        if p["protein_id"] == primary_id:
            continue
        if ref_species and p.get("species_id") and p["species_id"] != ref_species:
            continue
        reason = f"{p['source_label']} isoform"
        if plen is not None and p["length_aa"] is not None:
            if p["length_aa"] > plen and rule == "refseq_curated_over_predicted":
                reason = (f"{p['source_label']}; longer ({p['length_aa']} aa > {plen} aa) but not "
                          f"curated — curated RefSeq primary was preferred")
            elif p["length_aa"] > plen:
                reason = f"{p['source_label']}; longer ({p['length_aa']} aa) but not selected"
        alternatives.append({
            "protein_id": p["protein_id"], "transcript_id": p["transcript_id"],
            "length_aa": p["length_aa"], "source": p["source"], "reason": reason,
        })

    # Human-readable explanation of the selected primary.
    if primary:
        if rule == "refseq_curated_over_predicted":
            longer = [a for a in alternatives if a["length_aa"] and plen and a["length_aa"] > plen]
            explanation = (
                f"{primary['protein_id']} ({primary['transcript_id']}) was selected as primary "
                f"because it is the curated RefSeq protein (NM/NP). Predicted XP isoforms are "
                f"retained as alternatives")
            if longer:
                explanation += (f"; {len(longer)} predicted isoform(s) are longer but were not "
                                f"preferred over the curated protein")
            explanation += "."
        elif rule == "longest_protein_fallback":
            explanation = (
                f"{primary['protein_id']} was selected as the longest protein because no MANE / "
                f"APPRIS / Ensembl-canonical / curated-RefSeq evidence was available for this "
                f"organism. This is a fallback and should be treated with lower confidence.")
        else:
            explanation = (f"{primary['protein_id']} was selected as primary by rule "
                           f"'{rule_label}'.")
    else:
        explanation = "No primary protein could be selected."

    report = {
        "gene_symbol": collection_report.get("gene_symbol", ""),
        "species_id": ref_species or collection_report.get("species_id", ""),
        "species_primaries": list(species_primaries.values()),
        "n_species": len(species_primaries) or 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "selection_method_requested": collection_report.get("selection_method_requested",
                                                             collection_report.get("selection_method", "")),
        "selection_rule": rule,
        "selection_rule_label": rule_label,
        "selection_source": primary["source_label"] if primary else "",
        "evidence_status": "accepted" if primary else "failed",
        "confidence": confidence,
        "primary_protein_id": primary_id,
        "primary_transcript_id": primary["transcript_id"] if primary else "",
        "primary_length_aa": plen,
        "explanation": explanation,
        "hierarchy": [rid for rid, _, _ in SELECTION_RULES],
        "sources_available": sources_available,
        "tags_available": any_tags,
        "n_isoforms": len(proteins),
        "alternatives_considered": alternatives,
        "proteins": proteins,
    }
    return report


# ``species_id`` and ``species_primary`` exist because a multi-species run has one
# primary protein *per species*, not one per run. ``selected_primary`` remains the
# run-level reference protein, which is what single-species consumers have always
# meant by it; readers that must not collapse two species onto one record use
# ``species_primary`` instead. Dropping either column would silently reintroduce the
# defect where a species without the run-level primary gets an arbitrary protein.
_TSV_FIELDS = [
    "species_id", "protein_id", "transcript_id", "selected_primary", "species_primary",
    "selection_rule", "selection_source", "evidence_status", "confidence", "reason",
    "alternatives_considered",
]


def write_selection_evidence(report: Dict[str, Any], tsv_path: Path, json_path: Path) -> None:
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    alt_summary = "; ".join(
        f"{a['protein_id']}({a['length_aa']}aa)" for a in report.get("alternatives_considered", []))
    rule = report.get("selection_rule", "")
    src = report.get("selection_source", "")
    conf = report.get("confidence", "")
    primary_id = report.get("primary_protein_id", "")
    # One primary per species. The report already computes this; writing only the
    # run-level primary is what left every Mus row false in the two-species run and
    # forced the coordinate model into an alphabetical guess.
    per_species = {
        str(v.get("species_id") or ""): str(v.get("primary_protein_id") or "")
        for v in (report.get("species_primaries") or [])
    }
    with open(tsv_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=_TSV_FIELDS, delimiter="\t")
        w.writeheader()
        for p in report.get("proteins", []):
            is_primary = p["protein_id"] == primary_id
            sid = str(p.get("species_id") or "")
            is_species_primary = bool(sid) and per_species.get(sid) == p["protein_id"]
            if is_primary:
                reason = report.get("explanation", "")
                r_rule, r_src, r_conf = rule, src, conf
            else:
                # find the matching alternative reason
                alt = next((a for a in report.get("alternatives_considered", [])
                            if a["protein_id"] == p["protein_id"]), {})
                reason = alt.get("reason", "")
                r_rule, r_src, r_conf = "", p.get("source_label", ""), ""
            w.writerow({
                "species_id": sid,
                "protein_id": p["protein_id"],
                "transcript_id": p.get("transcript_id", ""),
                "selected_primary": "true" if is_primary else "false",
                "species_primary": "true" if is_species_primary else "false",
                "selection_rule": r_rule,
                "selection_source": r_src,
                "evidence_status": report.get("evidence_status", "") if is_primary else "alternative",
                "confidence": r_conf,
                "reason": reason,
                "alternatives_considered": alt_summary if is_primary else "",
            })
    json_path.write_text(__import__("json").dumps(report, indent=2) + "\n", encoding="utf-8")
