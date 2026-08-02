#!/usr/bin/env python3
"""
build_species_qc_master_pre_interpro.py  (Task 8)

Build ``species_qc_master_pre_interpro.tsv`` -- the single, pre-InterProScan
source of truth for per-species display and QC. It joins the already-validated
per-step outputs (registry, Step 4 direction calibration, Step 6b paralog screen,
Step 5b protein QC, Step 9 paper-ready QC, Step 10 pair-level resolver QC, Step 7
InterPro preparation) into one tidy, stable table.

This script performs **no** biological recomputation. It only aggregates and
classifies existing, upstream-validated evidence so that figures, reports and the
manuscript all read display classes from one place. InterProScan has not been run,
so no domain calls are present; ``interpro_status`` only reflects input readiness.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

SCRIPT_NAME = "build_species_qc_master_pre_interpro.py"
SCRIPT_VERSION = "1.0"

COLUMNS = [
    "species", "display_species_name", "taxon_group", "fgfr2_ortholog_status",
    "paralog_screen_status", "both_isoforms_detected", "direction_validation_status",
    "direction_confidence", "protein_validation_summary", "protein_validation_review_count",
    "resolver_status_summary", "native_coordinate_sanity", "normalized_slot_sanity",
    "iii_region_similarity_class", "cds_boundary_precision_summary", "main_analysis_eligible",
    "final_display_class", "review_reason_short", "review_reason_long", "recommended_use",
    "interpro_status",
    # Sprint Part 2: integrated phylogenetic/taxonomic ordering.
    "taxid", "taxon_group_display", "major_clade", "phylo_order",
    "phylo_order_source", "phylo_order_confidence",
    # CDS-boundary sprint Part A: explainable codon-phase/boundary summary.
    "cds_boundary_explainability",
]

# Lowest-wins confidence ranking for aggregation.
_CONF_RANK = {"high": 3, "medium": 2, "moderate": 2, "low": 1, "unresolved": 0, "": 0}
# 5b validation statuses that should count as review (not clean validation).
_REVIEW_PROTEIN = {
    "protein_conflicts_expected_isoform",
    "protein_ambiguous_inconclusive",
    "unresolved_segment_extraction",
    "systematic_isoform_inversion_review",
    "protein_supports_expected_isoform_below_strict_threshold",
}


def read_tsv(path: Optional[Path]) -> List[Dict[str, str]]:
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return []
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path: Path, rows: List[Dict[str, object]], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def norm(s: str) -> str:
    return str(s or "").strip().lower().replace(" ", "_")


def display_name(species_id: str, scientific: str) -> str:
    if scientific:
        return scientific
    parts = str(species_id or "").split("_")
    if len(parts) >= 2:
        return parts[0].capitalize() + " " + " ".join(parts[1:])
    return species_id


def boolish(v: str) -> bool:
    return str(v).strip().lower() in ("1", "true", "yes", "y")


def main() -> int:
    ap = argparse.ArgumentParser(description="Build the pre-InterPro species QC master table (Task 8).")
    ap.add_argument("--registry", type=Path, required=True)
    ap.add_argument("--isoform_evidence", type=Path, required=True, help="Step 4 fgfr2_isoform_evidence.tsv")
    ap.add_argument("--paralog_summary", type=Path, default=None, help="Step 6b paralog screen species summary")
    ap.add_argument("--protein_validation_summary", type=Path, required=True, help="Step 5b final selected protein validation summary")
    ap.add_argument("--paper_ready_qc", type=Path, required=True, help="Step 9 paper-ready species QC")
    ap.add_argument("--pair_qc", type=Path, required=True, help="Step 10 pair-level QC summary")
    ap.add_argument("--interpro_id_mapping", type=Path, default=None, help="Step 7 InterPro id mapping")
    ap.add_argument("--orthology_summary", type=Path, default=None, help="Addendum B orthology species summary")
    ap.add_argument("--phylo_order", type=Path, default=None, help="Sprint Part 2 species_phylogenetic_order.tsv")
    ap.add_argument("--cds_audit", type=Path, default=None, help="Sprint Part A cds_phase_boundary_audit.tsv")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--prefix", default="fgfr2")
    args = ap.parse_args()

    registry = {norm(r.get("species_id") or r.get("ncbi_species") or r.get("ensembl_species")): r
                for r in read_tsv(args.registry)}
    # Also index registry by scientific name normalized for fallback joins.
    reg_by_sci = {norm(r.get("scientific_name")): r for r in read_tsv(args.registry)}

    paralog = {norm(r.get("species")): r for r in read_tsv(args.paralog_summary)}

    # Step 4: aggregate direction status/confidence per species.
    dir_status: Dict[str, Set[str]] = defaultdict(set)
    dir_conf: Dict[str, List[str]] = defaultdict(list)
    for r in read_tsv(args.isoform_evidence):
        sp = norm(r.get("species_canonical") or r.get("species"))
        if not sp:
            continue
        st = str(r.get("direction_validation_status", "")).strip()
        if st:
            dir_status[sp].add(st)
        cf = str(r.get("direction_confidence", "")).strip()
        if cf:
            dir_conf[sp].append(cf)

    # Step 5b: aggregate protein validation per species.
    prot_status: Dict[str, List[str]] = defaultdict(list)
    for r in read_tsv(args.protein_validation_summary):
        sp = norm(r.get("species"))
        if sp:
            prot_status[sp].append(str(r.get("validation_status", "")).strip())

    paper = {norm(r.get("species")): r for r in read_tsv(args.paper_ready_qc)}
    pairqc = {norm(r.get("species_canonical") or r.get("species")): r for r in read_tsv(args.pair_qc)}

    # Addendum B: dedicated orthology evidence (preferred over paralog-only derivation).
    orthology = {norm(r.get("species")): r for r in read_tsv(args.orthology_summary)}

    # Sprint Part 2: integrated phylogenetic/taxonomic order.
    phylo = {norm(r.get("species")): r for r in read_tsv(args.phylo_order)}

    # CDS-boundary sprint Part A: per-species codon-phase/boundary explainability.
    cds_by_sp: Dict[str, List[Dict[str, str]]] = {}
    if args.cds_audit and Path(args.cds_audit).exists():
        for r in read_tsv(args.cds_audit):
            cds_by_sp.setdefault(norm(r.get("species")), []).append(r)

    def cds_explain(sp: str) -> str:
        rows_sp = cds_by_sp.get(sp, [])
        if not rows_sp:
            return "no_cds_audit"
        unk = [r for r in rows_sp if str(r.get("reason_if_unknown")) != "not_unknown"]
        spl = [r for r in rows_sp if str(r.get("reason_if_split")) != "not_split"]
        if not unk and not spl:
            return "all_boundaries_known"
        parts = []
        if unk:
            reasons = sorted({str(r.get("reason_if_unknown")) for r in unk})
            parts.append(f"unknown_x{len(unk)}:{'/'.join(reasons)}")
        if spl:
            parts.append(f"split_x{len(spl)}")
        return "; ".join(parts)

    # Step 7: species with prepared InterPro input.
    interpro_species: Set[str] = set()
    interpro_available = bool(args.interpro_id_mapping and Path(args.interpro_id_mapping).exists())
    for r in read_tsv(args.interpro_id_mapping):
        sp = norm(r.get("species_canonical") or r.get("species"))
        if sp:
            interpro_species.add(sp)

    species_universe = sorted(set(paper) | set(pairqc) | set(dir_status) | set(prot_status))

    rows: List[Dict[str, object]] = []
    for sp in species_universe:
        reg = registry.get(sp) or reg_by_sci.get(sp.replace("_", " ")) or {}
        pa = paper.get(sp, {})
        pq = pairqc.get(sp, {})
        pl = paralog.get(sp, {})

        # Direction.
        d_set = sorted(dir_status.get(sp, set()))
        d_status = d_set[0] if len(d_set) == 1 else (";".join(d_set) if d_set else "direction_unresolved_no_sequence")
        confs = dir_conf.get(sp, [])
        d_conf = min(confs, key=lambda c: _CONF_RANK.get(norm(c), 0)) if confs else ""

        # Protein QC.
        pv = [s for s in prot_status.get(sp, []) if s]
        review_count = sum(1 for s in pv if s in _REVIEW_PROTEIN)
        if not pv:
            prot_summary = "no_protein_validation_records"
        elif review_count == 0:
            prot_summary = "all_selected_proteins_validated_or_transcript_only"
        elif any(s == "protein_conflicts_expected_isoform" for s in pv):
            prot_summary = "protein_conflict_present_review"
        else:
            prot_summary = "protein_ambiguous_or_below_threshold_review"

        # Paralog / ortholog.
        paralog_status = str(pl.get("species_fgfr2_screen_status", "")).strip() or "not_paralog_screened"
        orth = orthology.get(sp, {})
        if orth.get("orthology_status_species"):
            # Prefer the dedicated multi-evidence orthology layer (Addendum B).
            ortholog_status = str(orth.get("orthology_status_species"))
        elif "high_confidence" in paralog_status:
            ortholog_status = "fgfr2_ortholog_high_confidence"
        elif "probable" in paralog_status:
            ortholog_status = "fgfr2_ortholog_probable"
        elif paralog_status == "not_paralog_screened":
            ortholog_status = "fgfr2_ortholog_assumed_not_screened"
        else:
            ortholog_status = "fgfr2_ortholog_review"

        both_iso = boolish(pq.get("has_both_isoforms", "")) if pq else False
        main_eligible = boolish(pa.get("main_analysis_eligible", "0"))

        native_sanity = str(pq.get("native_coordinate_sanity", "")).strip() or "unresolved"
        slot_sanity = str(pq.get("iii_slot_coordinate_sanity", "")).strip() or "normalized_III_slot_unresolved"
        sim_class = (str(pq.get("iii_region_similarity_class", "")).strip()
                     or str(pa.get("regional_local_identity", "")).strip() or "unresolved")
        cds_summary = str(pq.get("cds_boundary_precision_summary", "")).strip() or "unknown_codon_phase"
        resolver_summary = str(pq.get("resolver_status_refined_set", "")).strip() or "unresolved"

        # InterPro readiness.
        if not interpro_available:
            interpro_status = "interpro_pending"
        elif sp in interpro_species:
            interpro_status = "interpro_ready_input_prepared"
        else:
            interpro_status = "interpro_input_missing"

        # Derive display class + review reasons.
        reasons: List[str] = []
        if not both_iso:
            reasons.append("both_IIIb_IIIc_isoforms_not_jointly_resolved")
        if review_count > 0:
            reasons.append(f"protein_QC_review_x{review_count}")
        if native_sanity.endswith("review"):
            reasons.append(f"native_coordinate:{native_sanity}")
        if "review" in sim_class:
            reasons.append(f"iii_similarity:{sim_class}")
        if cds_summary in ("codon_split_both_sides", "unknown_codon_phase"):
            reasons.append(f"cds_boundary:{cds_summary}")
        if d_status not in ("sequence_calibrated_consistent_with_order_rule",
                            "sequence_calibrated_inverted_from_order_rule") and "unresolved" in d_status:
            reasons.append(f"direction:{d_status}")

        if main_eligible and not reasons:
            display_class = "main_analysis_high_confidence"
            recommended = "main_text_primary_claim"
        elif main_eligible:
            display_class = "main_analysis_with_minor_review"
            recommended = "main_text_with_footnote"
        else:
            display_class = "supplementary_review_not_primary_claim"
            recommended = "supplementary_only"

        review_short = "ok" if not reasons else reasons[0].split(":")[0]
        review_long = "; ".join(reasons) if reasons else "no_review_flags"

        # Sprint Part 2: phylogenetic/taxonomic order; broad taxon_group overrides raw clade.
        ph = phylo.get(sp, {})
        broad_group = str(ph.get("taxon_group", "")).strip() or (str(reg.get("clade", "")).strip() or "unknown")

        rows.append({
            "species": sp,
            "display_species_name": display_name(sp, str(reg.get("scientific_name", ""))),
            "taxon_group": broad_group,
            "fgfr2_ortholog_status": ortholog_status,
            "paralog_screen_status": paralog_status,
            "both_isoforms_detected": "true" if both_iso else "false",
            "direction_validation_status": d_status,
            "direction_confidence": d_conf or "unknown",
            "protein_validation_summary": prot_summary,
            "protein_validation_review_count": review_count,
            "resolver_status_summary": resolver_summary,
            "native_coordinate_sanity": native_sanity,
            "normalized_slot_sanity": slot_sanity,
            "iii_region_similarity_class": sim_class,
            "cds_boundary_precision_summary": cds_summary,
            "main_analysis_eligible": "true" if main_eligible else "false",
            "final_display_class": display_class,
            "review_reason_short": review_short,
            "review_reason_long": review_long,
            "recommended_use": recommended,
            "interpro_status": interpro_status,
            "taxid": str(ph.get("taxid", "")).strip() or str(reg.get("taxid", "")).strip(),
            "taxon_group_display": str(ph.get("taxon_group_display", "")).strip() or broad_group,
            "major_clade": str(ph.get("major_clade", "")).strip(),
            "phylo_order": str(ph.get("phylo_order", "")).strip(),
            "phylo_order_source": str(ph.get("order_source", "")).strip(),
            "phylo_order_confidence": str(ph.get("order_confidence", "")).strip(),
            "cds_boundary_explainability": cds_explain(sp),
        })

    # Canonical row order follows reproducible phylogenetic order when available.
    def _po(r: Dict[str, str]) -> Tuple[int, str]:
        v = str(r.get("phylo_order", "")).strip()
        return (int(v), r["species"]) if v.isdigit() else (10 ** 6, r["species"])

    rows.sort(key=_po)
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    # Addendum E: species_qc_master.tsv is the CANONICAL display/QC table.
    # species_qc_master_pre_interpro.tsv is an identical copy / documented alias
    # for the pre-InterProScan stage (no InterPro domain calls present yet).
    canonical_path = out / "species_qc_master.tsv"
    master_path = out / "species_qc_master_pre_interpro.tsv"
    write_tsv(canonical_path, rows, COLUMNS)
    write_tsv(master_path, rows, COLUMNS)

    from collections import Counter
    meta = {
        "script_name": SCRIPT_NAME, "script_version": SCRIPT_VERSION,
        "run_datetime_utc": datetime.now(timezone.utc).isoformat(),
        "n_species": len(rows),
        "canonical_master": "species_qc_master.tsv",
        "pre_interpro_alias": "species_qc_master_pre_interpro.tsv",
        "alias_note": "species_qc_master_pre_interpro.tsv is an identical copy of the canonical species_qc_master.tsv for the pre-InterProScan stage",
        "final_display_class_counts": dict(Counter(r["final_display_class"] for r in rows)),
        "main_analysis_eligible_counts": dict(Counter(r["main_analysis_eligible"] for r in rows)),
        "interpro_status_counts": dict(Counter(r["interpro_status"] for r in rows)),
        "both_isoforms_detected_counts": dict(Counter(r["both_isoforms_detected"] for r in rows)),
        "review_species": [r["species"] for r in rows if r["final_display_class"] != "main_analysis_high_confidence"],
        "inputs": {k: str(v) for k, v in vars(args).items() if isinstance(v, Path)},
    }
    (out / "species_qc_master_pre_interpro_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"[OK] wrote canonical {canonical_path} and alias {master_path} ({len(rows)} species)")
    print("final_display_class:", meta["final_display_class_counts"])
    print("interpro_status:", meta["interpro_status_counts"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
