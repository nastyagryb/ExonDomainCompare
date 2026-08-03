#!/usr/bin/env python3
"""
Build stable pre-InterPro plotting tables.

These tables do not require InterProScan results. Plotting reads them verbatim
and does not recompute the upstream biological QC.

Outputs (in --outdir):
  figure1_framework_counts_pre_interpro.tsv
  figure2_exon_to_protein_tracks_pre_interpro.tsv
  figure3_species_evidence_matrix_pre_interpro.tsv
  figure4_native_vs_normalized_qc_pre_interpro.tsv
  figure_review_cases_pre_interpro.tsv
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional



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


def g(d: Dict[str, str], key: str, default: str = "") -> str:
    v = d.get(key, default)
    return v if v not in (None, "") else default


def build_figure1(master: List[Dict[str, str]]) -> List[Dict[str, object]]:
    """Framework / category counts as a tidy (category, level, count) table."""
    rows: List[Dict[str, object]] = []
    dims = [
        ("final_display_class", "final_display_class"),
        ("taxon_group", "taxon_group"),
        ("main_analysis_eligible", "main_analysis_eligible"),
        ("both_isoforms_detected", "both_isoforms_detected"),
        ("fgfr2_ortholog_status", "fgfr2_ortholog_status"),
        ("paralog_screen_status", "paralog_screen_status"),
        ("direction_validation_status", "direction_validation_status"),
        ("cds_boundary_precision_summary", "cds_boundary_precision_summary"),
        ("iii_region_similarity_class", "iii_region_similarity_class"),
        ("interpro_status", "interpro_status"),
    ]
    for cat, col in dims:
        counts = Counter(g(m, col, "unknown") for m in master)
        for level, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
            rows.append({"category": cat, "level": level, "count": n,
                         "fraction": round(n / max(1, len(master)), 4)})
    return rows


def build_figure2(pairqc: List[Dict[str, str]], master_by_sp: Dict[str, Dict[str, str]]) -> List[Dict[str, object]]:
    """Exon-to-protein tracks: one tidy row per (species, isoform)."""
    rows: List[Dict[str, object]] = []
    for r in pairqc:
        sp = g(r, "species_canonical") or g(r, "species")
        m = master_by_sp.get(sp, {})
        for iso in ("IIIb", "IIIc"):
            rows.append({
                "species": sp,
                "display_species_name": g(m, "display_species_name", sp),
                "taxon_group": g(m, "taxon_group", "unknown"),
                "isoform": iso,
                "native_start_aa": g(r, f"{iso}_native_start_aa"),
                "native_end_aa": g(r, f"{iso}_native_end_aa"),
                "native_center_aa": g(r, f"{iso}_native_center_aa"),
                "native_length_aa": g(r, f"{iso}_native_len_aa"),
                "iii_slot_center_aa": g(r, f"{iso}_iii_slot_center_aa"),
                "iii_slot_length_aa": g(r, f"{iso}_iii_slot_len_aa"),
                "resolver_evidence_level": g(r, f"{iso}_evidence_level"),
                "resolver_status_refined": g(r, f"{iso}_refined"),
                "cds_boundary_precision_refined": g(r, f"{iso}_cds_boundary_precision_refined"),
                "cds_left_boundary_precision": g(r, f"{iso}_cds_left_boundary_precision"),
                "cds_right_boundary_precision": g(r, f"{iso}_cds_right_boundary_precision"),
                "transcript_id": g(r, f"{iso}_transcript_id_source"),
                "protein_id": g(r, f"{iso}_protein_id"),
                "final_display_class": g(m, "final_display_class", "unknown"),
                "main_analysis_eligible": g(m, "main_analysis_eligible", "false"),
                "interpro_domain_track": "InterProScan_pending",
            })
    return rows


def build_figure3(master: List[Dict[str, str]]) -> List[Dict[str, object]]:
    """Species x evidence matrix in tidy long form (species, dimension, status, ok_flag)."""
    rows: List[Dict[str, object]] = []
    dims = [
        ("fgfr2_ortholog", "fgfr2_ortholog_status", lambda v: "high_confidence" in v or "probable" in v),
        ("paralog_screen", "paralog_screen_status", lambda v: "high_confidence" in v),
        ("both_isoforms", "both_isoforms_detected", lambda v: v == "true"),
        ("direction", "direction_validation_status", lambda v: "unresolved" not in v and "ambiguous" not in v),
        ("protein_validation", "protein_validation_summary", lambda v: "review" not in v and "conflict" not in v),
        ("native_coordinate", "native_coordinate_sanity", lambda v: not v.endswith("review") and v != "unresolved"),
        ("normalized_slot", "normalized_slot_sanity", lambda v: "review" not in v and "unresolved" not in v),
        ("iii_similarity", "iii_region_similarity_class", lambda v: "review" not in v and v != "unresolved"),
        ("cds_boundary", "cds_boundary_precision_summary", lambda v: v == "exact" or v == "codon_split_one_side"),
        ("interpro_input", "interpro_status", lambda v: v == "interpro_ready_input_prepared"),
    ]
    for m in master:
        sp = g(m, "species")
        for dim_name, col, ok in dims:
            val = g(m, col, "unknown")
            rows.append({
                "species": sp,
                "display_species_name": g(m, "display_species_name", sp),
                "taxon_group": g(m, "taxon_group", "unknown"),
                "evidence_dimension": dim_name,
                "status": val,
                "ok_flag": 1 if ok(val) else 0,
                "final_display_class": g(m, "final_display_class", "unknown"),
            })
    return rows


def build_figure4(pairqc: List[Dict[str, str]], master_by_sp: Dict[str, Dict[str, str]]) -> List[Dict[str, object]]:
    """Native vs normalized III-slot pair QC per species."""
    rows: List[Dict[str, object]] = []
    for r in pairqc:
        sp = g(r, "species_canonical") or g(r, "species")
        m = master_by_sp.get(sp, {})
        rows.append({
            "species": sp,
            "display_species_name": g(m, "display_species_name", sp),
            "taxon_group": g(m, "taxon_group", "unknown"),
            "native_pair_center_distance_aa": g(r, "native_pair_center_distance_aa"),
            "iii_slot_pair_center_distance_aa": g(r, "iii_slot_pair_center_distance_aa"),
            "native_coordinate_sanity": g(r, "native_coordinate_sanity"),
            "normalized_slot_sanity": g(r, "iii_slot_coordinate_sanity"),
            "iii_region_similarity_class": g(r, "iii_region_similarity_class"),
            "final_display_class": g(m, "final_display_class", "unknown"),
            "main_analysis_eligible": g(m, "main_analysis_eligible", "false"),
        })
    return rows


def build_review(master: List[Dict[str, str]]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for m in master:
        if g(m, "final_display_class") == "main_analysis_high_confidence":
            continue
        rows.append({
            "species": g(m, "species"),
            "display_species_name": g(m, "display_species_name"),
            "taxon_group": g(m, "taxon_group", "unknown"),
            "final_display_class": g(m, "final_display_class"),
            "review_reason_short": g(m, "review_reason_short"),
            "review_reason_long": g(m, "review_reason_long"),
            "recommended_use": g(m, "recommended_use"),
            "direction_validation_status": g(m, "direction_validation_status"),
            "protein_validation_summary": g(m, "protein_validation_summary"),
            "native_coordinate_sanity": g(m, "native_coordinate_sanity"),
            "cds_boundary_precision_summary": g(m, "cds_boundary_precision_summary"),
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description="Build pre-InterPro figure input tables.")
    ap.add_argument("--master", type=Path, required=True, help="species_qc_master_pre_interpro.tsv")
    ap.add_argument("--pair_qc", type=Path, required=True, help="fgfr2_pair_level_qc_summary.tsv")
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    master = read_tsv(args.master)
    pairqc = read_tsv(args.pair_qc)
    master_by_sp = {g(m, "species"): m for m in master}

    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)

    f1 = build_figure1(master)
    write_tsv(out / "figure1_framework_counts_pre_interpro.tsv", f1, ["category", "level", "count", "fraction"])

    f2 = build_figure2(pairqc, master_by_sp)
    write_tsv(out / "figure2_exon_to_protein_tracks_pre_interpro.tsv", f2, list(f2[0].keys()) if f2 else
              ["species", "isoform"])

    f3 = build_figure3(master)
    write_tsv(out / "figure3_species_evidence_matrix_pre_interpro.tsv", f3,
              ["species", "display_species_name", "taxon_group", "evidence_dimension", "status", "ok_flag", "final_display_class"])

    f4 = build_figure4(pairqc, master_by_sp)
    write_tsv(out / "figure4_native_vs_normalized_qc_pre_interpro.tsv", f4, list(f4[0].keys()) if f4 else
              ["species"])

    fr = build_review(master)
    write_tsv(out / "figure_review_cases_pre_interpro.tsv", fr,
              ["species", "display_species_name", "taxon_group", "final_display_class", "review_reason_short",
               "review_reason_long", "recommended_use", "direction_validation_status", "protein_validation_summary",
               "native_coordinate_sanity", "cds_boundary_precision_summary"])

    print(f"[OK] figure1 rows={len(f1)} figure2 rows={len(f2)} figure3 rows={len(f3)} figure4 rows={len(f4)} review rows={len(fr)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
