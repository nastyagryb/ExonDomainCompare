#!/usr/bin/env python3
"""
final_pre_interpro_validation.py  (Sprint Part 1)

Single validation gate used by all final plotting/report scripts.

Before any figure or report is generated, this module verifies that every
required final pre-InterPro file exists, is non-empty, and contains the
required columns. If anything is missing or inconsistent, it fails with a
clear error message instead of silently producing misleading figures.

Importable API:
    locate_files(base) -> dict
    run_validation(base, outdir=None) -> (ok: bool, rows: list[dict], summary: dict)
    validate_or_raise(base) -> dict   # raises RuntimeError on hard failure

Outputs (when --outdir / metadata dir is provided):
    final_pre_interpro_validation_report.tsv
    final_pre_interpro_validation_report.json
"""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

SCRIPT_VERSION = "1.0"

# logical filename -> (required columns, dir hint substring, kind)
# kind: "tsv" | "fasta" | "text"
REQUIRED_SPEC: List[Tuple[str, List[str], str, str]] = [
    ("species_qc_master.tsv",
     ["species", "final_display_class", "recommended_use", "phylo_order", "taxon_group", "fgfr2_ortholog_status"],
     "11_pre_interpro_master", "tsv"),
    ("species_qc_master_pre_interpro.tsv",
     ["species", "final_display_class", "phylo_order"], "11_pre_interpro_master", "tsv"),
    ("fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv",
     ["species_canonical", "inferred_isoform"], "", "tsv"),
    ("fgfr2_pair_level_qc_summary.tsv",
     ["species_canonical", "has_both_isoforms", "pair_audit_status"], "", "tsv"),
    ("fgfr2_resolved_IIIb_IIIc_exon_CDS_mapping.tsv",
     ["species_canonical", "inferred_isoform", "cds_boundary_precision_refined"], "", "tsv"),
    ("fgfr2_resolved_IIIb_IIIc_candidate_scores.tsv",
     ["species_canonical", "transcript_id_source"], "", "tsv"),
    ("fgfr2_III_pair_audit.tsv",
     ["species", "has_IIIb", "has_IIIc"], "", "tsv"),
    ("fgfr2_orthology_evidence.tsv",
     ["species", "isoform", "orthology_status"], "", "tsv"),
    ("fgfr2_orthology_species_summary.tsv",
     ["species", "orthology_status_species", "recommended_use"], "", "tsv"),
    ("cds_phase_boundary_audit.tsv",
     ["species", "isoform", "cds_boundary_precision_refined", "cds_boundary_confidence",
      "reason_if_unknown", "reason_if_split", "transcript_cds_reconstruction_status",
      "protein_translation_check_status"], "", "tsv"),
    ("cds_phase_boundary_legacy_vs_refined_counts.tsv",
     ["precision_category", "legacy_count", "refined_count"], "", "tsv"),
    ("cds_phase_boundary_explainability_summary.tsv",
     ["dimension", "category", "count"], "", "tsv"),
    ("fgfr2_cassette_cds_block_map.tsv",
     ["species", "isoform", "matched_cds_rank", "matched_protein_start_aa",
      "cassette_overlap_status"], "", "tsv"),
    ("fgfr2_cassette_coordinate_sanity_audit.tsv",
     ["species", "isoform", "cassette_start_aa", "coordinate_sanity_status"], "", "tsv"),
    ("fgfr2_transcript_cds_reconstruction_audit.tsv",
     ["species", "isoform", "reconstruction_status", "translation_matches_selected_protein"], "", "tsv"),
    ("fgfr2_refined_uncertainty_classes.tsv",
     ["species", "isoform", "coordinate_resolution_state", "boundary_precision_state",
      "protein_evidence_state", "annotation_review_state", "display_uncertainty_class",
      "plot_visibility_level", "uncertainty_explanation_short", "uncertainty_explanation_long"],
     "", "tsv"),
    ("cds_phase_rescue_audit.tsv",
     ["species", "isoform", "original_boundary_precision", "rescue_attempted",
      "rescued_boundary_precision", "rescue_status"], "", "tsv"),
    ("fgfr2_ncbi_cds_boundary_patch_report.tsv",
     ["species", "isoform", "issue_type", "patch_attempted", "patch_status",
      "patch_used_in_final", "translation_validation_status"], "", "tsv"),
    ("fgfr2_interpro_clean_unique.fasta", [], "", "fasta"),
    ("fgfr2_interpro_unique_mapping.tsv",
     ["unique_id", "sequence_hash", "representative_clean_id"], "", "tsv"),
    ("fgfr2_interpro_id_mapping.tsv",
     ["clean_id", "unique_id", "is_unique_representative"], "", "tsv"),
    ("fgfr2_interpro_prepare_summary.tsv", ["metric", "value"], "", "tsv"),
    ("fgfr2_interpro_prepare_warnings.tsv", [], "", "tsv"),
    ("interproscan_input_manifest.tsv", ["file", "role"], "", "tsv"),
    ("interproscan_run_instructions.md", [], "", "text"),
]


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def locate(base: Path, name: str, hint: str) -> Optional[Path]:
    matches = sorted(base.rglob(name))
    if not matches:
        return None
    if hint:
        for m in matches:
            if hint in str(m):
                return m
    # prefer shallowest path otherwise
    return sorted(matches, key=lambda p: len(p.parts))[0]


def _read_header(path: Path) -> List[str]:
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        line = fh.readline().rstrip("\n").rstrip("\r")
    return line.split("\t") if line else []


def _count_tsv_rows(path: Path) -> int:
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return max(0, sum(1 for _ in fh) - 1)


def _count_fasta(path: Path) -> int:
    n = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def _distinct_species(path: Path, col: str) -> int:
    seen = set()
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v = (r.get(col) or "").strip()
            if v:
                seen.add(v.lower())
    return len(seen)


def run_validation(base: Path, outdir: Optional[Path] = None) -> Tuple[bool, List[Dict[str, object]], Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    hard_fail = False
    paths: Dict[str, Optional[Path]] = {}

    for name, req_cols, hint, kind in REQUIRED_SPEC:
        p = locate(base, name, hint)
        paths[name] = p
        exists = p is not None and p.exists()
        n = ""
        cols_present = ""
        status = "ok"
        warning = ""

        if not exists:
            status, warning, hard_fail = "fail", "file_missing", True
        elif p.stat().st_size == 0:
            status, warning, hard_fail = "fail", "file_empty", True
        else:
            if kind == "fasta":
                n = _count_fasta(p)
                if n == 0:
                    status, warning, hard_fail = "fail", "no_fasta_records", True
            elif kind == "text":
                n = sum(1 for _ in open(p, encoding="utf-8", errors="replace"))
            else:
                n = _count_tsv_rows(p)
                header = _read_header(p)
                missing = [c for c in req_cols if c not in header]
                if missing:
                    status, warning, hard_fail = "fail", f"missing_columns:{','.join(missing)}", True
                    cols_present = "false"
                else:
                    cols_present = "true"
                if n == 0 and kind == "tsv" and name not in ("fgfr2_interpro_prepare_warnings.tsv",):
                    if status == "ok":
                        status, warning = "warning", "zero_data_rows"

        rows.append({
            "file": name,
            "path": str(p) if p else "",
            "exists": "true" if exists else "false",
            "n_rows_or_records": n,
            "required_columns_present": cols_present,
            "status": status,
            "warning": warning,
        })

    # ---- cross-file consistency checks (do not silently pass) ----
    extra: List[Dict[str, object]] = []

    master = paths.get("species_qc_master.tsv")
    alias = paths.get("species_qc_master_pre_interpro.tsv")
    if master and alias and master.exists() and alias.exists():
        identical = master.read_bytes() == alias.read_bytes()
        extra.append({
            "file": "CHECK:master_alias_identical", "path": str(alias),
            "exists": "true", "n_rows_or_records": "",
            "required_columns_present": "", "status": "ok" if identical else "warning",
            "warning": "" if identical else "alias_not_byte_identical_documented_as_alias",
        })

    # orthology/paralog integrated into master
    if master and master.exists():
        hdr = _read_header(master)
        integrated = all(c in hdr for c in ("fgfr2_ortholog_status", "paralog_screen_status"))
        extra.append({
            "file": "CHECK:orthology_paralog_integrated_in_master", "path": str(master),
            "exists": "true", "n_rows_or_records": "",
            "required_columns_present": "true" if integrated else "false",
            "status": "ok" if integrated else "fail",
            "warning": "" if integrated else "orthology_paralog_not_integrated",
        })
        if not integrated:
            hard_fail = True

    # final resolved IIIb/IIIc mapping must remain 60 rows / 30 species
    mp = paths.get("fgfr2_resolved_IIIb_IIIc_exon_CDS_mapping.tsv")
    if mp and mp.exists():
        nrows = _count_tsv_rows(mp)
        nsp = _distinct_species(mp, "species_canonical")
        ok_map = (nrows == 60 and nsp == 30)
        extra.append({
            "file": "CHECK:resolved_mapping_60rows_30species", "path": str(mp),
            "exists": "true", "n_rows_or_records": f"{nrows}rows/{nsp}species",
            "required_columns_present": "", "status": "ok" if ok_map else "warning",
            "warning": "" if ok_map else f"expected_60rows_30species_got_{nrows}_{nsp}_justify_in_report",
        })

    # orthology evidence > 60 rows must be documented
    oe = paths.get("fgfr2_orthology_evidence.tsv")
    if oe and oe.exists():
        nrows = _count_tsv_rows(oe)
        extra.append({
            "file": "CHECK:orthology_evidence_row_count", "path": str(oe),
            "exists": "true", "n_rows_or_records": nrows,
            "required_columns_present": "", "status": "ok",
            "warning": "" if nrows <= 60 else
            "more_than_60_records_retained_ambiguous_review_documented_in_report",
        })

    # cassette coordinate sanity: no implausible NON-review cassette may remain
    sa = paths.get("fgfr2_cassette_coordinate_sanity_audit.tsv")
    if sa and sa.exists():
        bad = []
        with open(sa, encoding="utf-8", errors="replace", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                st = str(row.get("coordinate_sanity_status", ""))
                if ("implausible" in st or "first_cds_block" in st) and not st.endswith("_review_excluded"):
                    bad.append(f"{row.get('species')}/{row.get('isoform')}:{st}")
        extra.append({
            "file": "CHECK:cassette_coordinate_sanity", "path": str(sa),
            "exists": "true", "n_rows_or_records": len(bad),
            "required_columns_present": "", "status": "ok" if not bad else "fail",
            "warning": "" if not bad else f"implausible_non_review_cassettes:{bad[:8]}",
        })
        if bad:
            hard_fail = True

    # cassette map: no main-analysis cassette mapped to protein_start_aa==1
    cmap = paths.get("fgfr2_cassette_cds_block_map.tsv")
    if cmap and cmap.exists():
        n1 = 0
        with open(cmap, encoding="utf-8", errors="replace", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if str(row.get("matched_protein_start_aa", "")).strip() == "1":
                    n1 += 1
        extra.append({
            "file": "CHECK:no_cassette_at_protein_start_1", "path": str(cmap),
            "exists": "true", "n_rows_or_records": n1,
            "required_columns_present": "", "status": "ok" if n1 == 0 else "fail",
            "warning": "" if n1 == 0 else f"{n1}_cassettes_at_protein_start_aa_1_mapping_bug",
        })
        if n1 > 0:
            hard_fail = True

    # refined uncertainty classes: every row must be explainable (AC1) and carry a
    # known plot_visibility_level (AC8 relies on it)
    refp = paths.get("fgfr2_refined_uncertainty_classes.tsv")
    if refp and refp.exists():
        valid_vis = {"none", "subtle_symbol", "supplement_only", "main_warning", "hard_fail"}
        bad = []
        with open(refp, encoding="utf-8", errors="replace", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if (not str(row.get("display_uncertainty_class", "")).strip()
                        or not str(row.get("uncertainty_explanation_short", "")).strip()
                        or str(row.get("plot_visibility_level", "")).strip() not in valid_vis):
                    bad.append(f"{row.get('species')}/{row.get('isoform')}")
        extra.append({
            "file": "CHECK:refined_uncertainty_explainable", "path": str(refp),
            "exists": "true", "n_rows_or_records": len(bad),
            "required_columns_present": "", "status": "ok" if not bad else "fail",
            "warning": "" if not bad else f"unexplained_or_bad_visibility:{bad[:8]}",
        })
        if bad:
            hard_fail = True

    # phase rescue attempted for every phase-unavailable case (AC2)
    rescp = paths.get("cds_phase_rescue_audit.tsv")
    if rescp and rescp.exists():
        missing = []
        with open(rescp, encoding="utf-8", errors="replace", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                orig = str(row.get("original_boundary_precision", "")).lower()
                reason = str(row.get("original_reason_if_unknown", ""))
                needs = ("unknown" in orig or reason in (
                    "phase_not_propagated_from_source", "missing_gff3_phase",
                    "nucleotide_sequence_unavailable"))
                if needs and str(row.get("rescue_attempted", "")).strip().lower() != "true":
                    missing.append(f"{row.get('species')}/{row.get('isoform')}")
        extra.append({
            "file": "CHECK:phase_rescue_attempted_for_unavailable", "path": str(rescp),
            "exists": "true", "n_rows_or_records": len(missing),
            "required_columns_present": "", "status": "ok" if not missing else "fail",
            "warning": "" if not missing else f"phase_unavailable_without_rescue_attempt:{missing[:8]}",
        })
        if missing:
            hard_fail = True

    # no fake InterPro domain coordinates: ensure no domain-coordinate file present
    fake = list(base.rglob("*interpro_domain*coord*")) + list(base.rglob("*interpro*domains*.tsv"))
    extra.append({
        "file": "CHECK:no_fake_interpro_domain_coordinates", "path": "",
        "exists": "true", "n_rows_or_records": len(fake),
        "required_columns_present": "", "status": "ok" if not fake else "fail",
        "warning": "" if not fake else f"unexpected_domain_coordinate_files:{[str(f) for f in fake]}",
    })
    if fake:
        hard_fail = True

    rows.extend(extra)

    n_ok = sum(1 for r in rows if r["status"] == "ok")
    n_warn = sum(1 for r in rows if r["status"] == "warning")
    n_fail = sum(1 for r in rows if r["status"] == "fail")
    summary = {
        "script_version": SCRIPT_VERSION, "base": str(base), "generated_at": _now(),
        "n_checks": len(rows), "n_ok": n_ok, "n_warning": n_warn, "n_fail": n_fail,
        "overall_status": "fail" if hard_fail else ("warning" if n_warn else "pass"),
    }

    if outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        fields = ["file", "path", "exists", "n_rows_or_records",
                  "required_columns_present", "status", "warning"]
        with open(outdir / "final_pre_interpro_validation_report.tsv", "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
            w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, "") for k in fields})
        with open(outdir / "final_pre_interpro_validation_report.json", "w", encoding="utf-8") as fh:
            json.dump({"summary": summary, "checks": rows}, fh, indent=2)

    return (not hard_fail), rows, summary


def validate_or_raise(base: Path) -> Dict[str, object]:
    ok, rows, summary = run_validation(base, outdir=None)
    if not ok:
        failed = [f"{r['file']}: {r['warning']}" for r in rows if r["status"] == "fail"]
        raise RuntimeError(
            "Final pre-InterPro validation FAILED. Refusing to generate figures.\n  - "
            + "\n  - ".join(failed)
        )
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description="Final pre-InterPro validation gate (Part 1).")
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, default=None,
                    help="metadata dir for report (default: <base>/11_publication_figures_pre_interpro/metadata)")
    args = ap.parse_args()
    outdir = args.outdir or (args.base / "11_publication_figures_pre_interpro" / "metadata")
    ok, rows, summary = run_validation(args.base, outdir=outdir)
    print(f"[validation] {summary['overall_status'].upper()} "
          f"(ok={summary['n_ok']}, warn={summary['n_warning']}, fail={summary['n_fail']})")
    for r in rows:
        if r["status"] != "ok":
            print(f"  [{r['status']}] {r['file']} :: {r['warning']}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
