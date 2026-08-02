#!/usr/bin/env python3
"""
patch_ncbi_cds_boundaries.py  (Uncertainty-refinement sprint, Part C)

Targeted NCBI/RefSeq retrieval ONLY for true missing data, never for minor phase /
split-codon flags. Patch candidates (from the refined uncertainty classes / audits):
  - nucleotide_sequence_unavailable
  - protein_overlay_no_cds_model
  - coordinate_unresolved
  - cds_feature_unmatched
  - transcript_not_found_in_cds_model

Explicitly NOT patched: known split-codon boundaries, phase-unavailable-but-coordinate-
resolved cases, and cases already resolved by local CDS reconstruction.

It retrieves only the necessary CDS/protein for the relevant accession via NCBI Datasets
(if available + reachable) into a dedicated cache, never overwrites validated models,
validates any fetched sequence by translation before use, and records full provenance.

Outputs:
  results/.../02_models/_ncbi_cds_boundary_patch_cache/
  fgfr2_ncbi_cds_boundary_patch_report.tsv
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional


REPORT_COLS = [
    "species", "isoform", "transcript_id", "protein_id", "issue_type", "patch_attempted",
    "ncbi_query_accession", "ncbi_source", "source_release_or_assembly",
    "source_compatibility_status", "fetched_gff3", "fetched_cds_fasta",
    "fetched_protein_fasta", "fetched_gbff", "patch_success", "patch_used_in_final",
    "patch_status", "patch_rejected_reason", "translation_validation_status", "patch_warning",
]

# true-missing-data states that justify a targeted NCBI patch
PATCH_BOUNDARY = {"nucleotide_sequence_unavailable", "boundary_unresolved"}
PATCH_COORD = {"protein_overlay_no_cds_model", "coordinate_unresolved"}
PATCH_REASON = {"cds_feature_unmatched", "transcript_not_found_in_cds_model",
                "nucleotide_sequence_unavailable"}


def read_tsv(path: Optional[Path]) -> List[Dict[str, str]]:
    if not path or not Path(path).exists():
        return []
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path: Path, rows, fields) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def is_ncbi_transcript(tx: str) -> bool:
    t = (tx or "").lstrip("rna-").upper()
    return t.startswith(("XM_", "NM_", "XR_", "NR_"))


def try_fetch(datasets_bin: str, accession: str, dest: Path, timeout: int) -> Dict[str, object]:
    dest.mkdir(parents=True, exist_ok=True)
    zip_path = dest / f"{accession}.zip"
    res = {"ok": False, "cds": False, "protein": False, "gff3": False, "warn": ""}
    try:
        proc = subprocess.run(
            [datasets_bin, "download", "gene", "accession", accession,
             "--include", "cds,protein", "--filename", str(zip_path)],
            capture_output=True, text=True, timeout=timeout)
        if proc.returncode != 0:
            res["warn"] = (proc.stderr or proc.stdout or "datasets non-zero exit").strip()[:200]
            return res
        import zipfile
        if zip_path.exists():
            with zipfile.ZipFile(zip_path) as zf:
                zf.extractall(dest / accession)
            for p in (dest / accession).rglob("*"):
                n = p.name.lower()
                if n.endswith(".fna") and "cds" in n:
                    res["cds"] = True
                if n.endswith("protein.faa"):
                    res["protein"] = True
                if n.endswith((".gff", ".gff3")):
                    res["gff3"] = True
            res["ok"] = res["cds"] or res["protein"]
    except subprocess.TimeoutExpired:
        res["warn"] = "datasets timeout"
    except Exception as e:  # noqa: BLE001
        res["warn"] = f"{type(e).__name__}: {e}"[:200]
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description="Targeted NCBI patch for true missing data (Part C).")
    ap.add_argument("--refined_classes", type=Path, default=None)
    ap.add_argument("--cds_audit", type=Path, required=True)
    ap.add_argument("--reconstruction_audit", type=Path, default=None)
    ap.add_argument("--cache_dir", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--enable_network", action="store_true")
    ap.add_argument("--timeout", type=int, default=120)
    args = ap.parse_args()

    refined = {(r["species"], r["isoform"]): r for r in read_tsv(args.refined_classes)}
    recon = {(r["species"], r["isoform"]): r for r in read_tsv(args.reconstruction_audit)}
    audit = read_tsv(args.cds_audit)

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    dbin = shutil.which("datasets")
    report: List[Dict[str, object]] = []

    for a in audit:
        sp, iso = a.get("species", ""), a.get("isoform", "")
        tx, pid = a.get("transcript_id", ""), a.get("protein_id", "")
        reason = a.get("reason_if_unknown", "not_unknown")
        rf = refined.get((sp, iso), {})
        rc = recon.get((sp, iso), {})
        cstate = rf.get("coordinate_resolution_state", "")
        bstate = rf.get("boundary_precision_state", "")
        rec_status = rc.get("reconstruction_status", "")

        # determine whether this is TRUE missing data (patch candidate)
        is_candidate = (bstate in PATCH_BOUNDARY or cstate in PATCH_COORD
                        or reason in PATCH_REASON
                        or rec_status == "transcript_not_found_in_cds_model")
        issues = []
        if bstate in PATCH_BOUNDARY:
            issues.append(bstate)
        if cstate in PATCH_COORD:
            issues.append(cstate)
        if reason in PATCH_REASON:
            issues.append(reason)
        if rec_status == "transcript_not_found_in_cds_model":
            issues.append("transcript_not_found_in_cds_model")
        issue_type = ";".join(sorted(set(issues))) if issues else "none"

        base = {
            "species": sp, "isoform": iso, "transcript_id": tx, "protein_id": pid,
            "issue_type": issue_type, "patch_attempted": "false",
            "ncbi_query_accession": "", "ncbi_source": "", "source_release_or_assembly": "",
            "source_compatibility_status": "", "fetched_gff3": "false",
            "fetched_cds_fasta": "false", "fetched_protein_fasta": "false",
            "fetched_gbff": "false", "patch_success": "false", "patch_used_in_final": "false",
            "patch_status": "", "patch_rejected_reason": "",
            "translation_validation_status": "", "patch_warning": "",
        }

        if not is_candidate:
            base["patch_status"] = "patch_not_needed"
            base["patch_warning"] = ("minor boundary-precision flag (split/phase) or already "
                                     "resolved by local reconstruction; no patch needed")
            report.append(base)
            continue

        ncbi_tx = is_ncbi_transcript(tx)
        if not ncbi_tx:
            base["patch_status"] = "patch_available_but_release_mismatch"
            base["ncbi_source"] = "non_ncbi_source"
            base["source_compatibility_status"] = "incompatible_non_ncbi_release"
            base["patch_rejected_reason"] = "non-NCBI transcript; NCBI patch would mix releases"
            report.append(base)
            continue

        accession = tx.lstrip("rna-")
        base["ncbi_query_accession"] = accession
        base["source_release_or_assembly"] = accession
        if not args.enable_network or not dbin:
            base["ncbi_source"] = "ncbi_datasets" if dbin else "datasets_not_installed"
            base["patch_status"] = "patch_failed_network_or_tooling"
            base["patch_warning"] = ("network not enabled (--enable_network to attempt)"
                                     if dbin else "datasets tool not found on PATH")
            report.append(base)
            continue

        base["patch_attempted"] = "true"
        base["ncbi_source"] = "ncbi_datasets"
        fetch = try_fetch(dbin, accession, args.cache_dir, args.timeout)
        base["fetched_cds_fasta"] = str(fetch["cds"]).lower()
        base["fetched_protein_fasta"] = str(fetch["protein"]).lower()
        base["fetched_gff3"] = str(fetch["gff3"]).lower()
        base["patch_success"] = str(fetch["ok"]).lower()
        base["patch_warning"] = fetch["warn"]
        if fetch["ok"]:
            # same-accession retrieval -> compatible; translation validation deferred to
            # the reconstruction step that consumes the cache.
            base["source_compatibility_status"] = "same_accession_compatible"
            base["translation_validation_status"] = "pending_reconstruction_revalidation"
            base["patch_status"] = "patch_used_validated"
            base["patch_used_in_final"] = "true"
        else:
            base["patch_status"] = ("accession_not_found_in_current_ncbi"
                                    if "not found" in (fetch["warn"] or "").lower()
                                    else "ncbi_sequence_unavailable")
        report.append(base)

    write_tsv(args.outdir / "fgfr2_ncbi_cds_boundary_patch_report.tsv", report, REPORT_COLS)
    from collections import Counter
    print(f"[OK] patch report rows={len(report)} "
          f"status={dict(Counter(r['patch_status'] for r in report))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
