#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from exondomaincompare.adapters.fgfr2_to_generic_indices import DatasetSource, PROJECT_ROOT, FREEZE_ROOT, write_json, _projection_coords, _int_or_none
from exondomaincompare.framework.gene_config import GeneConfig, detector_for_analysis

DETECTOR_NAME = "fgfr2_iiib_iiic"
DETECTOR_VERSION = "1.0"
CONTRACT_VERSION = 1
SOURCE_LABEL = "legacy_fgfr2_adapter"


def _cassette_agreement(src: DatasetSource) -> Dict[tuple, Dict[str, Any]]:
    cassette = src.idx("cassette_residue_index.json", {}) or {}
    out: Dict[tuple, Dict[str, Any]] = {}
    for sp in cassette.get("species", []) or []:
        sp_id = sp.get("species", "")
        for label, panel in (sp.get("panels", {}) or {}).items():
            positions = panel.get("positions", []) or []
            disc = [p for p in positions if p.get("is_discriminating")]
            ident = [p for p in disc if p.get("agreement_class") == "identical_to_human"
                     or p.get("cls") == "identical"]
            out[(sp_id, label)] = {
                "disc_total": len(disc),
                "disc_identical": len(ident),
                "available": bool(panel.get("available")),
            }
    return out


def _tsv_write(path: Path, columns: List[str], rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in columns})


def build_detector_outputs(src: DatasetSource, cfg: GeneConfig, out_dir: Path) -> Dict[str, Any]:
    species_index = src.idx("species_index.json", []) or []
    coords = _projection_coords(src.closure)
    agreement = _cassette_agreement(src)
    marker_by_label = {lab["id"]: lab.get("marker_sequence", "")
                       for lab in cfg.event_labels}
    proj_rel = "results/13_final_pre_interpro_closure/MSA/final_cassette_msa_boundary_projection.tsv"

    candidates: List[Dict[str, Any]] = []
    region_rows: List[Dict[str, Any]] = []
    evidence_rows: List[Dict[str, Any]] = []

    n_species = 0
    n_candidates = 0
    n_primary = 0
    warnings: List[str] = []

    base = {"analysis_id": cfg.analysis_id, "gene_symbol": cfg.gene_symbol,
            "event_id": cfg.event_id}

    for sp in species_index if isinstance(species_index, list) else []:
        n_species += 1
        sp_id = sp.get("species", "")
        for iso in sp.get("isoforms", []) or []:
            label = iso.get("final_isoform_label") or iso.get("isoform") or ""
            protein_id = iso.get("protein_id", "")
            n_candidates += 1
            included = str(iso.get("interpro_included", "")).lower()
            if included == "primary":
                n_primary += 1

            candidates.append({
                **base,
                "species_id": sp_id,
                "event_label": label,
                "protein_id": protein_id,
                "transcript_id": iso.get("transcript_id", ""),
                "gene_id": iso.get("gene_id", ""),
                "sequence_id": iso.get("sequence_md5", ""),
                "protein_length": _int_or_none(iso.get("protein_length")),
                "candidate_status": iso.get("final_claim_status_after_rescue", ""),
                "evidence_level": iso.get("readiness_class", "") or included,
                "source": SOURCE_LABEL,
                "notes": "",
            })

            c = coords.get((sp_id, iso.get("isoform", "")), {})
            start = c.get("region_start_aa")
            end = c.get("region_end_aa")
            length = (end - start + 1) if (isinstance(start, int) and isinstance(end, int)) else None
            if start is None or end is None:
                warnings.append(f"missing region coordinates for {sp_id}/{label}")
            region_rows.append({
                **base,
                "species_id": sp_id,
                "event_label": label,
                "protein_id": protein_id,
                "region_start_aa": start,
                "region_end_aa": end,
                "region_length_aa": length,
                "coordinate_basis": "native_protein_aa",
                "coordinate_confidence": "high" if length else "unknown",
                "source_table": proj_rel,
                "notes": "",
            })

            ag = agreement.get((sp_id, label), {})
            disc_total = ag.get("disc_total", 0)
            disc_identical = ag.get("disc_identical", 0)
            ref_sim = (round(disc_identical / disc_total, 4) if disc_total else None)
            is_reference = (sp_id == cfg.reference_species)
            evidence_rows.append({
                **base,
                "species_id": sp_id,
                "event_label": label,
                "protein_id": protein_id,
                "marker_sequence": marker_by_label.get(label, ""),
                "marker_status": ("reference" if is_reference
                                  else ("evaluated" if disc_total else "not_evaluated")),
                "sequence_evidence": iso.get("final_claim_status_after_rescue", ""),
                "reference_similarity": ref_sim,
                "validated_label": iso.get("final_isoform_label", label),
                "confidence": iso.get("readiness_class", ""),
                "source": SOURCE_LABEL,
                "notes": "",
            })

    if not species_index:
        warnings.append("species_index.json is empty or missing; no candidates projected.")

    cand_cols = ["analysis_id", "gene_symbol", "species_id", "event_id", "event_label",
                 "protein_id", "transcript_id", "gene_id", "sequence_id", "protein_length",
                 "candidate_status", "evidence_level", "source", "notes"]
    region_cols = ["analysis_id", "gene_symbol", "species_id", "event_id", "event_label",
                   "protein_id", "region_start_aa", "region_end_aa", "region_length_aa",
                   "coordinate_basis", "coordinate_confidence", "source_table", "notes"]
    ev_cols = ["analysis_id", "gene_symbol", "species_id", "event_id", "event_label",
               "protein_id", "marker_sequence", "marker_status", "sequence_evidence",
               "reference_similarity", "validated_label", "confidence", "source", "notes"]

    _tsv_write(out_dir / "event_isoform_candidates.tsv", cand_cols, candidates)
    _tsv_write(out_dir / "event_region_coordinates.tsv", region_cols, region_rows)
    _tsv_write(out_dir / "event_label_evidence.tsv", ev_cols, evidence_rows)

    cassette = src.idx("cassette_residue_index.json", {}) or {}
    reference_used = bool(cassette.get("human_reference"))

    report = {
        "detector_name": DETECTOR_NAME,
        "detector_version": DETECTOR_VERSION,
        "contract_version": CONTRACT_VERSION,
        "analysis_id": cfg.analysis_id,
        "gene_config": cfg.source_path,
        "dataset_id": src.dataset_id,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "input_files": [
            "website_indices/species_index.json",
            proj_rel,
            "website_indices/cassette_residue_index.json",
        ],
        "output_files": [
            "event_isoform_candidates.tsv",
            "event_region_coordinates.tsv",
            "event_label_evidence.tsv",
            "event_detector_report.json",
        ],
        "n_species": n_species,
        "n_candidate_proteins": n_candidates,
        "n_accepted_primary_proteins": n_primary,
        "warnings": sorted(set(warnings)),
        "failures": [],
        "reference_control_used": reference_used,
        "note": "Projection of existing validated FGFR2 outputs; no biology recomputed.",
    }
    write_json(out_dir / "event_detector_report.json", report)
    return report


def _resolve_out_dir(src: DatasetSource, out_arg: Optional[str]) -> Path:
    if out_arg:
        out = Path(out_arg)
    elif src.kind == "run":
        out = src.run_root / "results" / "generic_event_detector"
    else:
        out = PROJECT_ROOT / "artifacts" / "generic_event_detector" / "example"
    if not out.is_absolute():
        out = PROJECT_ROOT / out
    if str(out.resolve()).startswith(str(FREEZE_ROOT.resolve())):
        raise SystemExit(
            f"Refusing to write detector outputs inside the example freeze: {out}. "
            "Use --out to choose a safe location (e.g. artifacts/generic_event_detector/example).")
    return out


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Project validated FGFR2 outputs into the generic event-detector contract.")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--run-id", help="Custom run id under runs/.")
    g.add_argument("--example", action="store_true", help="Read the example dataset (read-only).")
    ap.add_argument("--config", help="Override gene_config.yaml path (default: resolve/FGFR2).")
    ap.add_argument("--out", help="Output directory for detector contract outputs.")
    args = ap.parse_args(argv)

    src = DatasetSource.for_example() if args.example else DatasetSource.for_run(args.run_id)
    cfg = src.gene_config(args.config)

    # Guard: only run the FGFR2 detector for the analysis it supports.
    det = detector_for_analysis(cfg.analysis_id)
    if det is None or det.get("detector") != DETECTOR_NAME:
        print(f"NOTE: analysis '{cfg.analysis_id}' has no supported '{DETECTOR_NAME}' detector; "
              "projecting with the FGFR2 adapter anyway (dataset is FGFR2-derived).",
              file=sys.stderr)

    out_dir = _resolve_out_dir(src, args.out)
    report = build_detector_outputs(src, cfg, out_dir)

    print(f"OK  detector={DETECTOR_NAME}  dataset={src.dataset_id}  analysis={cfg.analysis_id}")
    print(f"    species={report['n_species']}  candidates={report['n_candidate_proteins']}  "
          f"primary={report['n_accepted_primary_proteins']}  ref_control={report['reference_control_used']}")
    print(f"    out: {out_dir}")
    for f in report["output_files"]:
        print(f"      - {f}")
    if report["warnings"]:
        print(f"    warnings: {len(report['warnings'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
