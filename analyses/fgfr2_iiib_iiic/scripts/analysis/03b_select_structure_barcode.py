from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select the structure-mapped FGFR2 barcode from FDR-significant alignment positions."
    )
    parser.add_argument("--jsd", required=True)
    parser.add_argument("--alignment-evidence", required=True)
    parser.add_argument("--structure-map", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--minimum-discriminating-score", type=float, default=0.70)
    parser.add_argument("--expected-selected", type=int, default=17)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    jsd_path = Path(args.jsd)
    evidence_path = Path(args.alignment_evidence)
    mapping_path = Path(args.structure_map)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    jsd = pd.read_csv(jsd_path, sep="\t")
    evidence = pd.read_csv(evidence_path, sep="\t")
    mapping = pd.read_csv(mapping_path, sep="\t")

    significant = jsd[jsd["significant_fdr"].astype(bool)].copy()
    audit = significant.merge(
        evidence,
        left_on="alignment_column_1based",
        right_on="alignment_col",
        how="left",
        validate="one_to_one",
        suffixes=("_jsd", "_evidence"),
    )

    audit["criterion_different_major_residues"] = (
        audit["between_isoform_difference"].eq(1)
        & audit["IIIb_major_aa_evidence"].notna()
        & audit["IIIc_major_aa_evidence"].notna()
        & audit["IIIb_major_aa_evidence"].ne("-")
        & audit["IIIc_major_aa_evidence"].ne("-")
        & audit["IIIb_major_aa_evidence"].ne(audit["IIIc_major_aa_evidence"])
    )
    audit["recomputed_discriminating_score"] = (
        audit[["IIIb_major_aa_fraction", "IIIc_major_aa_fraction"]].min(axis=1)
        * (1.0 - audit[["IIIb_gap_fraction", "IIIc_gap_fraction"]].max(axis=1))
        * audit["criterion_different_major_residues"].astype(int)
    )
    if not np.allclose(
        audit["recomputed_discriminating_score"], audit["discriminating_score"]
    ):
        raise ValueError("Stored discriminating scores do not reproduce from the recorded fractions.")
    audit["criterion_score"] = audit["recomputed_discriminating_score"].ge(
        args.minimum_discriminating_score
    )
    audit["criterion_informative"] = audit["informative_column"].fillna(False).astype(bool)
    audit["criterion_not_gap_rich"] = ~audit[
        "gap_rich_excluded_from_main_plot"
    ].fillna(True).astype(bool)
    audit["selected_for_structure"] = audit[
        [
            "criterion_different_major_residues",
            "criterion_score",
            "criterion_informative",
            "criterion_not_gap_rich",
        ]
    ].all(axis=1)

    selected = audit[audit["selected_for_structure"]].copy()
    if len(selected) != args.expected_selected:
        raise ValueError(
            f"Expected {args.expected_selected} selected positions, observed {len(selected)}."
        )

    selected_columns = selected["alignment_column_1based"].astype(int).tolist()
    mapped_columns = mapping["alignment_col"].astype(int).tolist()
    if set(selected_columns) != set(mapped_columns):
        missing = sorted(set(selected_columns) - set(mapped_columns))
        extra = sorted(set(mapped_columns) - set(selected_columns))
        raise ValueError(f"Structure-map mismatch; missing={missing}, extra={extra}")

    required_mapping = ["IIIb_PDB_1NUN_resi", "IIIc_PDB_2FDB_resi_P"]
    if mapping[required_mapping].isna().any().any():
        raise ValueError("At least one selected barcode position lacks a required PDB mapping.")

    audit.to_csv(outdir / "barcode_selection_audit.tsv", sep="\t", index=False)
    mapping.sort_values("alignment_col").to_csv(
        outdir / "structure_mapping_17_discriminating_positions.tsv",
        sep="\t",
        index=False,
    )

    summary = {
        "selection_contract": {
            "starting_set": "positions significant after Benjamini-Hochberg FDR correction",
            "different_major_residues": True,
            "minimum_discriminating_score": args.minimum_discriminating_score,
            "discriminating_score_definition": "min(IIIb major-residue fraction, IIIc major-residue fraction) * (1 - max(IIIb gap fraction, IIIc gap fraction))",
            "informative_alignment_column": True,
            "gap_rich_columns_excluded": True,
            "structural_outcomes_used_for_selection": False,
        },
        "counts": {
            "alignment_positions": int(len(jsd)),
            "fdr_significant_positions": int(len(significant)),
            "selected_structure_barcode_positions": int(len(selected)),
        },
        "selected_alignment_columns_1based": selected_columns,
        "required_structure_mappings": required_mapping,
        "input_sha256": {
            "weighted_jsd_positions": sha256(jsd_path),
            "alignment_evidence": sha256(evidence_path),
            "structure_map": sha256(mapping_path),
        },
    }
    (outdir / "barcode_selection_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
