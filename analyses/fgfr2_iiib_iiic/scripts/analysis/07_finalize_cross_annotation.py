
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


KNOWN_IG_SIGNATURES = {
    "G3DSA:2.60.40.10",
    "PS50835",
    "SSF48726",
    "SM00409",
    "PF07679",
    "SM00408",
    "PF13927",
    "SM00406",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Integrate InterProScan output with NCBI/Ensembl FGFR2 cassette coordinates and finalize cross-annotation boundary replication."
    )
    p.add_argument("--selected", required=True, help="selected_cross_annotation_candidates.tsv")
    p.add_argument("--interpro-tsv", required=True, help="InterProScan TSV output for selected proteins")
    p.add_argument("--outdir", required=True)
    p.add_argument("--distance-threshold", type=float, default=15.0)
    p.add_argument("--consensus", type=float, default=0.80)
    return p.parse_args()


def read_interpro_tsv(path: str) -> pd.DataFrame:
    rows = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError(f"InterProScan TSV line {line_number} has only {len(fields)} columns")
            fields += [""] * (15 - len(fields))
            rows.append(
                {
                    "sequence_id": fields[0],
                    "sequence_md5": fields[1],
                    "sequence_length": int(fields[2]),
                    "member_database": fields[3],
                    "signature_accession": fields[4],
                    "signature_description": fields[5],
                    "domain_start_aa": int(fields[6]),
                    "domain_end_aa": int(fields[7]),
                    "score": fields[8],
                    "status": fields[9],
                    "date": fields[10],
                    "interpro_accession": fields[11] if len(fields) > 11 else "",
                    "interpro_description": fields[12] if len(fields) > 12 else "",
                    "go_terms": fields[13] if len(fields) > 13 else "",
                    "pathways": fields[14] if len(fields) > 14 else "",
                }
            )
    return pd.DataFrame(rows)


def parse_sequence_id(sequence_id: str) -> tuple[str, str, str, str]:
    parts = sequence_id.split("|", 3)
    if len(parts) != 4:
        raise ValueError(
            f"Expected FASTA ID species|source|isoform|protein_id, received: {sequence_id}"
        )
    safe_species, source, isoform, protein_id = parts
    return safe_species, source.lower(), isoform, protein_id


def is_ig_call(row: pd.Series) -> bool:
    if str(row["signature_accession"]) in KNOWN_IG_SIGNATURES:
        return True
    text = " ".join(
        [
            str(row.get("signature_description", "")),
            str(row.get("interpro_description", "")),
        ]
    ).lower()
    return bool(re.search(r"immunoglobulin|\big[- ]?like\b|\bigc?\d?\b", text))


def main() -> None:
    args = parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    selected = pd.read_csv(args.selected, sep="\t")
    required = {
        "source",
        "species",
        "protein_id",
        "predicted_isoform",
        "cassette_start_aa",
        "cassette_end_aa",
    }
    missing = required - set(selected.columns)
    if missing:
        raise ValueError(f"Selected-candidate table is missing columns: {sorted(missing)}")

    interpro = read_interpro_tsv(args.interpro_tsv)
    parsed = interpro["sequence_id"].map(parse_sequence_id)
    interpro[["safe_species", "source", "predicted_isoform", "protein_id"]] = pd.DataFrame(
        parsed.tolist(), index=interpro.index
    )

    merged = interpro.merge(
        selected[
            [
                "source",
                "species",
                "protein_id",
                "predicted_isoform",
                "cassette_start_aa",
                "cassette_end_aa",
                "protein_length",
            ]
        ],
        on=["source", "protein_id", "predicted_isoform"],
        how="inner",
        validate="many_to_one",
    )
    if merged.empty:
        raise ValueError("No InterProScan rows matched selected candidates. Check FASTA IDs and input files.")

    merged["is_ig_call"] = merged.apply(is_ig_call, axis=1)
    merged["overlap_start"] = merged[["domain_start_aa", "cassette_start_aa"]].max(axis=1)
    merged["overlap_end"] = merged[["domain_end_aa", "cassette_end_aa"]].min(axis=1)
    merged["overlap_aa"] = (merged["overlap_end"] - merged["overlap_start"] + 1).clip(lower=0)
    calls = merged[merged["is_ig_call"] & (merged["overlap_aa"] > 0)].copy()
    if calls.empty:
        raise ValueError("No cassette-overlapping immunoglobulin-like InterPro calls were detected.")

    calls["end_signed_offset"] = calls["domain_end_aa"] - calls["cassette_end_aa"]
    calls["end_abs_offset"] = calls["end_signed_offset"].abs()
    calls["start_inside"] = (
        (calls["cassette_start_aa"] >= calls["domain_start_aa"])
        & (calls["cassette_start_aa"] <= calls["domain_end_aa"])
    )
    calls["start_distance_nearest_edge"] = np.minimum(
        (calls["cassette_start_aa"] - calls["domain_start_aa"]).abs(),
        (calls["domain_end_aa"] - calls["cassette_start_aa"]).abs(),
    )
    calls.to_csv(out / "cross_annotation_ig_overlap_calls.tsv", sep="\t", index=False)

    keys = ["species", "source", "predicted_isoform", "protein_id"]
    db = (
        calls.groupby(keys + ["member_database"], as_index=False)
        .agg(
            domain_start_aa=("domain_start_aa", "median"),
            domain_end_aa=("domain_end_aa", "median"),
            cassette_start_aa=("cassette_start_aa", "first"),
            cassette_end_aa=("cassette_end_aa", "first"),
            n_signatures=("signature_accession", "size"),
        )
    )
    db["end_signed_offset"] = db["domain_end_aa"] - db["cassette_end_aa"]
    db["end_abs_offset"] = db["end_signed_offset"].abs()
    db["end_support"] = db["end_abs_offset"] <= args.distance_threshold
    db["start_inside"] = (
        (db["cassette_start_aa"] >= db["domain_start_aa"])
        & (db["cassette_start_aa"] <= db["domain_end_aa"])
    )
    db["start_distance_nearest_edge"] = np.minimum(
        (db["cassette_start_aa"] - db["domain_start_aa"]).abs(),
        (db["domain_end_aa"] - db["cassette_start_aa"]).abs(),
    )
    db["start_support"] = db["start_inside"] & (
        db["start_distance_nearest_edge"] > args.distance_threshold
    )
    db.to_csv(out / "cross_annotation_member_database_calls.tsv", sep="\t", index=False)

    per_protein = (
        db.groupby(keys, as_index=False)
        .agg(
            n_member_databases=("member_database", "nunique"),
            cassette_start_aa=("cassette_start_aa", "first"),
            cassette_end_aa=("cassette_end_aa", "first"),
            median_domain_end_aa=("domain_end_aa", "median"),
            median_end_signed_offset=("end_signed_offset", "median"),
            median_end_abs_offset=("end_abs_offset", "median"),
            start_consensus=("start_support", "mean"),
            end_consensus=("end_support", "mean"),
        )
    )
    per_protein["start_pass"] = per_protein["start_consensus"] >= args.consensus
    per_protein["end_pass"] = per_protein["end_consensus"] >= args.consensus
    per_protein["topology_pass"] = per_protein["start_pass"] & per_protein["end_pass"]
    per_protein.to_csv(out / "cross_annotation_topology_by_protein.tsv", sep="\t", index=False)

    comparisons = []
    for (species, isoform), group in per_protein.groupby(["species", "predicted_isoform"]):
        if set(group["source"]) >= {"ncbi", "ensembl"}:
            n = group[group["source"] == "ncbi"].iloc[0]
            e = group[group["source"] == "ensembl"].iloc[0]
            comparisons.append(
                {
                    "species": species,
                    "isoform": isoform,
                    "ncbi_protein_id": n["protein_id"],
                    "ensembl_protein_id": e["protein_id"],
                    "cassette_start_delta_ncbi_minus_ensembl": n["cassette_start_aa"] - e["cassette_start_aa"],
                    "cassette_end_delta_ncbi_minus_ensembl": n["cassette_end_aa"] - e["cassette_end_aa"],
                    "median_d3_end_delta_ncbi_minus_ensembl": n["median_domain_end_aa"] - e["median_domain_end_aa"],
                    "same_topology_class": bool(n["topology_pass"] == e["topology_pass"]),
                    "ncbi_topology_pass": bool(n["topology_pass"]),
                    "ensembl_topology_pass": bool(e["topology_pass"]),
                    "minimum_start_consensus": min(n["start_consensus"], e["start_consensus"]),
                    "minimum_end_consensus": min(n["end_consensus"], e["end_consensus"]),
                }
            )
    comparison = pd.DataFrame(comparisons)
    comparison.to_csv(out / "cross_annotation_boundary_replication.tsv", sep="\t", index=False)

    if not comparison.empty:
        plot = comparison.sort_values(["species", "isoform"]).copy()
        labels = plot["species"].str.replace("_", " ") + " " + plot["isoform"]
        fig, ax = plt.subplots(figsize=(10, max(4, 0.35 * len(plot))))
        y = np.arange(len(plot))
        ax.scatter(plot["cassette_end_delta_ncbi_minus_ensembl"], y, label="Cassette end")
        ax.scatter(plot["median_d3_end_delta_ncbi_minus_ensembl"], y, marker="x", label="Median D3 end")
        ax.axvline(0, linewidth=1)
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_xlabel("NCBI minus Ensembl coordinate (aa)")
        ax.set_title("Cross-annotation coordinate differences with topology replication")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out / "figure_cross_annotation_replication.png", dpi=300)
        fig.savefig(out / "figure_cross_annotation_replication.svg")
        plt.close(fig)

    summary = {
        "distance_threshold_aa": args.distance_threshold,
        "consensus_threshold": args.consensus,
        "n_selected_proteins_with_ig_calls": int(len(per_protein)),
        "n_proteins_passing_topology": int(per_protein["topology_pass"].sum()),
        "n_ncbi_ensembl_pairs": int(len(comparison)),
        "n_pairs_same_topology_class": int(comparison["same_topology_class"].sum()) if not comparison.empty else 0,
        "cross_annotation_topology_agreement": float(comparison["same_topology_class"].mean()) if not comparison.empty else None,
    }
    (out / "cross_annotation_boundary_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
