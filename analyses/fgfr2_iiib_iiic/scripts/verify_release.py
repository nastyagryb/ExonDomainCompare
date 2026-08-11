from __future__ import annotations

import csv
import io
import json
import re
import tokenize
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_json(path: str) -> dict:
    return json.loads((ROOT / path).read_text())


def read_tsv(path: str) -> list[dict[str, str]]:
    with (ROOT / path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def close(value: float, expected: float, tolerance: float = 1e-10) -> bool:
    return abs(value - expected) <= tolerance


def main() -> None:
    required = [
        "data/framework_snapshot/cohort/fgfr2_post_rescue_final_truth_table.tsv",
        "data/framework_snapshot/alignments/final_fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa",
        "data/framework_snapshot/exon_protein_mapping/figure3C_exon_to_protein_cassette_coordinate_map.tsv",
        "data/framework_snapshot/protein_annotations/interpro_ensemble_coordinate_support_calls.tsv",
        "results/00_framework_outputs/website_figures/cmp_primary_msa_overview.svg",
        "results/00_framework_outputs/website_figures/cmp_boundary_matrix.tsv",
        "results/05_iqtree/parsed_AU_topology_test.tsv",
        "results/08_loco_asr_synthesis/analysis_summary.json",
    ]
    missing = [path for path in required if not (ROOT / path).is_file()]
    if missing:
        raise SystemExit("Missing required files: " + ", ".join(missing))

    truth = read_tsv(required[0])
    species = {row["species"] for row in truth}
    labels = {label: sum(row["final_isoform_label"] == label for row in truth) for label in ("IIIb", "IIIc")}
    if len(truth) != 60 or len(species) != 30 or labels != {"IIIb": 30, "IIIc": 30}:
        raise SystemExit("The final cohort does not match the frozen 30-species design")

    caller = read_json("results/01_domain_caller_bias/domain_caller_bias_summary.json")
    robust = read_json("results/02_robustness_surface/robustness_surface_summary.json")
    jsd = read_json("results/03_weighted_jsd/weighted_jsd_summary.json")
    cross = read_json("results/07_cross_annotation_final/cross_annotation_boundary_summary.json")
    final = read_json("results/08_loco_asr_synthesis/analysis_summary.json")

    assertions = [
        caller["n_proteins"] == 58,
        caller["n_member_databases"] == 5,
        close(caller["database_partial_r2_controlling_species_and_isoform"], 0.6221791428557758),
        robust["n_grid_points"] == 1326,
        close(robust["minimum_distance_for_100pct_at_80pct_consensus"], 12.0),
        jsd["n_species_pairs"] == 28,
        jsd["n_fdr_significant_positions"] == 25,
        close(jsd["global_paired_permutation_p"], 9.999000099990002e-05),
        cross["n_ncbi_ensembl_pairs"] == 16,
        cross["n_pairs_same_topology_class"] == 16,
        final["loco"]["correct_sequences"] == 56,
        final["loco"]["total_sequences"] == 56,
        final["asr"]["modern_barcode_sites_ancestrally_different"] == 15,
        final["asr"]["high_confidence_all_jackknife_stable_ancestral_core"] == 11,
    ]
    if not all(assertions):
        raise SystemExit("One or more frozen result values failed validation")

    au = read_tsv("results/05_iqtree/parsed_AU_topology_test.tsv")
    if len(au) != 3:
        raise SystemExit("The AU table must contain exactly three tested topologies")
    expected = [
        (1, -583.6588558, 0.0, 1.0),
        (2, -583.6591338, 0.00027794, 1.0),
        (3, -1497.465063, 913.81, 6.36e-08),
    ]
    for row, values in zip(au, expected):
        observed = (int(row["tree"]), float(row["log_likelihood"]), float(row["delta_log_likelihood"]), float(row["p_AU"]))
        if observed[0] != values[0] or any(not close(a, b, 1e-8) for a, b in zip(observed[1:], values[1:])):
            raise SystemExit("The parsed AU values do not match the IQ-TREE report")

    forbidden = []
    comments = []
    local_prefixes = ("/" + "Users" + "/", "/mnt" + "/data")
    for path in sorted((ROOT / "scripts").rglob("*.py")):
        text = path.read_text()
        if any(prefix in text for prefix in local_prefixes):
            forbidden.append(str(path.relative_to(ROOT)))
        for token in tokenize.generate_tokens(io.StringIO(text).readline):
            if token.type == tokenize.COMMENT:
                comments.append(f"{path.relative_to(ROOT)}:{token.start[0]}")
    if forbidden:
        raise SystemExit("Local absolute paths found in scripts: " + ", ".join(forbidden))
    if comments:
        raise SystemExit("Comments found in published scripts: " + ", ".join(comments))

    unwanted = []
    too_large = []
    for path in ROOT.rglob("*"):
        if path.is_dir():
            if path.name in {".venv", "__pycache__"}:
                unwanted.append(str(path.relative_to(ROOT)))
            continue
        if path.name == ".DS_Store" or path.suffix == ".pyc":
            unwanted.append(str(path.relative_to(ROOT)))
        if path.stat().st_size >= 100_000_000:
            too_large.append(str(path.relative_to(ROOT)))
    if unwanted or too_large:
        raise SystemExit("Release hygiene check failed: " + ", ".join(unwanted + too_large))

    report = {
        "status": "passed",
        "files_checked": sum(path.is_file() for path in ROOT.rglob("*")),
        "cohort_rows": len(truth),
        "species": len(species),
        "au_topologies": len(au),
        "published_python_comments": 0,
    }
    (ROOT / "checks/verification_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
