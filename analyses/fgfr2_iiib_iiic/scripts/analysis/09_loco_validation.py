from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO
from scipy.special import expit

from publication_style import COLORS, apply_publication_style, clean_axis, save_figure, title_block


CLADES = {
    "Mammals": {
        "bos_taurus", "callithrix_jacchus", "equus_caballus", "felis_catus",
        "gorilla_gorilla_gorilla", "homo_sapiens", "macaca_mulatta",
        "monodelphis_domestica", "mus_musculus", "ornithorhynchus_anatinus",
        "oryctolagus_cuniculus", "ovis_aries", "pan_troglodytes",
        "rattus_norvegicus", "sus_scrofa",
    },
    "Birds": {"gallus_gallus", "meleagris_gallopavo", "taeniopygia_guttata"},
    "Reptiles": {"alligator_mississippiensis", "anolis_carolinensis", "chrysemys_picta_bellii"},
    "Amphibians": {"ambystoma_mexicanum", "xenopus_tropicalis"},
    "Teleosts": {"danio_rerio", "gasterosteus_aculeatus", "oreochromis_niloticus", "oryzias_latipes", "takifugu_rubripes"},
}

EXCLUDED_REVIEW_SPECIES = {"canis_lupus_familiaris", "pongo_abelii"}
CANONICAL_BARCODE_1BASED = np.array([1, 2, 4, 6, 7, 8, 11, 16, 21, 28, 30, 31, 32, 46, 50, 51, 52])


def parse_args() -> argparse.Namespace:
    base = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Leave-one-clade-out portability analysis for the FGFR2 IIIb/IIIc sequence partition."
    )
    parser.add_argument(
        "--alignment",
        type=Path,
        default=base / "data/framework_snapshot/alignments/final_fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa",
    )
    parser.add_argument(
        "--outdir", type=Path, default=base / "results/reproduced/08_loco_validation"
    )
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--bootstraps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260804)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_alignment(path: Path):
    alignment = AlignIO.read(path, "fasta")
    species = sorted(set().union(*CLADES.values()))
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY-")
    state = {aa: i for i, aa in enumerate(amino_acids)}
    array = np.full(
        (len(species), 2, alignment.get_alignment_length()), state["-"], dtype=np.int16
    )
    present = set()
    species_index = {name: i for i, name in enumerate(species)}
    for record in alignment:
        species_name, isoform, *_ = record.id.split("|")
        if species_name in EXCLUDED_REVIEW_SPECIES or species_name not in species_index:
            continue
        isoform_index = 0 if isoform == "IIIb" else 1
        array[species_index[species_name], isoform_index] = [state[x] for x in str(record.seq)]
        present.add((species_name, isoform))
    expected = {(name, isoform) for name in species for isoform in ("IIIb", "IIIc")}
    if present != expected:
        raise ValueError(f"Incomplete LOCO cohort: missing={sorted(expected - present)}")
    return species, array, amino_acids


def discover_positions(training: np.ndarray, threshold: float = 0.70):
    positions = []
    details = []
    n_species = training.shape[0]
    for column in range(training.shape[2]):
        isoform_statistics = []
        for isoform_index in (0, 1):
            counts = np.bincount(training[:, isoform_index, column], minlength=21)
            gap_fraction = counts[20] / n_species
            non_gap = n_species - counts[20]
            if non_gap == 0:
                major_index, major_fraction = 20, 0.0
            else:
                major_index = int(np.argmax(counts[:20]))
                major_fraction = counts[major_index] / non_gap
            isoform_statistics.append((major_index, major_fraction, gap_fraction))
        score = min(isoform_statistics[0][1], isoform_statistics[1][1]) * (
            1.0 - max(isoform_statistics[0][2], isoform_statistics[1][2])
        )
        selected = isoform_statistics[0][0] != isoform_statistics[1][0] and score >= threshold
        if selected:
            positions.append(column)
        details.append((column, *isoform_statistics[0], *isoform_statistics[1], score, selected))
    return np.asarray(positions, dtype=int), details


def fit_predict(training: np.ndarray, test: np.ndarray, positions: np.ndarray, alpha: float = 0.5):
    if len(positions) == 0:
        return np.zeros((test.shape[0], 2), dtype=float)
    log_probability = np.zeros((2, len(positions), 21), dtype=float)
    for isoform_index in (0, 1):
        for site_index, column in enumerate(positions):
            counts = np.bincount(training[:, isoform_index, column], minlength=21)
            log_probability[isoform_index, site_index] = np.log(
                (counts + alpha) / (counts.sum() + alpha * 21)
            )
    log_likelihood_ratio = np.zeros((test.shape[0], 2), dtype=float)
    indices = np.arange(len(positions))
    for species_index in range(test.shape[0]):
        for isoform_index in (0, 1):
            states = test[species_index, isoform_index, positions]
            log_likelihood_ratio[species_index, isoform_index] = (
                log_probability[0, indices, states].sum()
                - log_probability[1, indices, states].sum()
            )
    return log_likelihood_ratio


def accuracy_for_partition(array: np.ndarray, clade_indices: dict[str, np.ndarray]) -> float:
    correct = 0.0
    total = 0
    all_indices = np.arange(array.shape[0])
    for held_out in clade_indices.values():
        training = np.setdiff1d(all_indices, held_out)
        positions, _ = discover_positions(array[training])
        scores = fit_predict(array[training], array[held_out], positions)
        correct += np.sum(scores[:, 0] > 0) + 0.5 * np.sum(scores[:, 0] == 0)
        correct += np.sum(scores[:, 1] < 0) + 0.5 * np.sum(scores[:, 1] == 0)
        total += 2 * len(held_out)
    return correct / total


def run_loco(species, array, amino_acids, n_permutations, n_bootstraps, seed):
    species_index = {name: i for i, name in enumerate(species)}
    clade_indices = {
        clade: np.array([species_index[name] for name in sorted(members)], dtype=int)
        for clade, members in CLADES.items()
    }
    canonical = set(CANONICAL_BARCODE_1BASED - 1)
    all_indices = np.arange(len(species))
    fold_rows = []
    prediction_rows = []
    site_rows = []

    for clade, held_out in clade_indices.items():
        training = np.setdiff1d(all_indices, held_out)
        positions, details = discover_positions(array[training])
        scores = fit_predict(array[training], array[held_out], positions)
        margins = []
        for test_index, held_species_index in enumerate(held_out):
            for isoform_index, isoform in enumerate(("IIIb", "IIIc")):
                margin = scores[test_index, isoform_index] if isoform_index == 0 else -scores[test_index, isoform_index]
                margins.append(margin)
                prediction_rows.append({
                    "held_out_clade": clade,
                    "species": species[held_species_index],
                    "true_isoform": isoform,
                    "predicted_isoform": "IIIb" if scores[test_index, isoform_index] > 0 else "IIIc",
                    "correct": bool(margin > 0),
                    "signed_log_likelihood_margin": float(margin),
                    "posterior_for_true_isoform": float(expit(margin)),
                    "n_training_species": len(training),
                    "n_selected_positions": len(positions),
                })
        selected = set(positions)
        fold_rows.append({
            "held_out_clade": clade,
            "n_training_species": len(training),
            "n_test_species": len(held_out),
            "n_test_sequences": 2 * len(held_out),
            "n_de_novo_positions": len(positions),
            "canonical_positions_recovered": len(selected & canonical),
            "canonical_positions_total": len(canonical),
            "canonical_recovery_fraction": len(selected & canonical) / len(canonical),
            "additional_positions": ",".join(str(x + 1) for x in sorted(selected - canonical)),
            "canonical_positions_missed": ",".join(str(x + 1) for x in sorted(canonical - selected)),
            "accuracy": float(np.mean(np.asarray(margins) > 0)),
            "minimum_signed_margin": float(np.min(margins)),
            "median_signed_margin": float(np.median(margins)),
        })
        for detail in details:
            column, b_index, b_fraction, b_gap, c_index, c_fraction, c_gap, score, selected_flag = detail
            site_rows.append({
                "held_out_clade": clade,
                "alignment_column": column + 1,
                "IIIb_major": amino_acids[b_index],
                "IIIc_major": amino_acids[c_index],
                "IIIb_major_fraction": b_fraction,
                "IIIc_major_fraction": c_fraction,
                "IIIb_gap_fraction": b_gap,
                "IIIc_gap_fraction": c_gap,
                "discriminating_score": score,
                "selected": selected_flag,
                "canonical_barcode": column in canonical,
            })

    fold = pd.DataFrame(fold_rows)
    prediction = pd.DataFrame(prediction_rows)
    site = pd.DataFrame(site_rows)
    observed_accuracy = float(prediction["correct"].mean())

    rng = np.random.default_rng(seed)
    null = np.empty(n_permutations, dtype=float)
    for permutation_index in range(n_permutations):
        swaps = rng.integers(0, 2, size=len(species)).astype(bool)
        permuted = array.copy()
        permuted[swaps] = permuted[swaps][:, ::-1].copy()
        null[permutation_index] = accuracy_for_partition(permuted, clade_indices)
    permutation_p = float((np.sum(null >= observed_accuracy) + 1) / (n_permutations + 1))

    rng = np.random.default_rng(seed + 1)
    bootstrap_fold_rows = []
    bootstrap_site_rows = []
    for clade, held_out in clade_indices.items():
        training = np.setdiff1d(all_indices, held_out)
        accuracies = []
        counts = np.zeros(array.shape[2], dtype=int)
        selected_counts = []
        all_correct = []
        for _ in range(n_bootstraps):
            sampled = rng.choice(training, size=len(training), replace=True)
            positions, _ = discover_positions(array[sampled])
            counts[positions] += 1
            scores = fit_predict(array[sampled], array[held_out], positions)
            correct = np.concatenate((scores[:, 0] > 0, scores[:, 1] < 0))
            accuracies.append(float(correct.mean()))
            selected_counts.append(len(positions))
            all_correct.append(bool(correct.all()))
        bootstrap_fold_rows.append({
            "held_out_clade": clade,
            "n_bootstraps": n_bootstraps,
            "median_accuracy": float(np.median(accuracies)),
            "accuracy_2.5pct": float(np.quantile(accuracies, 0.025)),
            "accuracy_97.5pct": float(np.quantile(accuracies, 0.975)),
            "fraction_bootstraps_with_100pct_accuracy": float(np.mean(all_correct)),
            "median_selected_positions": float(np.median(selected_counts)),
            "selected_positions_2.5pct": float(np.quantile(selected_counts, 0.025)),
            "selected_positions_97.5pct": float(np.quantile(selected_counts, 0.975)),
        })
        for column, count in enumerate(counts):
            bootstrap_site_rows.append({
                "held_out_clade": clade,
                "alignment_column": column + 1,
                "selection_frequency": count / n_bootstraps,
                "canonical_barcode": column in canonical,
            })

    return {
        "fold": fold,
        "prediction": prediction,
        "site": site,
        "permutation": pd.DataFrame({"permutation_index": np.arange(1, n_permutations + 1), "overall_accuracy": null}),
        "permutation_p": permutation_p,
        "bootstrap_fold": pd.DataFrame(bootstrap_fold_rows),
        "bootstrap_site": pd.DataFrame(bootstrap_site_rows),
    }


def plot_loco(result, outdir: Path) -> None:
    apply_publication_style()
    fold = result["fold"].set_index("held_out_clade").loc[list(CLADES)].reset_index()
    prediction = result["prediction"]
    null = result["permutation"]["overall_accuracy"].to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 6.8), gridspec_kw={"width_ratios": [1.0, 1.25]})
    fig.subplots_adjust(top=0.77, left=0.10, right=0.96, bottom=0.16, wspace=0.32)

    y = np.arange(len(fold))
    axes[0].barh(y, 100 * fold["accuracy"], color=COLORS["teal"])
    axes[0].set_yticks(y, fold["held_out_clade"], fontweight="bold")
    axes[0].invert_yaxis()
    axes[0].set_xlim(0, 112)
    axes[0].set_xlabel("Held-out accuracy (%)")
    for index, row in fold.iterrows():
        axes[0].text(101, index, f"{int(row.n_test_sequences)}/{int(row.n_test_sequences)}", va="center", fontweight="bold")
    clean_axis(axes[0], "x")

    axes[1].hist(null, bins=np.linspace(0.38, 0.85, 28), color=COLORS["sky"], edgecolor="white")
    axes[1].axvline(1.0, color=COLORS["red"], linewidth=2.5)
    axes[1].set_xlim(0.38, 1.03)
    axes[1].set_xlabel("LOCO accuracy after paired label swaps")
    axes[1].set_ylabel("Permutation count")
    clean_axis(axes[1], "y")
    axes[1].text(
        0.03, 0.94,
        f"Observed = {prediction.correct.mean():.2f}\nNull maximum = {null.max():.3f}\nEmpirical p = {result['permutation_p']:.4g}",
        transform=axes[1].transAxes, va="top", fontweight="bold",
    )
    title_block(
        fig,
        "The sequence-calibrated IIIb/IIIc partition transfers across held-out vertebrate clades",
        "Candidate sites are re-selected from training clades only; this tests portability within the analysed cohort, not independent truth of the isoform labels.",
    )
    save_figure(fig, outdir, "leave_one_clade_out_validation")


def main() -> None:
    args = parse_args()
    if args.outdir.exists():
        shutil.rmtree(args.outdir)
    (args.outdir / "tables").mkdir(parents=True)
    (args.outdir / "figures").mkdir()

    species, array, amino_acids = load_alignment(args.alignment)
    result = run_loco(
        species, array, amino_acids, args.permutations, args.bootstraps, args.seed
    )
    result["fold"].to_csv(args.outdir / "tables/loco_fold_summary.tsv", sep="\t", index=False)
    result["prediction"].to_csv(args.outdir / "tables/loco_held_out_predictions.tsv", sep="\t", index=False)
    result["site"].to_csv(args.outdir / "tables/loco_de_novo_site_selection.tsv", sep="\t", index=False)
    result["permutation"].to_csv(args.outdir / "tables/loco_paired_label_permutation_null.tsv", sep="\t", index=False)
    result["bootstrap_fold"].to_csv(args.outdir / "tables/loco_training_bootstrap_summary.tsv", sep="\t", index=False)
    result["bootstrap_site"].to_csv(args.outdir / "tables/loco_training_bootstrap_site_stability.tsv", sep="\t", index=False)
    plot_loco(result, args.outdir / "figures")

    summary = {
        "analysis": "leave-one-clade-out portability of the sequence-calibrated IIIb/IIIc partition",
        "interpretive_limit": "The labels and feature definition originate within the analysed FGFR2 cohort; the result is cross-clade portability, not independent validation of isoform truth.",
        "n_species_pairs": len(species),
        "n_held_out_proteins": int(len(result["prediction"])),
        "n_correct": int(result["prediction"]["correct"].sum()),
        "observed_accuracy": float(result["prediction"]["correct"].mean()),
        "paired_label_permutations": args.permutations,
        "paired_label_permutation_p": result["permutation_p"],
        "training_bootstraps_per_fold": args.bootstraps,
        "all_bootstraps_perfect_in_every_fold": bool(
            result["bootstrap_fold"]["fraction_bootstraps_with_100pct_accuracy"].eq(1.0).all()
        ),
        "seed": args.seed,
        "bootstrap_seed": args.seed + 1,
        "alignment_sha256": sha256(args.alignment),
        "ancestral_reconstruction_performed": False,
    }
    (args.outdir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    (args.outdir / "README.md").write_text(
        "# Leave-one-clade-out validation\n\n"
        "This directory contains the valid LOCO portability analysis only. Candidate positions are rediscovered from each training partition before classifying the excluded clade. The result must not be described as independent validation of the source isoform labels. An earlier ancestral reconstruction was separated and quarantined because the target ancestral split was not identifiable from the unrooted topology used.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
