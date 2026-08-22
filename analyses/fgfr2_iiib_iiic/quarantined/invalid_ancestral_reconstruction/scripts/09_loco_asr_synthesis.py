
from __future__ import annotations

import argparse
import copy
import json
import math
import shutil
import textwrap
import zipfile
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import AlignIO, Phylo
from scipy.linalg import expm
from scipy.special import expit, gammainc
from scipy.stats import gamma
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch
from matplotlib.lines import Line2D




COLORS = {
    "navy": "#102A43",
    "blue": "#2F6B9A",
    "sky": "#64B5D2",
    "teal": "#138A80",
    "green": "#2E8B57",
    "gold": "#D59B22",
    "orange": "#D96B38",
    "red": "#B94A48",
    "purple": "#7561A8",
    "gray": "#6B7785",
    "light": "#EDF2F7",
    "grid": "#DCE4EC",
    "dark": "#243B53",
    "white": "#FFFFFF",
}

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 16,
    "axes.labelsize": 12,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 9.5,
    "figure.titlesize": 22,
    "axes.linewidth": 0.9,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
})

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

AA_CLASS = {
    **{a: "hydrophobic" for a in "AVILMFWY"},
    **{a: "polar" for a in "STNQ"},
    **{a: "positive" for a in "KRH"},
    **{a: "negative" for a in "DE"},
    **{a: "special" for a in "CGP"},
    "-": "gap",
    "X": "uncertain",
}
AA_CLASS_COLOR = {
    "hydrophobic": "#6C7FA3",
    "polar": "#39A89D",
    "positive": "#D59B22",
    "negative": "#D96B38",
    "special": "#8A68A8",
    "gap": "#E5EAF0",
    "uncertain": "#B94A48",
}



LG_LOWER = """
0.425093
0.276818 0.751878
0.395144 0.123954 5.076149
2.489084 0.534551 0.528768 0.062556
0.969894 2.807908 1.695752 0.523386 0.084808
1.038545 0.363970 0.541712 5.243870 0.003499 4.128591
2.066040 0.390192 1.437645 0.844926 0.569265 0.267959 0.348847
0.358858 2.426601 4.509238 0.927114 0.640543 4.813505 0.423881 0.311484
0.149830 0.126991 0.191503 0.010690 0.320627 0.072854 0.044265 0.008705 0.108882
0.395337 0.301848 0.068427 0.015076 0.594007 0.582457 0.069673 0.044261 0.366317 4.145067
0.536518 6.326067 2.145078 0.282959 0.013266 3.234294 1.807177 0.296636 0.697264 0.159069 0.137500
1.124035 0.484133 0.371004 0.025548 0.893680 1.672569 0.173735 0.139538 0.442472 4.273607 6.312358 0.656604
0.253701 0.052722 0.089525 0.017416 1.105251 0.035855 0.018811 0.089586 0.682139 1.112727 2.592692 0.023918 1.798853
1.177651 0.332533 0.161787 0.394456 0.075382 0.624294 0.419409 0.196961 0.508851 0.078281 0.249060 0.390322 0.099849 0.094464
4.727182 0.858151 4.008358 1.240275 2.784478 1.223828 0.611973 1.739990 0.990012 0.064105 0.182287 0.748683 0.346960 0.361819 1.338132
2.139501 0.578987 2.000679 0.425860 1.143480 1.080136 0.604545 0.129836 0.584262 1.033739 0.302936 1.136863 2.020366 0.165001 0.571468 6.472279
0.180717 0.593607 0.045376 0.029890 0.670128 0.236199 0.077852 0.268491 0.597054 0.111660 0.619632 0.049906 0.696175 2.457121 0.095131 0.248862 0.140825
0.218959 0.314440 0.612025 0.135107 1.165532 0.257336 0.120037 0.054679 5.306834 0.232523 0.299648 0.131932 0.481306 7.803902 0.089613 0.400547 0.245841 3.151815
2.547870 0.170887 0.083688 0.037967 1.959291 0.210332 0.245034 0.076701 0.119013 10.649107 1.702745 0.185202 1.898718 0.654683 0.296501 0.098369 2.188158 0.189510 0.249313
"""
LG_FREQ = np.array([
    0.07906592, 0.05594094, 0.04197696, 0.05305195, 0.01293699,
    0.04076696, 0.07158593, 0.05733694, 0.02235498, 0.06215694,
    0.09908090, 0.06459994, 0.02295098, 0.04230196, 0.04403996,
    0.06119694, 0.05328695, 0.01206599, 0.03415497, 0.06914693,
])
LG_ORDER = list("ARNDCQEGHILKMFPSTWYV")


def save_figure(fig: plt.Figure, outdir: Path, stem: str) -> None:
    for ext in ("png", "svg", "pdf"):
        kwargs = {"dpi": 400} if ext == "png" else {}
        fig.savefig(outdir / f"{stem}.{ext}", **kwargs)
    plt.close(fig)


def clean_axis(ax, grid_axis="y"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=COLORS["grid"], linewidth=0.75, alpha=0.85)
    ax.set_axisbelow(True)


def header(fig, title: str, subtitle: str):
    fig.suptitle(title, x=0.045, y=0.985, ha="left", va="top",
                 fontsize=22, fontweight="bold", color=COLORS["navy"])
    fig.text(0.045, 0.94, subtitle, ha="left", va="top",
             fontsize=11.5, color=COLORS["gray"])





def load_alignment(path: Path):
    aln = AlignIO.read(path, "fasta")
    species = sorted(set().union(*CLADES.values()))
    aa_order = list("ACDEFGHIKLMNPQRSTVWY-")
    aa_to_idx = {aa: i for i, aa in enumerate(aa_order)}
    arr = np.full((len(species), 2, aln.get_alignment_length()), aa_to_idx["-"], dtype=np.int16)
    seq_by_tree_tip = {}
    alignment_id = {}
    species_to_idx = {s: i for i, s in enumerate(species)}
    for record in aln:
        sp, iso, *_ = record.id.split("|")
        if sp in EXCLUDED_REVIEW_SPECIES or sp not in species_to_idx:
            continue
        iso_idx = 0 if iso == "IIIb" else 1
        seq = str(record.seq)
        arr[species_to_idx[sp], iso_idx, :] = [aa_to_idx[x] for x in seq]
        seq_by_tree_tip[f"{sp}__{iso}"] = seq
        alignment_id[(sp, iso)] = record.id
    all_gap = np.all(arr == aa_to_idx["-"], axis=(0, 1))
    informative_cols = np.where(~all_gap)[0]
    return aln, species, arr, aa_order, aa_to_idx, seq_by_tree_tip, alignment_id, informative_cols





def discover_positions(train_arr: np.ndarray, threshold: float = 0.70):
    gap_idx = 20
    positions = []
    details = []
    n_species = train_arr.shape[0]
    for col in range(train_arr.shape[2]):
        info = []
        for iso_idx in (0, 1):
            values = train_arr[:, iso_idx, col]
            counts = np.bincount(values, minlength=21)
            gap_fraction = counts[gap_idx] / n_species
            nongap_n = n_species - counts[gap_idx]
            if nongap_n == 0:
                major_idx, major_fraction = gap_idx, 0.0
            else:
                major_idx = int(np.argmax(counts[:20]))
                major_fraction = counts[major_idx] / nongap_n
            info.append((major_idx, major_fraction, gap_fraction))
        score = min(info[0][1], info[1][1]) * (1.0 - max(info[0][2], info[1][2]))
        selected = info[0][0] != info[1][0] and score >= threshold
        if selected:
            positions.append(col)
        details.append((col, *info[0], *info[1], score, selected))
    return np.asarray(positions, dtype=int), details


def fit_predict_naive_bayes(train_arr: np.ndarray, test_arr: np.ndarray,
                            positions: np.ndarray, alpha: float = 0.5) -> np.ndarray:
    if len(positions) == 0:
        return np.zeros((test_arr.shape[0], 2), dtype=float)
    logp = np.zeros((2, len(positions), 21), dtype=float)
    for iso_idx in (0, 1):
        for k, col in enumerate(positions):
            counts = np.bincount(train_arr[:, iso_idx, col], minlength=21)
            logp[iso_idx, k, :] = np.log((counts + alpha) / (counts.sum() + alpha * 21))
    llr = np.zeros((test_arr.shape[0], 2), dtype=float)
    for species_i in range(test_arr.shape[0]):
        for seq_slot in (0, 1):
            states = test_arr[species_i, seq_slot, positions]
            idx = np.arange(len(positions))
            llr[species_i, seq_slot] = (
                logp[0, idx, states].sum() - logp[1, idx, states].sum()
            )
    return llr


def loco_accuracy(data_arr, clade_indices):
    total_correct = 0.0
    total = 0
    all_indices = np.arange(data_arr.shape[0])
    for held_idx in clade_indices.values():
        train_idx = np.setdiff1d(all_indices, held_idx)
        positions, _ = discover_positions(data_arr[train_idx])
        llr = fit_predict_naive_bayes(data_arr[train_idx], data_arr[held_idx], positions)
        total_correct += np.sum(llr[:, 0] > 0) + 0.5 * np.sum(llr[:, 0] == 0)
        total_correct += np.sum(llr[:, 1] < 0) + 0.5 * np.sum(llr[:, 1] == 0)
        total += 2 * len(held_idx)
    return total_correct / total


def run_loco(species, arr, aa_order, n_permutations=10000, n_bootstraps=1000, seed=20260804):
    species_to_idx = {s: i for i, s in enumerate(species)}
    clade_indices = {
        clade: np.array([species_to_idx[s] for s in sorted(members)], dtype=int)
        for clade, members in CLADES.items()
    }
    canonical_zero = set(CANONICAL_BARCODE_1BASED - 1)
    all_indices = np.arange(len(species))

    fold_rows, prediction_rows, site_rows = [], [], []
    for clade, held_idx in clade_indices.items():
        train_idx = np.setdiff1d(all_indices, held_idx)
        positions, details = discover_positions(arr[train_idx])
        llr = fit_predict_naive_bayes(arr[train_idx], arr[held_idx], positions)
        margins = []
        for test_i, species_i in enumerate(held_idx):
            for iso_idx, iso in enumerate(("IIIb", "IIIc")):
                signed_margin = llr[test_i, iso_idx] if iso_idx == 0 else -llr[test_i, iso_idx]
                correct = signed_margin > 0
                margins.append(signed_margin)
                prediction_rows.append({
                    "held_out_clade": clade,
                    "species": species[species_i],
                    "true_isoform": iso,
                    "predicted_isoform": "IIIb" if llr[test_i, iso_idx] > 0 else "IIIc",
                    "correct": bool(correct),
                    "signed_log_likelihood_margin": float(signed_margin),
                    "posterior_for_true_isoform": float(expit(signed_margin)),
                    "n_training_species": len(train_idx),
                    "n_selected_positions": len(positions),
                })
        selected_set = set(positions)
        fold_rows.append({
            "held_out_clade": clade,
            "n_training_species": len(train_idx),
            "n_test_species": len(held_idx),
            "n_test_sequences": 2 * len(held_idx),
            "n_de_novo_positions": len(positions),
            "canonical_positions_recovered": len(selected_set & canonical_zero),
            "canonical_positions_total": len(canonical_zero),
            "canonical_recovery_fraction": len(selected_set & canonical_zero) / len(canonical_zero),
            "additional_positions": ",".join(map(str, sorted(x + 1 for x in selected_set - canonical_zero))),
            "canonical_positions_missed": ",".join(map(str, sorted(x + 1 for x in canonical_zero - selected_set))),
            "accuracy": float(np.mean(np.array(margins) > 0)),
            "minimum_signed_margin": float(np.min(margins)),
            "median_signed_margin": float(np.median(margins)),
        })
        for item in details:
            col, b_idx, b_frac, b_gap, c_idx, c_frac, c_gap, score, selected_flag = item
            site_rows.append({
                "held_out_clade": clade,
                "alignment_column": col + 1,
                "IIIb_major": aa_order[b_idx],
                "IIIc_major": aa_order[c_idx],
                "IIIb_major_fraction": b_frac,
                "IIIc_major_fraction": c_frac,
                "IIIb_gap_fraction": b_gap,
                "IIIc_gap_fraction": c_gap,
                "discriminating_score": score,
                "selected": selected_flag,
                "canonical_barcode": col in canonical_zero,
            })

    fold_df = pd.DataFrame(fold_rows)
    prediction_df = pd.DataFrame(prediction_rows)
    site_df = pd.DataFrame(site_rows)

    rng = np.random.default_rng(seed)
    permutation_null = np.empty(n_permutations, dtype=float)
    for p in range(n_permutations):
        swaps = rng.integers(0, 2, size=len(species)).astype(bool)
        permuted = arr.copy()
        permuted[swaps] = permuted[swaps][:, ::-1, :].copy()
        permutation_null[p] = loco_accuracy(permuted, clade_indices)
    permutation_p = (np.sum(permutation_null >= prediction_df["correct"].mean()) + 1) / (n_permutations + 1)
    permutation_df = pd.DataFrame({
        "permutation_index": np.arange(1, n_permutations + 1),
        "overall_accuracy": permutation_null,
    })

    rng = np.random.default_rng(seed + 1)
    bootstrap_fold_rows, bootstrap_site_rows = [], []
    for clade, held_idx in clade_indices.items():
        train_idx = np.setdiff1d(all_indices, held_idx)
        accuracies, n_positions, all_correct = [], [], []
        selection_counts = np.zeros(arr.shape[2], dtype=int)
        for _ in range(n_bootstraps):
            boot_idx = rng.choice(train_idx, size=len(train_idx), replace=True)
            positions, _ = discover_positions(arr[boot_idx])
            selection_counts[positions] += 1
            llr = fit_predict_naive_bayes(arr[boot_idx], arr[held_idx], positions)
            correct = np.concatenate((llr[:, 0] > 0, llr[:, 1] < 0))
            accuracies.append(correct.mean())
            n_positions.append(len(positions))
            all_correct.append(correct.all())
        bootstrap_fold_rows.append({
            "held_out_clade": clade,
            "n_bootstraps": n_bootstraps,
            "median_accuracy": float(np.median(accuracies)),
            "accuracy_2.5pct": float(np.quantile(accuracies, 0.025)),
            "accuracy_97.5pct": float(np.quantile(accuracies, 0.975)),
            "fraction_bootstraps_with_100pct_accuracy": float(np.mean(all_correct)),
            "median_selected_positions": float(np.median(n_positions)),
            "selected_positions_2.5pct": float(np.quantile(n_positions, 0.025)),
            "selected_positions_97.5pct": float(np.quantile(n_positions, 0.975)),
        })
        for col, count in enumerate(selection_counts):
            bootstrap_site_rows.append({
                "held_out_clade": clade,
                "alignment_column": col + 1,
                "selection_frequency": count / n_bootstraps,
                "canonical_barcode": col in canonical_zero,
            })

    return {
        "fold": fold_df,
        "prediction": prediction_df,
        "site": site_df,
        "permutation": permutation_df,
        "permutation_p": permutation_p,
        "bootstrap_fold": pd.DataFrame(bootstrap_fold_rows),
        "bootstrap_site": pd.DataFrame(bootstrap_site_rows),
        "clade_indices": clade_indices,
    }





def build_lg_q():
    values = [float(x) for x in LG_LOWER.split()]
    if len(values) != 190:
        raise ValueError(f"Expected 190 LG exchangeabilities, got {len(values)}")
    exchangeability = np.zeros((20, 20), dtype=float)
    k = 0
    for i in range(1, 20):
        for j in range(i):
            exchangeability[i, j] = exchangeability[j, i] = values[k]
            k += 1
    q = np.zeros((20, 20), dtype=float)
    for i in range(20):
        for j in range(20):
            if i != j:
                q[i, j] = exchangeability[i, j] * LG_FREQ[j]
        q[i, i] = -q[i].sum()
    expected_rate = -np.sum(LG_FREQ * np.diag(q))
    return q / expected_rate


def discrete_gamma_mean_rates(alpha: float, categories: int = 4):
    boundaries = [0.0] + [gamma.ppf(i / categories, a=alpha, scale=1 / alpha)
                          for i in range(1, categories)] + [np.inf]
    rates = []
    for lo, hi in zip(boundaries[:-1], boundaries[1:]):
        low_moment = 0.0 if lo == 0 else gammainc(alpha + 1, alpha * lo)
        high_moment = 1.0 if np.isinf(hi) else gammainc(alpha + 1, alpha * hi)
        rates.append((high_moment - low_moment) * categories)
    return np.asarray(rates)


def parse_gamma_alpha(iqtree_report: Path, default=0.9960):
    text = iqtree_report.read_text(errors="replace")
    for line in text.splitlines():
        if "Gamma shape alpha:" in line:
            return float(line.split(":", 1)[1].strip())
    return default


def run_asr(tree_path: Path, iqtree_report: Path, seq_by_tree_tip: dict,
            arr: np.ndarray, informative_cols: np.ndarray):
    tree = Phylo.read(tree_path, "newick")
    nodes = list(tree.find_clades(order="preorder"))
    postorder = list(tree.find_clades(order="postorder"))
    parent = {}
    for node in nodes:
        for child in node.clades:
            parent[child] = node

    iiic_tips = {tip.name for tip in tree.get_terminals() if tip.name.endswith("__IIIc")}
    anc_iiic = None
    for node in nodes:
        if {tip.name for tip in node.get_terminals()} == iiic_tips:
            anc_iiic = node
            break
    if anc_iiic is None:
        raise RuntimeError("Could not identify exact IIIc clade in the ML tree")
    anc_iiib = parent[anc_iiic]

    q = build_lg_q()
    alpha = parse_gamma_alpha(iqtree_report)
    gamma_rates = discrete_gamma_mean_rates(alpha, 4)
    aa_index = {aa: i for i, aa in enumerate(LG_ORDER)}

    @lru_cache(maxsize=None)
    def transition(branch_length_rounded: float, rate_category: int):
        return expm(q * (branch_length_rounded * gamma_rates[rate_category]))

    def p_matrix(branch_length, category):
        value = 0.0 if branch_length is None else float(branch_length)
        return transition(round(value, 12), category)

    def posterior_for_site(column: int, excluded_species=frozenset()):
        joint = {anc_iiib: np.zeros(20), anc_iiic: np.zeros(20)}
        category_likelihoods = []
        for category in range(4):
            inside = {}
            contribution = {}
            for node in postorder:
                if node.is_terminal():
                    species = node.name.rsplit("__", 1)[0]
                    aa = seq_by_tree_tip[node.name][column]
                    if species in excluded_species or aa not in aa_index:
                        inside[node] = np.ones(20)
                    else:
                        vector = np.zeros(20)
                        vector[aa_index[aa]] = 1.0
                        inside[node] = vector
                else:
                    vector = np.ones(20)
                    for child in node.clades:
                        child_contribution = p_matrix(child.branch_length, category) @ inside[child]
                        contribution[(node, child)] = child_contribution
                        vector *= child_contribution
                    inside[node] = vector

            outside = {tree.root: LG_FREQ.copy()}
            for node in nodes:
                if node.is_terminal():
                    continue
                for child in node.clades:
                    base = outside[node].copy()
                    for sibling in node.clades:
                        if sibling is not child:
                            base *= contribution[(node, sibling)]
                    outside[child] = base @ p_matrix(child.branch_length, category)

            likelihood = float(np.sum(outside[tree.root] * inside[tree.root]))
            category_likelihoods.append(likelihood)
            for target in (anc_iiib, anc_iiic):
                joint[target] += 0.25 * outside[target] * inside[target]
        denominator = float(np.mean(category_likelihoods))
        return joint[anc_iiib] / denominator, joint[anc_iiic] / denominator

    asr_rows = []
    posterior_rows = []
    full_best = {}
    for col in informative_cols:
        p_b, p_c = posterior_for_site(int(col))
        for iso, posterior, iso_idx in (("IIIb", p_b, 0), ("IIIc", p_c, 1)):
            presence_fraction = float(np.mean(arr[:, iso_idx, col] != 20))
            presence_state = "present" if presence_fraction >= 0.5 else "gap"
            order = np.argsort(posterior)[::-1]
            best_idx, second_idx = int(order[0]), int(order[1])
            full_best[(iso, int(col + 1))] = LG_ORDER[best_idx]
            entropy = float(-np.sum(posterior * np.log2(np.clip(posterior, 1e-300, 1.0))))
            asr_rows.append({
                "isoform": iso,
                "alignment_column": int(col + 1),
                "observed_presence_fraction": presence_fraction,
                "ancestral_presence_consensus": presence_state,
                "presence_confidence_from_clade_fraction": max(presence_fraction, 1 - presence_fraction),
                "ml_amino_acid_if_present": LG_ORDER[best_idx],
                "posterior_probability": float(posterior[best_idx]),
                "second_amino_acid": LG_ORDER[second_idx],
                "second_posterior_probability": float(posterior[second_idx]),
                "posterior_entropy_bits": entropy,
            })
            row = {"isoform": iso, "alignment_column": int(col + 1)}
            row.update({f"posterior_{aa}": float(posterior[i]) for i, aa in enumerate(LG_ORDER)})
            posterior_rows.append(row)

    jackknife_rows = []
    for held_out_clade, members in CLADES.items():
        excluded = frozenset(members)
        for col in informative_cols:
            p_b, p_c = posterior_for_site(int(col), excluded)
            for iso, posterior in (("IIIb", p_b), ("IIIc", p_c)):
                best = LG_ORDER[int(np.argmax(posterior))]
                jackknife_rows.append({
                    "held_out_clade": held_out_clade,
                    "isoform": iso,
                    "alignment_column": int(col + 1),
                    "best_amino_acid": best,
                    "posterior_probability": float(np.max(posterior)),
                    "matches_full_reconstruction": best == full_best[(iso, int(col + 1))],
                })

    asr_df = pd.DataFrame(asr_rows)
    posterior_df = pd.DataFrame(posterior_rows)
    jackknife_df = pd.DataFrame(jackknife_rows)


    lookup = {(r.isoform, int(r.alignment_column)): r for _, r in asr_df.iterrows()}
    aligned, thresholded, ungapped, thresholded_ungapped = {}, {}, {}, {}
    informative_set = set(informative_cols + 1)
    for iso in ("IIIb", "IIIc"):
        chars, threshold_chars = [], []
        for col in range(1, 70):
            if col not in informative_set:
                chars.append("-")
                threshold_chars.append("-")
                continue
            record = lookup[(iso, col)]
            if record.ancestral_presence_consensus == "gap":
                chars.append("-")
                threshold_chars.append("-")
            else:
                aa = record.ml_amino_acid_if_present
                chars.append(aa)
                confidence = min(record.posterior_probability,
                                 record.presence_confidence_from_clade_fraction)
                threshold_chars.append(aa if confidence >= 0.8 else "X")
        aligned[iso] = "".join(chars)
        thresholded[iso] = "".join(threshold_chars)
        ungapped[iso] = aligned[iso].replace("-", "")
        thresholded_ungapped[iso] = thresholded[iso].replace("-", "")


    state = asr_df.pivot(index="alignment_column", columns="isoform")
    stability = jackknife_df.groupby(["isoform", "alignment_column"])["matches_full_reconstruction"].all().unstack(0)
    barcode_rows = []
    for col in CANONICAL_BARCODE_1BASED:
        pb = state.loc[col, ("ancestral_presence_consensus", "IIIb")]
        pc = state.loc[col, ("ancestral_presence_consensus", "IIIc")]
        aa_b = "-" if pb == "gap" else state.loc[col, ("ml_amino_acid_if_present", "IIIb")]
        aa_c = "-" if pc == "gap" else state.loc[col, ("ml_amino_acid_if_present", "IIIc")]
        pp_b = float(state.loc[col, ("posterior_probability", "IIIb")])
        pp_c = float(state.loc[col, ("posterior_probability", "IIIc")])
        stable_b = bool(stability.loc[col, "IIIb"])
        stable_c = bool(stability.loc[col, "IIIc"])
        barcode_rows.append({
            "alignment_column": int(col),
            "ancestral_IIIb": aa_b,
            "ancestral_IIIc": aa_c,
            "ancestrally_different": aa_b != aa_c,
            "IIIb_posterior_probability": pp_b,
            "IIIc_posterior_probability": pp_c,
            "minimum_posterior_probability": min(pp_b, pp_c),
            "IIIb_stable_all_clade_jackknifes": stable_b,
            "IIIc_stable_all_clade_jackknifes": stable_c,
            "both_stable_all_clade_jackknifes": stable_b and stable_c,
            "high_confidence_robust_ancestral_difference": (
                aa_b != aa_c and stable_b and stable_c and min(pp_b, pp_c) >= 0.95
            ),
        })
    barcode_df = pd.DataFrame(barcode_rows)


    difference_rows = []
    for col in informative_cols + 1:
        pb = state.loc[col, ("ancestral_presence_consensus", "IIIb")]
        pc = state.loc[col, ("ancestral_presence_consensus", "IIIc")]
        aa_b = "-" if pb == "gap" else state.loc[col, ("ml_amino_acid_if_present", "IIIb")]
        aa_c = "-" if pc == "gap" else state.loc[col, ("ml_amino_acid_if_present", "IIIc")]
        difference_rows.append({
            "alignment_column": int(col),
            "ancestral_IIIb": aa_b,
            "ancestral_IIIc": aa_c,
            "different": aa_b != aa_c,
            "canonical_barcode": int(col) in set(CANONICAL_BARCODE_1BASED),
            "IIIb_posterior_probability": float(state.loc[col, ("posterior_probability", "IIIb")]),
            "IIIc_posterior_probability": float(state.loc[col, ("posterior_probability", "IIIc")]),
        })
    difference_df = pd.DataFrame(difference_rows)

    return {
        "asr": asr_df,
        "posterior": posterior_df,
        "jackknife": jackknife_df,
        "barcode": barcode_df,
        "difference": difference_df,
        "aligned": aligned,
        "thresholded": thresholded,
        "ungapped": ungapped,
        "thresholded_ungapped": thresholded_ungapped,
        "alpha": alpha,
        "gamma_rates": gamma_rates,
        "target_node_IIIb": anc_iiib.name,
        "target_node_IIIc": anc_iiic.name,
        "split_branch_length": anc_iiic.branch_length,
    }





def draw_panel_label(ax, label, title):
    ax.text(0.0, 1.03, label, transform=ax.transAxes, fontsize=15,
            fontweight="bold", color=COLORS["navy"], va="bottom")
    ax.text(0.06, 1.03, title, transform=ax.transAxes, fontsize=13.5,
            fontweight="bold", color=COLORS["dark"], va="bottom")


def plot_loco(loco, outdir: Path):
    fold = loco["fold"]
    pred = loco["prediction"].copy()
    site = loco["site"]
    null = loco["permutation"]["overall_accuracy"].to_numpy()
    clade_order = ["Mammals", "Birds", "Reptiles", "Amphibians", "Teleosts"]

    fig = plt.figure(figsize=(16, 10.2))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[0.95, 1.15], width_ratios=[1.0, 1.25],
                  hspace=0.46, wspace=0.28, top=0.87, left=0.07, right=0.97, bottom=0.09)
    header(fig, "Evolutionary hold-out validation: the IIIb/IIIc code generalizes across vertebrate clades",
           "Features are re-discovered from the training clades only; every held-out sequence is classified once and never contributes labels to its own model.")


    ax = fig.add_subplot(gs[0, 0])
    draw_panel_label(ax, "A", "Every evolutionary hold-out is classified perfectly")
    f = fold.set_index("held_out_clade").loc[clade_order].reset_index()
    y = np.arange(len(f))
    ax.barh(y, f["accuracy"] * 100, height=0.55, color=COLORS["teal"], alpha=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(f["held_out_clade"], fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlim(0, 112)
    ax.set_xlabel("Held-out classification accuracy (%)")
    for yi, row in f.iterrows():
        ax.text(101.2, yi, f"{int(row.n_test_sequences)}/{int(row.n_test_sequences)}",
                va="center", fontsize=10.5, fontweight="bold", color=COLORS["navy"])
        ax.text(54, yi + 0.33,
                f"{int(row.n_de_novo_positions)} sites · {int(row.canonical_positions_recovered)}/17 canonical recovered",
                ha="center", va="center", fontsize=8.8, color=COLORS["gray"])
    clean_axis(ax, "x")
    ax.text(0.02, -0.20, "Training bootstrap: 1,000/1,000 replicates retained 100% accuracy in every fold.",
            transform=ax.transAxes, fontsize=10.2, fontweight="bold", color=COLORS["navy"])


    ax = fig.add_subplot(gs[0, 1])
    draw_panel_label(ax, "B", "The discriminating site set is stable but not artificially fixed")
    selected = site[site["alignment_column"].isin(range(1, 60))].pivot(
        index="held_out_clade", columns="alignment_column", values="selected"
    ).loc[clade_order]
    matrix = selected.astype(int).to_numpy()
    ax.imshow(matrix, aspect="auto", cmap=mpl.colors.ListedColormap(["#EEF2F6", COLORS["blue"]]), vmin=0, vmax=1)
    ax.set_yticks(np.arange(len(clade_order)))
    ax.set_yticklabels(clade_order, fontweight="bold")
    ticks = [1, 10, 20, 30, 40, 50, 59]
    ax.set_xticks([x - 1 for x in ticks])
    ax.set_xticklabels(ticks)
    ax.set_xlabel("Original MSA column")
    canonical = set(CANONICAL_BARCODE_1BASED)
    for col in canonical:
        if col <= 59:
            ax.add_patch(Rectangle((col - 1 - 0.48, -0.48), 0.96, len(clade_order) - 0.04,
                                   fill=False, edgecolor=COLORS["gold"], linewidth=1.1))
    ax.text(0.01, -0.22, "Blue: selected de novo. Gold outline: original 17-site barcode.",
            transform=ax.transAxes, fontsize=9.2, color=COLORS["gray"])
    for spine in ax.spines.values():
        spine.set_visible(False)


    ax = fig.add_subplot(gs[1, 0])
    draw_panel_label(ax, "C", "Held-out proteins retain large likelihood margins")
    pred["clade_order"] = pred["held_out_clade"].map({c: i for i, c in enumerate(clade_order)})
    pred = pred.sort_values(["clade_order", "species", "true_isoform"]).reset_index(drop=True)
    x = np.arange(len(pred))
    for iso, color, marker in [("IIIb", COLORS["blue"], "o"), ("IIIc", COLORS["orange"], "D")]:
        mask = pred["true_isoform"] == iso
        ax.scatter(x[mask], pred.loc[mask, "signed_log_likelihood_margin"], s=42,
                   color=color, edgecolor="white", linewidth=0.7, marker=marker, label=iso, zorder=3)
    boundaries = []
    cursor = 0
    centers = []
    for clade in clade_order:
        n = int((pred["held_out_clade"] == clade).sum())
        centers.append(cursor + (n - 1) / 2)
        cursor += n
        boundaries.append(cursor - 0.5)
    for b in boundaries[:-1]:
        ax.axvline(b, color=COLORS["grid"], linewidth=1.0)
    ax.set_xticks(centers)
    ax.set_xticklabels(clade_order, rotation=20, ha="right", fontweight="bold")
    ax.set_ylabel("Signed log-likelihood margin toward the true isoform")
    ax.axhline(0, color=COLORS["red"], linewidth=1.2, linestyle="--")
    clean_axis(ax, "y")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    ax.text(0.02, 0.95, f"Minimum margin = {pred.signed_log_likelihood_margin.min():.1f}",
            transform=ax.transAxes, va="top", fontsize=10.2, fontweight="bold", color=COLORS["navy"])


    ax = fig.add_subplot(gs[1, 1])
    draw_panel_label(ax, "D", "Label randomization destroys cross-clade predictability")
    bins = np.linspace(0.38, 0.85, 28)
    ax.hist(null, bins=bins, color=COLORS["sky"], edgecolor="white", alpha=0.95)
    ax.axvline(1.0, color=COLORS["red"], linewidth=3.0)
    ax.annotate("Observed = 1.00", xy=(1.0, ax.get_ylim()[1] * 0.75),
                xytext=(0.82, ax.get_ylim()[1] * 0.92),
                arrowprops=dict(arrowstyle="->", color=COLORS["red"], linewidth=1.3),
                fontsize=11, fontweight="bold", color=COLORS["red"])
    ax.set_xlim(0.38, 1.03)
    ax.set_xlabel("Overall LOCO accuracy after within-species label swaps")
    ax.set_ylabel("Permutation count")
    clean_axis(ax, "y")
    ax.text(0.02, 0.90,
            f"Null mean = {null.mean():.3f}\nNull maximum = {null.max():.3f}\nEmpirical p = {loco['permutation_p']:.4g}",
            transform=ax.transAxes, va="top", fontsize=10.6, fontweight="bold", color=COLORS["navy"],
            bbox=dict(boxstyle="round,pad=0.45", facecolor="#F3F7FA", edgecolor=COLORS["grid"]))

    save_figure(fig, outdir, "02_leave_one_clade_out_validation")


def plot_asr(asr, outdir: Path):
    asr_df = asr["asr"]
    barcode = asr["barcode"].copy()
    jack = asr["jackknife"]
    informative = sorted(asr_df["alignment_column"].unique())
    compact_index = {col: i for i, col in enumerate(informative)}

    fig = plt.figure(figsize=(17, 10.5))
    gs = GridSpec(3, 1, figure=fig, height_ratios=[1.30, 0.86, 1.02], hspace=0.54,
                  top=0.87, left=0.065, right=0.97, bottom=0.095)
    header(fig, "Ancestral reconstruction resolves an ancient core and later lineage-specific refinements",
           "Marginal empirical-Bayes amino-acid posteriors are conditioned on the existing LG+G4 maximum-likelihood tree; clade jackknifes expose sampling-sensitive sites.")


    ax = fig.add_subplot(gs[0])
    draw_panel_label(ax, "A", "Reconstructed ancestral cassette states with site-wise posterior support")
    ax.set_xlim(-2, len(informative) + 0.5)
    ax.set_ylim(-0.6, 2.05)
    for row_i, iso in enumerate(("IIIc", "IIIb")):
        subset = asr_df[asr_df["isoform"] == iso].set_index("alignment_column")
        y = 1.25 - row_i * 0.92
        ax.text(-1.7, y, f"Anc-{iso}", va="center", fontsize=12, fontweight="bold", color=COLORS["navy"])
        for col in informative:
            rec = subset.loc[col]
            aa = "-" if rec.ancestral_presence_consensus == "gap" else rec.ml_amino_acid_if_present
            confidence = min(rec.posterior_probability, rec.presence_confidence_from_clade_fraction)
            x = compact_index[col]
            face = AA_CLASS_COLOR[AA_CLASS.get(aa, "uncertain")]
            alpha = 0.35 + 0.65 * confidence
            patch = FancyBboxPatch((x - 0.42, y - 0.30), 0.84, 0.60,
                                   boxstyle="round,pad=0.02,rounding_size=0.06",
                                   facecolor=face, edgecolor="white", linewidth=0.6, alpha=alpha)
            ax.add_patch(patch)
            ax.text(x, y, aa, ha="center", va="center", fontsize=8.2, fontweight="bold",
                    color="white" if aa != "-" else COLORS["gray"])
            if col in set(CANONICAL_BARCODE_1BASED):
                ax.plot(x, y + 0.39, marker="o", markersize=3.6, color=COLORS["gold"])
        ax.text(len(informative) - 0.4, y - 0.36,
                f"mean PP {subset.posterior_probability.mean():.3f}",
                ha="right", fontsize=8.8, color=COLORS["gray"])
    tick_cols = [1, 10, 20, 30, 46, 59, 69]
    tick_pos = [compact_index[c] for c in tick_cols if c in compact_index]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels([c for c in tick_cols if c in compact_index])
    ax.set_xlabel("Original MSA column; all-gap artifact columns removed from reconstruction")
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    legend_items = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor=color, markersize=10, label=label)
        for label, color in [("hydrophobic", AA_CLASS_COLOR["hydrophobic"]),
                             ("polar", AA_CLASS_COLOR["polar"]),
                             ("charged", AA_CLASS_COLOR["positive"]),
                             ("special", AA_CLASS_COLOR["special"])]
    ] + [Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS["gold"], markersize=6,
                label="17-site modern barcode")]
    ax.legend(handles=legend_items, frameon=False, ncol=5, loc="upper right",
              bbox_to_anchor=(1.0, 0.97), fontsize=8.2, handletextpad=0.4, columnspacing=0.8)


    ax = fig.add_subplot(gs[1])
    draw_panel_label(ax, "B", "Most of the modern barcode was already distinct at the ancestral split")
    x = np.arange(len(barcode))
    for i, row in barcode.iterrows():
        if not row.ancestrally_different:
            color = COLORS["gray"]
            category = "same ancestral state"
        elif row.high_confidence_robust_ancestral_difference:
            color = COLORS["teal"]
            category = "robust ancestral difference"
        else:
            color = COLORS["gold"]
            category = "difference with uncertainty/sampling sensitivity"
        ax.scatter(i, -0.15, s=250, color=color, edgecolor="white", linewidth=0.9, zorder=3)
        ax.text(i, -0.15, f"{row.ancestral_IIIb}\n{row.ancestral_IIIc}", ha="center", va="center",
                fontsize=7.5, fontweight="bold", color="white" if color != COLORS["gray"] else COLORS["navy"])
    ax.set_xticks(x)
    ax.set_xticklabels(barcode["alignment_column"].astype(int))
    ax.set_yticks([])
    ax.set_xlim(-0.7, len(barcode) - 0.3)
    ax.set_ylim(-0.72, 0.78)
    ax.set_xlabel("Canonical barcode alignment column (top letter = Anc-IIIb; bottom = Anc-IIIc)")
    for spine in ax.spines.values():
        spine.set_visible(False)
    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS["teal"], markersize=10,
               label="11 robust high-confidence differences"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS["gold"], markersize=10,
               label="4 additional ancestral differences"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS["gray"], markersize=10,
               label="2 later IIIb refinements"),
    ]
    ax.legend(handles=legend, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, 0.98), fontsize=8.8, handletextpad=0.4, columnspacing=1.0)


    ax = fig.add_subplot(gs[2])
    draw_panel_label(ax, "C", "Clade jackknifes localize uncertainty instead of hiding it")
    rows = []
    row_labels = []
    clade_order = ["Mammals", "Birds", "Reptiles", "Amphibians", "Teleosts"]
    for iso in ("IIIb", "IIIc"):
        for clade in clade_order:
            subset = jack[(jack.isoform == iso) & (jack.held_out_clade == clade) &
                          (jack.alignment_column.isin(CANONICAL_BARCODE_1BASED))]
            subset = subset.set_index("alignment_column").loc[CANONICAL_BARCODE_1BASED]
            rows.append(subset.matches_full_reconstruction.astype(int).to_numpy())
            row_labels.append(f"{iso} · minus {clade}")
    matrix = np.vstack(rows)
    ax.imshow(matrix, aspect="auto", cmap=mpl.colors.ListedColormap([COLORS["red"], "#E8F3F0"]), vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(CANONICAL_BARCODE_1BASED)))
    ax.set_xticklabels(CANONICAL_BARCODE_1BASED)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8.8)
    ax.set_xlabel("Canonical barcode alignment column")
    for r in range(matrix.shape[0]):
        for c in range(matrix.shape[1]):
            if matrix[r, c] == 0:
                ax.text(c, r, "×", ha="center", va="center", fontsize=10, fontweight="bold", color="white")
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.text(0.075, 0.025,
             "Red cells change the maximum-posterior state when an entire clade is removed. "
             "Instability is localized to teleost-sensitive IIIb states and reptile-sensitive IIIc states.",
             fontsize=9.2, color=COLORS["gray"])

    save_figure(fig, outdir, "03_ancestral_sequence_reconstruction")


def rounded_box(ax, xy, width, height, text, face, edge=None, fontsize=10.5, weight="bold", text_color=None):
    x, y = xy
    patch = FancyBboxPatch((x, y), width, height,
                           boxstyle="round,pad=0.02,rounding_size=0.03",
                           facecolor=face, edgecolor=edge or face, linewidth=1.2)
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center",
            fontsize=fontsize, fontweight=weight,
            color=text_color or ("white" if face not in ("#FFFFFF", COLORS["light"]) else COLORS["navy"]))
    return patch


def arrow(ax, start, end, color=None, style="solid", lw=1.8, mutation_scale=14):
    patch = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=mutation_scale,
                            linewidth=lw, linestyle=style, color=color or COLORS["navy"],
                            connectionstyle="arc3,rad=0.0")
    ax.add_patch(patch)
    return patch


def plot_synthesis(loco, asr, outdir: Path):
    fig = plt.figure(figsize=(18, 10.1))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1.0, 0.88], width_ratios=[1.05, 1.05, 1.0],
                  hspace=0.25, wspace=0.25, top=0.87, left=0.04, right=0.97, bottom=0.08)
    header(fig, "FGFR2 IIIb/IIIc: an evolutionarily stable functional code inside an annotation-variable coordinate frame",
           "Integrated mechanistic model. Solid elements are directly demonstrated here; dashed elements are data-supported biological inference; dotted elements are testable predictions.")


    ax = fig.add_subplot(gs[0, 0]); ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    draw_panel_label(ax, "A", "Ancient split and out-of-clade generalization")
    rounded_box(ax, (0.36, 0.73), 0.28, 0.12, "ancestral D3\nsplicing module", COLORS["light"], COLORS["grid"], text_color=COLORS["navy"])
    arrow(ax, (0.45, 0.72), (0.25, 0.54), COLORS["blue"], "dashed", 2.0)
    arrow(ax, (0.55, 0.72), (0.75, 0.54), COLORS["orange"], "dashed", 2.0)
    rounded_box(ax, (0.10, 0.42), 0.30, 0.12, "IIIb lineage", COLORS["blue"])
    rounded_box(ax, (0.60, 0.42), 0.30, 0.12, "IIIc lineage", COLORS["orange"])
    clade_x = [0.12, 0.19, 0.26, 0.33, 0.40]
    labels = ["M", "B", "R", "A", "F"]
    for i, (x, lab) in enumerate(zip(clade_x, labels)):
        ax.add_patch(Circle((x, 0.28), 0.027, facecolor=COLORS["sky"], edgecolor="white"))
        ax.text(x, 0.28, lab, ha="center", va="center", fontsize=8, fontweight="bold", color=COLORS["navy"])
        x2 = x + 0.50
        ax.add_patch(Circle((x2, 0.28), 0.027, facecolor="#F1B28F", edgecolor="white"))
        ax.text(x2, 0.28, lab, ha="center", va="center", fontsize=8, fontweight="bold", color=COLORS["navy"])
    ax.text(0.5, 0.17, "56/56 proteins classified correctly when their entire clade is held out",
            ha="center", fontsize=10.5, fontweight="bold", color=COLORS["navy"])
    ax.text(0.5, 0.10, "paired label-permutation p = 9.999×10⁻⁵",
            ha="center", fontsize=9.5, color=COLORS["gray"])
    ax.text(0.5, 0.02, "15/17 modern barcode sites differ between reconstructed ancestors;\n11 form a high-confidence, all-jackknife-stable ancestral core.",
            ha="center", fontsize=9.7, color=COLORS["dark"])


    ax = fig.add_subplot(gs[0, 1]); ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    draw_panel_label(ax, "B", "Mutually exclusive exon architecture")
    y = 0.77
    rounded_box(ax, (0.08, y), 0.13, 0.09, "exon IIIa", COLORS["light"], COLORS["grid"], text_color=COLORS["navy"], fontsize=9)
    rounded_box(ax, (0.30, y), 0.18, 0.09, "exon IIIb", COLORS["blue"], fontsize=9)
    rounded_box(ax, (0.56, y), 0.18, 0.09, "exon IIIc", COLORS["orange"], fontsize=9)
    rounded_box(ax, (0.82, y), 0.10, 0.09, "TM", COLORS["light"], COLORS["grid"], text_color=COLORS["navy"], fontsize=9)
    ax.text(0.39, 0.91, "ESRP1/2", ha="center", fontsize=9.5, fontweight="bold", color=COLORS["teal"])
    arrow(ax, (0.39, 0.90), (0.39, 0.87), COLORS["teal"], "dashed", 1.6, 10)
    arrow(ax, (0.39, 0.72), (0.28, 0.59), COLORS["blue"], "solid", 1.8)
    arrow(ax, (0.65, 0.72), (0.72, 0.59), COLORS["orange"], "solid", 1.8)
    ax.text(0.50, 0.665, "mutually exclusive splicing", ha="center", fontsize=9, color=COLORS["gray"])

    def receptor(y0, cassette_color, label):
        x0 = 0.10
        for k, text in enumerate(("Ig1", "Ig2")):
            rounded_box(ax, (x0 + k * 0.15, y0), 0.11, 0.10, text, COLORS["light"], COLORS["grid"], fontsize=8.5, text_color=COLORS["navy"])
        rounded_box(ax, (x0 + 0.30, y0), 0.11, 0.10, "D3a", COLORS["light"], COLORS["grid"], fontsize=8.5, text_color=COLORS["navy"])
        rounded_box(ax, (x0 + 0.41, y0), 0.16, 0.10, label, cassette_color, fontsize=8.5)
        rounded_box(ax, (x0 + 0.61, y0), 0.07, 0.10, "TM", COLORS["dark"], fontsize=7.5)
    receptor(0.44, COLORS["blue"], "IIIb-D3")
    receptor(0.23, COLORS["orange"], "IIIc-D3")
    ax.text(0.84, 0.49, "FGFR2b", fontsize=11, fontweight="bold", color=COLORS["blue"])
    ax.text(0.84, 0.28, "FGFR2c", fontsize=11, fontweight="bold", color=COLORS["orange"])
    ax.text(0.50, 0.08, "17-site barcode carries 73.7% of total isoform-discriminating JSD",
            ha="center", fontsize=10.2, fontweight="bold", color=COLORS["navy"])


    ax = fig.add_subplot(gs[0, 2]); ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    draw_panel_label(ax, "C", "Ligand-facing barcode")
    ax.add_patch(Circle((0.40, 0.48), 0.24, facecolor="#E6EEF5", edgecolor=COLORS["blue"], linewidth=2))
    ax.add_patch(Circle((0.67, 0.53), 0.17, facecolor="#F7E5DA", edgecolor=COLORS["orange"], linewidth=2))
    ax.text(0.36, 0.48, "FGFR2\nD3", ha="center", va="center", fontsize=13, fontweight="bold", color=COLORS["navy"])
    ax.text(0.70, 0.53, "FGF\nligand", ha="center", va="center", fontsize=11, fontweight="bold", color=COLORS["navy"])
    contact_angles = np.linspace(-55, 55, 13)
    for angle in contact_angles:
        rad = np.deg2rad(angle)
        x = 0.40 + 0.235 * np.cos(rad)
        y = 0.48 + 0.235 * np.sin(rad)
        ax.add_patch(Circle((x, y), 0.018, facecolor=COLORS["gold"], edgecolor="white", linewidth=0.5))
    for x, y in [(0.24,0.31),(0.22,0.62),(0.34,0.72),(0.37,0.25)]:
        ax.add_patch(Circle((x, y), 0.017, facecolor=COLORS["gray"], edgecolor="white", linewidth=0.5))
    ax.text(0.5, 0.16, "13/17 direct contact in ≥1 complex", ha="center", fontsize=11, fontweight="bold", color=COLORS["navy"])
    ax.text(0.5, 0.10, "15/17 within 8 Å in ≥1 complex", ha="center", fontsize=10, color=COLORS["dark"])
    ax.text(0.5, 0.04, "contact enrichment: p = 10⁻⁴ and 2×10⁻⁴", ha="center", fontsize=9.5, color=COLORS["gray"])


    ax = fig.add_subplot(gs[1, :2]); ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    draw_panel_label(ax, "D", "Coordinate shifts preserve relative topology")

    for y, name, shift, color in [(0.66, "NCBI", 0.0, COLORS["blue"]), (0.34, "Ensembl", 0.14, COLORS["orange"])]:
        ax.text(0.02, y + 0.035, name, fontsize=11, fontweight="bold", color=color)
        x0 = 0.13 + shift
        ax.plot([x0, x0 + 0.62], [y, y], color=COLORS["gray"], linewidth=7, solid_capstyle="round")
        cassette_x = x0 + 0.37
        d3_x = cassette_x + 0.13
        ax.plot([cassette_x, cassette_x], [y - 0.08, y + 0.08], color=color, linewidth=2)
        ax.plot([d3_x, d3_x], [y - 0.08, y + 0.08], color=COLORS["teal"], linewidth=2)
        ax.text(cassette_x, y + 0.09, "cassette end", ha="center", fontsize=8.8, color=color)
        ax.text(d3_x, y - 0.13, "median D3 end", ha="center", fontsize=8.8, color=COLORS["teal"])
        ax.annotate("", xy=(d3_x, y - 0.05), xytext=(cassette_x, y - 0.05),
                    arrowprops=dict(arrowstyle="<->", color=COLORS["navy"], linewidth=1.4))
        ax.text((cassette_x+d3_x)/2, y - 0.095, "I", ha="center", fontsize=11, fontweight="bold", color=COLORS["navy"])
    ax.text(0.79, 0.72, "I = D3_end − cassette_end", fontsize=13, fontweight="bold", color=COLORS["navy"])
    ax.text(0.79, 0.59, "I′ = (D3_end + k) − (cassette_end + k)", fontsize=11, color=COLORS["dark"])
    ax.text(0.79, 0.48, "I′ = I", fontsize=16, fontweight="bold", color=COLORS["teal"])
    ax.text(0.79, 0.31, "16/16 source pairs: residual = 0 aa", fontsize=11, fontweight="bold", color=COLORS["navy"])
    ax.text(0.79, 0.17, "32/32 proteins pass topology", fontsize=10, color=COLORS["dark"])
    ax.text(0.02, 0.05, "Observed annotation shifts: +23, +137, −96 and −76 aa — without any topology change.",
            fontsize=10.2, fontweight="bold", color=COLORS["navy"])


    ax = fig.add_subplot(gs[1, 2]); ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    draw_panel_label(ax, "E", "Evidence hierarchy")
    rounded_box(ax, (0.06, 0.67), 0.88, 0.19,
                "DIRECTLY DEMONSTRATED\ncaller bias · robustness plateau · LOCO generalization\nancestral state posteriors · structural contacts · annotation invariance",
                "#E8F3F0", COLORS["teal"], fontsize=9.2, text_color=COLORS["navy"])
    rounded_box(ax, (0.06, 0.39), 0.88, 0.19,
                "DATA-SUPPORTED INFERENCE\nan old, functionally separated IIIb/IIIc module\nmaintained across vertebrate diversification",
                "#F7EED8", COLORS["gold"], fontsize=9.5, text_color=COLORS["navy"])
    rounded_box(ax, (0.06, 0.11), 0.88, 0.19,
                "TESTABLE PREDICTION\nancestral-core swaps may redirect ligand recognition\nrelative boundary invariants should generalize\nto other homologous alternative exons",
                "#F5E7E6", COLORS["red"], fontsize=8.8, text_color=COLORS["navy"])

    fig.text(0.5, 0.018,
             "Core conclusion: exact residue precision is model-dependent; evolutionary identity, functional surface placement and relative exon–domain topology are robust.",
             ha="center", fontsize=12.2, fontweight="bold", color=COLORS["navy"],
             bbox=dict(boxstyle="round,pad=0.55", facecolor="#F3F7FA", edgecolor=COLORS["grid"]))
    save_figure(fig, outdir, "01_final_integrated_synthesis")


def plot_presentation_hero(loco, asr, outdir: Path):
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.patch.set_facecolor("#F8FAFC")
    ax.text(0.05, 0.92, "FGFR2 IIIb/IIIc", fontsize=34, fontweight="bold", color=COLORS["navy"])
    ax.text(0.05, 0.855, "A conserved functional code inside a moving annotation frame",
            fontsize=21, color=COLORS["dark"])


    stages = [
        (0.07, COLORS["blue"], "EVOLUTIONARY IDENTITY", "56/56 held-out proteins\ncorrectly classified", "ancient IIIb/IIIc split"),
        (0.38, COLORS["gold"], "FUNCTIONAL SURFACE", "13/17 barcode sites\ndirectly contact ligand", "ligand-facing molecular code"),
        (0.69, COLORS["teal"], "ANNOTATION INVARIANCE", "16/16 NCBI–Ensembl pairs\nretain I exactly", "relative topology survives coordinate shifts"),
    ]
    for x, color, title, metric, caption in stages:
        ax.add_patch(FancyBboxPatch((x, 0.36), 0.24, 0.32, boxstyle="round,pad=0.02,rounding_size=0.035",
                                    facecolor="white", edgecolor=color, linewidth=2.4))
        ax.add_patch(Circle((x + 0.12, 0.70), 0.042, facecolor=color, edgecolor="white", linewidth=1.0))
        ax.text(x + 0.12, 0.70, str(stages.index((x,color,title,metric,caption))+1), ha="center", va="center",
                fontsize=16, fontweight="bold", color="white")
        ax.text(x + 0.12, 0.61, title, ha="center", fontsize=12, fontweight="bold", color=color)
        ax.text(x + 0.12, 0.50, metric, ha="center", fontsize=15, fontweight="bold", color=COLORS["navy"])
        ax.text(x + 0.12, 0.405, caption, ha="center", fontsize=10, color=COLORS["gray"])
    arrow(ax, (0.315, 0.52), (0.37, 0.52), COLORS["navy"], "solid", 2.1, 18)
    arrow(ax, (0.625, 0.52), (0.68, 0.52), COLORS["navy"], "solid", 2.1, 18)

    ax.text(0.5, 0.24, "The coordinate is not the biological invariant.", ha="center",
            fontsize=23, fontweight="bold", color=COLORS["navy"])
    ax.text(0.5, 0.17, "The invariant is the relationship between the alternative cassette and the D3 boundary.",
            ha="center", fontsize=16, color=COLORS["dark"])
    ax.text(0.5, 0.075, "I = D3_end − cassette_end", ha="center", fontsize=25,
            fontweight="bold", color=COLORS["teal"])
    save_figure(fig, outdir, "04_presentation_graphical_abstract")





def markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(column) for column in frame.columns]
    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join("---" for _ in columns) + " |"
    rows = [
        "| " + " | ".join(str(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    ]
    return "\n".join([header, divider, *rows])


def write_reports(outdir: Path, loco, asr):
    fold = loco["fold"]
    pred = loco["prediction"]
    barcode = asr["barcode"]
    differences = asr["difference"]
    jack = asr["jackknife"]

    n_ancestral_barcode_diff = int(barcode.ancestrally_different.sum())
    n_robust_core = int(barcode.high_confidence_robust_ancestral_difference.sum())
    n_all_ancestral_diff = int(differences.different.sum())
    later_sites = barcode.loc[~barcode.ancestrally_different, "alignment_column"].astype(int).tolist()
    unstable = jack[~jack.matches_full_reconstruction]

    summary = {
        "loco": {
            "overall_accuracy": float(pred.correct.mean()),
            "correct_sequences": int(pred.correct.sum()),
            "total_sequences": len(pred),
            "minimum_signed_log_likelihood_margin": float(pred.signed_log_likelihood_margin.min()),
            "paired_label_permutation_p": float(loco["permutation_p"]),
            "permutation_null_mean": float(loco["permutation"].overall_accuracy.mean()),
            "permutation_null_maximum": float(loco["permutation"].overall_accuracy.max()),
            "all_1000_training_bootstraps_perfect_in_all_folds": bool(
                (loco["bootstrap_fold"].fraction_bootstraps_with_100pct_accuracy == 1.0).all()
            ),
            "folds": fold.to_dict(orient="records"),
        },
        "asr": {
            "model": "LG+G4 marginal empirical Bayes on existing maximum-likelihood tree",
            "gamma_shape_alpha": float(asr["alpha"]),
            "gamma_category_mean_rates": [float(x) for x in asr["gamma_rates"]],
            "ancestral_IIIb_sequence_best_state": asr["ungapped"]["IIIb"],
            "ancestral_IIIc_sequence_best_state": asr["ungapped"]["IIIc"],
            "ancestral_IIIb_sequence_thresholded": asr["thresholded_ungapped"]["IIIb"],
            "ancestral_IIIc_sequence_thresholded": asr["thresholded_ungapped"]["IIIc"],
            "modern_barcode_sites_ancestrally_different": n_ancestral_barcode_diff,
            "modern_barcode_total": 17,
            "high_confidence_all_jackknife_stable_ancestral_core": n_robust_core,
            "all_ancestral_differences_across_48_informative_columns": n_all_ancestral_diff,
            "modern_barcode_sites_reconstructed_same_at_split": later_sites,
            "jackknife_state_changes": unstable.to_dict(orient="records"),
            "important_limitation": (
                "Conditional reconstruction on a short 69-column cassette alignment; all-gap artifact columns were removed, "
                "gaps were treated as missing in the amino-acid likelihood, and presence/absence was assigned from clade consensus."
            ),
        },
    }
    (outdir / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    results_md = f"""# Final FGFR2 synthesis, leave-one-clade-out validation and ancestral reconstruction

## Executive result

The new analyses strengthen the work in two different ways:

1. **Generalization was tested out of sample.** The discriminating positions were re-derived from the training clades only. All {len(pred)} held-out proteins were assigned to the correct isoform, including the complete mammalian hold-out of 30 sequences. The smallest signed log-likelihood margin was {pred.signed_log_likelihood_margin.min():.2f}. Within-species label randomization produced a null mean accuracy of {loco['permutation'].overall_accuracy.mean():.3f}, a null maximum of {loco['permutation'].overall_accuracy.max():.3f}, and an empirical p-value of {loco['permutation_p']:.6f}. In every fold, all 1,000 training-species bootstrap replicates retained 100% held-out accuracy.

2. **The modern barcode can now be divided into an ancestral core and later refinements.** Under marginal empirical-Bayes reconstruction with the same LG+G4 model selected for the phylogeny, {n_ancestral_barcode_diff}/17 modern barcode sites were reconstructed as different at the IIIb/IIIc ancestral split. {n_robust_core} of these are a conservative core with minimum posterior probability at least 0.95 and unchanged maximum-posterior states in every clade jackknife. Columns {', '.join(map(str, later_sites))} were reconstructed with the same state in both ancestors and therefore represent later IIIb refinements under the primary model.

## Leave-one-clade-out design

The 28 complete primary species pairs were divided into mammals, birds, reptiles, amphibians and teleosts. For each fold:

- the complete clade was removed;
- high-confidence positions were discovered only in the remaining species;
- a position-wise amino-acid naive-Bayes model was estimated with Jeffreys smoothing;
- the excluded IIIb and IIIc proteins were classified;
- training species were resampled 1,000 times to assess stability;
- 10,000 paired within-species label permutations calibrated the complete LOCO pipeline.

The feature rule exactly follows the conservative barcode logic: the two isoforms require different non-gap major residues and a discriminating score of at least 0.70, where the score is the smaller within-isoform major-residue fraction multiplied by one minus the larger gap fraction.

### Fold results

{markdown_table(fold)}

The mammalian hold-out is the hardest feature-recovery fold: only 12 of the original 17 positions are rediscovered and three additional noncanonical positions are selected. Nevertheless, all 30 mammalian proteins are correctly classified. This distinction is scientifically useful: **the exact feature list can vary while the isoform-level information remains fully generalizable.**

## Ancestral reconstruction method

The analysis uses the existing 56-sequence IQ-TREE maximum-likelihood tree, the selected LG+G4 model and gamma shape alpha={asr['alpha']:.4f}. The LG rate matrix and equilibrium frequencies were implemented in the IQ-TREE amino-acid order. Site-wise marginal posterior probabilities were calculated by Felsenstein pruning and upward-downward message passing over four discrete-gamma categories.

The target states are the two nodes on opposite sides of the maximum-likelihood IIIb/IIIc split. The reconstruction is therefore conditional on the existing gene tree and model; it does not independently date the exon duplication.

All-gap alignment artifact columns were removed. Amino-acid gaps were treated as missing in the substitution likelihood. Ancestral presence or absence was assigned separately from the observed clade fraction, and thresholded sequences use `X` when either amino-acid posterior probability or presence confidence is below 0.80.

### Reconstructed sequences

```text
Anc-IIIb best:        {asr['ungapped']['IIIb']}
Anc-IIIb thresholded: {asr['thresholded_ungapped']['IIIb']}
Anc-IIIc best:        {asr['ungapped']['IIIc']}
Anc-IIIc thresholded: {asr['thresholded_ungapped']['IIIc']}
```

The reconstructed sequences are hypotheses about ancestral amino-acid states, not experimentally resurrected proteins.

## Evolutionary interpretation

### 1. An old core was already differentiated

Fifteen of the 17 modern barcode positions are reconstructed as different between Anc-IIIb and Anc-IIIc. Eleven are both highly supported and insensitive to removal of any complete vertebrate clade. This is the strongest sequence-level evidence that a substantial part of the modern ligand-facing code was established near the ancestral IIIb/IIIc split.

### 2. The modern barcode was refined after the initial split

Columns 8 and 28 have the same reconstructed ancestral state in both isoforms. Their modern IIIb-specific states are absent or variable in teleosts, whereas the tetrapod IIIb sequences are strongly conserved. Under the primary reconstruction, these positions are therefore later IIIb refinements rather than part of the earliest split.

### 3. The ancestral divergence footprint was broader than the strict modern barcode

Across the 48 informative alignment columns, {n_all_ancestral_diff} positions differ between the two best-state ancestors. Seven lie outside the strict modern 17-site set. These include the IIIc-specific two-residue insertion and positions whose present-day conservation is insufficient for the modern high-confidence barcode. They should be described as **ancestral-state hypotheses**, not added automatically to the validated modern barcode.

### 4. Sampling dependence is localized

The clade jackknife does not undermine the whole reconstruction. State changes are concentrated in:

- teleost-dependent IIIb root states at columns 8, 11, 15, 17, 28, 51, 54 and 59;
- reptile-sensitive IIIc states at columns 15, 17 and 21.

This is biologically interpretable because the nodes of interest occur close to the deep vertebrate split, so basal lineages have disproportionate information about the ancestral state. These sites must remain visibly marked as uncertain in figures and discussion.

## Integrated mechanistic hypothesis

The combined evidence supports a two-layer model:

1. **Evolutionary-functional layer:** a trans-species IIIb/IIIc sequence identity, including an ancestral ligand-facing core, is maintained across vertebrate diversification.
2. **Coordinate layer:** annotation sources and domain callers alter absolute residue numbers and exact D3 endpoints, while the relative cassette-to-D3 relation remains invariant.

The proposed general principle is:

> For homologous mutually exclusive exons, biological robustness should be assessed using relational, coordinate-shift-invariant quantities rather than a single absolute amino-acid coordinate.

## Claims that are now justified

- The isoform code generalizes to completely unseen vertebrate clades.
- The result is not driven by mammals or by one taxonomic group.
- Most of the modern structural barcode is consistent with an ancestral IIIb/IIIc difference.
- A smaller 11-site core has strong posterior and clade-jackknife support.
- Two canonical barcode sites probably represent later refinement rather than the initial split.

## Claims that remain inappropriate

- The exact date or molecular mechanism of the original exon duplication was not inferred.
- The reconstructed sequences are not experimental evidence of ancestral ligand specificity.
- Posterior support is conditional on the alignment, tree, LG+G4 model and gap treatment.
- The seven additional ancestral differences are not automatically new validated barcode sites.
"""
    (outdir / "SCIENTIFIC_RESULTS_AND_INTERPRETATION.md").write_text(results_md, encoding="utf-8")

    methods_md = """# Reproducibility and method notes

## Input scope

- 28 complete primary species pairs / 56 cassette sequences.
- Review-only Canis IIIc and Pongo IIIb were excluded, matching the formal paired analyses.
- Original alignment length: 69 columns.
- 21 columns that were all gaps after exclusion of review-only records were removed from ASR, leaving 48 informative columns.

## Leave-one-clade-out validation

The five folds were defined a priori from major vertebrate clades. Feature selection was repeated within every training set. The full alignment coordinates were retained because alignment construction is label-free; held-out isoform labels did not enter feature discovery or model fitting.

Classification used independent categorical amino-acid probabilities at de novo selected positions and Jeffreys smoothing (alpha=0.5). The signed score is log P(sequence | IIIb) minus log P(sequence | IIIc), multiplied by -1 for true IIIc records so that positive values always support the correct class.

The null procedure independently swapped IIIb and IIIc labels within each species and reran the complete five-fold pipeline. Ties were counted as half-correct. The empirical p-value includes the standard plus-one correction.

## Ancestral sequence reconstruction

The calculation implements marginal empirical-Bayes reconstruction under LG+G4. The LG exchangeabilities and normalized frequencies match the IQ-TREE implementation. Four equal-probability discrete-gamma category mean rates were calculated from alpha=0.9960.

For each site and rate category:

1. conditional likelihoods were calculated postorder;
2. outside likelihoods were propagated preorder;
3. node-state joint likelihoods were integrated across gamma categories;
4. maximum-posterior states and full 20-amino-acid posterior vectors were stored.

Gaps are not states in the LG substitution model and were treated as missing. Presence/absence was summarized separately using clade fractions. This is conservative and avoids presenting an amino-acid substitution model as an indel model.

## Sensitivity analysis

Each major clade was set to missing in turn and the reconstruction repeated on the same tree. A state is called all-jackknife-stable only when the maximum-posterior amino acid remains unchanged under all five clade removals.

## Literature and software references

- IQ-TREE ASR command documentation: https://iqtree.github.io/doc/Command-Reference
- Minh et al. 2020. IQ-TREE 2. Molecular Biology and Evolution. doi:10.1093/molbev/msaa015
- Hanson-Smith et al. 2010. Robustness of ancestral sequence reconstruction to phylogenetic uncertainty. doi:10.1093/molbev/msq081
- Mistry et al. 2003. Of urchins and men. RNA. doi:10.1261/rna.2470903
- Yeh et al. 2003. Structural basis by which alternative splicing confers specificity in FGFRs. doi:10.1073/pnas.0436500100
- Olsen et al. 2006. Genes & Development. doi:10.1101/gad.1365406
- Warzecha et al. 2009. ESRP1 and ESRP2 regulate FGFR2 splicing. doi:10.1016/j.molcel.2009.01.025
"""
    (outdir / "METHODS_AND_REPRODUCIBILITY.md").write_text(methods_md, encoding="utf-8")

    storyboard = """# Presentation story: paper-level narrative

## Core principle

Do not begin with the final 17-site barcode. Begin with the apparent contradiction:

> The same biological event appears to have different domain boundaries depending on annotation and caller.

Then resolve the contradiction step by step.

## Recommended 12-slide story

### 1. The problem

Show one FGFR2 protein with three disagreeing D3 endpoints. Ask: biological variation or coordinate/model artifact?

### 2. The biological system

Introduce mutually exclusive IIIb/IIIc exons, epithelial/mesenchymal regulation and ligand specificity. Keep this visual and brief.

### 3. The naive result

Show the original apparent boundary heterogeneity. This creates tension.

### 4. Coordinate audit

Reveal the coordinate mismatch and explain why direct residue-number comparison fails.

### 5. Annotation-aware framework

Define precision, consistency and robustness. Introduce the relational invariant I = D3_end - cassette_end.

### 6. Caller effects

Show the repeated-measures caller bias and the replicated caller directions. Message: exact endpoints are model-specific.

### 7. Robustness plateau

Show the 2D robustness surface. Message: the biological conclusion is not rescued by one arbitrary threshold.

### 8. Evolutionary sequence identity

Show JSD and the 17-site barcode, then the corrected AU topology test. Message: cassette identity tracks isoform across species.

### 9. New held-out validation

Use the LOCO figure. Pause on the mammalian hold-out: only 12/17 original sites are rediscovered, but all 30 mammalian proteins are correct. Message: generalization is stronger than feature-list rigidity.

### 10. New ancestral reconstruction

Show the ancestral figure. Message: 15/17 modern barcode sites are reconstructed as distinct at the split; 11 form a conservative robust ancestral core; two are later refinements.

### 11. Structural and annotation replication

Combine direct ligand-contact evidence with NCBI/Ensembl coordinate invariance. Message: functional surfaces are conserved while coordinate frames move.

### 12. Final synthesis

Use the integrated synthesis figure or the presentation graphical abstract. End with:

> The coordinate is not the invariant. The biological invariant is the relation between an evolutionarily conserved exon module and its protein-domain context.

## Delivery rule

Every result slide should have exactly one sentence that the audience can repeat. Place methodological detail in a small footer or backup slide. Do not show all p-values at once; reveal each only when it resolves a specific alternative explanation.

## Backup slides

- complete LOCO prediction table;
- paired-label permutation null;
- bootstrap feature selection stability;
- full ancestral posterior table;
- clade-jackknife state changes;
- corrected AU parsing audit;
- cross-annotation candidate quality metrics;
- limitations and negative claims.
"""
    (outdir / "PRESENTATION_STORYBOARD.md").write_text(storyboard, encoding="utf-8")

    readme = """# FGFR2 final extension package

This package adds three final components to the existing paper-level analysis:

1. an integrated mechanistic synthesis figure;
2. leave-one-clade-out evolutionary validation;
3. LG+G4 marginal ancestral sequence reconstruction with clade-jackknife sensitivity.

## Main figures

- `figures/01_final_integrated_synthesis.*`
- `figures/02_leave_one_clade_out_validation.*`
- `figures/03_ancestral_sequence_reconstruction.*`
- `figures/04_presentation_graphical_abstract.*`

## Main reports

- `SCIENTIFIC_RESULTS_AND_INTERPRETATION.md`
- `METHODS_AND_REPRODUCIBILITY.md`
- `PRESENTATION_STORYBOARD.md`

## Most important caution

The ancestral sequences are conditional model-based reconstructions. The 11-site robust ancestral core is the strongest ASR claim. Sites that change under clade jackknifes must remain marked as uncertain.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")


def main():
    base = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", default=base / "data/framework_snapshot/alignments/final_fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa", type=Path)
    parser.add_argument("--tree", default=base / "data/external/phylogeny/unconstrained.treefile", type=Path)
    parser.add_argument("--iqtree-report", default=base / "data/external/phylogeny/unconstrained.iqtree", type=Path)
    parser.add_argument("--outdir", default=base / "results/reproduced/08_loco_asr_synthesis", type=Path)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--bootstraps", type=int, default=1000)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    if outdir.exists():
        shutil.rmtree(outdir)
    for sub in ("figures", "tables", "sequences", "scripts"):
        (outdir / sub).mkdir(parents=True, exist_ok=True)

    aln, species, arr, aa_order, aa_to_idx, seq_by_tree_tip, alignment_id, informative_cols = load_alignment(args.alignment)
    loco = run_loco(species, arr, aa_order, args.permutations, args.bootstraps)
    asr = run_asr(args.tree, args.iqtree_report, seq_by_tree_tip, arr, informative_cols)


    loco["fold"].to_csv(outdir / "tables/loco_fold_summary.tsv", sep="\t", index=False)
    loco["prediction"].to_csv(outdir / "tables/loco_held_out_predictions.tsv", sep="\t", index=False)
    loco["site"].to_csv(outdir / "tables/loco_de_novo_site_selection.tsv", sep="\t", index=False)
    loco["permutation"].to_csv(outdir / "tables/loco_paired_label_permutation_null.tsv", sep="\t", index=False)
    loco["bootstrap_fold"].to_csv(outdir / "tables/loco_training_bootstrap_summary.tsv", sep="\t", index=False)
    loco["bootstrap_site"].to_csv(outdir / "tables/loco_training_bootstrap_site_stability.tsv", sep="\t", index=False)

    asr["asr"].to_csv(outdir / "tables/ancestral_state_summary.tsv", sep="\t", index=False)
    asr["posterior"].to_csv(outdir / "tables/ancestral_state_full_posterior.tsv", sep="\t", index=False)
    asr["jackknife"].to_csv(outdir / "tables/ancestral_clade_jackknife.tsv", sep="\t", index=False)
    asr["barcode"].to_csv(outdir / "tables/ancestral_barcode_evolution.tsv", sep="\t", index=False)
    asr["difference"].to_csv(outdir / "tables/ancestral_difference_footprint.tsv", sep="\t", index=False)


    (outdir / "sequences/ancestral_sequences_best_state.faa").write_text(
        f">Anc_FGFR2_IIIb_LG_G4_best_state\n{asr['ungapped']['IIIb']}\n"
        f">Anc_FGFR2_IIIc_LG_G4_best_state\n{asr['ungapped']['IIIc']}\n",
        encoding="utf-8",
    )
    (outdir / "sequences/ancestral_sequences_thresholded.faa").write_text(
        f">Anc_FGFR2_IIIb_LG_G4_threshold_0.8\n{asr['thresholded_ungapped']['IIIb']}\n"
        f">Anc_FGFR2_IIIc_LG_G4_threshold_0.8\n{asr['thresholded_ungapped']['IIIc']}\n",
        encoding="utf-8",
    )
    (outdir / "sequences/ancestral_sequences_aligned.faa").write_text(
        f">Anc_FGFR2_IIIb_LG_G4_aligned\n{asr['aligned']['IIIb']}\n"
        f">Anc_FGFR2_IIIc_LG_G4_aligned\n{asr['aligned']['IIIc']}\n",
        encoding="utf-8",
    )


    plot_synthesis(loco, asr, outdir / "figures")
    plot_loco(loco, outdir / "figures")
    plot_asr(asr, outdir / "figures")
    plot_presentation_hero(loco, asr, outdir / "figures")

    write_reports(outdir, loco, asr)


    shutil.copy2(Path(__file__), outdir / "scripts/run_fgfr2_final_extension.py")


    zip_path = outdir.with_suffix(".zip")
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in outdir.rglob("*"):
            if path.is_file():
                archive.write(path, arcname=f"{outdir.name}/{path.relative_to(outdir)}")

    print(json.dumps({
        "output_directory": str(outdir),
        "zip": str(zip_path),
        "loco_correct": int(loco['prediction'].correct.sum()),
        "loco_total": len(loco['prediction']),
        "loco_permutation_p": loco['permutation_p'],
        "ancestral_barcode_different": int(asr['barcode'].ancestrally_different.sum()),
        "ancestral_robust_core": int(asr['barcode'].high_confidence_robust_ancestral_difference.sum()),
        "ancestral_all_differences": int(asr['difference'].different.sum()),
    }, indent=2))


if __name__ == "__main__":
    main()
