
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO
from statsmodels.stats.multitest import multipletests

from common import ensure_dir, read_tsv, setup_logging, write_json

AA = np.array(list("ARNDCQEGHILKMFPSTWYV"))
AA_TO_IDX = {a: i for i, a in enumerate(AA)}
GAPS = {"-", ".", "?", "X"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute sequence-weighted IIIb/IIIc Jensen-Shannon divergence and paired within-species label-permutation tests."
    )
    p.add_argument("--alignment", required=True, help="Combined IIIb/IIIc cassette protein MSA (FASTA)")
    p.add_argument("--truth-table", help="Optional truth table used to retain main_analysis proteins")
    p.add_argument("--outdir", required=True)
    p.add_argument("--weighting", choices=["henikoff", "uniform"], default="henikoff")
    p.add_argument("--permutations", type=int, default=10000)
    p.add_argument("--seed", type=int, default=20260804)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def parse_record_id(record_id: str) -> tuple[str, str, str]:
    parts = record_id.split("|")
    if len(parts) < 3:
        raise ValueError(f"Expected species|isoform|protein_id FASTA ID, got: {record_id}")
    return parts[0], parts[1], parts[2]


def henikoff_weights(chars: np.ndarray) -> np.ndarray:
    n, length = chars.shape
    w = np.zeros(n, dtype=float)
    informative = 0
    for j in range(length):
        column = chars[:, j]
        valid = np.array([x not in GAPS and x in AA_TO_IDX for x in column])
        residues = column[valid]
        if len(residues) == 0:
            continue
        values, counts = np.unique(residues, return_counts=True)
        r = len(values)
        count_map = dict(zip(values, counts))
        for i in np.where(valid)[0]:
            w[i] += 1.0 / (r * count_map[column[i]])
        informative += 1
    if informative == 0 or w.sum() == 0:
        return np.ones(n)
    w /= informative
    return w * n / w.sum()


def distributions(onehot: np.ndarray, weights: np.ndarray, labels: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mask_b = labels == 0
    mask_c = labels == 1
    weighted = onehot * weights[:, None, None]
    num_b = weighted[mask_b].sum(axis=0)
    num_c = weighted[mask_c].sum(axis=0)
    den_b = num_b.sum(axis=1, keepdims=True)
    den_c = num_c.sum(axis=1, keepdims=True)
    p = np.divide(num_b, den_b, out=np.zeros_like(num_b), where=den_b > 0)
    q = np.divide(num_c, den_c, out=np.zeros_like(num_c), where=den_c > 0)
    gap_b = 1 - (onehot[mask_b].sum(axis=2) * weights[mask_b, None]).sum(axis=0) / weights[mask_b].sum()
    gap_c = 1 - (onehot[mask_c].sum(axis=2) * weights[mask_c, None]).sum(axis=0) / weights[mask_c].sum()
    return p, q, gap_b, gap_c


def entropy(p: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(p > 0, p * np.log2(p), 0.0)
    return -terms.sum(axis=1)


def jsd(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    m = 0.5 * (p + q)
    with np.errstate(divide="ignore", invalid="ignore"):
        klp = np.where(p > 0, p * np.log2(np.divide(p, m, out=np.ones_like(p), where=m > 0)), 0.0).sum(axis=1)
        klq = np.where(q > 0, q * np.log2(np.divide(q, m, out=np.ones_like(q), where=m > 0)), 0.0).sum(axis=1)
    return 0.5 * (klp + klq)


def main() -> None:
    args = parse_args()
    setup_logging(args.verbose)
    out = ensure_dir(args.outdir)
    aln = AlignIO.read(args.alignment, "fasta")
    records = []
    for rec in aln:
        species, isoform, protein_id = parse_record_id(rec.id)
        records.append({"record": rec, "species": species, "isoform": isoform, "protein_id": protein_id})

    if args.truth_table:
        truth = read_tsv(args.truth_table, required=["species", "final_isoform_label", "protein_id"])
        if "recommended_use_post_rescue" in truth.columns:
            truth = truth[truth["recommended_use_post_rescue"] == "main_analysis"]
        keep = set(zip(truth["species"], truth["final_isoform_label"], truth["protein_id"]))
        records = [r for r in records if (r["species"], r["isoform"], r["protein_id"]) in keep]


    counts = pd.DataFrame(records).groupby(["species", "isoform"]).size().unstack(fill_value=0)
    complete_species = counts.index[(counts.get("IIIb", 0) == 1) & (counts.get("IIIc", 0) == 1)]
    records = [r for r in records if r["species"] in set(complete_species)]
    records.sort(key=lambda r: (r["species"], r["isoform"]))
    if not records:
        raise ValueError("No complete IIIb/IIIc species pairs remain after filtering")

    chars = np.array([list(str(r["record"].seq).upper()) for r in records])
    n_seq, n_col = chars.shape
    weights = henikoff_weights(chars) if args.weighting == "henikoff" else np.ones(n_seq)
    labels = np.array([0 if r["isoform"] == "IIIb" else 1 for r in records], dtype=int)
    species = np.array([r["species"] for r in records])
    protein_ids = np.array([r["protein_id"] for r in records])

    onehot = np.zeros((n_seq, n_col, len(AA)), dtype=float)
    for i in range(n_seq):
        for j in range(n_col):
            idx = AA_TO_IDX.get(chars[i, j])
            if idx is not None:
                onehot[i, j, idx] = 1.0

    p, q, gap_b, gap_c = distributions(onehot, weights, labels)
    obs_jsd = jsd(p, q)
    h_b = entropy(p)
    h_c = entropy(q)
    conservation_b = 1 - h_b / np.log2(len(AA))
    conservation_c = 1 - h_c / np.log2(len(AA))
    gap_penalty = 1 - np.maximum(gap_b, gap_c)
    idi = obs_jsd * conservation_b * conservation_c * gap_penalty

    pairs = []
    for sp in complete_species:
        idx = np.where(species == sp)[0]
        if len(idx) != 2:
            raise ValueError(f"Expected exactly two sequences for {sp}, found {len(idx)}")
        pairs.append(idx)

    rng = np.random.default_rng(args.seed)
    exceed = np.zeros(n_col, dtype=int)
    global_obs = float(obs_jsd.sum())
    global_exceed = 0
    for _ in range(args.permutations):
        perm_labels = labels.copy()
        swaps = rng.integers(0, 2, size=len(pairs), endpoint=False)
        for swap, idx in zip(swaps, pairs):
            if swap:
                perm_labels[idx] = perm_labels[idx[::-1]]
        pp, qq, _, _ = distributions(onehot, weights, perm_labels)
        perm_jsd = jsd(pp, qq)
        exceed += perm_jsd >= (obs_jsd - 1e-12)
        global_exceed += float(perm_jsd.sum()) >= global_obs - 1e-12

    p_emp = (exceed + 1) / (args.permutations + 1)
    q_value = multipletests(p_emp, alpha=args.fdr, method="fdr_bh")[1]
    global_p = (global_exceed + 1) / (args.permutations + 1)

    def major_aa(dist: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        idx = dist.argmax(axis=1)
        frac = dist[np.arange(len(dist)), idx]
        names = AA[idx].astype(object)
        names[dist.sum(axis=1) == 0] = "-"
        return names, frac

    major_b, frac_b = major_aa(p)
    major_c, frac_c = major_aa(q)
    table = pd.DataFrame(
        {
            "alignment_column_1based": np.arange(1, n_col + 1),
            "IIIb_major_aa": major_b,
            "IIIc_major_aa": major_c,
            "IIIb_major_fraction_weighted": frac_b,
            "IIIc_major_fraction_weighted": frac_c,
            "IIIb_entropy_bits": h_b,
            "IIIc_entropy_bits": h_c,
            "IIIb_gap_fraction_weighted": gap_b,
            "IIIc_gap_fraction_weighted": gap_c,
            "weighted_jsd_bits": obs_jsd,
            "isoform_discrimination_index": idi,
            "paired_permutation_p": p_emp,
            "bh_fdr_q": q_value,
            "significant_fdr": q_value <= args.fdr,
        }
    )
    table.to_csv(out / "weighted_jsd_positions.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "species": species,
            "isoform": [r["isoform"] for r in records],
            "protein_id": protein_ids,
            "sequence_weight": weights,
        }
    ).to_csv(out / "sequence_weights.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = table["alignment_column_1based"]
    ax.plot(x, table["weighted_jsd_bits"], label="Weighted JSD")
    ax.plot(x, table["isoform_discrimination_index"], label="IDI")
    sig = table["significant_fdr"]
    ax.scatter(x[sig], table.loc[sig, "weighted_jsd_bits"], marker="o", label=f"FDR≤{args.fdr:g}")
    ax.set_xlabel("Cassette MSA column")
    ax.set_ylabel("Score")
    ax.set_title("Sequence-weighted isoform divergence with paired label permutation")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "figure_weighted_jsd_permutation.png", dpi=300)
    fig.savefig(out / "figure_weighted_jsd_permutation.svg")
    plt.close(fig)

    write_json(
        {
            "n_species_pairs": int(len(complete_species)),
            "n_sequences": int(n_seq),
            "n_alignment_columns": int(n_col),
            "weighting": args.weighting,
            "n_permutations": args.permutations,
            "global_sum_jsd": global_obs,
            "global_paired_permutation_p": global_p,
            "n_fdr_significant_positions": int(table["significant_fdr"].sum()),
            "fdr_threshold": args.fdr,
        },
        out / "weighted_jsd_summary.json",
    )


if __name__ == "__main__":
    main()
