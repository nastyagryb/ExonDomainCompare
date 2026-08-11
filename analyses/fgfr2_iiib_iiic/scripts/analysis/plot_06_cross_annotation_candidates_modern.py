
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from publication_style import (
    COLORS, apply_publication_style, clean_axis,
    save_figure, takehome, title_block,
)

SOURCE_COLORS = {"ncbi": COLORS["blue"], "ensembl": COLORS["orange"]}

def parse_args():
    p = argparse.ArgumentParser(description="Create modern cross-annotation candidate plots.")
    p.add_argument("--results-dir", required=True)
    return p.parse_args()

def main():
    args = parse_args()
    apply_publication_style()
    out = Path(args.results_dir)
    selected_path = out / "selected_cross_annotation_candidates.tsv"
    if not selected_path.exists():
        raise FileNotFoundError(f"Missing {selected_path}; run cross-annotation acquisition first.")
    selected = pd.read_csv(selected_path, sep="\t")


    plot = selected.sort_values(["species", "predicted_isoform", "source"]).copy()
    plot["label"] = (
        plot["species"].str.replace("_", " ", regex=False)
        + " · " + plot["predicted_isoform"]
        + " · " + plot["source"].str.upper()
    )
    fig, ax = plt.subplots(figsize=(12.8, max(7.0, 0.36*len(plot))))
    fig.subplots_adjust(top=0.84, left=0.29, right=0.94, bottom=0.12)
    y = np.arange(len(plot))
    for source, group in plot.groupby("source"):
        idx = group.index
        positions = [plot.index.get_loc(i) for i in idx]
        ax.scatter(group["score_margin"], positions,
                   s=70 + 110*group["reference_coverage"],
                   color=SOURCE_COLORS.get(source, COLORS["gray"]),
                   edgecolor="white", linewidth=0.8,
                   label=source.upper())
    ax.axvline(0, color=COLORS["navy"], lw=1.2, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(plot["label"])
    ax.invert_yaxis()
    ax.set_xlabel("Best-versus-second-best cassette alignment score margin")
    clean_axis(ax, "x")
    ax.legend(frameon=False)
    title_block(
        fig,
        "NCBI and Ensembl candidates are classified by independent alignment evidence",
        "Point size represents coverage of the curated human IIIb/IIIc cassette reference."
    )
    takehome(
        ax,
        "Large positive margins and high reference coverage indicate unambiguous isoform classification."
    )
    save_figure(fig, out, "modern_cross_annotation_candidate_confidence")


    matrix = selected.pivot_table(
        index=["species", "predicted_isoform"],
        columns="source",
        values="reference_coverage",
        aggfunc="first"
    ).sort_index()
    fig, ax = plt.subplots(figsize=(9.5, max(6.5, 0.40*len(matrix))))
    fig.subplots_adjust(top=0.84, left=0.30, right=0.88, bottom=0.13)
    image = ax.imshow(matrix.to_numpy()*100, aspect="auto", vmin=0, vmax=100,
                      cmap="YlGnBu")
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels([x.upper() for x in matrix.columns], fontweight="bold")
    labels = [f"{sp.replace('_',' ')} · {iso}" for sp, iso in matrix.index]
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iloc[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value*100:.0f}%",
                        ha="center", va="center",
                        fontsize=10, fontweight="bold",
                        color="white" if value > 0.60 else COLORS["navy"])
    cbar = fig.colorbar(image, ax=ax, pad=0.025)
    cbar.set_label("Cassette reference coverage (%)")
    title_block(
        fig,
        "Reference-cassette coverage exposes source-specific annotation incompleteness",
        "Rows compare the selected NCBI and Ensembl candidate for every species–isoform combination."
    )
    save_figure(fig, out, "modern_cross_annotation_reference_coverage")


    comparison_path = out / "ncbi_vs_ensembl_coordinate_comparison.tsv"
    if comparison_path.exists():
        comparison = pd.read_csv(comparison_path, sep="\t")
        if not comparison.empty:
            comparison = comparison.sort_values(["species", "isoform"]).copy()
            labels = (
                comparison["species"].str.replace("_", " ", regex=False)
                + " · " + comparison["isoform"]
            )
            fig, ax = plt.subplots(figsize=(12.2, max(6.5, 0.40*len(comparison))))
            fig.subplots_adjust(top=0.84, left=0.29, right=0.94, bottom=0.13)
            y = np.arange(len(comparison))
            ax.scatter(comparison["cassette_start_delta_ncbi_minus_ensembl"],
                       y, s=75, color=COLORS["blue"],
                       edgecolor="white", linewidth=0.7,
                       label="Cassette start")
            ax.scatter(comparison["cassette_end_delta_ncbi_minus_ensembl"],
                       y, s=75, marker="D", color=COLORS["orange"],
                       edgecolor="white", linewidth=0.7,
                       label="Cassette end")
            ax.axvline(0, color=COLORS["navy"], lw=1.2, ls="--")
            ax.set_yticks(y)
            ax.set_yticklabels(labels)
            ax.invert_yaxis()
            ax.set_xlabel("NCBI minus Ensembl cassette coordinate (aa)")
            clean_axis(ax, "x")
            ax.legend(frameon=False, ncol=2)
            title_block(
                fig,
                "Independent annotation sources can disagree at residue level before domain integration",
                "The next InterProScan step determines whether coordinate differences alter the final exon–domain topology."
            )
            save_figure(fig, out, "modern_cross_annotation_sequence_coordinate_deltas")

if __name__ == "__main__":
    main()
