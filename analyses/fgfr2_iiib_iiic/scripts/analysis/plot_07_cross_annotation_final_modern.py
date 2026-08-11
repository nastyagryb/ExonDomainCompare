
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

def parse_args():
    p = argparse.ArgumentParser(description="Create modern final cross-annotation plots.")
    p.add_argument("--results-dir", required=True)
    p.add_argument("--consensus", type=float, default=0.80)
    return p.parse_args()

def main():
    args = parse_args()
    apply_publication_style()
    out = Path(args.results_dir)

    comparison_path = out / "cross_annotation_boundary_replication.tsv"
    protein_path = out / "cross_annotation_topology_by_protein.tsv"
    if not comparison_path.exists() or not protein_path.exists():
        raise FileNotFoundError("Run final cross-annotation integration before plotting.")

    comparison = pd.read_csv(comparison_path, sep="\t")
    per_protein = pd.read_csv(protein_path, sep="\t")


    if not comparison.empty:
        plot = comparison.sort_values(["species", "isoform"]).copy()
        labels = plot["species"].str.replace("_", " ", regex=False) + " · " + plot["isoform"]
        y = np.arange(len(plot))
        fig, ax = plt.subplots(figsize=(12.6, max(7.0, 0.43*len(plot))))
        fig.subplots_adjust(top=0.84, left=0.29, right=0.94, bottom=0.13)
        ax.scatter(plot["cassette_end_delta_ncbi_minus_ensembl"], y,
                   s=80, color=COLORS["blue"], edgecolor="white",
                   linewidth=0.7, label="Cassette end")
        ax.scatter(plot["median_d3_end_delta_ncbi_minus_ensembl"], y,
                   s=90, marker="D", color=COLORS["orange"],
                   edgecolor="white", linewidth=0.7, label="Median D3 end")
        for yi, (_, row) in zip(y, plot.iterrows()):
            ax.plot(
                [row["cassette_end_delta_ncbi_minus_ensembl"],
                 row["median_d3_end_delta_ncbi_minus_ensembl"]],
                [yi, yi], color=COLORS["grid"], linewidth=1.4, zorder=0
            )
        ax.axvline(0, color=COLORS["navy"], lw=1.2, ls="--")
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.invert_yaxis()
        ax.set_xlabel("NCBI minus Ensembl coordinate (aa)")
        clean_axis(ax, "x")
        ax.legend(frameon=False, ncol=2)
        title_block(
            fig,
            "Residue-level annotation differences are separated from topology-level replication",
            "Cassette coordinates and median D3 ends are compared independently for each species–isoform pair."
        )
        agreement = 100 * plot["same_topology_class"].mean()
        takehome(
            ax,
            f"Cross-annotation topology agreement: {agreement:.1f}% "
            f"({int(plot['same_topology_class'].sum())}/{len(plot)} pairs)."
        )
        save_figure(fig, out, "modern_cross_annotation_boundary_deltas")


    plot = per_protein.copy()
    fig, ax = plt.subplots(figsize=(10.4, 8.0))
    fig.subplots_adjust(top=0.82, left=0.13, right=0.94, bottom=0.13)
    colors = np.where(plot["topology_pass"], COLORS["teal"], COLORS["red"])
    markers = {"ncbi": "o", "ensembl": "D"}
    for source, group in plot.groupby("source"):
        ax.scatter(group["start_consensus"]*100,
                   group["end_consensus"]*100,
                   s=85, marker=markers.get(source, "o"),
                   color=np.where(group["topology_pass"], COLORS["teal"], COLORS["red"]),
                   edgecolor="white", linewidth=0.8,
                   label=source.upper())
    threshold = args.consensus * 100
    ax.axvline(threshold, color=COLORS["orange"], lw=1.4, ls="--")
    ax.axhline(threshold, color=COLORS["orange"], lw=1.4, ls="--")
    ax.axvspan(threshold, 101, ymin=threshold/101, ymax=1,
               color=COLORS["teal"], alpha=0.06)
    ax.set_xlim(0, 102)
    ax.set_ylim(0, 102)
    ax.set_xlabel("Cassette-start topology consensus (%)")
    ax.set_ylabel("Cassette-end boundary consensus (%)")
    clean_axis(ax, "both")
    ax.legend(frameon=False)
    title_block(
        fig,
        "Start and end evidence jointly determine cross-annotation topology support",
        "The upper-right quadrant contains proteins passing both annotation-aware consensus criteria."
    )
    takehome(
        ax,
        f"{int(plot['topology_pass'].sum())}/{len(plot)} selected proteins satisfy both topology criteria."
    )
    save_figure(fig, out, "modern_cross_annotation_consensus_scatter")


    if not comparison.empty:
        plot = comparison.sort_values(["species", "isoform"]).copy()
        plot["label"] = plot["species"].str.replace("_", " ", regex=False) + " · " + plot["isoform"]
        values = np.column_stack([
            plot["ncbi_topology_pass"].astype(int),
            plot["ensembl_topology_pass"].astype(int),
            plot["same_topology_class"].astype(int),
        ])
        fig, ax = plt.subplots(figsize=(9.8, max(6.8, 0.43*len(plot))))
        fig.subplots_adjust(top=0.84, left=0.31, right=0.90, bottom=0.13)
        cmap = plt.matplotlib.colors.ListedColormap(["#E8EEF4", COLORS["teal"]])
        image = ax.imshow(values, aspect="auto", vmin=0, vmax=1, cmap=cmap)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["NCBI passes", "Ensembl passes", "Same class"],
                           fontweight="bold")
        ax.set_yticks(range(len(plot)))
        ax.set_yticklabels(plot["label"])
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                ax.text(j, i, "✓" if values[i, j] else "×",
                        ha="center", va="center",
                        fontsize=14, fontweight="bold",
                        color="white" if values[i, j] else COLORS["gray"])
        title_block(
            fig,
            "A transparent evidence matrix shows where annotation sources agree or diverge",
            "Each row records NCBI support, Ensembl support and equality of the final topology class."
        )
        save_figure(fig, out, "modern_cross_annotation_topology_matrix")

if __name__ == "__main__":
    main()
