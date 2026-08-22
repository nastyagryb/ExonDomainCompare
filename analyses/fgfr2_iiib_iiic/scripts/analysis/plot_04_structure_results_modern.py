
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from publication_style import (
    COLORS, apply_publication_style, clean_axis,
    save_figure, takehome, title_block,
)

def parse_args():
    p = argparse.ArgumentParser(description="Create modern individual structure-interface plots.")
    p.add_argument("--results-dir", required=True)
    return p.parse_args()

def p_text(value):
    if value is None or not np.isfinite(value):
        return "not estimable"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3f}"

def sensitivity_text(primary, sensitivity):
    return (
        f"control-quartile p = {p_text(primary)}; "
        f"universe-quartile p = {p_text(sensitivity)}"
    )

def main():
    args = parse_args()
    apply_publication_style()
    out = Path(args.results_dir)
    summary_path = out / "structure_interface_enrichment_summary.tsv"
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing {summary_path}; run structure mapping first.")
    summary = pd.read_csv(summary_path, sep="\t")

    for _, srow in summary.iterrows():
        pdb_id = str(srow["pdb_id"])
        table = pd.read_csv(out / f"{pdb_id}_residue_interface_metrics.tsv", sep="\t")
        target = table["is_discriminating"].astype(bool)


        fig, ax = plt.subplots(figsize=(12.4, 7.4))
        fig.subplots_adjust(top=0.79, left=0.11, right=0.95, bottom=0.25)
        ax.plot(table["resi"], table["min_ligand_distance_A"],
                color=COLORS["blue"], linewidth=2.0, alpha=0.75)
        ax.scatter(table.loc[~target, "resi"],
                   table.loc[~target, "min_ligand_distance_A"],
                   s=34, color="#C8D1DB", edgecolor="white", linewidth=0.4,
                   label="Other cassette residues")
        ax.scatter(table.loc[target, "resi"],
                   table.loc[target, "min_ligand_distance_A"],
                   s=95, marker="*", color=COLORS["gold"],
                   edgecolor=COLORS["navy"], linewidth=0.7,
                   label="Isoform-discriminating residues", zorder=4)
        ax.axhspan(0, float(srow["direct_distance_A"]),
                   color=COLORS["teal"], alpha=0.12,
                   label=f'Direct-contact zone ≤{srow["direct_distance_A"]:g} Å')
        ax.axhspan(float(srow["direct_distance_A"]),
                   float(srow["near_distance_A"]),
                   color=COLORS["gold"], alpha=0.10,
                   label=f'Near-interface zone ≤{srow["near_distance_A"]:g} Å')
        ax.set_xlabel("PDB receptor residue number")
        ax.set_ylabel("Minimum heavy-atom distance to ligand (Å)")
        clean_axis(ax, "both")
        ax.legend(frameon=False, ncol=2, loc="upper right")
        title_block(
            fig,
            f"{pdb_id}: isoform-discriminating residues at the ligand interface",
            "Minimum heavy-atom distances quantify direct and near-interface positioning across the alternative D3 segment."
        )
        takehome(
            ax,
            f'{int(srow["n_discriminating_direct_contacts"])} direct contacts and '
            f'{int(srow["n_discriminating_near_interface"])} near-interface residues among '
            f'{int(srow["n_discriminating_residues_mapped"])} mapped barcode positions.\n'
            "Contact enrichment is retained under both common-bin specifications: "
            f'{sensitivity_text(srow["matched_permutation_p_direct_contacts"], srow["sensitivity_p_direct_contacts_universe_quartiles"])}.',
            loc=(0.98, -0.28),
        )
        save_figure(fig, out, f"modern_{pdb_id}_ligand_distance_profile")


        fig, ax = plt.subplots(figsize=(12.4, 7.4))
        fig.subplots_adjust(top=0.79, left=0.11, right=0.95, bottom=0.25)
        positive = table["delta_sasa_A2"].clip(lower=0)
        ax.bar(table["resi"], positive, width=0.85,
               color=np.where(target, COLORS["gold"], "#BFD0DF"),
               edgecolor="white", linewidth=0.25)
        ax.scatter(table.loc[target, "resi"],
                   table.loc[target, "delta_sasa_A2"],
                   s=85, marker="*", color=COLORS["orange"],
                   edgecolor=COLORS["navy"], linewidth=0.7, zorder=4)
        ax.set_xlabel("PDB receptor residue number")
        ax.set_ylabel("Buried solvent-accessible surface area, ΔSASA (Å²)")
        clean_axis(ax, "y")
        title_block(
            fig,
            f"{pdb_id}: ligand binding buries selected barcode residues",
            "ΔSASA is calculated as receptor-alone SASA minus complex SASA for each observed cassette residue."
        )
        takehome(
            ax,
            f'Summed barcode ΔSASA = {srow["sum_discriminating_delta_sasa_A2"]:.1f} Å².\n'
            "Exposure-matched result: "
            f'{sensitivity_text(srow["matched_permutation_p_sum_delta_sasa"], srow["sensitivity_p_sum_delta_sasa_universe_quartiles"])}; '
            "interpret specification-sensitive results cautiously.",
            loc=(0.98, -0.28),
        )
        save_figure(fig, out, f"modern_{pdb_id}_delta_SASA_profile")


    fig, ax = plt.subplots(figsize=(10.8, 7.0))
    fig.subplots_adjust(top=0.79, left=0.16, right=0.95, bottom=0.16)
    labels = summary["pdb_id"].astype(str) + " · " + summary["isoform"].astype(str)
    y = np.arange(len(summary))
    mapped = summary["n_discriminating_residues_mapped"].to_numpy(float)
    direct_pct = 100 * summary["n_discriminating_direct_contacts"].to_numpy(float) / mapped
    near_pct = 100 * summary["n_discriminating_near_interface"].to_numpy(float) / mapped
    ax.barh(y, near_pct, color=COLORS["gold"], alpha=0.55, label="Near interface")
    ax.barh(y, direct_pct, color=COLORS["teal"], alpha=0.95, label="Direct contact")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontweight="bold")
    ax.set_xlabel("Mapped barcode residues (%)")
    ax.set_xlim(0, 105)
    clean_axis(ax, "x")
    ax.legend(frameon=False)
    title_block(
        fig,
        "Structure complexes independently test whether the molecular barcode localizes to ligand-contacting surfaces",
        "Direct-contact enrichment persists under both shared-bin matching schemes; ΔSASA sensitivity is reported separately."
    )
    save_figure(fig, out, "modern_structure_interface_summary")

if __name__ == "__main__":
    main()
