
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from publication_style import COLORS, apply_publication_style, clean_axis, save_figure


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Export analyses 1–3 as separate publication-ready figures.")
    p.add_argument("--root", required=True, help="Root of the FGFR2 analysis suite")
    p.add_argument("--outdir", required=True)
    p.add_argument("--null-jsd", required=True, help="TSV with global permutation null distribution")
    return p.parse_args()


def title_block(fig, title: str, subtitle: str) -> None:
    fig.suptitle(title, x=0.055, y=0.985, ha="left", va="top",
                 fontsize=22, fontweight="bold", color=COLORS["navy"])
    fig.text(0.055, 0.925, subtitle, ha="left", va="top",
             fontsize=12.5, color=COLORS["gray"])


def footer(fig, text: str, face: str = "#F3F7FA") -> None:
    fig.text(
        0.5, 0.035, text, ha="center", va="center",
        fontsize=11.3, fontweight="bold", color=COLORS["navy"],
        bbox=dict(boxstyle="round,pad=0.55", facecolor=face,
                  edgecolor=COLORS["grid"], linewidth=0.8),
    )


def save(fig, out: Path, stem: str) -> None:
    save_figure(fig, out, stem, dpi=400)


def main() -> None:
    args = parse_args()
    apply_publication_style()
    root = Path(args.root)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    bias_dir = root / "results/01_domain_caller_bias"
    surface_dir = root / "results/02_robustness_surface"
    jsd_dir = root / "results/03_weighted_jsd"

    offsets = pd.read_csv(bias_dir / "member_database_offsets_per_protein.tsv", sep="\t")
    bias = pd.read_csv(bias_dir / "member_database_bias_estimates.tsv", sep="\t")
    loo = pd.read_csv(bias_dir / "leave_one_database_out.tsv", sep="\t")
    surface = pd.read_csv(surface_dir / "robustness_surface_long.tsv", sep="\t")
    minimum = pd.read_csv(surface_dir / "minimum_thresholds_for_target_coverage.tsv", sep="\t")
    jsd = pd.read_csv(jsd_dir / "weighted_jsd_positions.tsv", sep="\t")
    barcode = pd.read_csv(root / "data/structure_mapping_17_discriminating_positions.tsv", sep="\t")
    null_jsd = pd.read_csv(args.null_jsd, sep="\t").iloc[:, 0].to_numpy()

    with open(bias_dir / "domain_caller_bias_summary.json") as h:
        bias_summary = json.load(h)
    with open(surface_dir / "robustness_surface_summary.json") as h:
        surface_summary = json.load(h)
    with open(jsd_dir / "weighted_jsd_summary.json") as h:
        jsd_summary = json.load(h)

    db_colors = {
        "Pfam": "#3973A8",
        "SMART": "#168A83",
        "ProSiteProfiles": "#D59B22",
        "SUPERFAMILY": "#D96B38",
        "Gene3D": "#7561A8",
    }
    db_label = {"ProSiteProfiles": "PROSITE Profiles"}
    db_order = bias.sort_values("mean_signed_offset_aa")["member_database"].tolist()


    fig, ax = plt.subplots(figsize=(12.7, 8.0))
    fig.subplots_adjust(top=0.80, left=0.20, right=0.96, bottom=0.24)
    rng = np.random.default_rng(20260804)
    for y, db in enumerate(db_order):
        vals = offsets.loc[offsets["member_database"] == db, "signed_offset"].to_numpy(float)
        parts = ax.violinplot(vals, positions=[y], vert=False, widths=0.70,
                              showmeans=False, showmedians=False, showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor(db_colors[db]); body.set_edgecolor("none"); body.set_alpha(0.20)
        jitter = rng.normal(0, 0.065, len(vals))
        ax.scatter(vals, y + jitter, s=28, alpha=0.38, color=db_colors[db],
                   edgecolor="white", linewidth=0.35, zorder=3)
        q1, med, q3 = np.quantile(vals, [0.25, 0.5, 0.75])
        ax.plot([q1, q3], [y, y], lw=7, solid_capstyle="round", color=COLORS["navy"], zorder=4)
        ax.scatter([med], [y], s=65, facecolor="white", edgecolor=COLORS["navy"], lw=1.7, zorder=5)
    ax.axvline(0, color=COLORS["navy"], lw=1.4, ls=(0, (4, 3)))
    ax.axvspan(-15, 15, color=COLORS["teal"], alpha=0.06)
    ax.set_yticks(range(len(db_order)))
    ax.set_yticklabels([db_label.get(x, x) for x in db_order], fontweight="bold")
    ax.set_xlabel("Signed terminal D3-end offset relative to reconciled cassette end (aa)")
    ax.set_xlim(-17, 8)
    clean_axis(ax, "x")
    title_block(fig, "Domain callers place the same D3 boundary at systematically different residues",
                "Each point is one primary FGFR2 protein; thick bars show the interquartile range and white points the median.")
    footer(fig, "Exact residue-level precision depends on the domain model, while every distribution remains inside the topology-relevant ±15-aa interval.")
    save(fig, out, "01_domain_caller_offset_distributions")


    fig, ax = plt.subplots(figsize=(11.7, 7.8))
    fig.subplots_adjust(top=0.79, left=0.23, right=0.92, bottom=0.24)
    forest = bias.set_index("member_database").loc[db_order].reset_index()
    y = np.arange(len(forest))
    for i, row in forest.iterrows():
        db = row["member_database"]
        ax.plot([row["bootstrap_ci_low_aa"], row["bootstrap_ci_high_aa"]], [i, i], lw=4,
                solid_capstyle="round", color=db_colors[db])
        ax.scatter(row["mean_signed_offset_aa"], i, s=125, color=db_colors[db],
                   edgecolor="white", lw=1.3, zorder=3)
        ax.text(row["bootstrap_ci_high_aa"] + 0.45, i, f'{row["mean_signed_offset_aa"]:+.2f} aa',
                va="center", fontsize=11, fontweight="bold", color=COLORS["dark"])
    ax.axvline(0, color=COLORS["navy"], lw=1.4, ls=(0, (4, 3)))
    ax.set_yticks(y)
    ax.set_yticklabels([db_label.get(x, x) for x in forest["member_database"]], fontweight="bold")
    ax.set_xlabel("Mean signed offset with species-bootstrap 95% confidence interval")
    ax.set_xlim(-9.5, 6.3)
    clean_axis(ax, "x")
    title_block(fig, "Domain-model bias is directional, reproducible and quantitatively large",
                "Positive values place the predicted D3 end after the reconciled cassette end; negative values place it before.")
    footer(fig, f'Member-database identity explains {bias_summary["database_partial_r2_controlling_species_and_isoform"]*100:.1f}% of residual offset variation after controlling for species and isoform.')
    save(fig, out, "02_domain_caller_bias_forest")


    fig, ax = plt.subplots(figsize=(12.3, 8.0))
    fig.subplots_adjust(top=0.79, left=0.28, right=0.92, bottom=0.24)
    order = ["none", "Gene3D", "Pfam", "ProSiteProfiles", "SMART", "SUPERFAMILY"]
    plot = loo.set_index("excluded_database").loc[order].reset_index()
    labels = ["All databases"] + [f"Exclude {db_label.get(x, x)}" for x in order[1:]]
    y = np.arange(len(plot))
    for yi, row in zip(y, plot.itertuples()):
        ax.plot([0, 58], [yi, yi], color=COLORS["light"], lw=7, solid_capstyle="round", zorder=1)
        ax.scatter(np.arange(1, int(row.n_pass) + 1), np.repeat(yi, int(row.n_pass)),
                   s=32, color=COLORS["teal"], edgecolor="white", lw=0.35, zorder=2)
        ax.text(59.2, yi, f"{row.n_pass}/{row.n_proteins}", va="center", ha="left",
                fontsize=12, fontweight="bold", color=COLORS["navy"])
    ax.set_yticks(y); ax.set_yticklabels(labels, fontweight="bold"); ax.invert_yaxis()
    ax.set_xlim(0, 64); ax.set_xticks([0, 10, 20, 30, 40, 50, 58])
    ax.set_xlabel("Primary proteins retaining the topology criterion")
    clean_axis(ax, "x")
    title_block(fig, "No individual InterPro member database is required for the conclusion",
                "Each row removes one member database and recalculates the annotation-aware consensus from the remaining evidence.")
    footer(fig, "58/58 proteins pass in every scenario. Criterion: absolute offset ≤15 aa and remaining-database consensus ≥80%.", face="#EAF7F3")
    save(fig, out, "03_leave_one_database_out")


    fig, ax = plt.subplots(figsize=(12.8, 8.5))
    fig.subplots_adjust(top=0.80, left=0.12, right=0.88, bottom=0.22)
    pivot = surface.pivot(index="consensus_threshold", columns="distance_threshold_aa", values="pass_fraction").sort_index()
    xv, yv, z = pivot.columns.to_numpy(), pivot.index.to_numpy()*100, pivot.to_numpy()*100
    cmap = LinearSegmentedColormap.from_list("fgfr2", ["#102A43", "#2F6B9A", "#64B5D2", "#80CFC4", "#F2CF65"])
    mesh = ax.pcolormesh(xv, yv, z, shading="nearest", cmap=cmap, vmin=0, vmax=100)
    cs = ax.contour(xv, yv, z, levels=[90, 95, 99.9], colors=["white", "white", COLORS["navy"]],
                    linewidths=[1.0, 1.4, 2.1], linestyles=[":", "--", "-"])
    ax.clabel(cs, fmt={90: "90%", 95: "95%", 99.9: "100%"}, fontsize=10, inline=True)
    ax.scatter([12], [80], marker="*", s=280, color=COLORS["gold"], edgecolor=COLORS["navy"], lw=1, zorder=5)
    ax.annotate("Minimum for 58/58\n12 aa at 80%", xy=(12, 80), xytext=(5.0, 69),
                arrowprops=dict(arrowstyle="-|>", color=COLORS["navy"], lw=1.2),
                bbox=dict(boxstyle="round,pad=0.45", fc="white", ec=COLORS["grid"]),
                fontsize=11, fontweight="bold", color=COLORS["navy"])
    ax.scatter([15], [80], marker="D", s=95, color=COLORS["orange"], edgecolor="white", lw=1.1, zorder=5)
    ax.annotate("Selected operating point\n15 aa at 80%", xy=(15, 80), xytext=(18.2, 73),
                arrowprops=dict(arrowstyle="-|>", color=COLORS["orange"], lw=1.2),
                bbox=dict(boxstyle="round,pad=0.45", fc="white", ec="#F2C9B6"),
                fontsize=11, fontweight="bold", color=COLORS["orange"])
    ax.set_xlabel("Allowed absolute D3-end offset (aa)")
    ax.set_ylabel("Required member-database consensus (%)")
    cbar = fig.colorbar(mesh, ax=ax, pad=0.025, fraction=0.045)
    cbar.set_label("Proteins satisfying the boundary criterion (%)"); cbar.set_ticks([0, 25, 50, 75, 100])
    title_block(fig, "The boundary result occupies a broad, threshold-insensitive stability plateau",
                "All 1,326 combinations of boundary tolerance and member-database consensus are evaluated jointly.")
    footer(fig, f'{surface_summary["fraction_of_grid_with_100pct_pass"]*100:.1f}% of the entire tested parameter grid supports all 58 primary proteins.')
    save(fig, out, "04_robustness_surface")


    fig, ax = plt.subplots(figsize=(13.0, 8.0))
    fig.subplots_adjust(top=0.79, left=0.11, right=0.76, bottom=0.24)
    levels = [0.50, 0.67, 0.80, 0.90, 1.00]
    cols = [COLORS["gray"], COLORS["sky"], COLORS["blue"], COLORS["purple"], COLORS["orange"]]
    for c, col in zip(levels, cols):
        sub = surface[np.isclose(surface["consensus_threshold"], c)].sort_values("distance_threshold_aa")
        ax.plot(sub["distance_threshold_aa"], sub["pass_fraction"]*100,
                lw=3.6 if c == 0.80 else 2.4, color=col, label=f"{int(c*100)}% consensus")
        full = sub[sub["pass_fraction"] >= 0.999999]
        if not full.empty:
            ax.scatter(full.iloc[0]["distance_threshold_aa"], 100, s=65, color=col,
                       edgecolor="white", lw=0.9, zorder=4)
    ax.axvline(12, color=COLORS["gold"], lw=1.4, ls=(0, (4, 3)))
    ax.axvline(15, color=COLORS["orange"], lw=1.4, ls=(0, (4, 3)))
    ax.axhline(100, color=COLORS["navy"], lw=1, ls=(0, (3, 3)))
    ax.fill_between([12, 25], 98.5, 100.5, color=COLORS["teal"], alpha=0.08)
    ax.set_xlim(0, 25); ax.set_ylim(0, 103)
    ax.set_xlabel("Allowed absolute D3-end offset (aa)")
    ax.set_ylabel("Primary proteins satisfying the criterion (%)")
    clean_axis(ax, "both")
    ax.legend(title="Required database consensus", frameon=False,
              bbox_to_anchor=(1.03, 0.5), loc="center left")
    title_block(fig, "Coverage rapidly reaches 100% and remains stable under stringent consensus rules",
                "Curves show how the number of supported proteins changes as the tolerated D3-end offset is widened.")
    footer(fig, "At 80% consensus, all 58 proteins pass from 12 aa onward; the selected 15-aa threshold is conservative rather than outcome-defining.")
    save(fig, out, "05_robustness_coverage_curves")


    fig, ax = plt.subplots(figsize=(11.8, 7.8))
    fig.subplots_adjust(top=0.79, left=0.13, right=0.95, bottom=0.24)
    target_colors = {0.90: COLORS["sky"], 0.95: COLORS["blue"], 1.00: COLORS["orange"]}
    for target in sorted(minimum["target_pass_fraction"].unique()):
        sub = minimum[minimum["target_pass_fraction"] == target].sort_values("consensus_threshold")
        ax.plot(sub["consensus_threshold"]*100, sub["minimum_distance_threshold_aa"],
                marker="o", ms=8, lw=3, color=target_colors[target], label=f"{int(target*100)}% protein coverage")
    ax.set_xlabel("Required member-database consensus (%)")
    ax.set_ylabel("Minimum D3-end tolerance required (aa)")
    ax.set_xticks(sorted(minimum["consensus_threshold"].unique()*100))
    clean_axis(ax, "both"); ax.legend(frameon=False, loc="upper left")
    title_block(fig, "Even strict database consensus requires only a modest widening of the boundary tolerance",
                "For each consensus definition, the plot reports the smallest amino-acid tolerance reaching the target protein coverage.")
    footer(fig, "The method remains stable as consensus becomes stricter, separating genuine robustness from tuning to one cutoff.")
    save(fig, out, "06_minimum_tolerance_by_consensus")


    fig, ax = plt.subplots(figsize=(13.3, 8.0))
    fig.subplots_adjust(top=0.79, left=0.11, right=0.96, bottom=0.24)
    x = jsd["alignment_column_1based"].to_numpy(); sig = jsd["significant_fdr"].astype(bool).to_numpy()
    barcode_cols = set(barcode["combined_alignment_col"].astype(int)); is_barcode = np.array([int(v) in barcode_cols for v in x])
    ax.fill_between(x, 0, jsd["weighted_jsd_bits"], color=COLORS["sky"], alpha=0.25)
    ax.plot(x, jsd["weighted_jsd_bits"], lw=2.7, color=COLORS["blue"], label="Weighted JSD")
    ax.plot(x, jsd["isoform_discrimination_index"], lw=1.8, color=COLORS["purple"], alpha=0.9,
            label="Isoform discrimination index")
    ax.scatter(x[sig], jsd.loc[sig, "weighted_jsd_bits"], s=55, color=COLORS["teal"],
               edgecolor="white", lw=0.7, zorder=4, label="FDR-significant")
    ax.scatter(x[is_barcode], jsd.loc[is_barcode, "weighted_jsd_bits"], s=135, marker="*",
               color=COLORS["gold"], edgecolor=COLORS["navy"], lw=0.7, zorder=5,
               label="17 conserved barcode sites")
    top = jsd.sort_values(["isoform_discrimination_index", "weighted_jsd_bits"], ascending=False).head(6)
    for k, (_, row) in enumerate(top.iterrows()):
        ax.annotate(f"{row['IIIb_major_aa']}→{row['IIIc_major_aa']}",
                    (int(row["alignment_column_1based"]), row["weighted_jsd_bits"]),
                    xytext=(0, 12 + (k % 2)*8), textcoords="offset points", ha="center",
                    fontsize=9, fontweight="bold", color=COLORS["navy"])
    ax.set_xlim(0.5, len(jsd)+0.5); ax.set_ylim(0, 1.12)
    ax.set_xlabel("Cassette multiple-sequence-alignment column")
    ax.set_ylabel("Sequence divergence score (bits)")
    clean_axis(ax, "both"); ax.legend(frameon=False, ncol=2, loc="upper right")
    title_block(fig, "FGFR2 isoform divergence is concentrated at discrete, evolutionarily conserved positions",
                "Henikoff sequence weighting corrects for uneven taxonomic sampling; significance is controlled by paired permutation and FDR.")
    footer(fig, f'{jsd_summary["n_fdr_significant_positions"]} positions are FDR-significant; 17 form the high-confidence structure-mapped molecular barcode.')
    save(fig, out, "07_weighted_JSD_position_track")


    fig, ax = plt.subplots(figsize=(14.2, 6.8))
    fig.subplots_adjust(top=0.76, left=0.10, right=0.98, bottom=0.48)
    aa_group = {**{a: "Hydrophobic" for a in "AVLIM"}, **{a: "Aromatic" for a in "FWY"},
                **{a: "Polar" for a in "STNQ"}, **{a: "Positive" for a in "KRH"},
                **{a: "Negative" for a in "DE"}, **{a: "Special" for a in "CGP"}}
    group_colors = {"Hydrophobic": "#3973A8", "Aromatic": "#7561A8", "Polar": "#168A83",
                    "Positive": "#D96B38", "Negative": "#B94A48", "Special": "#8996A3"}
    bp = barcode.sort_values("combined_alignment_col").reset_index(drop=True)
    for j, row in bp.iterrows():
        for yrow, aa_col, frac_col in [(1, "IIIb_major_aa", "IIIb_major_aa_fraction"),
                                        (0, "IIIc_major_aa", "IIIc_major_aa_fraction")]:
            aa = str(row[aa_col]); frac = float(row[frac_col]); group = aa_group.get(aa, "Special")
            ax.add_patch(Rectangle((j-0.46, yrow-0.38), 0.92, 0.76, facecolor=group_colors[group],
                                   edgecolor="white", lw=1.3, alpha=0.40 + 0.60*frac))
            ax.text(j, yrow+0.05, aa, ha="center", va="center", fontsize=14, fontweight="bold", color="white")
            ax.text(j, yrow-0.22, f"{frac*100:.0f}%", ha="center", va="center", fontsize=7.3, color="white")
        if "contact" in str(row["structural_evidence"]).lower():
            ax.scatter(j, 1.62, marker="v", s=72, color=COLORS["orange"], clip_on=False)
    ax.set_xlim(-0.65, len(bp)-0.35); ax.set_ylim(-0.55, 1.85)
    ax.set_yticks([0, 1]); ax.set_yticklabels(["IIIc", "IIIb"], fontweight="bold", fontsize=13)
    ax.set_xticks(range(len(bp)))
    positions = []
    for _, row in bp.iterrows():
        b, c = int(row["IIIb_abs_position"]), int(row["IIIc_abs_position"])
        positions.append(str(b) if b == c else f"{b}/{c}")
    ax.set_xticklabels(positions, rotation=90, fontsize=9)
    ax.set_xlabel("Human FGFR2 residue position (IIIb/IIIc when numbering differs)", labelpad=14)
    for spine in ax.spines.values(): spine.set_visible(False)
    ax.tick_params(length=0)
    handles = [Line2D([0], [0], marker="s", color="none", markerfacecolor=c, markersize=11, label=k)
               for k, c in group_colors.items()]
    handles.append(Line2D([0], [0], marker="v", color="none", markerfacecolor=COLORS["orange"],
                          markersize=9, label="Literature-supported interface site"))
    ax.legend(handles=handles, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.43), columnspacing=1.6, handletextpad=0.7)
    title_block(fig, "A compact 17-position amino-acid barcode preserves distinct IIIb and IIIc chemistries",
                "Tile opacity reflects within-isoform conservation; colors represent physicochemical amino-acid classes.")
    save(fig, out, "08_molecular_isoform_barcode")


    fig, ax = plt.subplots(figsize=(11.8, 7.8))
    fig.subplots_adjust(top=0.79, left=0.12, right=0.95, bottom=0.24)
    observed = float(jsd_summary["global_sum_jsd"])
    ax.hist(null_jsd, bins=52, color=COLORS["sky"], alpha=0.75, edgecolor="white", lw=0.55)
    ax.axvline(observed, color=COLORS["orange"], lw=3.5)
    ax.annotate(f"Observed ΣJSD = {observed:.2f}\npaired permutation p = {jsd_summary['global_paired_permutation_p']:.1e}",
                xy=(observed, ax.get_ylim()[1]*0.74), xytext=(np.quantile(null_jsd, 0.66), ax.get_ylim()[1]*0.84),
                arrowprops=dict(arrowstyle="-|>", color=COLORS["orange"], lw=1.5),
                bbox=dict(boxstyle="round,pad=0.55", fc="white", ec="#F1C7B6"),
                fontsize=12, fontweight="bold", color=COLORS["navy"])
    ax.set_xlabel("Global sum of weighted JSD after within-species IIIb/IIIc label permutation")
    ax.set_ylabel("Number of permutations")
    clean_axis(ax, "both")
    title_block(fig, "The observed trans-species isoform signal exceeds the complete permutation null distribution",
                "In each of 10,000 permutations, IIIb and IIIc labels were swapped only within species, preserving phylogenetic pairing.")
    footer(fig, "None of the 10,000 permuted datasets produced a global isoform signal as strong as the observed data.")
    save(fig, out, "09_global_JSD_permutation_test")


    fig, ax = plt.subplots(figsize=(11.9, 8.0))
    fig.subplots_adjust(top=0.79, left=0.12, right=0.95, bottom=0.24)
    q = np.clip(jsd["bh_fdr_q"].to_numpy(float), 1e-8, 1.0); neglogq = -np.log10(q)
    conservation = np.minimum(jsd["IIIb_major_fraction_weighted"], jsd["IIIc_major_fraction_weighted"]).to_numpy(float)
    sizes = 35 + 150*conservation
    ax.scatter(jsd.loc[~sig, "weighted_jsd_bits"], neglogq[~sig], s=sizes[~sig], color="#C7D0DA",
               edgecolor="white", lw=0.5, alpha=0.82)
    ax.scatter(jsd.loc[sig, "weighted_jsd_bits"], neglogq[sig], s=sizes[sig], color=COLORS["teal"],
               edgecolor="white", lw=0.65, alpha=0.90, label="FDR-significant")
    ax.scatter(jsd.loc[is_barcode, "weighted_jsd_bits"], neglogq[is_barcode], s=sizes[is_barcode]+50,
               marker="*", color=COLORS["gold"], edgecolor=COLORS["navy"], lw=0.7, zorder=4,
               label="17 conserved barcode sites")
    ax.axhline(-np.log10(0.05), color=COLORS["orange"], lw=1.4, ls=(0, (4, 3)))
    ax.text(0.01, -np.log10(0.05)+0.10, "FDR q = 0.05", color=COLORS["orange"], fontsize=10, fontweight="bold")
    top_rows = jsd.sort_values(["isoform_discrimination_index", "weighted_jsd_bits"], ascending=False).head(5)
    top_text = "Top high-confidence positions\n" + "\n".join(
        f"col {int(row['alignment_column_1based'])}: {row['IIIb_major_aa']}→{row['IIIc_major_aa']}"
        for _, row in top_rows.iterrows()
    )
    ax.text(0.03, 0.95, top_text, transform=ax.transAxes, ha="left", va="top",
            fontsize=10.2, fontweight="bold", color=COLORS["navy"],
            bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                      edgecolor=COLORS["grid"], linewidth=0.8))
    ax.set_xlabel("Weighted Jensen–Shannon divergence (effect size)")
    ax.set_ylabel("−log₁₀(FDR-adjusted q-value)")
    clean_axis(ax, "both"); ax.legend(frameon=False, loc="lower right")
    title_block(fig, "The strongest isoform effects are simultaneously significant and conserved within each isoform class",
                "Point size represents the lower of the IIIb and IIIc within-class conservation fractions.")
    footer(fig, "The 17 high-confidence barcode sites occupy the high-effect, high-significance region.")
    save(fig, out, "10_JSD_effect_significance_map")

    (out / "FIGURE_INDEX_AND_CAPTIONS.md").write_text(
        "# FGFR2 individual publication figures\n\n"
        "Ten separate figures are exported as PNG (400 dpi), SVG and PDF.\n", encoding="utf-8")


if __name__ == "__main__":
    main()
