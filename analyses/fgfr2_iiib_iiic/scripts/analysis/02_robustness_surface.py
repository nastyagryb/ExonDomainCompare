
from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import aggregate_member_calls, ensure_dir, read_tsv, setup_logging, write_json


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute a two-dimensional boundary robustness surface.")
    p.add_argument("--calls", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--distance-min", type=float, default=0)
    p.add_argument("--distance-max", type=float, default=25)
    p.add_argument("--distance-step", type=float, default=1)
    p.add_argument("--consensus-min", type=float, default=0.50)
    p.add_argument("--consensus-max", type=float, default=1.00)
    p.add_argument("--consensus-step", type=float, default=0.01)
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    setup_logging(args.verbose)
    out = ensure_dir(args.outdir)
    calls = read_tsv(
        args.calls,
        required=["species", "isoform", "member_database", "end_signed_offset"],
    )
    agg = aggregate_member_calls(calls)
    distances = np.round(
        np.arange(args.distance_min, args.distance_max + args.distance_step / 2, args.distance_step), 6
    )
    consensuses = np.round(
        np.arange(args.consensus_min, args.consensus_max + args.consensus_step / 2, args.consensus_step), 6
    )
    proteins = agg[["species", "isoform"]].drop_duplicates()
    n_proteins = len(proteins)

    rows = []
    support_by_distance: dict[float, pd.DataFrame] = {}
    for distance in distances:
        per = (
            agg.assign(supported=agg["abs_offset"] <= distance)
            .groupby(["species", "isoform"], as_index=False)
            .agg(
                support_fraction=("supported", "mean"),
                n_databases=("member_database", "nunique"),
            )
        )
        support_by_distance[float(distance)] = per
        for consensus in consensuses:
            n_pass = int((per["support_fraction"] >= consensus).sum())
            rows.append(
                {
                    "distance_threshold_aa": distance,
                    "consensus_threshold": consensus,
                    "n_pass": n_pass,
                    "n_proteins": n_proteins,
                    "pass_fraction": n_pass / n_proteins,
                }
            )
    surface = pd.DataFrame(rows)
    surface.to_csv(out / "robustness_surface_long.tsv", sep="\t", index=False)
    pivot = surface.pivot(
        index="consensus_threshold", columns="distance_threshold_aa", values="pass_fraction"
    )
    pivot.to_csv(out / "robustness_surface_matrix.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(10, 6.5))
    im = ax.imshow(
        100 * pivot.values,
        origin="lower",
        aspect="auto",
        extent=[distances.min(), distances.max(), consensuses.min(), consensuses.max()],
        vmin=0,
        vmax=100,
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Proteins passing (%)")
    ax.set_xlabel("Allowed absolute D3-end offset (aa)")
    ax.set_ylabel("Required member-database consensus")
    ax.set_title("Annotation-aware boundary robustness surface")
    fig.tight_layout()
    fig.savefig(out / "figure_robustness_surface.png", dpi=300)
    fig.savefig(out / "figure_robustness_surface.svg")
    plt.close(fig)

    selected_consensus = [0.50, 0.67, 0.80, 0.90, 1.00]
    threshold_rows = []
    for c in selected_consensus:
        sub = surface[np.isclose(surface["consensus_threshold"], c)]
        for target in [0.90, 0.95, 1.00]:
            eligible = sub[sub["pass_fraction"] >= target]
            threshold_rows.append(
                {
                    "consensus_threshold": c,
                    "target_pass_fraction": target,
                    "minimum_distance_threshold_aa": (
                        eligible["distance_threshold_aa"].min() if not eligible.empty else np.nan
                    ),
                }
            )
    thresholds = pd.DataFrame(threshold_rows)
    thresholds.to_csv(out / "minimum_thresholds_for_target_coverage.tsv", sep="\t", index=False)

    plateau = surface[(surface["pass_fraction"] == 1.0)]
    write_json(
        {
            "n_proteins": n_proteins,
            "n_grid_points": len(surface),
            "fraction_of_grid_with_100pct_pass": float(len(plateau) / len(surface)),
            "minimum_distance_for_100pct_at_80pct_consensus": float(
                thresholds.loc[
                    np.isclose(thresholds["consensus_threshold"], 0.80)
                    & np.isclose(thresholds["target_pass_fraction"], 1.00),
                    "minimum_distance_threshold_aa",
                ].iloc[0]
            ),
        },
        out / "robustness_surface_summary.json",
    )


if __name__ == "__main__":
    main()
