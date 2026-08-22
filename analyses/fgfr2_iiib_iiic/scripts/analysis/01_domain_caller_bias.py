
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from common import (
    aggregate_member_calls,
    ensure_dir,
    read_tsv,
    setup_logging,
    write_json,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Estimate systematic InterPro member-database boundary offsets and perform leave-one-database-out robustness analysis."
    )
    p.add_argument("--calls", required=True, help="interpro_ensemble_coordinate_support_calls.tsv")
    p.add_argument("--outdir", required=True)
    p.add_argument("--distance-threshold", type=float, default=15.0)
    p.add_argument("--consensus", type=float, default=0.80)
    p.add_argument("--bootstrap", type=int, default=5000)
    p.add_argument("--seed", type=int, default=20260804)
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def partial_r2_database(agg: pd.DataFrame) -> tuple[float, float, float]:



    reduced = smf.ols("signed_offset ~ C(species) + C(isoform)", data=agg).fit()
    full = smf.ols(
        "signed_offset ~ C(species) + C(isoform) + C(member_database)", data=agg
    ).fit()
    sse_reduced = float(np.sum(reduced.resid**2))
    sse_full = float(np.sum(full.resid**2))
    partial = (sse_reduced - sse_full) / sse_reduced if sse_reduced else np.nan
    return float(partial), float(reduced.rsquared), float(full.rsquared)


def main() -> None:
    args = parse_args()
    setup_logging(args.verbose)
    out = ensure_dir(args.outdir)
    calls = read_tsv(
        args.calls,
        required=["species", "isoform", "member_database", "end_signed_offset"],
    )
    agg = aggregate_member_calls(calls)
    agg.to_csv(out / "member_database_offsets_per_protein.tsv", sep="\t", index=False)

    rows = []
    for db, sub in agg.groupby("member_database", sort=True):


        species_means = sub.groupby("species")["signed_offset"].mean().to_numpy(float)
        mean = float(species_means.mean())
        rng = np.random.default_rng(args.seed)
        sampled = rng.choice(species_means, size=(args.bootstrap, len(species_means)), replace=True)
        lo, hi = np.quantile(sampled.mean(axis=1), [0.025, 0.975])
        rows.append(
            {
                "member_database": db,
                "n_proteins": sub["protein_key"].nunique(),
                "n_species": sub["species"].nunique(),
                "mean_signed_offset_aa": mean,
                "bootstrap_ci_low_aa": lo,
                "bootstrap_ci_high_aa": hi,
                "median_signed_offset_aa": sub["signed_offset"].median(),
                "median_absolute_offset_aa": sub["abs_offset"].median(),
                "fraction_within_threshold": (sub["abs_offset"] <= args.distance_threshold).mean(),
            }
        )
    bias = pd.DataFrame(rows).sort_values("mean_signed_offset_aa")
    bias.to_csv(out / "member_database_bias_estimates.tsv", sep="\t", index=False)

    partial, reduced_r2, full_r2 = partial_r2_database(agg)

    databases = sorted(agg["member_database"].unique())
    loo_rows = []
    scenarios: list[tuple[str, str | None]] = [("none", None)] + [(db, db) for db in databases]
    for label, excluded in scenarios:
        x = agg if excluded is None else agg[agg["member_database"] != excluded]
        per_protein = (
            x.assign(supported=x["abs_offset"] <= args.distance_threshold)
            .groupby(["species", "isoform"], as_index=False)
            .agg(
                n_databases=("member_database", "nunique"),
                support_fraction=("supported", "mean"),
            )
        )
        per_protein["pass"] = per_protein["support_fraction"] >= args.consensus
        n = len(per_protein)
        k = int(per_protein["pass"].sum())
        per_species = per_protein.groupby("species", as_index=False).agg(
            n_models=("isoform", "size"),
            all_observed_models_pass=("pass", "all"),
        )
        n_species = len(per_species)
        k_species = int(per_species["all_observed_models_pass"].sum())
        loo_rows.append(
            {
                "excluded_database": label,
                "n_proteins": n,
                "n_pass": k,
                "pass_fraction": k / n if n else np.nan,
                "n_species": n_species,
                "n_species_all_observed_models_pass": k_species,
                "species_complete_pass_fraction": k_species / n_species if n_species else np.nan,
                "distance_threshold_aa": args.distance_threshold,
                "consensus_threshold": args.consensus,
                "min_remaining_databases": per_protein["n_databases"].min() if n else np.nan,
            }
        )
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(out / "leave_one_database_out.tsv", sep="\t", index=False)



    model = smf.ols(
        "signed_offset ~ C(member_database) + C(isoform)", data=agg
    ).fit(cov_type="cluster", cov_kwds={"groups": agg["species"]})
    regression = pd.DataFrame(
        {
            "term": model.params.index,
            "estimate": model.params.values,
            "cluster_robust_se": model.bse.values,
            "p_value": model.pvalues.values,
            "ci_low": model.conf_int()[0].values,
            "ci_high": model.conf_int()[1].values,
        }
    )
    regression.to_csv(out / "caller_bias_cluster_robust_regression.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(8, max(4, 0.55 * len(bias))))
    y = np.arange(len(bias))
    xerr = np.vstack(
        [
            bias["mean_signed_offset_aa"] - bias["bootstrap_ci_low_aa"],
            bias["bootstrap_ci_high_aa"] - bias["mean_signed_offset_aa"],
        ]
    )
    ax.errorbar(bias["mean_signed_offset_aa"], y, xerr=xerr, fmt="o", capsize=3)
    ax.axvline(0, linewidth=1)
    ax.set_yticks(y, bias["member_database"])
    ax.set_xlabel("Signed D3-end offset relative to cassette end (aa)")
    ax.set_title("Systematic InterPro member-database boundary offsets")
    fig.tight_layout()
    fig.savefig(out / "figure_domain_caller_bias.png", dpi=300)
    fig.savefig(out / "figure_domain_caller_bias.svg")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.bar(loo["excluded_database"], 100 * loo["pass_fraction"])
    ax.set_ylim(0, 105)
    ax.set_ylabel("Proteins passing consensus (%)")
    ax.set_xlabel("Excluded member database")
    ax.set_title(
        f"Leave-one-database-out robustness (|offset|≤{args.distance_threshold:g} aa, consensus≥{args.consensus:.0%})"
    )
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(out / "figure_leave_one_database_out.png", dpi=300)
    fig.savefig(out / "figure_leave_one_database_out.svg")
    plt.close(fig)

    write_json(
        {
            "n_aggregated_protein_database_observations": len(agg),
            "n_proteins": int(agg["protein_key"].nunique()),
            "n_species": int(agg["species"].nunique()),
            "n_member_databases": int(agg["member_database"].nunique()),
            "database_partial_r2_controlling_species_and_isoform": partial,
            "reduced_model_r2": reduced_r2,
            "full_model_r2": full_r2,
            "all_leave_one_out_scenarios_pass_all_proteins": bool((loo["pass_fraction"] == 1).all()),
            "all_leave_one_out_scenarios_pass_all_species": bool((loo["species_complete_pass_fraction"] == 1).all()),
            "success_fractions_are_descriptive": True,
            "protein_level_binomial_interval_reported": False,
        },
        out / "domain_caller_bias_summary.json",
    )


if __name__ == "__main__":
    main()
