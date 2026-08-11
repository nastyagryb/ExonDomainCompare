from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

LOG = logging.getLogger("fgfr2-analysis")


def setup_logging(verbose: bool = False) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def read_tsv(path: str | Path, required: Iterable[str] = ()) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Input file not found: {p}")
    df = pd.read_csv(p, sep="\t")
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{p} is missing required columns: {missing}")
    return df


def protein_key(df: pd.DataFrame) -> pd.Series:
    return df["species"].astype(str) + "|" + df["isoform"].astype(str)


def aggregate_member_calls(
    df: pd.DataFrame,
    offset_col: str = "end_signed_offset",
) -> pd.DataFrame:
    """One observation per protein and InterPro member database.

    Some member databases can emit more than one overlapping signature. Taking
    a within-protein/database median prevents databases with duplicate models
    from receiving extra statistical weight.
    """
    required = ["species", "isoform", "member_database", offset_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns for aggregation: {missing}")
    x = df.dropna(subset=required).copy()
    x[offset_col] = pd.to_numeric(x[offset_col], errors="coerce")
    x = x.dropna(subset=[offset_col])
    agg = (
        x.groupby(["species", "isoform", "member_database"], as_index=False)
        .agg(
            signed_offset=(offset_col, "median"),
            n_signature_calls=(offset_col, "size"),
        )
    )
    agg["abs_offset"] = agg["signed_offset"].abs()
    agg["protein_key"] = protein_key(agg)
    return agg


def bootstrap_species_statistic(
    df: pd.DataFrame,
    statistic,
    n_boot: int = 5000,
    seed: int = 20260804,
) -> tuple[float, float, float]:
    """Cluster bootstrap where species, not individual domain calls, are sampled."""
    species = np.array(sorted(df["species"].dropna().unique()))
    if len(species) == 0:
        return np.nan, np.nan, np.nan
    observed = float(statistic(df))
    rng = np.random.default_rng(seed)
    values = np.empty(n_boot, dtype=float)
    groups = {s: df[df["species"] == s] for s in species}
    for b in range(n_boot):
        sampled = rng.choice(species, size=len(species), replace=True)
        pieces = []
        for j, s in enumerate(sampled):
            part = groups[s].copy()
            part["_bootstrap_species"] = f"{s}__{j}"
            pieces.append(part)
        boot_df = pd.concat(pieces, ignore_index=True)
        values[b] = statistic(boot_df)
    lo, hi = np.nanquantile(values, [0.025, 0.975])
    return observed, float(lo), float(hi)


def write_json(data: dict, path: str | Path) -> None:
    Path(path).write_text(json.dumps(data, indent=2, default=str) + "\n")
