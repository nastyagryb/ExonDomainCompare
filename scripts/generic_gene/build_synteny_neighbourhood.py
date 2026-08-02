"""Generic synteny neighbourhood table (gene-agnostic).

Reads standardized ``synteny_neighbors.tsv`` and writes canonical
``synteny_neighbourhood.tsv``. Synteny supports locus / orthology of the target
gene; it does not by itself assign isoform identity. No FGFR2 anchor wording.
"""
from __future__ import annotations

import argparse
from typing import Any, Dict

from .common import GenericContext, load_context, read_tsv, write_tsv

COLUMNS = [
    "species_id", "gene_symbol", "neighbor_symbol", "side", "order",
    "orientation", "classification", "status", "source",
]


def build(ctx: GenericContext) -> Dict[str, Any]:
    ctx.assert_not_freeze()
    rows_in = read_tsv(ctx.core("synteny_neighbors.tsv"))
    rows = [{
        "species_id": r.get("species_id", ""),
        "gene_symbol": r.get("gene_symbol", ctx.gene_symbol),
        "neighbor_symbol": r.get("neighbor_symbol", ""),
        "side": r.get("side", ""),
        "order": r.get("order", ""),
        "orientation": r.get("orientation", ""),
        "classification": r.get("classification", ""),
        "status": r.get("status", ""),
        "source": r.get("source", ""),
    } for r in rows_in]
    write_tsv(ctx.out("synteny_neighbourhood.tsv"), rows, COLUMNS)
    resolved = [r for r in rows if r["status"] == "resolved"]
    return {
        "synteny_neighbourhood.tsv": len(rows),
        "n_neighbors": len(rows),
        "n_resolved_neighbors": len(resolved),
        "available": bool(rows),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()
    ctx = load_context(args.run_id)
    res = build(ctx)
    print(f"OK synteny_neighbourhood  neighbours={res['n_neighbors']}  resolved={res['n_resolved_neighbors']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
