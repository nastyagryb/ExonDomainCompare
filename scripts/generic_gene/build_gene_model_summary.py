"""Generic gene/transcript/protein model summary (gene-agnostic).

Reads the standardized ``gene_model_index.tsv`` + ``protein_isoform_index.tsv``
and writes canonical ``gene_model_summary.tsv`` and ``protein_isoform_summary.tsv``.
No FGFR2 / event assumptions.
"""
from __future__ import annotations

import argparse
from typing import Any, Dict

from .common import GenericContext, load_context, read_tsv, write_tsv

GENE_MODEL_COLUMNS = [
    "analysis_id", "gene_symbol", "species_id", "gene_id", "transcript_id",
    "protein_id", "source", "protein_length", "model_status", "notes",
]
ISOFORM_COLUMNS = [
    "species_id", "protein_id", "transcript_id", "protein_length",
    "primary_status", "source_kind", "notes",
]


def _source_kind(protein_id: str, transcript_id: str) -> str:
    pid = (protein_id or "").upper()
    tid = (transcript_id or "").upper()
    if pid.startswith("NP_") or tid.startswith("NM_"):
        return "refseq_curated"
    if pid.startswith("XP_") or tid.startswith("XM_"):
        return "refseq_predicted"
    if pid.startswith("ENSP") or tid.startswith("ENST"):
        return "ensembl"
    return "other"


def build(ctx: GenericContext) -> Dict[str, Any]:
    ctx.assert_not_freeze()
    gmi = read_tsv(ctx.core("gene_model_index.tsv"))
    iso = read_tsv(ctx.core("protein_isoform_index.tsv"))

    write_tsv(ctx.out("gene_model_summary.tsv"), gmi, GENE_MODEL_COLUMNS)

    iso_rows = []
    for r in iso:
        iso_rows.append({
            "species_id": r.get("species_id", ""),
            "protein_id": r.get("protein_id", ""),
            "transcript_id": r.get("transcript_id", ""),
            "protein_length": r.get("protein_length", ""),
            "primary_status": r.get("primary_status", ""),
            "source_kind": _source_kind(r.get("protein_id", ""), r.get("transcript_id", "")),
            "notes": r.get("notes", ""),
        })
    write_tsv(ctx.out("protein_isoform_summary.tsv"), iso_rows, ISOFORM_COLUMNS)

    species = sorted({r.get("species_id", "") for r in gmi if r.get("species_id")})
    return {
        "gene_model_summary.tsv": len(gmi),
        "protein_isoform_summary.tsv": len(iso_rows),
        "n_species": len(species),
        "n_gene_models": len(gmi),
        "n_protein_isoforms": len(iso_rows),
        "species": species,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()
    ctx = load_context(args.run_id)
    res = build(ctx)
    print(f"OK gene_model_summary  gene={ctx.gene_symbol}  models={res['n_gene_models']}  isoforms={res['n_protein_isoforms']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
