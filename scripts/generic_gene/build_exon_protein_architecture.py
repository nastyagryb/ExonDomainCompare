from __future__ import annotations

import argparse
from typing import Any, Dict

from exondomaincompare.generic_gene.common import GenericContext, load_context, read_tsv, write_tsv

COLUMNS = [
    "species_id", "protein_id", "transcript_id", "exon_number", "exon_id",
    "protein_start_aa", "protein_end_aa", "length_aa", "cds_start", "cds_end",
    "phase", "confidence", "source",
]


def build(ctx: GenericContext) -> Dict[str, Any]:
    ctx.assert_not_freeze()
    rows_in = read_tsv(ctx.core("exon_protein_map.tsv"))
    rows = []
    for r in rows_in:
        try:
            start = int(r.get("protein_start_aa") or 0)
            end = int(r.get("protein_end_aa") or 0)
        except ValueError:
            start = end = 0
        rows.append({
            "species_id": r.get("species_id", ""),
            "protein_id": r.get("protein_id", ""),
            "transcript_id": r.get("transcript_id", ""),
            "exon_number": r.get("exon_number", ""),
            "exon_id": r.get("exon_id", ""),
            "protein_start_aa": start,
            "protein_end_aa": end,
            "length_aa": max(0, end - start + 1) if end and start else "",
            "cds_start": r.get("cds_start", ""),
            "cds_end": r.get("cds_end", ""),
            "phase": r.get("phase", ""),
            "confidence": r.get("confidence", ""),
            "source": r.get("source", ""),
        })
    write_tsv(ctx.out("exon_protein_architecture.tsv"), rows, COLUMNS)
    proteins = sorted({r["protein_id"] for r in rows if r["protein_id"]})
    return {
        "exon_protein_architecture.tsv": len(rows),
        "n_proteins": len(proteins),
        "n_exon_blocks": len(rows),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description='Generic exon/protein architecture table (gene-agnostic).')
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()
    ctx = load_context(args.run_id)
    res = build(ctx)
    print(f"OK exon_protein_architecture  blocks={res['n_exon_blocks']}  proteins={res['n_proteins']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
