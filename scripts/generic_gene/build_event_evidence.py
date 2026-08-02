"""Generic exploratory event-evidence (gene-agnostic).

Normalizes the standardized isoform-difference outputs into canonical
``event_region_evidence.tsv`` + ``event_region_candidate_clusters.tsv``. These are
EXPLORATORY isoform-difference candidates, never validated events. No marker /
cassette / IIIb-IIIc logic.
"""
from __future__ import annotations

import argparse
from typing import Any, Dict

from .common import GenericContext, load_context, read_tsv, write_tsv

EVID_COLUMNS = [
    "analysis_id", "gene_symbol", "species_id", "event_candidate_id",
    "evidence_source", "evidence_status", "protein_a", "protein_b",
    "region_start_aa", "region_end_aa", "region_length_aa", "event_type_candidate",
    "exon_aligned", "confidence", "confidence_reason", "notes",
]
CLUSTER_COLUMNS = [
    "candidate_cluster_id", "gene_symbol", "species_id", "representative_start_aa",
    "representative_end_aa", "representative_length_aa", "support_count",
    "proteins_involved", "transcripts_involved", "evidence_sources",
    "exon_aligned_support", "confidence", "confidence_reason", "notes",
]


def build(ctx: GenericContext) -> Dict[str, Any]:
    ctx.assert_not_freeze()
    evid = read_tsv(ctx.core("event_region_evidence.tsv"))
    clusters = read_tsv(ctx.core("event_region_candidate_clusters.tsv"))
    write_tsv(ctx.out("event_region_evidence.tsv"), evid, EVID_COLUMNS)
    write_tsv(ctx.out("event_region_candidate_clusters.tsv"), clusters, CLUSTER_COLUMNS)
    return {
        "event_region_evidence.tsv": len(evid),
        "event_region_candidate_clusters.tsv": len(clusters),
        "n_evidence": len(evid),
        "n_clusters": len(clusters),
        "available": bool(clusters or evid),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()
    ctx = load_context(args.run_id)
    res = build(ctx)
    print(f"OK event_evidence  clusters={res['n_clusters']}  evidence={res['n_evidence']} (exploratory only)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
