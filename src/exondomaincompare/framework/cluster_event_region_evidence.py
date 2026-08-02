#!/usr/bin/env python3
"""Cluster / rank the exploratory event-region evidence into candidate regions.

Raw pairwise isoform comparisons are noisy: the same biological region shows up
once per isoform pair. This groups overlapping / nearby regions (within the same
gene + species + similar candidate type) into a small number of candidate
CLUSTERS, and records how many isoform pairs support each one.

The output stays EXPLORATORY: nothing here validates an event region or turns on
any event-specific analysis. It only summarises evidence so a user sees a few
meaningful candidate regions instead of dozens of raw rows.

Input (from a core run's results/core_gene_analysis/):
  * event_region_evidence.tsv    (produced by build_event_region_evidence.py)

Output (same directory):
  * event_region_candidate_clusters.tsv   (always written; header-only if empty)

Usage:
  python scripts/framework/cluster_event_region_evidence.py --run-id <run_id>
  python scripts/framework/cluster_event_region_evidence.py --core-dir <dir> [--gap 5]
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any, Dict, List, Optional

from exondomaincompare.config import discover_repository_root, load_config
from exondomaincompare.runs.registry import RegistryError, resolve_run_record

PROJECT_ROOT = discover_repository_root(__file__)
RUNTIME_CONFIG = load_config(repository_root=PROJECT_ROOT)

EVIDENCE_FILE = "event_region_evidence.tsv"
CLUSTERS_FILE = "event_region_candidate_clusters.tsv"

CLUSTER_COLUMNS = [
    "candidate_cluster_id", "gene_symbol", "species_id",
    "representative_start_aa", "representative_end_aa", "representative_length_aa",
    "support_count", "proteins_involved", "transcripts_involved",
    "evidence_sources", "exon_aligned_support",
    "confidence", "confidence_reason", "notes",
]

# Regions whose type family differs are never merged (an insertion is not a
# substitution). Exon-aligned variants map to the same family as their base type.
_TYPE_FAMILY = {
    "insertion": "indel", "deletion": "indel",
    "substitution_block": "substitution",
}

# Conservative confidence ranking; a cluster is at most as confident as its
# strongest supporting row, but we never exceed "medium" for exploratory data.
_CONF_RANK = {"": 0, "low": 1, "medium": 2, "high": 3}
_MAX_EXPLORATORY_CONF = "medium"


def read_tsv(p: Path) -> List[Dict[str, str]]:
    if not Path(p).is_file():
        return []
    with open(p, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, columns: List[str], rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in columns})


def _int(v: Any) -> Optional[int]:
    try:
        if v is None or v == "":
            return None
        return int(float(v))
    except (TypeError, ValueError):
        return None


def _family(event_type: str) -> str:
    return _TYPE_FAMILY.get((event_type or "").strip(), (event_type or "").strip() or "other")


def _overlap_or_near(a_start: int, a_end: int, b_start: int, b_end: int, gap: int) -> bool:
    """True if [a] and [b] overlap or are within `gap` aa of each other."""
    return (a_start - gap) <= b_end and (b_start - gap) <= a_end


def cluster_rows(evidence: List[Dict[str, str]], gap: int) -> List[Dict[str, Any]]:
    """Greedy 1-D interval clustering per (gene, species, type family)."""
    # bucket by (gene, species, family)
    buckets: Dict[Any, List[Dict[str, Any]]] = {}
    for r in evidence:
        s, e = _int(r.get("region_start_aa")), _int(r.get("region_end_aa"))
        if s is None or e is None:
            continue
        key = (r.get("gene_symbol", ""), r.get("species_id", ""),
               _family(r.get("event_type_candidate", "")))
        buckets.setdefault(key, []).append({**r, "_s": s, "_e": e})

    clusters: List[Dict[str, Any]] = []
    _cluster_n = 0
    for (gene, species, family), rows in buckets.items():
        rows.sort(key=lambda x: (x["_s"], x["_e"]))
        # greedy sweep merging nearby/overlapping intervals
        groups: List[List[Dict[str, Any]]] = []
        for row in rows:
            placed = False
            for grp in groups:
                gs = min(x["_s"] for x in grp)
                ge = max(x["_e"] for x in grp)
                if _overlap_or_near(gs, ge, row["_s"], row["_e"], gap):
                    grp.append(row)
                    placed = True
                    break
            if not placed:
                groups.append([row])

        for grp in groups:
            _cluster_n += 1
            starts = [x["_s"] for x in grp]
            ends = [x["_e"] for x in grp]
            # Representative region: the span reported by most isoform pairs. Ties are
            # broken by the longer span, then by the earlier start. Without that
            # tie-break a group that mixes a short block with a long one a few residues
            # away — every span observed once — would be represented by whichever row
            # happened to sort first, and the longer block would silently vanish from
            # every view built on the clusters while staying in the raw evidence.
            span_counts: Dict[Any, int] = {}
            for x in grp:
                span_counts[(x["_s"], x["_e"])] = span_counts.get((x["_s"], x["_e"]), 0) + 1
            rep_start, rep_end = max(
                span_counts.items(),
                key=lambda kv: (kv[1], kv[0][1] - kv[0][0], -kv[0][0]))[0]
            proteins = sorted({p for x in grp for p in (x.get("protein_a", ""), x.get("protein_b", "")) if p})
            transcripts = sorted({t for x in grp for t in (x.get("transcript_a", ""), x.get("transcript_b", "")) if t})
            sources = sorted({x.get("evidence_source", "") for x in grp if x.get("evidence_source")})
            exon_aligned_support = sum(1 for x in grp if str(x.get("exon_aligned", "")).lower() == "yes")

            best_conf = ""
            for x in grp:
                if _CONF_RANK.get(x.get("confidence", ""), 0) > _CONF_RANK.get(best_conf, 0):
                    best_conf = x.get("confidence", "")
            if _CONF_RANK.get(best_conf, 0) > _CONF_RANK[_MAX_EXPLORATORY_CONF]:
                best_conf = _MAX_EXPLORATORY_CONF

            support = len(grp)
            reason_bits = [f"Supported by {support} isoform pair(s)"]
            if exon_aligned_support:
                reason_bits.append(f"{exon_aligned_support} exon-aligned")
            reason_bits.append("exploratory only, not a validated event")
            clusters.append({
                "candidate_cluster_id": f"{gene}:{species}:{family}:{rep_start}-{rep_end}",
                "gene_symbol": gene,
                "species_id": species,
                "representative_start_aa": rep_start,
                "representative_end_aa": rep_end,
                "representative_length_aa": rep_end - rep_start + 1,
                "support_count": support,
                "proteins_involved": ";".join(proteins),
                "transcripts_involved": ";".join(transcripts),
                "evidence_sources": ";".join(sources),
                "exon_aligned_support": exon_aligned_support,
                "confidence": best_conf or "low",
                "confidence_reason": "; ".join(reason_bits) + ".",
                "notes": "Exploratory candidate region cluster; not a validated event region.",
                "_span": (min(starts), max(ends)),
            })

    # rank: exon-aligned support first, then support count, then longer region
    clusters.sort(key=lambda c: (c["exon_aligned_support"], c["support_count"],
                                 c["representative_length_aa"]), reverse=True)
    for c in clusters:
        c.pop("_span", None)
    return clusters


def build_clusters(core_dir: Path, gap: int) -> Dict[str, Any]:
    evidence = read_tsv(core_dir / EVIDENCE_FILE)
    clusters = cluster_rows(evidence, gap)
    out = core_dir / CLUSTERS_FILE
    write_tsv(out, CLUSTER_COLUMNS, clusters)
    return {"out": out, "n_clusters": len(clusters), "n_evidence": len(evidence)}


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Cluster exploratory event-region evidence into candidate regions.")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--run-id", help="Run id under runs/.")
    g.add_argument("--core-dir", help="Explicit results/core_gene_analysis directory.")
    ap.add_argument("--gap", type=int, default=5,
                    help="Max aa gap between regions to still merge them (default 5).")
    args = ap.parse_args(argv)

    if args.run_id:
        try:
            record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
        except RegistryError as exc:
            raise SystemExit(str(exc)) from exc
        if record is None:
            raise SystemExit(f"Run not found: {args.run_id}")
        if record.read_only:
            raise SystemExit(
                "Run is registered read-only; copy it before writing derived evidence.")
        run_dir = record.path
        core_dir = run_dir / "results" / "core_gene_analysis"
    else:
        core_dir = Path(args.core_dir)
        if not core_dir.is_absolute():
            core_dir = PROJECT_ROOT / core_dir

    core_dir.mkdir(parents=True, exist_ok=True)
    res = build_clusters(core_dir, args.gap)
    print(f"OK  candidate clusters: {res['n_clusters']} cluster(s) "
          f"from {res['n_evidence']} evidence row(s).")
    print(f"    wrote: {res['out']}")
    print("    Clusters are exploratory summaries; no validated event region is asserted.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
