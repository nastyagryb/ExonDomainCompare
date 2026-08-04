#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import difflib
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from exondomaincompare.config import discover_repository_root, load_config

PROJECT_ROOT = discover_repository_root(__file__)
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))
from exondomaincompare.runs.registry import RegistryError, resolve_run_record
RUNTIME_CONFIG = load_config(repository_root=PROJECT_ROOT)

COLUMNS = [
    "species_id", "gene_symbol", "transcript_a", "transcript_b",
    "protein_a", "protein_b", "coordinate_reference_protein",
    "candidate_start_aa", "candidate_end_aa",
    "candidate_length_aa", "candidate_type", "evidence", "confidence", "notes",
]


def read_tsv(p: Path) -> List[Dict[str, str]]:
    if not Path(p).is_file():
        return []
    with open(p, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_fasta(p: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    if not Path(p).is_file():
        return seqs
    cur, buf = None, []
    with open(p, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0].strip()
                buf = []
            else:
                buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


def _diff_regions(a: str, b: str) -> List[Tuple[str, int, int]]:
    out: List[Tuple[str, int, int]] = []
    sm = difflib.SequenceMatcher(a=a, b=b, autojunk=False)
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "insert":
            out.append(("insertion", j1 + 1, j2))
        elif tag == "delete":
            out.append(("deletion", i1 + 1, i2))
        elif tag == "replace":
            out.append(("substitution_block", i1 + 1, i2))
    return out


# Length at which a block is long enough to be a plausible coding cassette on its own.
# Below it a block still needs a second, independent signal to rise above "low".
_SUBSTANTIAL_AA = 6


def _confidence(length_aa: int, n_comparisons: int, exon_aligned: bool) -> str:
    score = (2 if exon_aligned else 0)
    score += 1 if n_comparisons >= 2 else 0
    score += 1 if length_aa >= _SUBSTANTIAL_AA else 0
    return "high" if score >= 3 else "medium" if score >= 2 else "low"


def scan(run_dir: Path) -> List[Dict[str, Any]]:
    core = run_dir / "results" / "core_gene_analysis"
    rc = json.loads((run_dir / "run_config.json").read_text(encoding="utf-8")) \
        if (run_dir / "run_config.json").is_file() else {}
    gene_symbol = rc.get("gene_symbol", "")

    iso = read_tsv(core / "protein_isoform_index.tsv")
    exon_map = read_tsv(core / "exon_protein_map.tsv")

    # collect sequences from any FASTA present in the core dir
    seqs: Dict[str, str] = {}
    for faa in core.glob("*.faa"):
        seqs.update(load_fasta(faa))

    # exon-set per protein for exon-aligned candidates
    exons_by_prot: Dict[str, List[Tuple[int, int]]] = {}
    for e in exon_map:
        pid = e.get("protein_id", "")
        try:
            s, t = int(e.get("protein_start_aa")), int(e.get("protein_end_aa"))
        except (TypeError, ValueError):
            continue
        exons_by_prot.setdefault(pid, []).append((s, t))

    # group isoforms by species
    by_species: Dict[str, List[Dict[str, str]]] = {}
    for r in iso:
        by_species.setdefault(r.get("species_id", ""), []).append(r)

    rows: List[Dict[str, Any]] = []
    for species, isos in by_species.items():
        # pairwise compare isoforms we actually have sequences for
        have = [r for r in isos if seqs.get(r.get("protein_id", ""))]
        for i in range(len(have)):
            for j in range(i + 1, len(have)):
                ra, rb = have[i], have[j]
                pa, pb = ra["protein_id"], rb["protein_id"]
                sa, sb = seqs[pa], seqs[pb]
                if sa == sb:
                    continue
                for ctype, s, e in _diff_regions(sa, sb):
                    length = e - s + 1
                    # Insertions use protein B coordinates; other blocks use protein A.
                    reference = pb if ctype == "insertion" else pa
                    # exon-aligned if the region matches an exon block in either protein
                    exon_aligned = any(abs(s - es) <= 2 and abs(e - et) <= 2
                                       for es, et in exons_by_prot.get(pa, []) + exons_by_prot.get(pb, []))
                    rows.append({
                        "species_id": species, "gene_symbol": gene_symbol,
                        "transcript_a": ra.get("transcript_id", ""),
                        "transcript_b": rb.get("transcript_id", ""),
                        "protein_a": pa, "protein_b": pb,
                        "coordinate_reference_protein": reference,
                        "candidate_start_aa": s, "candidate_end_aa": e,
                        "candidate_length_aa": length,
                        "candidate_type": ("exon_aligned_" + ctype) if exon_aligned else ctype,
                        "evidence": "isoform_sequence_diff"
                                    + ("+exon_map" if exon_aligned else ""),
                        "confidence": "",  # filled in below, once support is known
                        "notes": "Exploratory candidate; not a validated event region.",
                        "_exon_aligned": exon_aligned,
                    })

    # How many independent isoform comparisons report the very same block. This is a
    # scoring input, so it can only be known after every pair has been compared.
    support: Dict[Tuple[str, str, int, int], int] = {}
    for r in rows:
        support[(r["species_id"], r["candidate_type"],
                 r["candidate_start_aa"], r["candidate_end_aa"])] = 0
    for r in rows:
        support[(r["species_id"], r["candidate_type"],
                 r["candidate_start_aa"], r["candidate_end_aa"])] += 1
    for r in rows:
        key = (r["species_id"], r["candidate_type"],
               r["candidate_start_aa"], r["candidate_end_aa"])
        r["confidence"] = _confidence(int(r["candidate_length_aa"]), support[key],
                                      bool(r.pop("_exon_aligned")))
    return rows


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Exploratory isoform-specific region scan (optional, non-blocking).")
    ap.add_argument("--run-id", required=True, help="Run id under runs/.")
    args = ap.parse_args(argv)

    try:
        record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
    except RegistryError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    if record is None:
        print(f"ERROR: run not found: {args.run_id}", file=sys.stderr)
        return 2
    if record.read_only:
        print("ERROR: run is registered read-only; copy it before retrying.",
              file=sys.stderr)
        return 2
    run_dir = record.path

    rows = scan(run_dir)
    out = run_dir / "results" / "core_gene_analysis" / "event_candidate_regions.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    if rows:
        print(f"Found {len(rows)} exploratory isoform-specific candidate region(s).")
        print("These are candidates only — NOT validated event regions.")
    else:
        print("No protein-isoform difference block was detected "
              "(needs >=2 isoforms with sequences that are not identical).")
        print("Core gene-level analysis is unaffected.")
    print(f"Wrote: {out.relative_to(PROJECT_ROOT)}")

    # Always (re)build the summarised evidence + candidate-cluster layers so the
    # webapp can prioritise a small number of clustered candidate regions over
    # the raw pairwise rows. Both are non-blocking and stay exploratory.
    try:
        from exondomaincompare.framework.build_event_region_evidence import build_evidence  # noqa: E402
        from exondomaincompare.framework.cluster_event_region_evidence import build_clusters  # noqa: E402
        rc = json.loads((run_dir / "run_config.json").read_text(encoding="utf-8")) \
            if (run_dir / "run_config.json").is_file() else {}
        ev = build_evidence(out.parent, rc.get("analysis_id", ""), rc.get("gene_symbol", ""))
        cl = build_clusters(out.parent, gap=5)
        print(f"Wrote: {ev['out'].relative_to(PROJECT_ROOT)} ({ev['n_evidence']} evidence row(s))")
        print(f"Wrote: {cl['out'].relative_to(PROJECT_ROOT)} ({cl['n_clusters']} candidate cluster(s))")
    except Exception as exc:  # noqa: BLE001 - evidence layer is best-effort
        print(f"WARN could not build evidence/cluster layer (non-blocking): {exc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
