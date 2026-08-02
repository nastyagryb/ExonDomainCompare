#!/usr/bin/env python3
"""Do the returned cluster results still describe the run's proteins?

A run's domain and topology layers are produced off-machine, so the only honest
way to decide whether they may be used is to compare what the cluster scored
against what the run currently asks about. File presence cannot answer that: a
repaired coordinate model can change a protein sequence while every output file
stays exactly where it was, and the run would then present domain calls for a
sequence it no longer analyses.

InterProScan records the MD5 of each scored sequence, and pyTMHMM output is keyed
by the same protein ids, so both can be checked against the primary FASTA
directly. The result is one of:

``fresh``
    every primary protein was scored, and the scored sequence is byte-identical
    to the current one.
``stale``
    at least one protein was scored under a different sequence. The outputs stay
    on disk for diagnostics but must not be used.
``incomplete``
    a primary protein was never scored.
``missing``
    no cluster output exists yet.

Usage::

    python -m shared_gene_analysis.cluster_output_freshness --run-id <run_id>
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from exondomaincompare.config import discover_repository_root, load_config
from exondomaincompare.runs.registry import discover_run

ROOT = discover_repository_root(__file__)
RUNTIME_CONFIG = load_config(repository_root=ROOT)
ROOT = RUNTIME_CONFIG.repository_root
RUNS_ROOT = RUNTIME_CONFIG.runs_root

FRESH = "fresh"
STALE = "stale"
INCOMPLETE = "incomplete"
MISSING = "missing"

# Where the two pipelines keep the primary InterProScan and pyTMHMM payloads.
_INTERPRO_DIRS = (
    "results/14_interproscan/primary/output",
    "results/interproscan/primary/output",
)
_TMHMM_FILES = (
    "results/15_exon_domain_boundary_post_interpro/pytmhmm_primary/output/"
    "pytmhmm_summary_all.tsv",
    "results/15_exon_domain_boundary_post_interpro/pytmhmm_primary/output/"
    "pytmhmm_transmembrane_hits.tsv",
)
_FASTA_CANDIDATES = (
    "results/14_interproscan/primary/input/final_pre_interpro_proteins_primary.faa",
    "results/13_final_pre_interpro_closure/fasta/final_pre_interpro_proteins_primary.faa",
    "results/core_gene_analysis/primary_proteins.faa",
)


def _repo_relative(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def sequence_md5(seq: str) -> str:
    """InterProScan's sequence digest: uppercase residues, no whitespace."""
    return hashlib.md5("".join(seq.split()).upper().encode()).hexdigest()


def read_fasta(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    name: Optional[str] = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = ""
        elif name is not None:
            out[name] += line.strip()
    return out


def primary_fasta(run_dir: Path) -> Optional[Path]:
    for rel in _FASTA_CANDIDATES:
        fp = run_dir / rel
        if fp.is_file():
            return fp
    hits = sorted(run_dir.glob("results/**/*primary*.faa"))
    return hits[0] if hits else None


def _interpro_scored(run_dir: Path) -> Dict[str, str]:
    """protein id -> MD5 of the sequence InterProScan actually scored."""
    scored: Dict[str, str] = {}
    for rel in _INTERPRO_DIRS:
        for fp in sorted((run_dir / rel).glob("*.json")) if (run_dir / rel).is_dir() else []:
            try:
                doc = json.loads(fp.read_text())
            except (OSError, ValueError):
                continue
            for res in doc.get("results") or []:
                md5 = res.get("md5")
                seq = res.get("sequence") or ""
                digest = md5 or (sequence_md5(seq) if seq else "")
                for xref in res.get("xref") or []:
                    pid = xref.get("id")
                    if pid and digest:
                        scored[pid] = digest
    return scored


def _tmhmm_scored(run_dir: Path) -> List[str]:
    """Protein ids present in the returned topology tables."""
    ids: List[str] = []
    for rel in _TMHMM_FILES:
        fp = run_dir / rel
        if not fp.is_file():
            continue
        lines = fp.read_text().splitlines()
        if not lines:
            continue
        header = lines[0].split("\t")
        try:
            col = next(i for i, h in enumerate(header)
                       if h.strip() in ("protein_id", "id", "sequence_id", "query"))
        except StopIteration:
            col = 0
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) > col and parts[col].strip():
                ids.append(parts[col].strip())
    return ids


def evaluate(run_dir: Path) -> Dict[str, Any]:
    """Compare the returned cluster outputs against the current primary FASTA."""
    fasta = primary_fasta(run_dir)
    report: Dict[str, Any] = {
        "run_id": run_dir.name,
        "primary_fasta": _repo_relative(fasta) if fasta else "",
        "n_primary_proteins": 0,
        "interproscan": {"status": MISSING, "n_scored": 0, "mismatched": [],
                         "unscored": [], "extra": []},
        "pytmhmm": {"status": MISSING, "n_scored": 0, "unscored": []},
        "status": MISSING,
    }
    if not fasta:
        report["reason"] = "no primary protein FASTA in this run"
        return report

    proteins = read_fasta(fasta)
    digests = {pid: sequence_md5(seq) for pid, seq in proteins.items()}
    report["n_primary_proteins"] = len(proteins)

    scored = _interpro_scored(run_dir)
    ip = report["interproscan"]
    ip["n_scored"] = len(scored)
    if scored:
        ip["mismatched"] = sorted(pid for pid, d in digests.items()
                                  if pid in scored and scored[pid] != d)
        ip["unscored"] = sorted(pid for pid in digests if pid not in scored)
        ip["extra"] = sorted(pid for pid in scored if pid not in digests)
        ip["status"] = (STALE if ip["mismatched"]
                        else INCOMPLETE if ip["unscored"] else FRESH)

    tm_ids = set(_tmhmm_scored(run_dir))
    tm = report["pytmhmm"]
    tm["n_scored"] = len(tm_ids)
    if tm_ids:
        tm["unscored"] = sorted(pid for pid in digests if pid not in tm_ids)
        tm["status"] = INCOMPLETE if tm["unscored"] else FRESH

    # The run-level verdict is the worst of the two layers: a stale sequence is
    # the strongest signal, and a layer that was never returned still blocks a
    # results_ready claim on any analysis that depends on it.
    order = {FRESH: 0, INCOMPLETE: 1, MISSING: 2, STALE: 3}
    report["status"] = max((ip["status"], tm["status"]), key=lambda s: order[s])
    report["usable"] = report["status"] == FRESH
    return report


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args(argv)
    run_dir = discover_run(RUNTIME_CONFIG, args.run_id) or (
        RUNS_ROOT / args.run_id).resolve()
    if not run_dir.is_dir():
        ap.error(f"no such run: {run_dir}")
    report = evaluate(run_dir)
    print(json.dumps(report, indent=2))
    return 0 if report["status"] in (FRESH, MISSING) else 1


if __name__ == "__main__":
    raise SystemExit(main())
