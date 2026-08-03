#!/usr/bin/env python3
"""Validate a Core Gene Analysis run against the required milestones.

Prints completed milestones, missing required/optional files, the inferred
status and the suggested next action. Uses the exact same milestone logic as the
webapp backend (src/exondomaincompare/framework/core_run_milestones.py), so a partial or empty
run is classified identically in both places.

Usage:
  python -m exondomaincompare.framework.validate_core_gene_run --run-id <run_id>
  python -m exondomaincompare.framework.validate_core_gene_run --run-id <run_id> --json
"""
from __future__ import annotations

import argparse
import json
import sys

from exondomaincompare.framework.core_run_milestones import evaluate_core_run  # noqa: E402
from exondomaincompare.runs.registry import RegistryError, resolve_run_record  # noqa: E402
from exondomaincompare.config import discover_repository_root, load_config  # noqa: E402

PROJECT_ROOT = discover_repository_root(__file__)
RUNTIME_CONFIG = load_config(repository_root=PROJECT_ROOT)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Validate a Core Gene Analysis run's milestones.")
    ap.add_argument("--run-id", required=True, help="Run id under runs/.")
    ap.add_argument("--json", action="store_true", help="Emit the full report as JSON.")
    args = ap.parse_args(argv)

    try:
        record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
    except RegistryError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    if record is None:
        print(f"ERROR: run not found: {args.run_id}", file=sys.stderr)
        return 2
    run_dir = record.path

    rep = evaluate_core_run(run_dir)
    if args.json:
        print(json.dumps(rep, indent=2))
        return 0

    print(f"Core Gene Analysis run: {rep['run_id']}")
    print(f"  analysis_id     : {rep['analysis_id'] or '(unknown)'}")
    print(f"  gene_symbol     : {rep['gene_symbol'] or '(unknown)'}")
    print(f"  core-only       : {rep['is_core_only']}  (has_event={rep['has_event']})")
    print(f"  inferred status : {rep['inferred_status']}")
    print(f"  next action     : {rep['suggested_next_action']}")
    print()
    print("Milestones:")
    for ms in rep["milestones"]:
        mark = "OK " if ms["complete"] else ("!! " if ms["required"] else ".. ")
        req = "required" if ms["required"] else "optional"
        line = f"  [{mark}] {ms['name']} ({req})"
        print(line)
        if not ms["complete"]:
            if ms["missing_files"]:
                print(f"         missing: {', '.join(ms['missing_files'])}")
            if ms["reason"]:
                print(f"         reason : {ms['reason']}")
    print()
    c = rep["counts"]
    print("Counts:")
    print(f"  gene_models={c['gene_models']}  isoforms={c['protein_isoforms']}  "
          f"primary_proteins={c['primary_proteins']}")
    print(f"  exon_map_rows={c['exon_map_rows']}  synteny_neighbors={c['synteny_neighbors']}")
    print(f"  domain_features={c['domain_features']}  tm_features={c['tm_features']}  "
          f"boundary_rows={c['boundary_rows']}")
    print()
    if rep["missing_required"]:
        print("Missing REQUIRED outputs:")
        for f in rep["missing_required"]:
            print(f"  - {f}")
    else:
        print("All required outputs present for the current phase.")
    if rep["missing_optional"]:
        print("Missing optional outputs:")
        for f in rep["missing_optional"]:
            print(f"  - {f}")
    if rep["inferred_status"] == "cluster_required":
        print(f"\nNext: {rep['cluster_command']}")
    if rep["logs"]:
        print("\nLogs:")
        for lg in rep["logs"]:
            print(f"  - {lg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
