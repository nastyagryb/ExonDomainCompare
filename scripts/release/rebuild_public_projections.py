#!/usr/bin/env python3
"""Refresh portable website/download projections without recomputing biology."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from exondomaincompare.shared_gene_analysis.public_paths import (
    rebuild_existing_public_projections,
    sanitize_public_payload,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", action="append", default=[])
    parser.add_argument(
        "--path",
        action="append",
        default=[],
        help="Repository-relative derived/public JSON projection to sanitize in place.",
    )
    args = parser.parse_args()
    if not args.run_id and not args.path:
        parser.error("at least one --run-id or --path is required")
    for run_id in args.run_id:
        run_dir = ROOT / "runs" / run_id
        if not run_dir.is_dir():
            parser.error(f"run not found: {run_id}")
        written = rebuild_existing_public_projections(run_dir)
        print(f"{run_id}: {len(written)} portable projection(s) written")
    for relative in args.path:
        path = (ROOT / relative).resolve()
        try:
            path.relative_to(ROOT.resolve())
        except ValueError:
            parser.error(f"path escapes repository: {relative}")
        if path.suffix.lower() != ".json" or not path.is_file():
            parser.error(f"public JSON projection not found: {relative}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        projected = sanitize_public_payload(payload)
        if projected != payload:
            path.write_text(
                json.dumps(projected, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            print(f"{relative}: portable projection written")
        else:
            print(f"{relative}: already portable")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
