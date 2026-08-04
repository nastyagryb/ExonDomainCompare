from __future__ import annotations

import argparse
from typing import Any, Dict

from exondomaincompare.generic_gene.common import GenericContext, load_context, read_json, read_tsv

try:
    from exondomaincompare.framework import primary_selection as _ps
except Exception:
    import sys
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
    from exondomaincompare.framework import primary_selection as _ps


def build(ctx: GenericContext) -> Dict[str, Any]:
    ctx.assert_not_freeze()
    iso = read_tsv(ctx.core("protein_isoform_index.tsv"))
    collection = read_json(ctx.core("core_model_collection_report.json"), {}) or {}
    report = _ps.build_primary_selection(iso, collection_report=collection)
    tsv_path = ctx.out("primary_selection_evidence.tsv")
    json_path = ctx.out("primary_selection_report.json")
    _ps.write_selection_evidence(report, tsv_path, json_path)
    return {
        "primary_selection_evidence.tsv": len(report.get("proteins", [])),
        "primary_protein_id": report.get("primary_protein_id", ""),
        "selection_rule": report.get("selection_rule", ""),
        "confidence": report.get("confidence", ""),
        "report": report,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description='Generic primary-protein selection (gene-agnostic).')
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()
    ctx = load_context(args.run_id)
    res = build(ctx)
    print(f"OK primary_selection  primary={res['primary_protein_id']}  rule={res['selection_rule']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
