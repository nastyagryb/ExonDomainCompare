#!/usr/bin/env python3
"""The one place a run's status.json is decided.

Two pipelines finish a run, and each used to write the run-level status from a
field it happened to have at hand. The post-cluster runner derived the overall
status from the persisted ``pre_interpro_status`` — so a pre-cluster stage that
had failed once and was then repaired kept the run out of ``results_ready``
forever, even though every scientific view had current data and the roundtrip's
own readiness evaluation said so. That is the FGFR2 human/cat readiness
regression: a status field describing a superseded attempt outranked the run's
artifacts.

This module answers the question from the artifacts instead, for both layouts:

* runs with the event-pipeline closure are judged by
  ``shared_gene_analysis.run_availability.readiness``;
* generic core runs are judged by ``framework.species_completion``.

Both are gated on the cluster results actually describing the run's current
proteins, so returned annotations for a superseded sequence can never carry a run
to ``results_ready``.

Usage::

    python -m exondomaincompare.shared_gene_analysis.finalize_run_status --run-id <run_id>
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from exondomaincompare.config import discover_repository_root, load_config
from exondomaincompare.contracts import stamp_payload
from exondomaincompare.runs.registry import resolve_run_record

RUNTIME_CONFIG = load_config(repository_root=discover_repository_root(__file__))
ROOT = RUNTIME_CONFIG.repository_root
RUNS_ROOT = RUNTIME_CONFIG.runs_root

RESULTS_READY = "results_ready"
POST_CLUSTER_PARTIAL = "post_cluster_partial"
CLUSTER_REQUIRED = "cluster_required"
FAILED = "failed"

# Failure notes describe one attempt. Once the artifacts of that stage are valid
# again, keeping them would contradict the state the run is published with.
_STALE_FAILURE_FIELDS = ("failed_reason", "failed_step", "failed_species",
                         "error", "detail")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read_json(path: Path, default: Any = None) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return default


def _cluster_state(run_dir: Path) -> Tuple[str, Dict[str, Any]]:
    from .cluster_output_freshness import evaluate
    report = evaluate(run_dir)
    return report["status"], report


def _event_pipeline_verdict(run_dir: Path, status: Dict[str, Any],
                            ) -> Optional[Tuple[str, str, List[str]]]:
    """Verdict for a run with the event-pipeline closure, or None."""
    from .run_availability import models_run, readiness
    if not models_run(run_dir):
        return None
    run_config = _read_json(run_dir / "run_config.json", {}) or {}
    gene = str(run_config.get("gene_symbol") or "FGFR2").upper()
    verdict = readiness(run_dir,
                        n_species=int(status.get("species_count") or 0),
                        has_event=gene == "FGFR2")
    if verdict.ready:
        return RESULTS_READY, verdict.reason, []
    blocking = [f"{v.view}={v.state}" for v in verdict.blocking]
    return POST_CLUSTER_PARTIAL, verdict.reason, blocking


def _core_verdict(run_dir: Path) -> Optional[Tuple[str, str, List[str]]]:
    """Verdict for a generic core run, or None when the layout does not apply."""
    from exondomaincompare.framework.species_completion import (
        aggregate_run_status, build_species_completion,
    )
    completion = build_species_completion(run_dir)
    if not completion:
        return None
    status_doc = _read_json(run_dir / "status.json", {}) or {}
    cluster_complete = str(status_doc.get("cluster_status")
                           or status_doc.get("cluster_analysis_status")
                           or "").lower() == "complete"
    state, reasons = aggregate_run_status(completion, cluster_complete=cluster_complete)
    mapped = {"results_ready": RESULTS_READY,
              "post_cluster_partial": POST_CLUSTER_PARTIAL,
              "pre_cluster_ready": CLUSTER_REQUIRED,
              "failed": FAILED}.get(state, POST_CLUSTER_PARTIAL)
    reason = ("every species reached the requested analysis stage"
              if mapped == RESULTS_READY else "; ".join(reasons))
    return mapped, reason, reasons


def evaluate_run(run_dir: Path) -> Dict[str, Any]:
    """Derive the canonical run status without writing anything."""
    status = _read_json(run_dir / "status.json", {}) or {}
    cluster_status, cluster_report = _cluster_state(run_dir)

    verdict = _event_pipeline_verdict(run_dir, status) or _core_verdict(run_dir)
    if verdict is None:
        return {"run_id": run_dir.name, "status": status.get("status") or "",
                "reason": "no availability contract describes this run layout",
                "blocking": [], "cluster_outputs": cluster_status,
                "cluster_report": cluster_report, "decided": False}

    state, reason, blocking = verdict
    # Returned annotations that no longer describe the run's proteins cannot carry
    # it to a ready state, whatever the view-level evaluation says about files.
    if state == RESULTS_READY and cluster_status == "stale":
        state = POST_CLUSTER_PARTIAL
        reason = ("returned cluster results were scored for superseded protein "
                  "sequences and must be recomputed")
        blocking = [f"cluster_outputs={cluster_status}"]

    return {"run_id": run_dir.name, "status": state, "reason": reason,
            "blocking": blocking, "cluster_outputs": cluster_status,
            "cluster_report": cluster_report, "decided": True}


def finalize(run_dir: Path, dry_run: bool = False) -> Dict[str, Any]:
    """Write the canonical run status derived from the run's artifacts."""
    report = evaluate_run(run_dir)
    if not report["decided"] or dry_run:
        return report

    path = run_dir / "status.json"
    status = _read_json(path, {}) or {}
    ready = report["status"] == RESULTS_READY
    status.update({
        "status": report["status"],
        "run_status": report["status"],
        "current_step": report["status"],
        "readiness_reason": report["reason"],
        "explorable": ready,
        "next_action": "open_results" if ready else status.get("next_action") or "",
        "cluster_output_status": report["cluster_outputs"],
        "status_source": "shared_gene_analysis.finalize_run_status",
        "last_updated": _now(),
    })
    if ready:
        # The stage's artifacts are valid, so a note from an earlier attempt would
        # contradict the state the run is now published with.
        status.update({"pre_interpro_status": "complete",
                       "post_interpro_status": "complete",
                       "website_indices_status": "complete",
                       "cluster_analysis_status": "complete",
                       "cluster_fetch_status": "complete",
                       "interproscan_status": "complete",
                       "pytmhmm_status": "complete",
                       "next_actions": []})
        for field in _STALE_FAILURE_FIELDS:
            status.pop(field, None)
        status.pop("blocking_analyses", None)
        roundtrip = status.get("cluster_roundtrip")
        if isinstance(roundtrip, dict):
            roundtrip.update({"phase": "complete", "reason": "",
                              "updated_at": _now()})
    else:
        status["blocking_analyses"] = report["blocking"]
    status = stamp_payload(
        status, payload_type="status", run_id=run_dir.name,
        dataset_id=run_dir.name, profile=RUNTIME_CONFIG.public_identity(),
        generator="src/exondomaincompare/shared_gene_analysis/finalize_run_status.py",
    )
    path.write_text(json.dumps(status, indent=2))
    report["written"] = True
    return report


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
    run_dir = record.path if record else (RUNS_ROOT / args.run_id).resolve()
    if not run_dir.is_dir():
        ap.error(f"no such run: {run_dir}")
    if record and record.read_only and not args.dry_run:
        ap.error("registered legacy run is read-only; copy it before finalizing")
    report = finalize(run_dir, dry_run=args.dry_run)
    report.pop("cluster_report", None)
    print(json.dumps(report, indent=2))
    return 0 if report.get("status") == RESULTS_READY else 1


if __name__ == "__main__":
    raise SystemExit(main())
