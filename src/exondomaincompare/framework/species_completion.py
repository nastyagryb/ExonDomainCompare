#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

STATE_AVAILABLE = "available"
STATE_PENDING = "pending"
STATE_MISSING = "missing"
STATE_PARTIAL = "partial"
STATE_UNAVAILABLE = "unavailable"
STATE_NOT_APPLICABLE = "not_applicable"
#: Returned cluster results that no longer describe the run's current proteins.
STATE_STALE = "stale"

#: States that are settled answers and therefore never block a run's completion.
RESOLVED_STATES = (STATE_AVAILABLE, STATE_UNAVAILABLE, STATE_NOT_APPLICABLE)

RUN_READY = "results_ready"
RUN_PARTIAL = "post_cluster_partial"
RUN_PRE_CLUSTER = "pre_cluster_ready"
RUN_FAILED = "failed"

# The analyses a species must have to be considered complete post-cluster. Ordered as
# the pipeline produces them, so the first non-available entry is the blocking one.
REQUIRED_ANALYSES = (
    "cluster_outputs",
    "interproscan",
    "pytmhmm",
    "representative_domains",
    "candidate_domain_context",
    "boundary_analysis",
    "protein_coordinate_model",
    "website_indices",
)


def _read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        return []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _read_json(path: Path, default: Any = None) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return default


def _rows_for(rows: List[Dict[str, str]], species_id: str,
              protein_id: Optional[str] = None) -> List[Dict[str, str]]:
    out = [r for r in rows if str(r.get("species_id") or "") == species_id]
    if protein_id:
        out = [r for r in out if str(r.get("protein_id") or "") == protein_id]
    return out


def build_species_completion(run_dir: Path) -> Dict[str, Dict[str, Any]]:
    run_dir = Path(run_dir)
    core = run_dir / "results" / "core_gene_analysis"
    status_doc = _read_json(run_dir / "status.json", {}) or {}
    cluster_complete = str(status_doc.get("cluster_status")
                           or status_doc.get("cluster_analysis_status") or "").lower() == "complete"

    isoforms = _read_tsv(core / "protein_isoform_index.tsv")
    species = sorted({str(r.get("species_id") or "") for r in isoforms
                      if r.get("species_id")})

    domain_rows = _read_tsv(core / "domain_features.tsv")
    tm_rows = _read_tsv(core / "tm_features.tsv")
    boundary_rows = _read_tsv(core / "exon_domain_boundary_distances.tsv")
    candidate_rows = _read_tsv(core / "event_candidate_regions.tsv")
    interpro_rows = _read_tsv(core / "interpro_annotations.tsv")

    model = _read_json(
        run_dir / "website_indices" / "generic" / "protein_coordinate_model.json", {}) or {}
    models_by_species = {str(m.get("species_id") or ""): m
                         for m in (model.get("models") or [])}

    try:
        from exondomaincompare.framework.primary_resolution import resolve_primaries
        primaries = resolve_primaries(core)
    except Exception:
        primaries = {}

    # Whether the returned cluster results still describe the proteins this run
    # analyses. A repaired coordinate model can change a sequence while every
    # output file stays in place, and domain calls for a superseded sequence must
    # not reach a reader, so this is checked by sequence digest rather than by
    # file presence.
    cluster_freshness: Dict[str, Any] = {}
    stale_proteins: set = set()
    try:
        from exondomaincompare.shared_gene_analysis.cluster_output_freshness import (
            FRESH, INCOMPLETE, MISSING, STALE, evaluate,
        )
        cluster_freshness = evaluate(run_dir)
        stale_proteins = set(cluster_freshness["interproscan"]["mismatched"])
        # Only a positive sequence mismatch is a defect. A payload without a
        # digest cannot be checked, so the run's own cluster status decides —
        # calling "unverifiable" a stale result would withhold current findings.
        _CLUSTER_STATE = {
            FRESH: STATE_AVAILABLE,
            STALE: STATE_STALE,
            INCOMPLETE: STATE_AVAILABLE if cluster_complete else STATE_PENDING,
            MISSING: STATE_AVAILABLE if cluster_complete else STATE_PENDING,
        }
        cluster_state = _CLUSTER_STATE.get(cluster_freshness["status"], STATE_PENDING)
    except Exception:
        cluster_state = STATE_AVAILABLE if cluster_complete else STATE_PENDING

    # How many internal coding-exon boundaries each species' primary model actually has.
    # Derived from the coding-exon count, never from the total exon count.
    internal_boundaries: Dict[str, int] = {}
    try:
        from exondomaincompare.shared_gene_analysis.analysis_availability import has_core_tables, prerequisites
        if has_core_tables(run_dir):
            for sid in species:
                internal_boundaries[sid] = \
                    prerequisites(run_dir, sid).internal_coding_exon_boundary_count
    except Exception:
        internal_boundaries = {}

    out: Dict[str, Dict[str, Any]] = {}
    for sid in species:
        primary = primaries.get(sid) or {}
        pid = primary.get("protein_id") or ""
        m = models_by_species.get(sid) or {}

        def _feature_state(rows: List[Dict[str, str]]) -> str:
            if _rows_for(rows, sid, pid):
                return STATE_AVAILABLE
            # No rows for this species' primary. Before the cluster has returned this
            # is genuinely pending; afterwards it is a gap.
            return STATE_MISSING if cluster_complete else STATE_PENDING

        interproscan = _feature_state(interpro_rows or domain_rows)

        # "pyTMHMM found no transmembrane region" is a result, not a missing analysis.
        # A soluble protein such as TP53 legitimately has zero TM regions, and counting
        # that as incomplete would mark a finished run partial forever. The model's
        # tm_analysis block already separates "ran" from "found something".
        tm_analysis = m.get("tm_analysis") or {}
        if tm_analysis.get("pending"):
            pytmhmm = STATE_PENDING
        elif tm_analysis.get("performed"):
            pytmhmm = STATE_AVAILABLE
        else:
            pytmhmm = _feature_state(tm_rows)
        rep_domains = (STATE_AVAILABLE if (m.get("representative_domains") or [])
                       else _feature_state(domain_rows))
        candidates = (STATE_AVAILABLE if _rows_for(candidate_rows, sid)
                      else STATE_UNAVAILABLE)
        boundaries = (STATE_AVAILABLE if (m.get("exon_boundaries") or [])
                      else _feature_state(boundary_rows))
        # Zero internal coding-exon boundaries settles the question before the domain
        # layer is consulted: there is nothing to relate a domain to, so neither a missing
        # boundary table nor a missing domain table is a gap.
        if boundaries != STATE_AVAILABLE and internal_boundaries.get(sid, 1) == 0:
            boundaries = STATE_NOT_APPLICABLE
        elif boundaries != STATE_AVAILABLE and rep_domains != STATE_AVAILABLE:
            boundaries = "blocked_missing_domains"

        model_state = (
            STATE_AVAILABLE if str(m.get("status") or "") == STATE_AVAILABLE
            else STATE_PENDING if str(m.get("status") or "").startswith("pending")
            else STATE_MISSING if not m else STATE_PARTIAL
        )
        indices = (STATE_AVAILABLE
                   if (run_dir / "website_indices" / "generic").is_dir() and m
                   else STATE_MISSING)

        record = {
            "species_id": sid,
            "primary_protein": pid,
            "primary_transcript": primary.get("transcript_id") or "",
            "cluster_outputs": (STATE_STALE if pid in stale_proteins
                                else cluster_state),
            "interproscan": interproscan,
            "pytmhmm": pytmhmm,
            "representative_domains": rep_domains,
            "candidate_domain_context": candidates,
            "boundary_analysis": boundaries,
            "protein_coordinate_model": model_state,
            "website_indices": indices,
        }
        record["domain_architecture"] = rep_domains
        record["boundary"] = boundaries
        record["cluster_output_status"] = cluster_freshness.get("status") or ""
        blocking = [k for k in REQUIRED_ANALYSES
                    if record.get(k) not in RESOLVED_STATES]
        record["complete"] = not blocking
        record["blocking_analyses"] = blocking
        out[sid] = record
    return out


def aggregate_run_status(completion: Dict[str, Dict[str, Any]],
                         cluster_complete: bool = True,
                         failed_reasons: Optional[List[str]] = None,
                         ) -> Tuple[str, List[str]]:
    reasons = list(failed_reasons or [])
    if reasons:
        return RUN_FAILED, reasons

    if not completion:
        return RUN_PRE_CLUSTER, ["no species completion records"]

    complete = [s for s, r in completion.items() if r.get("complete")]
    incomplete = [s for s, r in completion.items() if not r.get("complete")]

    for sid in sorted(incomplete):
        blocking = completion[sid].get("blocking_analyses") or []
        reasons.append(f"{sid}: " + ", ".join(
            f"{k}={completion[sid].get(k)}" for k in blocking) if blocking
            else f"{sid}: incomplete")

    if not complete:
        return (RUN_PRE_CLUSTER if not cluster_complete else RUN_PARTIAL), reasons
    if incomplete:
        return RUN_PARTIAL, reasons
    return RUN_READY, []


def species_status_summary(completion: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    return {
        sid: {"domain_architecture": r.get("representative_domains"),
              "boundary": r.get("boundary_analysis"),
              "protein_coordinate_model": r.get("protein_coordinate_model"),
              "primary_protein": r.get("primary_protein"),
              "complete": r.get("complete")}
        for sid, r in sorted(completion.items())
    }
