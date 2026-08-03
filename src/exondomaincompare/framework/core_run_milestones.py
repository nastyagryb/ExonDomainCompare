#!/usr/bin/env python3
"""Canonical, file-based milestone checks for a Core Gene Analysis run.

This is the SINGLE source of truth for classifying a core-only run's state,
used by both:
  * src/exondomaincompare/framework/validate_core_gene_run.py (CLI), and
  * the webapp backend (webapp/backend/main.py) so a dashboard refresh uses the
    exact same milestone logic.

A run must NEVER look analysis-ready if required core outputs are missing. The
inferred status is derived purely from on-disk artefacts (not from a possibly
stale status.json), so empty/partial runs are always classified honestly.

Inferred status values (core-only):
  created_not_started        - run folder scaffolded, no core collection yet
  running                    - pre-cluster pipeline actively collecting/retrieving
  core_model_collection_failed - collection ran but produced no gene models
  incomplete                 - some required core outputs missing unexpectedly
  cluster_required           - core collection complete, cluster annotation pending
  cluster_running            - cluster jobs submitted / running
  post_interpro_incomplete   - cluster outputs present but domain build missing
  results_ready              - domain architecture + generic indices built
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Tuple

from exondomaincompare.config import discover_repository_root, load_config
from exondomaincompare.contracts import resolve_path_reference

PROJECT_ROOT = discover_repository_root(__file__)

RUNTIME_CONFIG = load_config(repository_root=PROJECT_ROOT)
PROJECT_ROOT = RUNTIME_CONFIG.repository_root
CORE_REL = "results/core_gene_analysis"


def _read_json(p: Path, default: Any = None) -> Any:
    try:
        return json.loads(Path(p).read_text(encoding="utf-8"))
    except Exception:
        return default


def _run_config(run_dir: Path) -> Dict[str, Any]:
    """Read the canonical run record with legacy compatibility fields as fallback."""
    compatibility = _read_json(run_dir / "run_config.json", {}) or {}
    canonical = _read_json(run_dir / "run.json", {}) or {}
    return {**compatibility, **canonical}


def _tsv_rows(p: Path) -> int:
    """Number of data rows (excluding header). 0 if missing/empty/header-only."""
    if not Path(p).is_file():
        return -1  # sentinel: file missing
    try:
        with open(p, "r", encoding="utf-8", newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            rows = list(reader)
        if not rows:
            return 0
        return max(0, len(rows) - 1)
    except Exception:
        return 0


def _fasta_records(p: Path) -> int:
    if not Path(p).is_file():
        return -1
    try:
        n = 0
        with open(p, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    n += 1
        return n
    except Exception:
        return 0


def _exists(p: Path) -> bool:
    return Path(p).exists()


def _boundary_applicability(run_dir: Path) -> Tuple[bool, str]:
    """Whether an exon–domain boundary analysis can be performed at all.

    Returns ``(applicable, reason)``. Falls back to applicable when the shared module
    cannot be imported, which keeps the previous, stricter behaviour rather than declaring
    something not applicable on the strength of an import error.
    """
    try:
        from exondomaincompare.shared_gene_analysis import analysis_availability as aa
    except Exception:
        return True, ""
    try:
        if not aa.has_core_tables(Path(run_dir)):
            return True, ""
        pre = aa.prerequisites(Path(run_dir))
    except Exception:
        return True, ""
    if pre.internal_coding_exon_boundary_count > 0:
        return True, ""
    reason = (aa.SINGLE_CODING_EXON if pre.coding_exon_count == 1
              else aa.NO_INTERNAL_BOUNDARIES)
    return False, aa.MESSAGES.get(reason, "")


def is_core_only_run(run_dir: Path) -> bool:
    rc = _run_config(Path(run_dir))
    if str(rc.get("run_mode", "")).lower() == "core_only_pilot":
        return True
    if rc.get("has_event") is False and rc.get("analysis_id"):
        return True
    st = _read_json(run_dir / "status.json", {}) or {}
    return str(st.get("run_mode", "")).lower() == "core_only_pilot"


def evaluate_core_run(run_dir: Path) -> Dict[str, Any]:
    """Return a structured milestone report for a core-only run directory."""
    run_dir = Path(run_dir)
    run_id = run_dir.name
    core = run_dir / CORE_REL
    rc = _run_config(run_dir)
    st = _read_json(run_dir / "status.json", {}) or {}
    report = _read_json(core / "core_gene_report.json", {}) or {}
    collection_report = _read_json(core / "core_model_collection_report.json", {}) or {}

    milestones: List[Dict[str, Any]] = []

    def add(mid: str, name: str, required: bool, present: List[str],
            missing: List[str], complete: bool, reason: str = "") -> None:
        milestones.append({
            "id": mid, "name": name, "required": required,
            "complete": complete, "present_files": present,
            "missing_files": missing, "reason": reason,
        })

    # -- Milestone 1: run setup ------------------------------------------- #
    setup_files = {
        "run_config.json": _exists(run_dir / "run_config.json"),
        "gene_config.yaml": _exists(run_dir / "gene_config.yaml"),
        "species_list.txt": _exists(run_dir / "species_list.txt"),
        "status.json": _exists(run_dir / "status.json"),
    }
    setup_missing = [k for k, v in setup_files.items() if not v]
    add("run_setup", "Run setup", True,
        [k for k, v in setup_files.items() if v], setup_missing, not setup_missing)

    # -- Milestone 2: model collection ------------------------------------ #
    gm = _tsv_rows(core / "gene_model_index.tsv")
    iso = _tsv_rows(core / "protein_isoform_index.tsv")
    m2_present, m2_missing = [], []
    (m2_present if gm > 0 else m2_missing).append("gene_model_index.tsv")
    (m2_present if iso > 0 else m2_missing).append("protein_isoform_index.tsv")
    m2_complete = gm > 0 and iso > 0
    add("model_collection", "Gene / protein model collection", True,
        m2_present, m2_missing, m2_complete,
        "" if m2_complete else "No gene models / protein isoforms collected.")

    # -- Milestone 3: primary protein FASTA ------------------------------- #
    faa = _fasta_records(core / "proteins_primary.faa")
    fasta_reason = ""
    if faa <= 0:
        fasta_reason = (report.get("fasta_reason")
                        or collection_report.get("fasta_reason")
                        or "no_primary_protein_selected")
    add("protein_fasta", "Primary protein FASTA", True,
        ["proteins_primary.faa"] if faa > 0 else [],
        [] if faa > 0 else ["proteins_primary.faa"], faa > 0, fasta_reason)

    # -- Milestone 4: exon map (optional; reason allowed) ----------------- #
    exon = _tsv_rows(core / "exon_protein_map.tsv")
    exon_reason = "" if exon > 0 else (report.get("exon_map_reason")
                                       or collection_report.get("exon_map_reason")
                                       or "exon_map_unavailable")
    add("exon_map", "Exon → protein map", False,
        ["exon_protein_map.tsv"] if exon > 0 else [],
        [] if exon > 0 else ["exon_protein_map.tsv"], exon > 0, exon_reason)

    # -- Milestone 5: synteny (optional; reason allowed) ------------------ #
    syn = _tsv_rows(core / "synteny_neighbors.tsv")
    syn_reason = "" if syn > 0 else (report.get("synteny_reason") or "synteny_unavailable")
    add("synteny", "Synteny neighbours", False,
        ["synteny_neighbors.tsv"] if syn > 0 else [],
        [] if syn > 0 else ["synteny_neighbors.tsv"], syn > 0, syn_reason)

    # -- Milestone 6: cluster input --------------------------------------- #
    cluster_fasta_rel = rc.get("cluster_input_fasta") or rc.get("primary_fasta_path") or ""
    cluster_fasta_ok = False
    if cluster_fasta_rel:
        cp = resolve_path_reference(
            str(cluster_fasta_rel), repository_root=PROJECT_ROOT, run_root=run_dir)
        cluster_fasta_ok = cp.is_file() and _fasta_records(cp) > 0
    cluster_command = f".venv/bin/edc cluster roundtrip --run-id {run_id}"
    add("cluster_input", "Cluster input FASTA", True,
        ["cluster_input_fasta"] if cluster_fasta_ok else [],
        [] if cluster_fasta_ok else ["cluster_input_fasta"], cluster_fasta_ok,
        "" if cluster_fasta_ok else "cluster_input_fasta not set or FASTA empty.")

    # -- Milestone 7: post-domain ----------------------------------------- #
    # Only *applicable* outputs are required. Demanding a boundary row unconditionally
    # made a single-coding-exon gene look like a broken pipeline: chicken MC1R has one
    # coding exon, hence zero internal boundaries, hence nothing a boundary table could
    # contain — and the run was reported as post_cluster_partial for having the exon
    # structure it has. Applicability is decided from the real coding-exon count.
    dom = _tsv_rows(core / "domain_features.tsv")
    tm = _tsv_rows(core / "tm_features.tsv")
    bnd = _tsv_rows(core / "exon_domain_boundary_distances.tsv")
    generic_views = _exists(run_dir / "website_indices" / "generic" / "available_views.json")
    dom_arch = _read_json(run_dir / "website_indices" / "generic" / "domain_architecture_index.json", {}) or {}
    dom_arch_available = bool(dom_arch.get("available"))
    boundary_applicable, boundary_note = _boundary_applicability(run_dir)
    m7_present, m7_missing = [], []
    (m7_present if dom > 0 else m7_missing).append("domain_features.tsv")
    if boundary_applicable:
        (m7_present if bnd > 0 else m7_missing).append("exon_domain_boundary_distances.tsv")
    (m7_present if _exists(core / "core_gene_report.json") else m7_missing).append("core_gene_report.json")
    (m7_present if generic_views else m7_missing).append("generic_website_indices")
    # tm_features may legitimately be empty (soluble proteins); informational only
    m7_complete = (dom > 0 and dom_arch_available
                   and (bnd > 0 or not boundary_applicable))
    if m7_complete:
        m7_reason = ""
    elif not boundary_applicable and not (dom > 0 and dom_arch_available):
        m7_reason = ("Domain architecture not built yet (requires cluster "
                     "InterProScan/pyTMHMM outputs).")
    else:
        m7_reason = ("Domain architecture / boundaries not built yet (requires cluster "
                     "InterProScan/pyTMHMM outputs).")
    add("post_domain", "Post-domain analysis", True,
        m7_present, m7_missing, m7_complete, m7_reason)
    if not boundary_applicable and boundary_note:
        milestones[-1]["not_applicable"] = ["exon_domain_boundary_distances.tsv"]
        milestones[-1]["not_applicable_reason"] = boundary_note

    # -- infer status ----------------------------------------------------- #
    cluster_status = str(st.get("cluster_analysis_status") or "").lower()
    cluster_jobs = st.get("cluster_jobs") or {}
    ips_out = list((run_dir / "results" / "14_interproscan" / "primary" / "output").glob("*.tsv")) \
        if (run_dir / "results" / "14_interproscan" / "primary" / "output").is_dir() else []
    ips_out = [p for p in ips_out if "input" not in p.name or p.stat().st_size > 0]

    m = {ms["id"]: ms for ms in milestones}
    missing_required = [f"{ms['id']}:{f}" for ms in milestones if ms["required"]
                        for f in ms["missing_files"]]
    missing_optional = [f"{ms['id']}:{f}" for ms in milestones if not ms["required"]
                        for f in ms["missing_files"]]

    # Explicit signals written by the runner (single source of truth for the
    # pre-cluster phase). A run that recorded a terminal failure or is actively
    # running must be classified honestly BEFORE the artefact heuristics — a
    # blocked run creates no core artefacts, so the "no files => not started"
    # heuristic would otherwise mislabel a genuine failure as "created".
    raw_status = str(st.get("status") or "").lower()
    raw_step = str(st.get("current_step") or "").lower()
    pre_status = str(st.get("pre_interpro_status") or "").lower()
    failed_step = str(st.get("failed_step") or "")
    failed_reason = str(st.get("failed_reason") or st.get("error") or "")
    terminal_failed = bool(
        failed_step
        or raw_status in ("failed", "error", "incomplete", "core_model_collection_failed")
        or raw_status.endswith("_failed")
        or pre_status in ("failed", "error")
        or "blocked" in raw_step
    )
    is_running = (not terminal_failed) and (
        raw_status == "running" or pre_status == "running"
        or raw_step in ("collecting_models", "retrieving_models", "core_collection_running")
    )

    if not m["run_setup"]["complete"]:
        inferred = "created_not_started"
    elif not m["model_collection"]["complete"]:
        if terminal_failed:
            # distinguish a FASTA-stage incompleteness from a collection failure.
            inferred = "incomplete" if failed_step == "core_primary_fasta" else "core_model_collection_failed"
        elif is_running:
            inferred = "running"
        else:
            # setup exists but no models and no explicit signal: distinguish
            # "never started" vs "started but produced nothing".
            started = bool(collection_report) or bool(report) or _exists(core)
            inferred = "core_model_collection_failed" if started else "created_not_started"
    elif not m["protein_fasta"]["complete"] or not m["cluster_input"]["complete"]:
        inferred = "incomplete"
    elif m["post_domain"]["complete"]:
        inferred = "results_ready"
    elif ips_out or cluster_status in ("fetched", "complete", "completed"):
        inferred = "post_interpro_incomplete"
    elif cluster_jobs or cluster_status in ("submitted", "running", "queued", "pending"):
        inferred = "cluster_running"
    else:
        inferred = "cluster_required"

    next_action = {
        "created_not_started": "run_core_collection",
        "running": "wait_pre_interpro",
        "core_model_collection_failed": "inspect_logs",
        "incomplete": "inspect_logs",
        "cluster_required": "run_cluster_roundtrip_command",
        "cluster_running": "wait_cluster",
        "post_interpro_incomplete": "run_core_post",
        "results_ready": "open_results",
    }[inferred]

    completed = [ms["id"] for ms in milestones if ms["complete"]]
    logs = []
    log_dir = run_dir / "logs"
    if log_dir.is_dir():
        for path in sorted(log_dir.glob("*.log")):
            try:
                logs.append(str(path.relative_to(PROJECT_ROOT)))
            except ValueError:
                logs.append("run:" + str(path.relative_to(run_dir)))

    return {
        "run_id": run_id,
        "is_core_only": is_core_only_run(run_dir),
        "run_mode": rc.get("run_mode") or st.get("run_mode") or "",
        "analysis_id": rc.get("analysis_id") or report.get("analysis_id") or "",
        "gene_symbol": rc.get("gene_symbol") or report.get("gene_symbol") or "",
        "has_event": bool(rc.get("has_event", False)),
        "milestones": milestones,
        "completed_milestones": completed,
        "missing_required": missing_required,
        "missing_optional": missing_optional,
        "inferred_status": inferred,
        "suggested_next_action": next_action,
        "cluster_command": cluster_command,
        "failed_stage": failed_step,
        "failed_reason": failed_reason,
        "failed_species": str(st.get("failed_species") or ""),
        "run_status_raw": raw_status,
        "detail": str(st.get("detail") or ""),
        "logs": logs,
        "counts": {
            "gene_models": max(gm, 0), "protein_isoforms": max(iso, 0),
            "primary_proteins": max(faa, 0), "exon_map_rows": max(exon, 0),
            "synteny_neighbors": max(syn, 0), "domain_features": max(dom, 0),
            "tm_features": max(tm, 0), "boundary_rows": max(bnd, 0),
        },
    }
