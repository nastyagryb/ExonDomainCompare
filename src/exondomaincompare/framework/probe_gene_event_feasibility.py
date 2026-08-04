#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from exondomaincompare.framework.gene_config import (  # noqa: E402
    GeneConfig, GeneConfigError, PROJECT_ROOT, REQUIRES_EVENT_DETECTOR,
    GENE_CONFIG_DIR, GENE_DRAFT_DIR, capability_summary, load_gene_config_lenient,
    detector_for_analysis, load_event_detector_contract, load_core_analysis_contract,
)

# Overall feasibility levels (event-oriented, kept for backward compatibility).
SUPPORTED = "supported"
PARTIAL = "partially_supported"
NEEDS_DETECTOR = "requires_gene_specific_event_detector"
NOT_YET = "not_supported_yet"

# Core-vs-event classification (new model: core analysis is separate & optional-event).
EVENT_CONFIGURED = "event_configured"
EVENT_AUTO_POSSIBLE = "event_auto_possible"
EVENT_USER_DEFINED_POSSIBLE = "event_user_defined_possible"
EVENT_REQUIRES_USER_DEFINITION = "event_requires_user_definition"
EVENT_NOT_AVAILABLE = "event_not_available"

# Views core analysis alone can provide vs views that require an event region.
_DEFAULT_CORE_VIEWS = ["overview", "gene_models", "isoforms", "msa",
                       "domain_architecture", "synteny", "exon_domain_boundaries"]
_EVENT_VIEWS = ["event_region", "event_specific_comparison",
                "event_discriminating_columns", "event_specific_boundary_relation"]


def _core_views() -> List[str]:
    c = load_core_analysis_contract()
    v = c.get("core_views")
    return list(v) if isinstance(v, list) and v else list(_DEFAULT_CORE_VIEWS)


def _iso_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _find_configs_for_gene(gene: str) -> List[GeneConfig]:
    out: List[GeneConfig] = []
    g = gene.strip().lower()
    for d in (GENE_CONFIG_DIR, GENE_DRAFT_DIR):
        base = PROJECT_ROOT / d
        if not base.is_dir():
            continue
        for p in sorted(base.glob("*.yaml")):
            if p.name.startswith("TEMPLATE"):
                continue
            try:
                cfg = load_gene_config_lenient(p)
            except GeneConfigError:
                continue
            if cfg.gene_symbol.strip().lower() == g or g in cfg.analysis_id.lower():
                out.append(cfg)
    return out


def _event_state(cfg: Optional[GeneConfig]) -> Dict[str, Any]:
    if cfg is None:
        return {
            "event_configured": False,
            "event_type": None,
            "markers_provided": False,
            "labels": [],
            "detail": "No gene/event config exists for this gene yet.",
        }
    labels = cfg.event_labels
    markers = [lab for lab in labels if lab.get("marker_sequence")]
    event_defined = bool(cfg.event_type) and cfg.event_type != REQUIRES_EVENT_DETECTOR
    return {
        "event_configured": event_defined,
        "event_type": cfg.event_type or None,
        "markers_provided": bool(markers),
        "labels": [lab["id"] for lab in labels if lab["id"]],
        "detail": ("Event type and discriminating markers are configured."
                   if (event_defined and markers)
                   else "Event region / markers are not fully configured; a "
                        "gene-specific event detector is required."),
    }


def _detector_state(cfg: Optional[GeneConfig]) -> Dict[str, Any]:
    contract = load_event_detector_contract()
    required = sorted((contract.get("outputs", {}) or {}).keys())
    if cfg is None:
        return {
            "supported_detector_exists": False,
            "detector": None,
            "implementation": None,
            "contract_can_be_satisfied": False,
            "required_outputs": required,
            "detail": "No config, so no detector can be associated yet.",
        }
    det = detector_for_analysis(cfg.analysis_id)
    exists = det is not None
    return {
        "supported_detector_exists": exists,
        "detector": (det or {}).get("detector"),
        "implementation": (det or {}).get("implementation"),
        # The contract is satisfiable when a supported detector is registered that
        # projects/produces the required outputs. For FGFR2 this is the legacy adapter.
        "contract_can_be_satisfied": exists,
        "required_outputs": required,
        "detail": ("A supported event detector is registered and produces the "
                   "contract outputs."
                   if exists else
                   "No supported event detector is registered for this analysis; a "
                   "gene-specific detector must be implemented and registered in "
                   "configs/framework/event_detectors.yaml."),
    }


def classify_core_event(cfg: Optional[GeneConfig], event: Dict[str, Any],
                        detector: Dict[str, Any]) -> Dict[str, Any]:
    core_disabled = bool(cfg is not None and not cfg.core_gene_analysis_enabled)
    core = "disabled" if core_disabled else "potentially_supported"

    has_event = bool(cfg is not None and cfg.has_event)
    has_detector = bool(detector["supported_detector_exists"])
    mode = (cfg.event_analysis_mode if cfg is not None else "not_configured")
    status = (cfg.event_status if cfg is not None else "not_configured")

    # Does the config carry a concrete user-defined region (start+end aa)?
    user_region = False
    if cfg is not None:
        b = cfg.event_region_bounds
        user_region = b.get("region_start_aa") is not None and b.get("region_end_aa") is not None

    if has_event and has_detector:
        event_analysis = "configured"
        event_level = EVENT_CONFIGURED
    elif status == "user_defined":
        event_analysis = "user_defined"
        # A user-supplied region makes event views possible even without a
        # gene-specific detector; a bare user_defined status still needs a region.
        if has_detector:
            event_level = EVENT_CONFIGURED
        elif user_region:
            event_level = EVENT_USER_DEFINED_POSSIBLE
        else:
            event_level = EVENT_REQUIRES_USER_DEFINITION
    elif user_region:
        event_analysis = "user_defined"
        event_level = EVENT_USER_DEFINED_POSSIBLE
    elif mode == "auto":
        event_analysis = "auto_possible"
        event_level = EVENT_AUTO_POSSIBLE
    elif mode == "disabled":
        event_analysis = "disabled"
        event_level = EVENT_NOT_AVAILABLE
    else:
        event_analysis = "not_configured"
        event_level = EVENT_REQUIRES_USER_DEFINITION

    return {
        "core_gene_analysis": core,
        "event_analysis": event_analysis,
        "event_level": event_level,
        "available_if_core_runs": _core_views(),
        "requires_event": list(_EVENT_VIEWS),
    }


def classify(gene: str, cfg: Optional[GeneConfig], event: Dict[str, Any],
             detector: Dict[str, Any]) -> str:
    has_detector = detector["supported_detector_exists"]
    fully_configured = bool(cfg is not None and event["event_configured"]
                            and event["markers_provided"])
    if cfg is not None and cfg.runnable and fully_configured and has_detector:
        return SUPPORTED
    if fully_configured or has_detector:
        # configured but not enabled/runnable, or missing one requirement
        return PARTIAL
    if cfg is not None:
        # a config/draft exists but neither markers nor a supported detector are ready
        return NEEDS_DETECTOR
    return NOT_YET


def _recommended_action(level: str, gene: str) -> str:
    return {
        SUPPORTED: f"{gene} is configured and runnable. No further action needed.",
        PARTIAL: (f"{gene} has a full event/marker config but is not enabled as "
                  "runnable yet. Add adapter coverage and validate on a run, then "
                  "enable in the UI."),
        NEEDS_DETECTOR: (f"Define the {gene} event region and implement a "
                         "gene-specific event detector (and discriminating markers), "
                         "then register it in configs/framework/event_detectors.yaml "
                         "as 'supported' before this gene can be configured for analysis."),
        NOT_YET: (f"Create a draft gene config for {gene} "
                  f"(configs/genes/drafts/{gene}_draft.yaml), then re-run this probe. "
                  "An event detector will still be required."),
    }.get(level, "Review the report and decide next steps.")


def build_report(gene: str, reference_species: str) -> Dict[str, Any]:
    configs = _find_configs_for_gene(gene)
    primary = configs[0] if configs else None
    event = _event_state(primary)
    detector = _detector_state(primary)
    caps = capability_summary()
    level = classify(gene, primary, event, detector)
    core_event = classify_core_event(primary, event, detector)

    reusable = caps["reusable_modules"]
    gene_specific = caps["gene_specific_modules"]
    partial = caps["partial_or_event_specific_modules"]

    return {
        "schema_version": 1,
        "generated_at": _iso_now(),
        "gene_symbol": gene,
        "reference_species": reference_species,
        "feasibility": level,
        "runnable_now": bool(primary is not None and primary.runnable
                             and level == SUPPORTED and detector["supported_detector_exists"]),
        # New core-vs-event model: core analysis is separate and does not require
        # a configured event region.
        "core_gene_analysis": core_event["core_gene_analysis"],
        "event_analysis": core_event["event_analysis"],
        "event_analysis_level": core_event["event_level"],
        "available_if_core_runs": core_event["available_if_core_runs"],
        "requires_event": core_event["requires_event"],
        "event_detector": detector,
        # candidate identifiers / isoforms are NOT fetched here (offline probe);
        # a full probe would query NCBI/Ensembl. Report conservatively.
        "candidate_identifiers": {
            "status": "not_inspected_offline",
            "note": "This offline probe does not query external databases. A full "
                    "probe would resolve gene ids / RefSeq accessions per species.",
        },
        "known_protein_isoforms": {
            "status": "known" if event["labels"] else "not_inspected_offline",
            "isoform_labels": event["labels"],
        },
        "event_region": {
            "configured": event["event_configured"],
            "event_type": event["event_type"],
            "markers_provided": event["markers_provided"],
            "detail": event["detail"],
        },
        "domain_annotation": {
            "possible": True,
            "basis": "external_tool",
            "note": "InterProScan is gene-agnostic; domain architecture is reusable "
                    "once the gene's ortholog proteins are available.",
        },
        "synteny": {
            "possible": "in_principle",
            "note": "Reusable module, but requires the gene locus and per-species "
                    "genomic annotation sources; not gene-specific logic.",
        },
        "reusable_modules": reusable,
        "gene_specific_modules_needed": gene_specific,
        "partial_or_event_specific_modules": partial,
        "existing_configs": [
            {"analysis_id": c.analysis_id, "config": c.source_path,
             "status": c.status, "runnable": c.runnable}
            for c in configs
        ],
        "recommended_next_action": _recommended_action(level, gene),
        "disclaimer": "Core gene-level analysis (gene models, domains, synteny, all-exon "
                      "boundary distances) is generic and does NOT require an event region. "
                      "Event-specific views require a configured/detected event region; the "
                      "biological event detection is still FGFR2-specific.",
    }


def render_markdown(rep: Dict[str, Any]) -> str:
    def yesno(v: Any) -> str:
        return "yes" if v is True else ("no" if v is False else str(v))

    lines: List[str] = []
    lines.append(f"# Gene/event feasibility — {rep['gene_symbol']}")
    lines.append("")
    lines.append(f"- **Generated:** {rep['generated_at']}")
    lines.append(f"- **Gene symbol:** {rep['gene_symbol']}")
    lines.append(f"- **Reference species:** {rep['reference_species']}")
    lines.append(f"- **Feasibility:** `{rep['feasibility']}`")
    lines.append(f"- **Runnable now:** {yesno(rep['runnable_now'])}")
    lines.append("")
    lines.append("> " + rep["disclaimer"])
    lines.append("")

    lines.append("## Event region")
    er = rep["event_region"]
    lines.append(f"- Configured: {yesno(er['configured'])}")
    lines.append(f"- Event type: {er['event_type']}")
    lines.append(f"- Markers provided: {yesno(er['markers_provided'])}")
    lines.append(f"- {er['detail']}")
    lines.append("")

    lines.append("## Core vs event analysis")
    lines.append(f"- Core gene analysis: `{rep['core_gene_analysis']}` "
                 "(generic; no event region required)")
    lines.append(f"- Event analysis: `{rep['event_analysis']}` "
                 f"(level: `{rep['event_analysis_level']}`)")
    lines.append("- Available if core runs: " + ", ".join(rep["available_if_core_runs"]))
    lines.append("- Requires an event region: " + ", ".join(rep["requires_event"]))
    lines.append("")

    lines.append("## Event detector")
    ed = rep["event_detector"]
    lines.append(f"- Supported detector exists: {yesno(ed['supported_detector_exists'])}")
    lines.append(f"- Detector: {ed['detector'] or 'none'} "
                 f"(implementation: {ed['implementation'] or 'n/a'})")
    lines.append(f"- Contract outputs can be satisfied: {yesno(ed['contract_can_be_satisfied'])}")
    lines.append(f"- Required contract outputs: {', '.join(ed['required_outputs']) or 'none'}")
    lines.append(f"- {ed['detail']}")
    lines.append("")

    lines.append("## Identifiers & isoforms")
    ci = rep["candidate_identifiers"]
    ki = rep["known_protein_isoforms"]
    lines.append(f"- Candidate identifiers: `{ci['status']}` — {ci['note']}")
    lines.append(f"- Known protein isoforms: `{ki['status']}` "
                 f"({', '.join(ki['isoform_labels']) or 'none'})")
    lines.append("")

    lines.append("## Annotation feasibility")
    lines.append(f"- Domain annotation: {yesno(rep['domain_annotation']['possible'])} "
                 f"({rep['domain_annotation']['basis']}) — {rep['domain_annotation']['note']}")
    lines.append(f"- Synteny: {rep['synteny']['possible']} — {rep['synteny']['note']}")
    lines.append("")

    lines.append("## Module reuse")
    lines.append("**Reusable now (gene-agnostic / external tools):**")
    for m in rep["reusable_modules"]:
        lines.append(f"- {m}")
    lines.append("")
    lines.append("**Partial / event-type-specific (need per-gene inputs or mapping):**")
    for m in rep["partial_or_event_specific_modules"]:
        lines.append(f"- {m}")
    lines.append("")
    lines.append("**Gene-specific (must be implemented for a new gene):**")
    for m in rep["gene_specific_modules_needed"]:
        lines.append(f"- {m}")
    lines.append("")

    if rep["existing_configs"]:
        lines.append("## Existing configs for this gene")
        for c in rep["existing_configs"]:
            lines.append(f"- `{c['analysis_id']}` ({c['config']}) — "
                         f"status={c['status']}, runnable={yesno(c['runnable'])}")
        lines.append("")

    lines.append("## Recommended next action")
    lines.append(rep["recommended_next_action"])
    lines.append("")
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Probe gene/event feasibility (offline, read-only).")
    ap.add_argument("--gene", required=True, help="Gene symbol, e.g. FGFR1.")
    ap.add_argument("--reference-species", default="homo_sapiens",
                    help="Reference/control species (default: homo_sapiens).")
    ap.add_argument("--outdir", help="Output directory (default: artifacts/gene_feasibility/<GENE>).")
    args = ap.parse_args(argv)

    gene = args.gene.strip()
    outdir = Path(args.outdir) if args.outdir else (
        PROJECT_ROOT / "artifacts" / "gene_feasibility" / gene)
    if not outdir.is_absolute():
        outdir = PROJECT_ROOT / outdir
    # Safety: never write into the example freeze.
    freeze = PROJECT_ROOT / "results" / "final_30_until_interpro_prepare"
    if str(outdir.resolve()).startswith(str(freeze.resolve())):
        print("Refusing to write inside the example freeze.", file=sys.stderr)
        return 2

    rep = build_report(gene, args.reference_species)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "feasibility_report.json").write_text(
        json.dumps(rep, indent=2) + "\n", encoding="utf-8")
    (outdir / "feasibility_report.md").write_text(render_markdown(rep), encoding="utf-8")

    print(f"OK  gene={gene}  feasibility={rep['feasibility']}  runnable_now={rep['runnable_now']}")
    print(f"    out: {outdir}")
    print(f"      - feasibility_report.json")
    print(f"      - feasibility_report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
