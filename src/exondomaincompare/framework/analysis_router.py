from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional

try:  # gene_config is the config-driven registry we build on top of.
    from . import gene_config as _gc  # type: ignore
except Exception:  # pragma: no cover - allow standalone import
    try:
        from exondomaincompare.framework import gene_config as _gc  # type: ignore
    except Exception:  # pragma: no cover
        _gc = None


WORKFLOW_VALIDATED = "validated_event_analysis"
WORKFLOW_SHARED = "shared_exploratory"

EVENT_LAYER_VALIDATED = "validated"
EVENT_LAYER_EXPLORATORY = "exploratory"

# The one validated gene that ships a frozen pipeline. This is intentionally an
# explicit constant *and* cross-checked against the config registry below, so a
# missing/broken config can never silently downgrade FGFR2 to exploratory.
_CANONICAL_VALIDATED = "FGFR2"


@dataclass
class GeneWorkflow:
    gene_symbol: str
    workflow: str                     # WORKFLOW_VALIDATED | WORKFLOW_SHARED
    event_layer: str                  # EVENT_LAYER_VALIDATED | EVENT_LAYER_EXPLORATORY
    is_validated: bool
    has_event: bool
    support_level: str
    case_study: str
    analysis_id: str
    creator: str                      # which script creates + drives the run
    mode: str = "auto"
    reason: str = ""
    ui_note: str = ""
    analysis_id_candidates: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def normalize_gene_symbol(symbol: Optional[str]) -> str:
    return str(symbol or "").strip().upper()


def _validated_registry() -> Dict[str, Dict[str, Any]]:
    registry: Dict[str, Dict[str, Any]] = {}
    if _gc is not None:
        try:
            discovered = _gc.discover_analyses()
            for detail in discovered.get("supported_detail", []) or []:
                sym = normalize_gene_symbol(detail.get("gene_symbol"))
                if sym:
                    registry[sym] = detail
        except Exception:
            registry = {}
    # Safety net: FGFR2 is always validated even if config discovery failed.
    if _CANONICAL_VALIDATED not in registry:
        registry[_CANONICAL_VALIDATED] = {
            "analysis_id": "FGFR2_IIIb_IIIc",
            "gene_symbol": _CANONICAL_VALIDATED,
            "support_level": "validated_event_analysis",
            "has_event": True,
        }
    return registry




def list_validated_genes() -> List[str]:
    return sorted(_validated_registry().keys())


def resolve_gene_workflow(
    gene_symbol: Optional[str],
    mode: str = "auto",
    species: Optional[List[str]] = None,
) -> GeneWorkflow:
    sym = normalize_gene_symbol(gene_symbol)
    registry = _validated_registry()

    if sym in registry:
        detail = registry[sym]
        analysis_id = str(detail.get("analysis_id") or f"{sym}_validated")
        return GeneWorkflow(
            gene_symbol=sym,
            workflow=WORKFLOW_VALIDATED,
            event_layer=EVENT_LAYER_VALIDATED,
            is_validated=True,
            has_event=bool(detail.get("has_event", True)),
            support_level=str(detail.get("support_level") or "validated_event_analysis"),
            case_study="FGFR2_IIIb_IIIc" if sym == _CANONICAL_VALIDATED else f"{sym}_validated",
            analysis_id=analysis_id,
            creator="run_pre_interpro_for_run.py",
            mode=mode,
            reason=f"{sym} has an active, runnable validated config with a supported event detector.",
            ui_note=(
                "Routes to the frozen, validated FGFR2 IIIb/IIIc workflow "
                "(events, cassette, human comparison)."
                if sym == _CANONICAL_VALIDATED
                else f"Routes to the validated workflow configured for {sym}."
            ),
            analysis_id_candidates=[analysis_id],
        )

    # Everything else → shared, gene-agnostic exploratory workflow.
    return GeneWorkflow(
        gene_symbol=sym,
        workflow=WORKFLOW_SHARED,
        event_layer=EVENT_LAYER_EXPLORATORY,
        is_validated=False,
        has_event=False,
        support_level="core_only_pilot",
        case_study=f"{sym}_core_only_pilot" if sym else "core_only_pilot",
        analysis_id=f"{sym}_core_only_pilot" if sym else "core_only_pilot",
        creator="run_core_gene_analysis.py",
        mode=mode,
        reason=(
            f"No active validated config for {sym}; routed to the shared exploratory workflow."
            if sym else "No gene symbol provided; defaulting to the shared exploratory workflow."
        ),
        ui_note=(
            "Routes to the shared exploratory gene workflow: same run stages, but an "
            "exploratory isoform-difference evidence layer instead of validated events."
        ),
        analysis_id_candidates=[
            f"{sym}_core_only_pilot" if sym else "core_only_pilot",
            f"{sym}_draft" if sym else "",
        ],
    )


def _main(argv: Optional[List[str]] = None) -> int:  # pragma: no cover - CLI helper
    import argparse
    import json

    ap = argparse.ArgumentParser(description="Resolve the analysis workflow for a gene symbol.")
    ap.add_argument("gene", help="Gene symbol, e.g. FGFR2 or FGFR1.")
    ap.add_argument("--mode", default="auto")
    args = ap.parse_args(argv)
    print(json.dumps(resolve_gene_workflow(args.gene, mode=args.mode).to_dict(), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
