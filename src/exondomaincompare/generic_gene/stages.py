"""Canonical shared pipeline-stage contract (single source of truth).

Every gene run uses the SAME stage folders and the SAME conceptual per-stage
outputs. Only the *event layer* implementation differs:

  * FGFR2 -> validated_fgfr2_iiib_iiic  (IIIb/IIIc, markers, cassette, human comp.)
  * other -> exploratory_event_evidence (isoform-difference candidates only)

FGFR2's historical folders remain unchanged to preserve the validated pipeline
and freeze.
"""
from __future__ import annotations

# Canonical, gene-agnostic stage folders (order matters for display).
STAGES = [
    "00_run_setup",
    "01_input_validation",
    "02_gene_models",
    "03_orthology_or_model_selection",
    "04_event_evidence",
    "05_event_region_detection",
    "06_coordinate_mapping",
    "07_msa",
    "08_synteny",
    "09_qc",
    "10_figures_pre_domain",
    "14_interproscan",
    "15_domain_architecture",
    "16_final_analyses",
]

# Extra shared outputs that are not numbered stages.
WEBSITE_INDICES = "website_indices"

# Event-layer routing constants.
EVENT_LAYER_VALIDATED = "validated_fgfr2_iiib_iiic"
EVENT_LAYER_EXPLORATORY = "exploratory_event_evidence"
PIPELINE_SHARED = "shared_gene_pipeline"


def event_layer_for_gene(gene_symbol: str) -> dict:
    """Routing: same pipeline for all; only the event layer differs."""
    if (gene_symbol or "").upper() == "FGFR2":
        return {
            "pipeline_type": PIPELINE_SHARED,
            "event_layer_type": EVENT_LAYER_VALIDATED,
            "event_status": "validated",
            "event_analysis_enabled": True,
            "support_level": "validated_event_analysis",
        }
    return {
        "pipeline_type": PIPELINE_SHARED,
        "event_layer_type": EVENT_LAYER_EXPLORATORY,
        "event_status": "exploratory_candidates_only",
        "event_analysis_enabled": False,
        "support_level": "generic_gene_analysis_experimental",
    }
