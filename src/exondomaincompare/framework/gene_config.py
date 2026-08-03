#!/usr/bin/env python3
"""Gene/event analysis configuration loader.

Loads a gene_config.yaml (e.g. configs/genes/FGFR2_IIIb_IIIc.yaml) and exposes a
small, stable, gene/event-agnostic view of it. This is a DESCRIPTION layer only:
it does not run or change the FGFR2 pipeline and never touches the example freeze.

Usage (validation CLI):

    python -m exondomaincompare.framework.gene_config \
        --config configs/genes/FGFR2_IIIb_IIIc.yaml --validate

Programmatic:

    from scripts.framework.gene_config import load_gene_config, DEFAULT_GENE_CONFIG
    cfg = load_gene_config("configs/genes/FGFR2_IIIb_IIIc.yaml")
    cfg.analysis_id      # "FGFR2_IIIb_IIIc"
    cfg.event_labels     # [{"id": "IIIb", ...}, {"id": "IIIc", ...}]
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from exondomaincompare.config import discover_repository_root

try:
    import yaml  # PyYAML
except Exception as exc:  # pragma: no cover - import guard
    yaml = None
    _YAML_IMPORT_ERROR = exc
else:
    _YAML_IMPORT_ERROR = None

PROJECT_ROOT = discover_repository_root(__file__)

# The canonical, validated analysis. Used as the backward-compatible default for
# runs / datasets that predate the gene/event config layer.
DEFAULT_GENE_CONFIG = "configs/genes/FGFR2_IIIb_IIIc.yaml"
GENE_CONFIG_DIR = "configs/genes"
GENE_DRAFT_DIR = "configs/genes/drafts"
MODULE_CAPABILITIES_PATH = "configs/framework/module_capabilities.yaml"
EVENT_DETECTORS_PATH = "configs/framework/event_detectors.yaml"
EVENT_DETECTOR_CONTRACT_PATH = "configs/framework/event_detector_contract.yaml"
CORE_ANALYSIS_CONTRACT_PATH = "configs/framework/core_gene_analysis_contract.yaml"

# Marker used by draft/undefined-event configs to state that the event region
# still needs a gene-specific detector before the gene can be configured/run.
REQUIRES_EVENT_DETECTOR = "requires_gene_specific_event_detector"

# Event-analysis lifecycle status (event.status). "configured" = a known event
# region with labels/markers (FGFR2). "user_defined" = coordinates supplied by a
# user. "not_configured"/"optional"/"disabled" = no event-specific analysis.
EVENT_STATUS_CONFIGURED = "configured"
EVENT_STATUS_USER_DEFINED = "user_defined"
EVENT_STATUS_NOT_CONFIGURED = "not_configured"
# event-region-bearing statuses (i.e. an event region actually exists)
_EVENT_PRESENT_STATUSES = {EVENT_STATUS_CONFIGURED, EVENT_STATUS_USER_DEFINED}
# analysis_modes.event_analysis vocabulary
_EVENT_MODE_VALUES = {"configured", "auto", "user_defined", "disabled", "optional", "not_configured"}

# Safe fallback UI labels (identical to the current FGFR2 wording) so any
# consumer keeps working even if a config omits ui_labels entirely.
_DEFAULT_UI_LABELS: Dict[str, str] = {
    "gene_explorer": "Gene Explorer",
    "event_region": "Cassette",
    "event_region_full": "IIIb/IIIc cassette",
    "event_discriminating_columns": "IIIb/IIIc-discriminating columns",
    "boundary_relation": "Boundary Consistency",
    "domain_relation_description": "Cassette-to-domain boundary consistency",
    "reference_comparison": "Human comparison",
}


class GeneConfigError(ValueError):
    """Raised when a gene_config file is missing required fields or is malformed."""


@dataclass
class GeneConfig:
    """A validated, gene/event-agnostic view over a gene_config.yaml."""

    raw: Dict[str, Any]
    source_path: Optional[str] = None

    # --- analysis ---------------------------------------------------------- #
    @property
    def schema_version(self) -> int:
        return int(self.raw.get("schema_version", 1) or 1)

    @property
    def analysis_id(self) -> str:
        return str(self.raw.get("analysis", {}).get("id", "")).strip()

    @property
    def status(self) -> str:
        """Lifecycle status: 'active' (runnable) or e.g. 'draft_not_runnable'."""
        return str(self.raw.get("status", "active") or "active").strip()

    @property
    def is_draft(self) -> bool:
        return self.status != "active" or str(self.status).startswith("draft")

    @property
    def is_core_only_pilot(self) -> bool:
        """A Core-only proof-of-concept: core gene analysis, no configured event."""
        return self.status == "core_only_pilot"

    @property
    def experimental(self) -> bool:
        """Experimental analyses are not offered as normal runnable choices, but may
        be exposed behind an explicit 'experimental' label (e.g. core-only pilots)."""
        r = self.raw.get("runnable")
        if isinstance(r, str) and r.strip().lower() in ("experimental", "false_or_experimental"):
            return True
        return self.is_core_only_pilot

    @property
    def runnable(self) -> bool:
        """Whether this analysis may be offered as a *normal* runnable choice.

        Only 'active' configs default to runnable. Draft / pilot / experimental
        configs stay non-runnable unless a real boolean 'runnable: true' is set.
        A string 'runnable: experimental' means experimental, i.e. NOT normally runnable.
        """
        if "runnable" in self.raw:
            r = self.raw.get("runnable")
            if isinstance(r, str):
                return False  # e.g. "experimental" / "false_or_experimental"
            return bool(r)
        return self.status == "active"

    @property
    def support_level(self) -> str:
        """Coarse support level used across docs, probe and API.

        validated_event_analysis | core_only_pilot | draft_not_runnable
        """
        if self.is_core_only_pilot:
            return "core_only_pilot"
        if self.status == "active" and self.runnable and self.has_event:
            return "validated_event_analysis"
        return "draft_not_runnable"

    @property
    def analysis_display_name(self) -> str:
        a = self.raw.get("analysis", {})
        return str(a.get("display_name") or self.analysis_id).strip()

    @property
    def analysis_description(self) -> str:
        return str(self.raw.get("analysis", {}).get("description", "") or "").strip()

    # --- gene -------------------------------------------------------------- #
    @property
    def gene_symbol(self) -> str:
        return str(self.raw.get("gene", {}).get("symbol", "")).strip()

    @property
    def gene_display_name(self) -> str:
        g = self.raw.get("gene", {})
        return str(g.get("display_name") or self.gene_symbol).strip()

    @property
    def reference_species(self) -> str:
        return str(self.raw.get("gene", {}).get("reference_species", "") or "").strip()

    # --- analysis modes ---------------------------------------------------- #
    @property
    def event_status(self) -> str:
        """Event lifecycle status. Explicit event.status wins; otherwise inferred.

        Backward compatible: FGFR2 (has a real event.type + labels) infers
        'configured'; a draft with type == requires_gene_specific_event_detector
        (or no event type) infers 'not_configured'.
        """
        e = self.raw.get("event", {}) or {}
        st = str(e.get("status", "") or "").strip().lower()
        if st:
            return st
        t = str(e.get("type", "") or "").strip()
        if not t or t == REQUIRES_EVENT_DETECTOR:
            return EVENT_STATUS_NOT_CONFIGURED
        return EVENT_STATUS_CONFIGURED

    @property
    def analysis_modes(self) -> Dict[str, Any]:
        """The core/event analysis-mode block, with backward-compatible defaults.

        core_gene_analysis defaults to True (every supported gene can run core
        analysis). event_analysis is inferred from event.status when not given.
        """
        am = self.raw.get("analysis_modes", {}) or {}
        core = am.get("core_gene_analysis", True)
        ev = am.get("event_analysis")
        if ev is None:
            st = self.event_status
            ev = ("configured" if st == EVENT_STATUS_CONFIGURED
                  else "user_defined" if st == EVENT_STATUS_USER_DEFINED
                  else "optional")
        ev = str(ev).strip().lower()
        if ev not in _EVENT_MODE_VALUES:
            ev = "optional"
        return {"core_gene_analysis": bool(core), "event_analysis": ev}

    @property
    def core_gene_analysis_enabled(self) -> bool:
        return bool(self.analysis_modes["core_gene_analysis"])

    @property
    def event_analysis_mode(self) -> str:
        return self.analysis_modes["event_analysis"]

    @property
    def has_event(self) -> bool:
        """True iff a usable event region is configured (configured/user_defined
        with a concrete event type). Core-only genes return False."""
        if self.event_analysis_mode == "disabled":
            return False
        if self.event_status not in _EVENT_PRESENT_STATUSES:
            return False
        t = self.event_type
        return bool(t) and t != REQUIRES_EVENT_DETECTOR

    @property
    def event_region_bounds(self) -> Dict[str, Any]:
        """User-defined event region coordinates, if provided (else empty)."""
        e = self.raw.get("event", {}) or {}
        return {
            "reference_protein": str(e.get("reference_protein", "") or "").strip(),
            "region_start_aa": e.get("region_start_aa"),
            "region_end_aa": e.get("region_end_aa"),
        }

    # --- event ------------------------------------------------------------- #
    @property
    def event_id(self) -> str:
        return str(self.raw.get("event", {}).get("id", "")).strip()

    @property
    def event_type(self) -> str:
        return str(self.raw.get("event", {}).get("type", "")).strip()

    @property
    def event_display_name(self) -> str:
        e = self.raw.get("event", {})
        return str(e.get("display_name") or self.event_id).strip()

    @property
    def event_generic_label(self) -> str:
        return str(self.raw.get("event", {}).get("generic_label", "event_region") or "event_region").strip()

    @property
    def event_labels(self) -> List[Dict[str, Any]]:
        labels = self.raw.get("event", {}).get("labels", []) or []
        out: List[Dict[str, Any]] = []
        for lab in labels:
            if not isinstance(lab, dict):
                continue
            out.append({
                "id": str(lab.get("id", "")).strip(),
                "display_name": str(lab.get("display_name") or lab.get("id") or "").strip(),
                "marker_sequence": str(lab.get("marker_sequence", "") or "").strip(),
                "reference_role": str(lab.get("reference_role", "") or "").strip(),
            })
        return out

    @property
    def event_label_ids(self) -> List[str]:
        return [lab["id"] for lab in self.event_labels if lab["id"]]

    # --- reference control ------------------------------------------------- #
    @property
    def reference_control(self) -> Dict[str, Any]:
        rc = self.raw.get("reference_control", {}) or {}
        return {
            "enabled": bool(rc.get("enabled", False)),
            "species": str(rc.get("species", self.reference_species) or "").strip(),
            "source": str(rc.get("source", "") or "").strip(),
            "role": str(rc.get("role", "reference_control_only") or "reference_control_only").strip(),
            "note": str(rc.get("note", "") or "").strip(),
        }

    # --- views / labels ---------------------------------------------------- #
    @property
    def views(self) -> Dict[str, bool]:
        views = self.raw.get("views", {}) or {}
        return {str(k): bool(v) for k, v in views.items()}

    @property
    def enabled_views(self) -> List[str]:
        return [k for k, v in self.views.items() if v]

    @property
    def ui_labels(self) -> Dict[str, str]:
        labels = dict(_DEFAULT_UI_LABELS)
        for k, v in (self.raw.get("ui_labels", {}) or {}).items():
            if v:
                labels[str(k)] = str(v)
        return labels

    @property
    def canonical_outputs(self) -> Dict[str, str]:
        return {str(k): str(v) for k, v in (self.raw.get("canonical_outputs", {}) or {}).items()}

    # --- serialisation ----------------------------------------------------- #

    def summary_dict(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "analysis_id": self.analysis_id,
            "analysis_display_name": self.analysis_display_name,
            "gene_symbol": self.gene_symbol,
            "gene_display_name": self.gene_display_name,
            "analysis_modes": self.analysis_modes,
            "core_gene_analysis_enabled": self.core_gene_analysis_enabled,
            "event_analysis_mode": self.event_analysis_mode,
            "event_status": self.event_status,
            "has_event": self.has_event,
            "support_level": self.support_level,
            "experimental": self.experimental,
            "event_id": self.event_id,
            "event_type": self.event_type,
            "event_display_name": self.event_display_name,
            "event_labels": self.event_labels,
            "reference_species": self.reference_species,
            "reference_control": self.reference_control,
            "enabled_views": self.enabled_views,
            "ui_labels": self.ui_labels,
        }


# --------------------------------------------------------------------------- #
# Loading + validation
# --------------------------------------------------------------------------- #
# Always required (core gene-level analysis is possible for any protein-coding
# gene). Event fields are required ONLY when the config declares a configured or
# user-defined event region (see validate_gene_config).
_REQUIRED = [
    ("analysis.id", lambda c: c.analysis_id),
    ("gene.symbol", lambda c: c.gene_symbol),
]


def _resolve_path(path: str | Path) -> Path:
    p = Path(path)
    if not p.is_absolute():
        p = (PROJECT_ROOT / p)
    return p


def load_gene_config(path: str | Path) -> GeneConfig:
    """Load and validate a gene_config.yaml. Raises GeneConfigError on problems."""
    if yaml is None:  # pragma: no cover
        raise GeneConfigError(
            f"PyYAML is required to read gene configs but could not be imported: {_YAML_IMPORT_ERROR}. "
            "Install it with: python -m pip install pyyaml")

    p = _resolve_path(path)
    if not p.is_file():
        raise GeneConfigError(f"Gene config not found: {p}")
    try:
        raw = yaml.safe_load(p.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise GeneConfigError(f"Gene config is not valid YAML ({p}): {exc}") from exc
    if not isinstance(raw, dict):
        raise GeneConfigError(f"Gene config must be a YAML mapping at the top level ({p}).")

    # Store the path relative to the project root when possible (portable).
    try:
        src = str(p.relative_to(PROJECT_ROOT))
    except ValueError:
        src = str(p)
    cfg = GeneConfig(raw=raw, source_path=src)
    validate_gene_config(cfg)
    return cfg


def validate_gene_config(cfg: GeneConfig) -> List[str]:
    """Validate required fields. Returns a list of warnings; raises on hard errors."""
    problems: List[str] = []
    for field_name, getter in _REQUIRED:
        try:
            value = getter(cfg)
        except Exception as exc:  # pragma: no cover - defensive
            problems.append(f"{field_name}: could not read ({exc})")
            continue
        if not value:
            problems.append(f"{field_name}: missing or empty (required)")
    # Event fields are only mandatory when an event region is actually declared.
    if cfg.has_event:
        if not cfg.event_id:
            problems.append("event.id: missing or empty (required when an event is configured)")
        if not cfg.event_type:
            problems.append("event.type: missing or empty (required when an event is configured)")

    if problems:
        raise GeneConfigError(
            "Invalid gene config:\n  - " + "\n  - ".join(problems))

    warnings: List[str] = []
    if not cfg.core_gene_analysis_enabled:
        warnings.append("analysis_modes.core_gene_analysis is false; no analysis mode is enabled.")
    if cfg.has_event:
        if len(cfg.event_labels) < 1:
            warnings.append("event.labels is empty; event-region views will have no labels.")
        for i, lab in enumerate(cfg.event_labels):
            if not lab["id"]:
                warnings.append(f"event.labels[{i}].id is empty.")
        if cfg.event_status == EVENT_STATUS_USER_DEFINED:
            b = cfg.event_region_bounds
            if b["region_start_aa"] is None or b["region_end_aa"] is None:
                warnings.append("event.status is user_defined but region_start_aa/region_end_aa are not both set.")
    else:
        warnings.append(
            "No event region is configured (event_analysis="
            f"{cfg.event_analysis_mode}); only core gene-level views will be available.")
    if not cfg.reference_species:
        warnings.append("gene.reference_species is empty; reference comparison may be disabled.")
    if not cfg.views:
        warnings.append("views section is empty; all views default to disabled.")
    return warnings


def default_gene_config() -> GeneConfig:
    """Load the canonical FGFR2 config; the backward-compatible default."""
    return load_gene_config(DEFAULT_GENE_CONFIG)


def load_gene_config_lenient(path: str | Path) -> GeneConfig:
    """Load a config WITHOUT hard validation (for drafts with undefined events)."""
    if yaml is None:  # pragma: no cover
        raise GeneConfigError("PyYAML is required to read gene configs.")
    p = _resolve_path(path)
    if not p.is_file():
        raise GeneConfigError(f"Gene config not found: {p}")
    raw = yaml.safe_load(p.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise GeneConfigError(f"Gene config must be a YAML mapping ({p}).")
    try:
        src = str(p.relative_to(PROJECT_ROOT))
    except ValueError:
        src = str(p)
    return GeneConfig(raw=raw, source_path=src)


def build_generic_gene_config(
    gene_symbol: str,
    *,
    generated_by: str = "dynamic_run_creation",
    reference_species: str = "",
    extra_provenance: Optional[Dict[str, Any]] = None,
) -> GeneConfig:
    """Synthesize a generic *core-only* gene config from just a gene symbol.

    This is what makes the framework truly generic: a user can pick ANY protein-
    coding gene symbol and we build a valid, no-event core-only analysis config
    in memory — no pre-existing ``configs/genes/**`` YAML is required. FGFR2 keeps
    its validated, hand-authored specialization; every other gene routes here.

    The resulting config has ``source_path=None`` (it did not come from a file);
    callers that want a run-local copy should serialize it with
    :func:`gene_config_to_yaml`.
    """
    sym = str(gene_symbol or "").strip().upper()
    if not sym:
        raise GeneConfigError("A gene symbol is required to build a generic gene config.")

    provenance: Dict[str, Any] = {
        "generated": True,
        "generated_by": generated_by,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": "dynamic_generic_run_creation",
        "note": (
            "Auto-generated generic core-only config. No gene-specific YAML was "
            "required or edited. Only FGFR2 uses a hand-authored validated config."
        ),
    }
    if extra_provenance:
        provenance.update(extra_provenance)

    raw: Dict[str, Any] = {
        "schema_version": 1,
        # Lifecycle: an experimental, exploratory core-only analysis (no event).
        "status": "core_only_pilot",
        "runnable": "experimental",
        "analysis": {
            "id": f"{sym}_core_only_pilot",
            "display_name": f"{sym} core-only exploratory analysis",
            "description": (
                f"Generic exploratory analysis for {sym}: gene-level comparative "
                "exon–protein mapping and post-InterPro domain/boundary analysis. "
                "No validated event region is configured."
            ),
        },
        "gene": {
            "symbol": sym,
            "display_name": sym,
            "reference_species": str(reference_species or "").strip(),
        },
        # Core gene analysis on; no event analysis (optional/undefined).
        "analysis_modes": {"core_gene_analysis": True, "event_analysis": "optional"},
        "event": {
            "status": EVENT_STATUS_NOT_CONFIGURED,
            "type": REQUIRES_EVENT_DETECTOR,
            "labels": [],
        },
        "reference_control": {"enabled": False},
        # Only core/shared views; event-specific views stay disabled.
        "views": {},
        "ui_labels": {},
        "canonical_outputs": {},
        "provenance": provenance,
    }
    return GeneConfig(raw=raw, source_path=None)


def gene_config_to_yaml(cfg: GeneConfig) -> str:
    """Serialize a GeneConfig's raw mapping to YAML text (with a provenance header).

    Used to drop a run-local ``gene_config.yaml`` for generated configs so a run
    is fully self-describing and reproducible without any repo-level YAML.
    """
    if yaml is None:  # pragma: no cover
        raise GeneConfigError("PyYAML is required to serialize gene configs.")
    prov = cfg.raw.get("provenance", {}) if isinstance(cfg.raw, dict) else {}
    header = [
        "# Auto-generated generic gene config — no manual editing required.",
        f"# gene: {cfg.gene_symbol}",
    ]
    if prov.get("generated_at"):
        header.append(f"# generated_at: {prov.get('generated_at')}")
    if prov.get("generated_by"):
        header.append(f"# generated_by: {prov.get('generated_by')}")
    header.append("")
    body = yaml.safe_dump(cfg.raw, sort_keys=False, default_flow_style=False, allow_unicode=True)
    return "\n".join(header) + body


def load_module_capabilities() -> Dict[str, Any]:
    """Load configs/framework/module_capabilities.yaml (best-effort, {} on failure)."""
    if yaml is None:
        return {}
    p = _resolve_path(MODULE_CAPABILITIES_PATH)
    if not p.is_file():
        return {}
    try:
        data = yaml.safe_load(p.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def capability_summary() -> Dict[str, Any]:
    """Reusability roll-up derived from the module capability inventory."""
    caps = load_module_capabilities()
    modules = caps.get("modules", {}) if isinstance(caps, dict) else {}
    reusable, gene_specific, needs_work = [], [], []
    for name, info in modules.items():
        if not isinstance(info, dict):
            continue
        flag = info.get("reusable_for_future_genes")
        scope = info.get("scope", "")
        if flag is True:
            reusable.append(name)
        elif flag is False:
            gene_specific.append(name)
        else:  # "partial" or unknown
            needs_work.append(name)
        _ = scope
    return {
        "reusable_modules": sorted(reusable),
        "gene_specific_modules": sorted(gene_specific),
        "partial_or_event_specific_modules": sorted(needs_work),
        "modules": modules,
    }


def load_event_detectors() -> Dict[str, Any]:
    """Load configs/framework/event_detectors.yaml (best-effort, {} on failure)."""
    if yaml is None:
        return {}
    p = _resolve_path(EVENT_DETECTORS_PATH)
    if not p.is_file():
        return {}
    try:
        data = yaml.safe_load(p.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def supported_detectors() -> Dict[str, Dict[str, Any]]:
    """Detector name -> spec, for detectors with status == 'supported'."""
    dets = (load_event_detectors().get("detectors", {}) or {})
    return {k: v for k, v in dets.items()
            if isinstance(v, dict) and str(v.get("status", "")).lower() == "supported"}


def planned_detectors() -> Dict[str, Dict[str, Any]]:
    """Detector name -> spec, for detectors that are not yet supported."""
    dets = (load_event_detectors().get("detectors", {}) or {})
    return {k: v for k, v in dets.items()
            if isinstance(v, dict) and str(v.get("status", "")).lower() != "supported"}


def detector_for_analysis(analysis_id: str) -> Optional[Dict[str, Any]]:
    """Return the detector spec (with its registry name) for an analysis_id, or None.

    Prefers a supported detector; only returns a supported one so that callers can
    treat "has detector" as "is runnable".
    """
    if not analysis_id:
        return None
    for name, spec in supported_detectors().items():
        if str(spec.get("analysis_id", "")) == analysis_id:
            return {"detector": name, **spec}
    return None


def load_event_detector_contract() -> Dict[str, Any]:
    """Load the machine-readable detector contract (best-effort, {} on failure)."""
    if yaml is None:
        return {}
    p = _resolve_path(EVENT_DETECTOR_CONTRACT_PATH)
    if not p.is_file():
        return {}
    try:
        data = yaml.safe_load(p.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def load_core_analysis_contract() -> Dict[str, Any]:
    """Load the machine-readable core gene-analysis contract (best-effort)."""
    if yaml is None:
        return {}
    p = _resolve_path(CORE_ANALYSIS_CONTRACT_PATH)
    if not p.is_file():
        return {}
    try:
        data = yaml.safe_load(p.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def discover_analyses() -> Dict[str, Any]:
    """List configured analyses: supported (active/runnable) vs drafts."""
    supported: List[Dict[str, Any]] = []
    drafts: List[Dict[str, Any]] = []

    pilots: List[Dict[str, Any]] = []

    def _summarise(cfg: GeneConfig) -> Dict[str, Any]:
        det = detector_for_analysis(cfg.analysis_id)
        return {
            "analysis_id": cfg.analysis_id,
            "gene_symbol": cfg.gene_symbol,
            "event_id": cfg.event_id,
            "event_type": cfg.event_type,
            "status": cfg.status,
            "runnable": cfg.runnable,
            "experimental": cfg.experimental,
            "support_level": cfg.support_level,
            "core_gene_analysis_enabled": cfg.core_gene_analysis_enabled,
            "has_event": cfg.has_event,
            "event_analysis_mode": cfg.event_analysis_mode,
            "has_supported_detector": det is not None,
            "detector": (det or {}).get("detector"),
            "config": cfg.source_path,
        }

    def _is_runnable(summary: Dict[str, Any]) -> bool:
        # A config must be active/runnable and enable core gene analysis.
        if not (summary["runnable"] and summary["core_gene_analysis_enabled"]):
            return False
        # If it declares a configured event region, a supported detector is
        # required. Core-only genes are runnable on core analysis alone.
        if summary["has_event"]:
            return bool(summary["has_supported_detector"])
        return True

    gene_dir = _resolve_path(GENE_CONFIG_DIR)
    if gene_dir.is_dir():
        for p in sorted(gene_dir.glob("*.yaml")):
            if p.name.startswith("TEMPLATE"):
                continue
            try:
                cfg = load_gene_config_lenient(p)
            except GeneConfigError:
                continue
            s = _summarise(cfg)
            (supported if _is_runnable(s) else drafts).append(s)

    draft_dir = _resolve_path(GENE_DRAFT_DIR)
    if draft_dir.is_dir():
        for p in sorted(draft_dir.glob("*.yaml")):
            if p.name.startswith("TEMPLATE"):
                continue
            try:
                cfg = load_gene_config_lenient(p)
            except GeneConfigError:
                continue
            s = _summarise(cfg)
            (pilots if cfg.is_core_only_pilot else drafts).append(s)

    return {
        "supported": [a["analysis_id"] for a in supported],
        "core_only_pilots": [a["analysis_id"] for a in pilots],
        "drafts": [a["analysis_id"] for a in drafts],
        "supported_detail": supported,
        "core_only_pilots_detail": pilots,
        "drafts_detail": drafts,
        "note": "Support levels: validated_event_analysis (runnable, FGFR2), "
                "core_only_pilot (experimental core-only, no event region), "
                "draft_not_runnable (config only). An analysis is normally runnable "
                "only with an active config AND a supported event detector.",
    }


def resolve_run_config_path(run_config: Dict[str, Any], run_dir: Optional[Path] = None) -> str:
    """Determine the gene_config path for a run, defaulting to FGFR2 for old runs.

    Order of preference:
      1. runs/<run_id>/gene_config.yaml (run-local copy), if present
      2. run_config.json["gene_config"]
      3. DEFAULT_GENE_CONFIG (FGFR2 IIIb/IIIc)
    """
    if run_dir is not None:
        local = Path(run_dir) / "gene_config.yaml"
        if local.is_file():
            try:
                return str(local.relative_to(PROJECT_ROOT))
            except ValueError:
                return str(local)
    cfg_path = str((run_config or {}).get("gene_config") or "").strip()
    if cfg_path and _resolve_path(cfg_path).is_file():
        return cfg_path
    return DEFAULT_GENE_CONFIG


def resolve_run_analysis(run_config: Dict[str, Any], run_dir: Optional[Path] = None) -> GeneConfig:
    """Return the GeneConfig for a run, defaulting to FGFR2 for backward compat."""
    try:
        return load_gene_config(resolve_run_config_path(run_config, run_dir))
    except GeneConfigError:
        # Never fail a read because of a bad/missing config: fall back to FGFR2.
        return default_gene_config()


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Load / validate a gene_config.yaml.")
    ap.add_argument("--config", default=DEFAULT_GENE_CONFIG,
                    help="Path to the gene config YAML (default: FGFR2 IIIb/IIIc).")
    ap.add_argument("--validate", action="store_true",
                    help="Validate required fields and print a summary.")
    ap.add_argument("--json", action="store_true",
                    help="Print the resolved config summary as JSON.")
    args = ap.parse_args(argv)

    try:
        cfg = load_gene_config(args.config)
        warnings = validate_gene_config(cfg)
    except GeneConfigError as exc:
        print(f"INVALID: {exc}", file=sys.stderr)
        return 2

    if args.json:
        import json
        print(json.dumps(cfg.summary_dict(), indent=2))
        return 0

    print(f"OK: {args.config}")
    print(f"  analysis_id      : {cfg.analysis_id}")
    print(f"  gene_symbol      : {cfg.gene_symbol}")
    print(f"  core_analysis    : {cfg.core_gene_analysis_enabled}")
    print(f"  event_analysis   : {cfg.event_analysis_mode} (status={cfg.event_status}, has_event={cfg.has_event})")
    print(f"  event_id         : {cfg.event_id or '(none)'}")
    print(f"  event_type       : {cfg.event_type or '(none)'}")
    print(f"  event_labels     : {', '.join(cfg.event_label_ids) or '(none)'}")
    print(f"  reference_species: {cfg.reference_species or '(none)'}")
    print(f"  enabled_views    : {', '.join(cfg.enabled_views) or '(none)'}")
    if warnings:
        print("  warnings:")
        for w in warnings:
            print(f"    - {w}")
    else:
        print("  warnings         : none")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
