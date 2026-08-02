#!/usr/bin/env python3
"""
build_website_indices.py

Phase-1 website data generator for ExonDomainCompare.

Reads a *finished* FGFR2 pre-InterPro closure run directory (the
``13_final_pre_interpro_closure`` folder, which carries the truth table, the
evidence stack, the freeze, figures, tables, reports, the consistency gate and
the run/step metadata) and emits six compact JSON indices that the web frontend
consumes:

    run_index.json        run id, mode, gate, KPI counts, step timeline
    species_index.json    per-species IIIb/IIIc final status (truth table)
    evidence_stack.json   normalized per-row evidence layers (heatmap source)
    figure_index.json     curated, grouped figure catalogue with captions
    download_index.json    key downloadable artefacts
    freeze_index.json     reproducibility / InterPro-ready export view

Hard rules honoured:
  * FINAL biological status comes ONLY from final_pre_interpro_truth_table.tsv
    (and the closure evidence stack derived from it). Provenance/legacy columns
    are surfaced as detail, never as decision logic.
  * No invented data: every value is read from the run folder; missing inputs
    degrade gracefully (empty/absent rather than fabricated).
  * Rescued-and-validated rows are reported as accepted; only review/excluded
    rows are flagged.

The indices use project-root-relative POSIX file paths so the backend can serve
them through its sandboxed /api/download endpoint.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from fgfr2 import human_reference_control
from shared_gene_analysis import run_availability as ra
from shared_gene_analysis import species_order as _species_order
from shared_gene_analysis import synteny_contract as sc
from shared_gene_analysis.public_paths import (
    sanitize_public_payload,
    write_public_download_projections,
)
from framework.data_contract import write_freshness_contract, write_payload_contracts

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# --------------------------------------------------------------------------- #
# small IO helpers
# --------------------------------------------------------------------------- #
def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path or not Path(path).exists():
        return []
    with Path(path).open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            return []
        return [dict(r) for r in reader]


def read_json(path: Path, fallback: Any = None) -> Any:
    if not path or not Path(path).exists():
        return fallback
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return fallback


def rel(path: Path) -> str:
    """project-root-relative POSIX path (for the /api/download endpoint)."""
    try:
        return Path(path).resolve().relative_to(PROJECT_ROOT).as_posix()
    except Exception:
        return Path(path).as_posix()


def file_meta(path: Path) -> Optional[Dict[str, Any]]:
    p = Path(path)
    if not p.exists() or not p.is_file():
        return None
    size = p.stat().st_size
    return {"name": p.name, "path": rel(p), "size_bytes": size, "size_human": human_size(size)}


def _availability_block(label: str, source: Path, has_payload: bool,
                        extra_sources: Sequence[Path] = ()) -> Dict[str, Any]:
    """Why this index is or is not usable, in the index itself.

    A bare ``available: false`` forces the reader to guess between "this species has no
    such result" and "this file was never written". The frontend guessed biology and told
    users a cassette was absent when the cassette table had merely not been rebuilt. The
    distinction is knowable here — the source table is either missing, present but empty,
    or present with rows — so it is recorded here.
    """
    sources = [Path(source), *[Path(p) for p in extra_sources]]
    missing = [rel(p) for p in sources if not p.is_file()]
    if has_payload:
        return {"state": ra.AVAILABLE, "reason": f"{label.capitalize()} is available.",
                "source_tables": [rel(p) for p in sources], "missing_inputs": []}
    if missing:
        return {"state": ra.TECHNICALLY_MISSING,
                "reason": (f"Expected {label} outputs were not generated. "
                           "Retry local analysis."),
                "source_tables": [rel(p) for p in sources],
                "missing_inputs": missing}
    return {"state": ra.SCIENTIFICALLY_UNAVAILABLE,
            "reason": (f"The {label} source table was written but contains no supported "
                       "row for this run's models."),
            "source_tables": [rel(p) for p in sources], "missing_inputs": []}


def human_size(n: int) -> str:
    f = float(n)
    for unit in ("B", "KB", "MB", "GB"):
        if f < 1024 or unit == "GB":
            return f"{f:.0f} {unit}" if unit == "B" else f"{f:.1f} {unit}"
        f /= 1024
    return f"{n} B"


def count_fasta(path: Path) -> int:
    p = Path(path)
    if not p.exists():
        return 0
    n = 0
    with p.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def to_int(v: Any) -> Optional[int]:
    try:
        if v in (None, "") or str(v).lower() == "nan":
            return None
        return int(float(str(v)))
    except Exception:
        return None


# --------------------------------------------------------------------------- #
# status semantics (single place so frontend + backend agree)
# --------------------------------------------------------------------------- #
# Normalized evidence-stack cell tokens -> colour class.
TOKEN_CLASS = {
    "pass": "accepted", "strong": "accepted", "robust": "accepted",
    "rescued_ok": "accepted", "high": "accepted", "true_ok": "accepted",
    "minor": "minor", "moderate": "minor", "supported": "minor",
    "review": "review", "supplement": "review", "manual_review": "review",
    "unresolved": "review", "warn": "minor",
    "fail": "excluded", "excluded": "excluded", "missing": "unknown", "": "unknown",
}


def token_class(token: str) -> str:
    t = (token or "").strip().lower()
    return TOKEN_CLASS.get(t, "neutral")


def readiness_class(value: str) -> str:
    """Row-level readiness -> {accepted, minor, review, excluded}."""
    v = (value or "").lower()
    if not v:
        return "unknown"
    if "excluded" in v:
        return "excluded"
    if "supplement" in v or "review" in v or "unresolved" in v:
        return "review"
    if "with_minor_flags" in v or "minor" in v:
        return "minor"
    if "ready" in v or "primary" in v or "supported" in v:
        return "accepted"
    return "neutral"


def claim_class(value: str) -> str:
    v = (value or "").lower()
    if "excluded" in v:
        return "excluded"
    if "supplement" in v or "review" in v:
        return "review"
    if "minor" in v:
        return "minor"
    if "primary" in v or "supported" in v:
        return "accepted"
    return "neutral"


# --------------------------------------------------------------------------- #
# Final-truth-driven UI status semantics
# --------------------------------------------------------------------------- #
# Every evidence layer is expressed as a small UI cell:
#   {"class": colour-family, "tone": optional refinement, "value": short label,
#    "note": one-sentence explanation}
# Colour family is one of {accepted, minor, review, excluded, neutral}; "minor"
# is an ACCEPTED-family colour (accepted with minor flags), never a warning.
# The tone drives a small provenance/offset icon (corrected / rescued / offset).
#
# HARD RULE: a native/reference coordinate offset or an upstream-label swap is a
# framework CORRECTION, not biological uncertainty. Only genuine hard failures
# (AA1/unsafe mapping, missing protein link) or a final review/supplement row
# may downgrade a layer to review/excluded.

# Known status tokens -> (colour family, short label, one-sentence note).
SHORT_STATUS: Dict[str, tuple] = {
    "fgfr2_ortholog_high_confidence": ("accepted", "high-confidence ortholog",
                                       "FGFR2 orthology confidently established."),
    "fgfr2_ortholog_supported_with_warnings": ("minor", "ortholog · minor warnings",
                                                "FGFR2 orthology supported with minor warnings."),
    "full_length_msa_pass": ("accepted", "conserved",
                             "Full-length alignment passes quality control."),
    "full_length_msa_pass_with_minor_flags": ("minor", "conserved · minor flags",
                                              "Full-length alignment passes with minor flags."),
    "robust_boundary": ("accepted", "robust boundary",
                        "Cassette boundary is robust across the alignment."),
    "supported_boundary_with_minor_flags": ("minor", "supported boundary",
                                            "Cassette boundary supported with minor flags."),
    "synteny_strong": ("accepted", "strong synteny",
                       "Local gene order strongly conserved around FGFR2."),
    "synteny_supported_with_minor_rearrangement": ("minor", "synteny · minor rearrangement",
                                                   "Local synteny supported with a minor rearrangement."),
    "synteny_partial_blast_supported": ("minor", "partial synteny (BLAST)",
                                        "Synteny partially supported via BLAST homology."),
    "protein_integrity_pass": ("accepted", "intact",
                               "Protein length / integrity within expected range."),
    "protein_length_outlier_review": ("minor", "length outlier",
                                       "Protein length is an outlier (annotation note, final row accepted)."),
    "pass": ("accepted", "pass", ""),
}


def _humanize(token: str) -> str:
    return (token or "").replace("_", " ").strip() or "—"


def status_from_token(raw: str) -> Dict[str, Any]:
    """Generic evidence layer (orthology / MSA / synteny / integrity)."""
    t = (raw or "").strip().lower()
    if not t:
        return {"class": "neutral", "tone": "", "value": "—", "note": ""}
    if t in SHORT_STATUS:
        cls, val, note = SHORT_STATUS[t]
        return {"class": cls, "tone": "", "value": val, "note": note}
    if "fail" in t or "excluded" in t:
        return {"class": "excluded", "tone": "", "value": _humanize(raw), "note": ""}
    if "review" in t or "unresolved" in t or "supplement" in t:
        return {"class": "review", "tone": "", "value": _humanize(raw), "note": ""}
    if any(k in t for k in ("minor", "warning", "partial", "rearrangement", "supported_with", "with_flags")):
        return {"class": "minor", "tone": "", "value": _humanize(raw), "note": ""}
    if any(k in t for k in ("pass", "strong", "robust", "high", "confirmed", "ortholog", "supported")):
        return {"class": "accepted", "tone": "", "value": _humanize(raw), "note": ""}
    return {"class": statusify(raw), "tone": "", "value": _humanize(raw), "note": ""}


def label_assignment_status(consistency_raw: str, readiness_cls: str) -> Dict[str, Any]:
    """Final isoform assignment / annotation reconciliation (NOT upstream label)."""
    v = (consistency_raw or "").lower()
    if "swapped" in v:
        return {"class": "accepted", "tone": "corrected", "value": "corrected",
                "note": "Upstream annotation was corrected by sequence evidence."}
    if "consistent" in v or "agree" in v:
        return {"class": "accepted", "tone": "", "value": "consistent",
                "note": "Source and sequence-calibrated labels agree."}
    if "ambiguous" in v:
        cls = "review" if readiness_cls in ("review", "excluded") else "minor"
        return {"class": cls, "tone": "", "value": "ambiguous label",
                "note": "Upstream label was ambiguous; final label set from sequence evidence."}
    if "unresolved" in v:
        cls = "review" if readiness_cls in ("review", "excluded") else "minor"
        return {"class": cls, "tone": "", "value": "no upstream sequence",
                "note": "No upstream sequence to reconcile; final label taken from framework evidence."}
    return {"class": statusify(consistency_raw), "tone": "", "value": _humanize(consistency_raw), "note": ""}


def coordinate_status(coord_raw: str, readiness_cls: str, integrity_raw: str = "") -> Dict[str, Any]:
    """Coordinate mapping: a native/reference offset is a note, not a failure."""
    v = (coord_raw or "").lower()
    # hard mapping failures dominate and are NEVER weakened
    if any(k in v for k in ("aa1", "unsafe", "mapping_bug", "no_protein",
                            "no_native_mapping", "unmapped", "missing")):
        return {"class": "excluded", "tone": "", "value": "mapping failed",
                "note": "Coordinate mapping failed a hard integrity gate."}
    if "same_native" in v or "validated" in v or "exact" in v or v == "pass":
        return {"class": "accepted", "tone": "", "value": "validated",
                "note": "Native reference coordinate validated."}
    if "offset" in v:
        if readiness_cls in ("review", "excluded"):
            return {"class": "review", "tone": "offset", "value": "offset · review row",
                    "note": "Native coordinate offset on a review row."}
        sev = "major" if "major" in v else "moderate" if "moderate" in v else "native"
        return {"class": "minor", "tone": "offset", "value": f"validated · {sev} offset",
                "note": (f"Final coordinate validated; {sev} native offset versus the upstream "
                         "annotation is a note, not a biological failure.")}
    if "fail" in v:
        return {"class": "excluded", "tone": "", "value": "fail",
                "note": "Coordinate validation failed."}
    if "review" in v:
        cls = "review" if readiness_cls in ("review", "excluded") else "minor"
        return {"class": cls, "tone": "", "value": "flagged",
                "note": ("Coordinate flagged; final row accepted." if cls == "minor"
                         else "Coordinate flagged for review.")}
    return {"class": statusify(coord_raw), "tone": "", "value": _humanize(coord_raw), "note": ""}


def rescue_status_ui(rescue_raw: str, readiness_cls: str) -> Dict[str, Any]:
    v = (rescue_raw or "").lower()
    if "rescued" in v:
        return {"class": "accepted", "tone": "rescued", "value": "rescued & validated",
                "note": "Rescued with a source-compatible validated candidate; accepted with provenance retained."}
    if v in ("pass", "") or "no_rescue" in v or "not_required" in v or "confirmed" in v:
        return {"class": "accepted", "tone": "", "value": "no rescue needed",
                "note": "Current candidate confirmed; no external rescue required."}
    if readiness_cls in ("review", "excluded"):
        return {"class": "review", "tone": "", "value": _humanize(rescue_raw),
                "note": "Kept as review / supplement with provenance."}
    return {"class": statusify(rescue_raw), "tone": "", "value": _humanize(rescue_raw), "note": ""}


def readiness_status_ui(readiness_raw: str) -> Dict[str, Any]:
    v = (readiness_raw or "").lower()
    if "excluded" in v:
        return {"class": "excluded", "tone": "", "value": "excluded", "included": "excluded",
                "note": "Excluded from the InterPro input."}
    if "supplement" in v or "review" in v:
        return {"class": "review", "tone": "", "value": "supplement / review", "included": "supplement",
                "note": "Kept in the review-included FASTA (optional/supplementary); not in the primary set."}
    if "with_minor_flags" in v or "minor" in v:
        return {"class": "minor", "tone": "", "value": "primary · minor flags", "included": "primary",
                "note": "Accepted into the primary InterPro FASTA (with minor annotation flags)."}
    if "ready" in v or "primary" in v:
        return {"class": "accepted", "tone": "", "value": "primary", "included": "primary",
                "note": "Accepted into the primary InterPro FASTA."}
    return {"class": "neutral", "tone": "", "value": _humanize(readiness_raw), "included": "", "note": ""}


# User-facing evidence layers (scientific decision layers, not internal enums).
EVIDENCE_LAYERS = [
    ("label", "Final isoform assignment"),
    ("rescue", "Rescue / provenance"),
    ("orthology", "Orthology / paralog"),
    ("coordinates", "Coordinate mapping"),
    ("full_msa", "Full-length MSA"),
    ("cassette_msa", "Cassette MSA / residue agreement"),
    ("synteny", "Synteny / locus"),
    ("readiness", "InterPro input"),
]


def _evidence_sources(run_dir: Path) -> Dict[str, str]:
    return {
        "truth": rel(run_dir / "final_pre_interpro_truth_table.tsv"),
        "review": rel(run_dir / "tables" / "final_review_case_explanation.tsv"),
        "coord": rel(run_dir / "tables" / "figure3C_exon_to_protein_cassette_coordinate_map.tsv"),
        "msa": rel(run_dir / "MSA" / "final_full_length_msa_conservation_summary.tsv"),
        "cassette": rel(run_dir / "tables" / "figure6B_species_resolved_IIIb_IIIc_cassette_residue_map.tsv"),
        "synteny": rel(run_dir / "final_pre_interpro_truth_table.tsv"),
    }


def compute_layers(r: Dict[str, str], ref: Dict[str, str], sources: Dict[str, str]) -> Dict[str, Any]:
    """Single source of truth for the per-row evidence layers (stack + summary).

    ``r``   = truth-table row (authoritative final columns + provenance)
    ``ref`` = matching framework-evidence-stack row (carries the *_raw tokens and
              the numeric reference_agreement); may be empty.
    """
    rcls = readiness_class(r.get("pre_interpro_readiness_class", ""))

    def raw(name: str) -> str:
        return ref.get(f"{name}_raw") or ref.get(name) or r.get(name, "")

    label = label_assignment_status(r.get("label_consistency_status", ""), rcls)
    rescue = rescue_status_ui(raw("rescue_status") or r.get("rescue_decision", ""), rcls)
    orth = status_from_token(raw("orthology_status"))
    coord = coordinate_status(r.get("coordinate_validation_status", ""), rcls,
                              r.get("protein_integrity_status", ""))
    fmsa = status_from_token(raw("MSA_full_length_status"))
    boundary_raw = raw("boundary_robustness_class")
    cass = dict(status_from_token(boundary_raw))
    ref_ag = to_float(ref.get("reference_agreement_raw") or ref.get("reference_agreement"))
    if ref_ag is not None:
        pct = round(ref_ag * 100)
        cass["value"] = f"{cass['value']} · {pct}% ref"
        cass["note"] = (cass["note"] + f" Human-referenced residue agreement {pct}%.").strip()
        if pct < 60 and cass["class"] == "accepted":
            cass["class"] = "minor"
    syn = status_from_token(raw("synteny_validation_class"))
    ready = readiness_status_ui(r.get("pre_interpro_readiness_class", ""))

    spec = {
        "label": (label, r.get("label_consistency_status", ""), sources["truth"]),
        "rescue": (rescue, raw("rescue_status") or r.get("rescue_decision", ""), sources["review"]),
        "orthology": (orth, raw("orthology_status"), sources["truth"]),
        "coordinates": (coord, r.get("coordinate_validation_status", ""), sources["coord"]),
        "full_msa": (fmsa, raw("MSA_full_length_status"), sources["msa"]),
        "cassette_msa": (cass, boundary_raw, sources["cassette"]),
        "synteny": (syn, raw("synteny_validation_class"), sources["synteny"]),
        "readiness": (ready, r.get("pre_interpro_readiness_class", ""), sources["truth"]),
    }
    out: Dict[str, Any] = {}
    for key, (st, rawval, src) in spec.items():
        out[key] = {
            "class": st["class"], "tone": st.get("tone", ""), "value": st["value"],
            "note": st.get("note", ""), "raw": rawval, "source_table": src,
        }
    return out


# --- Phase-2 residue / synteny semantics ----------------------------------- #
# Neutral residue-agreement colour classes (no functional-effect claims).
def residue_class(value: str) -> str:
    """agreement_class / substitution_class -> {identical, conservative, nonconservative, gap}."""
    v = (value or "").strip().lower()
    if not v:
        return "gap"
    if "gap" in v or "missing" in v:
        return "gap"
    if "identical" in v:
        return "identical"
    if "nonconservative" in v or "non_conservative" in v or "non-conservative" in v:
        return "nonconservative"
    if "conservative" in v or "semi" in v:  # conservative + semi_conservative
        return "conservative"
    return "nonconservative"


# orthology-resolution method -> styling class for the synteny viewer.
SYNTENY_METHOD_CLASS = {
    "anchor": "anchor",
    "exact_symbol_match": "exact",
    "curated_ortholog": "curated",
    "broad_proteome_best_hit": "rbh",
    "broad_proteome_weak_best_hit": "weak",
    "raw_annotation_only": "unresolved",
    "": "unresolved",
}


def synteny_method_class(method: str) -> str:
    return SYNTENY_METHOD_CLASS.get((method or "").strip().lower(), "unresolved")


def to_float(v: Any) -> Optional[float]:
    try:
        if v in (None, "") or str(v).lower() == "nan":
            return None
        return round(float(str(v)), 4)
    except Exception:
        return None


def truthy(v: Any) -> bool:
    return str(v).strip().lower() in ("true", "1", "yes")


def _synteny_source_dir(run_dir: Path) -> Optional[Path]:
    """Locate the sibling step-12 MSA/synteny directory (synteny tables live there)."""
    parent = Path(run_dir).parent
    for cand in sorted(parent.glob("12_*")):
        if (cand / "synteny").exists() or (cand / "tables").exists():
            return cand
    return None


def _review_lookup(run_dir: Path) -> Dict[tuple, str]:
    """(species, isoform) -> readiness colour class, the single biological truth."""
    out: Dict[tuple, str] = {}
    for r in _truth(run_dir):
        out[(r.get("species", ""), r.get("isoform", ""))] = \
            readiness_class(r.get("pre_interpro_readiness_class", ""))
    return out


def _parse_alignment(path: Path) -> Dict[str, Any]:
    """Read a wrapped FASTA alignment (>species|isoform|protein|tag) into rows."""
    p = Path(path)
    if not p.exists():
        return {}
    rows: List[Dict[str, Any]] = []
    header: Optional[str] = None
    seq_parts: List[str] = []

    def flush():
        if header is None:
            return
        parts = header.split("|")
        species = parts[0] if parts else header
        rows.append({
            "species": species,
            "isoform": parts[1] if len(parts) > 1 else "",
            "protein_id": parts[2] if len(parts) > 2 else "",
            "tag": parts[3] if len(parts) > 3 else "",
            "seq": "".join(seq_parts),
        })

    with p.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                seq_parts = []
            elif header is not None:
                seq_parts.append(line.strip())
    flush()
    n_cols = max((len(r["seq"]) for r in rows), default=0)
    return {"rows": rows, "n_columns": n_cols}


# --------------------------------------------------------------------------- #
# index builders
# --------------------------------------------------------------------------- #
def _truth(run_dir: Path) -> List[Dict[str, str]]:
    return read_tsv(run_dir / "final_pre_interpro_truth_table.tsv")


def build_run_index(run_dir: Path) -> Dict[str, Any]:
    run_mode = read_json(run_dir / "final_pre_interpro_run_mode.json", {}) or {}
    gate = read_json(run_dir / "gates" / "final_pre_interpro_cross_table_consistency_gate.json", {}) or {}
    gate_checks = gate.get("checks", []) if isinstance(gate, dict) else []
    gate_fail = [c for c in gate_checks if str(c.get("status", "")).lower() not in ("pass", "ok", "")]

    steps_raw = read_tsv(run_dir / "final_pre_interpro_step_status.tsv")
    # Friendly phase labels for the canonical A-step ids (never fabricate granularity).
    phase_label = {
        "A0": "Resolve species list",
        "A1": "Collect models · select · classify · coordinates (steps 1–11)",
        "A2": "MSA · reconcile · rescue · synteny",
        "A2b": "Publication figures",
        "A3": "Synteny paper figures",
        "A4": "Framework evidence figure",
        "A5": "Freeze & closure",
        "A6": "Run-mode summary",
    }
    steps = []
    for s in steps_raw:
        sid = s.get("step_id", "")
        steps.append({
            "step_id": sid,
            "step_name": s.get("step_name", ""),
            "label": phase_label.get(sid, (s.get("step_name", "") or sid).replace("_", " ")),
            "status": s.get("status", ""),
            "return_code": to_int(s.get("return_code")),
            "runtime_seconds": to_int(s.get("runtime_seconds")),
            "command": s.get("command", ""),
            "output_files": s.get("output_files", ""),
            "warning_summary": s.get("warning_summary", ""),
        })

    rows = _truth(run_dir)
    species = sorted({r.get("species", "") for r in rows if r.get("species")})
    primary = [r for r in rows if readiness_class(r.get("pre_interpro_readiness_class", "")) in ("accepted", "minor")]
    review = [r for r in rows if readiness_class(r.get("pre_interpro_readiness_class", "")) == "review"]
    excluded = [r for r in rows if readiness_class(r.get("pre_interpro_readiness_class", "")) == "excluded"]
    rescued = [r for r in rows if "rescued" in (r.get("rescue_decision", "") or "").lower()
               or "rescued" in (r.get("final_label_source", "") or "").lower()]

    primary_fasta = run_dir / "freeze" / "final_pre_interpro_proteins_primary.faa"
    review_fasta = run_dir / "freeze" / "final_pre_interpro_proteins_all_review_included.faa"

    return {
        "run_id": run_mode.get("run_id") or "unknown",
        "case_study": "FGFR2 IIIb/IIIc",
        "generated_at": now_iso(),
        "run_mode": run_mode,
        "full_clean_run_completed": bool(run_mode.get("full_clean_run_completed")),
        "used_cached_v3_outputs": bool(run_mode.get("used_cached_v3_outputs")),
        "used_cached_msa_outputs": bool(run_mode.get("used_cached_msa_outputs")),
        "species_list_resolved": run_mode.get("species_list_resolved", ""),
        "gate_status": "pass" if gate_checks and not gate_fail else ("fail" if gate_fail else "unknown"),
        "gate_checks_total": len(gate_checks),
        "gate_checks_failed": len(gate_fail),
        "gate_failed_checks": [c.get("check") for c in gate_fail],
        "kpi": {
            "species": len(species),
            "isoform_rows": len(rows),
            "primary_ready": len(primary),
            "review_only": len(review),
            "excluded": len(excluded),
            "rescued_validated": len(rescued),
            "primary_fasta_sequences": count_fasta(primary_fasta),
            "review_fasta_sequences": count_fasta(review_fasta),
        },
        "freeze_ready": primary_fasta.exists(),
        "steps": steps,
    }


def build_species_index(run_dir: Path) -> List[Dict[str, Any]]:
    rows = _truth(run_dir)
    review_expl = {(r.get("species", ""), r.get("isoform", "")): r
                   for r in read_tsv(run_dir / "tables" / "final_review_case_explanation.tsv")}
    stack = read_tsv(run_dir / "tables" / "figure_final_framework_evidence_stack.tsv")
    ref_by = {(x.get("species", ""), x.get("isoform", "")): x for x in stack}
    sources = _evidence_sources(run_dir)
    by_species: Dict[str, Dict[str, Any]] = {}
    for r in rows:
        sp = r.get("species", "")
        if not sp:
            continue
        iso = r.get("isoform", "")
        readiness = r.get("pre_interpro_readiness_class", "")
        claim = r.get("final_claim_status_after_rescue", "") or r.get("final_claim_status", "")
        expl = review_expl.get((sp, iso), {})
        ready = readiness_status_ui(readiness)
        layers = compute_layers(r, ref_by.get((sp, iso), {}), sources)
        iso_entry = {
            "isoform": iso,
            "final_isoform_label": r.get("final_isoform_label", ""),
            "transcript_id": r.get("transcript_id", ""),
            "protein_id": r.get("protein_id", ""),
            "gene_id": r.get("gene_id", ""),
            "protein_length": to_int(r.get("protein_length")),
            "sequence_md5": r.get("sequence_md5", ""),
            "final_claim_status_after_rescue": claim,
            "claim_class": claim_class(claim),
            "pre_interpro_readiness_class": readiness,
            "readiness_class": readiness_class(readiness),
            # user-facing final status (single source of truth for Summary/Isoform)
            "interpro_included": ready["included"],
            "readiness_status": ready,
            "layers": layers,
            "orthology_status": r.get("orthology_status", ""),
            "synteny_validation_class": r.get("combined_synteny_validation_class") or r.get("synteny_validation_class", ""),
            "MSA_full_length_status": r.get("MSA_full_length_status", ""),
            "MSA_cassette_status": r.get("MSA_cassette_status", ""),
            "coordinate_validation_status": r.get("coordinate_validation_status", ""),
            "boundary_robustness_class": r.get("boundary_robustness_class", ""),
            "protein_integrity_status": r.get("protein_integrity_status", ""),
            "rescue_decision": r.get("rescue_decision", ""),
            "final_label_source": r.get("final_label_source", ""),
            "pre_interpro_warning": r.get("pre_interpro_warning", ""),
            "unresolved_reason_if_any": r.get("unresolved_reason_if_any", ""),
            # provenance (detail only, never overrides final status)
            "provenance": {
                "upstream_label": r.get("upstream_label", ""),
                "legacy_label": r.get("legacy_label", ""),
                "previous_pipeline_label": r.get("previous_pipeline_label", ""),
                "validated_exon_type": r.get("validated_exon_type", ""),
                "label_consistency_status": r.get("label_consistency_status", ""),
                "rescue_required": r.get("rescue_required", ""),
                "final_label_source": r.get("final_label_source", ""),
            },
            "upstream_corrected": "swapped" in (r.get("label_consistency_status", "") or "").lower(),
            "review_explanation": expl.get("final_interpretation", ""),
        }
        node = by_species.setdefault(sp, {
            "species": sp,
            "display_species_name": r.get("display_species_name", "") or sp.replace("_", " ").title(),
            "taxon_group": r.get("taxon_group", ""),
            "isoforms": [],
        })
        node["isoforms"].append(iso_entry)

    out: List[Dict[str, Any]] = []
    for sp, node in by_species.items():
        node["isoforms"].sort(key=lambda x: x["isoform"])
        classes = [i["readiness_class"] for i in node["isoforms"]]
        if any(c == "excluded" for c in classes):
            overall = "excluded"
        elif any(c == "review" for c in classes):
            overall = "review"
        elif any(c == "minor" for c in classes):
            overall = "minor"
        else:
            overall = "accepted"
        node["overall_readiness"] = overall
        node["synteny_summary"] = node["isoforms"][0].get("synteny_validation_class", "") if node["isoforms"] else ""
        node["msa_warning"] = any(token_class(i.get("MSA_full_length_status", "")) not in ("accepted",)
                                  or i.get("cassette_MSA_warning") for i in node["isoforms"])
        out.append(node)
    out.sort(key=lambda x: (x["taxon_group"], x["display_species_name"]))
    return out


def build_evidence_stack(run_dir: Path) -> Dict[str, Any]:
    """Per-row biological decision layers.

    Driven by the FINAL truth table (authoritative status + provenance) and
    enriched with the framework evidence stack's raw tokens + numeric reference
    agreement. Upstream/legacy labels are provenance only and never set status.
    """
    truth = _truth(run_dir)
    stack = read_tsv(run_dir / "tables" / "figure_final_framework_evidence_stack.tsv")
    ref_by = {(x.get("species", ""), x.get("isoform", "")): x for x in stack}
    sources = _evidence_sources(run_dir)
    columns = [{"key": k, "label": lbl} for k, lbl in EVIDENCE_LAYERS]

    rows = []
    for r in truth:
        sp, iso = r.get("species", ""), r.get("isoform", "")
        if not sp or not iso:
            continue
        ref = ref_by.get((sp, iso), {})
        rcls = readiness_class(r.get("pre_interpro_readiness_class", ""))
        ready = readiness_status_ui(r.get("pre_interpro_readiness_class", ""))
        consistency = (r.get("label_consistency_status", "") or "").lower()
        rows.append({
            "species": sp,
            "display_species_name": r.get("display_species_name", "") or sp.replace("_", " ").title(),
            "taxon_group": r.get("taxon_group", ""),
            "isoform": iso,
            "final_isoform_label": r.get("final_isoform_label", ""),
            "validated_exon_type": r.get("validated_exon_type", ""),
            "final_claim_status_after_rescue": r.get("final_claim_status_after_rescue", ""),
            "readiness_class": rcls,
            "row_class": rcls,
            "interpro_included": ready["included"],
            "upstream_corrected": "swapped" in consistency,
            "upstream_consistent": ("consistent" in consistency or "agree" in consistency),
            "visual_review_flag": rcls == "review",
            "provenance": {
                "upstream_label": r.get("upstream_label", ""),
                "legacy_label": r.get("legacy_label", ""),
                "previous_pipeline_label": r.get("previous_pipeline_label", ""),
                "final_label_source": r.get("final_label_source", ""),
                "label_consistency_status": r.get("label_consistency_status", ""),
                "rescue_decision": r.get("rescue_decision", ""),
                "validated_exon_type": r.get("validated_exon_type", ""),
            },
            "cells": compute_layers(r, ref, sources),
        })
    return {"columns": columns, "rows": rows}


# --- figures ---------------------------------------------------------------- #
def _parse_captions(run_dir: Path) -> Dict[str, str]:
    md = run_dir / "reports" / "final_pre_interpro_figure_captions.md"
    out: Dict[str, str] = {}
    if not md.exists():
        return out
    text = md.read_text(encoding="utf-8", errors="replace")
    # join soft-wrapped caption lines into single paragraphs per bold marker
    for m in re.finditer(r"\*\*(Figure\s+([0-9A-Za-z]+)|Supplement)\*\*\s*[—-]\s*(.+?)(?=\n\*\*|\Z)",
                         text, flags=re.DOTALL):
        num = (m.group(2) or "SUPP").upper()
        cap = " ".join(m.group(3).split())
        out[num] = cap
    return out


def _figure_group(stem: str) -> str:
    s = stem.lower()
    if s.startswith("supplement"):
        return "Supplement"
    if "synteny" in s or re.search(r"figure_9", s):
        return "Synteny"
    if ("msa" in s or "residue" in s or "discriminating" in s
            or re.search(r"figure_(5|6|7)\b", s) or "cassette_residue" in s):
        return "Sequence & MSA"
    return "Framework"


def _figure_number(stem: str) -> Optional[str]:
    m = re.match(r"figure_([0-9]+[A-Za-z]?)_", stem, flags=re.IGNORECASE)
    return m.group(1).upper() if m else None


def _withdrawn_figure_stems() -> set:
    """Figures the FGFR2 catalogue has withdrawn from every visible Gallery.

    Imported lazily from the catalogue module so the list lives in exactly one
    place; if it cannot be imported the index is built unfiltered rather than
    failing, because a missing filter is a cosmetic regression and a failed index
    build takes the whole dataset offline.
    """
    try:
        from fgfr2.gallery_catalogue import withdrawn_figure_stems
        return withdrawn_figure_stems()
    except Exception:
        return set()


def build_figure_index(run_dir: Path) -> Dict[str, Any]:
    figdir = run_dir / "figures"
    tabledir = run_dir / "tables"
    captions = _parse_captions(run_dir)
    table_stems = {p.stem.lower(): p for p in tabledir.glob("*.tsv")} if tabledir.exists() else {}

    by_stem: Dict[str, Dict[str, Any]] = {}
    if figdir.exists():
        for p in sorted(figdir.iterdir()):
            if p.suffix.lower() not in (".png", ".svg", ".pdf"):
                continue
            stem = p.stem
            entry = by_stem.setdefault(stem, {"stem": stem, "formats": {}})
            entry["formats"][p.suffix.lower().lstrip(".")] = rel(p)

    withdrawn = _withdrawn_figure_stems()

    figures: List[Dict[str, Any]] = []
    for stem, entry in by_stem.items():
        # A withdrawn figure is still on disk and still readable as validated output;
        # it simply is not a Gallery card any more. Skipping it here rather than
        # deleting the file is what keeps downloads and QC intact.
        if stem in withdrawn:
            continue
        num = _figure_number(stem)
        group = _figure_group(stem)
        is_supp = stem.lower().startswith("supplement")
        title = _humanise_figure_title(stem, num)
        caption = captions.get(num or "", "") if num else ""
        if not caption and is_supp:
            caption = captions.get("SUPP", "")
        # match a source table by figure number prefix (figure3B..., figure6...)
        src_table = ""
        if num:
            key = f"figure{num.lower()}_"
            for tstem, tpath in table_stems.items():
                if tstem.startswith(key):
                    src_table = rel(tpath)
                    break
        figures.append({
            "id": stem,
            "title": title,
            "number": num or "",
            "group": group,
            "kind": "supplement" if is_supp else "main",
            "caption": caption,
            "formats": entry["formats"],
            "thumbnail": entry["formats"].get("png", ""),
            "source_table": src_table,
        })

    # stable, reader-friendly ordering
    def sort_key(f):
        n = f["number"]
        m = re.match(r"(\d+)([A-Za-z]?)", n)
        major = int(m.group(1)) if m else 999
        minor = m.group(2) if m else ""
        return (0 if f["kind"] == "main" else 1, major, minor, f["title"])

    # Post-InterPro / pyTMHMM domain architecture figures (sibling step-15 folder)
    arch_group = "Domain & exon-boundary"
    has_arch = _append_architecture_figures(run_dir, figures, arch_group)

    # Module 1 — exon-domain boundary consistency (sibling step-16 folder)
    bc_group = "Boundary consistency"
    has_bc = _append_boundary_consistency_figures(run_dir, figures, bc_group)

    figures.sort(key=sort_key)
    groups = ["Framework", "Sequence & MSA", "Synteny"]
    if has_arch:
        groups.append(arch_group)
    if has_bc:
        groups.append(bc_group)
    groups.append("Supplement")
    return {
        "groups": groups,
        "figures": figures,
        "curated_ids": [f["id"] for f in figures if f["kind"] == "main"],
    }


def _append_architecture_figures(run_dir: Path, figures: List[Dict[str, Any]],
                                 group: str) -> bool:
    """Add step-15 domain/exon-boundary overview + per-species plots to the gallery.

    Overview panels are 'main' (always visible); per-species plots are
    'supplement' (shown when the gallery's supplement toggle is on). Static
    figures are downloadable previews; the interactive Domain architecture tab
    is the primary view (UI rule).
    """
    arch = _post_interpro_dir(run_dir)
    if not arch:
        return False
    feat_tbl = arch / "tables" / "exon_domain_architecture_features.tsv"
    src_table = rel(feat_tbl) if feat_tbl.exists() else ""

    def _formats(base_dir: Path, stem: str) -> Dict[str, str]:
        out: Dict[str, str] = {}
        for ext in ("png", "svg", "pdf"):
            p = base_dir / f"{stem}.{ext}"
            if p.exists():
                out[ext] = rel(p)
        return out

    overview_defs = [
        ("Figure_10_all_species_FGFR2_exon_domain_architecture_primary",
         "Figure 10 — All species · FGFR2 exon–domain architecture"),
        ("Figure_10A_IIIb_exon_domain_architecture_primary",
         "Figure 10A — IIIb · exon–domain architecture"),
        ("Figure_10B_IIIc_exon_domain_architecture_primary",
         "Figure 10B — IIIc · exon–domain architecture"),
        ("Figure_10C_mammals_exon_domain_architecture_primary",
         "Figure 10C — Mammals · exon–domain architecture"),
        ("Figure_10D_nonmammals_exon_domain_architecture_primary",
         "Figure 10D — Non-mammals · exon–domain architecture"),
    ]
    ov_dir = arch / "figures" / "overview"
    added = False
    for stem, title in overview_defs:
        fmts = _formats(ov_dir, stem)
        if not fmts:
            continue
        figures.append({
            "id": stem, "title": title, "number": "10", "group": group,
            "kind": "main",
            "caption": ("Post-InterProScan exon–domain architecture: InterProScan domains, "
                        "pyTMHMM transmembrane helix, numbered coding exons and the IIIb/IIIc "
                        "cassette slot. Labels come from the final truth table; pyTMHMM is the "
                        "TM layer. Interactive view: Gene Explorer → Domain architecture."),
            "formats": fmts, "thumbnail": fmts.get("png", ""), "source_table": src_table,
        })
        added = True

    ps_dir = arch / "figures" / "per_species"
    if ps_dir.exists():
        stems = sorted({p.stem for p in ps_dir.glob("*.png")})
        for stem in stems:
            fmts = _formats(ps_dir, stem)
            if not fmts:
                continue
            pretty = stem.replace("_exon_domain_architecture", "").replace("_", " ").strip()
            pretty = pretty[:1].upper() + pretty[1:]
            figures.append({
                "id": stem, "title": f"{pretty} — exon–domain architecture",
                "number": "10", "group": group, "kind": "supplement",
                "caption": ("Per-protein exon–domain architecture (InterProScan domains, pyTMHMM "
                            "TM, numbered coding exons, IIIb/IIIc cassette slot)."),
                "formats": fmts, "thumbnail": fmts.get("png", ""), "source_table": src_table,
            })
            added = True
    return added


def _append_boundary_consistency_figures(run_dir: Path, figures: List[Dict[str, Any]],
                                          group: str) -> bool:
    """Add the Module 1 boundary-consistency static figures (heatmap + distance
    distribution) to the gallery as downloadable thesis previews. The interactive
    Boundary Consistency explorer is the primary web experience (UI rule)."""
    bc = _boundary_consistency_dir(run_dir)
    if not bc:
        return False
    src_table = bc / "tables" / "exon_domain_boundary_distances.tsv"
    src = rel(src_table) if src_table.exists() else ""

    def _formats(stem: str) -> Dict[str, str]:
        out: Dict[str, str] = {}
        for ext in ("png", "svg", "pdf"):
            p = bc / "figures" / f"{stem}.{ext}"
            if p.exists():
                out[ext] = rel(p)
        return out

    defs = [
        ("Figure_11_exon_domain_boundary_consistency_heatmap", "11",
         "Figure 11 — FGFR2 IIIb/IIIc exon\u2013domain boundary consistency",
         "Cassette start/end boundary class (aligned / near / within / between / "
         "missing) relative to the nearest protein-domain boundary, per species and "
         "isoform. Rows are ordered taxonomically. Static thesis figure — the "
         "interactive Boundary Consistency explorer is the primary view."),
        ("Figure_12_boundary_distance_distribution", "12",
         "Figure 12 — Cassette boundary distance to nearest domain boundary",
         "Distribution of the absolute distance from cassette start/end boundaries "
         "to the nearest protein-domain boundary, by isoform. Static thesis figure — "
         "the interactive Boundary Consistency explorer is the primary view."),
    ]
    added = False
    for stem, num, title, caption in defs:
        fmts = _formats(stem)
        if not fmts:
            continue
        figures.append({
            "id": stem, "title": title, "number": num, "group": group,
            "kind": "main", "caption": caption,
            "formats": fmts, "thumbnail": fmts.get("png", ""), "source_table": src,
        })
        added = True
    return added


def _humanise_figure_title(stem: str, num: Optional[str]) -> str:
    s = stem
    s = re.sub(r"^(Figure_|Supplement_Figure_|Supplement_)", "", s)
    s = re.sub(r"^([0-9]+[A-Za-z]?)_", "", s)
    s = s.replace("_pre_interpro", "").replace("_paper", "")
    s = s.replace("_", " ").strip()
    s = s[:1].upper() + s[1:] if s else stem
    if num:
        prefix = "Suppl." if stem.lower().startswith("supplement") else "Figure"
        return f"{prefix} {num} — {s}"
    if stem.lower().startswith("supplement"):
        return f"Supplement — {s}"
    return s


def build_download_index(run_dir: Path) -> List[Dict[str, Any]]:
    items: List[Dict[str, Any]] = []

    def add(group: str, label: str, path: Path, fmt: str = ""):
        meta = file_meta(path)
        if meta:
            items.append({"group": group, "label": label, "format": fmt or Path(path).suffix.lstrip("."), **meta})

    add("Freeze", "Primary protein FASTA (InterPro input)", run_dir / "freeze" / "final_pre_interpro_proteins_primary.faa")
    add("Freeze", "Review-included protein FASTA", run_dir / "freeze" / "final_pre_interpro_proteins_all_review_included.faa")
    add("Freeze", "Sequence manifest (TSV)", run_dir / "freeze" / "final_pre_interpro_sequence_manifest.tsv")
    add("Freeze", "File checksums", run_dir / "freeze" / "final_pre_interpro_file_checksums.tsv")
    add("Tables", "Final truth table", run_dir / "final_pre_interpro_truth_table.tsv")
    add("Tables", "Evidence stack", run_dir / "tables" / "figure_final_framework_evidence_stack.tsv")
    add("Tables", "Review-case explanations", run_dir / "tables" / "final_review_case_explanation.tsv")
    add("Reports", "QC report", run_dir / "reports" / "final_pre_interpro_QC_report.md")
    add("Reports", "Results summary", run_dir / "reports" / "final_pre_interpro_results_summary.md")
    add("Reports", "Methods summary", run_dir / "reports" / "final_pre_interpro_methods_summary.md")
    add("Reports", "Figure captions", run_dir / "reports" / "final_pre_interpro_figure_captions.md")
    add("Gate", "Cross-table consistency gate", run_dir / "gates" / "final_pre_interpro_cross_table_consistency_gate.tsv")

    # The scientific tables and alignments the views are drawn from. Offering only the
    # closure summary meant a finished run exposed its conclusion and none of the
    # evidence: no gene models, no cassette map, no alignment, no domain annotation.
    add("Tables", "Cassette coordinate map (exon → protein)",
        run_dir / "tables" / "figure3C_exon_to_protein_cassette_coordinate_map.tsv")
    add("Tables", "Cassette amino-acid motif map",
        run_dir / "tables" / "figure3B_IIIb_IIIc_cassette_amino_acid_motif_map.tsv")
    add("Tables", "Cassette residue map (species-resolved)",
        run_dir / "tables" / "figure6B_species_resolved_IIIb_IIIc_cassette_residue_map.tsv")
    add("Tables", "Exon → protein architecture",
        run_dir / "tables" / "figure2_final_exon_to_protein_architecture_pre_interpro.tsv")
    add("Tables", "Label reconciliation and rescue",
        run_dir / "tables" / "figure4_label_reconciliation_and_rescue_summary.tsv")
    for key, label in (("final_fgfr2_full_length_protein_msa.aln.faa",
                        "Full-length protein alignment (FASTA)"),
                       ("final_fgfr2_IIIb_cassette_msa.aln.faa",
                        "IIIb cassette alignment (FASTA)"),
                       ("final_fgfr2_IIIc_cassette_msa.aln.faa",
                        "IIIc cassette alignment (FASTA)"),
                       ("final_fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa",
                        "Combined IIIb/IIIc cassette alignment (FASTA)"),
                       ("final_full_length_msa_conservation_summary.tsv",
                        "Alignment conservation summary"),
                       ("final_isoform_discriminating_residues.tsv",
                        "Isoform-discriminating residues")):
        add("Alignment", label, run_dir / "MSA" / key)

    # Stages outside the closure directory. The closure is a freeze of conclusions; the
    # run's own model tables and returned cluster annotation live beside it.
    results = _results_root(run_dir)
    if results is not None:
        for rel_path, group, label in (
            ("02_models/genes.tsv", "Gene models", "Genes"),
            ("02_models/transcripts.tsv", "Gene models", "Transcripts"),
            ("02_models/exons.tsv", "Gene models", "Coding exons"),
            ("02_models/cds_features.tsv", "Gene models", "CDS features"),
            ("02_models/gene_candidates.tsv", "Gene models", "Gene candidate inventory"),
            ("02_models/ncbi_assembly_selection.tsv", "Gene models",
             "Assembly candidates and decisions"),
            ("02_models/internal_consistency_checks.tsv", "Gene models",
             "Internal consistency checks"),
            ("05b_selection_with_isoforms_v2_7_marker_validated/selected_transcripts.tsv",
             "Selected models", "Selected transcripts"),
            ("05b_selection_with_isoforms_v2_7_marker_validated/"
             "fgfr2_III_final_selected_protein_validation_summary.tsv",
             "Selected models", "Selected-protein validation"),
            # The header of this file carries the *upstream* selection role, assigned
            # before marker-based reconciliation. In this dataset that role is
            # systematically inverted relative to the final label, so the file is
            # offered as provenance and the name says so. The final label lives in the
            # truth table and in the alignment headers.
            ("06_protein_export_v2_7_marker_validated/selected_fgfr2_proteins.faa",
             "Selected models",
             "Exported proteins with upstream selection role (provenance — final "
             "isoform label is in the truth table)"),
            ("06b_paralog_screen_v2_7_marker_validated/fgfr2_paralog_screen_detailed.tsv",
             "Selected models", "Paralog screen"),
            ("14_interproscan/primary/output/input.fasta.tsv", "Domain annotation",
             "InterProScan results (TSV)"),
            ("14_interproscan/primary/output/input.fasta.gff3", "Domain annotation",
             "InterProScan results (GFF3)"),
            ("15_exon_domain_boundary_post_interpro/pytmhmm_primary/output/"
             "pytmhmm_transmembrane_hits.tsv", "Domain annotation",
             "pyTMHMM transmembrane hits"),
            ("15_exon_domain_boundary_post_interpro/tables/"
             "exon_domain_architecture_features.tsv", "Domain annotation",
             "Exon–domain architecture features"),
            ("15_exon_domain_boundary_post_interpro/tables/"
             "interpro_domain_features_normalized.tsv", "Domain annotation",
             "Normalized InterPro domain features"),
            ("16_final_thesis_analyses/exon_domain_boundary_consistency/tables",
             "Boundary analysis", "Exon–domain boundary"),
            ("16_final_thesis_analyses/final_audit/tables",
             "Boundary analysis", "Final audit"),
        ):
            target = results / rel_path
            if target.is_dir():
                for found in sorted(target.rglob("*.tsv")):
                    add(group, f"{label} — {found.stem.replace('_', ' ')}", found)
            else:
                add(group, label, target)

    # Figures and the source table behind each, so a figure can be checked rather than
    # only looked at.
    figures = run_dir / "figures"
    if figures.is_dir():
        for svg in sorted(figures.glob("*.svg")):
            add("Figures", svg.stem.replace("_", " "), svg)

    # archive (timestamped)
    arch = run_dir / "archive"
    if arch.exists():
        for z in sorted(arch.glob("*.zip")):
            add("Archive", "Reproducibility freeze archive", z)
    return items


def _results_root(run_dir: Path) -> Optional[Path]:
    """The ``results/`` directory a closure directory belongs to.

    Indices are built from the closure directory, but a run's evidence is spread across
    its sibling stage folders. Deriving the root from the closure keeps the builder
    callable with a single argument, exactly as before.
    """
    parent = Path(run_dir).parent
    return parent if (parent / "02_models").is_dir() else None


def build_freeze_index(run_dir: Path) -> Dict[str, Any]:
    run_mode = read_json(run_dir / "final_pre_interpro_run_mode.json", {}) or {}
    primary = run_dir / "freeze" / "final_pre_interpro_proteins_primary.faa"
    review = run_dir / "freeze" / "final_pre_interpro_proteins_all_review_included.faa"
    manifest = run_dir / "freeze" / "final_pre_interpro_sequence_manifest.tsv"
    checksums = run_dir / "freeze" / "final_pre_interpro_file_checksums.tsv"

    checksum_rows = read_tsv(checksums)
    gate = read_json(run_dir / "gates" / "final_pre_interpro_cross_table_consistency_gate.json", {}) or {}
    gate_checks = gate.get("checks", []) if isinstance(gate, dict) else []
    md5_check = next((c for c in gate_checks if "md5" in str(c.get("check", "")).lower()), None)

    if run_mode.get("full_clean_run_completed"):
        run_mode_label = "full clean run"
    elif run_mode.get("used_cached_v3_outputs") or run_mode.get("used_cached_msa_outputs"):
        run_mode_label = "cached / freeze-based"
    else:
        run_mode_label = "unknown"

    cards = []

    def card(key, title, path, role, count=None):
        meta = file_meta(path)
        cards.append({
            "key": key, "title": title, "role": role,
            "available": meta is not None,
            "sequences": count,
            **(meta or {"name": Path(path).name, "path": rel(path), "size_human": "—"}),
        })

    n_primary = count_fasta(primary)
    n_review = count_fasta(review)
    card("primary_fasta", "Primary InterPro FASTA", primary,
         "Recommended main InterProScan input — final primary-ready sequences.", n_primary)
    card("review_fasta", "Review-included FASTA", review,
         "Optional / supplementary input — adds supplement/review cases; not equal main evidence.", n_review)
    card("manifest", "Sequence manifest", manifest, "Per-sequence provenance + readiness")
    card("checksums", "Checksums", checksums, "Integrity of frozen artefacts")
    return {
        "run_mode_label": run_mode_label,
        "run_mode": run_mode,
        "checksum_count": len(checksum_rows),
        "checksum_gate": (md5_check or {}).get("status", "n/a"),
        "cards": cards,
        "interpro_policy": {
            "primary_count": n_primary,
            "review_included_count": n_review,
            "primary_role": "main recommended input",
            "review_role": "optional / supplementary input",
            "notes": [
                f"Primary FASTA contains the {n_primary} final primary-ready sequences and is the recommended main InterProScan analysis.",
                f"Review-included FASTA ({n_review} sequences) adds supplement/review cases and should not be treated as equal main evidence.",
                "InterProScan annotates protein domains; it does not validate IIIb/IIIc final claim status.",
                "Excluding the two review rows from the primary set is by design, not an error.",
            ],
        },
    }


# --------------------------------------------------------------------------- #
# Phase-2 interactive indices (all derived from final run outputs)
# --------------------------------------------------------------------------- #
def build_cassette_residue_index(run_dir: Path) -> Dict[str, Any]:
    """Module 1 source: reproduce Figure 6B interactively + human-reference + discriminating."""
    fig6b = read_tsv(run_dir / "tables" / "figure6B_species_resolved_IIIb_IIIc_cassette_residue_map.tsv")
    agree = read_tsv(run_dir / "MSA" / "final_human_referenced_residue_agreement.tsv")
    motif = read_tsv(run_dir / "tables" / "figure3B_IIIb_IIIc_cassette_amino_acid_motif_map.tsv")
    review = _review_lookup(run_dir)

    # enrichment join key: (species, isoform, human_reference_residue_index)
    agree_idx: Dict[tuple, Dict[str, str]] = {}
    for r in agree:
        agree_idx[(r.get("species", ""), r.get("isoform", ""),
                   r.get("human_reference_residue_index", ""))] = r

    # Human reference per panel. IIIb is 46 cassette residues and IIIc is 48, and each
    # numbers its own cassette. Deriving both from one shared index (one `seen_ref` set
    # over the motif map's single index column) forced them onto the same axis: IIIc was
    # truncated to IIIb's length, which destroyed its GVNTTDKEI marker, and positions
    # where one panel has an alignment gap lost the other panel's residue. The reference
    # therefore comes from the validated control, which stores the two panels separately.
    human_ref: Dict[str, List[Dict[str, Any]]] = {"IIIb": [], "IIIc": []}
    reference_status = "unavailable"
    try:
        control = human_reference_control.load()
        for panel in ("IIIb", "IIIc"):
            human_ref[panel] = [dict(res) for res in
                                human_reference_control.panel_residues(control, panel)]
        reference_status = "validated_control"
    except human_reference_control.ReferenceControlError:
        control = None

    # The discriminating layer is a IIIb-vs-IIIc comparison, so it lives on the combined
    # cassette alignment column — the one axis both panels genuinely share.
    discriminating: Dict[str, Dict[str, Any]] = {}
    seen_ref: set = set()
    _has_human_index = any(to_int(r.get("human_reference_residue_index")) is not None for r in motif)
    residue_index_basis = "combined_cassette_alignment_column"
    # A motif map that predates the two per-panel index columns (the read-only freeze,
    # and every run whose closure reused it) has them recovered from the combined
    # alignment first, so both the per-column payload below and the per-panel marker
    # sets speak the same coordinate systems.
    motif_indexed = motif if human_reference_control.has_panel_indices(motif) \
        else human_reference_control.panel_indices_from_combined_alignment(
            motif, column_fields=("MSA_column", "combined_alignment_col", "alignment_col"),
            aa_fields={"IIIb": ("human_IIIb_aa_one_letter", "human_IIIb_aa"),
                       "IIIc": ("human_IIIc_aa_one_letter", "human_IIIc_aa")})
    disc_by_panel = human_reference_control.discriminating_positions_by_panel(motif_indexed)
    for r in motif_indexed:
        column = to_int(r.get("MSA_column"))
        if column is None or column in seen_ref:
            continue
        seen_ref.add(column)
        iiib_aa = r.get("human_IIIb_aa_one_letter", "")
        iiic_aa = r.get("human_IIIc_aa_one_letter", "")
        if truthy(r.get("is_isoform_discriminating")):
            discriminating[str(column)] = {
                "i": column,
                "msa_column": column,
                "human_IIIb_reference_index": to_int(r.get("human_IIIb_reference_index")),
                "human_IIIc_reference_index": to_int(r.get("human_IIIc_reference_index")),
                "IIIb_aa": iiib_aa,
                "IIIc_aa": iiic_aa,
                "IIIb_property": r.get("IIIb_residue_property_class", ""),
                "IIIc_property": r.get("IIIc_residue_property_class", ""),
                "substitution_class": r.get("substitution_class_IIIb_vs_IIIc", ""),
                "position_class": r.get("position_class", ""),
                "discriminating_score": to_float(r.get("discriminating_score")),
            }
    human_ref["IIIb"].sort(key=lambda x: x["i"])
    human_ref["IIIc"].sort(key=lambda x: x["i"])
    # Species cells are keyed by each panel's OWN reference index, so the marker set
    # (`disc_by_panel`, derived above) is per panel too. Taking the union marked IIIb at
    # 16 and 17 — combined columns that exist only in IIIc, because of its two-residue
    # insertion — and shifted every IIIc mark behind it. The same rule serves Figure 6B.
    species: Dict[str, Dict[str, Any]] = {}
    for r in fig6b:
        sp = r.get("species", "")
        iso = r.get("isoform", "")
        if not sp or not iso:
            continue
        hri = to_int(r.get("human_reference_residue_index"))
        enr = agree_idx.get((sp, iso, r.get("human_reference_residue_index", "")), {})
        # The panel a row is drawn on decides which coordinate system its position is
        # in — the final label, since that is the biology the cell reports. The motif
        # map is the authority whenever it carries the per-panel indices; the table's
        # own flag is only a fallback for a run whose motif map predates them, and it
        # must not override the per-panel answer.
        panel = r.get("final_isoform_label") or iso
        is_disc = (hri in disc_by_panel.get(panel, set())) if any(disc_by_panel.values()) \
            else truthy(r.get("is_discriminating_position"))
        pos = {
            "i": hri,
            "h_aa": r.get("human_reference_aa", ""),
            "sp_aa": r.get("species_aa", "") or "-",
            "agreement_class": r.get("agreement_class", ""),
            "cls": residue_class(r.get("agreement_class", "")),
            "is_discriminating": is_disc,
            "substitution_class": enr.get("substitution_class", ""),
            "blosum": to_float(enr.get("blosum62_score_if_available")),
            "msa_column": to_int(enr.get("alignment_col")),
        }
        node = species.setdefault(sp, {
            "species": sp,
            "display_species_name": sp.replace("_", " ").title(),
            "panels": {},
        })
        panel = node["panels"].setdefault(iso, {
            "isoform": iso,
            "available": True,
            "final_isoform_label": r.get("final_isoform_label", ""),
            "validated_exon_type": r.get("validated_exon_type", ""),
            "final_claim_status_after_rescue": r.get("final_claim_status_after_rescue", ""),
            "claim_class": claim_class(r.get("final_claim_status_after_rescue", "")),
            "readiness_class": review.get((sp, iso), "neutral"),
            "is_review": review.get((sp, iso)) == "review",
            "source_table": rel(run_dir / "tables" / "figure6B_species_resolved_IIIb_IIIc_cassette_residue_map.tsv"),
            "positions": [],
        })
        panel["positions"].append(pos)

    out_species = []
    for sp, node in species.items():
        for panel in node["panels"].values():
            panel["positions"].sort(key=lambda x: (x["i"] is None, x["i"]))
            idxs = [p["i"] for p in panel["positions"] if p["i"] is not None]
            panel["cassette_start"] = min(idxs) if idxs else None
            panel["cassette_end"] = max(idxs) if idxs else None
            panel["n_identical"] = sum(1 for p in panel["positions"] if p["cls"] == "identical")
            panel["n_diff"] = sum(1 for p in panel["positions"] if p["cls"] in ("conservative", "nonconservative"))
            panel["n_gap"] = sum(1 for p in panel["positions"] if p["cls"] == "gap")
        node["panels"] = {k: node["panels"][k] for k in sorted(node["panels"])}
        out_species.append(node)
    out_species.sort(key=lambda x: x["display_species_name"])

    # Sequence / marker level cassette evidence. Always built so that runs which have
    # no residue-level human-referenced agreement table (figure6B empty — typical for
    # small custom runs) still expose cassette coordinates, length and marker status.
    sequence_evidence = _build_cassette_sequence_evidence(run_dir, review)

    human_reference_role = None
    comparison_source = ""
    note = ""
    if fig6b:
        # Full run with its own residue-level human-referenced agreement (example dataset):
        # keep the exact original behaviour.
        evidence_level = "residue_map"
    else:
        # Custom / single-species run: try to build a human-vs-species residue comparison by
        # aligning each run-local species cassette to the validated human FGFR2 reference
        # (allowed reference/control layer). This renders with the same Human-comparison,
        # heatmap and discriminating UI as the example dataset.
        cassette_ref = _human_reference_cassette(run_dir)
        comparison = _build_cassette_human_comparison(run_dir, cassette_ref, review) if cassette_ref else None
        if comparison:
            evidence_level = "human_reference_msa_comparison"
            residue_index_basis = "msa_column"
            human_reference_role = "reference_control_only"
            human_ref = comparison["human_reference"]
            discriminating = comparison["discriminating"]
            out_species = comparison["species"]
            comparison_source = comparison.get("source_table", "")
            note = ("Custom run comparison is indexed by cassette MSA column; "
                    "human reference is shown for orientation only.")
        elif sequence_evidence:
            evidence_level = "sequence_marker"
            note = ("Cassette evidence available from sequence/MSA marker layer; "
                    "exact residue-level human-referenced agreement is not available for this run. "
                    "The human FGFR2 IIIb/IIIc reference is shown for orientation only.")
        else:
            evidence_level = "none"

    def _last_i(rows: List[Dict[str, Any]]):
        idxs = [r["i"] for r in rows if r.get("i") is not None]
        return max(idxs) if idxs else None

    _cassette_available = (evidence_level in ("residue_map",
                                              "human_reference_msa_comparison")
                           or bool(sequence_evidence))
    return {
        "available": _cassette_available,
        # The cassette layer has three source tables and degrades between them, so its
        # diagnosis names the coordinate map: the one table that must exist for any
        # cassette evidence at all, residue-level or marker-level.
        "availability": _availability_block(
            "cassette",
            run_dir / "tables" / "figure3C_exon_to_protein_cassette_coordinate_map.tsv",
            _cassette_available),
        "evidence_level": evidence_level,
        "residue_index_basis": residue_index_basis,
        "human_reference_role": human_reference_role,
        "human_reference_status": reference_status,
        "note": note,
        "panels": ["IIIb", "IIIc"],
        "human_reference": human_ref,
        "discriminating": discriminating,
        "cassette_length": {p: _last_i(hr) for p, hr in human_ref.items()},
        "species": out_species,
        "sequence_evidence": sequence_evidence,
        "source_tables": {
            "residue_map": rel(run_dir / "tables" / "figure6B_species_resolved_IIIb_IIIc_cassette_residue_map.tsv"),
            "agreement": rel(run_dir / "MSA" / "final_human_referenced_residue_agreement.tsv"),
            "motif_map": rel(run_dir / "tables" / "figure3B_IIIb_IIIc_cassette_amino_acid_motif_map.tsv"),
            "coordinate_map": rel(run_dir / "tables" / "figure3C_exon_to_protein_cassette_coordinate_map.tsv"),
            "cassette_zoom": rel(run_dir / "tables" / "figure3_final_IIIb_IIIc_cassette_zoom_pre_interpro.tsv"),
            "comparison": comparison_source or rel(run_dir / "MSA" / "final_fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa"),
        },
    }


def _build_cassette_sequence_evidence(run_dir: Path, review: Dict[tuple, str]) -> List[Dict[str, Any]]:
    """Per species/isoform cassette evidence from the sequence/marker layer.

    Sources (all run-local, no cross-run reuse):
      * figure3C exon→protein cassette coordinate map — cassette AA span + protein length
      * figure3 cassette zoom — cassette length (aa) and MSA column span
      * final pre-InterPro truth table — transcript / protein IDs + validation status
    """
    fig3c = read_tsv(run_dir / "tables" / "figure3C_exon_to_protein_cassette_coordinate_map.tsv")
    zoom = read_tsv(run_dir / "tables" / "figure3_final_IIIb_IIIc_cassette_zoom_pre_interpro.tsv")
    truth = _truth(run_dir)

    zoom_idx = {(r.get("species", ""), r.get("isoform", "")): r for r in zoom}
    truth_idx = {(r.get("species", ""), r.get("isoform", "")): r for r in truth}

    coord_src = rel(run_dir / "tables" / "figure3C_exon_to_protein_cassette_coordinate_map.tsv")
    by: Dict[tuple, Dict[str, Any]] = {}
    for r in fig3c:
        sp, iso = r.get("species", ""), r.get("isoform", "")
        if not sp or not iso:
            continue
        key = (sp, iso)
        node = by.get(key)
        if node is None:
            z = zoom_idx.get(key, {})
            t = truth_idx.get(key, {})
            cass_start = to_int(r.get("cassette_start_aa"))
            cass_end = to_int(r.get("cassette_end_aa"))
            node = {
                "species": sp,
                "isoform": iso,
                "display_species_name": sp.replace("_", " ").title(),
                "final_isoform_label": r.get("final_isoform_label", "") or iso,
                "validated_exon_type": r.get("validated_exon_type", ""),
                "protein_length": to_int(r.get("protein_length")),
                "transcript_id": t.get("transcript_id", ""),
                "protein_id": t.get("protein_id", ""),
                "cassette_start_aa": cass_start,
                "cassette_end_aa": cass_end,
                "cassette_available": cass_start is not None and cass_end is not None and cass_end > 0,
                "cassette_length_aa": to_int(z.get("cassette_length_aa")),
                "cassette_msa_start_col": to_int(z.get("cassette_msa_start_col")),
                "cassette_msa_end_col": to_int(z.get("cassette_msa_end_col")),
                "final_claim_status_after_rescue": r.get("final_claim_status_after_rescue", ""),
                "claim_class": claim_class(r.get("final_claim_status_after_rescue", "")),
                "final_plot_status": r.get("final_plot_status", ""),
                "visual_review_flag": r.get("visual_review_flag", ""),
                "readiness_class": review.get(key, "neutral"),
                "is_review": review.get(key) == "review",
                "cassette_exons": [],
                "source_table": coord_src,
            }
            by[key] = node
        if truthy(r.get("is_IIIb_cassette")) or truthy(r.get("is_IIIc_cassette")):
            node["cassette_exons"].append({
                "id": r.get("exon_or_cds_id", ""),
                "label": r.get("feature_label", ""),
                "start": to_int(r.get("block_start_aa")),
                "end": to_int(r.get("block_end_aa")),
            })

    out = list(by.values())
    out.sort(key=lambda x: (x["display_species_name"], x["isoform"]))
    return out


# --------------------------------------------------------------------------- #
# Run-local human-vs-species cassette comparison.
#
# Custom / single-species runs do not carry the residue-level human-referenced
# agreement table (figure6B is empty) because the human reference is NOT part of
# the run's cassette MSA. To still render the same Human-comparison UI as the
# validated example, we align each species cassette to the *validated human
# FGFR2 IIIb/IIIc reference* (the only allowed reused layer) with BLOSUM62 and
# emit the exact schema the frontend already consumes. No non-human example
# result rows are used; the comparison is computed fresh from run-local cassettes.
# --------------------------------------------------------------------------- #
_ALIGNER_CACHE: Dict[str, Any] = {}


def _cassette_aligner():
    if "a" in _ALIGNER_CACHE:
        return _ALIGNER_CACHE["a"]
    try:
        from Bio.Align import PairwiseAligner, substitution_matrices
        mat = substitution_matrices.load("BLOSUM62")
        al = PairwiseAligner()
        al.substitution_matrix = mat
        al.open_gap_score = -10.0
        al.extend_gap_score = -0.5
        al.mode = "global"
        res = (al, mat)
    except Exception:  # noqa: BLE001 — Biopython missing → graceful fallback to cards
        res = (None, None)
    _ALIGNER_CACHE["a"] = res
    return res


def _read_fasta(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    p = Path(path)
    if not p.exists():
        return out
    name = None
    buf: List[str] = []
    for line in p.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(buf)
            name = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf)
    return out


def _human_reference_cassette(run_dir: Path) -> Optional[Dict[str, Any]]:
    """Locate the curated human IIIb/IIIc cassette reference layer (control only)."""
    p = Path(run_dir).resolve()
    for base in [p, *p.parents]:
        cand = base / "results" / "web_state" / "human_reference_control.json"
        if cand.exists():
            try:
                data = json.loads(cand.read_text(encoding="utf-8"))
            except Exception:  # noqa: BLE001
                return None
            cass = data.get("cassette") or {}
            hr = cass.get("human_reference") or {}
            if hr.get("IIIb") or hr.get("IIIc"):
                return cass
            return None
    return None


def _classify_cassette_residue(h_aa: str, sp_aa: str, blosum_fn) -> tuple:
    """(agreement_class, cls, blosum, substitution_class) — mirrors the pipeline classifier."""
    if not sp_aa or sp_aa in ("-", ".", "*"):
        return "gap_or_missing", "gap", None, "gap"
    if not h_aa or h_aa in ("-", ".", "*"):
        return "unmapped_review", "nonconservative", None, "unknown"
    sc = blosum_fn(h_aa, sp_aa)
    if sp_aa == h_aa:
        return "identical_to_human", "identical", sc, "identical"
    if sc is None:
        return "nonconservative_substitution", "nonconservative", None, "nonconservative"
    if sc > 0:
        return "conservative_substitution", "conservative", sc, "conservative"
    if sc == 0:
        return "nonconservative_substitution", "conservative", sc, "semi_conservative"
    return "nonconservative_substitution", "nonconservative", sc, "nonconservative"


def _pairwise_positions(aligner, ref_seq: str, ref_ng: List[tuple], sp_seq: str,
                        discriminating: Dict[str, Any], blosum_fn) -> List[Dict[str, Any]]:
    """Align one species cassette to the human reference cassette; emit per-position rows
    keyed by the human reference residue index (same coordinate as the example UI)."""
    try:
        aln = aligner.align(ref_seq, sp_seq)[0]
        ref_aln, sp_aln = str(aln[0]), str(aln[1])
    except Exception:  # noqa: BLE001
        return []
    positions: List[Dict[str, Any]] = []
    ri = 0
    for rc, sc in zip(ref_aln, sp_aln):
        if rc in ("-", ".", " "):
            continue  # insertion in species relative to human reference: not a ref column
        if ri >= len(ref_ng):
            break
        i, h_aa = ref_ng[ri]
        ri += 1
        sp_aa = "" if sc in ("-", ".", " ") else sc
        agreement_class, cls, blo, subs = _classify_cassette_residue(h_aa, sp_aa, blosum_fn)
        positions.append({
            "i": i,
            "h_aa": h_aa,
            "sp_aa": sp_aa or "-",
            "agreement_class": agreement_class,
            "cls": cls,
            "is_discriminating": str(i) in discriminating,
            "substitution_class": subs,
            "blosum": blo,
            "msa_column": None,
        })
    return positions


def _build_cassette_human_comparison(run_dir: Path, cassette_ref: Dict[str, Any],
                                     review: Dict[tuple, str]) -> Optional[Dict[str, Any]]:
    aligner, blosum_mat = _cassette_aligner()
    if aligner is None:
        return None

    def blosum_fn(a: str, b: str) -> Optional[float]:
        if blosum_mat is None:
            return None
        for k in ((a, b), (b, a)):
            try:
                return float(blosum_mat[k])
            except Exception:  # noqa: BLE001
                continue
        return None

    human_ref_in = cassette_ref.get("human_reference") or {}
    discriminating = cassette_ref.get("discriminating") or {}
    ref_positions: Dict[str, List[tuple]] = {}
    for panel in ("IIIb", "IIIc"):
        rows = sorted((r for r in human_ref_in.get(panel, []) if r.get("i") is not None),
                      key=lambda r: r["i"])
        ref_positions[panel] = [(r["i"], (r.get("aa") or "").strip().upper()) for r in rows]

    def nongap(panel: str) -> List[tuple]:
        return [(i, aa) for (i, aa) in ref_positions[panel] if aa and aa not in ("-", ".", "*")]

    if not nongap("IIIb") and not nongap("IIIc"):
        return None

    truth_idx = {(r.get("species", ""), r.get("isoform", "")): r for r in _truth(run_dir)}
    species_nodes: Dict[str, Dict[str, Any]] = {}
    n_aligned = 0
    for panel in ("IIIb", "IIIc"):
        ref_ng = nongap(panel)
        if not ref_ng:
            continue
        ref_seq = "".join(aa for _, aa in ref_ng)
        msa = _read_fasta(run_dir / "MSA" / f"final_fgfr2_{panel}_cassette_msa.aln.faa")
        for header, aln_seq in msa.items():
            sp = header.split("|")[0].strip()
            if not sp or sp == "homo_sapiens":  # human never enters the analysed-species panel
                continue
            sp_seq = "".join(ch for ch in aln_seq.upper() if ch.isalpha())
            if not sp_seq:
                continue
            positions = _pairwise_positions(aligner, ref_seq, ref_ng, sp_seq, discriminating, blosum_fn)
            if not positions:
                continue
            n_aligned += 1
            t = truth_idx.get((sp, panel), {})
            idxs = [p["i"] for p in positions if p["i"] is not None]
            node = species_nodes.setdefault(sp, {
                "species": sp,
                "display_species_name": sp.replace("_", " ").title(),
                "panels": {},
            })
            node["panels"][panel] = {
                "isoform": panel,
                "available": True,
                "final_isoform_label": t.get("final_isoform_label", "") or panel,
                "validated_exon_type": t.get("validated_exon_type", ""),
                "final_claim_status_after_rescue": t.get("final_claim_status_after_rescue", ""),
                "claim_class": claim_class(t.get("final_claim_status_after_rescue", "")),
                "readiness_class": review.get((sp, panel), "neutral"),
                "is_review": review.get((sp, panel)) == "review",
                "source_table": rel(run_dir / "MSA" / f"final_fgfr2_{panel}_cassette_msa.aln.faa"),
                "positions": positions,
                "cassette_start": min(idxs) if idxs else None,
                "cassette_end": max(idxs) if idxs else None,
                "n_identical": sum(1 for p in positions if p["cls"] == "identical"),
                "n_diff": sum(1 for p in positions if p["cls"] in ("conservative", "nonconservative")),
                "n_gap": sum(1 for p in positions if p["cls"] == "gap"),
            }
    if n_aligned == 0:
        return None

    # pass the curated human reference rows through unchanged (keeps aa + property for tooltips)
    human_reference = {
        p: sorted((r for r in human_ref_in.get(p, []) if r.get("i") is not None), key=lambda r: r["i"])
        for p in ("IIIb", "IIIc")
    }
    out_species = sorted(species_nodes.values(), key=lambda x: x["display_species_name"])
    return {
        "human_reference": human_reference,
        "discriminating": discriminating,
        "species": out_species,
        "source_table": cassette_ref.get("source_table", ""),
    }


def build_coordinate_track_index(run_dir: Path) -> Dict[str, Any]:
    """Module 2 source: exon/CDS blocks + cassette span on protein AA coordinates (Figure 3C)."""
    rows = read_tsv(run_dir / "tables" / "figure3C_exon_to_protein_cassette_coordinate_map.tsv")
    review = _review_lookup(run_dir)
    by: Dict[tuple, Dict[str, Any]] = {}
    for r in rows:
        sp, iso = r.get("species", ""), r.get("isoform", "")
        if not sp or not iso:
            continue
        key = (sp, iso)
        cass_start = to_int(r.get("cassette_start_aa"))
        cass_end = to_int(r.get("cassette_end_aa"))
        node = by.setdefault(key, {
            "species": sp,
            "isoform": iso,
            "display_species_name": sp.replace("_", " ").title(),
            "final_isoform_label": r.get("final_isoform_label", ""),
            "protein_length": to_int(r.get("protein_length")),
            # never default an unresolved cassette to AA1: keep None when absent
            "cassette_start_aa": cass_start,
            "cassette_end_aa": cass_end,
            "cassette_available": cass_start is not None and cass_end is not None and cass_end > 0,
            "boundary_left_precision": r.get("boundary_left_precision", ""),
            "boundary_right_precision": r.get("boundary_right_precision", ""),
            "final_plot_status": r.get("final_plot_status", ""),
            "final_claim_status_after_rescue": r.get("final_claim_status_after_rescue", ""),
            "claim_class": claim_class(r.get("final_claim_status_after_rescue", "")),
            "readiness_class": review.get(key, "neutral"),
            "is_review": review.get(key) == "review",
            "blocks": [],
            "source_table": rel(run_dir / "tables" / "figure3C_exon_to_protein_cassette_coordinate_map.tsv"),
        })
        start = to_int(r.get("block_start_aa"))
        end = to_int(r.get("block_end_aa"))
        node["blocks"].append({
            "feature_type": r.get("feature_type", ""),
            "id": r.get("exon_or_cds_id", ""),
            "label": r.get("feature_label", ""),
            "start": start,
            "end": end,
            "is_iiib_cassette": truthy(r.get("is_IIIb_cassette")),
            "is_iiic_cassette": truthy(r.get("is_IIIc_cassette")),
            "in_cassette": (start is not None and end is not None and cass_start is not None and cass_end is not None
                            and end >= cass_start and start <= cass_end),
        })

    species: Dict[str, Dict[str, Any]] = {}
    for (sp, iso), node in by.items():
        node["blocks"].sort(key=lambda b: (b["start"] is None, b["start"]))
        s = species.setdefault(sp, {"species": sp,
                                    "display_species_name": node["display_species_name"],
                                    "panels": {}})
        s["panels"][iso] = node
    out = sorted(species.values(), key=lambda x: x["display_species_name"])
    return {"available": bool(rows),
            "domain_layer": "pending_interproscan",
            "availability": _availability_block(
                "exon map",
                run_dir / "tables" / "figure3C_exon_to_protein_cassette_coordinate_map.tsv",
                bool(rows)),
            "species": out}


def build_msa_index(run_dir: Path) -> Dict[str, Any]:
    """Module 3 source: parsed alignments + conservation summary + discriminating columns."""
    msa = run_dir / "MSA"
    files = {
        "full_length": ("Full-length FGFR2", msa / "final_fgfr2_full_length_protein_msa.aln.faa"),
        "iiib_cassette": ("IIIb cassette", msa / "final_fgfr2_IIIb_cassette_msa.aln.faa"),
        "iiic_cassette": ("IIIc cassette", msa / "final_fgfr2_IIIc_cassette_msa.aln.faa"),
        "combined_cassette": ("Combined IIIb/IIIc cassette", msa / "final_fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa"),
    }
    truth = _truth(run_dir)
    taxon = {r.get("species", ""): r.get("taxon_group", "") for r in truth}
    review = _review_lookup(run_dir)

    alignments: Dict[str, Any] = {}
    for key, (label, path) in files.items():
        parsed = _parse_alignment(path)
        if not parsed:
            alignments[key] = {"available": False, "label": label}
            continue
        for row in parsed["rows"]:
            sp = row["species"]
            row["display_species_name"] = sp.replace("_", " ").title()
            row["taxon_group"] = taxon.get(sp, "")
            row["is_human"] = sp == "homo_sapiens"
            row["is_review"] = review.get((sp, row["isoform"])) == "review"
        parsed["rows"].sort(key=lambda r: (r["isoform"], r["display_species_name"]))
        alignments[key] = {"available": True, "label": label, "file": rel(path),
                           "n_columns": parsed["n_columns"], "rows": parsed["rows"]}

    disc_rows = read_tsv(run_dir / "MSA" / "final_isoform_discriminating_residues.tsv")
    discriminating = []
    disc_cols_combined = []
    for r in disc_rows:
        if not truthy(r.get("informative_column")):
            continue
        is_disc = "discriminating" in (r.get("position_class", "") or "")
        col = to_int(r.get("combined_alignment_col"))
        entry = {
            "human_reference_residue_index": to_int(r.get("human_reference_residue_index")),
            "combined_alignment_col": col,
            "alignment_col": to_int(r.get("alignment_col")),
            "IIIb_major_aa": r.get("IIIb_major_aa", ""),
            "IIIc_major_aa": r.get("IIIc_major_aa", ""),
            "position_class": r.get("position_class", ""),
            "discriminating_score": to_float(r.get("discriminating_score")),
            "is_discriminating": is_disc,
        }
        discriminating.append(entry)
        if is_disc and col is not None:
            disc_cols_combined.append(col)

    conservation = []
    for r in read_tsv(run_dir / "MSA" / "final_full_length_msa_conservation_summary.tsv"):
        conservation.append({
            "species": r.get("species", ""),
            "isoform": r.get("isoform", ""),
            "region_type": r.get("region_type", ""),
            "mean_conservation_score": to_float(r.get("mean_conservation_score")),
            "mean_gap_fraction": to_float(r.get("mean_gap_fraction")),
            "conservation_status": r.get("conservation_status", ""),
        })

    return {
        "available": any(a.get("available") for a in alignments.values()),
        "availability": _availability_block(
            "alignment", msa / "final_fgfr2_full_length_protein_msa.aln.faa",
            any(a.get("available") for a in alignments.values())),
        "alignments": alignments,
        "discriminating": discriminating,
        "discriminating_columns_combined": sorted(set(disc_cols_combined)),
        "conservation": conservation,
        "tabs": [
            {"key": "full_length", "label": "Full-length FGFR2"},
            {"key": "iiib_cassette", "label": "IIIb cassette"},
            {"key": "iiic_cassette", "label": "IIIc cassette"},
            {"key": "combined_cassette", "label": "Combined cassette"},
            {"key": "discriminating", "label": "Discriminating residues"},
        ],
    }


SHARED_HUMAN_SYNTENY_REFERENCE = (
    PROJECT_ROOT / "references" / "synteny" / "human_fgfr2_10neighbor_reference.tsv")


def _shared_human_synteny_reference_node() -> Optional[Dict[str, Any]]:
    """Build a Homo sapiens FGFR2-neighborhood node from the curated shared reference.

    This is the SAME allowed human reference/control layer used elsewhere: it lets the
    "Compare to human" synteny view work for ANY custom run (where homo_sapiens is not an
    analysed species) without copying non-human Example results. Returns None if the shared
    reference file is unavailable, so the caller can signal "human reference not available"
    rather than render an empty view.
    """
    rows = read_tsv(SHARED_HUMAN_SYNTENY_REFERENCE)
    if not rows:
        return None

    neighbors: List[Dict[str, Any]] = []
    for r in rows:
        rank = to_int(r.get("human_neighbor_rank"))
        side = (r.get("human_neighbor_side", "") or "").strip()
        sym = (r.get("human_gene_symbol", "") or "").strip()
        if not rank or side not in ("upstream", "downstream") or not sym:
            continue
        neighbors.append(sc.neighbour_locus(
            side=side, rank=rank, source_symbol=sym, resolved_symbol=sym,
            gene_id=(r.get("human_gene_id", "") or "").strip(),
            protein_id=(r.get("human_protein_id", "") or "").strip(),
            strand=(r.get("human_strand", "") or "").strip(),
            orthology_class="curated", identity_status="human_reference_control",
            method="curated_human_reference"))
    if not neighbors:
        return None

    row = sc.species_row(
        "homo_sapiens", gene_symbol="FGFR2",
        target=sc.target_locus(gene_symbol="FGFR2", strand="-"),
        neighbours=neighbors, display_name="Homo sapiens",
        synteny_status="synteny_supported", taxon_group="Primates",
        is_human_reference_control=True,
        extra={"synteny_class": "reference_control",
               "synteny_status_class": "neutral",
               "has_resolved": True})
    return _with_legacy_keys(row)


def _with_legacy_keys(row: Dict[str, Any]) -> Dict[str, Any]:
    """Add the historical neighbour keys next to the canonical contract fields."""
    legacy = sc.legacy_nodes(row)
    row.update({
        "n_resolved": sum(1 for n in legacy if n["resolved"] and not n["is_anchor"]),
        "n_neighbors": row["displayed_flanking_count"],
        "neighbors5": legacy,
        "neighbors10": legacy,
    })
    row.setdefault("has_resolved", row["n_resolved"] > 0)
    return row


def build_synteny_locus_index(run_dir: Path) -> Dict[str, Any]:
    """Module 4 source: local FGFR2 gene neighborhood (resolved 5-neighbor + raw 10-neighbor)."""
    sdir = _synteny_source_dir(run_dir)
    if sdir is None:
        return {
            "schema_version": 3,
            "contract": "shared_synteny_v1",
            "available": False,
            "synteny_status": "not_computed",
            "synteny_reason": ("No local FGFR2 gene-neighborhood table was produced for this run; "
                               "synteny supports FGFR2 locus/orthology, not IIIb/IIIc identity."),
            "n_resolved_neighbors": 0,
            "n_flanking_loci": 0,
            "orthology_classes_present": [],
            "requested_neighbour_count": sc.REQUESTED_NEIGHBOUR_COUNT,
            "gene_symbol": "FGFR2",
            "target_symbol": "FGFR2",
            "extraction_warning": "",
            "has_10neighbor": False,
            "human_reference": None,
            "human_reference_available": False,
            "species": [],
            "source_tables": {},
        }
    fig9a = read_tsv(sdir / "tables" / "figure9A_fgfr2_local_gene_neighborhood_5neighbors.tsv")
    supp10 = read_tsv(sdir / "synteny" / "fgfr2_local_gene_neighborhood_10neighbors_supplement.tsv")
    truth = _truth(run_dir)
    taxon = {r.get("species", ""): r.get("taxon_group", "") for r in truth}
    synteny_truth = {r.get("species", ""): (r.get("combined_synteny_validation_class")
                                            or r.get("synteny_validation_class", "")) for r in truth}
    review = _review_lookup(run_dir)

    # Curated 5-neighbour panel: the resolved orthology assignment per slot.
    curated: Dict[str, Dict[Tuple[str, int], Dict[str, str]]] = {}
    synteny_class: Dict[str, str] = {}
    target_strand: Dict[str, str] = {}
    target_gene_id: Dict[str, str] = {}
    for r in fig9a:
        sp = r.get("species", "")
        if not sp:
            continue
        synteny_class.setdefault(sp, r.get("synteny_validation_class", ""))
        side, rank = r.get("side", ""), to_int(r.get("rank")) or 0
        if side == "fgfr2":
            target_strand.setdefault(sp, r.get("strand", ""))
            continue
        curated.setdefault(sp, {})[(side, rank)] = r

    # Raw 10-neighbour supplement: the genomic truth (gene ids, strands, spans).
    raw: Dict[str, List[Dict[str, str]]] = {}
    for r in supp10:
        sp = r.get("species", "")
        if not sp:
            continue
        target_strand.setdefault(sp, r.get("fgfr2_strand", ""))
        target_gene_id.setdefault(sp, r.get("fgfr2_gene_id", ""))
        raw.setdefault(sp, []).append(r)

    all_species = sorted(set(list(curated) + list(raw)))

    def _merged_neighbours(sp: str) -> List[Dict[str, Any]]:
        """One locus per real neighbour, raw genomics plus the curated symbol.

        The supplement holds the assembly facts and the curated panel holds the
        orthology assignment for the innermost five slots. Merging them means the
        displayed neighbourhood carries curated symbols and confidence classes
        instead of collapsing to raw identifiers.
        """
        out: List[Dict[str, Any]] = []
        seen: set = set()
        for r in raw.get(sp, []):
            side = r.get("neighbor_side", "")
            rank = to_int(r.get("neighbor_rank")) or 0
            if side not in ("upstream", "downstream") or not rank:
                continue
            seen.add((side, rank))
            cur = (curated.get(sp) or {}).get((side, rank), {})
            raw_symbol = (r.get("neighbor_symbol_raw", "") or "").strip()
            out.append(sc.neighbour_locus(
                side=side, rank=rank,
                source_symbol=raw_symbol,
                resolved_symbol=(cur.get("normalized_symbol", "") or "").strip(),
                gene_id=r.get("neighbor_gene_id", ""),
                protein_id=r.get("neighbor_protein_id", ""),
                strand=r.get("neighbor_strand", ""),
                orthology_class=(synteny_method_class(cur.get("identity_method", ""))
                                 if cur else ""),
                identity_status=cur.get("identity_status", ""),
                percent_identity=to_float(cur.get("broad_homology_percent_identity")),
                coverage=to_float(cur.get("broad_homology_query_coverage")),
                distance=to_int(r.get("distance_to_fgfr2")),
                seqid=r.get("seqid", ""),
                start=to_int(r.get("neighbor_start")), end=to_int(r.get("neighbor_end")),
                method=cur.get("identity_method", "") or "raw_annotation_only"))
        # A curated slot with no supplement row is still a real locus.
        for (side, rank), cur in sorted((curated.get(sp) or {}).items()):
            if (side, rank) in seen or side not in ("upstream", "downstream"):
                continue
            out.append(sc.neighbour_locus(
                side=side, rank=rank,
                source_symbol=(cur.get("raw_symbol", "") or "").strip(),
                resolved_symbol=(cur.get("normalized_symbol", "") or "").strip(),
                strand=cur.get("strand", ""),
                orthology_class=synteny_method_class(cur.get("identity_method", "")),
                identity_status=cur.get("identity_status", ""),
                percent_identity=to_float(cur.get("broad_homology_percent_identity")),
                coverage=to_float(cur.get("broad_homology_query_coverage")),
                method=cur.get("identity_method", "")))
        return out

    def assemble(sp: str) -> Dict[str, Any]:
        neighbours = _merged_neighbours(sp)
        first_raw = next(iter(raw.get(sp, [])), {})
        cls = synteny_class.get(sp) or synteny_truth.get(sp, "")
        row = sc.species_row(
            sp, gene_symbol="FGFR2",
            target=sc.target_locus(
                gene_symbol="FGFR2", gene_id=target_gene_id.get(sp, ""),
                strand=target_strand.get(sp, ""), seqid=first_raw.get("seqid", ""),
                start=to_int(first_raw.get("fgfr2_start")),
                end=to_int(first_raw.get("fgfr2_end")),
                protein_id=first_raw.get("fgfr2_protein_id", "")),
            neighbours=neighbours,
            display_name=sc.display_binomial(sp),
            taxon_group=taxon.get(sp, ""),
            synteny_status=cls or "local_neighbourhood",
            is_review=(review.get((sp, "IIIb")) == "review"
                       or review.get((sp, "IIIc")) == "review"),
            extra={"synteny_class": cls,
                   "synteny_status_class": readiness_class(cls),
                   "has_resolved": sp in curated})
        return _with_legacy_keys(row)

    # The canonical taxonomic order, shared with every other comparative view.
    # Sorting by (taxon_group, name) put Amphibians before Birds before Other
    # mammals — alphabetical group names, not a taxonomic sequence.
    species = [assemble(sp) for sp in _species_order.order_species(all_species)]

    # Count neighbours that actually resolved to a gene. Small custom runs
    # frequently produce a synteny table where every neighbour is `missing`
    # because the local genome-neighbourhood annotation was not available —
    # that is "not computed", not a failure.
    totals = sc.summarise(species)
    n_resolved = totals["n_resolved_neighbours"]
    extraction_statuses = {(r.get("extraction_status", "") or "").strip() for r in supp10}
    extraction_warnings = [w for w in ((r.get("extraction_warning", "") or "").strip() for r in supp10) if w]
    # exact run-local raw input the synteny extractor needs (genome-wide NCBI Datasets GFF).
    # Surfaced in the UI Technical details when synteny could not be computed. This is a
    # RAW INPUT path, never an example result — custom runs must fetch their own annotation.
    source_files = [s for s in ((r.get("source_file", "") or "").strip() for r in supp10) if s]
    models_cache = rel(run_dir.parent / "02_models" / "_ncbi_datasets_cache")
    synteny_missing_source = "" if source_files else (
        f"{models_cache}/ncbi_<taxid>/<assembly_accession>/unzipped/ncbi_dataset/data/"
        "<assembly_accession>/genomic.gff")

    if not species:
        synteny_status = "not_computed"
        synteny_reason = "No local FGFR2 gene-neighborhood table was produced for this run."
    elif n_resolved > 0:
        synteny_status = "available"
        synteny_reason = ""
    elif any("unavailable" in s or "missing" in s for s in extraction_statuses) or extraction_warnings:
        synteny_status = "not_computed"
        synteny_reason = ("Synteny not computed because the required genomic annotation could not be "
                          "obtained for this run (the FGFR2 genomic neighborhood GFF was not available "
                          "in this run's local genome cache), so synteny neighbors could not be "
                          "resolved. The run remains usable: IIIb/IIIc identity rests on sequence and "
                          "marker evidence; synteny supports FGFR2 locus/orthology, not isoform identity.")
    else:
        synteny_status = "not_applicable"
        synteny_reason = ("Synteny was not resolved for this custom run. Synteny supports FGFR2 "
                          "locus/orthology and is not required for IIIb/IIIc identity, which is "
                          "based on sequence and marker evidence.")

    # Human reference/control row for the "Compare to human" view. If homo_sapiens is an
    # analysed species (Example dataset) use its own run-local node; otherwise fall back to
    # the curated shared human FGFR2 neighborhood (the only allowed reusable reference layer),
    # so custom runs can compare-to-human without copying non-human Example results.
    human_ref = next((s for s in species if s["species"] == "homo_sapiens"), None)
    human_ref_role = "analysed_species" if human_ref else ""
    if human_ref is None:
        human_ref = _shared_human_synteny_reference_node()
        if human_ref is not None:
            human_ref_role = "human_reference_control"

    for row in species:
        row["comparison_available"] = (human_ref is not None
                                       and row["species_id"] != "homo_sapiens")

    return {
        "schema_version": 3,
        "contract": "shared_synteny_v1",
        "available": n_resolved > 0,
        "synteny_status": synteny_status,
        "synteny_reason": synteny_reason,
        "n_resolved_neighbors": n_resolved,
        "n_flanking_loci": totals["n_flanking_loci"],
        "orthology_classes_present": totals["classes_present"],
        "requested_neighbour_count": sc.REQUESTED_NEIGHBOUR_COUNT,
        "gene_symbol": "FGFR2",
        "target_symbol": "FGFR2",
        "scope": ("comparative_synteny" if len(species) > 1
                  else "single_species_local_neighbourhood"),
        "extraction_warning": extraction_warnings[0] if extraction_warnings else "",
        "synteny_missing_source": synteny_missing_source,
        "synteny_source_file": source_files[0] if source_files else "",
        "has_10neighbor": bool(raw),
        "human_reference": human_ref,
        "human_reference_available": human_ref is not None,
        "human_reference_role": human_ref_role,
        "species_order": _species_order.build_species_order(
            [s["species_id"] for s in species]),
        "species": species,
        "source_tables": {
            "resolved_5neighbor": rel(sdir / "tables" / "figure9A_fgfr2_local_gene_neighborhood_5neighbors.tsv"),
            "raw_10neighbor": rel(sdir / "synteny" / "fgfr2_local_gene_neighborhood_10neighbors_supplement.tsv"),
        },
    }


def build_species_story_index(run_dir: Path) -> Dict[str, Any]:
    """Module 5 source: per species/isoform evidence story (truth table + review + projection)."""
    truth = _truth(run_dir)
    review_expl = {(r.get("species", ""), r.get("isoform", "")): r
                   for r in read_tsv(run_dir / "tables" / "final_review_case_explanation.tsv")}
    truth_rel = rel(run_dir / "final_pre_interpro_truth_table.tsv")
    proj_rel = rel(run_dir / "MSA" / "final_cassette_msa_boundary_projection.tsv")
    review_rel = rel(run_dir / "tables" / "final_review_case_explanation.tsv")
    manifest_rel = rel(run_dir / "freeze" / "final_pre_interpro_sequence_manifest.tsv")

    species: Dict[str, Dict[str, Any]] = {}
    for r in truth:
        sp, iso = r.get("species", ""), r.get("isoform", "")
        if not sp or not iso:
            continue
        readiness = readiness_class(r.get("pre_interpro_readiness_class", ""))
        is_review = readiness == "review"
        expl = review_expl.get((sp, iso), {})

        def step(key, title, status, text, source, ids=None, figure=None):
            return {"key": key, "title": title, "class": status, "text": text,
                    "source_table": source, "ids": ids or {}, "figure": figure}

        upstream = r.get("upstream_label", "")
        validated = r.get("validated_exon_type", "") or r.get("final_isoform_label", "")
        consistency = r.get("label_consistency_status", "")
        rescue_dec = r.get("rescue_decision", "")
        rescue_required = truthy(r.get("rescue_required"))
        final_src = r.get("final_label_source", "")

        steps = [
            step("input", "Input / upstream label", "neutral",
                 f"Source annotation labelled this isoform {upstream or '—'}; sequence-calibrated type is {validated or '—'}.",
                 truth_rel,
                 {"upstream_label": upstream, "validated_exon_type": validated}),
            step("reconciliation", "Label reconciliation",
                 "minor" if "swapped" in consistency.lower() else statusify(consistency),
                 reconciliation_text(consistency, upstream, r.get("final_isoform_label", "")),
                 truth_rel, {"final_isoform_label": r.get("final_isoform_label", "")},
                 "Figure 4"),
            step("rescue", "Rescue / provenance",
                 rescue_step_class(rescue_dec, final_src, is_review),
                 rescue_text(rescue_dec, final_src, rescue_required, expl),
                 review_rel if expl else truth_rel,
                 {"final_label_source": final_src}),
            step("coordinates", "Coordinate validation",
                 statusify(r.get("coordinate_validation_status", "")),
                 f"Coordinates: {r.get('coordinate_validation_status','—')}; protein integrity: {r.get('protein_integrity_status','—')}.",
                 truth_rel,
                 {"protein_length": r.get("protein_length", ""),
                  "transcript_id": r.get("transcript_id", ""),
                  "protein_id": r.get("protein_id", "")},
                 "Figure 3C"),
            step("msa", "MSA support",
                 statusify(r.get("MSA_full_length_status", "")),
                 f"Full-length MSA: {r.get('MSA_full_length_status','—')}; cassette MSA: {r.get('MSA_cassette_status','—')}; boundary: {r.get('boundary_robustness_class','—')}.",
                 proj_rel, {}, "Figure 5 / 6"),
            step("synteny", "Synteny support",
                 statusify(r.get("combined_synteny_validation_class") or r.get("synteny_validation_class", "")),
                 f"Synteny: {r.get('combined_synteny_validation_class') or r.get('synteny_validation_class','—')}; orthology: {r.get('orthology_status','—')}.",
                 truth_rel, {}, "Figure 9A"),
            step("readiness", "Final pre-InterPro readiness",
                 readiness,
                 (expl.get("final_interpretation") if is_review and expl
                  else f"Final status: {r.get('pre_interpro_readiness_class','—')}."),
                 truth_rel,
                 {"pre_interpro_readiness_class": r.get("pre_interpro_readiness_class", "")}),
            step("files", "Files / exports", "neutral",
                 "Included in primary freeze FASTA + manifest." if not is_review
                 else "Kept in review/supplement FASTA with provenance; excluded from the primary set.",
                 manifest_rel,
                 {"sequence_md5": r.get("sequence_md5", "")}),
        ]

        node = species.setdefault(sp, {
            "species": sp,
            "display_species_name": r.get("display_species_name", "") or sp.replace("_", " ").title(),
            "taxon_group": r.get("taxon_group", ""),
            "panels": {},
        })
        node["panels"][iso] = {
            "isoform": iso,
            "overall": readiness,
            "is_review": is_review,
            "final_claim_status_after_rescue": r.get("final_claim_status_after_rescue", ""),
            "final_isoform_label": r.get("final_isoform_label", ""),
            "rescued": ("rescued" in rescue_dec.lower() or "rescued" in final_src.lower()),
            "review_explanation": expl.get("final_interpretation", ""),
            "ids": {
                "transcript_id": r.get("transcript_id", ""),
                "protein_id": r.get("protein_id", ""),
                "gene_id": r.get("gene_id", ""),
                "protein_length": to_int(r.get("protein_length")),
                "sequence_md5": r.get("sequence_md5", ""),
            },
            "steps": steps,
        }

    out = sorted(species.values(), key=lambda x: (x["taxon_group"], x["display_species_name"]))
    return {"available": bool(truth), "species": out}


def statusify(value: str) -> str:
    """Map an arbitrary status token to a colour class (re-uses readiness_class + token_class)."""
    rc = readiness_class(value)
    if rc not in ("neutral", "unknown"):
        return rc
    tc = token_class(value)
    return tc if tc not in ("neutral", "unknown") else \
        ("accepted" if any(k in (value or "").lower()
                           for k in ("pass", "confirmed", "supported", "strong", "high", "robust", "ortholog"))
         else "neutral")


def reconciliation_text(consistency: str, upstream: str, final: str) -> str:
    c = (consistency or "").lower()
    if "swapped" in c:
        return f"Sequence reconciliation re-assigned the label relative to the source ({upstream or '—'} → {final or '—'})."
    if "consistent" in c or "agree" in c:
        return "Source label and sequence-calibrated label agree."
    return f"Label consistency: {consistency or '—'}."


def rescue_step_class(rescue_dec: str, final_src: str, is_review: bool) -> str:
    d = (rescue_dec + " " + final_src).lower()
    if "rescued" in d or "external" in d and "validated" in d:
        return "accepted"
    if is_review or "unresolved" in d or "supplement" in d or "sequence_support_only" in d:
        return "review"
    if "confirmed" in d:
        return "accepted"
    return "neutral"


def rescue_text(rescue_dec: str, final_src: str, rescue_required: bool, expl: Dict[str, str]) -> str:
    d = (rescue_dec + " " + final_src).lower()
    if "rescued" in d:
        return "Rescued with an external source-compatible validated candidate; accepted as a primary row with provenance retained."
    if "sequence_support_only" in d or "unresolved" in d or "supplement" in (expl.get("final_claim_status_after_rescue", "") or "").lower():
        return "No source-compatible externally validated candidate found (sequence support only); kept as supplement/review with provenance."
    if "confirmed" in d:
        return "Current candidate confirmed after exhaustive screen / sequence reconciliation; no external rescue required."
    return rescue_dec or "—"


# --------------------------------------------------------------------------- #
# orchestration
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Post-InterPro / pyTMHMM domain-architecture indices (step 15)
#
# These are derived from the fixed post-InterProScan analysis folder
# (15_exon_domain_boundary_post_interpro), a *sibling* of the closure run dir.
# The post-InterPro step is not part of the interactive web pipeline, so these
# indices are only populated for a run that has the folder next to it (the
# example freeze); for fresh web runs they degrade to {"available": False}.
#
# Hard rules honoured here:
#   * IIIb/IIIc labels ALWAYS come from the final truth table
#     (final_isoform_label); InterProScan and pyTMHMM never relabel a cassette.
#   * pyTMHMM is the transmembrane layer because InterProScan does not annotate
#     TM helices for these proteins.
#   * The 3 review_unusual_domain_order cases were audited
#     (post_interpro_qc_review_case_audit.tsv) and shown to be low-confidence
#     coordinate artifacts, not biological anomalies: their web *display* status
#     is softened to a minor flag with an explicit note, but the raw QC status
#     and primary/review membership are never mutated.
# --------------------------------------------------------------------------- #
_ARCH_STEP_NAME = "15_exon_domain_boundary_post_interpro"


def _post_interpro_dir(run_dir: Path) -> Optional[Path]:
    """Locate the sibling step-15 post-InterPro/pyTMHMM architecture folder."""
    parent = Path(run_dir).parent
    cand = parent / _ARCH_STEP_NAME
    if (cand / "tables" / "fgfr2_domain_architecture_qc.tsv").exists():
        return cand
    for c in sorted(parent.glob("15_*")):
        if (c / "tables" / "fgfr2_domain_architecture_qc.tsv").exists():
            return c
    return None


_EXON_NUM_RE = re.compile(r"exon\s+(\d+)", re.IGNORECASE)


def _exon_number(label: str) -> Optional[int]:
    m = _EXON_NUM_RE.search(label or "")
    return int(m.group(1)) if m else None


def _arch_qc_class(status: str) -> str:
    """final_qc_status -> {accepted, minor, review}."""
    s = (status or "").lower()
    if s == "architecture_supported":
        return "accepted"
    if "review" in s:
        return "review"
    if "minor" in s or "partial" in s or "supported" in s:
        return "minor"
    return "neutral"


def _split_warnings(raw: str) -> List[str]:
    return [w.strip() for w in (raw or "").split(";") if w.strip()]


def _arch_audit_lookup(arch: Path) -> Dict[tuple, Dict[str, str]]:
    """(species, isoform) -> review-case audit row, if the audit table exists."""
    rows = read_tsv(arch / "tables" / "post_interpro_qc_review_case_audit.tsv")
    out: Dict[tuple, Dict[str, str]] = {}
    for r in rows:
        out[(r.get("species", ""), r.get("isoform", ""))] = r
    return out


def _arch_recon_lookup(arch: Path) -> Dict[tuple, Dict[str, str]]:
    """(species, isoform) -> exon-block reconstruction audit row, if present."""
    rows = read_tsv(arch / "tables" / "exon_block_coordinate_reconstruction_audit.tsv")
    out: Dict[tuple, Dict[str, str]] = {}
    for r in rows:
        # audit isoform column carries the display label (e.g. "FGFR2 IIIb");
        # normalise back to the pipeline key.
        iso = r.get("isoform", "")
        iso_key = "IIIb" if "IIIb" in iso else ("IIIc" if "IIIc" in iso else iso)
        out[(r.get("species", ""), iso_key)] = r
    return out


_RECON_NOTES = {
    "native_exon_blocks_reconstructed": (
        "Coding-exon blocks were rebuilt from the native local CDS coordinates of this "
        "exact transcript (they previously came from a different resolver transcript). "
        "The cassette slot uses the species-specific validated reference coordinate, "
        "which maps upstream of the TM and kinase. IIIb/IIIc label and primary/review "
        "membership are unchanged."),
    "cassette_only_high_confidence": (
        "Native CDS coordinates for this final RefSeq protein are not available locally, "
        "so the earlier (misleading) template exon blocks are hidden and only the "
        "validated cassette slot is shown. The cassette maps upstream of the TM and "
        "kinase. IIIb/IIIc label and primary/review membership are unchanged."),
    "minor_length_clamped": (
        "The final coding-exon block extended 1-2 aa past the protein length "
        "(codon-boundary rounding) and was clamped to the protein length for display. "
        "Coordinates are otherwise unchanged; domain architecture and IIIb/IIIc label "
        "are unchanged."),
    "exon_blocks_hidden_untrusted": (
        "Coding-exon block coordinates could not be reconstructed from native CDS data "
        "for this transcript, so the untrusted blocks are hidden. Domains, TM helix, "
        "kinase and (where available) the cassette slot remain shown. IIIb/IIIc label "
        "and primary/review membership are unchanged."),
}


def _arch_display_status(final_status: str, audit: Optional[Dict[str, str]],
                         exon_block_status: str = "") -> tuple:
    """(display_status, display_class, display_note).

    The exon-block reconstruction status (from the generator) is the primary source
    of the explanatory note. The legacy review-case audit softening is retained as a
    fallback for any run whose raw QC still carries review_unusual_domain_order.
    """
    note = _RECON_NOTES.get(exon_block_status or "", "")
    if final_status == "review_unusual_domain_order":
        soft_issues = {"fallback_coordinate_used", "coordinate_join_error"}
        if audit and audit.get("likely_issue", "") in soft_issues:
            note = note or (
                "Flagged only by a low-confidence exon-block cassette coordinate "
                "(native codon-phase offset); the species-specific reference cassette "
                "maps upstream of the TM/kinase as expected. IIIb/IIIc label and "
                "primary/review membership are unchanged.")
            return "architecture_supported_with_minor_flags", "minor", note
    return final_status, _arch_qc_class(final_status), note


def _load_arch_features(arch: Path) -> Dict[tuple, Dict[str, Any]]:
    """Assemble per-(species, isoform) domains / TM / exons / cassette / warnings."""
    rows = read_tsv(arch / "tables" / "exon_domain_architecture_features.tsv")
    feat_src = rel(arch / "tables" / "exon_domain_architecture_features.tsv")
    by: Dict[tuple, Dict[str, Any]] = {}
    for r in rows:
        sp, iso = r.get("species", ""), r.get("isoform", "")
        if not sp or not iso:
            continue
        key = (sp, iso)
        node = by.setdefault(key, {
            "species": sp,
            "isoform": iso,
            "display_species_name": sp.replace("_", " ").title(),
            "transcript_id": r.get("transcript_id", ""),
            "protein_id": r.get("protein_id", ""),
            "protein_length": to_int(r.get("protein_length")),
            "domains": [], "tm": [], "exons": [], "cassette": None,
            "warnings": [],
            "source_table": feat_src,
        })
        ftype = r.get("feature_type", "")
        start, end = to_int(r.get("start_aa")), to_int(r.get("end_aa"))
        label = r.get("feature_label", "")
        source = r.get("source", "")
        status = r.get("status", "")
        if ftype == "transmembrane_pytmhmm":
            node["tm"].append({"label": label, "start": start, "end": end,
                               "status": status, "source": source})
        elif ftype in ("ig_like_domain", "kinase_domain", "other_domain",
                       "signal_peptide"):
            node["domains"].append({"class": ftype, "label": label,
                                    "start": start, "end": end, "source": source})
        elif ftype == "coding_exon":
            node["exons"].append({"number": _exon_number(label), "label": label,
                                  "start": start, "end": end, "is_cassette": False})
        elif ftype in ("IIIb_slot", "IIIc_slot"):
            num = _exon_number(label)
            node["exons"].append({"number": num, "label": label,
                                  "start": start, "end": end, "is_cassette": True})
            node["cassette"] = {"slot_type": ftype, "label": label,
                                "number": num, "start": start, "end": end}
        elif ftype == "warning":
            node["warnings"].append(label)
    for node in by.values():
        node["exons"].sort(key=lambda e: (e["start"] is None, e["start"]))
        node["domains"].sort(key=lambda d: (d["start"] is None, d["start"]))
        node["tm"].sort(key=lambda t: (t["start"] is None, t["start"]))
    return by


def _load_arch_qc(arch: Path) -> Dict[tuple, Dict[str, str]]:
    rows = read_tsv(arch / "tables" / "fgfr2_domain_architecture_qc.tsv")
    return {(r.get("species", ""), r.get("isoform", "")): r for r in rows}


def _arch_figure_paths(arch: Path, species: str, isoform: str) -> Dict[str, Any]:
    stem = f"{species}_{isoform}_exon_domain_architecture"
    base = arch / "figures" / "per_species"
    out: Dict[str, Any] = {}
    for ext in ("png", "pdf", "svg"):
        p = base / f"{stem}.{ext}"
        if p.exists():
            out[ext] = rel(p)
    return out


def build_species_domain_architecture(run_dir: Path) -> Dict[str, Any]:
    """Per (species, isoform) protein domain architecture for the interactive tab."""
    arch = _post_interpro_dir(run_dir)
    if not arch:
        return {"available": False, "tm_layer": "pyTMHMM",
                "reason": "post-InterPro architecture folder not present for this run.",
                "species": []}
    feats = _load_arch_features(arch)
    qc = _load_arch_qc(arch)
    audit = _arch_audit_lookup(arch)
    recon = _arch_recon_lookup(arch)
    truth = {(r.get("species", ""), r.get("isoform", "")): r for r in _truth(run_dir)}

    species: Dict[str, Dict[str, Any]] = {}
    for key, node in feats.items():
        sp, iso = key
        q = qc.get(key, {})
        au = audit.get(key)
        rc = recon.get(key)
        exon_block_status = q.get("exon_block_display_status", "")
        final_status = q.get("final_qc_status", "")
        disp_status, disp_class, disp_note = _arch_display_status(
            final_status, au, exon_block_status)
        tr = truth.get(key, {})
        entry = {
            **node,
            "final_isoform_label": tr.get("final_isoform_label", q.get("isoform", iso)),
            "axis": {"start": 1, "end": node.get("protein_length") or 0},
            "qc": {
                "final_qc_status": final_status,
                "final_qc_class": _arch_qc_class(final_status),
                "display_qc_status": disp_status,
                "display_qc_class": disp_class,
                "display_note": disp_note,
                "exon_block_display_status": exon_block_status,
                "domain_order_status": q.get("domain_order_status", ""),
                "cassette_slot_position_status": q.get("cassette_slot_position_status", ""),
                "exon_domain_mapping_status": q.get("exon_domain_mapping_status", ""),
                "tm_agreement": q.get("tm_agreement", ""),
                "pytmhmm_tm_found": truthy(q.get("pytmhmm_tm_found")),
                "interpro_tm_found": truthy(q.get("interpro_tm_found")),
                "kinase_found": truthy(q.get("kinase_found")),
                "signal_region_supported": truthy(q.get("signal_region_supported")),
                "expected_ig_like_domain_count": to_int(q.get("expected_ig_like_domain_count")),
                "warnings": _split_warnings(q.get("warnings", "")),
            },
            "audited": bool(au) or bool(rc),
            "audit": None,
            "reconstruction": None,
            "figures": _arch_figure_paths(arch, sp, iso),
        }
        if rc:
            entry["reconstruction"] = {
                "final_display_status": rc.get("final_display_status", exon_block_status),
                "reconstruction_success": rc.get("reconstruction_success", ""),
                "number_of_exon_blocks": to_int(rc.get("number_of_exon_blocks")),
                "max_exon_end_aa": to_int(rc.get("max_exon_end_aa")),
                "exceeds_protein_length": rc.get("exceeds_protein_length", ""),
                "reconstructed_exon_block_source": rc.get("reconstructed_exon_block_source", ""),
                "notes": rc.get("notes", ""),
            }
        if au:
            entry["audit"] = {
                "likely_issue": au.get("likely_issue", ""),
                "overlaps_kinase": au.get("overlaps_kinase", ""),
                "action": au.get("action", ""),
                "final_interpretation": au.get("final_interpretation", ""),
                "cassette_start_aa_used": to_int(au.get("cassette_start_aa_used")),
                "cassette_end_aa_used": to_int(au.get("cassette_end_aa_used")),
                "kinase_start_aa": to_int(au.get("kinase_start_aa")),
                "kinase_end_aa": to_int(au.get("kinase_end_aa")),
                "source_coordinate_file": au.get("source_coordinate_file", ""),
            }
        s = species.setdefault(sp, {"species": sp,
                                    "display_species_name": node["display_species_name"],
                                    "panels": {}})
        s["panels"][iso] = entry

    out = sorted(species.values(), key=lambda x: x["display_species_name"])
    return {"available": True, "tm_layer": "pyTMHMM",
            "label_source": "final_pre_interpro_truth_table.tsv",
            "species": out}


def build_domain_architecture_qc(run_dir: Path) -> Dict[str, Any]:
    """Flat per-protein QC view (raw status + audited display status)."""
    arch = _post_interpro_dir(run_dir)
    if not arch:
        return {"available": False, "rows": []}
    qc = _load_arch_qc(arch)
    audit = _arch_audit_lookup(arch)
    recon = _arch_recon_lookup(arch)
    rows_out: List[Dict[str, Any]] = []
    for (sp, iso), q in sorted(qc.items()):
        au = audit.get((sp, iso))
        rc = recon.get((sp, iso))
        exon_block_status = q.get("exon_block_display_status", "")
        final_status = q.get("final_qc_status", "")
        disp_status, disp_class, disp_note = _arch_display_status(
            final_status, au, exon_block_status)
        rows_out.append({
            "species": sp,
            "display_species_name": sp.replace("_", " ").title(),
            "isoform": iso,
            "protein_id": q.get("protein_id", ""),
            "final_qc_status": final_status,
            "final_qc_class": _arch_qc_class(final_status),
            "display_qc_status": disp_status,
            "display_qc_class": disp_class,
            "display_note": disp_note,
            "exon_block_display_status": exon_block_status,
            "domain_order_status": q.get("domain_order_status", ""),
            "cassette_slot_position_status": q.get("cassette_slot_position_status", ""),
            "exon_domain_mapping_status": q.get("exon_domain_mapping_status", ""),
            "tm_agreement": q.get("tm_agreement", ""),
            "receptor_tm_start_aa": to_int(q.get("receptor_tm_start_aa")),
            "receptor_tm_end_aa": to_int(q.get("receptor_tm_end_aa")),
            "pytmhmm_tm_found": truthy(q.get("pytmhmm_tm_found")),
            "interpro_tm_found": truthy(q.get("interpro_tm_found")),
            "kinase_found": truthy(q.get("kinase_found")),
            "signal_region_supported": truthy(q.get("signal_region_supported")),
            "warnings": _split_warnings(q.get("warnings", "")),
            "audited": bool(au) or bool(rc),
            "likely_issue": au.get("likely_issue", "") if au else "",
            "final_interpretation": au.get("final_interpretation", "") if au else "",
            "reconstruction_success": rc.get("reconstruction_success", "") if rc else "",
            "reconstruction_notes": rc.get("notes", "") if rc else "",
        })
    recon_tbl = arch / "tables" / "exon_block_coordinate_reconstruction_audit.tsv"
    return {
        "available": True,
        "source_table": rel(arch / "tables" / "fgfr2_domain_architecture_qc.tsv"),
        "audit_table": rel(arch / "tables" / "post_interpro_qc_review_case_audit.tsv")
        if (arch / "tables" / "post_interpro_qc_review_case_audit.tsv").exists() else None,
        "reconstruction_audit_table": rel(recon_tbl) if recon_tbl.exists() else None,
        "rows": rows_out,
    }


def build_domain_architecture_summary(run_dir: Path) -> Dict[str, Any]:
    """Compact counts for the Overview 'Post-InterPro architecture' card."""
    arch = _post_interpro_dir(run_dir)
    if not arch:
        return {"available": False}
    qc = _load_arch_qc(arch)
    tm_rows = read_tsv(arch / "tables" / "pytmhmm_tm_features_normalized.tsv")
    interpro_rows = read_tsv(arch / "tables" / "interpro_domain_features_normalized.tsv")
    _audit = _arch_audit_lookup(arch)

    status_counts: Dict[str, int] = {}
    kinase = tm_pred = sp_supported = 0
    proteins = set()
    for (sp, iso), q in qc.items():
        proteins.add((sp, iso))
        st = q.get("final_qc_status", "")
        status_counts[st] = status_counts.get(st, 0) + 1
        if truthy(q.get("kinase_found")):
            kinase += 1
        if truthy(q.get("pytmhmm_tm_found")):
            tm_pred += 1
        if truthy(q.get("signal_region_supported")):
            sp_supported += 1

    with_interpro = len({(r.get("species"), r.get("isoform")) for r in interpro_rows})
    with_tm = len({(r.get("species"), r.get("isoform")) for r in tm_rows
                   if r.get("status") == "receptor_tm"})

    # coordinate-artifact cases: previously review_unusual_domain_order, now
    # resolved via native exon-block reconstruction / validated cassette slot.
    recon_rows = read_tsv(arch / "tables" / "exon_block_coordinate_reconstruction_audit.tsv")
    resolved_cases = [{
        "species": r.get("species", ""),
        "display_species_name": r.get("species", "").replace("_", " ").title(),
        "isoform": r.get("isoform", ""),
        "final_display_status": r.get("final_display_status", ""),
        "reconstruction_success": r.get("reconstruction_success", ""),
        "exceeds_protein_length": r.get("exceeds_protein_length", ""),
        "final_interpretation": r.get("notes", ""),
    } for r in recon_rows]

    return {
        "available": True,
        "generated_at": now_iso(),
        "step": _ARCH_STEP_NAME,
        "tm_layer": "pyTMHMM (InterProScan did not annotate transmembrane helices)",
        "proteins_annotated": len(proteins),
        "with_interpro_hits": with_interpro,
        "with_kinase": kinase,
        "with_tm_pytmhmm": with_tm if with_tm else tm_pred,
        "with_signal_region_supported": sp_supported,
        "qc_status_counts": status_counts,
        "supported": status_counts.get("architecture_supported", 0),
        "minor_flags": status_counts.get("architecture_supported_with_minor_flags", 0),
        "review_warnings": status_counts.get("review_unusual_domain_order", 0),
        # coordinate-artifact cases resolved by exon-block reconstruction
        "coordinate_artifact_cases_resolved": len(resolved_cases),
        "review_cases_audited": len(resolved_cases),
        "review_all_coordinate_artifacts": bool(resolved_cases),
        "review_cases": resolved_cases,
        "resolved_coordinate_artifacts": resolved_cases,
        "report": rel(arch / "reports" / "post_interpro_exon_domain_architecture_summary.md")
        if (arch / "reports" / "post_interpro_exon_domain_architecture_summary.md").exists() else None,
    }


def build_domain_architecture_index(run_dir: Path) -> Dict[str, Any]:
    """Top-level catalogue for the Domain architecture tab + Figure Gallery group."""
    arch = _post_interpro_dir(run_dir)
    if not arch:
        return {"available": False, "species": [], "overview_figures": [], "tables": []}
    qc = _load_arch_qc(arch)
    audit = _arch_audit_lookup(arch)
    truth = {(r.get("species", ""), r.get("isoform", "")): r for r in _truth(run_dir)}
    taxon = {r.get("species", ""): r.get("taxon_group", "") for r in _truth(run_dir)}

    species: Dict[str, Dict[str, Any]] = {}
    for (sp, iso), q in qc.items():
        au = audit.get((sp, iso))
        final_status = q.get("final_qc_status", "")
        disp_status, disp_class, _n = _arch_display_status(
            final_status, au, q.get("exon_block_display_status", ""))
        node = species.setdefault(sp, {
            "species": sp,
            "display_species_name": sp.replace("_", " ").title(),
            "taxon_group": taxon.get(sp, ""),
            "isoforms": {},
        })
        node["isoforms"][iso] = {
            "final_isoform_label": truth.get((sp, iso), {}).get("final_isoform_label", iso),
            "protein_id": q.get("protein_id", ""),
            "final_qc_status": final_status,
            "display_qc_status": disp_status,
            "display_qc_class": disp_class,
            "audited": bool(au),
            "figure_png": (rel(arch / "figures" / "per_species" /
                               f"{sp}_{iso}_exon_domain_architecture.png")
                           if (arch / "figures" / "per_species" /
                               f"{sp}_{iso}_exon_domain_architecture.png").exists() else None),
        }

    ov = arch / "figures" / "overview"
    overview_defs = [
        ("Figure_10_all_species_FGFR2_exon_domain_architecture_primary",
         "All species — FGFR2 exon–domain architecture", "all_species"),
        ("Figure_10A_IIIb_exon_domain_architecture_primary",
         "IIIb isoforms — exon–domain architecture", "iiib"),
        ("Figure_10B_IIIc_exon_domain_architecture_primary",
         "IIIc isoforms — exon–domain architecture", "iiic"),
        ("Figure_10C_mammals_exon_domain_architecture_primary",
         "Mammals — exon–domain architecture", "mammals"),
        ("Figure_10D_nonmammals_exon_domain_architecture_primary",
         "Non-mammals — exon–domain architecture", "nonmammals"),
    ]
    overview_figures = []
    for stem, title, group in overview_defs:
        entry: Dict[str, Any] = {"id": stem, "title": title, "group": group}
        present = False
        for ext in ("png", "pdf", "svg"):
            p = ov / f"{stem}.{ext}"
            if p.exists():
                entry[ext] = rel(p)
                present = True
        if present:
            overview_figures.append(entry)

    def _tbl(name: str, desc: str) -> Optional[Dict[str, Any]]:
        p = arch / "tables" / name
        return {"name": name, "path": rel(p), "description": desc} if p.exists() else None

    tables = [t for t in [
        _tbl("interpro_domain_features_normalized.tsv",
             "InterProScan domain matches (normalized, AA coordinates)."),
        _tbl("pytmhmm_tm_features_normalized.tsv",
             "pyTMHMM transmembrane predictions (receptor TM + N-terminal anchors)."),
        _tbl("exon_domain_architecture_features.tsv",
             "Representative domains, TM, numbered coding exons and cassette slots per protein."),
        _tbl("fgfr2_domain_architecture_qc.tsv",
             "Per-protein architecture QC (domain order, cassette position, TM agreement)."),
        _tbl("post_interpro_qc_review_case_audit.tsv",
             "Audit of the review_unusual_domain_order cases."),
        _tbl("exon_block_coordinate_reconstruction_audit.tsv",
             "Native exon-block reconstruction audit for the coordinate-artifact cases "
             "(canis/gorilla/xenopus IIIb)."),
    ] if t]

    return {
        "available": True,
        "generated_at": now_iso(),
        "step": _ARCH_STEP_NAME,
        "tm_layer": "pyTMHMM",
        "label_source": "final_pre_interpro_truth_table.tsv",
        "counts": {
            "proteins": len(qc),
            "species": len(species),
            "iiib": sum(1 for (_s, i) in qc if i == "IIIb"),
            "iiic": sum(1 for (_s, i) in qc if i == "IIIc"),
        },
        "ui_rules": {
            "interactive_from_tables": True,
            "static_figures_are_downloadable_previews": True,
            "labels_from_truth_table": True,
            "tm_layer_is_pytmhmm": True,
            "qc_warnings_visible_not_exaggerated": True,
        },
        "species": sorted(species.values(), key=lambda x: x["display_species_name"]),
        "overview_figures": overview_figures,
        "tables": tables,
        "report": rel(arch / "reports" / "post_interpro_exon_domain_architecture_summary.md")
        if (arch / "reports" / "post_interpro_exon_domain_architecture_summary.md").exists() else None,
    }


# --------------------------------------------------------------------------- #
# Module 1 — exon-domain boundary consistency (step 16 final thesis analysis)
# --------------------------------------------------------------------------- #
_BOUNDARY_CLASS_UI = {
    "aligned_to_domain_boundary": {"label": "Aligned", "color": "#1B7837",
        "tooltip": "Boundary coincides with a protein-domain boundary (0-3 aa)."},
    "near_domain_boundary": {"label": "Near boundary", "color": "#A6DBA0",
        "tooltip": "Boundary is close to a protein-domain boundary (4-15 aa)."},
    "within_domain": {"label": "Within domain", "color": "#FDB863",
        "tooltip": "Boundary lies inside a protein domain (>15 aa from its edges)."},
    "between_domains": {"label": "Between domains", "color": "#B2ABD2",
        "tooltip": "Boundary lies in a linker region between two domains."},
    "review_or_missing": {"label": "Missing / review", "color": "#D9D9D9",
        "tooltip": "No cassette / domain coordinate available for this protein."},
}
_BOUNDARY_CLASS_ORDER = ["aligned_to_domain_boundary", "near_domain_boundary",
                         "within_domain", "between_domains", "review_or_missing"]
_TAXON_ORDER = _species_order.TAXON_GROUP_ORDER


def _boundary_consistency_dir(run_dir: Path) -> Optional[Path]:
    """Locate the sibling step-16 exon-domain boundary consistency folder."""
    parent = Path(run_dir).parent
    cand = parent / "16_final_thesis_analyses" / "exon_domain_boundary_consistency"
    if (cand / "tables" / "exon_domain_boundary_distances.tsv").exists():
        return cand
    for c in sorted(parent.glob("*final_thesis*")):
        d = c / "exon_domain_boundary_consistency"
        if (d / "tables" / "exon_domain_boundary_distances.tsv").exists():
            return d
    return None


def _bc_species_order() -> Dict[str, int]:
    """The one canonical species order, shared with every other view."""
    return _species_order.reference_panel_order()


def _bc_figure_links(bc: Path) -> Dict[str, Any]:
    def _fmts(stem: str) -> Dict[str, str]:
        out: Dict[str, str] = {}
        for ext in ("png", "svg", "pdf"):
            p = bc / "figures" / f"{stem}.{ext}"
            if p.exists():
                out[ext] = rel(p)
        return out
    return {
        "heatmap": _fmts("Figure_11_exon_domain_boundary_consistency_heatmap"),
        "distance_distribution": _fmts("Figure_12_boundary_distance_distribution"),
    }


def _bc_interpretation(exon_block_status: str) -> str:
    s = (exon_block_status or "").lower()
    if s == "native_exon_blocks_reconstructed":
        return ("Coding-exon blocks were rebuilt from native CDS coordinates; the "
                "validated cassette slot is unchanged. Inspection case (display "
                "confidence), not a biological failure.")
    if s == "cassette_only_high_confidence":
        return ("Only the validated cassette slot is shown (native exon coordinates "
                "unavailable). Inspection case (display confidence), not a biological "
                "failure.")
    if s == "minor_length_clamped":
        return ("A coding-exon block was clamped 1-2 aa to protein length "
                "(codon-boundary rounding). Inspection case, not a biological failure.")
    return "Inspection case flagged for display-coordinate confidence."


def build_boundary_consistency_summary(run_dir: Path) -> Dict[str, Any]:
    bc = _boundary_consistency_dir(run_dir)
    if not bc:
        return {"available": False, "reason": "missing module 1 outputs"}
    rows = read_tsv(bc / "tables" / "exon_domain_boundary_consistency_summary.tsv")
    overall = next((r for r in rows if r.get("scope") == "overall"), {})
    counts = {c: to_int(overall.get(c)) or 0 for c in _BOUNDARY_CLASS_ORDER}
    iso_summary: Dict[str, Any] = {}
    taxon_summary: List[Dict[str, Any]] = []
    for r in rows:
        if r.get("scope") == "isoform":
            iso_summary[r.get("level", "")] = {
                "n_proteins": to_int(r.get("n_proteins")),
                "boundary_class_counts": {c: to_int(r.get(c)) or 0 for c in _BOUNDARY_CLASS_ORDER},
                "median_distance": to_float(r.get("median_distance_to_nearest_domain_boundary")),
                "mean_distance": to_float(r.get("mean_distance_to_nearest_domain_boundary")),
            }
        elif r.get("scope") == "taxon_group":
            taxon_summary.append({
                "taxon_group": r.get("level", ""),
                "n_proteins": to_int(r.get("n_proteins")),
                "boundary_class_counts": {c: to_int(r.get(c)) or 0 for c in _BOUNDARY_CLASS_ORDER},
                "median_distance": to_float(r.get("median_distance_to_nearest_domain_boundary")),
                "mean_distance": to_float(r.get("mean_distance_to_nearest_domain_boundary")),
            })
    taxon_summary.sort(key=lambda x: _TAXON_ORDER.index(x["taxon_group"])
                       if x["taxon_group"] in _TAXON_ORDER else 99)
    return {
        "available": True,
        "generated_at": now_iso(),
        "total_primary_proteins": to_int(overall.get("total_primary_proteins")),
        "proteins_with_cassette_data": to_int(overall.get("proteins_with_cassette_boundary_data")),
        "proteins_with_interpro_domain_data": to_int(overall.get("proteins_with_interpro_domain_data")),
        "proteins_with_tm_data": to_int(overall.get("proteins_with_tm_data")),
        "n_cassette_boundaries": to_int(overall.get("n_cassette_boundaries")),
        "boundary_class_counts": counts,
        "median_distance": to_float(overall.get("median_distance_to_nearest_domain_boundary")),
        "mean_distance": to_float(overall.get("mean_distance_to_nearest_domain_boundary")),
        "isoform_summary": iso_summary,
        "taxon_summary": taxon_summary,
        "key_interpretation": (
            "Cassette-end boundaries are consistently aligned/near the Ig-like domain "
            "boundary, while cassette-start boundaries usually lie within the Ig-like "
            "region. This is consistent with robust exon\u2013domain boundary "
            "identification across the primary vertebrate FGFR2 dataset."),
        "source_files": [
            rel(bc / "tables" / "exon_domain_boundary_consistency_summary.tsv"),
            rel(bc / "tables" / "exon_domain_boundary_distances.tsv"),
        ],
        "figure_links": _bc_figure_links(bc),
    }


def build_boundary_consistency_matrix(run_dir: Path) -> Dict[str, Any]:
    bc = _boundary_consistency_dir(run_dir)
    if not bc:
        return {"available": False, "reason": "missing module 1 outputs", "rows": []}
    dist = read_tsv(bc / "tables" / "exon_domain_boundary_distances.tsv")
    taxon = {r.get("species", ""): r.get("taxon_group", "") for r in _truth(run_dir)}
    order = _bc_species_order()

    def _cell(r: Dict[str, str]) -> Dict[str, Any]:
        return {
            "boundary_type": r.get("boundary_type", ""),
            "boundary_label": r.get("boundary_label", ""),
            "boundary_aa": to_int(r.get("boundary_aa")),
            "nearest_domain_label": r.get("nearest_domain_label", ""),
            "nearest_domain_class": r.get("nearest_domain_class", ""),
            "nearest_domain_start_aa": to_int(r.get("nearest_domain_start_aa")),
            "nearest_domain_end_aa": to_int(r.get("nearest_domain_end_aa")),
            "distance_to_nearest_domain_boundary": to_int(r.get("distance_to_nearest_domain_boundary")),
            "boundary_class": r.get("boundary_class", "review_or_missing"),
            "notes": r.get("notes", ""),
        }

    by_key: Dict[tuple, Dict[str, Any]] = {}
    for r in dist:
        bt = r.get("boundary_type", "")
        if bt not in ("cassette_start", "cassette_end"):
            continue
        sp, iso = r.get("species", ""), r.get("isoform", "")
        key = (sp, iso)
        node = by_key.setdefault(key, {
            "species": sp,
            "display_species_name": sp.replace("_", " ").title(),
            "isoform": iso,
            "taxon_group": taxon.get(sp, ""),
            "display_label": f"{sp.replace('_', ' ').title()} · {r.get('final_isoform_label', iso)}",
            "protein_id": r.get("protein_id", ""),
            "transcript_id": r.get("transcript_id", ""),
            "protein_length": to_int(r.get("protein_length")),
            "architecture_qc_status": r.get("architecture_qc_status", ""),
            "exon_block_display_status": r.get("exon_block_display_status", ""),
            "source_coordinate_status": r.get("source_coordinate_status", ""),
            "cells": {},
        })
        node["cells"][bt] = _cell(r)

    def _row_sort(node: Dict[str, Any]) -> tuple:
        ti = _TAXON_ORDER.index(node["taxon_group"]) if node["taxon_group"] in _TAXON_ORDER else 99
        si = order.get(node["species"], 999)
        return (ti, si, node["species"], node["isoform"])

    rows = sorted(by_key.values(), key=_row_sort)
    # normalise cells into an ordered list per row (cassette_start, cassette_end)
    for node in rows:
        node["cells"] = [node["cells"].get(bt) for bt in ("cassette_start", "cassette_end")
                         if node["cells"].get(bt)]
    return {
        "available": True,
        "generated_at": now_iso(),
        "columns": [
            {"key": "cassette_start", "label": "Cassette start"},
            {"key": "cassette_end", "label": "Cassette end"},
        ],
        "boundary_classes": _BOUNDARY_CLASS_UI,
        "row_order": "taxonomic",
        "source_table": rel(bc / "tables" / "exon_domain_boundary_distances.tsv"),
        "rows": rows,
    }


def build_boundary_consistency_outliers(run_dir: Path) -> Dict[str, Any]:
    bc = _boundary_consistency_dir(run_dir)
    if not bc:
        return {"available": False, "reason": "missing module 1 outputs", "outliers": []}
    rows = read_tsv(bc / "tables" / "exon_domain_boundary_outliers.tsv")
    out: List[Dict[str, Any]] = []
    for r in rows:
        eb = r.get("exon_block_display_status", "")
        out.append({
            "species": r.get("species", ""),
            "display_species_name": r.get("species", "").replace("_", " ").title(),
            "isoform": r.get("isoform", ""),
            "final_isoform_label": r.get("final_isoform_label", r.get("isoform", "")),
            "boundary_type": r.get("boundary_type", ""),
            "boundary_aa": to_int(r.get("boundary_aa")),
            "nearest_domain_label": r.get("nearest_domain_label", ""),
            "distance": to_int(r.get("distance_to_nearest_domain_boundary")),
            "boundary_class": r.get("boundary_class", ""),
            "reason": r.get("outlier_reason", ""),
            "interpretation": _bc_interpretation(eb),
            "exon_block_display_status": eb,
            "architecture_qc_status": r.get("architecture_qc_status", ""),
            "link_target": {"species": r.get("species", ""),
                            "isoform": r.get("isoform", ""), "tab": "boundary"},
        })
    return {
        "available": True,
        "generated_at": now_iso(),
        "count": len(out),
        "title": "Inspection cases",
        "note": ("Cases flagged for display-coordinate confidence (reconstructed or "
                 "cassette-only exon blocks). These are inspection cases, not biological "
                 "failures; IIIb/IIIc labels and primary/review membership are unchanged."),
        "source_table": rel(bc / "tables" / "exon_domain_boundary_outliers.tsv"),
        "outliers": out,
    }


def build_boundary_consistency_index(run_dir: Path) -> Dict[str, Any]:
    bc = _boundary_consistency_dir(run_dir)
    if not bc:
        return {"available": False, "reason": "missing module 1 outputs"}
    report = bc / "reports" / "exon_domain_boundary_consistency_report.md"
    return {
        "available": True,
        "created_at": now_iso(),
        "title": "Boundary Consistency Explorer",
        "step": "16_final_thesis_analyses/exon_domain_boundary_consistency",
        "source_tables": [
            rel(bc / "tables" / "exon_domain_boundary_distances.tsv"),
            rel(bc / "tables" / "exon_domain_boundary_consistency_summary.tsv"),
            rel(bc / "tables" / "exon_domain_boundary_outliers.tsv"),
        ],
        "source_figures": _bc_figure_links(bc),
        "report": rel(report) if report.exists() else None,
        "endpoints": {
            "summary": "/api/runs/current/boundary-consistency/summary",
            "matrix": "/api/runs/current/boundary-consistency/matrix",
            "outliers": "/api/runs/current/boundary-consistency/outliers",
            "index": "/api/runs/current/boundary-consistency",
        },
        "boundary_classes": _BOUNDARY_CLASS_UI,
        "boundary_class_order": _BOUNDARY_CLASS_ORDER,
        "ui_copy": {
            "outlier_section_title": "Inspection cases",
            "outlier_semantics": "inspection, not biological failure",
            "thresholds": {"aligned_to_domain_boundary": "0-3 aa",
                           "near_domain_boundary": "4-15 aa"},
            "key_interpretation": (
                "Cassette-end boundaries are consistently near Ig-like domain "
                "boundaries, while cassette-start boundaries are usually located "
                "within the Ig-like region."),
        },
    }


INDEX_BUILDERS = {
    "run_index.json": build_run_index,
    "species_index.json": build_species_index,
    "evidence_stack.json": build_evidence_stack,
    "figure_index.json": build_figure_index,
    "download_index.json": build_download_index,
    "freeze_index.json": build_freeze_index,
    # Phase-2 interactive viewers (graceful when optional inputs are absent)
    "cassette_residue_index.json": build_cassette_residue_index,
    "coordinate_track_index.json": build_coordinate_track_index,
    "msa_index.json": build_msa_index,
    "synteny_locus_index.json": build_synteny_locus_index,
    "species_story_index.json": build_species_story_index,
    # Post-InterPro / pyTMHMM domain-architecture (step 15; example freeze only)
    "domain_architecture_index.json": build_domain_architecture_index,
    "domain_architecture_summary.json": build_domain_architecture_summary,
    "species_domain_architecture.json": build_species_domain_architecture,
    "domain_architecture_qc.json": build_domain_architecture_qc,
    # Module 1 — exon-domain boundary consistency (step 16; example freeze)
    "boundary_consistency_index.json": build_boundary_consistency_index,
    "boundary_consistency_summary.json": build_boundary_consistency_summary,
    "boundary_consistency_matrix.json": build_boundary_consistency_matrix,
    "boundary_consistency_outliers.json": build_boundary_consistency_outliers,
}




def write_all(run_dir: Path, outdir: Optional[Path] = None) -> Path:
    """Build and write every website index.

    Each index is built independently and tolerant of missing inputs: a builder
    that fails (e.g. post-InterPro domain/boundary indices for a run that has only
    completed pre-InterPro) is skipped instead of aborting the whole build. This
    makes partial exploration possible after pre-InterPro, before cluster
    annotation, while a full closure run still writes all indices.
    """
    run_dir = Path(run_dir)
    outdir = Path(outdir) if outdir else (run_dir / "website_indices")
    outdir.mkdir(parents=True, exist_ok=True)
    errors: Dict[str, str] = {}
    for name, builder in INDEX_BUILDERS.items():
        try:
            payload = builder(run_dir)
        except Exception as exc:  # partial run: skip indices whose inputs are absent
            errors[name] = f"{type(exc).__name__}: {exc}"
            continue
        payload = sanitize_public_payload(payload)
        (outdir / name).write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
    # The curated FGFR2 Gallery, written last so it replaces the per-file listing above.
    # This is the production path every FGFR2 run goes through, which is what makes the
    # modern catalogue automatic rather than something a migration command has to apply.
    gallery = _write_fgfr2_gallery(run_dir, outdir)
    if gallery:
        errors["figure_index.json"] = gallery

    errors_path = outdir / "_index_build_errors.json"
    if errors:
        errors_path.write_text(json.dumps(errors, indent=2, ensure_ascii=False), encoding="utf-8")
    elif errors_path.exists():
        errors_path.unlink()
    # Raw run metadata remains a runtime record. The website receives only
    # portable projections of the two JSON files offered in its Files view.
    run_base = run_dir.parent.parent
    write_public_download_projections(run_base)
    write_freshness_contract(
        run_base, outdir, generator="scripts/build_website_indices.py")
    write_payload_contracts(
        outdir, run_id=run_base.name, dataset_id=run_base.name,
        generator="scripts/build_website_indices.py")
    return outdir


def _write_fgfr2_gallery(closure_dir: Path, outdir: Path) -> str:
    """Overwrite ``figure_index.json`` with the curated FGFR2 catalogue.

    Returns "" when the catalogue was written or does not apply, and the cause when writing
    it failed. A failure is reported and the per-file index is left in place: a plain Gallery
    is a cosmetic regression, whereas aborting the build takes the whole dataset offline.

    ``closure_dir`` is the closure, so the run base is two levels up. The validated dataset
    is not built through here — its catalogue is produced into the derived overlay by
    ``rebuild_dataset_indices`` and the freeze itself is never written to.
    """
    run_dir = Path(closure_dir).parent.parent
    if not (run_dir / "run_config.json").is_file():
        return ""
    try:
        from fgfr2 import run_gallery
        if not run_gallery.is_fgfr2_closure_run(run_dir):
            return ""
        run_gallery.write(run_dir, outdir)
    except Exception as exc:
        return f"{type(exc).__name__}: {exc}"
    return ""


def main() -> int:
    ap = argparse.ArgumentParser(description="Build Phase-1 website JSON indices from a closure run dir.")
    ap.add_argument("--run-dir", type=Path, required=True,
                    help="path to a *_final_pre_interpro_closure directory")
    ap.add_argument("--outdir", type=Path, default=None,
                    help="output dir (default: <run-dir>/website_indices)")
    args = ap.parse_args()
    if not args.run_dir.exists():
        print(f"[FAIL] run dir not found: {args.run_dir}")
        return 2
    out = write_all(args.run_dir, args.outdir)
    print(f"[OK] wrote {len(INDEX_BUILDERS)} website indices -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
