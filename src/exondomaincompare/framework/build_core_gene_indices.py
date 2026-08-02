#!/usr/bin/env python3
"""Build generic website indices from CORE gene-analysis contract outputs.

This is the gene/event-AGNOSTIC index builder for the Core-only path. It reads
the core contract TSVs (see docs/core_gene_analysis_contract.md) produced by a
core runner/adapter and writes the generic website indices — WITHOUT requiring
any event region, cassette, IIIb/IIIc logic, or marker reconciliation.

When the gene config has no configured event region:
  * event_region_index.json is written with available=false,
    reason="no_event_configured"
  * event-specific "Boundary Consistency" is disabled in available_views
  * the generic "Exon–Domain Boundaries" view is enabled when exon+domain data exist

Inputs (one of):
  --run-id <id>     read runs/<id>/results/core_gene_analysis/, config from run
  --core-dir <dir>  read core contract TSVs from an explicit directory (e.g. a mock)

Output:
  --out <dir>       default: runs/<id>/website_indices/generic (for --run-id)

Examples:
  python scripts/framework/build_core_gene_indices.py --run-id <run_id>
  python scripts/framework/build_core_gene_indices.py \
      --core-dir artifacts/core_gene_analysis/mock \
      --config configs/genes/drafts/TPM1_core_only_pilot.yaml \
      --dataset-id mock:tpm1_core_only \
      --out artifacts/generic_indices/mock_core_only
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from exondomaincompare.framework.gene_config import (  # noqa: E402
    GeneConfig, GeneConfigError, PROJECT_ROOT, load_gene_config, load_gene_config_lenient,
    resolve_run_analysis, default_gene_config,
)
from exondomaincompare.framework.primary_selection import (  # noqa: E402
    build_primary_selection, write_selection_evidence,
)
from exondomaincompare.config import load_config  # noqa: E402

SCHEMA_VERSION = 1
NEAR_EDGE_THRESHOLD_AA = 5
FREEZE_ROOT = PROJECT_ROOT / "results" / "final_30_until_interpro_prepare"
RUNTIME_CONFIG = load_config(repository_root=PROJECT_ROOT)
LOCAL_RUNS_ROOT = RUNTIME_CONFIG.runs_root
from exondomaincompare.runs.registry import RegistryError, resolve_run_record  # noqa: E402
from exondomaincompare.runs.legacy import LegacyRunAdapter  # noqa: E402


# --------------------------------------------------------------------------- #
# IO helpers (kept local so the framework builder is self-contained)
# --------------------------------------------------------------------------- #
def read_tsv(p: Path) -> List[Dict[str, str]]:
    if not Path(p).is_file():
        return []
    with open(p, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_json(p: Path, default: Any = None) -> Any:
    try:
        return json.loads(Path(p).read_text(encoding="utf-8"))
    except Exception:
        return default


def write_json(p: Path, data: Any) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def _int_or_none(v: Any) -> Optional[int]:
    try:
        if v is None or v == "":
            return None
        return int(float(v))
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------------- #
# core-contract source
# --------------------------------------------------------------------------- #
class CoreSource:
    def __init__(self, core_dir: Path, dataset_id: str, cfg: GeneConfig):
        self.core_dir = core_dir
        self.dataset_id = dataset_id
        self.cfg = cfg
        self.report = read_json(core_dir / "core_gene_report.json", {}) or {}

    def tsv(self, name: str) -> List[Dict[str, str]]:
        return read_tsv(self.core_dir / name)


# --------------------------------------------------------------------------- #
# index builders (generic; no event assumptions)
# --------------------------------------------------------------------------- #
def _gene_identity_for(src: "CoreSource") -> Dict[str, Any]:
    """The requested-versus-source symbol record, read from the run's own config.

    ``core_dir`` is ``runs/<id>/results/core_gene_analysis``, so the run root is two levels
    up. Absent for runs created before the identity record existed, and for genes whose
    annotation symbol already matches the request.
    """
    run_dir = src.core_dir.parent.parent
    data = read_json(run_dir / "run_config.json", {}) or {}
    identity = data.get("gene_identity")
    if not isinstance(identity, dict) or not identity:
        return {}
    per_species = data.get("gene_identity_by_species")
    return {**identity,
            "by_species": per_species if isinstance(per_species, dict) else {}}


def _display_species(sid: str) -> str:
    return (sid or "").replace("_", " ").strip().title()


def _domain_index(src: CoreSource) -> Dict[str, Any]:
    # domain_features.tsv is the curated layered model (layer = domain|family|feature);
    # interpro_annotations.tsv is the raw-signature provenance layer.
    curated = src.tsv("domain_features.tsv")
    raw_sigs = src.tsv("interpro_annotations.tsv")
    tms = src.tsv("tm_features.tsv")
    exon_map = src.tsv("exon_protein_map.tsv")
    iso = src.tsv("protein_isoform_index.tsv")

    meta = {r.get("protein_id", ""): {
        "length_aa": _int_or_none(r.get("protein_length")),
        "role": "primary" if str(r.get("primary_status", "")).lower() == "primary" else "alternative",
        "transcript_id": r.get("transcript_id", ""),
    } for r in iso if r.get("protein_id")}

    proteins: Dict[str, Dict[str, Any]] = {}

    def prot(sp: str, pid: str) -> Dict[str, Any]:
        key = f"{sp}||{pid}"
        m = meta.get(pid, {})
        return proteins.setdefault(key, {
            "species_id": sp, "protein_id": pid,
            "length_aa": m.get("length_aa"),
            "role": m.get("role", "alternative"),
            "transcript_id": m.get("transcript_id", ""),
            "domains": [], "families": [], "features": [],
            "raw_signatures": [], "tm": [], "exons": [],
        })

    def _annot(d: Dict[str, Any]) -> Dict[str, Any]:
        return {
            "interpro_accession": d.get("interpro_accession", ""),
            "interpro_name": d.get("interpro_name", "") or d.get("domain_name", ""),
            "interpro_type": d.get("interpro_type", ""),
            "member_databases": d.get("member_databases", "") or d.get("domain_source", ""),
            "supporting_interpro": d.get("supporting_interpro", ""),
            "n_signatures": _int_or_none(d.get("n_signatures")),
            "representative_signature": d.get("representative_signature", ""),
            "start_aa": _int_or_none(d.get("start_aa")),
            "end_aa": _int_or_none(d.get("end_aa")),
            # legacy keys kept so older consumers/tooltips still resolve
            "domain_id": d.get("domain_id", "") or d.get("interpro_accession", ""),
            "domain_name": d.get("domain_name", "") or d.get("interpro_name", ""),
            "domain_source": d.get("member_databases", "") or d.get("domain_source", ""),
        }

    for d in curated:
        p = prot(d.get("species_id", ""), d.get("protein_id", ""))
        layer = (d.get("layer") or "").strip()
        if layer == "domain":
            p["domains"].append(_annot(d))
        elif layer == "family":
            p["families"].append(_annot(d))
        elif layer == "feature":
            p["features"].append(_annot(d))
    for r in raw_sigs:
        p = prot(r.get("species_id", ""), r.get("protein_id", ""))
        p["raw_signatures"].append({
            "member_database": r.get("member_database", ""),
            "signature_accession": r.get("signature_accession", ""),
            "signature_name": r.get("signature_name", ""),
            "interpro_accession": r.get("interpro_accession", ""),
            "interpro_name": r.get("interpro_name", ""),
            "interpro_type": r.get("interpro_type", ""),
            "layer": r.get("layer", ""),
            "is_integrated": str(r.get("is_integrated", "")) in ("1", "True", "true"),
            "start_aa": _int_or_none(r.get("start_aa")),
            "end_aa": _int_or_none(r.get("end_aa")),
        })
    for t in tms:
        p = prot(t.get("species_id", ""), t.get("protein_id", ""))
        p["tm"].append({
            "start_aa": _int_or_none(t.get("start_aa")),
            "end_aa": _int_or_none(t.get("end_aa")),
            "source": t.get("source", ""),
            "topology": t.get("topology", "transmembrane"),
        })
    for e in exon_map:
        p = prot(e.get("species_id", ""), e.get("protein_id", ""))
        p["exons"].append({
            "exon_id": e.get("exon_id", ""),
            "exon_number": _int_or_none(e.get("exon_number")),
            "protein_start_aa": _int_or_none(e.get("protein_start_aa")),
            "protein_end_aa": _int_or_none(e.get("protein_end_aa")),
        })

    def _dedup(items: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        seen, out = set(), []
        for d in sorted(items, key=lambda x: (x.get("start_aa") or 0, x.get("end_aa") or 0)):
            k = (d.get("interpro_accession"), d.get("start_aa"), d.get("end_aa"))
            if k not in seen:
                seen.add(k)
                out.append(d)
        return out

    for p in proteins.values():
        p["domains"] = _dedup(p["domains"])
        p["families"] = _dedup(p["families"])
        p["features"] = _dedup(p["features"])
        seen_t, uniq_t = set(), []
        for t in sorted(p["tm"], key=lambda x: (x["start_aa"] or 0, x["end_aa"] or 0)):
            k = (t["start_aa"], t["end_aa"])
            if k not in seen_t:
                seen_t.add(k)
                uniq_t.append(t)
        p["tm"] = uniq_t
        p["raw_signatures"].sort(key=lambda x: (x.get("start_aa") or 0, x.get("end_aa") or 0))
        p["exons"].sort(key=lambda x: (x["protein_start_aa"] or 0))
        if not p.get("length_aa"):
            p["length_aa"] = max((e["protein_end_aa"] or 0) for e in p["exons"]) if p["exons"] else None
        # a protein is "annotated" once it carries any real InterPro/pyTMHMM result
        p["annotated"] = bool(p["domains"] or p["families"] or p["features"] or p["tm"])
        p["n_representative_domains"] = len(p["domains"])

    # group by species (primary first, annotated before unannotated)
    by_species: Dict[str, Dict[str, Any]] = {}
    for p in proteins.values():
        node = by_species.setdefault(p["species_id"], {
            "species": p["species_id"],
            "display_species_name": _display_species(p["species_id"]),
            "proteins": [],
        })
        node["proteins"].append(p)
    for node in by_species.values():
        node["proteins"].sort(key=lambda p: (p["role"] != "primary", not p["annotated"], p["protein_id"]))
        node["n_annotated_proteins"] = sum(1 for p in node["proteins"] if p["annotated"])

    # dynamic legend: distinct InterPro entry types actually present in the data
    type_set = sorted({d.get("interpro_type", "") for p in proteins.values()
                       for d in (p["domains"] + p["families"] + p["features"])
                       if d.get("interpro_type")})
    has_domains = any(p["domains"] for p in proteins.values())
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": src.cfg.analysis_id,
        "gene_symbol": src.cfg.gene_symbol,
        "mode": "generic",
        "available": has_domains or any(p["families"] or p["features"] or p["tm"]
                                        for p in proteins.values()),
        "domain_layer": "representative_domain",
        "entry_types_present": type_set,
        "relation_description": "Protein domain architecture (core, gene-agnostic; "
                                "representative InterPro domain layer + separate family / "
                                "feature / topology layers).",
        "species": list(by_species.values()),
        "source": "core_gene_analysis",
    }


def _synteny_index(src: CoreSource) -> Dict[str, Any]:
    rows = src.tsv("synteny_neighbors.tsv")
    by_species: Dict[str, Dict[str, Any]] = {}
    n_resolved = 0
    for r in rows:
        sp = r.get("species_id", "")
        node = by_species.setdefault(sp, {"species": sp, "neighbors": []})
        resolved = str(r.get("status", "")).lower() == "resolved"
        n_resolved += 1 if resolved else 0
        node["neighbors"].append({
            "symbol": r.get("neighbor_symbol", ""),
            "side": r.get("side", ""),
            "order": _int_or_none(r.get("order")),
            "orientation": r.get("orientation", ""),
            "classification": r.get("classification", ""),
            "status": r.get("status", ""),
        })
    available = n_resolved > 0
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": src.cfg.analysis_id,
        "available": available,
        "synteny_status": "computed" if available else "not_computed",
        "synteny_reason": "" if available else
        "No resolved synteny neighbours in core outputs (annotation may be unavailable).",
        "n_resolved_neighbors": n_resolved,
        "species": list(by_species.values()),
        "source": "core_gene_analysis",
    }


def _exon_domain_boundary_index(src: CoreSource) -> Dict[str, Any]:
    rows = src.tsv("exon_domain_boundary_distances.tsv")
    counts: Dict[str, int] = {}
    by_protein: Dict[str, Dict[str, Any]] = {}
    for r in rows:
        cat = r.get("category", "")
        counts[cat] = counts.get(cat, 0) + 1
        key = f"{r.get('species_id','')}||{r.get('protein_id','')}"
        node = by_protein.setdefault(key, {
            "species_id": r.get("species_id", ""),
            "protein_id": r.get("protein_id", ""),
            "boundaries": [],
        })
        node["boundaries"].append({
            "analysis_id": r.get("analysis_id", ""),
            "gene_symbol": r.get("gene_symbol", ""),
            "transcript_id": r.get("transcript_id", ""),
            "exon_boundary_id": r.get("exon_boundary_id", ""),
            "boundary_position_aa": _int_or_none(r.get("boundary_position_aa")),
            "nearest_domain_accession": r.get("nearest_domain_accession", "") or r.get("nearest_domain_id", ""),
            "nearest_domain_id": r.get("nearest_domain_id", ""),
            "nearest_domain_name": r.get("nearest_domain_name", ""),
            "nearest_domain_type": r.get("nearest_domain_type", ""),
            "nearest_edge": r.get("nearest_edge", "") or r.get("domain_edge_type", ""),
            "nearest_domain_boundary_type": r.get("nearest_domain_boundary_type", ""),
            "domain_edge_type": r.get("domain_edge_type", ""),
            "signed_distance_aa": _int_or_none(r.get("signed_distance_aa")),
            "absolute_distance_aa": _int_or_none(r.get("absolute_distance_aa")),
            "distance_aa": _int_or_none(r.get("distance_aa")),
            "classification": r.get("classification", "") or cat,
            "category": cat,
            "domain_layer": r.get("domain_layer", "representative_domain"),
            "source": r.get("source", ""),
        })
    proteins_list = list(by_protein.values())
    species_scope = sorted({p["species_id"] for p in proteins_list if p.get("species_id")})
    # normalise counts to the explicit generic classification vocabulary
    ordered_counts = {c: counts.get(c, 0) for c in
                      ("exact_edge", "near_edge", "inside_domain", "outside_domain", "unknown")}
    for k, v in counts.items():
        ordered_counts.setdefault(k, v)
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": src.cfg.analysis_id,
        "gene_symbol": src.cfg.gene_symbol,
        "available": bool(rows),
        "relation_label": src.cfg.ui_labels.get("boundary_relation", "Exon–Domain Boundaries"),
        "scope": "internal_coding_exon_boundaries",
        "protein_scope": "primary_only",
        "isoform_scope": "primary_only",
        "boundary_scope": "internal_coding_exon_boundaries",
        "domain_layer": "representative_domain",
        "near_edge_threshold_aa": NEAR_EDGE_THRESHOLD_AA,
        "species_scope": species_scope,
        "n_species": len(species_scope),
        "n_proteins": len(proteins_list),
        "category_counts": ordered_counts,
        "n_boundaries": len(rows),
        "proteins": proteins_list,
        "source": "core_gene_analysis",
    }


def _event_region_index(src: CoreSource, has_event: bool) -> Dict[str, Any]:
    if not has_event:
        return {
            "schema_version": SCHEMA_VERSION,
            "dataset_id": src.dataset_id,
            "analysis_id": src.cfg.analysis_id,
            "available": False,
            "reason": "no_event_configured",
            "message": "No event region is configured for this gene. "
                       "Core gene-level analysis is available.",
            "event_labels": [],
            "species": [],
        }
    # Event configured: core builder does not compute event specifics; point to
    # the event-detector/adapter path instead.
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": src.cfg.analysis_id,
        "available": False,
        "reason": "built_by_event_adapter",
        "message": "Event region is configured; build event indices with the "
                   "event-detector adapter, not the core builder.",
        "event_labels": src.cfg.event_label_ids,
        "species": [],
    }


def _primary_proteins(iso: List[Dict[str, str]]) -> List[Dict[str, str]]:
    return [r for r in iso if str(r.get("primary_status", "")).lower() == "primary"]


def _gene_analysis_index(src: CoreSource) -> Dict[str, Any]:
    gm = src.tsv("gene_model_index.tsv")
    iso = src.tsv("protein_isoform_index.tsv")
    exon_map = src.tsv("exon_protein_map.tsv")
    species = sorted({r.get("species_id", "") for r in gm if r.get("species_id")})
    primary = _primary_proteins(iso)
    # selection method / rule are documented in the core report (gene-agnostic).
    selection_method = (src.report.get("selection_method")
                        or src.report.get("selection_rule_applied") or "")
    isoforms = [{
        "species_id": r.get("species_id", ""),
        "protein_id": r.get("protein_id", ""),
        "transcript_id": r.get("transcript_id", ""),
        "protein_length": _int_or_none(r.get("protein_length")),
        "primary_status": r.get("primary_status", ""),
    } for r in iso]
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": src.cfg.analysis_id,
        "analysis_display_name": src.cfg.analysis_display_name,
        "gene": {
            "symbol": src.cfg.gene_symbol,
            "display_name": src.cfg.gene_display_name,
            "reference_species": src.cfg.reference_species,
        },
        "analysis_modes": src.cfg.analysis_modes,
        "has_event": src.cfg.has_event,
        "event": {
            "status": src.cfg.event_status,
            "type": src.cfg.event_type,
            "labels": src.cfg.event_label_ids,
        },
        "selection_method": selection_method,
        "species": species,
        "n_species": len(species),
        "n_gene_models": len(gm),
        # n_proteins kept for backward compatibility (== isoform count).
        "n_proteins": len(iso),
        "n_protein_isoforms": len(iso),
        "n_primary_proteins": len(primary),
        "primary_protein_ids": [r.get("protein_id", "") for r in primary],
        "exon_map_available": bool(exon_map),
        "isoforms": isoforms,
        "source": "core_gene_analysis",
    }


def _count_tsv_rows(p: Path) -> int:
    if not Path(p).is_file():
        return -1
    try:
        with open(p, "r", encoding="utf-8", newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        return max(0, len(rows) - 1) if rows else 0
    except Exception:
        return 0


def _capability_report(src: CoreSource, gene_idx: Dict[str, Any],
                       domain_idx: Dict[str, Any], synteny_idx: Dict[str, Any],
                       boundary_idx: Dict[str, Any],
                       n_evidence: int, n_clusters: int) -> Dict[str, Any]:
    """PART 8 capability report — a single, honest status object for the UI.

    Domain-dependent milestones are reported as ``pending`` (not ``failed``)
    before the cluster InterProScan/pyTMHMM step, and only become ``available``
    once real cluster outputs exist.
    """
    cfg = src.cfg
    report = src.report or {}
    n_models = gene_idx["n_gene_models"]
    core_ok = n_models > 0 and gene_idx["n_primary_proteins"] > 0
    domain_available = bool(domain_idx.get("available"))
    boundary_available = bool(boundary_idx.get("available"))
    domain_status = str(report.get("domain_status", "") or "").lower()

    # cluster/domain status: complete once real domain outputs exist; otherwise
    # pending (never "failed" just because the cluster step has not run yet).
    if domain_available:
        cluster_status = "complete"
    elif domain_status in ("failed", "error"):
        cluster_status = "failed"
    else:
        cluster_status = "pending"

    def _pending_or(available: bool) -> str:
        if available:
            return "available"
        return "failed" if cluster_status == "failed" else "pending"

    has_event = cfg.has_event
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": cfg.analysis_id,
        "gene_symbol": cfg.gene_symbol,
        "generated_at": report.get("generated_at", ""),
        "core_model_collection": "available" if core_ok else "failed",
        "protein_isoforms_count": gene_idx["n_protein_isoforms"],
        "primary_proteins_count": gene_idx["n_primary_proteins"],
        "gene_models_count": n_models,
        "species_count": gene_idx["n_species"],
        "exon_map": "available" if gene_idx.get("exon_map_available") else "unavailable",
        "synteny": "available" if synteny_idx.get("available") else "unavailable",
        "synteny_neighbours_count": synteny_idx.get("n_resolved_neighbors", 0),
        "cluster_status": cluster_status,
        "domain_architecture": _pending_or(domain_available),
        "exon_domain_boundaries": _pending_or(boundary_available),
        "event_configured": bool(has_event),
        "exploratory_event_evidence": "available" if n_evidence > 0 else "none",
        "candidate_clusters_count": max(n_clusters, 0),
        "event_analysis_enabled": False if not has_event else bool(has_event),
        "support_level": cfg.support_level,
        "selection_method": gene_idx.get("selection_method", ""),
        "cluster_command": (
            f".venv/bin/python scripts/edc.py cluster roundtrip "
            f"--run-id {src.dataset_id[4:]}" if src.dataset_id.startswith("run:") else ""),
    }


def _ensure_primary_selection(src: CoreSource) -> Dict[str, Any]:
    """Load (or build+persist) the primary-selection evidence report.

    If primary_selection_report.json is absent it is derived from the isoform
    index + collection report using the documented selection hierarchy and
    written next to the other core outputs (so the UI + TSV stay in sync).
    """
    report = read_json(src.core_dir / "primary_selection_report.json", None)
    if isinstance(report, dict) and report.get("proteins"):
        return report
    iso = src.tsv("protein_isoform_index.tsv")
    coll = read_json(src.core_dir / "core_model_collection_report.json", {}) or {}
    report = build_primary_selection(iso, collection_report=coll)
    try:
        write_selection_evidence(report,
                                 src.core_dir / "primary_selection_evidence.tsv",
                                 src.core_dir / "primary_selection_report.json")
    except Exception:
        pass
    return report


def _rel(core_dir: Path, name: str) -> str:
    p = core_dir / name
    if not p.is_file():
        return ""
    try:
        return str(p.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(p)


def _evidence_stack_index(src: CoreSource, gene_idx: Dict[str, Any],
                          synteny_idx: Dict[str, Any], domain_idx: Dict[str, Any],
                          sel: Dict[str, Any], n_evidence: int, n_clusters: int,
                          uniprot: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Ordered evidence-stack items (source / status / confidence / explanation /
    file), mirroring the FGFR2 evidence-stack concept for a generic gene."""
    cd = src.core_dir
    lengths = [i.get("protein_length") for i in gene_idx.get("isoforms", []) if i.get("protein_length")]
    len_span = f"{min(lengths)}–{max(lengths)} aa" if lengths else "—"
    species = ", ".join(gene_idx.get("species", [])) or "—"
    coll = read_json(cd / "core_model_collection_report.json", {}) or {}
    assembly = (coll.get("source", {}) or {}).get("assembly_accession", "")

    items: List[Dict[str, Any]] = []
    items.append({
        "id": "gene_model", "title": "Gene model evidence",
        "source": f"NCBI RefSeq annotation{f' · {assembly}' if assembly else ''}",
        "status": "accepted" if gene_idx["n_gene_models"] > 0 else "failed",
        "confidence": "high",
        "explanation": (f"{gene_idx['n_gene_models']} transcript/gene-model record(s) collected for "
                        f"{gene_idx['gene']['symbol']} in {species}."),
        "file": _rel(cd, "gene_model_index.tsv"),
    })
    items.append({
        "id": "protein_isoform", "title": "Protein isoform evidence",
        "source": "RefSeq protein set",
        "status": "accepted" if gene_idx["n_protein_isoforms"] > 0 else "failed",
        "confidence": "high",
        "explanation": (f"{gene_idx['n_protein_isoforms']} protein-coding isoform(s) "
                        f"(lengths {len_span})."),
        "file": _rel(cd, "protein_isoform_index.tsv"),
    })
    items.append({
        "id": "primary_selection", "title": "Primary protein selection",
        "source": sel.get("selection_source", "") or sel.get("selection_rule_label", ""),
        "status": sel.get("evidence_status", "accepted"),
        "confidence": sel.get("confidence", "medium"),
        "explanation": sel.get("explanation", ""),
        "rule": sel.get("selection_rule", ""),
        "rule_label": sel.get("selection_rule_label", ""),
        "file": _rel(cd, "primary_selection_evidence.tsv"),
    })
    items.append({
        "id": "synteny", "title": "Synteny evidence",
        "source": "Neighbouring-gene order (RefSeq annotation)",
        "status": "accepted" if synteny_idx.get("available") else "unavailable",
        "confidence": "medium" if synteny_idx.get("available") else "n/a",
        "explanation": (f"{synteny_idx.get('n_resolved_neighbors', 0)} synteny neighbour(s) resolved "
                        f"around the target locus." if synteny_idx.get("available") else
                        synteny_idx.get("synteny_reason", "Synteny not available.")),
        "file": _rel(cd, "synteny_neighbors.tsv"),
    })
    items.append({
        "id": "exploratory_event", "title": "Exploratory event evidence",
        "source": "Isoform-difference scan",
        "status": "exploratory" if n_clusters > 0 else "none",
        "confidence": "exploratory",
        "explanation": (f"{n_clusters} candidate cluster(s) from {n_evidence} isoform-difference "
                        f"comparison(s). Exploratory only — not validated event regions."
                        if n_clusters > 0 else
                        # The scan applies no length threshold, so "none" means the
                        # isoforms encode identical proteins — say that, don't imply
                        # something was filtered away.
                        "No protein-isoform difference block was detected."),
        "file": _rel(cd, "event_region_candidate_clusters.tsv"),
    })
    up_status, up_expl = "unavailable", "External curated evidence not collected."
    if isinstance(uniprot, dict):
        if uniprot.get("status") == "uniprot_evidence_appended":
            up_status = "available"
            up_expl = f"{uniprot.get('n_curated_rows_appended', 0)} curated UniProt annotation(s) mapped."
        else:
            up_status = "not_found"
            up_expl = uniprot.get("reason", "No curated UniProt evidence found.")
    items.append({
        "id": "external_uniprot", "title": "External curated evidence (UniProt)",
        "source": "UniProt", "status": up_status, "confidence": "medium" if up_status == "available" else "n/a",
        "explanation": up_expl, "file": _rel(cd, "uniprot_event_evidence_report.json"),
    })
    dom_available = bool(domain_idx.get("available"))
    items.append({
        "id": "domain", "title": "Domain evidence",
        "source": "InterProScan / pyTMHMM (cluster)",
        "status": "available" if dom_available else "pending",
        "confidence": "high" if dom_available else "n/a",
        "explanation": ("Domain architecture available." if dom_available else
                        "Domain architecture pending — requires the cluster InterProScan / pyTMHMM step."),
        "file": _rel(cd, "domain_features.tsv"),
    })
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": src.cfg.analysis_id,
        "gene_symbol": src.cfg.gene_symbol,
        "items": items,
    }


def _protein_architecture_index(src: CoreSource, gene_idx: Dict[str, Any]) -> Dict[str, Any]:
    """Exon/protein track data (generic). Candidate regions are exploratory
    overlays on the primary protein; domains/TM appear after the cluster step."""
    gm = src.tsv("gene_model_index.tsv")
    iso = src.tsv("protein_isoform_index.tsv")
    exon_map = src.tsv("exon_protein_map.tsv")
    domain_rows = src.tsv("domain_features.tsv")
    tm_rows = src.tsv("tm_features.tsv")
    clusters = src.tsv("event_region_candidate_clusters.tsv")

    meta = {r.get("protein_id", ""): {
        "transcript_id": r.get("transcript_id", ""),
        "length_aa": _int_or_none(r.get("protein_length")),
        "role": "primary" if str(r.get("primary_status", "")).lower() == "primary" else "alternative",
    } for r in iso if r.get("protein_id")}

    def exons_for(pid):
        out = [{"exon_id": e.get("exon_id", ""), "exon_number": _int_or_none(e.get("exon_number")),
                "protein_start_aa": _int_or_none(e.get("protein_start_aa")),
                "protein_end_aa": _int_or_none(e.get("protein_end_aa")),
                "confidence": e.get("confidence", "")}
               for e in exon_map if e.get("protein_id") == pid]
        out.sort(key=lambda x: (x["protein_start_aa"] or 0))
        return out

    def domains_for(pid):
        return [{"domain_source": d.get("domain_source", ""), "domain_id": d.get("domain_id", ""),
                 "domain_name": d.get("domain_name", ""),
                 "start_aa": _int_or_none(d.get("start_aa")), "end_aa": _int_or_none(d.get("end_aa"))}
                for d in domain_rows if d.get("protein_id") == pid]

    def tms_for(pid):
        return [{"start_aa": _int_or_none(t.get("start_aa")), "end_aa": _int_or_none(t.get("end_aa")),
                 "source": t.get("source", "")} for t in tm_rows if t.get("protein_id") == pid]

    def clusters_for(is_primary):
        if not is_primary:
            return []
        return [{"candidate_cluster_id": c.get("candidate_cluster_id", ""),
                 "start_aa": _int_or_none(c.get("representative_start_aa")),
                 "end_aa": _int_or_none(c.get("representative_end_aa")),
                 "length_aa": _int_or_none(c.get("representative_length_aa")),
                 "support_count": _int_or_none(c.get("support_count")),
                 "confidence": c.get("confidence", ""),
                 "exon_aligned": (_int_or_none(c.get("exon_aligned_support")) or 0) > 0,
                 "evidence_status": "exploratory"} for c in clusters]

    by_species: Dict[str, Dict[str, Any]] = {}
    seen = set()
    for r in gm:
        pid = r.get("protein_id", "")
        if not pid or r.get("model_status") != "protein_coding" or pid in seen:
            continue
        seen.add(pid)
        sp = r.get("species_id", "")
        node = by_species.setdefault(sp, {"species_id": sp, "proteins": []})
        m = meta.get(pid, {})
        is_primary = m.get("role") == "primary"
        node["proteins"].append({
            "protein_id": pid, "transcript_id": m.get("transcript_id", r.get("transcript_id", "")),
            "length_aa": m.get("length_aa") or _int_or_none(r.get("protein_length")),
            "role": m.get("role", "alternative"),
            "exons": exons_for(pid), "candidate_regions": clusters_for(is_primary),
            "domains": domains_for(pid), "tm_regions": tms_for(pid),
        })
    for node in by_species.values():
        node["proteins"].sort(key=lambda p: (p["role"] != "primary", p["protein_id"]))
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": src.cfg.analysis_id,
        "gene_symbol": gene_idx["gene"]["symbol"],
        "domain_status": "available" if domain_rows else "pending_cluster",
        "selection_method": gene_idx.get("selection_method", ""),
        "species": list(by_species.values()),
    }


def _event_evidence_index(src: CoreSource, domain_available: bool,
                          uniprot: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Source-by-source exploratory candidate model (Part 5). Each cluster gets an
    isoform-difference / exon-alignment / external-curated / domain evidence
    breakdown plus a cautious interpretation. NOT a validated event."""
    clusters = src.tsv("event_region_candidate_clusters.tsv")
    evidence = src.tsv("event_region_evidence.tsv")
    gene = src.cfg.gene_symbol

    def rows_for(c):
        s = _int_or_none(c.get("representative_start_aa"))
        e = _int_or_none(c.get("representative_end_aa"))
        out = []
        for r in evidence:
            if (r.get("species_id", "") != c.get("species_id", "")):
                continue
            rs, re_ = _int_or_none(r.get("region_start_aa")), _int_or_none(r.get("region_end_aa"))
            if rs is None or re_ is None or s is None or e is None:
                continue
            if (rs - 5) <= e and (s - 5) <= re_:
                out.append(r)
        return out

    up_found = isinstance(uniprot, dict) and uniprot.get("status") == "uniprot_evidence_appended"
    out_clusters = []
    for c in clusters:
        proteins = [p for p in (c.get("proteins_involved", "") or "").split(";") if p]
        transcripts = [t for t in (c.get("transcripts_involved", "") or "").split(";") if t]
        exon_aligned_n = _int_or_none(c.get("exon_aligned_support")) or 0
        support = _int_or_none(c.get("support_count")) or 0
        raw = rows_for(c)
        conf = c.get("confidence", "")
        interp = (f"This region differs between {len(proteins) or 'multiple'} {gene} protein "
                  f"isoform(s)"
                  + (" and aligns with exon boundaries" if exon_aligned_n > 0 else "")
                  + ". It is an exploratory candidate region, not a validated event. It may become "
                  "more interpretable after domain annotation or external curated evidence.")
        out_clusters.append({
            "candidate_cluster_id": c.get("candidate_cluster_id", ""),
            "region": f"{c.get('representative_start_aa','')}–{c.get('representative_end_aa','')}",
            "start_aa": _int_or_none(c.get("representative_start_aa")),
            "end_aa": _int_or_none(c.get("representative_end_aa")),
            "length_aa": _int_or_none(c.get("representative_length_aa")),
            "confidence": conf,
            "status": "exploratory",
            "evidence": {
                "isoform_difference": {
                    "n_isoform_pairs": support,
                    "proteins": proteins, "transcripts": transcripts,
                    "confidence": conf,
                    "explanation": (f"Supported by {support} isoform pair comparison(s) across "
                                    f"{len(proteins)} isoform(s)."),
                },
                "exon_alignment": {
                    "aligned": exon_aligned_n > 0,
                    "n_exon_aligned": exon_aligned_n,
                    "confidence": "medium" if exon_aligned_n > 0 else "low",
                    "explanation": (f"{exon_aligned_n} supporting comparison(s) align with exon "
                                    f"boundaries." if exon_aligned_n > 0 else
                                    "Region does not align cleanly with exon boundaries."),
                },
                "external_curated": {
                    "uniprot": "found" if up_found else ("unavailable" if uniprot is None else "not_found"),
                    "confidence": "medium" if up_found else "n/a",
                    "explanation": ("Curated UniProt alternative-sequence evidence mapped." if up_found
                                    else "No curated UniProt alternative-sequence evidence for this region."),
                },
                "domain": {
                    "status": "available" if domain_available else "pending",
                    "explanation": ("Domain overlap can be evaluated." if domain_available else
                                    "Pending until cluster domain annotation; later: overlaps domain / outside domain."),
                },
            },
            "interpretation": interp,
            "raw_support_rows": raw,
        })
    # strongest by support
    out_clusters.sort(key=lambda c: -(c["evidence"]["isoform_difference"]["n_isoform_pairs"]))
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": src.cfg.analysis_id,
        "gene_symbol": gene,
        "available": bool(out_clusters),
        "evidence_status": "exploratory",
        "n_clusters": len(out_clusters),
        "n_raw_evidence": len(evidence),
        "clusters": out_clusters,
        "explanation": ("Candidate regions are sequence differences between protein isoforms of the "
                        "same gene. They are exploratory and NOT validated event regions."),
    }


def _figures_index(src: CoreSource, arch_idx: Dict[str, Any], synteny_idx: Dict[str, Any],
                   event_idx: Dict[str, Any], domain_available: bool) -> Dict[str, Any]:
    """Stage-aware generic figures (Part 7). Pre-cluster plots are generic and
    computable locally; domain/boundary plots are pending until the cluster step."""
    has_arch = any((s.get("proteins") for s in arch_idx.get("species", [])))
    available = []
    if has_arch:
        available.append({"id": "protein_architecture", "title": "Gene / protein architecture",
                          "kind": "exon_protein_track", "available": True,
                          "caption": "Primary protein exon architecture with exploratory candidate overlays."})
    if synteny_idx.get("available"):
        available.append({"id": "synteny_neighbourhood", "title": "Synteny neighbourhood",
                          "kind": "synteny_strip", "available": True,
                          "caption": "Upstream/downstream gene neighbourhood around the target locus."})
    if event_idx.get("n_clusters", 0) > 0:
        available.append({"id": "event_evidence", "title": "Exploratory candidate evidence",
                          "kind": "candidate_overview", "available": True,
                          "caption": "Exploratory isoform-difference candidate regions (not validated events)."})
    pending = [
        {"id": "domain_architecture", "title": "Domain architecture", "requires": "cluster",
         "available": domain_available},
        {"id": "exon_domain_boundary_distribution", "title": "Exon–domain boundary distribution",
         "requires": "cluster", "available": domain_available},
    ]
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": src.cfg.analysis_id,
        "gene_symbol": src.cfg.gene_symbol,
        "kind": "core_only",
        "available": available,
        "pending": pending,
    }


def build_core_indices(src: CoreSource, out_dir: Path) -> Dict[str, Any]:
    cfg = src.cfg
    has_event = cfg.has_event

    domain_idx = _domain_index(src)
    synteny_idx = _synteny_index(src)
    boundary_idx = _exon_domain_boundary_index(src)
    gene_idx = _gene_analysis_index(src)

    # Exploratory event-evidence + candidate-cluster counts (Part 4/5). These
    # files are optional; missing == none.
    n_evidence = max(_count_tsv_rows(src.core_dir / "event_region_evidence.tsv"), 0)
    n_clusters = max(_count_tsv_rows(src.core_dir / "event_region_candidate_clusters.tsv"), 0)
    capability = _capability_report(src, gene_idx, domain_idx, synteny_idx,
                                    boundary_idx, n_evidence, n_clusters)

    # Scientific primary-selection evidence (Part 3) + generic FGFR2-parity
    # indices (Part 9). Domain-dependent indices carry an explicit
    # status=pending_cluster before InterProScan so nothing looks precomputed.
    domain_available = bool(domain_idx.get("available"))
    uniprot = read_json(src.core_dir / "uniprot_event_evidence_report.json", None)
    sel = _ensure_primary_selection(src)
    gene_idx["selection_rule"] = sel.get("selection_rule", "")
    gene_idx["selection_rule_label"] = sel.get("selection_rule_label", "")
    gene_idx["primary_selection"] = {
        "primary_protein_id": sel.get("primary_protein_id", ""),
        "selection_rule": sel.get("selection_rule", ""),
        "selection_source": sel.get("selection_source", ""),
        "confidence": sel.get("confidence", ""),
        "explanation": sel.get("explanation", ""),
    }
    evidence_stack = _evidence_stack_index(src, gene_idx, synteny_idx, domain_idx, sel,
                                           n_evidence, n_clusters, uniprot)
    arch_idx = _protein_architecture_index(src, gene_idx)
    event_evidence_idx = _event_evidence_index(src, domain_available, uniprot)
    figures_idx = _figures_index(src, arch_idx, synteny_idx, event_evidence_idx, domain_available)

    # explicit pending_cluster status on domain-dependent indices
    domain_idx["status"] = "available" if domain_available else "pending_cluster"
    boundary_idx["status"] = "available" if boundary_idx.get("available") else "pending_cluster"

    # available_views: core (gene-agnostic) always allowed if data exists; event
    # views only if an event region is configured.
    views_cfg = cfg.views
    views = {
        "overview": True,
        "gene_explorer": gene_idx["n_gene_models"] > 0,
        "gene_models": gene_idx["n_gene_models"] > 0,
        "domain_architecture": domain_idx["available"],
        "exon_domain_boundaries": boundary_idx["available"],
        "synteny": synteny_idx["available"],
        "figure_gallery": bool(figures_idx.get("figures")),
        # event-specific
        "event_region": has_event,
        "boundary_relation": has_event,
    }
    for k in list(views.keys()):
        if k in views_cfg and not views_cfg[k]:
            views[k] = False

    dataset_summary = {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "dataset_kind": "core_only",
        "analysis_id": cfg.analysis_id,
        "gene_symbol": cfg.gene_symbol,
        "support_level": cfg.support_level,
        "experimental": cfg.experimental,
        "analysis_modes": cfg.analysis_modes,
        "event_analysis_mode": cfg.event_analysis_mode,
        "has_event": has_event,
        "event_message": ("" if has_event else
                          "No event region is configured for this gene. "
                          "Core gene-level analysis is available."),
        "event_id": cfg.event_id,
        "event_type": cfg.event_type,
        "analysed_species_count": gene_idx["n_species"],
        # protein_isoform_count == all coding isoforms; primary_protein_count ==
        # the isoforms selected as primary for annotation (Part 2 cards).
        "protein_isoform_count": gene_idx["n_protein_isoforms"],
        "primary_protein_count": gene_idx["n_primary_proteins"],
        "gene_model_count": gene_idx["n_gene_models"],
        "synteny_neighbours_count": synteny_idx.get("n_resolved_neighbors", 0),
        "exploratory_event_evidence": "available" if n_evidence > 0 else "none",
        "candidate_clusters_count": n_clusters,
        "capability": capability,
        "status": "results_ready" if boundary_idx["available"] or domain_idx["available"] else "in_progress",
        "available_views": views,
        "ui_labels": cfg.ui_labels,
    }

    # overview_index.json (Part 8): the Overview page reads a coherent, explicit
    # index (KPIs + pending flags) rather than assembling ad-hoc from raw tables.
    overview_index = {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": cfg.analysis_id,
        "gene_symbol": cfg.gene_symbol,
        # Which assembly locus the requested symbol reached, and how. Present only when
        # the annotation names the gene differently, e.g. HBA -> LOC122209636.
        "gene_identity": _gene_identity_for(src),
        "kpis": {
            "species_analysed": gene_idx["n_species"],
            "gene_models": gene_idx["n_gene_models"],
            "protein_isoforms": gene_idx["n_protein_isoforms"],
            "primary_proteins": gene_idx["n_primary_proteins"],
            "synteny_neighbours": synteny_idx.get("n_resolved_neighbors", 0),
            "exploratory_event_candidates": n_clusters,
        },
        "domain_annotation": "available" if domain_available else "pending_cluster",
        "exon_domain_boundaries": "available" if boundary_idx.get("available") else "pending_cluster",
        "evidence_summary": [
            {"id": it.get("id"), "title": it.get("title"), "status": it.get("status"),
             "confidence": it.get("confidence")}
            for it in (evidence_stack.get("items") or [])
        ],
        "has_event": has_event,
        "available_views": views,
    }

    outputs = {
        "overview_index.json": overview_index,
        "dataset_summary.json": dataset_summary,
        "gene_analysis_index.json": gene_idx,
        "gene_event_index.json": gene_idx,   # alias for consumers expecting this name
        "gene_explorer_index.json": gene_idx,  # FGFR2-parity name (Part 9)
        "evidence_stack.json": evidence_stack,
        "primary_selection_index.json": sel,
        "protein_architecture_index.json": arch_idx,
        "event_evidence_index.json": event_evidence_idx,
        "figures_index.json": figures_idx,
        "domain_architecture_index.json": domain_idx,
        "synteny_index.json": synteny_idx,
        "exon_domain_boundary_index.json": boundary_idx,
        "exon_domain_boundaries_index.json": boundary_idx,  # FGFR2-parity name (Part 9)
        "event_region_index.json": _event_region_index(src, has_event),
        "gene_capability_report.json": capability,
        "available_views.json": {
            "schema_version": SCHEMA_VERSION,
            "dataset_id": src.dataset_id,
            "analysis_id": cfg.analysis_id,
            "has_event": has_event,
            "available_views": views,
        },
    }
    for name, data in outputs.items():
        write_json(out_dir / name, data)
    return {"out_dir": str(out_dir), "files": sorted(outputs.keys()),
            "analysis_id": cfg.analysis_id, "dataset_id": src.dataset_id,
            "has_event": has_event, "support_level": cfg.support_level,
            "available_views": views}


def _assert_not_in_freeze(out_dir: Path) -> None:
    if str(out_dir.resolve()).startswith(str(FREEZE_ROOT.resolve())):
        raise SystemExit(f"Refusing to write inside the example freeze: {out_dir}.")


def _load_cfg(config_arg: Optional[str], run_dir: Optional[Path]) -> GeneConfig:
    if config_arg:
        try:
            return load_gene_config(config_arg)
        except GeneConfigError:
            return load_gene_config_lenient(config_arg)
    if run_dir is not None:
        rc = read_json(run_dir / "run_config.json", {}) or {}
        return resolve_run_analysis(rc, run_dir)
    return default_gene_config()


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Build generic website indices from CORE gene-analysis contract outputs.")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--run-id", help="Run id under runs/ (reads results/core_gene_analysis/).")
    g.add_argument("--core-dir", help="Explicit directory with core contract TSVs (e.g. a mock).")
    ap.add_argument("--config", help="Gene config path (required for --core-dir; default resolves).")
    ap.add_argument("--dataset-id", help="Dataset id label for the indices.")
    ap.add_argument("--out", help="Output directory for generic indices.")
    args = ap.parse_args(argv)

    if args.run_id:
        try:
            record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
        except RegistryError as exc:
            raise SystemExit(str(exc)) from exc
        if record is None:
            raise SystemExit(f"Run not found: {args.run_id}")
        if record.read_only:
            raise SystemExit(
                "Run is registered read-only; copy it before rebuilding indices.")
        run_dir = record.path
        LegacyRunAdapter(run_dir, expected_run_id=args.run_id)
        core_dir = run_dir / "results" / "core_gene_analysis"
        cfg = _load_cfg(args.config, run_dir)
        dataset_id = args.dataset_id or f"run:{args.run_id}"
        out_dir = Path(args.out) if args.out else (run_dir / "website_indices" / "generic")
    else:
        core_dir = Path(args.core_dir)
        if not core_dir.is_absolute():
            core_dir = PROJECT_ROOT / core_dir
        cfg = _load_cfg(args.config, None)
        dataset_id = args.dataset_id or f"core:{core_dir.name}"
        if not args.out:
            raise SystemExit("--out is required with --core-dir.")
        out_dir = Path(args.out)

    if not out_dir.is_absolute():
        out_dir = PROJECT_ROOT / out_dir
    _assert_not_in_freeze(out_dir)
    if not core_dir.is_dir():
        raise SystemExit(f"Core outputs directory not found: {core_dir}. "
                         "Run a core adapter/runner first.")

    src = CoreSource(core_dir, dataset_id, cfg)
    result = build_core_indices(src, out_dir)

    print(f"OK  core indices  dataset={dataset_id}  analysis={cfg.analysis_id}  "
          f"support_level={cfg.support_level}")
    print(f"    has_event={result['has_event']}")
    print(f"    out: {out_dir}")
    for f in result["files"]:
        print(f"      - {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
