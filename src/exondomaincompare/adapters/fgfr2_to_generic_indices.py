#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from exondomaincompare.framework.gene_config import (  # noqa: E402
    GeneConfig, load_gene_config, resolve_run_analysis, default_gene_config,
)
from exondomaincompare.runs.legacy import LegacyRunAdapter  # noqa: E402
from exondomaincompare.runs.registry import RegistryError, resolve_run_record  # noqa: E402
from exondomaincompare.config import discover_repository_root, load_config  # noqa: E402

PROJECT_ROOT = discover_repository_root(__file__)
RUNTIME_CONFIG = load_config(repository_root=PROJECT_ROOT)
FREEZE_ROOT = PROJECT_ROOT / "results" / "final_30_until_interpro_prepare"
EXAMPLE_CLOSURE = FREEZE_ROOT / "13_final_pre_interpro_closure"

SCHEMA_VERSION = 1


# --------------------------------------------------------------------------- #
# small IO helpers
# --------------------------------------------------------------------------- #
def read_json(p: Path, default: Any = None) -> Any:
    try:
        return json.loads(Path(p).read_text(encoding="utf-8"))
    except Exception:
        return default


def read_tsv(p: Path) -> List[Dict[str, str]]:
    if not Path(p).is_file():
        return []
    with open(p, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


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
# dataset location resolution
# --------------------------------------------------------------------------- #
class DatasetSource:
    def __init__(self, *, kind: str, dataset_id: str, closure: Path,
                 website_indices: Path, run_config: Dict[str, Any],
                 status: Dict[str, Any], run_root: Optional[Path]):
        self.kind = kind                      # "run" | "example"
        self.dataset_id = dataset_id
        self.closure = closure
        self.website_indices = website_indices
        self.run_config = run_config
        self.status = status
        self.run_root = run_root

    @classmethod
    def for_run(cls, run_id: str) -> "DatasetSource":
        try:
            record = resolve_run_record(RUNTIME_CONFIG, run_id)
        except RegistryError as exc:
            raise SystemExit(str(exc)) from exc
        if record is None:
            raise SystemExit(f"Run not found: {run_id}")
        run_root = record.path
        adapter = LegacyRunAdapter(run_root, expected_run_id=run_id)
        return cls(
            kind="run",
            dataset_id=f"run:{run_id}",
            closure=run_root / "results" / "13_final_pre_interpro_closure",
            website_indices=run_root / "website_indices",
            run_config=adapter.config(),
            status=adapter.status(),
            run_root=run_root,
        )

    @classmethod
    def for_example(cls) -> "DatasetSource":
        return cls(
            kind="example",
            dataset_id="example",
            closure=EXAMPLE_CLOSURE,
            website_indices=EXAMPLE_CLOSURE / "website_indices",
            run_config={},   # example has no run_config.json -> FGFR2 default
            status={"status": "results_ready"},
            run_root=None,
        )

    def gene_config(self, override: Optional[str]) -> GeneConfig:
        if override:
            return load_gene_config(override)
        if self.kind == "run":
            return resolve_run_analysis(self.run_config, self.run_root)
        return default_gene_config()

    def idx(self, name: str, default: Any = None) -> Any:
        return read_json(self.website_indices / name, default)


# --------------------------------------------------------------------------- #
# generic index builders
# --------------------------------------------------------------------------- #
def _projection_coords(closure: Path) -> Dict[tuple, Dict[str, Optional[int]]]:
    proj = closure / "MSA" / "final_cassette_msa_boundary_projection.tsv"
    out: Dict[tuple, Dict[str, Optional[int]]] = {}
    for row in read_tsv(proj):
        sp = (row.get("species") or "").strip()
        iso = (row.get("isoform") or "").strip()
        if not sp or not iso:
            continue
        start = _int_or_none(row.get("native_cassette_start_aa"))
        end = _int_or_none(row.get("native_cassette_end_aa"))
        out[(sp, iso)] = {"region_start_aa": start, "region_end_aa": end}
    return out


def _available_views(src: DatasetSource, cfg: GeneConfig) -> Dict[str, bool]:
    def idx_available(name: str) -> bool:
        data = src.idx(name)
        if data is None:
            return False
        if isinstance(data, dict) and "available" in data:
            return bool(data["available"])
        return True

    run_index = src.idx("run_index.json")
    truth_table = (src.closure / "final_pre_interpro_truth_table.tsv").exists()
    views_cfg = cfg.views
    has_event = cfg.has_event
    domain_ok = idx_available("species_domain_architecture.json")
    computed = {
        # core (gene-agnostic) views — never require a configured event region
        "overview": bool(run_index) or truth_table,
        "gene_explorer": bool(src.idx("species_index.json")),
        "domain_architecture": domain_ok,
        "exon_domain_boundaries": domain_ok,
        "synteny": idx_available("synteny_locus_index.json"),
        # event-specific views — only when an event region is configured
        "event_region": has_event and idx_available("cassette_residue_index.json"),
        "boundary_relation": has_event and idx_available("boundary_consistency_summary.json"),
    }
    # A view is shown only if the gene config enables it AND the data exists.
    # (event_region/boundary_relation map from config keys of the same name.)
    for k in list(computed.keys()):
        if k in views_cfg and not views_cfg[k]:
            computed[k] = False
    return computed


def build_dataset_summary(src: DatasetSource, cfg: GeneConfig,
                          views: Dict[str, bool]) -> Dict[str, Any]:
    species = src.idx("species_index.json", []) or []
    analysed_count = len(species) if isinstance(species, list) else 0
    primary_count = _int_or_none(src.run_config.get("primary_fasta_count"))
    if primary_count is None:
        # fall back to counting primary isoforms in the species index
        primary_count = 0
        for sp in species if isinstance(species, list) else []:
            for iso in sp.get("isoforms", []) or []:
                if str(iso.get("interpro_included", "")).lower() == "primary":
                    primary_count += 1
    # Prefer the real, view-derived readiness over a possibly-stale raw status
    # field: once the boundary_relation (post-InterPro) view exists, the dataset
    # is results_ready regardless of leftover intermediate status text.
    raw_status = str(src.status.get("status") or "").strip()
    if views.get("boundary_relation"):
        status = "results_ready"
    elif raw_status:
        status = raw_status
    else:
        status = "in_progress"
    rc = cfg.reference_control
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "dataset_kind": src.kind,
        "analysis_id": cfg.analysis_id,
        "gene_symbol": cfg.gene_symbol,
        "analysis_modes": cfg.analysis_modes,
        "event_analysis_mode": cfg.event_analysis_mode,
        "has_event": cfg.has_event,
        "event_id": cfg.event_id,
        "event_type": cfg.event_type,
        "event_display_name": cfg.event_display_name,
        "analysed_species_count": analysed_count,
        "primary_protein_count": primary_count,
        "reference_control": {
            "enabled": rc["enabled"],
            "species": rc["species"],
            "role": rc["role"],
        },
        "status": status,
        "available_views": views,
        "ui_labels": cfg.ui_labels,
    }


def build_gene_event_index(src: DatasetSource, cfg: GeneConfig) -> Dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": cfg.analysis_id,
        "analysis_display_name": cfg.analysis_display_name,
        "analysis_description": cfg.analysis_description,
        "gene": {
            "symbol": cfg.gene_symbol,
            "display_name": cfg.gene_display_name,
            "reference_species": cfg.reference_species,
        },
        "event": {
            "id": cfg.event_id,
            "type": cfg.event_type,
            "display_name": cfg.event_display_name,
            "generic_label": cfg.event_generic_label,
            "labels": cfg.event_labels,
        },
        "reference_control": cfg.reference_control,
    }


def _detector_dir(src: DatasetSource) -> Optional[Path]:
    if src.run_root is not None:
        d = src.run_root / "results" / "generic_event_detector"
        if (d / "event_region_coordinates.tsv").is_file():
            return d
    return None


def _event_region_from_detector(src: DatasetSource, cfg: GeneConfig,
                                det_dir: Path) -> List[Dict[str, Any]]:
    regions = read_tsv(det_dir / "event_region_coordinates.tsv")
    candidates = read_tsv(det_dir / "event_isoform_candidates.tsv")
    cand_by_key = {(r.get("species_id", ""), r.get("protein_id", ""), r.get("event_label", "")): r
                   for r in candidates}
    display = {}
    for sp in (src.idx("species_index.json", []) or []):
        display[sp.get("species", "")] = sp.get("display_species_name", sp.get("species", ""))

    by_species: Dict[str, Dict[str, Any]] = {}
    for r in regions:
        sp_id = r.get("species_id", "")
        if not sp_id:
            continue
        label = r.get("event_label", "")
        cand = cand_by_key.get((sp_id, r.get("protein_id", ""), label), {})
        node = by_species.setdefault(sp_id, {
            "species_id": sp_id,
            "display_name": display.get(sp_id, sp_id),
            "records": [],
        })
        node["records"].append({
            "event_label": label,
            "protein_id": r.get("protein_id", ""),
            "transcript_id": cand.get("transcript_id", ""),
            "protein_length": _int_or_none(cand.get("protein_length")),
            "region_start_aa": _int_or_none(r.get("region_start_aa")),
            "region_end_aa": _int_or_none(r.get("region_end_aa")),
            "region_length_aa": _int_or_none(r.get("region_length_aa")),
            "validation_status": cand.get("candidate_status", ""),
        })
    return list(by_species.values())


def build_event_region_index(src: DatasetSource, cfg: GeneConfig) -> Dict[str, Any]:
    label_ids = cfg.event_label_ids

    # Preferred path: build from generic event-detector contract outputs if present.
    det_dir = _detector_dir(src)
    if det_dir is not None:
        species_out = _event_region_from_detector(src, cfg, det_dir)
        if species_out:
            return _event_region_payload(src, cfg, label_ids, species_out,
                                         source="generic_event_detector")

    # Fallback path: read the FGFR2-specific website indices + closure projection.
    species_index = src.idx("species_index.json", []) or []
    coords = _projection_coords(src.closure)

    species_out: List[Dict[str, Any]] = []
    for sp in species_index if isinstance(species_index, list) else []:
        sp_id = sp.get("species", "")
        records: List[Dict[str, Any]] = []
        for iso in sp.get("isoforms", []) or []:
            iso_label = iso.get("final_isoform_label") or iso.get("isoform") or ""
            c = coords.get((sp_id, iso.get("isoform", "")), {})
            start = c.get("region_start_aa")
            end = c.get("region_end_aa")
            length = (end - start + 1) if (isinstance(start, int) and isinstance(end, int)) else None
            records.append({
                "event_label": iso_label,
                "protein_id": iso.get("protein_id", ""),
                "transcript_id": iso.get("transcript_id", ""),
                "protein_length": _int_or_none(iso.get("protein_length")),
                "region_start_aa": start,
                "region_end_aa": end,
                "region_length_aa": length,
                "validation_status": iso.get("final_claim_status_after_rescue")
                                     or iso.get("readiness_status", {}).get("value", ""),
            })
        species_out.append({
            "species_id": sp_id,
            "display_name": sp.get("display_species_name", sp_id),
            "records": records,
        })

    return _event_region_payload(src, cfg, label_ids, species_out,
                                 source="fgfr2_website_indices")


def _event_region_payload(src: DatasetSource, cfg: GeneConfig, label_ids: List[str],
                          species_out: List[Dict[str, Any]], source: str) -> Dict[str, Any]:
    # Reference comparison: human curated IIIb/IIIc, taken from the cassette index's
    # human_reference block if present (reference/control only, never analysed).
    cassette = src.idx("cassette_residue_index.json", {}) or {}
    human_ref = cassette.get("human_reference") or {}
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "event_id": cfg.event_id,
        "event_type": cfg.event_type,
        "event_labels": label_ids,
        "region_label": cfg.event_display_name,
        "source": source,
        "species": species_out,
        "reference_comparison": {
            "enabled": bool(cfg.reference_control["enabled"] and human_ref),
            "reference_species": cfg.reference_species,
            "basis": "msa_column",
            "role": cfg.reference_control["role"],
            "rows": [],
        },
    }


def build_domain_architecture_index(src: DatasetSource, cfg: GeneConfig) -> Dict[str, Any]:
    summary = src.idx("domain_architecture_summary.json", {}) or {}
    per_species = src.idx("species_domain_architecture.json", {}) or {}
    available = bool(summary.get("available")) if isinstance(summary, dict) else False
    if isinstance(per_species, dict) and "available" in per_species:
        available = available or bool(per_species.get("available"))
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": cfg.analysis_id,
        "available": available,
        "relation_description": cfg.ui_labels.get("domain_relation_description", ""),
        "summary": summary if isinstance(summary, dict) else {},
        "source_indices": ["domain_architecture_summary.json", "species_domain_architecture.json"],
    }


def build_synteny_index(src: DatasetSource, cfg: GeneConfig) -> Dict[str, Any]:
    syn = src.idx("synteny_locus_index.json", {}) or {}
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": cfg.analysis_id,
        "available": bool(syn.get("available")),
        "synteny_status": syn.get("synteny_status", "not_computed"),
        "synteny_reason": syn.get("synteny_reason", ""),
        "n_resolved_neighbors": syn.get("n_resolved_neighbors", 0),
        "reference_comparison": {
            "enabled": bool(syn.get("human_reference_available")),
            "reference_species": cfg.reference_species,
            "role": syn.get("human_reference_role", ""),
        },
        "species": syn.get("species", []),
        "source_indices": ["synteny_locus_index.json"],
    }


def build_boundary_relation_index(src: DatasetSource, cfg: GeneConfig) -> Dict[str, Any]:
    summary = src.idx("boundary_consistency_summary.json", {}) or {}
    available = bool(summary.get("available")) if isinstance(summary, dict) else False
    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": src.dataset_id,
        "analysis_id": cfg.analysis_id,
        "event_id": cfg.event_id,
        "available": available,
        "relation_description": cfg.ui_labels.get("domain_relation_description", ""),
        "relation_label": cfg.ui_labels.get("boundary_relation", "Boundary Consistency"),
        "boundary_class_counts": summary.get("boundary_class_counts", {}) if isinstance(summary, dict) else {},
        "median_distance": summary.get("median_distance") if isinstance(summary, dict) else None,
        "mean_distance": summary.get("mean_distance") if isinstance(summary, dict) else None,
        "n_cassette_boundaries": summary.get("n_cassette_boundaries") if isinstance(summary, dict) else None,
        "source_indices": ["boundary_consistency_summary.json", "boundary_consistency_index.json"],
    }


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def generate_generic_indices(src: DatasetSource, out_dir: Path,
                             config_override: Optional[str] = None) -> Dict[str, Any]:
    cfg = src.gene_config(config_override)
    views = _available_views(src, cfg)

    outputs = {
        "dataset_summary.json": build_dataset_summary(src, cfg, views),
        "gene_event_index.json": build_gene_event_index(src, cfg),
        "event_region_index.json": build_event_region_index(src, cfg),
        "domain_architecture_index.json": build_domain_architecture_index(src, cfg),
        "synteny_index.json": build_synteny_index(src, cfg),
        "boundary_relation_index.json": build_boundary_relation_index(src, cfg),
        "available_views.json": {
            "schema_version": SCHEMA_VERSION,
            "dataset_id": src.dataset_id,
            "analysis_id": cfg.analysis_id,
            "available_views": views,
        },
    }
    for name, data in outputs.items():
        write_json(out_dir / name, data)
    return {"out_dir": str(out_dir), "files": sorted(outputs.keys()),
            "analysis_id": cfg.analysis_id, "dataset_id": src.dataset_id}


def _assert_not_in_freeze(out_dir: Path) -> None:
    resolved = out_dir.resolve()
    if str(resolved).startswith(str(FREEZE_ROOT.resolve())):
        raise SystemExit(
            f"Refusing to write generic indices inside the example freeze: {resolved}. "
            "Use --out to choose a safe location (e.g. artifacts/generic_indices/example).")


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="FGFR2 -> generic website-index adapter (additive).")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--run-id", help="Custom run id under runs/.")
    g.add_argument("--example", action="store_true",
                   help="Read the example dataset (freeze) read-only.")
    ap.add_argument("--config", help="Override gene_config.yaml path (default: resolve/FGFR2).")
    ap.add_argument("--out", help="Output directory for generic indices.")
    args = ap.parse_args(argv)

    if args.example:
        src = DatasetSource.for_example()
        out_dir = Path(args.out) if args.out else (PROJECT_ROOT / "artifacts" / "generic_indices" / "example")
    else:
        src = DatasetSource.for_run(args.run_id)
        out_dir = Path(args.out) if args.out else (src.run_root / "website_indices" / "generic")
    if not out_dir.is_absolute():
        out_dir = PROJECT_ROOT / out_dir

    _assert_not_in_freeze(out_dir)
    if not src.website_indices.is_dir():
        print(f"WARNING: canonical website indices not found at {src.website_indices}. "
              "Generic indices will reflect only what exists.", file=sys.stderr)

    result = generate_generic_indices(src, out_dir, config_override=args.config)
    print(f"OK  dataset={result['dataset_id']}  analysis={result['analysis_id']}")
    print(f"    out: {out_dir}")
    for f in result["files"]:
        print(f"      - {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
