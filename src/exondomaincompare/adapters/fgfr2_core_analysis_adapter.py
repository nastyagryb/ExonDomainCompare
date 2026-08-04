#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from exondomaincompare.adapters.fgfr2_to_generic_indices import (  # noqa: E402
    DatasetSource, PROJECT_ROOT, FREEZE_ROOT, write_json, _int_or_none,
)
from exondomaincompare.framework.gene_config import (  # noqa: E402
    GeneConfig, load_core_analysis_contract,
)

SOURCE_LABEL = "fgfr2_core_adapter"
DEFAULT_NEAR_THRESHOLD_AA = 5


def _near_threshold() -> int:
    c = load_core_analysis_contract()
    v = _int_or_none(c.get("near_boundary_threshold_aa"))
    return v if v is not None else DEFAULT_NEAR_THRESHOLD_AA


def _tsv_write(path: Path, columns: List[str], rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in columns})


def _panels(src: DatasetSource) -> List[Dict[str, Any]]:
    sda = src.idx("species_domain_architecture.json", {}) or {}
    out: List[Dict[str, Any]] = []
    for sp in (sda.get("species", []) or []):
        sp_id = sp.get("species", "")
        panels = sp.get("panels", {}) or {}
        iterable = panels.values() if isinstance(panels, dict) else panels
        for panel in iterable:
            if isinstance(panel, dict):
                panel = dict(panel)
                panel.setdefault("species", sp_id)
                out.append(panel)
    return out


def _classify_boundary(pos: int, domains: List[Dict[str, Any]],
                       threshold: int) -> Tuple[Optional[Dict[str, Any]], str, str, Optional[int]]:
    if not domains:
        return None, "", "outside_annotated_domain", None
    best = None
    best_type = ""
    best_dist = None
    for d in domains:
        for kind, coord in (("domain_start", d.get("start")), ("domain_end", d.get("end"))):
            c = _int_or_none(coord)
            if c is None:
                continue
            dist = abs(pos - c)
            if best_dist is None or dist < best_dist:
                best_dist, best, best_type = dist, d, kind
    inside = any(
        _int_or_none(d.get("start")) is not None and _int_or_none(d.get("end")) is not None
        and _int_or_none(d.get("start")) < pos < _int_or_none(d.get("end"))
        for d in domains)
    if best_dist == 0:
        category = "exactly_aligned"
    elif best_dist is not None and best_dist <= threshold:
        category = "near_boundary"
    elif inside:
        category = "inside_domain"
    else:
        category = "outside_annotated_domain"
    return best, best_type, category, best_dist


def build_core_outputs(src: DatasetSource, cfg: GeneConfig, out_dir: Path) -> Dict[str, Any]:
    species_index = src.idx("species_index.json", []) or []
    panels = _panels(src)
    synteny = src.idx("synteny_locus_index.json", {}) or {}
    threshold = _near_threshold()

    gene_models: List[Dict[str, Any]] = []
    isoforms: List[Dict[str, Any]] = []
    exon_map: List[Dict[str, Any]] = []
    domain_rows: List[Dict[str, Any]] = []
    tm_rows: List[Dict[str, Any]] = []
    synteny_rows: List[Dict[str, Any]] = []
    boundary_rows: List[Dict[str, Any]] = []
    warnings: List[str] = []

    n_species = 0
    n_proteins = 0

    for sp in species_index if isinstance(species_index, list) else []:
        n_species += 1
        sp_id = sp.get("species", "")
        for iso in sp.get("isoforms", []) or []:
            n_proteins += 1
            protein_id = iso.get("protein_id", "")
            gene_models.append({
                "analysis_id": cfg.analysis_id, "gene_symbol": cfg.gene_symbol,
                "species_id": sp_id, "gene_id": iso.get("gene_id", ""),
                "transcript_id": iso.get("transcript_id", ""), "protein_id": protein_id,
                "source": SOURCE_LABEL,
                "protein_length": _int_or_none(iso.get("protein_length")),
                "model_status": iso.get("final_claim_status_after_rescue", ""),
                "notes": "",
            })
            isoforms.append({
                "species_id": sp_id, "protein_id": protein_id,
                "transcript_id": iso.get("transcript_id", ""),
                "isoform_label": iso.get("final_isoform_label", "") or iso.get("isoform", ""),
                "protein_length": _int_or_none(iso.get("protein_length")),
                "sequence_path": "",
                "primary_status": iso.get("interpro_included", ""),
                "notes": "",
            })

    # Exon map, domain/TM features, and generic boundary distances from panels.
    for panel in panels:
        sp_id = panel.get("species", "")
        protein_id = panel.get("protein_id", "")
        transcript_id = panel.get("transcript_id", "")
        domains = [d for d in (panel.get("domains", []) or []) if isinstance(d, dict)]
        for d in domains:
            domain_rows.append({
                "species_id": sp_id, "protein_id": protein_id,
                "domain_source": d.get("source", ""),
                "domain_id": d.get("class", ""), "domain_name": d.get("label", ""),
                "start_aa": _int_or_none(d.get("start")), "end_aa": _int_or_none(d.get("end")),
                "score": d.get("score", ""),
            })
        for t in (panel.get("tm", []) or []):
            if not isinstance(t, dict):
                continue
            tm_rows.append({
                "species_id": sp_id, "protein_id": protein_id,
                "start_aa": _int_or_none(t.get("start")), "end_aa": _int_or_none(t.get("end")),
                "source": t.get("source", ""),
            })
        exons = [e for e in (panel.get("exons", []) or []) if isinstance(e, dict)]
        for e in exons:
            num = e.get("number")
            start = _int_or_none(e.get("start"))
            end = _int_or_none(e.get("end"))
            exon_map.append({
                "species_id": sp_id, "protein_id": protein_id, "transcript_id": transcript_id,
                "exon_id": e.get("label", "") or (f"exon{num}" if num is not None else ""),
                "exon_number": num, "cds_start": "", "cds_end": "",
                "protein_start_aa": start, "protein_end_aa": end,
                "phase": "", "confidence": "", "source": "species_domain_architecture",
            })
            # Generic exon-domain boundary distance at the exon C-terminal boundary.
            if end is None:
                continue
            nearest, btype, category, dist = _classify_boundary(end, domains, threshold)
            boundary_rows.append({
                "analysis_id": cfg.analysis_id, "gene_symbol": cfg.gene_symbol,
                "species_id": sp_id, "protein_id": protein_id, "transcript_id": transcript_id,
                "exon_boundary_id": f"{protein_id}:exon{num}_end" if num is not None else f"{protein_id}:{end}",
                "boundary_position_aa": end,
                "nearest_domain_id": (nearest or {}).get("class", ""),
                "nearest_domain_name": (nearest or {}).get("label", ""),
                "nearest_domain_boundary_type": btype,
                "distance_aa": dist, "category": category, "source": SOURCE_LABEL,
            })

    if not panels:
        warnings.append("species_domain_architecture.json has no panels; exon/domain/TM "
                        "and boundary outputs are empty (domain annotation not available).")

    # Synteny neighbours (optional).
    n_resolved = 0
    for sp in (synteny.get("species", []) or []):
        sp_id = sp.get("species", "")
        neigh = sp.get("neighbors10") or sp.get("neighbors5") or sp.get("neighbors") or []
        for nb in neigh:
            if not isinstance(nb, dict):
                continue
            resolved = bool(nb.get("resolved"))
            n_resolved += 1 if resolved else 0
            synteny_rows.append({
                "species_id": sp_id, "gene_symbol": cfg.gene_symbol,
                "neighbor_symbol": nb.get("symbol", ""),
                "side": nb.get("side", ""), "order": nb.get("rank", ""),
                "orientation": nb.get("strand", ""),
                "classification": nb.get("identity_status", "") or nb.get("method_class", ""),
                "source": nb.get("method", "") or "synteny_locus_index",
                "status": "resolved" if resolved else "unresolved",
            })
    if not synteny.get("available"):
        warnings.append("synteny not available for this dataset; synteny_neighbors is empty.")

    _tsv_write(out_dir / "gene_model_index.tsv",
               ["analysis_id", "gene_symbol", "species_id", "gene_id", "transcript_id",
                "protein_id", "source", "protein_length", "model_status", "notes"], gene_models)
    _tsv_write(out_dir / "protein_isoform_index.tsv",
               ["species_id", "protein_id", "transcript_id", "isoform_label", "protein_length",
                "sequence_path", "primary_status", "notes"], isoforms)
    _tsv_write(out_dir / "exon_protein_map.tsv",
               ["species_id", "protein_id", "transcript_id", "exon_id", "exon_number",
                "cds_start", "cds_end", "protein_start_aa", "protein_end_aa", "phase",
                "confidence", "source"], exon_map)
    _tsv_write(out_dir / "domain_features.tsv",
               ["species_id", "protein_id", "domain_source", "domain_id", "domain_name",
                "start_aa", "end_aa", "score"], domain_rows)
    _tsv_write(out_dir / "tm_features.tsv",
               ["species_id", "protein_id", "start_aa", "end_aa", "source"], tm_rows)
    _tsv_write(out_dir / "synteny_neighbors.tsv",
               ["species_id", "gene_symbol", "neighbor_symbol", "side", "order", "orientation",
                "classification", "source", "status"], synteny_rows)
    _tsv_write(out_dir / "exon_domain_boundary_distances.tsv",
               ["analysis_id", "gene_symbol", "species_id", "protein_id", "transcript_id",
                "exon_boundary_id", "boundary_position_aa", "nearest_domain_id",
                "nearest_domain_name", "nearest_domain_boundary_type", "distance_aa",
                "category", "source"], boundary_rows)

    # Available views: core views are gene-agnostic; event views only if configured.
    has_domains = bool(domain_rows)
    has_exons = bool(exon_map)
    has_msa = (src.closure / "MSA" / "final_cassette_msa_boundary_projection.tsv").exists()
    available_views = {
        "overview": True,
        "gene_models": bool(gene_models),
        "isoforms": bool(isoforms),
        "msa": has_msa,
        "domain_architecture": has_domains,
        "synteny": bool(synteny.get("available")) and n_resolved > 0,
        "exon_domain_boundaries": has_domains and has_exons,
        # event-specific views require a configured event region
        "event_region": cfg.has_event,
        "event_specific_boundary_relation": cfg.has_event,
    }

    report = {
        "analysis_id": cfg.analysis_id,
        "gene_symbol": cfg.gene_symbol,
        "dataset_id": src.dataset_id,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "adapter": SOURCE_LABEL,
        "contract_version": int((load_core_analysis_contract().get("contract_version", 1)) or 1),
        "event_analysis_mode": cfg.event_analysis_mode,
        "has_event": cfg.has_event,
        "summary": {
            "n_species": n_species,
            "n_proteins": n_proteins,
            "n_domain_features": len(domain_rows),
            "n_tm_features": len(tm_rows),
            "n_exon_boundaries": len(boundary_rows),
            "n_synteny_neighbors": len(synteny_rows),
            "n_resolved_synteny_neighbors": n_resolved,
        },
        "inputs": [
            "website_indices/species_index.json",
            "website_indices/species_domain_architecture.json",
            "website_indices/synteny_locus_index.json",
        ],
        "outputs": [
            "gene_model_index.tsv", "protein_isoform_index.tsv", "exon_protein_map.tsv",
            "domain_features.tsv", "tm_features.tsv", "synteny_neighbors.tsv",
            "exon_domain_boundary_distances.tsv", "core_gene_report.json",
        ],
        "warnings": sorted(set(warnings)),
        "failures": [],
        "available_views": available_views,
        "note": "Projection of existing validated FGFR2 outputs into the gene-agnostic "
                "core contract; only boundary geometry is computed (no biology recomputed).",
    }
    write_json(out_dir / "core_gene_report.json", report)
    return report


def _resolve_out_dir(src: DatasetSource, out_arg: Optional[str]) -> Path:
    if out_arg:
        out = Path(out_arg)
    elif src.kind == "run":
        out = src.run_root / "results" / "core_gene_analysis"
    else:
        out = PROJECT_ROOT / "artifacts" / "core_gene_analysis" / "example"
    if not out.is_absolute():
        out = PROJECT_ROOT / out
    if str(out.resolve()).startswith(str(FREEZE_ROOT.resolve())):
        raise SystemExit(
            f"Refusing to write core outputs inside the example freeze: {out}. "
            "Use --out (e.g. artifacts/core_gene_analysis/example).")
    return out


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Project validated FGFR2 outputs into the generic CORE gene-analysis contract.")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--run-id", help="Custom run id under runs/.")
    g.add_argument("--example", action="store_true", help="Read the example dataset (read-only).")
    ap.add_argument("--config", help="Override gene_config.yaml path.")
    ap.add_argument("--out", help="Output directory for core contract outputs.")
    args = ap.parse_args(argv)

    src = DatasetSource.for_example() if args.example else DatasetSource.for_run(args.run_id)
    cfg = src.gene_config(args.config)
    out_dir = _resolve_out_dir(src, args.out)
    report = build_core_outputs(src, cfg, out_dir)

    s = report["summary"]
    print(f"OK  core adapter  dataset={src.dataset_id}  analysis={cfg.analysis_id}")
    print(f"    species={s['n_species']} proteins={s['n_proteins']} domains={s['n_domain_features']} "
          f"tm={s['n_tm_features']} exon_boundaries={s['n_exon_boundaries']} "
          f"synteny={s['n_resolved_synteny_neighbors']}/{s['n_synteny_neighbors']}")
    print(f"    has_event={report['has_event']} (mode={report['event_analysis_mode']})")
    print(f"    out: {out_dir}")
    for f in report["outputs"]:
        print(f"      - {f}")
    if report["warnings"]:
        print(f"    warnings: {len(report['warnings'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
