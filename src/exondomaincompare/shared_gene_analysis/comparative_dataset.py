#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from exondomaincompare.config import discover_repository_root, load_config
from exondomaincompare.runs.registry import RegistryError, resolve_run_record

ROOT = discover_repository_root(__file__)
from . import boundary_observations, species_order  # noqa: E402
from .msa_coordinates import (  # noqa: E402
    build_msa_coordinate_map, lookup_tables, read_aligned_fasta, parse_header,
)
RUNTIME_CONFIG = load_config(repository_root=ROOT)


def _rel(path: Path, root: Path = ROOT) -> str:
    try:
        return str(path.resolve().relative_to(root))
    except ValueError:
        return str(path)


def _write_tsv(path: Path, rows: Sequence[Dict[str, Any]],
               fieldnames: Optional[Sequence[str]] = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and not fieldnames:
        path.write_text("", encoding="utf-8")
        return
    names = list(fieldnames or rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=names, extrasaction="ignore",
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: "" if row.get(k) is None else row.get(k) for k in names})


def _load_json(path: Path) -> Any:
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _native_to_column(tables: Dict[str, Dict[int, int]], species_id: str,
                      native: Optional[int]) -> Optional[int]:
    if native is None:
        return None
    try:
        return tables.get(species_id, {}).get(int(native))
    except (TypeError, ValueError):
        return None


def _project_interval(tables: Dict[str, Dict[int, int]], species_id: str,
                      start: Optional[int], end: Optional[int]
                      ) -> Tuple[Optional[int], Optional[int], str]:
    c0 = _native_to_column(tables, species_id, start)
    c1 = _native_to_column(tables, species_id, end)
    if c0 is None and c1 is None:
        return None, None, "unmapped"
    if c0 is None or c1 is None:
        return c0, c1, "partial"
    return min(c0, c1), max(c0, c1), "mapped"


def _species_inventory(models: Sequence[Dict[str, Any]],
                       completion: Optional[Dict[str, Any]] = None
                       ) -> List[Dict[str, Any]]:
    by_status = {
        (r.get("species_id") or ""): r
        for r in (completion or {}).get("species") or []
    }
    rows = []
    for m in models:
        sid = m.get("species_id") or ""
        prov = m.get("provenance") or {}
        st = by_status.get(sid) or {}
        rows.append({
            "species_id": sid,
            "scientific_name": m.get("scientific_name") or sid,
            "taxonomic_group": prov.get("taxonomic_group") or "",
            "clade": prov.get("clade") or "",
            "protein_id": m.get("protein_id") or "",
            "transcript_id": m.get("transcript_id") or "",
            "protein_length": m.get("protein_length"),
            "n_exons": len(m.get("exons") or []),
            "n_domains": len(m.get("representative_domains") or []),
            "n_boundaries": len(m.get("exon_boundaries") or []),
            "n_transcript_models": m.get("n_transcript_models")
            or len(m.get("transcript_models") or []),
            "analysis_status": st.get("status") or m.get("status") or "unknown",
            "domain_status": st.get("domain_status") or (
                "available" if m.get("representative_domains") else "pending"),
            "boundary_status": st.get("boundary_status") or (
                "available" if m.get("exon_boundaries") else "pending"),
        })
    rows.sort(key=lambda r: (r.get("taxonomic_group") or "",
                             r.get("scientific_name") or "",
                             r.get("species_id") or ""))
    return rows


def _msa_aligned_exons(models: Sequence[Dict[str, Any]],
                       tables: Dict[str, Dict[int, int]]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for m in models:
        sid = m.get("species_id") or ""
        for ex in m.get("exons") or []:
            start = ex.get("start")
            end = ex.get("end")
            c0, c1, status = _project_interval(tables, sid, start, end)
            rows.append({
                "species_id": sid,
                "scientific_name": m.get("scientific_name") or sid,
                "protein_id": m.get("protein_id") or "",
                "exon_id": ex.get("id") or "",
                "exon_label": ex.get("label") or "",
                "native_start": start,
                "native_end": end,
                "msa_start_column": c0,
                "msa_end_column": c1,
                "msa_mapping_status": status,
            })
    return rows


def _msa_aligned_domains(models: Sequence[Dict[str, Any]],
                         tables: Dict[str, Dict[int, int]]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for m in models:
        sid = m.get("species_id") or ""
        for d in m.get("representative_domains") or []:
            start = d.get("start")
            end = d.get("end")
            c0, c1, status = _project_interval(tables, sid, start, end)
            rows.append({
                "species_id": sid,
                "scientific_name": m.get("scientific_name") or sid,
                "protein_id": m.get("protein_id") or "",
                "domain_instance_id": d.get("domain_instance_id") or d.get("id") or "",
                "interpro_accession": d.get("interpro_accession") or "",
                "label": d.get("short_label") or d.get("label") or "",
                "instance_number": d.get("instance_number"),
                "native_start": start,
                "native_end": end,
                "msa_start_column": c0,
                "msa_end_column": c1,
                "msa_mapping_status": status,
                "order_along_protein": d.get("display_order") or d.get("instance_number"),
            })
    return rows


def _domain_annotation_matrix(models: Sequence[Dict[str, Any]],
                              aligned_domains: Sequence[Dict[str, Any]]
                              ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    # Bucket instances by (accession, instance_number) when numbers agree across
    # species; otherwise fall back to MSA-interval proximity for the same accession.
    by_key: Dict[Tuple[str, int], List[Dict[str, Any]]] = defaultdict(list)
    for row in aligned_domains:
        acc = row.get("interpro_accession") or ""
        num = int(row.get("instance_number") or 0)
        if not acc or not num:
            continue
        by_key[(acc, num)].append(row)

    groups: List[Dict[str, Any]] = []
    matrix_rows: List[Dict[str, Any]] = []
    gid = 0
    for (acc, num), members in sorted(by_key.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        if len({m["species_id"] for m in members}) < 1:
            continue
        gid += 1
        group_id = f"CDG{gid}"
        label = next((m.get("label") for m in members if m.get("label")), acc)
        groups.append({
            "comparable_domain_group_id": group_id,
            "interpro_accession": acc,
            "instance_number": num,
            "label": label,
            "n_species": len({m["species_id"] for m in members}),
            "mapping_status": (
                "supported" if len({m["species_id"] for m in members}) >= 2
                else "single_species_only"
            ),
            "members": members,
        })
        by_sp = {m["species_id"]: m for m in members}
        for m in models:
            sid = m.get("species_id") or ""
            hit = by_sp.get(sid)
            if hit:
                state = "detected" if hit.get("msa_mapping_status") == "mapped" \
                    else "uncertain mapping"
                matrix_rows.append({
                    "species_id": sid,
                    "scientific_name": m.get("scientific_name") or sid,
                    "comparable_domain_group_id": group_id,
                    "interpro_accession": acc,
                    "label": label,
                    "state": state,
                    "domain_instance_id": hit.get("domain_instance_id") or "",
                    "native_start": hit.get("native_start"),
                    "native_end": hit.get("native_end"),
                    "msa_start_column": hit.get("msa_start_column"),
                    "msa_end_column": hit.get("msa_end_column"),
                })
            else:
                # Distinguishes annotation gap from pipeline incompleteness.
                status = m.get("status") or ""
                state = "pending" if status in ("pending", "pending_cluster") else (
                    "unavailable" if status == "failed" else "not detected")
                matrix_rows.append({
                    "species_id": sid,
                    "scientific_name": m.get("scientific_name") or sid,
                    "comparable_domain_group_id": group_id,
                    "interpro_accession": acc,
                    "label": label,
                    "state": state,
                    "domain_instance_id": "",
                    "native_start": None,
                    "native_end": None,
                    "msa_start_column": None,
                    "msa_end_column": None,
                })
    return groups, matrix_rows


def _pairwise_identity(alignment_path: Path) -> List[Dict[str, Any]]:
    if not alignment_path.is_file():
        return []
    records = read_aligned_fasta(alignment_path)
    parsed = []
    for header, seq in records:
        meta = parse_header(header)
        parsed.append((meta.get("species_id") or meta.get("protein_id") or header,
                       meta.get("protein_id") or "", seq))
    rows: List[Dict[str, Any]] = []
    for i, (sa, pa, a) in enumerate(parsed):
        for sb, pb, b in parsed[i:]:
            if len(a) != len(b):
                continue
            matches = 0
            compared = 0
            for ca, cb in zip(a, b):
                if ca in "-." or cb in "-.":
                    continue
                compared += 1
                if ca.upper() == cb.upper():
                    matches += 1
            identity = (matches / compared) if compared else None
            rows.append({
                "species_a": sa, "protein_a": pa,
                "species_b": sb, "protein_b": pb,
                "n_compared_columns": compared,
                "n_identical": matches,
                "percent_identity": None if identity is None else round(100 * identity, 2),
            })
    return rows


def _isoform_diversity(models: Sequence[Dict[str, Any]],
                       run_dir: Path) -> List[Dict[str, Any]]:
    isoform_index = _load_json(
        run_dir / "website_indices" / "isoform_alignment_index.json") or {}
    by_sp = {}
    for block in isoform_index.get("species") or []:
        sid = block.get("species_id") or block.get("species") or ""
        by_sp[sid] = block

    rows = []
    for m in models:
        sid = m.get("species_id") or ""
        block = by_sp.get(sid) or {}
        models_list = block.get("proteins") or block.get("isoforms") or m.get("transcript_models") or []
        lengths = []
        n_curated = 0
        n_predicted = 0
        for p in models_list:
            if isinstance(p, dict):
                ln = p.get("protein_length") or p.get("length")
                if ln:
                    lengths.append(int(ln))
                src = str(p.get("source") or p.get("biotype") or "").lower()
                if "predict" in src or "xp_" in str(p.get("protein_id") or "").lower():
                    n_predicted += 1
                else:
                    n_curated += 1
        primary_len = m.get("protein_length")
        if primary_len and lengths:
            max_diff = max(abs(l - int(primary_len)) for l in lengths)
        else:
            max_diff = None
        cand = m.get("candidate_regions") or []
        rows.append({
            "species_id": sid,
            "scientific_name": m.get("scientific_name") or sid,
            "primary_protein_id": m.get("protein_id") or "",
            "n_protein_models": len(models_list) or m.get("n_transcript_models") or 1,
            "n_curated_models": n_curated or None,
            "n_predicted_models": n_predicted or None,
            "primary_protein_length": primary_len,
            "protein_length_min": min(lengths) if lengths else primary_len,
            "protein_length_max": max(lengths) if lengths else primary_len,
            "n_variable_alignment_blocks": block.get("n_variable_blocks"),
            "n_exploratory_candidates": len(cand),
            "max_difference_from_primary_aa": max_diff,
        })
    return rows


def _synteny_blocks(run_dir: Path,
                    models: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows = (_load_json(run_dir / "website_indices" / "synteny_locus_index.json")
            or {}).get("species") or []
    name_by_id = {m.get("species_id"): m.get("scientific_name") for m in models}
    for row in rows:
        sid = row.get("species_id") or row.get("species") or ""
        row["scientific_name"] = name_by_id.get(sid) or row.get("display_species_name") or sid
    return rows


def _synteny_rows(blocks: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for block in blocks:
        sid = block.get("species_id") or block.get("species") or ""
        for locus in block.get("loci") or []:
            rows.append({
                "species_id": sid,
                "scientific_name": block.get("scientific_name") or sid,
                "target_gene": block.get("target_symbol") or "",
                "slot_x": locus.get("slot_x"),
                "is_target": locus.get("is_target"),
                "neighbour_symbol": "" if locus.get("is_target") else locus.get("symbol") or "",
                "side": locus.get("side") or "",
                "order": locus.get("rank"),
                "orientation": locus.get("strand") or "",
                "gene_id": locus.get("gene_id") or "",
                "placeholder": locus.get("placeholder"),
                "classification": locus.get("orthology_class") or "",
                "status": locus.get("identity_status") or "",
                "orthology_confidence": locus.get("mapping_confidence") or "",
                "distance_to_target": locus.get("distance"),
            })
    return rows


def _analysis_availability(inventory: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows = []
    for r in inventory:
        rows.append({
            "species_id": r["species_id"],
            "scientific_name": r["scientific_name"],
            "exon_structure": "available" if r.get("n_exons") else "unavailable",
            "domain_architecture": r.get("domain_status") or "pending",
            "exon_domain_boundaries": r.get("boundary_status") or "pending",
            "isoform_models": "available" if r.get("n_transcript_models") else "unavailable",
            "primary_protein": "available" if r.get("protein_id") else "unavailable",
        })
    return rows


_COMPARABLE_BOUNDARY_COLUMNS = [
    "comparable_boundary_group_id", "mapping_method", "mapping_status",
    "group_n_species", "group_species_coverage", "group_msa_column",
    "run_id", "species_id", "scientific_name", "protein_id", "transcript_id",
    "boundary_id", "left_exon_id", "right_exon_id", "exon_transition",
    "native_protein_position", "msa_column", "msa_mapping_status",
    "nearest_domain_instance_id", "nearest_domain_accession", "nearest_domain_label",
    "nearest_edge", "nearest_edge_position", "signed_distance", "absolute_distance",
    "boundary_class", "domain_annotation_available", "mapping_confidence",
    "review_reason",
]


def _comparable_boundary_rows(groups: Sequence[Dict[str, Any]],
                              run_id: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for g in groups:
        observations = g.get("per_species_native_positions") or []
        for o in observations:
            rows.append({
                "comparable_boundary_group_id": g.get("comparable_boundary_group_id"),
                "mapping_method": g.get("mapping_method"),
                "mapping_status": g.get("mapping_status"),
                "group_n_species": g.get("n_species") or len(observations),
                "group_species_coverage": g.get("species_coverage") or "",
                "group_msa_column": g.get("msa_column"),
                "run_id": run_id,
                "species_id": o.get("species_id"),
                "scientific_name": o.get("scientific_name"),
                "protein_id": o.get("protein_id"),
                "transcript_id": o.get("transcript_id"),
                "boundary_id": o.get("boundary_id"),
                "left_exon_id": o.get("left_exon_id"),
                "right_exon_id": o.get("right_exon_id"),
                "exon_transition": o.get("exon_transition"),
                "native_protein_position": o.get("native_position"),
                "msa_column": o.get("msa_column"),
                "msa_mapping_status": o.get("msa_mapping_status"),
                "nearest_domain_instance_id": o.get("nearest_domain_instance_id"),
                "nearest_domain_accession": o.get("nearest_domain_accession"),
                "nearest_domain_label": o.get("nearest_domain_label"),
                "nearest_edge": o.get("nearest_edge"),
                "nearest_edge_position": o.get("nearest_edge_position"),
                "signed_distance": o.get("signed_distance"),
                "absolute_distance": o.get("absolute_distance"),
                "boundary_class": o.get("boundary_class"),
                "domain_annotation_available": o.get("domain_annotation_available"),
                "mapping_confidence": o.get("mapping_confidence"),
                "review_reason": o.get("review_reason") or "",
            })
    return rows


def build_comparative_dataset(run_dir: Path,
                              coordinate_index: Optional[Dict[str, Any]] = None,
                              project_root: Optional[Path] = None) -> Dict[str, Any]:
    run_dir = Path(run_dir)
    root = Path(project_root or ROOT)
    pcm_path = run_dir / "website_indices" / "generic" / "protein_coordinate_model.json"
    index = coordinate_index or _load_json(pcm_path) or {}
    models = list(index.get("models") or [])
    gene = index.get("gene_symbol") or "GENE"
    run_id = index.get("run_id") or run_dir.name

    out_dir = run_dir / "results" / "generic_gene_analysis" / "comparative"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Rebuild the full native→column tables from the alignment file. The published
    # coordinate-model strips the pair list to keep the JSON small; figures and the
    # package builder still need the geometry.
    alignment_rel = (index.get("msa_coordinate_map") or {}).get("alignment_file") or ""
    alignment_path = root / alignment_rel if alignment_rel else (
        run_dir / "results/generic_gene_analysis/msa/primaries_msa.aln.faa")
    if not alignment_path.is_file():
        alignment_path = run_dir / "results/generic_gene_analysis/msa/primaries_msa.aln.faa"
    coord_map = build_msa_coordinate_map(alignment_path)
    tables = lookup_tables(coord_map) if coord_map.get("available") else {}

    # One canonical species order for every comparative view of this dataset, so
    # a reader scanning down a column meets the same species in the same places
    # in the matrix, the architecture figure and the exported tables.
    order_doc = species_order.build_species_order(
        [m.get("species_id") for m in models if m.get("species_id")])
    rank = {r["species_id"]: r["display_order"] for r in order_doc["species"]}
    models.sort(key=lambda m: rank.get(m.get("species_id"), len(rank)))

    completion = _load_json(run_dir / "status.json") or {}
    inventory = _species_inventory(models, completion.get("species_completion"))
    availability = _analysis_availability(inventory)
    msa_exons = _msa_aligned_exons(models, tables)
    msa_domains = _msa_aligned_domains(models, tables)
    domain_groups, domain_matrix = _domain_annotation_matrix(models, msa_domains)
    identity = _pairwise_identity(alignment_path)
    isoforms = _isoform_diversity(models, run_dir)
    synteny_blocks = _synteny_blocks(run_dir, models)
    synteny = _synteny_rows(synteny_blocks)

    multi = ((index.get("boundary_dashboard") or {}).get("multi_species") or {})
    comparable_boundaries = multi.get("comparable_boundary_groups") or []
    boundary_matrix = multi.get("boundary_matrix") or []
    distance_stats = multi.get("distance_statistics") or []
    inspection_cases = multi.get("inspection_cases") or []

    # One row per species-specific Boundary observation, built from the canonical
    # per-species Boundary index rather than from whichever member list a
    # comparable group happens to expose.
    boundary_long = boundary_observations.build_rows(run_dir)

    consistency_rows = []
    for s in distance_stats:
        consistency_rows.append({
            "comparable_boundary_group_id": s.get("comparable_boundary_group_id") or s.get("id"),
            "n_species_observed": s.get("n_species_observed") or s.get("species_coverage"),
            "mapping_coverage": s.get("mapping_coverage"),
            "exact_near_proportion": s.get("exact_near_proportion"),
            "raw_signed_distances": ",".join(
                str(x) for x in (s.get("raw_signed_distances") or s.get("signed_distances") or [])),
            "cross_species_difference": s.get("cross_species_difference")
            or s.get("pairwise_difference"),
            "distance_range": s.get("distance_range"),
            "dominant_class": s.get("dominant_class"),
            "domain_annotation_availability": s.get("domain_annotation_availability"),
            "mapping_confidence_distribution": json.dumps(
                s.get("mapping_confidence_distribution") or {}),
        })

    table_files = {
        "species_inventory": ("species_inventory.tsv", inventory,
                              list(inventory[0].keys()) if inventory else None),
        "analysis_availability": ("analysis_availability.tsv", availability,
                                  list(availability[0].keys()) if availability else None),
        "msa_aligned_exons": ("msa_aligned_exons.tsv", msa_exons,
                             list(msa_exons[0].keys()) if msa_exons else None),
        "msa_aligned_domains": ("msa_aligned_domains.tsv", msa_domains,
                               list(msa_domains[0].keys()) if msa_domains else None),
        "domain_annotation_matrix": ("domain_annotation_matrix.tsv", domain_matrix,
                                    list(domain_matrix[0].keys()) if domain_matrix else None),
        "comparable_domain_groups": (
            "comparable_domain_groups.tsv",
            [{k: g.get(k) for k in (
                "comparable_domain_group_id", "interpro_accession", "instance_number",
                "label", "n_species", "mapping_status")} for g in domain_groups],
            None),
        "species_order": ("species_order.tsv", order_doc["species"],
                          list(species_order.TSV_COLUMNS)),
        "pairwise_identity": ("pairwise_identity.tsv", identity,
                              list(identity[0].keys()) if identity else None),
        "isoform_diversity": ("isoform_diversity.tsv", isoforms,
                              list(isoforms[0].keys()) if isoforms else None),
        "comparative_synteny": ("comparative_synteny.tsv", synteny,
                                list(synteny[0].keys()) if synteny else None),
        # One row per species observation, not one row per group. A group summary alone
        # states that species agree on a boundary without saying which species, at which
        # residue, against which domain instance — so a reader cannot check the claim. The
        # group columns repeat on each of its rows, which keeps both readings available.
        "comparable_boundary_groups": (
            "comparable_boundary_groups.tsv",
            _comparable_boundary_rows(comparable_boundaries, run_id),
            None),
        "exon_domain_boundaries_long": ("exon_domain_boundaries_long.tsv", boundary_long,
                                        list(boundary_observations.COLUMNS)),
        "boundary_consistency_summary": ("boundary_consistency_summary.tsv", consistency_rows,
                                         list(consistency_rows[0].keys()) if consistency_rows else None),
        "inspection_cases": (
            "inspection_cases.tsv",
            [{
                "case_id": c.get("case_id") or c.get("id"),
                "case_type": c.get("case_type") or c.get("type"),
                "comparable_boundary_group_id": c.get("comparable_boundary_group_id"),
                "summary": c.get("summary") or c.get("title") or c.get("label"),
            } for c in inspection_cases],
            None),
    }

    artefacts: Dict[str, str] = {}
    for key, (fname, rows, fields) in table_files.items():
        path = out_dir / fname
        _write_tsv(path, rows, fields)
        artefacts[key] = _rel(path, root)

    # Primary proteins FASTA / MSA paths (referenced, not duplicated here).
    primary_faa = run_dir / "results/generic_gene_analysis/msa/primaries_msa_input.faa"
    primary_aln = alignment_path if alignment_path.is_file() else None

    comparative = {
        "schema_version": "1.0",
        "dataset_id": run_id,
        "gene_symbol": gene,
        "n_species": len(models),
        "available": len(models) >= 2,
        "reason": ("" if len(models) >= 2
                   else "comparative layer requires at least two species"),
        "msa": {
            "available": bool(coord_map.get("available")),
            "reason": coord_map.get("reason") or "",
            "n_columns": coord_map.get("n_columns") or 0,
            "alignment_file": _rel(alignment_path, root) if primary_aln else "",
            "primary_proteins_file": _rel(primary_faa, root) if primary_faa.is_file() else "",
        },
        "species_order": order_doc,
        "species_inventory": inventory,
        "analysis_availability": availability,
        "msa_aligned_exons": msa_exons,
        "msa_aligned_domains": msa_domains,
        "comparable_domain_groups": [
            {k: g[k] for k in g if k != "members"} for g in domain_groups
        ],
        "domain_annotation_matrix": domain_matrix,
        "pairwise_identity": identity,
        "isoform_diversity": isoforms,
        "synteny": synteny,
        "synteny_neighbourhood": synteny_blocks,
        "comparable_boundary_groups": comparable_boundaries,
        "boundary_matrix": boundary_matrix,
        "boundary_consistency": distance_stats,
        "inspection_cases": inspection_cases,
        "boundary_long": boundary_long,
        "artefacts": artefacts,
        "source_coordinate_model": _rel(pcm_path, root) if pcm_path.is_file() else "",
    }

    # Publish under website_indices/generic so the API and package builder share one
    # path, and keep a copy next to the TSV tables for the ZIP layout.
    for dest in (
        run_dir / "website_indices" / "generic" / "comparative_dataset_index.json",
        out_dir / "comparative_dataset_index.json",
    ):
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(comparative, indent=2), encoding="utf-8")

    return comparative


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse
    p = argparse.ArgumentParser(description='Shared comparative dataset layer for multi-species runs.')
    p.add_argument("--run-id", required=True)
    p.add_argument("--runs-root", type=Path, default=None)
    args = p.parse_args(argv)
    if args.runs_root is not None:
        run_dir = Path(args.runs_root).expanduser().resolve() / args.run_id
    else:
        try:
            record = resolve_run_record(RUNTIME_CONFIG, args.run_id)
        except RegistryError as exc:
            raise SystemExit(str(exc)) from exc
        if record is None:
            raise SystemExit(f"run not found: {args.run_id}")
        if record.read_only:
            raise SystemExit("run is registered read-only; copy it before rebuilding")
        run_dir = record.path
    if not run_dir.is_dir():
        raise SystemExit(f"run not found: {run_dir}")
    doc = build_comparative_dataset(run_dir)
    print(f"comparative dataset: available={doc['available']} "
          f"n_species={doc['n_species']} "
          f"msa_columns={doc['msa']['n_columns']} "
          f"domain_groups={len(doc['comparable_domain_groups'])} "
          f"boundary_groups={len(doc['comparable_boundary_groups'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
