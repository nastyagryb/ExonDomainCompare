from __future__ import annotations

"""Read-only adapters for the versioned, species-centred dataset API model."""

import json
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

_REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(_REPOSITORY_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(_REPOSITORY_ROOT / "scripts"))
from exondomaincompare.contracts import stamp_payload  # noqa: E402


CANONICAL_DATASET_MODEL_VERSION = "1.0.0"

_FGFR2_LEGACY_FILES = {
    "overview": "run_index.json",
    "species": "species_index.json",
    "evidence": "evidence_stack.json",
    "protein_architecture": "coordinate_track_index.json",
    "synteny": "synteny_locus_index.json",
    "event_evidence": "cassette_residue_index.json",
    "domain_architecture": "domain_architecture_index.json",
    "boundary": "boundary_consistency_index.json",
    "figures": "figure_index.json",
    "downloads": "download_index.json",
}

_FGFR2_VALIDATED_EVENT_FILES = {
    "cassette_residue_index": "cassette_residue_index.json",
    "coordinate_track_index": "coordinate_track_index.json",
    "msa_index": "msa_index.json",
    "species_story_index": "species_story_index.json",
    "boundary_consistency_index": "boundary_consistency_index.json",
    "boundary_consistency_summary": "boundary_consistency_summary.json",
    "boundary_consistency_matrix": "boundary_consistency_matrix.json",
    "boundary_consistency_outliers": "boundary_consistency_outliers.json",
}

_SHARED_FILES = {
    "overview": "overview_index.json",
    "evidence": "evidence_stack.json",
    "gene_explorer": "gene_explorer_index.json",
    "transcript_exon_structure": "transcript_exon_structure_index.json",
    "primary_selection_evidence": "primary_selection_index.json",
    "isoform_alignment": "isoform_alignment_index.json",
    "candidate_evidence": "event_candidate_evidence_index.json",
    "protein_architecture": "protein_architecture_index.json",
    "synteny": "synteny_index.json",
    "event_evidence": "event_evidence_index.json",
    "domain_architecture": "domain_architecture_index.json",
    "boundary": "exon_domain_boundaries_index.json",
    "figures": "figures_index.json",
    "downloads": "download_index.json",
    "available_views": "available_views.json",
    "dataset_summary": "dataset_summary.json",
    # FGFR2-compatible shapes for shared Gene Explorer modules (CoordinateTrack, MsaExplorer, SyntenyViewer).
    "coordinate_track_index": "coordinate_track_index.json",
    "msa_index": "msa_index.json",
    "synteny_locus_index": "synteny_locus_index.json",
    "protein_coordinate_model": "protein_coordinate_model.json",
}


def _read_json(path: Path, fallback: Any) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return fallback


def _index_reader(descriptor: Mapping[str, Any]):
    indices_dir = Path(descriptor["indices_dir"])
    derived = descriptor.get("derived_indices_dir")
    derived_dir = Path(derived) if derived else None

    def read(filename: str, fallback: Any) -> Any:
        if derived_dir is not None:
            found = _read_json(derived_dir / filename, None)
            if found is not None:
                return found
        return _read_json(indices_dir / filename, fallback)

    return read


def _read_shared(indices_dir: Path, filename: str, fallback: Any) -> Any:
    for path in (indices_dir / filename, indices_dir / "generic" / filename):
        if path.is_file():
            return _read_json(path, fallback)
    return fallback


# clade (species_registry.tsv) -> the same taxon-group vocabulary the FGFR2 UI uses.
_CLADE_TO_TAXON_GROUP = {
    "mammal": "Other mammals",
    "primate": "Primates",
    "bird": "Birds",
    "reptile": "Reptiles",
    "amphibian": "Amphibians",
    "fish": "Teleost fish",
}


def _read_tsv_rows(path: Path) -> List[Dict[str, str]]:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return []
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    rows: List[Dict[str, str]] = []
    for line in lines[1:]:
        parts = line.split("\t")
        rows.append({header[i]: (parts[i] if i < len(parts) else "") for i in range(len(header))})
    return rows


def _species_registry_meta(run_base: Optional[Path]) -> Dict[str, Dict[str, Any]]:
    if not run_base:
        return {}
    reg = Path(run_base) / "results" / "01_species_registry" / "species_registry.tsv"
    meta: Dict[str, Dict[str, Any]] = {}
    for row in _read_tsv_rows(reg):
        sid = row.get("species_id") or ""
        if not sid:
            continue
        clade = (row.get("clade") or "").strip().lower()
        scientific = row.get("scientific_name") or sid.replace("_", " ").capitalize()
        common = (row.get("common_name") or "").strip()
        meta[sid] = {
            "species_id": sid,
            "scientific_name": scientific,
            "common_name": common.capitalize() if common else "",
            "clade": clade,
            "taxonomic_group": _CLADE_TO_TAXON_GROUP.get(clade, "Analysed species"),
            "analysis_status": "analysed",
        }
    return meta


def _curation_status(protein: Mapping[str, Any]) -> str:
    curated = str(protein.get("curated") or "").strip().lower()
    if curated in {"yes", "true", "curated"}:
        return "curated"
    source = str(protein.get("source") or "").lower()
    if "curated" in source:
        return "curated"
    if "predicted" in source:
        return "predicted"
    return "unknown"


def _normalize_protein_models(
    primary_selection: Mapping[str, Any], species_id: str
) -> List[Dict[str, Any]]:
    proteins = primary_selection.get("proteins") if isinstance(primary_selection, dict) else None
    if not isinstance(proteins, list):
        return []
    # Per-species primary (multi-species aware). ``species_primaries`` maps each
    # species to its own primary; fall back to the global primary for single-species.
    sp_primaries = {
        str(sp.get("species_id")): str(sp.get("primary_protein_id") or "")
        for sp in (primary_selection.get("species_primaries") or [])
        if isinstance(sp, dict)
    }
    primary_id = sp_primaries.get(species_id) or str(primary_selection.get("primary_protein_id") or "")
    rule = str(primary_selection.get("selection_rule") or "")
    rule_label = str(primary_selection.get("selection_rule_label") or "")
    # If proteins carry species_id, restrict to the requested species so models
    # from other species never leak into a species' model list.
    has_species = any(isinstance(p, dict) and p.get("species_id") for p in proteins)
    models: List[Dict[str, Any]] = []
    for protein in proteins:
        if not isinstance(protein, dict):
            continue
        if has_species and species_id and str(protein.get("species_id") or "") != species_id:
            continue
        pid = str(protein.get("protein_id") or "")
        status = str(protein.get("primary_status") or "").lower()
        is_primary = (bool(primary_id) and pid == primary_id) or (not primary_id and status == "primary")
        models.append({
            "species_id": species_id,
            "protein_id": pid,
            "transcript_id": str(protein.get("transcript_id") or ""),
            "length_aa": protein.get("length_aa"),
            "protein_length": protein.get("length_aa"),
            "is_primary": is_primary,
            "primary_status": "primary" if is_primary else "alternative",
            "role": "primary" if is_primary else "alternative",
            "curation_status": _curation_status(protein),
            "source_label": str(protein.get("source_label") or protein.get("source") or ""),
            "selection_rule": rule,
            "selection_reason": (
                str(primary_selection.get("explanation") or rule_label)
                if is_primary else str(protein.get("source_label") or "Alternative isoform")
            ),
        })
    return models


def _dir_has_files(folder: Path, suffix: Optional[str] = None) -> bool:
    if not folder.is_dir():
        return False
    for p in folder.rglob("*"):
        if p.is_file() and (suffix is None or p.name.endswith(suffix)):
            return True
    return False


def _has_content(path: Path) -> bool:
    try:
        if path.stat().st_size == 0:
            return False
        suffix = path.suffix.lower()
        if suffix in {".tsv", ".csv", ".txt"}:
            with path.open("r", encoding="utf-8", errors="replace") as fh:
                rows = 0
                for line in fh:
                    if line.strip():
                        rows += 1
                    if rows > 1:
                        return True
            return False
        if suffix == ".json":
            payload = json.loads(path.read_text(encoding="utf-8"))
            return bool(payload)
    except (OSError, ValueError):
        return False
    return True




def _post_cluster_complete(base: Path) -> bool:
    model = base / "website_indices" / "generic" / "protein_coordinate_model.json"
    try:
        payload = json.loads(model.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return any(m.get("status") == "available" for m in payload.get("models", []))


# Provenance catalogue for shared exploratory runs. Each entry documents how an
# artefact was generated so the Files tab becomes a provenance hub. Columns:
# (category, label, path-relative-to-run, generated_by, stage, source_inputs)
_PROVENANCE_CATALOGUE = [
    ("Input and run configuration", "Run configuration (JSON)", "run_config.json",
     "run_core_gene_analysis.py", "run_setup", []),
    ("Input and run configuration", "Species registry (TSV)",
     "results/01_species_registry/species_registry.tsv", "build_species_registry", "run_setup", ["NCBI taxonomy"]),
    ("Gene models", "Gene model index (TSV)", "results/core_gene_analysis/gene_model_index.tsv",
     "parse_gene_models", "models", ["genomic.gff"]),
    ("Protein sequences", "All protein isoforms (FASTA)", "results/core_gene_analysis/proteins_all_isoforms.faa",
     "build_core_contract", "models", ["protein.faa"]),
    ("Protein sequences", "Primary protein (FASTA)", "results/core_gene_analysis/protein_primary.faa",
     "select_primary_protein", "models", ["proteins_all_isoforms.faa"]),
    ("Primary selection", "Primary selection report (JSON)", "results/core_gene_analysis/primary_selection_report.json",
     "select_primary_protein", "models", ["gene_model_index.tsv"]),
    ("Primary selection", "Primary selection evidence (TSV)", "results/core_gene_analysis/primary_selection_evidence.tsv",
     "select_primary_protein", "models", ["gene_model_index.tsv"]),
    ("Exon mapping", "Exon → protein map (TSV)", "results/core_gene_analysis/exon_protein_map.tsv",
     "map_exons_to_protein", "coordinate_mapping", ["genomic.gff", "protein_primary.faa"]),
    ("Coordinate evidence", "Coordinate evidence register (TSV)",
     "results/core_gene_analysis/evidence_register/coordinate_evidence_register.tsv",
     "build_coordinate_evidence_register", "coordinate_audit",
     ["gene_model_index.tsv", "exon_protein_map.tsv", "domain_features.tsv",
      "exon_domain_boundary_distances.tsv", "post_cluster_qc.json"]),
    ("Coordinate evidence", "Coordinate evidence register with provenance (JSON)",
     "results/core_gene_analysis/evidence_register/coordinate_evidence_register.json",
     "build_coordinate_evidence_register", "coordinate_audit",
     ["run_config.json", "core_model_collection_report.json", "core_gene_report.json",
      "post_cluster_qc.json"]),
    ("Isoform alignment", "Protein isoform alignment (FASTA)", "results/07_msa/protein_alignment.faa",
     "MAFFT (--auto)", "msa", ["proteins_all_isoforms.faa"]),
    ("Local neighbourhood", "Local neighbourhood (TSV)", "results/core_gene_analysis/synteny_neighbors.tsv",
     "extract_local_neighbourhood", "synteny", ["genomic.gff"]),
    ("Candidate evidence", "Exploratory candidate evidence (TSV)", "results/core_gene_analysis/event_candidate_evidence.tsv",
     "build_exploratory_evidence", "event_evidence", ["protein_alignment.faa", "exon_protein_map.tsv"]),
    ("Candidate evidence", "Candidate ranking (TSV)", "results/generic_gene_analysis/event_candidate_ranking.tsv",
     "rank_candidates", "event_evidence", ["event_candidate_evidence.tsv"]),
    ("Cluster input", "Cluster InterProScan input (FASTA)",
     "results/13_final_pre_interpro_closure/freeze/final_pre_interpro_proteins_primary.faa",
     "prepare_cluster_input", "cluster_input", ["protein_primary.faa"]),
    # Real post-cluster scientific outputs. These are listed individually so a
    # results-ready run shows the actual tables instead of a pending placeholder.
    ("Post-cluster analysis outputs", "Normalized InterPro domain table (TSV)",
     "results/core_gene_analysis/interpro_annotations.tsv",
     "normalize_interproscan_output", "post_cluster", ["input.fasta.tsv"]),
    ("Post-cluster analysis outputs", "Representative domain table (TSV)",
     "results/core_gene_analysis/domain_features.tsv",
     "select_representative_domains", "post_cluster", ["interpro_annotations.tsv"]),
    ("Post-cluster analysis outputs", "Domain architecture table (TSV)",
     "results/15_domain_architecture/domain_architecture.tsv",
     "build_domain_architecture", "post_cluster", ["domain_features.tsv"]),
    ("Post-cluster analysis outputs", "Normalized pyTMHMM topology table (TSV)",
     "results/15_domain_architecture/tm_features.tsv",
     "normalize_pytmhmm_output", "post_cluster", ["pytmhmm_transmembrane_hits.tsv"]),
    ("Post-cluster analysis outputs", "Candidate domain-context table (TSV)",
     "results/generic_gene_analysis/candidate_domain_context.tsv",
     "build_candidate_domain_context", "post_cluster",
     ["domain_features.tsv", "event_candidate_ranking.tsv"]),
    ("Post-cluster analysis outputs", "Exon–domain boundary table (TSV)",
     "results/16_final_analyses/exon_domain_boundary_analysis.tsv",
     "boundary_classification.py", "post_cluster", ["domain_features.tsv", "exon_protein_map.tsv"]),
    ("Post-cluster analysis outputs", "Exon–domain boundary distances (TSV)",
     "results/16_final_analyses/exon_domain_boundary_distances.tsv",
     "boundary_classification.py", "post_cluster", ["domain_features.tsv", "exon_protein_map.tsv"]),
    ("Post-cluster analysis outputs", "Exon–domain boundary summary (TSV)",
     "results/16_final_analyses/exon_domain_boundary_summary.tsv",
     "boundary_classification.py", "post_cluster", ["exon_domain_boundary_distances.tsv"]),
    ("Post-cluster analysis outputs", "Post-cluster QC report (JSON)",
     "results/15_domain_architecture/post_cluster_qc.json",
     "build_domain_architecture", "post_cluster", ["domain_architecture.tsv"]),
]

# Directories that can hold the real post-cluster analysis products. The pending
# row for "Post-cluster analysis" is only emitted when none of them has content.
_POST_CLUSTER_DIRS = (
    ("results", "15_domain_architecture"),
    ("results", "16_final_analyses"),
    ("results", "16_final_thesis_analyses"),
)

_PREVIEWABLE = {"tsv", "json", "txt", "csv", "faa", "fasta", "md"}
_PUBLIC_DOWNLOAD_PROJECTIONS = {
    "run_config.json": "website_indices/public/run_config.json",
    "results/15_domain_architecture/post_cluster_qc.json":
        "website_indices/public/post_cluster_qc.json",
}


def _shared_downloads(
    indices_dir: Path, run_base: Optional[Path], run_id: str
) -> Dict[str, Any]:
    raw = _read_shared(indices_dir, "download_index.json", None)
    if isinstance(raw, list):
        return {"available": bool(raw), "items": raw}
    if isinstance(raw, dict) and isinstance(raw.get("items"), list):
        return {"available": bool(raw.get("items")), "items": raw["items"]}
    if not run_base:
        return {"available": False, "items": [], "reason": "download_index_not_provided"}

    base = Path(run_base)
    post_complete = _post_cluster_complete(base)
    items: List[Dict[str, Any]] = []
    for group, label, rel, generated_by, stage, sources in _PROVENANCE_CATALOGUE:
        public_rel = _PUBLIC_DOWNLOAD_PROJECTIONS.get(rel, rel)
        fp = base / public_rel
        if not fp.is_file():
            continue
        # A post-cluster placeholder written before the cluster round-trip must
        # not be offered as a finished result. Once the coordinate model has been
        # rebuilt from real results, an empty table is a genuine zero finding.
        if stage == "post_cluster" and not post_complete and not _has_content(fp):
            continue
        try:
            size = fp.stat().st_size
        except OSError:
            continue
        fmt = fp.suffix.lstrip(".").lower()
        items.append({
            "group": group,
            "category": group,
            "label": label,
            "name": fp.name,
            "path": f"runs/{run_id}/{public_rel}",
            "format": fmt,
            "size_bytes": size,
            "size_human": _human_size(size),
            "generated_by": generated_by,
            "stage": stage,
            "source_inputs": sources,
            "status": "available",
            "previewable": fmt in _PREVIEWABLE and size <= 256 * 1024,
        })

    # Generated figures deliberately do NOT appear here: the Files page carries
    # scientific data files only, every figure lives in the Figure Gallery.

    # Cluster stages: shown as pending only while the real outputs are absent.
    ips_out = base / "results" / "14_interproscan" / "primary" / "output"
    tm_out = (base / "results" / "15_exon_domain_boundary_post_interpro"
              / "pytmhmm_primary" / "output")
    for cat, folder, gen in [
        ("InterProScan", ips_out, "InterProScan (cluster)"),
        ("pyTMHMM", tm_out, "pyTMHMM (cluster)"),
    ]:
        output_files = sorted(path for path in folder.rglob("*") if path.is_file())
        if not output_files:
            items.append({
                "group": cat, "category": cat,
                "label": f"{cat} output",
                "name": folder.name,
                "path": None, "format": "",
                "generated_by": gen, "stage": "post_cluster", "source_inputs": [],
                "status": "pending", "previewable": False,
            })
            continue
        for output_file in output_files:
            size = output_file.stat().st_size
            fmt = output_file.suffix.lstrip(".").lower()
            items.append({
                "group": cat, "category": cat,
                "label": f"{cat} output — {output_file.name}",
                "name": output_file.name,
                "path": f"runs/{run_id}/{output_file.relative_to(base).as_posix()}",
                "format": fmt, "size_bytes": size, "size_human": _human_size(size),
                "generated_by": gen, "stage": "post_cluster", "source_inputs": [],
                "status": "available",
                "previewable": fmt in _PREVIEWABLE and size <= 256 * 1024,
            })

    # Only non-empty artefacts count as produced; pending entries expose no path.
    produced = {d["label"] for d in items if d["stage"] == "post_cluster"}
    for group, label, rel, generated_by, stage, sources in _PROVENANCE_CATALOGUE:
        if stage != "post_cluster" or label in produced:
            continue
        items.append({
            "group": group, "category": group, "label": label,
            "name": Path(rel).name,
            "path": None, "format": Path(rel).suffix.lstrip(".").lower(),
            "generated_by": generated_by, "stage": "post_cluster",
            "source_inputs": sources, "status": "pending", "previewable": False,
        })

    return {"available": bool(items), "items": items}


def _human_size(num: int) -> str:
    value = float(num)
    for unit in ("B", "KB", "MB", "GB"):
        if value < 1024 or unit == "GB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    return f"{num} B"


def _analysis(
    analysis_id: str,
    gene_symbol: str,
    *,
    support_level: str,
    experimental: bool,
    has_event: bool,
    event_id: str = "",
    event_type: str = "",
    display_name: str = "",
) -> Dict[str, Any]:
    return {
        "id": analysis_id,
        "display_name": display_name or analysis_id,
        "gene_symbol": gene_symbol,
        "support_level": support_level,
        "experimental": experimental,
        "has_event": has_event,
        "event_id": event_id,
        "event_type": event_type,
    }


def _tabs(available_views: Mapping[str, Any], event_layer_type: str) -> List[Dict[str, Any]]:
    labels = {
        "overview": "Overview",
        "gene_explorer": "Gene Explorer",
        "gene_models": "Gene Models",
        "protein_architecture": "Protein Architecture",
        "event_region": "Event Region",
        "event_evidence": "Event Evidence",
        "domain_architecture": "Domain Architecture",
        "exon_domain_boundaries": "Exon–Domain Boundaries",
        "synteny": "Synteny",
        "boundary_consistency": "Boundary Consistency",
        "boundary_relation": "Boundary Relation",
        "figure_gallery": "Figures",
        "downloads": "Downloads",
    }
    result = []
    for tab_id, available in available_views.items():
        result.append({
            "id": tab_id,
            "label": labels.get(tab_id, tab_id.replace("_", " ").title()),
            "available": bool(available),
            "event_layer": event_layer_type if tab_id in {
                "event_region", "event_evidence", "boundary_consistency", "boundary_relation"
            } else "none",
        })
    return result


def _dataset_metadata(
    descriptor: Mapping[str, Any], *, adapter: str, source_layout: str
) -> Dict[str, Any]:
    kind = str(descriptor.get("kind") or "run")
    run_id = str(descriptor.get("run_id") or ("example" if kind == "example" else ""))
    dataset_id = "example" if kind == "example" else f"run:{run_id}"
    return {
        "id": dataset_id,
        "kind": kind,
        "run_id": run_id,
        "read_only": bool(descriptor.get("read_only", False)),
        "adapter": adapter,
        "adapter_read_only": True,
        "source_layout": source_layout,
    }


def adapt_fgfr2_legacy(descriptor: Mapping[str, Any]) -> Dict[str, Any]:
    _indices_dir = Path(descriptor["indices_dir"])
    run_base = descriptor.get("run_base")
    config = _read_json(Path(run_base) / "run_config.json", {}) if run_base else {}
    read_index = _index_reader(descriptor)

    blobs = {
        field: read_index(filename, [] if field == "species" else {})
        for field, filename in _FGFR2_LEGACY_FILES.items()
    }
    validated = {
        field: read_index(filename, {})
        for field, filename in _FGFR2_VALIDATED_EVENT_FILES.items()
    }
    overview = blobs["overview"] if isinstance(blobs["overview"], dict) else {}
    available_views = {
        "overview": bool(blobs["overview"]),
        "gene_explorer": bool(blobs["species"]),
        "protein_architecture": bool(blobs["protein_architecture"]),
        "event_region": bool(validated["cassette_residue_index"]),
        "event_evidence": bool(blobs["event_evidence"]),
        "domain_architecture": bool(blobs["domain_architecture"]),
        "exon_domain_boundaries": bool(blobs["boundary"]),
        "synteny": bool(blobs["synteny"]),
        "boundary_consistency": bool(blobs["boundary"]),
        "figure_gallery": bool(blobs["figures"]),
        "downloads": bool(blobs["downloads"]),
    }
    analysis_id = str(config.get("analysis_id") or "FGFR2_IIIb_IIIc")
    gene_symbol = str(config.get("gene_symbol") or "FGFR2")
    event_id = str(config.get("event_id") or "FGFR2_IIIb_IIIc_cassette")
    event_type = str(config.get("event_type") or "mutually_exclusive_cassette")

    return {
        "schema_version": CANONICAL_DATASET_MODEL_VERSION,
        "model_type": "CanonicalDatasetModel",
        "dataset": _dataset_metadata(
            descriptor, adapter="fgfr2_legacy", source_layout="legacy_website_indices"
        ),
        "analysis": _analysis(
            analysis_id,
            gene_symbol,
            support_level=str(config.get("support_level") or "validated_event_analysis"),
            experimental=bool(config.get("experimental", False)),
            has_event=True,
            event_id=event_id,
            event_type=event_type,
            display_name=str(overview.get("case_study") or "FGFR2 IIIb/IIIc"),
        ),
        "event_layer": {
            "type": "validated",
            "configured": True,
            "event_id": event_id,
            "event_type": event_type,
        },
        "available_views": available_views,
        "tabs": _tabs(available_views, "validated"),
        **blobs,
        # Exact legacy JSON payloads are retained under stable, explicit names.
        "legacy_fgfr2_indices": {
            Path(filename).stem: read_index(filename, {})
            for filename in sorted({
                *_FGFR2_LEGACY_FILES.values(),
                *_FGFR2_VALIDATED_EVENT_FILES.values(),
                "freeze_index.json",
                "domain_architecture_summary.json",
                "species_domain_architecture.json",
                "domain_architecture_qc.json",
            })
        },
        "validated_event_indices": validated,
    }


def _species_key(row: Mapping[str, Any]) -> str:
    return str(row.get("species_id") or row.get("species") or "")


def _species_rows(index: Any) -> Iterable[Mapping[str, Any]]:
    if isinstance(index, dict) and isinstance(index.get("species"), list):
        return (row for row in index["species"] if isinstance(row, dict))
    return ()


def _shared_species(
    gene_explorer: Mapping[str, Any],
    protein_architecture: Any,
    synteny: Any,
    event_evidence: Any,
    domain_architecture: Any,
    boundary: Any,
    species_meta: Optional[Mapping[str, Mapping[str, Any]]] = None,
    primary_selection: Optional[Mapping[str, Any]] = None,
) -> List[Dict[str, Any]]:
    species_meta = species_meta or {}
    primary_selection = primary_selection if isinstance(primary_selection, dict) else {}
    selected_primary_protein = str(primary_selection.get("primary_protein_id") or "")
    selected_primary_transcript = str(primary_selection.get("primary_transcript_id") or "")
    # Per-species primary (multi-species aware). Falls back to the global primary.
    sp_primaries = {
        str(sp.get("species_id")): sp
        for sp in (primary_selection.get("species_primaries") or [])
        if isinstance(sp, dict)
    }
    species_ids = {
        str(value) for value in gene_explorer.get("species", []) if value
    }
    for index in (protein_architecture, synteny, domain_architecture):
        species_ids.update(filter(None, (_species_key(row) for row in _species_rows(index))))
    species_ids.update(
        str(row.get("species_id"))
        for row in gene_explorer.get("isoforms", [])
        if isinstance(row, dict) and row.get("species_id")
    )

    def one(index: Any, species_id: str) -> Dict[str, Any]:
        return next(
            (dict(row) for row in _species_rows(index) if _species_key(row) == species_id),
            {},
        )

    clusters = event_evidence.get("clusters", []) if isinstance(event_evidence, dict) else []
    boundary_proteins = boundary.get("proteins", []) if isinstance(boundary, dict) else []
    result = []
    for species_id in sorted(species_ids):
        isoforms = [
            row for row in gene_explorer.get("isoforms", [])
            if isinstance(row, dict) and row.get("species_id") == species_id
        ]
        species_clusters = [
            cluster for cluster in clusters
            if isinstance(cluster, dict) and (
                not cluster.get("raw_support_rows")
                or any(
                    isinstance(row, dict) and row.get("species_id") == species_id
                    for row in cluster.get("raw_support_rows", [])
                )
            )
        ]
        species_boundaries = [
            row for row in boundary_proteins
            if isinstance(row, dict) and row.get("species_id") == species_id
        ]
        protein_node = one(protein_architecture, species_id)
        if not protein_node and isinstance(protein_architecture, dict) \
                and isinstance(protein_architecture.get("proteins"), list):
            protein_node = {"species_id": species_id,
                            "proteins": protein_architecture["proteins"],
                            "domain_status": protein_architecture.get("domain_annotation_status", "")}
        synteny_node = one(synteny, species_id)
        if not synteny_node and isinstance(synteny, dict) \
                and isinstance(synteny.get("neighbours"), list):
            synteny_node = {**synteny, "species_id": species_id,
                            "neighbors": synteny["neighbours"]}
        meta = species_meta.get(species_id, {})
        scientific = str(meta.get("scientific_name")
                         or species_id.replace("_", " ").capitalize())
        # Per-species primary protein/transcript (falls back to global).
        sp_prim = sp_primaries.get(species_id, {})
        sp_primary_protein = str(sp_prim.get("primary_protein_id") or "") or selected_primary_protein
        sp_primary_transcript = str(sp_prim.get("primary_transcript_id") or "") or selected_primary_transcript
        # One normalized protein-model list (single source of truth for primary).
        protein_models = _normalize_protein_models(primary_selection, species_id)
        result.append({
            "species_id": species_id,
            "species": species_id,
            "display_species_name": scientific,
            "scientific_name": scientific,
            "common_name": meta.get("common_name", ""),
            "taxonomic_group": meta.get("taxonomic_group", "Analysed species"),
            "taxon_group": meta.get("taxonomic_group", "Analysed species"),
            "clade": meta.get("clade", ""),
            "analysis_status": meta.get("analysis_status", "analysed"),
            "overall_status": "accepted",
            "selected_primary_protein": sp_primary_protein,
            "selected_primary_transcript": sp_primary_transcript,
            "proteins": protein_models,
            "protein_models": protein_models,
            "gene_explorer": {"isoforms": isoforms},
            "protein_architecture": protein_node,
            "synteny": synteny_node,
            "event_evidence": {"clusters": species_clusters},
            "domain_architecture": one(domain_architecture, species_id),
            "boundary": {"proteins": species_boundaries},
        })
    return result


def adapt_shared_run(descriptor: Mapping[str, Any]) -> Dict[str, Any]:
    indices_dir = Path(descriptor["indices_dir"])
    run_base = Path(descriptor["run_base"])
    config = _read_json(run_base / "run_config.json", {})
    run_status = _read_json(run_base / "status.json", {})
    shared = {
        field: _read_shared(indices_dir, filename, {})
        for field, filename in _SHARED_FILES.items()
    }
    summary = shared["dataset_summary"] if isinstance(shared["dataset_summary"], dict) else {}
    explorer = shared["gene_explorer"] if isinstance(shared["gene_explorer"], dict) else {}
    views_blob = shared["available_views"] if isinstance(shared["available_views"], dict) else {}
    available_views = dict(
        views_blob.get("available_views")
        or summary.get("available_views")
        or shared["overview"].get("available_views", {})
    )
    event_evidence = shared["event_evidence"]
    has_exploratory = bool(
        isinstance(event_evidence, dict)
        and event_evidence.get("available")
        and event_evidence.get("evidence_status") == "exploratory"
    )
    configured_event = bool(summary.get("has_event") or config.get("has_event"))
    event_layer_type = "validated" if configured_event else ("exploratory" if has_exploratory else "none")
    gene_symbol = str(
        summary.get("gene_symbol")
        or config.get("gene_symbol")
        or explorer.get("gene", {}).get("symbol")
        or ""
    )
    analysis_id = str(summary.get("analysis_id") or config.get("analysis_id") or "")
    run_id = str(descriptor.get("run_id") or run_base.name)
    downloads = _shared_downloads(indices_dir, run_base, run_id)

    # Single source of truth for the primary protein/transcript selection.
    primary_selection = shared["primary_selection_evidence"] \
        if isinstance(shared["primary_selection_evidence"], dict) else {}
    species_meta = _species_registry_meta(run_base)
    selected_primary_protein = str(primary_selection.get("primary_protein_id") or "")
    selected_primary_transcript = str(primary_selection.get("primary_transcript_id") or "")
    default_species = next(
        (str(v) for v in explorer.get("species", []) if v), ""
    ) or str(primary_selection.get("species_id") or "")
    protein_models = _normalize_protein_models(primary_selection, default_species)

    return {
        "schema_version": CANONICAL_DATASET_MODEL_VERSION,
        "model_type": "CanonicalDatasetModel",
        "dataset": _dataset_metadata(
            descriptor, adapter="shared_run", source_layout="website_indices_root_or_generic"
        ),
        "analysis": _analysis(
            analysis_id,
            gene_symbol,
            support_level=str(summary.get("support_level") or config.get("support_level") or ""),
            experimental=bool(summary.get("experimental", config.get("experimental", False))),
            has_event=configured_event,
            event_id=str(summary.get("event_id") or config.get("event_id") or ""),
            event_type=str(summary.get("event_type") or config.get("event_type") or ""),
            display_name=str(explorer.get("analysis_display_name") or analysis_id),
        ),
        "event_layer": {
            "type": event_layer_type,
            "configured": configured_event,
            "evidence_status": (
                event_evidence.get("evidence_status", "none")
                if isinstance(event_evidence, dict) else "none"
            ),
        },
        "available_views": available_views,
        "tabs": _tabs(available_views, event_layer_type),
        "overview": shared["overview"],
        "analysis_stage": {
            "status": run_status.get("status", "models_ready"),
            "current_step": run_status.get("current_step", ""),
            "cluster_status": run_status.get("cluster_analysis_status", "not_started"),
            "post_interpro_status": run_status.get("post_interpro_status", "not_started"),
            "reason": run_status.get("failed_reason") or run_status.get("error") or "",
            "cluster_command": run_status.get("cluster_command", ""),
        },
        "species": _shared_species(
            explorer,
            shared["protein_architecture"],
            shared["synteny"],
            shared["event_evidence"],
            shared["domain_architecture"],
            shared["boundary"],
            species_meta=species_meta,
            primary_selection=primary_selection,
        ),
        "selected_primary_protein": selected_primary_protein,
        "selected_primary_transcript": selected_primary_transcript,
        "evidence": shared["evidence"],
        "transcript_exon_structure": shared["transcript_exon_structure"],
        "protein_models": protein_models or explorer.get("isoforms", []),
        "primary_selection_evidence": shared["primary_selection_evidence"],
        "isoform_alignment": shared["isoform_alignment"],
        "candidate_evidence": shared["candidate_evidence"],
        "local_synteny": shared["synteny"],
        "domain_features": shared["domain_architecture"],
        "tm_features": shared["domain_architecture"],
        "candidate_domain_context": (
            shared["candidate_evidence"].get("candidate_domain_context", {})
            if isinstance(shared["candidate_evidence"], dict) else {}
        ),
        "exon_domain_boundaries": shared["boundary"],
        "protein_architecture": shared["protein_architecture"],
        "synteny": shared["synteny"],
        "event_evidence": shared["event_evidence"],
        "domain_architecture": shared["domain_architecture"],
        "boundary": shared["boundary"],
        "figures": shared["figures"],
        "downloads": downloads,
        # Validated protein-coordinate model (single source of truth for the Exon Map).
        "protein_coordinate_model": shared.get("protein_coordinate_model") or {},
        # Convenience aliases used by the Comparative Figure Gallery scope selector.
        "models": (shared.get("protein_coordinate_model") or {}).get("models") or [],
        "coordinate_models": (shared.get("protein_coordinate_model") or {}).get("models") or [],
        "comparative_dataset": _read_shared(
            indices_dir, "comparative_dataset_index.json", {}) or {},
        # Exact root/generic payloads remain available during frontend migration.
        "shared_indices": shared,
        # Same field names as validated_event_indices so the frontend reuses FGFR2 tab components.
        "shared_event_indices": {
            "coordinates": shared.get("coordinate_track_index") or {},
            "msa": shared.get("msa_index") or {},
            "synteny": shared.get("synteny_locus_index") or {},
            "domainArch": shared.get("domain_architecture") or {},
            "boundaryMatrix": shared.get("boundary") or {},
        },
    }


def _stamp_identity(model: Dict[str, Any], descriptor: Mapping[str, Any]) -> Dict[str, Any]:
    dataset = model.get("dataset") or {}
    run_id = dataset.get("run_id") or descriptor.get("run_id") or ""
    model["dataset_id"] = dataset.get("id") or descriptor.get("id") or ""
    model["run_id"] = run_id
    run_base = descriptor.get("run_base")
    version = ""
    if run_base:
        try:
            from exondomaincompare.shared_gene_analysis.analysis_availability import index_version
            version = index_version(Path(run_base))
        except Exception:
            version = ""
    model["index_version"] = version
    return stamp_payload(
        model,
        payload_type="canonical_dataset",
        run_id=str(run_id),
        dataset_id=str(model["dataset_id"]),
        generator="webapp/backend/canonical_dataset.py",
    )


def build_canonical_dataset_model(descriptor: Mapping[str, Any]) -> Dict[str, Any]:
    if descriptor.get("kind") == "example":
        return _stamp_identity(adapt_fgfr2_legacy(descriptor), descriptor)

    run_base: Optional[Path] = (
        Path(descriptor["run_base"]) if descriptor.get("run_base") else None
    )
    config = _read_json(run_base / "run_config.json", {}) if run_base else {}
    gene_symbol = str(config.get("gene_symbol") or "").upper()
    has_legacy_species = (Path(descriptor["indices_dir"]) / "species_index.json").is_file()
    if gene_symbol == "FGFR2" or has_legacy_species:
        model = adapt_fgfr2_legacy(descriptor)
    else:
        model = adapt_shared_run(descriptor)
    model = _stamp_identity(model, descriptor)
    # Availability is derived here rather than in each index builder so runs already on
    # disk report the canonical states without being rebuilt.
    if run_base:
        try:
            from exondomaincompare.shared_gene_analysis.analysis_availability import annotate_dataset_model
            model = annotate_dataset_model(model, run_base)
        except Exception:
            pass
    return stamp_payload(
        model, payload_type="canonical_dataset",
        run_id=str(model.get("run_id") or ""),
        dataset_id=str(model.get("dataset_id") or ""),
        generator="webapp/backend/canonical_dataset.py",
    )
