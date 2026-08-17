from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


SCHEMA_VERSION = "1.0"
REGISTER_DIR = Path("results/core_gene_analysis/evidence_register")
TSV_NAME = "coordinate_evidence_register.tsv"
JSON_NAME = "coordinate_evidence_register.json"

CORE_TABLES = (
    "gene_model_index.tsv",
    "protein_isoform_index.tsv",
    "exon_protein_map.tsv",
    "domain_features.tsv",
    "tm_features.tsv",
    "exon_domain_boundary_distances.tsv",
)

TSV_COLUMNS = [
    "schema_version", "run_id", "analysis_id", "dataset_id", "register_phase",
    "gene_symbol", "species_id", "transcript_id", "protein_id",
    "model_source", "model_status", "record_type", "record_id",
    "coordinate_system", "coordinate_start_aa", "coordinate_end_aa",
    "coordinate_position_aa", "coordinate_source", "mapping_confidence",
    "exon_id", "exon_number", "cds_start", "cds_end", "phase",
    "domain_id", "domain_name", "domain_source", "domain_start_aa",
    "domain_end_aa", "domain_score", "boundary_type", "signed_distance_aa",
    "absolute_distance_aa", "boundary_category", "qc_state", "qc_detail",
    "run_status", "interpro_qc_status", "pytmhmm_qc_status",
    "provenance_source_file", "provenance_source_row",
]


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read_json(path: Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return {}
    return value if isinstance(value, dict) else {}


def _read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        return []
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    except (OSError, csv.Error):
        return []


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _rel(path: Path, run_dir: Path) -> str:
    try:
        return path.resolve().relative_to(run_dir.resolve()).as_posix()
    except (OSError, ValueError):
        return str(path)


def _first(row: Mapping[str, Any], *names: str) -> str:
    for name in names:
        value = row.get(name)
        if value is not None and str(value).strip() != "":
            return str(value).strip()
    return ""


def _norm_accession(value: str) -> str:
    return (value or "").split(".", 1)[0]


def _record_id(*parts: Any) -> str:
    cleaned = [str(part).strip().replace("\t", " ") for part in parts
               if part is not None and str(part).strip()]
    return ":".join(cleaned)


def _model_index(rows: Sequence[Mapping[str, Any]]) -> Dict[Tuple[str, str], List[Dict[str, str]]]:
    index: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    seen: set[Tuple[str, str, str, str]] = set()
    for row in rows:
        species = _first(row, "species_id")
        protein = _first(row, "protein_id")
        if not protein:
            continue
        model = {
            "transcript_id": _first(row, "transcript_id"),
            "model_source": _first(row, "source"),
            "model_status": _first(row, "model_status", "primary_status"),
        }
        marker = (species, protein, model["transcript_id"], model["model_source"])
        if marker in seen:
            continue
        seen.add(marker)
        index[(species, protein)].append(model)
        norm = _norm_accession(protein)
        if norm and norm != protein:
            index[(species, norm)].append(model)
    return index


def _models_for(index: Mapping[Tuple[str, str], List[Dict[str, str]]],
                species: str, protein: str, transcript: str = "") -> List[Dict[str, str]]:
    models = list(index.get((species, protein), []))
    if not models:
        models = list(index.get((species, _norm_accession(protein)), []))
    if transcript:
        exact = [model for model in models if model.get("transcript_id") == transcript]
        if exact:
            return exact
        return [{"transcript_id": transcript, "model_source": "", "model_status": ""}]
    return models or [{"transcript_id": "", "model_source": "", "model_status": ""}]


def _qc_for_coordinate(*, start: str = "", end: str = "", position: str = "",
                       model_linked: bool = True, pending: bool = False,
                       applicable: bool = True) -> Tuple[str, str]:
    if not applicable:
        return "not_applicable", "This evidence type is not applicable to this coordinate record."
    if pending:
        return "pending_cluster", "Domain-dependent evidence is unavailable before cluster annotation."
    if not (position or (start and end)):
        return "incomplete_coordinate", "A required protein-space coordinate is missing."
    if not model_linked:
        return "missing_model_link", "No transcript/protein model link was available for this record."
    return "available", "Coordinate and its declared source context are available."


def _base_record(context: Mapping[str, str], model: Mapping[str, str],
                 *, species: str, protein: str, transcript: str = "") -> Dict[str, str]:
    return {column: "" for column in TSV_COLUMNS} | {
        "schema_version": SCHEMA_VERSION,
        "run_id": context["run_id"],
        "analysis_id": context["analysis_id"],
        "dataset_id": context["dataset_id"],
        "register_phase": context["register_phase"],
        "gene_symbol": context["gene_symbol"],
        "species_id": species,
        "transcript_id": transcript or model.get("transcript_id", ""),
        "protein_id": protein,
        "model_source": model.get("model_source", ""),
        "model_status": model.get("model_status", ""),
        "coordinate_system": "protein_aa_1_based_inclusive",
        "run_status": context["run_status"],
        "interpro_qc_status": context["interpro_qc_status"],
        "pytmhmm_qc_status": context["pytmhmm_qc_status"],
    }


def _boundary_index(rows: Sequence[Mapping[str, Any]]) -> Dict[Tuple[str, str, str, str], Mapping[str, Any]]:
    index: Dict[Tuple[str, str, str, str], Mapping[str, Any]] = {}
    for row in rows:
        key = (_first(row, "species_id"), _first(row, "protein_id"),
               _first(row, "transcript_id"), _first(row, "boundary_position_aa"))
        index.setdefault(key, row)
    return index


def _attach_domain(record: Dict[str, str], row: Mapping[str, Any]) -> None:
    record.update({
        "domain_id": _first(row, "nearest_domain_id", "nearest_domain_accession", "domain_id",
                            "interpro_accession", "signature_accession"),
        "domain_name": _first(row, "nearest_domain_name", "domain_name", "interpro_name",
                              "signature_name"),
        "domain_source": _first(row, "domain_source", "member_databases", "member_database"),
        "domain_start_aa": _first(row, "nearest_domain_start_aa", "domain_start_aa", "start_aa"),
        "domain_end_aa": _first(row, "nearest_domain_end_aa", "domain_end_aa", "end_aa"),
        "domain_score": _first(row, "score", "score_or_evalue"),
        "boundary_type": _first(row, "nearest_domain_boundary_type", "domain_edge_type",
                                "nearest_edge"),
        "signed_distance_aa": _first(row, "signed_distance_aa"),
        "absolute_distance_aa": _first(row, "absolute_distance_aa", "distance_aa"),
        "boundary_category": _first(row, "category", "classification"),
    })


def build_coordinate_evidence_register(run_dir: Path) -> Dict[str, Any]:
    run_dir = Path(run_dir).resolve()
    core = run_dir / "results" / "core_gene_analysis"
    if not core.is_dir():
        raise FileNotFoundError(f"Core gene-analysis directory not found: {core}")

    run_config = _read_json(run_dir / "run_config.json")
    status = _read_json(run_dir / "status.json")
    core_report = _read_json(core / "core_gene_report.json")
    collection_report = _read_json(core / "core_model_collection_report.json")
    qc_path = run_dir / "results" / "15_domain_architecture" / "post_cluster_qc.json"
    post_qc = _read_json(qc_path)

    tables = {name: _read_tsv(core / name) for name in CORE_TABLES}
    model_rows = tables["gene_model_index.tsv"] or tables["protein_isoform_index.tsv"]
    models = _model_index(model_rows)
    domains = tables["domain_features.tsv"]
    boundaries = tables["exon_domain_boundary_distances.tsv"]

    domain_status = _first(core_report, "domain_status")
    post_cluster = bool(domains) or domain_status in {"complete", "available"}
    context = {
        "run_id": _first(run_config, "run_id", "dataset_id") or run_dir.name,
        "analysis_id": _first(core_report, "analysis_id") or _first(run_config, "analysis_id", "case_study"),
        "dataset_id": _first(core_report, "dataset_id") or _first(run_config, "dataset_id", "run_id") or run_dir.name,
        "register_phase": "post_cluster" if post_cluster else "pre_cluster",
        "gene_symbol": _first(core_report, "gene_symbol") or _first(run_config, "gene_symbol", "gene"),
        "run_status": _first(status, "run_status", "status", "current_step"),
        "interpro_qc_status": _first(post_qc.get("interproscan", {}), "status") or (
            "pending_cluster" if not post_cluster else "not_reported"),
        "pytmhmm_qc_status": _first(post_qc.get("pytmhmm", {}), "status") or (
            "pending_cluster" if not post_cluster else "not_reported"),
    }

    records: List[Dict[str, str]] = []
    boundary_lookup = _boundary_index(boundaries)

    for source_row, row in enumerate(tables["exon_protein_map.tsv"], start=2):
        species, protein, transcript = (_first(row, "species_id"), _first(row, "protein_id"),
                                        _first(row, "transcript_id"))
        for model in _models_for(models, species, protein, transcript):
            start = _first(row, "protein_start_aa")
            end = _first(row, "protein_end_aa")
            record = _base_record(context, model, species=species, protein=protein,
                                  transcript=transcript)
            record.update({
                "record_type": "exon_protein_interval",
                "record_id": _record_id(protein, _first(row, "exon_id") or "exon",
                                        _first(row, "exon_number")),
                "coordinate_start_aa": start, "coordinate_end_aa": end,
                "coordinate_source": _first(row, "source"),
                "mapping_confidence": _first(row, "confidence"),
                "exon_id": _first(row, "exon_id"), "exon_number": _first(row, "exon_number"),
                "cds_start": _first(row, "cds_start"), "cds_end": _first(row, "cds_end"),
                "phase": _first(row, "phase"),
                "provenance_source_file": "results/core_gene_analysis/exon_protein_map.tsv",
                "provenance_source_row": str(source_row),
            })
            boundary = boundary_lookup.get((species, protein, transcript, end))
            if boundary:
                _attach_domain(record, boundary)
            qc_state, qc_detail = _qc_for_coordinate(
                start=start, end=end, model_linked=bool(record["transcript_id"]), pending=not post_cluster)
            record.update({"qc_state": qc_state, "qc_detail": qc_detail})
            records.append(record)

    for source_row, row in enumerate(domains, start=2):
        species, protein = _first(row, "species_id"), _first(row, "protein_id")
        linked = _models_for(models, species, protein)
        for model in linked:
            start, end = _first(row, "start_aa"), _first(row, "end_aa")
            record = _base_record(context, model, species=species, protein=protein)
            record.update({
                "record_type": "domain_interval",
                "record_id": _record_id(protein, _first(row, "domain_id", "interpro_accession"), start, end,
                                        model.get("transcript_id", "")),
                "coordinate_start_aa": start, "coordinate_end_aa": end,
                "coordinate_source": _first(row, "domain_source", "member_databases"),
                "provenance_source_file": "results/core_gene_analysis/domain_features.tsv",
                "provenance_source_row": str(source_row),
            })
            _attach_domain(record, row)
            qc_state, qc_detail = _qc_for_coordinate(
                start=start, end=end, model_linked=bool(record["transcript_id"]))
            record.update({"qc_state": qc_state, "qc_detail": qc_detail})
            records.append(record)

    for source_row, row in enumerate(tables["tm_features.tsv"], start=2):
        species, protein = _first(row, "species_id"), _first(row, "protein_id")
        for model in _models_for(models, species, protein):
            start, end = _first(row, "start_aa"), _first(row, "end_aa")
            record = _base_record(context, model, species=species, protein=protein)
            record.update({
                "record_type": "transmembrane_interval",
                "record_id": _record_id(protein, "tm", start, end, model.get("transcript_id", "")),
                "coordinate_start_aa": start, "coordinate_end_aa": end,
                "coordinate_source": _first(row, "source"),
                "domain_source": "not_applicable",
                "provenance_source_file": "results/core_gene_analysis/tm_features.tsv",
                "provenance_source_row": str(source_row),
            })
            qc_state, qc_detail = _qc_for_coordinate(
                start=start, end=end, model_linked=bool(record["transcript_id"]))
            record.update({"qc_state": qc_state, "qc_detail": qc_detail})
            records.append(record)

    for source_row, row in enumerate(boundaries, start=2):
        species, protein, transcript = (_first(row, "species_id"), _first(row, "protein_id"),
                                        _first(row, "transcript_id"))
        for model in _models_for(models, species, protein, transcript):
            position = _first(row, "boundary_position_aa")
            record = _base_record(context, model, species=species, protein=protein,
                                  transcript=transcript)
            record.update({
                "record_type": "exon_domain_boundary_relation",
                "record_id": _first(row, "exon_boundary_id") or _record_id(protein, "boundary", position),
                "coordinate_position_aa": position,
                "coordinate_source": _first(row, "source"),
                "exon_id": _first(row, "exon_boundary_id"),
                "provenance_source_file": "results/core_gene_analysis/exon_domain_boundary_distances.tsv",
                "provenance_source_row": str(source_row),
            })
            _attach_domain(record, row)
            domain_available = bool(record["domain_id"] or record["domain_start_aa"] or record["domain_end_aa"])
            qc_state, qc_detail = _qc_for_coordinate(
                position=position, model_linked=bool(record["transcript_id"]),
                pending=not post_cluster)
            if qc_state == "available" and not domain_available:
                qc_state = "no_domain_call"
                qc_detail = "Boundary coordinate is available, but no representative domain call was linked."
            record.update({"qc_state": qc_state, "qc_detail": qc_detail})
            records.append(record)

    output_dir = run_dir / REGISTER_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = output_dir / TSV_NAME
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in TSV_COLUMNS}
                         for row in records)

    source_files: List[Dict[str, Any]] = []
    for name in CORE_TABLES:
        path = core / name
        if path.is_file():
            source_files.append({
                "path": _rel(path, run_dir), "sha256": _sha256(path),
                "n_data_rows": len(tables[name]),
            })
    for path in (run_dir / "run_config.json", run_dir / "status.json",
                 core / "core_model_collection_report.json", core / "core_gene_report.json", qc_path):
        if path.is_file():
            source_files.append({"path": _rel(path, run_dir), "sha256": _sha256(path)})

    counts = dict(sorted(Counter(record["record_type"] for record in records).items()))
    payload: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "register_status": "available" if records else "empty_no_coordinate_records",
        "run_id": context["run_id"], "analysis_id": context["analysis_id"],
        "dataset_id": context["dataset_id"], "gene_symbol": context["gene_symbol"],
        "register_phase": context["register_phase"], "generated_at": _now_iso(),
        "coordinate_system": "protein_aa_1_based_inclusive",
        "description": (
            "Additive coordinate-level audit register. Scientific values are copied from the "
            "listed source artefacts; missing or inapplicable evidence is marked explicitly."
        ),
        "counts": {"total_records": len(records), "by_record_type": counts},
        "source_files": source_files,
        "run_provenance": {
            "run_config": run_config,
            "status": status,
            "core_model_collection_report": collection_report,
            "core_gene_report": core_report,
            "post_cluster_qc": post_qc,
        },
        "records": records,
    }
    json_path = output_dir / JSON_NAME
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return payload


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Build a run-local coordinate evidence register.")
    parser.add_argument("--run-dir", type=Path, help="Path to a run directory.")
    parser.add_argument("--run-id", help="Run id below --runs-root (defaults to ./runs).")
    parser.add_argument("--runs-root", type=Path, default=Path("runs"))
    args = parser.parse_args(argv)
    if bool(args.run_dir) == bool(args.run_id):
        parser.error("Provide exactly one of --run-dir or --run-id.")
    run_dir = args.run_dir if args.run_dir else args.runs_root / args.run_id
    payload = build_coordinate_evidence_register(run_dir)
    print(f"OK  coordinate evidence register: {payload['counts']['total_records']} records "
          f"({payload['register_phase']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
