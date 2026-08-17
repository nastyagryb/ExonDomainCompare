from __future__ import annotations

import csv
import json
from pathlib import Path

from exondomaincompare.framework.coordinate_evidence_register import (
    JSON_NAME,
    REGISTER_DIR,
    TSV_NAME,
    build_coordinate_evidence_register,
)


def _write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _run_fixture(tmp_path: Path, *, post_cluster: bool) -> Path:
    run = tmp_path / "runs" / "test_run"
    core = run / "results" / "core_gene_analysis"
    core.mkdir(parents=True)
    (run / "run_config.json").write_text(json.dumps({
        "run_id": "test_run", "dataset_id": "test_run", "analysis_id": "GENE_core",
        "gene_symbol": "GENE",
    }))
    (run / "status.json").write_text(json.dumps({
        "status": "results_ready" if post_cluster else "cluster_required",
    }))
    (core / "core_gene_report.json").write_text(json.dumps({
        "analysis_id": "GENE_core", "dataset_id": "run:test_run", "gene_symbol": "GENE",
        "domain_status": "complete" if post_cluster else "pending_cluster",
    }))
    (core / "core_model_collection_report.json").write_text(json.dumps({
        "source": {"assembly_accession": "GCF_TEST"},
    }))
    _write_tsv(core / "gene_model_index.tsv",
               ["analysis_id", "gene_symbol", "species_id", "transcript_id", "protein_id",
                "source", "model_status"], [{
                    "analysis_id": "GENE_core", "gene_symbol": "GENE", "species_id": "species_a",
                    "transcript_id": "TX1", "protein_id": "P1", "source": "ncbi_gff",
                    "model_status": "protein_coding",
                }])
    _write_tsv(core / "protein_isoform_index.tsv",
               ["species_id", "protein_id", "transcript_id"], [])
    _write_tsv(core / "exon_protein_map.tsv",
               ["species_id", "protein_id", "transcript_id", "exon_id", "exon_number",
                "cds_start", "cds_end", "protein_start_aa", "protein_end_aa", "phase",
                "confidence", "source"], [{
                    "species_id": "species_a", "protein_id": "P1", "transcript_id": "TX1",
                    "exon_id": "P1:cds1", "exon_number": 1, "cds_start": 100, "cds_end": 129,
                    "protein_start_aa": 1, "protein_end_aa": 10, "phase": 0,
                    "confidence": "gff_cds_derived", "source": "ncbi_gff",
                }])
    domain_rows = [{
        "species_id": "species_a", "protein_id": "P1", "domain_source": "PFAM",
        "domain_id": "IPR_TEST", "domain_name": "Test domain", "start_aa": 2,
        "end_aa": 10, "score": "1e-20",
    }] if post_cluster else []
    _write_tsv(core / "domain_features.tsv",
               ["species_id", "protein_id", "domain_source", "domain_id", "domain_name",
                "start_aa", "end_aa", "score"], domain_rows)
    _write_tsv(core / "tm_features.tsv",
               ["species_id", "protein_id", "start_aa", "end_aa", "source"], [])
    boundary_rows = [{
        "analysis_id": "GENE_core", "gene_symbol": "GENE", "species_id": "species_a",
        "protein_id": "P1", "transcript_id": "TX1", "exon_boundary_id": "P1:cds1_end",
        "boundary_position_aa": 10, "nearest_domain_id": "IPR_TEST",
        "nearest_domain_name": "Test domain", "nearest_domain_start_aa": 2,
        "nearest_domain_end_aa": 10, "nearest_domain_boundary_type": "domain_end",
        "signed_distance_aa": 0, "absolute_distance_aa": 0, "distance_aa": 0,
        "category": "exact_edge", "source": "core_post_interpro",
    }] if post_cluster else []
    _write_tsv(core / "exon_domain_boundary_distances.tsv",
               ["analysis_id", "gene_symbol", "species_id", "protein_id", "transcript_id",
                "exon_boundary_id", "boundary_position_aa", "nearest_domain_id",
                "nearest_domain_name", "nearest_domain_start_aa", "nearest_domain_end_aa",
                "nearest_domain_boundary_type", "signed_distance_aa", "absolute_distance_aa",
                "distance_aa", "category", "source"], boundary_rows)
    if post_cluster:
        qc = run / "results" / "15_domain_architecture" / "post_cluster_qc.json"
        qc.parent.mkdir(parents=True)
        qc.write_text(json.dumps({
            "interproscan": {"status": "valid"}, "pytmhmm": {"status": "valid"},
        }))
    return run


def test_precluster_register_marks_domain_evidence_pending(tmp_path: Path) -> None:
    run = _run_fixture(tmp_path, post_cluster=False)
    payload = build_coordinate_evidence_register(run)
    assert payload["register_phase"] == "pre_cluster"
    assert payload["counts"] == {
        "total_records": 1, "by_record_type": {"exon_protein_interval": 1},
    }
    row = payload["records"][0]
    assert row["transcript_id"] == "TX1"
    assert row["protein_id"] == "P1"
    assert row["mapping_confidence"] == "gff_cds_derived"
    assert row["qc_state"] == "pending_cluster"
    assert row["provenance_source_file"].endswith("exon_protein_map.tsv")
    assert (run / REGISTER_DIR / TSV_NAME).is_file()
    assert (run / REGISTER_DIR / JSON_NAME).is_file()


def test_postcluster_register_links_domain_boundary_qc_and_provenance(tmp_path: Path) -> None:
    run = _run_fixture(tmp_path, post_cluster=True)
    payload = build_coordinate_evidence_register(run)
    assert payload["register_phase"] == "post_cluster"
    assert payload["counts"]["by_record_type"] == {
        "domain_interval": 1,
        "exon_domain_boundary_relation": 1,
        "exon_protein_interval": 1,
    }
    assert all(row["transcript_id"] == "TX1" for row in payload["records"])
    assert all(row["protein_id"] == "P1" for row in payload["records"])
    assert all(row["qc_state"] == "available" for row in payload["records"])
    exon = next(row for row in payload["records"]
                if row["record_type"] == "exon_protein_interval")
    assert exon["domain_id"] == "IPR_TEST"
    relation = next(row for row in payload["records"]
                    if row["record_type"] == "exon_domain_boundary_relation")
    assert relation["boundary_category"] == "exact_edge"
    assert relation["absolute_distance_aa"] == "0"
    assert payload["run_provenance"]["post_cluster_qc"]["interproscan"]["status"] == "valid"
    assert all(source["sha256"] for source in payload["source_files"])
