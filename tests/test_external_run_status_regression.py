"""Regression coverage for completed generic runs stored outside the repository."""
from __future__ import annotations

import json
from pathlib import Path

from exondomaincompare.framework.core_run_milestones import evaluate_core_run
from webapp.backend import main as backend


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _completed_external_run(tmp_path: Path) -> Path:
    run = tmp_path / "Application Support" / "ExonDomainCompare" / "runs" / "ptpn11_external"
    core = run / "results" / "core_gene_analysis"
    generic = run / "website_indices" / "generic"

    canonical = {
        "run_id": run.name,
        "analysis_id": "PTPN11_core_only_pilot",
        "gene_symbol": "PTPN11",
        "has_event": False,
        "run_mode": "core_only_pilot",
        "cluster_input_fasta": "run:results/core_gene_analysis/proteins_primary.faa",
        "primary_fasta_path": "run:results/core_gene_analysis/proteins_primary.faa",
        "species_count": 1,
        "species_ids": ["rattus_norvegicus"],
    }
    compatibility = {
        "run_id": run.name,
        "analysis_id": canonical["analysis_id"],
        "gene_symbol": "PTPN11",
        "has_event": False,
    }
    status = {
        "run_id": run.name,
        "status": "results_ready",
        "post_interpro_status": "complete",
        "cluster_analysis_status": "complete",
        "cluster_fetch_status": "complete",
        "next_action": "open_results",
    }
    _write(run / "run.json", json.dumps(canonical))
    _write(run / "run_config.json", json.dumps(compatibility))
    _write(run / "status.json", json.dumps(status))
    _write(run / "gene_config.yaml", "gene_symbol: PTPN11\n")
    _write(run / "species_list.txt", "rattus_norvegicus\n")
    _write(core / "gene_model_index.tsv",
           "species_id\ttranscript_id\tprotein_id\nrattus_norvegicus\tT1\tP1\n")
    _write(core / "protein_isoform_index.tsv",
           "species_id\tprotein_id\tprimary_status\nrattus_norvegicus\tP1\tprimary\n")
    _write(core / "proteins_primary.faa", ">P1\nMPEPTIDE\n")
    _write(core / "exon_protein_map.tsv",
           "species_id\tprotein_id\texon_id\nrattus_norvegicus\tP1\tE1\n")
    _write(core / "synteny_neighbors.tsv", "species_id\tgene\nrattus_norvegicus\tPTPN11\n")
    _write(core / "domain_features.tsv", "protein_id\tstart\tend\nP1\t2\t7\n")
    _write(core / "core_gene_report.json", "{}")
    _write(generic / "domain_architecture_index.json", json.dumps({"available": True}))
    _write(generic / "available_views.json", json.dumps({
        "available_views": {
            "overview": True,
            "gene_explorer": True,
            "domain_architecture": True,
            "exon_domain_boundaries": True,
            "synteny": True,
        }
    }))
    _write(run / "logs" / "cluster_roundtrip.log", "completed\n")
    return run


def test_external_completed_run_uses_the_canonical_record_and_stays_ready(tmp_path: Path):
    run = _completed_external_run(tmp_path)

    report = evaluate_core_run(run)

    assert report["inferred_status"] == "results_ready"
    assert report["suggested_next_action"] == "open_results"
    assert report["cluster_command"] == (
        ".venv/bin/edc cluster roundtrip --run-id ptpn11_external"
    )
    assert report["logs"] == ["run:logs/cluster_roundtrip.log"]


def test_backend_badge_and_action_agree_with_the_completed_external_run(tmp_path: Path):
    model = backend.derive_status_model(_completed_external_run(tmp_path))

    assert model["status"] == "results_ready"
    assert model["status_label"] == "Results ready"
    assert model["post_interpro_status"] == "complete"
    assert model["next_action"] == "open_results"
    assert model["explorable"] is True
