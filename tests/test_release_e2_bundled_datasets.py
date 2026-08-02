"""Release E2 contracts for bundled read-only example datasets."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest

from exondomaincompare.config import RuntimeConfig
from exondomaincompare.runs.registry import required_roots


ROOT = Path(__file__).resolve().parents[1]
DATASETS = ROOT / "datasets"
PTPN11_RUN_ID = "2026-08-02_1910_ptpn11_3species"


def test_repository_bundled_run_root_is_always_read_only(tmp_path):
    (tmp_path / "datasets" / "runs").mkdir(parents=True)
    config = RuntimeConfig(
        repository_root=tmp_path,
        runs_root=tmp_path / "private-runs",
        local_profile_name="default",
        lrz_profile_name="lrz",
        local={},
        lrz={},
        config_source="test",
    )

    roots = {row["id"]: row for row in required_roots(config)}

    assert roots["configured-runs"]["read_only"] is False
    assert roots["bundled-release-datasets"]["read_only"] is True
    assert Path(roots["bundled-release-datasets"]["path"]) == tmp_path / "datasets" / "runs"


def _bundled_or_skip() -> None:
    if not DATASETS.is_dir():
        pytest.skip("bundled release datasets exist only in the clean release candidate")


def test_bundled_dataset_registry_has_exact_release_scope():
    _bundled_or_skip()
    registry = json.loads((DATASETS / "registry.json").read_text(encoding="utf-8"))

    assert registry["default_dataset"] == "example"
    assert registry["datasets"] == [
        {"id": "example", "path": "fgfr2_30_species", "read_only": True},
        {"id": f"run:{PTPN11_RUN_ID}", "path": f"runs/{PTPN11_RUN_ID}",
         "read_only": True},
    ]


def test_bundled_dataset_checksums_are_complete_and_current():
    _bundled_or_skip()
    lines = (DATASETS / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
    assert len(lines) > 100
    for line in lines:
        digest, relative = line.split("  ", 1)
        path = DATASETS / relative
        assert path.is_file(), relative
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, relative


def test_bundled_dataset_semantics_and_private_state_are_explicit():
    _bundled_or_skip()
    fgfr2 = json.loads((DATASETS / "fgfr2_30_species" / "dataset.json").read_text())
    ptpn11 = json.loads((DATASETS / "runs" / PTPN11_RUN_ID / "dataset.json").read_text())
    status = json.loads((DATASETS / "runs" / PTPN11_RUN_ID / "status.json").read_text())

    assert fgfr2["species_count"] == 30
    assert fgfr2["support_level"] == "validated_event_analysis"
    assert fgfr2["read_only"] is True
    assert ptpn11["species_count"] == 3
    assert ptpn11["scientific_semantics"] == "exploratory_not_validated"
    assert ptpn11["read_only"] is True
    assert status["read_only"] is True
    assert "cluster_jobs" not in status
    assert "cluster_status_detail" not in status


def test_bundled_datasets_do_not_expose_local_machine_paths():
    _bundled_or_skip()
    forbidden = ("/Users/", "/home/", "/private/", "/dss/", "/gpfs/")

    for path in DATASETS.rglob("*"):
        if path.is_file() and path.suffix.lower() in {
            ".csv", ".faa", ".fasta", ".fa", ".gff", ".gff3", ".json",
            ".md", ".svg", ".tsv", ".txt", ".yaml", ".yml",
        }:
            text = path.read_text(encoding="utf-8", errors="replace")
            assert not any(token in text for token in forbidden), path


def test_backend_discovers_both_bundled_datasets_without_a_live_run_root():
    _bundled_or_skip()
    from webapp.backend import main

    datasets = {row["id"]: row for row in main.list_datasets()["datasets"]}
    assert datasets["example"]["read_only"] is True
    assert datasets["example"]["gene_symbol"] == "FGFR2"
    ptpn11 = datasets[f"run:{PTPN11_RUN_ID}"]
    assert ptpn11["read_only"] is True
    assert ptpn11["bundled_example"] is True
    assert ptpn11["gene_symbol"] == "PTPN11"
    assert ptpn11["status"] == "results_ready"
    assert ptpn11["available_views"]["figure_gallery"] is True


def test_bundled_examples_do_not_appear_as_user_owned_runs():
    _bundled_or_skip()
    from webapp.backend import main

    run_ids = {row["run_id"] for row in main.local_runs()}

    assert PTPN11_RUN_ID not in run_ids


def test_homepage_excludes_bundled_examples_from_my_runs():
    start_page = (ROOT / "webapp" / "frontend" / "src" / "pages"
                  / "StartPage.jsx").read_text(encoding="utf-8")

    assert 'd.kind === "run" && !d.bundled_example' in start_page


def test_bundled_models_preserve_the_accepted_scientific_scope():
    _bundled_or_skip()
    from webapp.backend import main
    from webapp.backend.canonical_dataset import build_canonical_dataset_model

    fgfr2 = build_canonical_dataset_model(main.resolve_dataset("example"))
    ptpn11 = build_canonical_dataset_model(
        main.resolve_dataset(f"run:{PTPN11_RUN_ID}"))

    assert len(fgfr2["species"]) == 30
    assert all(fgfr2["available_views"].values())
    assert len(fgfr2["figures"]["figures"]) == 232
    assert len(fgfr2["downloads"]) == 12
    assert len(ptpn11["species"]) == 3
    assert len(ptpn11["figures"]["figures"]) == 60
    assert len(ptpn11["downloads"]["items"]) == 34
    candidates = [
        candidate
        for species in ptpn11["candidate_evidence"]["species"]
        for candidate in species["candidates"]
    ]
    observed = {
        (row["species_id"], row["candidate_id"].split(":")[-2],
         row["aa_start"], row["aa_end"], row["confidence_class"], row["status"])
        for row in candidates
    }
    assert observed == {
        ("homo_sapiens", "indel", 465, 592, "high", "exploratory_not_validated"),
        ("homo_sapiens", "indel", 46, 46, "medium", "exploratory_not_validated"),
        ("homo_sapiens", "indel", 408, 411, "medium", "exploratory_not_validated"),
        ("mus_musculus", "indel", 409, 412, "low", "exploratory_not_validated"),
        ("rattus_norvegicus", "indel", 409, 412, "low", "exploratory_not_validated"),
    }


def test_every_bundled_gallery_file_and_download_resolves():
    _bundled_or_skip()
    from webapp.backend import main
    from webapp.backend.canonical_dataset import build_canonical_dataset_model

    fgfr2 = build_canonical_dataset_model(main.resolve_dataset("example"))
    ptpn_descriptor = main.resolve_dataset(f"run:{PTPN11_RUN_ID}")
    ptpn11 = build_canonical_dataset_model(ptpn_descriptor)
    checked = 0

    for card in fgfr2["figures"]["figures"]:
        modes = card.get("modes") or [card]
        for mode in modes:
            for value in (mode.get("formats") or {}).values():
                assert main._resolve_public_file_path(value).is_file(), value
                checked += 1
            thumbnail = mode.get("thumbnail")
            if thumbnail:
                assert main._resolve_public_file_path(thumbnail).is_file(), thumbnail
                checked += 1
    for item in fgfr2["downloads"]:
        assert main._resolve_public_file_path(item["path"]).is_file(), item["path"]
        checked += 1

    run_base = Path(ptpn_descriptor["run_base"])
    for card in ptpn11["figures"]["figures"]:
        for key in ("png_url", "svg_url", "pdf_url", "table_url"):
            url = card.get(key)
            if not url:
                continue
            relative = parse_qs(urlparse(url).query).get("path", [""])[0]
            assert relative and (run_base / relative).is_file(), url
            checked += 1
    for item in ptpn11["downloads"]["items"]:
        assert main._resolve_public_file_path(item["path"]).is_file(), item["path"]
        checked += 1
    assert checked > 500
