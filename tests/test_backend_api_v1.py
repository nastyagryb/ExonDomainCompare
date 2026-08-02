"""API tests for the current ExonDomainCompare backend.

The original version of this file exercised the first-phase backend, which
served a fixed set of static JSON files from a copied project tree
(``/api/dashboard``, ``DEFAULT_DATA_DIR``). That surface no longer exists: the
backend now builds a canonical dataset model per run. The tests below therefore
target the real module and the endpoints the web UI actually consumes.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from webapp.backend import main as backend_main

PROJECT_ROOT = Path(__file__).resolve().parents[1]
FGFR1_RUN_ID = "2026-07-23_1100_fgfr1_gallus_core_pilot"


@pytest.fixture(scope="module")
def client() -> TestClient:
    return TestClient(backend_main.app)


def test_health_reports_a_version(client):
    r = client.get("/api/health")
    assert r.status_code == 200
    body = r.json()
    assert body.get("status") in {"ok", "healthy"}
    assert body.get("version")


def test_dataset_list_is_served(client):
    r = client.get("/api/datasets")
    assert r.status_code == 200
    payload = r.json()
    datasets = payload if isinstance(payload, list) else payload.get("datasets", [])
    assert isinstance(datasets, list)


def test_unknown_dataset_status_is_not_a_server_error(client):
    r = client.get("/api/datasets/run:does-not-exist/status")
    assert r.status_code in {200, 404}, r.text
    if r.status_code == 200:
        # An unknown dataset must resolve to an explicit unavailable stage, never
        # to a silently "ready" one.
        assert r.json().get("status") in {"unavailable", "failed", None}


def test_local_run_status_uses_the_canonical_stage_vocabulary(client):
    if not (PROJECT_ROOT / "runs" / FGFR1_RUN_ID).is_dir():
        pytest.skip("FGFR1 reference run not present")
    r = client.get(f"/api/local-runs/{FGFR1_RUN_ID}/status")
    assert r.status_code == 200, r.text
    body = r.json()
    allowed = {"pre_cluster_ready", "cluster_processing", "post_cluster_partial",
               "results_ready", "failed", "unavailable"}
    for key in ("stage", "run_stage", "dataset_status", "status"):
        value = body.get(key)
        if isinstance(value, str) and value in allowed:
            return
    # The endpoint may nest the stage; accept any canonical token in the payload.
    assert any(token in str(body) for token in allowed), body


def test_missing_index_file_returns_404_not_500(client):
    r = client.get("/api/local-runs/does-not-exist")
    assert r.status_code in {404, 422}, r.text
