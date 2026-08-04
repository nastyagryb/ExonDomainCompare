from __future__ import annotations

from concurrent.futures import Future, ThreadPoolExecutor
from threading import Event


def test_repeated_file_links_are_checked_once_per_model(monkeypatch, tmp_path):
    from webapp.backend import main

    available = tmp_path / "available.tsv"
    available.write_text("ok\n", encoding="utf-8")
    missing = tmp_path / "missing.tsv"
    calls: list[str] = []

    def resolve(path: str, *, dataset: str):
        calls.append(path)
        assert dataset == "example"
        return available if path.endswith("available.tsv") else missing

    monkeypatch.setattr(main, "_resolve_public_file_path", resolve)
    payload = {
        "source_table": "results/available.tsv",
        "source_tables": {
            "same": "results/available.tsv",
            "missing": "results/missing.tsv",
        },
        "figures": [
            {"format": "tsv", "path": "results/available.tsv"},
            {"format": "tsv", "path": "results/missing.tsv"},
        ],
    }

    cleaned = main._prune_missing_file_links(payload, "example")

    assert calls.count("results/available.tsv") == 1
    assert calls.count("results/missing.tsv") == 1
    assert cleaned["source_table"] == "results/available.tsv"
    assert cleaned["source_tables"] == {"same": "results/available.tsv"}
    assert cleaned["figures"] == [
        {"format": "tsv", "path": "results/available.tsv"}
    ]


def test_concurrent_identical_model_requests_share_one_build(monkeypatch):
    from webapp.backend import main

    owner_started = Event()
    waiter_started = Event()
    release_owner = Event()
    build_calls = 0

    class ObservableFuture(Future):
        def result(self, timeout=None):
            waiter_started.set()
            return super().result(timeout)

    def build(_resolved):
        nonlocal build_calls
        build_calls += 1
        owner_started.set()
        assert release_owner.wait(timeout=5)
        return {"dataset": "example"}

    monkeypatch.setattr(main, "Future", ObservableFuture)
    monkeypatch.setattr(main, "resolve_dataset", lambda selected: selected)
    monkeypatch.setattr(main, "build_canonical_dataset_model", build)
    monkeypatch.setattr(
        main, "_prune_missing_file_links", lambda model, selected: model
    )
    with main._DATASET_MODEL_INFLIGHT_LOCK:
        main._DATASET_MODEL_INFLIGHT.clear()

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            owner = pool.submit(main.current_dataset_model, "example")
            assert owner_started.wait(timeout=5)
            waiter = pool.submit(main.current_dataset_model, "example")
            assert waiter_started.wait(timeout=5)
            release_owner.set()
            assert owner.result(timeout=5) == {"dataset": "example"}
            assert waiter.result(timeout=5) == {"dataset": "example"}
    finally:
        release_owner.set()
        with main._DATASET_MODEL_INFLIGHT_LOCK:
            main._DATASET_MODEL_INFLIGHT.clear()

    assert build_calls == 1
    assert main._DATASET_MODEL_INFLIGHT == {}
