"""Run discovery precedence for bundled examples and user-owned runs."""

from pathlib import Path

import pytest

from exondomaincompare.config import RuntimeConfig
from exondomaincompare.runs.registry import RunCollisionError, resolve_run_record


RUN_ID = "2026-08-02_1910_ptpn11_3species"


def _config(root: Path) -> RuntimeConfig:
    return RuntimeConfig(
        repository_root=root,
        runs_root=root / "private-runs",
        local_profile_name="default",
        lrz_profile_name="lrz",
        local={},
        lrz={},
        config_source="test",
    )


def test_user_run_shadows_same_named_bundled_example(tmp_path):
    bundled = tmp_path / "datasets" / "runs" / RUN_ID
    local = tmp_path / "private-runs" / RUN_ID
    bundled.mkdir(parents=True)
    local.mkdir(parents=True)

    record = resolve_run_record(_config(tmp_path), RUN_ID)

    assert record is not None
    assert record.path == local.resolve()
    assert record.kind == "canonical"
    assert record.read_only is False


def test_two_user_roots_still_report_a_real_collision(tmp_path):
    (tmp_path / "datasets" / "runs" / RUN_ID).mkdir(parents=True)
    (tmp_path / "private-runs" / RUN_ID).mkdir(parents=True)
    (tmp_path / "runs" / RUN_ID).mkdir(parents=True)

    with pytest.raises(RunCollisionError) as error:
        resolve_run_record(_config(tmp_path), RUN_ID)

    assert {row.kind for row in error.value.candidates} == {
        "canonical", "repository_legacy",
    }
