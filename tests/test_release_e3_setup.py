"""Release E3 contracts for setup and remote tool onboarding."""
from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from exondomaincompare import cluster_setup
from exondomaincompare.cluster_setup import (
    INTERPROSCAN_VERSION,
    MANAGED_MARKER,
    PYTMHMM_VERSION,
    _installer_script,
    initial_user_config,
    remote_install_plan,
    remote_preflight,
    render_managed_config,
    write_managed_config,
)
from exondomaincompare.config import (
    ConfigurationError,
    RuntimeConfig,
    load_config,
    remote_shell_path,
    user_config_path,
)
from exondomaincompare.cli import main as cli_main
from exondomaincompare.runs.registry import _normalize


ROOT = Path(__file__).resolve().parents[1]


def test_roundtrip_stops_before_ssh_when_cluster_setup_is_missing():
    cluster_scripts = ROOT / "scripts" / "interpro_cluster"
    if str(cluster_scripts) not in sys.path:
        sys.path.insert(0, str(cluster_scripts))
    import run_cluster_roundtrip

    class MissingProfile:
        @staticmethod
        def require_cluster():
            raise ConfigurationError("missing profile fields: user, host")

    with pytest.raises(SystemExit, match="edc cluster configure"):
        run_cluster_roundtrip._require_cluster_profile(MissingProfile())

    source = (cluster_scripts / "run_cluster_roundtrip.py").read_text(encoding="utf-8")
    main = source.split("def main(", 1)[1]
    assert main.index("_require_cluster_profile(RUNTIME_CONFIG)") \
        < main.index("rt.open_ssh_master()")


def _config(tmp_path: Path) -> RuntimeConfig:
    return RuntimeConfig(
        repository_root=ROOT,
        runs_root=tmp_path / "runs",
        local_profile_name="default",
        lrz_profile_name="lrz",
        local={
            "python": "python", "node": "node", "npm": "npm",
            "datasets": "datasets", "mafft": "mafft",
            "ssh": "/usr/bin/ssh", "scp": "/usr/bin/scp",
        },
        lrz={
            "user": "TEST_USER", "host": "TEST_HOST",
            "remote_root": "~/ExonDomainCompare/runs",
            "partition": "TEST_PARTITION", "account": "",
            "managed_tool_root": "~/.local/share/ExonDomainCompare",
            "ssh_options": [], "scp_options": [],
            "interproscan": {
                "launcher": "~/.local/share/ExonDomainCompare/tools/"
                f"interproscan-{INTERPROSCAN_VERSION}/interproscan.sh",
                "environment": "~/.local/share/ExonDomainCompare/envs/interproscan",
            },
            "pytmhmm": {
                "launcher": "~/.local/share/ExonDomainCompare/envs/pytmhmm/bin/pyTMHMM",
                "environment": "~/.local/share/ExonDomainCompare/envs/pytmhmm",
                "python": "~/.local/share/ExonDomainCompare/envs/pytmhmm/bin/python",
            },
        },
        config_source="test",
    )


def test_setup_and_discovery_use_the_same_personal_path(tmp_path: Path):
    env = {"EDC_DATA_DIR": str(tmp_path / "data")}
    expected = tmp_path / "data" / "config" / "config.toml"
    assert user_config_path(env) == expected
    expected.parent.mkdir(parents=True)
    expected.write_text(initial_user_config(), encoding="utf-8")
    cfg = load_config(repository_root=ROOT, env=env)
    assert cfg.config_source == "explicit/user configuration"


def test_managed_configuration_is_private_portable_and_loadable(tmp_path: Path):
    cfg = _config(tmp_path)
    text = render_managed_config(
        cfg, user="AB12CD", host="cluster.example", partition="cpu")
    assert MANAGED_MARKER in text
    assert "password" not in text.lower()
    assert str(Path.home()) not in text
    target = write_managed_config(
        cfg, output=tmp_path / "config.toml",
        user="AB12CD", host="cluster.example", partition="cpu")
    assert target.stat().st_mode & 0o077 == 0
    loaded = load_config(repository_root=ROOT, config_path=target, env={})
    assert loaded.ssh_target == "AB12CD@cluster.example"
    assert loaded.lrz["remote_root"] == "~/ExonDomainCompare/runs"


def test_existing_remote_tools_can_be_configured_without_managed_install(tmp_path: Path):
    text = render_managed_config(
        _config(tmp_path), user="AB12CD", host="cluster.example", partition="cpu",
        interproscan_launcher="interproscan.sh", interproscan_module="ipr/5.78",
        pytmhmm_launcher="pyTMHMM", pytmhmm_module="pytmhmm/1.3.6",
        pytmhmm_python="python3")
    assert 'launcher = "interproscan.sh"' in text
    assert 'module = "ipr/5.78"' in text
    assert 'module = "pytmhmm/1.3.6"' in text
    assert 'python = "python3"' in text


def test_managed_registry_roots_follow_the_current_clone(tmp_path: Path):
    repository = tmp_path / "relocated"
    (repository / "runs").mkdir(parents=True)
    (repository / "datasets" / "runs").mkdir(parents=True)
    cfg = replace(_config(tmp_path), repository_root=repository)
    normalized = _normalize(cfg, {
        "registry_version": "2.0",
        "roots": [
            {"id": "configured-runs", "path": "/OLD/RUNS"},
            {"id": "repository-legacy-runs", "path": "/OLD/CLONE/runs"},
            {"id": "bundled-release-datasets", "path": "/OLD/CLONE/datasets/runs"},
            {"id": "user-reference", "path": "/USER/REFERENCE", "read_only": True},
        ],
        "runs": [],
    })
    roots = {row["id"]: row for row in normalized["roots"]}
    assert roots["configured-runs"]["path"] == str(cfg.paths.runs)
    assert roots["repository-legacy-runs"]["path"] == str(repository / "runs")
    assert roots["bundled-release-datasets"]["path"] == str(
        repository / "datasets" / "runs")
    assert roots["user-reference"]["path"] == "/USER/REFERENCE"


def test_manual_configuration_is_not_overwritten_without_confirmation(tmp_path: Path):
    target = tmp_path / "manual.toml"
    target.write_text('schema_version = "1.0"\n', encoding="utf-8")
    with pytest.raises(ConfigurationError, match="Refusing"):
        write_managed_config(
            _config(tmp_path), output=target,
            user="AB12CD", host="cluster.example", partition="cpu")


def test_remote_home_paths_expand_without_embedding_a_username():
    assert remote_shell_path("~/ExonDomainCompare/runs") == \
        '"$HOME"/' + "ExonDomainCompare/runs"
    assert "TEST_USER" not in remote_shell_path("~/ExonDomainCompare/runs")


def test_remote_preflight_is_read_only_and_submits_no_job(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    captured: dict[str, object] = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["script"] = kwargs["input"]
        output = "\n".join(
            f"EDC:{key}={value}" for key, value in {
                "linux": "yes", "home_writable": "yes",
                "remote_root_writable": "yes", "sbatch": "yes",
                "squeue": "yes", "sacct": "yes", "partition": "yes",
                "interproscan": "yes", "pytmhmm": "yes",
                "downloader": "yes", "archive_tools": "yes",
                "environment_manager": "yes", "free_kb": "50000000",
                "tool_root_writable": "yes", "tool_free_kb": "50000000",
                "home": "/remote/home/AB12CD",
            }.items())
        return subprocess.CompletedProcess(command, 0, stdout=output)

    monkeypatch.setattr(cluster_setup.subprocess, "run", fake_run)
    report = remote_preflight(_config(tmp_path), redact_paths=True)
    assert report["ready_for_cluster_runs"] is True
    assert report["job_submitted"] is False
    assert report["remote_mutation"] is False
    assert report["remote_home"] == "<REMOTE_HOME>"
    assert report["expected_versions"]["interproscan"] == INTERPROSCAN_VERSION
    assert report["expected_versions"]["pytmhmm"] == PYTMHMM_VERSION
    script = str(captured["script"])
    assert not any(
        line.lstrip().startswith("sbatch ") for line in script.splitlines())
    assert "mkdir -p" not in script
    assert not any(
        line.lstrip().startswith(("curl ", "wget "))
        for line in script.splitlines())


def test_remote_preflight_resolves_launchers_inside_configured_environments(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    captured: dict[str, str] = {}

    def fake_run(command, **kwargs):
        captured["script"] = kwargs["input"]
        values = {
            "linux": "yes", "home_writable": "yes",
            "remote_root_writable": "yes", "sbatch": "yes",
            "squeue": "yes", "sacct": "yes", "partition": "yes",
            "interproscan": "yes", "pytmhmm": "yes",
            "downloader": "yes", "archive_tools": "yes",
            "environment_manager": "yes", "free_kb": "50000000",
            "tool_root_writable": "yes", "tool_free_kb": "50000000",
            "home": "/remote/home/AB12CD",
        }
        output = "\n".join(f"EDC:{key}={value}" for key, value in values.items())
        return subprocess.CompletedProcess(command, 0, stdout=output)

    monkeypatch.setattr(cluster_setup.subprocess, "run", fake_run)
    remote_preflight(_config(tmp_path))
    script = captured["script"]
    assert 'PATH="$INTERPRO_ENV/bin:$PATH" command_ok "$INTERPRO"' in script
    assert 'PATH="$PYTMHMM_ENV/bin:$PATH" command_ok "$PYTMHMM"' in script


@pytest.mark.parametrize(
    ("missing", "expected", "forbidden"),
    [
        ({"partition"}, "Review the reported LRZ profile", "tools install"),
        ({"interproscan"}, "tools install --tool interproscan", "--confirm"),
        ({"interproscan", "pytmhmm"}, "tools install --tool all", "--confirm"),
    ],
)
def test_remote_doctor_recommends_installation_only_for_missing_tools(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
        missing: set[str], expected: str, forbidden: str):
    keys = {
        "linux", "home_writable", "remote_root_writable", "sbatch", "squeue",
        "sacct", "partition", "interproscan", "pytmhmm", "tool_root_writable",
        "downloader", "archive_tools", "environment_manager",
    }

    def fake_run(command, **_kwargs):
        values = {key: "no" if key in missing else "yes" for key in keys}
        values.update({"free_kb": "50000000", "tool_free_kb": "50000000",
                       "home": "/remote/home/AB12CD"})
        output = "\n".join(f"EDC:{key}={value}" for key, value in values.items())
        return subprocess.CompletedProcess(command, 0, stdout=output)

    monkeypatch.setattr(cluster_setup.subprocess, "run", fake_run)
    action = remote_preflight(_config(tmp_path))["next_action"]
    assert expected in action
    assert forbidden not in action


def test_remote_installer_is_explicit_pinned_and_idempotent(tmp_path: Path):
    cfg = _config(tmp_path)
    plan = remote_install_plan(cfg, "all")
    assert plan["confirmation_required"] is True
    assert plan["network_contacted"] is False
    assert plan["submits_cluster_jobs"] is False
    script = _installer_script(cfg, "all")
    assert INTERPROSCAN_VERSION in script
    assert f'pyTMHMM=={PYTMHMM_VERSION}' in script
    assert "md5sum -c" in script
    assert "30 GB free" in script
    assert "install -y -p" in script
    assert "version('pyTMHMM')" in script
    assert "sbatch" not in script
    assert ".bashrc" not in script and ".profile" not in script
    assert "rm -rf" not in script


def test_local_setup_is_reproducible_and_start_does_not_install_on_demand():
    setup = (ROOT / "scripts" / "setup_local.sh").read_text(encoding="utf-8")
    start = (ROOT / "scripts" / "start_local.sh").read_text(encoding="utf-8")
    backend = (ROOT / "webapp" / "start_backend.sh").read_text(encoding="utf-8")
    frontend = (ROOT / "webapp" / "start_frontend.sh").read_text(encoding="utf-8")
    assert "python3.13" in setup and "npm ci" in setup
    assert "pip install" in setup and "-e '.[test,render,synteny]'" in setup
    assert "constraints-py313-tested.txt" in setup
    assert "pip install" not in start
    assert "pip install" not in backend
    assert "npm install" not in frontend


def test_interrupted_remote_authentication_has_a_clean_exit(monkeypatch, capsys):
    monkeypatch.setattr(
        "exondomaincompare.cli.cluster_doctor",
        lambda _args: (_ for _ in ()).throw(KeyboardInterrupt()))
    assert cli_main(["cluster", "doctor"]) == 130
    assert "no credentials were stored" in capsys.readouterr().err
