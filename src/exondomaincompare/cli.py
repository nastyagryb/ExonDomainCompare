#!/usr/bin/env python3
"""Unified ExonDomainCompare command-line entry point."""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import os
from pathlib import Path

from exondomaincompare.config import (
    ConfigurationError,
    discover_repository_root,
    load_config,
    user_config_path,
)
from exondomaincompare.config_cli import doctor_report
from exondomaincompare.cluster_setup import (
    initial_user_config,
    install_remote_tools,
    remote_install_plan,
    remote_preflight,
    write_managed_config,
)
from exondomaincompare.runs.migration import MigrationError, MigrationService
from exondomaincompare.runs.registry import discover_runs, write_initial_registry

ROOT = discover_repository_root(Path(__file__))


def _common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config")
    parser.add_argument("--local-profile")
    parser.add_argument("--lrz-profile")


def _config(args: argparse.Namespace):
    return load_config(
        repository_root=ROOT,
        config_path=getattr(args, "config", None),
        local_profile=getattr(args, "local_profile", None),
        lrz_profile=getattr(args, "lrz_profile", None),
    )


def setup(args: argparse.Namespace) -> int:
    cfg = _config(args)
    for path in cfg.paths.setup_roots():
        path.mkdir(parents=True, exist_ok=True)
    template = cfg.paths.config / "config.example.toml"
    if not template.exists():
        shutil.copyfile(ROOT / "config" / "exondomain.example.toml", template)
    personal = Path(args.config).expanduser() if args.config else user_config_path()
    configuration_created = False
    if not personal.exists():
        personal.parent.mkdir(parents=True, exist_ok=True)
        personal.write_text(initial_user_config(), encoding="utf-8")
        os.chmod(personal, 0o600)
        configuration_created = True
    registry = write_initial_registry(cfg)
    print(json.dumps({
        "schema_version": "1.0",
        "created_or_verified_roots": len(cfg.paths.setup_roots()),
        "configuration_template": str(template),
        "configuration_file": str(personal),
        "configuration_created": configuration_created,
        "registry": str(registry),
        "runs_moved": False,
        "data_downloaded": False,
        "network_contacted": False,
    }, indent=2))
    return 0


def doctor(args: argparse.Namespace) -> int:
    print(json.dumps(doctor_report(_config(args), redact_paths=args.redact_paths), indent=2))
    return 0


def profile(args: argparse.Namespace) -> int:
    cfg = _config(args)
    print(json.dumps({
        "schema_version": "1.0",
        "local_profile": cfg.local_profile_name,
        "lrz_profile": cfg.lrz_profile_name,
        "identity": cfg.public_identity(),
    }, indent=2))
    return 0


def serve(args: argparse.Namespace) -> int:
    cfg = _config(args)
    runtime = cfg.local_python()
    command = [
        runtime.selected, "-m", "uvicorn", "webapp.backend.main:app",
        "--host", args.host, "--port", str(args.port),
    ]
    if args.dry_run:
        print(json.dumps({
            "command": cfg.command(["python", "-m", "uvicorn",
                                    "webapp.backend.main:app",
                                    "--host", args.host, "--port", str(args.port)]),
            "network_contacted": False,
            "process_started": False,
            "python_runtime": runtime.report(redact_paths=True),
        }, indent=2))
        return 0
    import subprocess
    return subprocess.run(command, cwd=cfg.repository_root).returncode


def cluster_roundtrip(args: argparse.Namespace) -> int:
    cluster_dir = ROOT / "scripts" / "interpro_cluster"
    if str(cluster_dir) not in sys.path:
        sys.path.insert(0, str(cluster_dir))
    from run_cluster_roundtrip import main as roundtrip_main
    forwarded = ["--run-id", args.run_id]
    for name, flag in (
        ("config", "--config"), ("local_profile", "--local-profile"),
        ("lrz_profile", "--lrz-profile"),
    ):
        value = getattr(args, name)
        if value:
            forwarded += [flag, value]
    if args.dry_run:
        forwarded.append("--dry-run")
    roundtrip_main(forwarded)
    return 0


def cluster_configure(args: argparse.Namespace) -> int:
    target = Path(args.config).expanduser() if args.config else user_config_path()
    cfg = load_config(
        repository_root=ROOT,
        config_path=target if target.is_file() else None,
        local_profile=args.local_profile,
        lrz_profile=args.lrz_profile,
    )
    written = write_managed_config(
        cfg, output=target, replace=args.replace_config,
        user=args.user, host=args.host, partition=args.partition,
        account=args.account, remote_root=args.remote_root,
        tool_root=args.tool_root,
        interproscan_launcher=args.interproscan_launcher,
        interproscan_module=args.interproscan_module,
        interproscan_environment=args.interproscan_environment,
        pytmhmm_launcher=args.pytmhmm_launcher,
        pytmhmm_module=args.pytmhmm_module,
        pytmhmm_environment=args.pytmhmm_environment,
        pytmhmm_python=args.pytmhmm_python,
    )
    print(json.dumps({
        "schema_version": "1.0",
        "configuration_file": str(written),
        "credentials_stored": False,
        "network_contacted": False,
        "next_command": "edc cluster doctor",
    }, indent=2))
    return 0


def cluster_doctor(args: argparse.Namespace) -> int:
    print(json.dumps(
        remote_preflight(_config(args), redact_paths=args.redact_paths), indent=2))
    return 0


def cluster_tools_install(args: argparse.Namespace) -> int:
    cfg = _config(args)
    if not args.confirm:
        report = remote_install_plan(cfg, args.tool)
        report["next_command"] = (
            f"edc cluster tools install --tool {args.tool} --confirm")
        print(json.dumps(report, indent=2))
        return 0
    return install_remote_tools(cfg, tool=args.tool, confirmed=True)


def cluster_init(args: argparse.Namespace) -> int:
    if not args.dry_run:
        raise SystemExit(
            "Remote initialization is reserved for Phase E; use --dry-run to inspect the plan.")
    cfg = _config(args)
    cfg.require_cluster()
    print(json.dumps(cfg.remote_directory_plan(args.run_id), indent=2))
    return 0


def runs_list(args: argparse.Namespace) -> int:
    cfg = _config(args)
    records, collisions = discover_runs(cfg)
    print(json.dumps({
        "schema_version": "1.0",
        "runs": [{
            "run_id": row.run_id,
            "dataset_id": row.dataset_id,
            "root_id": row.root_id,
            "kind": row.kind,
            "read_only": row.read_only,
            "explicit": row.explicit,
        } for row in records],
        "collisions": [{
            "run_id": run_id,
            "roots": [row.root_id for row in rows],
            "requires_explicit_binding": True,
        } for run_id, rows in sorted(collisions.items())],
        "private_paths_included": False,
    }, indent=2))
    return 0


def migrate_legacy_runs(args: argparse.Namespace) -> int:
    cfg = _config(args)
    service = MigrationService(cfg)
    if args.rollback_journal:
        result = service.rollback(args.rollback_journal)
    elif args.mode == "register":
        result = service.register(
            source_root=Path(args.source_root),
            run_ids=args.run_ids,
            read_only=not args.allow_registered_writes,
            dry_run=args.dry_run,
        )
    elif args.mode == "copy":
        result = service.copy(
            source_root=Path(args.source_root),
            run_ids=args.run_ids,
            destination_root=Path(args.destination).expanduser()
            if args.destination else None,
            dry_run=args.dry_run,
        )
    else:
        if not args.copy_journal:
            raise MigrationError("--copy-journal is required for move mode.")
        result = service.move(
            copy_journal_id=args.copy_journal,
            confirmed=args.confirm_move,
            dry_run=args.dry_run,
        )
    print(json.dumps({
        "schema_version": "1.0",
        "operation": "legacy-runs",
        "mode": args.mode,
        "dry_run": args.dry_run,
        "network_contacted": False,
        "scientific_tools_executed": False,
        "result": result,
    }, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="edc")
    commands = parser.add_subparsers(dest="command", required=True)
    setup_parser = commands.add_parser("setup")
    _common(setup_parser)
    setup_parser.set_defaults(handler=setup)
    doctor_parser = commands.add_parser("doctor")
    _common(doctor_parser)
    doctor_parser.add_argument("--redact-paths", action="store_true")
    doctor_parser.set_defaults(handler=doctor)
    profile_parser = commands.add_parser("profile")
    _common(profile_parser)
    profile_parser.set_defaults(handler=profile)
    serve_parser = commands.add_parser("serve")
    _common(serve_parser)
    serve_parser.add_argument("--host", default="127.0.0.1")
    serve_parser.add_argument("--port", type=int, default=8000)
    serve_parser.add_argument("--dry-run", action="store_true")
    serve_parser.set_defaults(handler=serve)
    cluster_parser = commands.add_parser("cluster")
    cluster_commands = cluster_parser.add_subparsers(dest="cluster_command", required=True)
    roundtrip = cluster_commands.add_parser("roundtrip")
    _common(roundtrip)
    roundtrip.add_argument("--run-id", required=True)
    roundtrip.add_argument("--dry-run", action="store_true")
    roundtrip.set_defaults(handler=cluster_roundtrip)
    configure = cluster_commands.add_parser("configure")
    _common(configure)
    configure.add_argument("--user", required=True)
    configure.add_argument("--host", required=True)
    configure.add_argument("--partition", required=True)
    configure.add_argument("--account", default="")
    configure.add_argument("--remote-root", default="~/ExonDomainCompare/runs")
    configure.add_argument("--tool-root", default="~/.local/share/ExonDomainCompare")
    configure.add_argument("--interproscan-launcher")
    configure.add_argument("--interproscan-module", default="")
    configure.add_argument("--interproscan-environment", default="")
    configure.add_argument("--pytmhmm-launcher")
    configure.add_argument("--pytmhmm-module", default="")
    configure.add_argument("--pytmhmm-environment", default="")
    configure.add_argument("--pytmhmm-python", default="")
    configure.add_argument("--replace-config", action="store_true")
    configure.set_defaults(handler=cluster_configure)
    cluster_doctor_parser = cluster_commands.add_parser("doctor")
    _common(cluster_doctor_parser)
    cluster_doctor_parser.add_argument("--redact-paths", action="store_true")
    cluster_doctor_parser.set_defaults(handler=cluster_doctor)
    tools_parser = cluster_commands.add_parser("tools")
    tools_commands = tools_parser.add_subparsers(dest="tools_command", required=True)
    install = tools_commands.add_parser("install")
    _common(install)
    install.add_argument(
        "--tool", choices=["all", "interproscan", "pytmhmm"], default="all")
    install.add_argument("--confirm", action="store_true")
    install.set_defaults(handler=cluster_tools_install)
    init = cluster_commands.add_parser("init")
    _common(init)
    init.add_argument("--run-id", required=True)
    init.add_argument("--dry-run", action="store_true")
    init.set_defaults(handler=cluster_init)
    runs_parser = commands.add_parser("runs")
    runs_commands = runs_parser.add_subparsers(dest="runs_command", required=True)
    runs_list_parser = runs_commands.add_parser("list")
    _common(runs_list_parser)
    runs_list_parser.set_defaults(handler=runs_list)
    migrate_parser = commands.add_parser("migrate")
    migrate_commands = migrate_parser.add_subparsers(
        dest="migrate_command", required=True)
    legacy = migrate_commands.add_parser("legacy-runs")
    _common(legacy)
    legacy.add_argument("--from", dest="source_root", required=True)
    legacy.add_argument("--mode", choices=["register", "copy", "move"], required=True)
    legacy.add_argument("--run-id", dest="run_ids", action="append")
    legacy.add_argument("--destination")
    legacy.add_argument("--dry-run", action="store_true")
    legacy.add_argument("--allow-registered-writes", action="store_true")
    legacy.add_argument("--copy-journal")
    legacy.add_argument("--confirm-move", action="store_true")
    legacy.add_argument("--rollback-journal")
    legacy.set_defaults(handler=migrate_legacy_runs)
    return parser


def main(argv: list[str] | None = None) -> int:
    try:
        args = build_parser().parse_args(argv)
        return int(args.handler(args))
    except (ConfigurationError, MigrationError) as exc:
        print(f"Configuration or migration error: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("Cancelled; no credentials were stored.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
