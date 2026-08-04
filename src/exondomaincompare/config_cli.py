#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys

from exondomaincompare.config import ConfigurationError, RuntimeConfig, load_config


def doctor_report(cfg: RuntimeConfig, *, redact_paths: bool = False) -> dict:
    paths = cfg.paths
    runtime = cfg.local_python()
    if runtime.matches_current:
        workbook_available = importlib.util.find_spec("openpyxl") is not None
    else:
        probe = subprocess.run(
            [runtime.selected, "-c", "import openpyxl"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=15,
            check=False,
        )
        workbook_available = probe.returncode == 0
    return {
        "schema_version": "1.0",
        "offline": True,
        "network_contacted": False,
        "repository_root": ("<REPOSITORY_ROOT>" if redact_paths
                            else str(cfg.repository_root)),
        "runs_root": "<RUNS_ROOT>" if redact_paths else str(cfg.runs_root),
        "application_roots": {
            name: (f"<{name.upper()}_ROOT>" if redact_paths else str(value))
            for name, value in (
                ("data", paths.data), ("config", paths.config),
                ("cache", paths.cache), ("logs", paths.logs),
                ("temp", paths.temp), ("packages", paths.packages),
                ("datasets", paths.datasets), ("registry", paths.registry),
                ("runs", paths.runs), ("legacy_runs", paths.legacy_runs),
            )
        },
        "local_profile": cfg.local_profile_name,
        "lrz_profile": cfg.lrz_profile_name,
        "config_source": cfg.config_source,
        "local_python_runtime": runtime.report(redact_paths=redact_paths),
        "workbook_capability": {
            "available": workbook_available,
            "checked_with_selected_interpreter": True,
            "requirement": "openpyxl",
        },
        "cluster_profile_complete": not any(
            not str(value or "").strip() for value in (
                cfg.lrz.get("user"), cfg.lrz.get("host"), cfg.lrz.get("remote_root"),
                cfg.lrz.get("partition"),
                cfg.lrz.get("interproscan", {}).get("launcher"),
                cfg.lrz.get("pytmhmm", {}).get("launcher"),
            )
        ),
        "capabilities": cfg.capabilities(redact_paths=redact_paths),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate ExonDomainCompare configuration offline.")
    parser.add_argument("command", choices=["doctor", "show-profiles"], nargs="?", default="doctor")
    parser.add_argument("--config")
    parser.add_argument("--local-profile")
    parser.add_argument("--lrz-profile")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--redact-paths", action="store_true")
    args = parser.parse_args(argv)
    try:
        cfg = load_config(
            config_path=args.config,
            local_profile=args.local_profile,
            lrz_profile=args.lrz_profile,
        )
    except ConfigurationError as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        return 2
    report = doctor_report(cfg, redact_paths=args.redact_paths)
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(f"Configuration schema: {report['schema_version']}")
        print(f"Local profile: {cfg.local_profile_name}")
        print(f"LRZ profile: {cfg.lrz_profile_name}")
        print(f"Runs root: {cfg.runs_root}")
        runtime = report["local_python_runtime"]
        print(f"Local Python selection: {runtime['selection_mode']}; "
              f"matches current: {runtime['matches_current']}")
        workbook = report["workbook_capability"]
        print(f"Workbook support: {'available' if workbook['available'] else 'missing'}")
        print("Offline doctor: no LRZ connection attempted")
        for row in report["capabilities"]:
            state = "configured (remote, not contacted)" if row["scope"] == "remote" \
                else ("available" if row["available"] else "missing")
            requirement = "required" if row["required"] else "optional"
            print(f"- {row['capability']}: {state}; {requirement}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
