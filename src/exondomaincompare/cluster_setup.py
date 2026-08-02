"""Safe LRZ configuration, remote preflight and managed tool provisioning."""
from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
from pathlib import Path
from typing import Any, Mapping

from exondomaincompare.config import (
    ConfigurationError,
    RuntimeConfig,
    remote_shell_path,
    user_config_path,
)


INTERPROSCAN_VERSION = "5.78-109.0"
PYTMHMM_VERSION = "1.3.6"
DEFAULT_REMOTE_ROOT = "~/ExonDomainCompare/runs"
DEFAULT_TOOL_ROOT = "~/.local/share/ExonDomainCompare"
MANAGED_MARKER = "# Managed by ExonDomainCompare setup."
INTERPROSCAN_ARCHIVE = f"interproscan-{INTERPROSCAN_VERSION}-64-bit.tar.gz"
INTERPROSCAN_URL = (
    "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/"
    f"{INTERPROSCAN_VERSION}/{INTERPROSCAN_ARCHIVE}"
)


def _toml_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, list):
        return "[" + ", ".join(_toml_value(item) for item in value) + "]"
    return json.dumps(str(value), ensure_ascii=False)


def _table_name(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_-]+", value):
        raise ConfigurationError(f"Unsafe profile name: {value!r}")
    return value


def _section(lines: list[str], name: str, values: Mapping[str, Any]) -> None:
    lines.extend(("", f"[{name}]"))
    for key, value in values.items():
        if isinstance(value, Mapping) or value is None:
            continue
        lines.append(f"{key} = {_toml_value(value)}")


def render_managed_config(
        cfg: RuntimeConfig, *, user: str, host: str, partition: str,
        account: str = "", remote_root: str = DEFAULT_REMOTE_ROOT,
        tool_root: str = DEFAULT_TOOL_ROOT,
        interproscan_launcher: str | None = None,
        interproscan_module: str = "", interproscan_environment: str = "",
        pytmhmm_launcher: str | None = None, pytmhmm_module: str = "",
        pytmhmm_environment: str = "", pytmhmm_python: str = "") -> str:
    local_name = _table_name(cfg.local_profile_name)
    lrz_name = _table_name(cfg.lrz_profile_name)
    remote_root = str(remote_root).strip().rstrip("/")
    tool_root = str(tool_root).strip().rstrip("/")
    if not user.strip() or not host.strip() or not partition.strip():
        raise ConfigurationError("LRZ user, host and partition are required.")
    if not remote_root or not tool_root:
        raise ConfigurationError("Remote run and tool roots must not be empty.")

    lrz = dict(cfg.lrz)
    lrz.update({
        "user": user.strip(),
        "host": host.strip(),
        "remote_root": remote_root,
        "remote_work_dir": "",
        "remote_temp_dir": "",
        "partition": partition.strip(),
        "account": account.strip(),
        "managed_tool_root": tool_root,
    })
    interpro = dict(lrz.pop("interproscan", {}) or {})
    pytmhmm = dict(lrz.pop("pytmhmm", {}) or {})
    if interproscan_launcher is None:
        interpro.update({
            "launcher": f"{tool_root}/tools/interproscan-{INTERPROSCAN_VERSION}/interproscan.sh",
            "module": "",
            "environment": f"{tool_root}/envs/interproscan",
        })
    else:
        interpro.update({
            "launcher": interproscan_launcher.strip(),
            "module": interproscan_module.strip(),
            "environment": interproscan_environment.strip(),
        })
    if pytmhmm_launcher is None:
        pytmhmm.update({
            "launcher": f"{tool_root}/envs/pytmhmm/bin/pyTMHMM",
            "module": "",
            "environment": f"{tool_root}/envs/pytmhmm",
            "python": f"{tool_root}/envs/pytmhmm/bin/python",
        })
    else:
        pytmhmm.update({
            "launcher": pytmhmm_launcher.strip(),
            "module": pytmhmm_module.strip(),
            "environment": pytmhmm_environment.strip(),
            "python": pytmhmm_python.strip() or "python",
        })

    lines = [
        MANAGED_MARKER,
        'schema_version = "1.0"',
        f"active_local_profile = {_toml_value(local_name)}",
        f"active_lrz_profile = {_toml_value(lrz_name)}",
    ]
    _section(lines, f"local_profiles.{local_name}", cfg.local)
    _section(lines, f"lrz_profiles.{lrz_name}", lrz)
    _section(lines, f"lrz_profiles.{lrz_name}.interproscan", interpro)
    _section(lines, f"lrz_profiles.{lrz_name}.pytmhmm", pytmhmm)
    return "\n".join(lines) + "\n"


def write_managed_config(
        cfg: RuntimeConfig, *, output: Path | None = None, replace: bool = False,
        **values: str) -> Path:
    target = Path(output or user_config_path()).expanduser()
    if target.exists() and MANAGED_MARKER not in target.read_text(
            encoding="utf-8", errors="replace") and not replace:
        raise ConfigurationError(
            f"Refusing to replace a manually maintained configuration: {target}. "
            "Use --replace-config only after reviewing it."
        )
    text = render_managed_config(cfg, **values)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.chmod(temporary, 0o600)
    temporary.replace(target)
    return target


def initial_user_config() -> str:
    return "\n".join((
        MANAGED_MARKER,
        'schema_version = "1.0"',
        'active_local_profile = "default"',
        'active_lrz_profile = "lrz"',
        "",
    ))


def _preflight_script(cfg: RuntimeConfig) -> str:
    lrz = cfg.lrz
    interpro = lrz.get("interproscan", {}) or {}
    pytmhmm = lrz.get("pytmhmm", {}) or {}
    values = {
        "REMOTE_ROOT": str(lrz.get("remote_root") or DEFAULT_REMOTE_ROOT),
        "TOOL_ROOT": str(lrz.get("managed_tool_root") or DEFAULT_TOOL_ROOT),
        "PARTITION": str(lrz.get("partition") or ""),
        "INTERPRO": str(interpro.get("launcher") or ""),
        "INTERPRO_MODULE": str(interpro.get("module") or ""),
        "INTERPRO_ENV": str(interpro.get("environment") or ""),
        "PYTMHMM": str(pytmhmm.get("launcher") or ""),
        "PYTMHMM_MODULE": str(pytmhmm.get("module") or ""),
        "PYTMHMM_ENV": str(pytmhmm.get("environment") or ""),
        "PYTMHMM_PYTHON": str(pytmhmm.get("python") or "python"),
    }
    assignments = "\n".join(
        f"{key}={remote_shell_path(value)}" for key, value in values.items())
    return f"""set -u
{assignments}
emit() {{ printf 'EDC:%s=%s\\n' "$1" "$2"; }}
yesno() {{ if "$@" >/dev/null 2>&1; then printf yes; else printf no; fi; }}
command_ok() {{
  candidate="$1"
  if [ -z "$candidate" ]; then return 1; fi
  case "$candidate" in */*) [ -x "$candidate" ];; *) command -v "$candidate" >/dev/null 2>&1;; esac
}}
load_configured_module() {{
  requested="$1"
  [ -z "$requested" ] && return 0
  if ! command -v module >/dev/null 2>&1 && [ -r /etc/profile.d/modules.sh ]; then
    . /etc/profile.d/modules.sh
  fi
  command -v module >/dev/null 2>&1 || return 1
  module load "$requested" >/dev/null 2>&1
}}
nearest_existing() {{
  candidate="$1"
  while [ ! -e "$candidate" ] && [ "$candidate" != "/" ]; do candidate=$(dirname "$candidate"); done
  printf '%s' "$candidate"
}}
interpro_ok() {{
  load_configured_module "$INTERPRO_MODULE" || return 1
  if [ -n "$INTERPRO_ENV" ]; then
    PATH="$INTERPRO_ENV/bin:$PATH" command_ok "$INTERPRO" || return 1
    PATH="$INTERPRO_ENV/bin:$PATH" "$INTERPRO" --version 2>&1 | grep -Fq '{INTERPROSCAN_VERSION}'
  else
    command_ok "$INTERPRO" || return 1
    "$INTERPRO" --version 2>&1 | grep -Fq '{INTERPROSCAN_VERSION}'
  fi
}}
pytmhmm_ok() {{
  load_configured_module "$PYTMHMM_MODULE" || return 1
  if [ -n "$PYTMHMM_ENV" ]; then
    PATH="$PYTMHMM_ENV/bin:$PATH" command_ok "$PYTMHMM" || return 1
    PATH="$PYTMHMM_ENV/bin:$PATH" "$PYTMHMM" -h >/dev/null 2>&1 || return 1
    PATH="$PYTMHMM_ENV/bin:$PATH" "$PYTMHMM_PYTHON" -c "from importlib.metadata import version; raise SystemExit(version('pyTMHMM') != '{PYTMHMM_VERSION}')" >/dev/null 2>&1
  else
    command_ok "$PYTMHMM" || return 1
    "$PYTMHMM" -h >/dev/null 2>&1 || return 1
    "$PYTMHMM_PYTHON" -c "from importlib.metadata import version; raise SystemExit(version('pyTMHMM') != '{PYTMHMM_VERSION}')" >/dev/null 2>&1
  fi
}}
existing=$(nearest_existing "$REMOTE_ROOT")
tool_existing=$(nearest_existing "$TOOL_ROOT")
emit linux "$(yesno test "$(uname -s 2>/dev/null)" = Linux)"
emit home_writable "$(yesno test -w "$HOME")"
emit remote_root_writable "$(yesno test -d "$existing" -a -w "$existing")"
emit tool_root_writable "$(yesno test -d "$tool_existing" -a -w "$tool_existing")"
emit sbatch "$(yesno command -v sbatch)"
emit squeue "$(yesno command -v squeue)"
emit sacct "$(yesno command -v sacct)"
if command -v sinfo >/dev/null 2>&1 && [ -n "$PARTITION" ]; then
  emit partition "$(sinfo -h -p "$PARTITION" 2>/dev/null | grep -q . && printf yes || printf no)"
else emit partition no; fi
emit interproscan "$(yesno interpro_ok)"
emit pytmhmm "$(yesno pytmhmm_ok)"
emit downloader "$(command -v curl >/dev/null 2>&1 || command -v wget >/dev/null 2>&1; [ $? -eq 0 ] && printf yes || printf no)"
emit archive_tools "$(yesno sh -c 'command -v tar >/dev/null && command -v md5sum >/dev/null')"
emit environment_manager "$(command -v micromamba >/dev/null 2>&1 || command -v mamba >/dev/null 2>&1 || command -v conda >/dev/null 2>&1 || command -v python3.11 >/dev/null 2>&1 || command -v python3 >/dev/null 2>&1; [ $? -eq 0 ] && printf yes || printf no)"
emit free_kb "$(df -Pk "$existing" 2>/dev/null | awk 'NR==2 {{print $4}}')"
emit tool_free_kb "$(df -Pk "$tool_existing" 2>/dev/null | awk 'NR==2 {{print $4}}')"
emit home "$HOME"
"""


def remote_preflight(cfg: RuntimeConfig, *, redact_paths: bool = False) -> dict[str, Any]:
    cfg.require_cluster_connection()
    process = subprocess.run(
        cfg.ssh_argv("bash -s"), input=_preflight_script(cfg), text=True,
        stdout=subprocess.PIPE, stderr=None, check=False,
    )
    if process.returncode != 0:
        raise ConfigurationError(
            f"Remote preflight failed with SSH exit code {process.returncode}."
        )
    raw: dict[str, str] = {}
    for line in process.stdout.splitlines():
        if line.startswith("EDC:") and "=" in line:
            key, value = line[4:].split("=", 1)
            raw[key] = value
    required = (
        "linux", "home_writable", "remote_root_writable", "sbatch", "squeue",
        "partition", "interproscan", "pytmhmm",
    )
    checks = {key: raw.get(key) == "yes" for key in required}
    missing = [key for key, available in checks.items() if not available]
    missing_tools = [key for key in ("interproscan", "pytmhmm") if key in missing]
    missing_infrastructure = [key for key in missing if key not in missing_tools]
    if not missing:
        next_action = ".venv/bin/edc cluster roundtrip --run-id <RUN_ID>"
    elif missing_infrastructure:
        next_action = (
            "Review the reported LRZ profile, destination and scheduler checks; "
            "do not install tools yet."
        )
    else:
        tool = missing_tools[0] if len(missing_tools) == 1 else "all"
        next_action = f".venv/bin/edc cluster tools install --tool {tool}"
    return {
        "schema_version": "1.0",
        "profile": cfg.lrz_profile_name,
        "ssh_connection": "available",
        "network_contacted": True,
        "job_submitted": False,
        "remote_mutation": False,
        "checks": checks,
        "supporting_capabilities": {
            key: raw.get(key) == "yes"
            for key in ("sacct", "tool_root_writable", "downloader",
                        "archive_tools", "environment_manager")
        },
        "remote_home": "<REMOTE_HOME>" if redact_paths else raw.get("home", ""),
        "free_kb": int(raw.get("free_kb") or 0),
        "tool_free_kb": int(raw.get("tool_free_kb") or 0),
        "expected_versions": {
            "interproscan": INTERPROSCAN_VERSION,
            "pytmhmm": PYTMHMM_VERSION,
        },
        "ready_for_cluster_runs": not missing,
        "missing": missing,
        "next_action": next_action,
    }


def remote_install_plan(cfg: RuntimeConfig, tool: str) -> dict[str, Any]:
    cfg.require_cluster_connection()
    return {
        "schema_version": "1.0",
        "profile": cfg.lrz_profile_name,
        "tool": tool,
        "versions": {
            "interproscan": INTERPROSCAN_VERSION,
            "pytmhmm": PYTMHMM_VERSION,
        },
        "interproscan_download": "about 5.5 GB compressed",
        "checksum_required": True,
        "idempotent": True,
        "modifies_shell_startup_files": False,
        "submits_cluster_jobs": False,
        "network_contacted": False,
        "remote_mutation": False,
        "confirmation_required": True,
    }


def _installer_script(cfg: RuntimeConfig, tool: str) -> str:
    root = str(cfg.lrz.get("managed_tool_root") or "").strip().rstrip("/")
    if not root:
        raise ConfigurationError(
            "The selected profile has no managed_tool_root; run cluster configure first."
        )
    install_interpro = tool in {"all", "interproscan"}
    install_pytmhmm = tool in {"all", "pytmhmm"}
    root_assignment = remote_shell_path(root)
    return f"""set -euo pipefail
TOOL_ROOT={root_assignment}
DOWNLOADS="$TOOL_ROOT/downloads"
TOOLS="$TOOL_ROOT/tools"
ENVS="$TOOL_ROOT/envs"
mkdir -p "$DOWNLOADS" "$TOOLS" "$ENVS"

find_manager() {{
  for candidate in micromamba mamba conda; do
    if command -v "$candidate" >/dev/null 2>&1; then command -v "$candidate"; return 0; fi
  done
  return 1
}}

ensure_conda_env() {{
  target="$1"; shift
  manager=$(find_manager) || return 1
  if [ -d "$target/conda-meta" ]; then
    "$manager" install -y -p "$target" "$@"
  elif [ -x "$target/bin/python" ]; then
    return 0
  else
    "$manager" create -y -p "$target" "$@"
  fi
}}

install_pytmhmm() {{
  envdir="$ENVS/pytmhmm"
  if [ ! -x "$envdir/bin/python" ]; then
    if ! ensure_conda_env "$envdir" python=3.11 pip; then
      python_cmd=""
      for candidate in python3.11 python3; do
        if command -v "$candidate" >/dev/null 2>&1; then python_cmd=$(command -v "$candidate"); break; fi
      done
      [ -n "$python_cmd" ] || {{ echo 'No suitable Python or Conda/Mamba installation was found.' >&2; exit 20; }}
      "$python_cmd" -m venv "$envdir"
    fi
  fi
  "$envdir/bin/python" -m pip install "pyTMHMM=={PYTMHMM_VERSION}"
  "$envdir/bin/pyTMHMM" -h >/dev/null
  "$envdir/bin/python" -c "from importlib.metadata import version; raise SystemExit(version('pyTMHMM') != '{PYTMHMM_VERSION}')"
  echo "EDC:installed_pytmhmm=yes"
}}

install_interproscan() {{
  launcher="$TOOLS/interproscan-{INTERPROSCAN_VERSION}/interproscan.sh"
  if [ ! -x "$launcher" ]; then
    free_kb=$(df -Pk "$TOOL_ROOT" | awk 'NR==2 {{print $4}}')
    [ "${{free_kb:-0}}" -ge 30000000 ] || {{ echo 'InterProScan installation requires at least 30 GB free.' >&2; exit 21; }}
    if manager=$(find_manager); then
      ensure_conda_env "$ENVS/interproscan" python=3.11 openjdk=11 perl
    else
      command -v java >/dev/null && command -v python3 >/dev/null && command -v perl >/dev/null || {{
        echo 'InterProScan needs Java 11, Python 3 and Perl, or Conda/Mamba.' >&2; exit 22;
      }}
      java_major=$(java -version 2>&1 | awk -F'[\".]' 'NR==1 {{if ($2 == "1") print $3; else print $2}}')
      [ "${{java_major:-0}}" -ge 11 ] || {{ echo 'InterProScan needs Java 11 or newer.' >&2; exit 22; }}
    fi
    cd "$DOWNLOADS"
    archive={shlex.quote(INTERPROSCAN_ARCHIVE)}
    url={shlex.quote(INTERPROSCAN_URL)}
    if [ ! -f "$archive" ]; then
      if command -v curl >/dev/null 2>&1; then curl -fL --retry 3 -C - -o "$archive" "$url"
      elif command -v wget >/dev/null 2>&1; then wget -c -O "$archive" "$url"
      else echo 'curl or wget is required.' >&2; exit 23; fi
    fi
    if [ ! -f "$archive.md5" ]; then
      if command -v curl >/dev/null 2>&1; then curl -fL --retry 3 -o "$archive.md5" "$url.md5"
      else wget -O "$archive.md5" "$url.md5"; fi
    fi
    md5sum -c "$archive.md5"
    tar -pxzf "$archive" -C "$TOOLS"
  fi
  PATH="$ENVS/interproscan/bin:$PATH" "$launcher" --version 2>&1 | grep -Fq '{INTERPROSCAN_VERSION}'
  echo "EDC:installed_interproscan=yes"
}}

{"install_interproscan" if install_interpro else ":"}
{"install_pytmhmm" if install_pytmhmm else ":"}
echo 'EDC:installation_complete=yes'
"""


def install_remote_tools(cfg: RuntimeConfig, *, tool: str, confirmed: bool) -> int:
    if tool not in {"all", "interproscan", "pytmhmm"}:
        raise ConfigurationError(f"Unsupported managed tool: {tool}")
    if not confirmed:
        raise ConfigurationError("Remote installation requires explicit --confirm.")
    return subprocess.run(
        cfg.ssh_argv("bash -s"), input=_installer_script(cfg, tool), text=True,
        check=False,
    ).returncode
