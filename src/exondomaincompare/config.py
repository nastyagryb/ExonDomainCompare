"""Versioned, shared configuration and capability contract."""
from __future__ import annotations

import copy
import os
import shlex
import shutil
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

try:
    from platformdirs import user_cache_dir, user_config_dir, user_data_dir, user_log_dir
except ImportError:  # fallback for installations without platformdirs
    def _fallback_user_dir(kind: str, app: str) -> str:
        home = Path.home()
        if sys.platform == "darwin":
            base = home / "Library" / (
                "Caches" if kind == "cache" else "Logs" if kind == "log"
                else "Application Support")
        elif os.name == "nt":
            base = Path(os.environ.get("LOCALAPPDATA") or home / "AppData" / "Local")
        else:
            variable = "XDG_CACHE_HOME" if kind == "cache" else (
                "XDG_CONFIG_HOME" if kind == "config" else "XDG_DATA_HOME")
            fallback = home / (".cache" if kind == "cache" else
                               ".config" if kind == "config" else ".local/share")
            base = Path(os.environ.get(variable) or fallback)
        return str(base / app)

    def user_cache_dir(app: str, _author: str) -> str:
        return _fallback_user_dir("cache", app)

    def user_config_dir(app: str, _author: str) -> str:
        return _fallback_user_dir("config", app)

    def user_data_dir(app: str, _author: str) -> str:
        return _fallback_user_dir("data", app)

    def user_log_dir(app: str, _author: str) -> str:
        return _fallback_user_dir("log", app)

SCHEMA_VERSION = "1.0"
CONFIG_ENV = "EXONDOMAIN_CONFIG"
REPO_ENV = "EXONDOMAIN_REPO_ROOT"
RUNS_ENV = "EXONDOMAIN_RUNS_ROOT"
DATA_ENV = "EDC_DATA_DIR"
LOCAL_PROFILE_ENV = "EXONDOMAIN_LOCAL_PROFILE"
LRZ_PROFILE_ENV = "EXONDOMAIN_LRZ_PROFILE"
CONTROL_PATH_ENV = "EXONDOMAIN_SSH_CONTROL_PATH"

ENV_FIELDS = {
    RUNS_ENV: ("local", "runs_root"),
    "EDC_RUNS_ROOT": ("local", "runs_root"),
    DATA_ENV: ("local", "data_dir"),
    "EXONDOMAIN_PYTHON": ("local", "python"),
    "EXONDOMAIN_NODE": ("local", "node"),
    "EXONDOMAIN_NPM": ("local", "npm"),
    "EXONDOMAIN_DATASETS": ("local", "datasets"),
    "EXONDOMAIN_MAFFT": ("local", "mafft"),
    "EXONDOMAIN_SSH": ("local", "ssh"),
    "EXONDOMAIN_SCP": ("local", "scp"),
    "EXONDOMAIN_RSYNC": ("local", "rsync"),
    "EXONDOMAIN_LRZ_USER": ("lrz", "user"),
    "EXONDOMAIN_LRZ_HOST": ("lrz", "host"),
    "EXONDOMAIN_LRZ_REMOTE_ROOT": ("lrz", "remote_root"),
    "EDC_LRZ_REMOTE_WORK_DIR": ("lrz", "remote_work_dir"),
    "EDC_LRZ_REMOTE_TEMP_DIR": ("lrz", "remote_temp_dir"),
    "EXONDOMAIN_LRZ_ACCOUNT": ("lrz", "account"),
    "EXONDOMAIN_LRZ_PARTITION": ("lrz", "partition"),
    "EXONDOMAIN_INTERPROSCAN_LAUNCHER": ("interproscan", "launcher"),
    "EXONDOMAIN_INTERPROSCAN_MODULE": ("interproscan", "module"),
    "EXONDOMAIN_INTERPROSCAN_ENV": ("interproscan", "environment"),
    "EXONDOMAIN_PYTMHMM_LAUNCHER": ("pytmhmm", "launcher"),
    "EXONDOMAIN_PYTMHMM_MODULE": ("pytmhmm", "module"),
    "EXONDOMAIN_PYTMHMM_ENV": ("pytmhmm", "environment"),
}

_FORBIDDEN_KEYS = {
    "password", "passwd", "passphrase", "private_key", "private_key_path",
    "access_token", "api_token", "secret", "client_secret",
}


class ConfigurationError(ValueError):
    """The selected configuration cannot safely drive the requested operation."""


_GENERIC_LOCAL_PYTHON_TOKENS = frozenset({"python", "python3"})


@dataclass(frozen=True)
class LocalPythonRuntime:
    """Validated interpreter selected for local application child processes."""

    configured_token: str
    selected: str
    current: str
    selection_mode: str
    matches_current: bool

    def report(self, *, redact_paths: bool = False) -> dict[str, Any]:
        configured = self.configured_token
        if redact_paths and ("/" in configured or "\\" in configured):
            configured = "<EXPLICIT_PYTHON>"
        if redact_paths:
            selected = ("<CURRENT_PYTHON>" if self.selection_mode == "current_interpreter"
                        else "<CONFIGURED_PYTHON>")
            current = "<CURRENT_PYTHON>"
        else:
            selected = self.selected
            current = self.current
        return {
            "configured": configured,
            "selection_mode": self.selection_mode,
            "selected": selected,
            "current": current,
            "matches_current": self.matches_current,
            "validated_executable": True,
            "mismatch_visible": not self.matches_current,
        }


@dataclass(frozen=True)
class AppPaths:
    """All mutable application roots, separate from repository source."""

    repository: Path
    data: Path
    config: Path
    cache: Path
    logs: Path
    temp: Path
    packages: Path
    datasets: Path
    registry: Path
    runs: Path
    legacy_runs: Path
    deleted: Path
    quarantine: Path
    migration_staging: Path

    def setup_roots(self) -> tuple[Path, ...]:
        return (
            self.data, self.config, self.cache, self.logs, self.temp,
            self.packages, self.datasets, self.registry, self.runs,
            self.deleted, self.quarantine, self.migration_staging,
        )


def discover_repository_root(anchor: Path | str | None = None,
                             env: Mapping[str, str] | None = None) -> Path:
    env = env or os.environ
    if env.get(REPO_ENV):
        candidate = Path(env[REPO_ENV]).expanduser().resolve()
        if not (candidate / "scripts").is_dir():
            raise ConfigurationError(f"{REPO_ENV} is not an ExonDomainCompare repository.")
        return candidate
    start = Path(anchor or __file__).expanduser().resolve()
    if start.is_file():
        start = start.parent
    for candidate in (start, *start.parents):
        if (candidate / "scripts").is_dir() and (candidate / "webapp").is_dir():
            return candidate
    raise ConfigurationError("Could not discover repository root from the executable location.")


def _load_toml(path: Path) -> dict[str, Any]:
    try:
        return tomllib.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise ConfigurationError(f"Configuration file does not exist: {path}") from None
    except tomllib.TOMLDecodeError as exc:
        raise ConfigurationError(f"Invalid TOML in {path}: {exc}") from None


def _merge(base: dict[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, Mapping) and isinstance(result.get(key), Mapping):
            result[key] = _merge(dict(result[key]), value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def _user_config_path(repo_root: Path, explicit: Path | str | None,
                      env: Mapping[str, str]) -> Path | None:
    if explicit:
        return Path(explicit).expanduser().resolve()
    if env.get(CONFIG_ENV):
        return Path(env[CONFIG_ENV]).expanduser().resolve()
    repository_local = repo_root / "config" / "exondomain.local.toml"
    if repository_local.is_file():
        return repository_local
    canonical = user_config_path(env)
    if canonical.is_file():
        return canonical
    legacy = Path(env.get("XDG_CONFIG_HOME") or (Path.home() / ".config")) \
        / "exondomaincompare" / "config.toml"
    return legacy if legacy.is_file() else None


def user_config_path(env: Mapping[str, str] | None = None) -> Path:
    """Return the personal configuration path used by setup and discovery."""
    env = env or os.environ
    if env.get(CONFIG_ENV):
        return Path(env[CONFIG_ENV]).expanduser()
    if env.get(DATA_ENV):
        return Path(env[DATA_ENV]).expanduser() / "config" / "config.toml"
    if env.get("XDG_CONFIG_HOME"):
        return Path(env["XDG_CONFIG_HOME"]).expanduser() \
            / "exondomaincompare" / "config.toml"
    return Path(user_config_dir("ExonDomainCompare", "ExonDomainCompare")) \
        / "config.toml"


def remote_shell_path(value: str) -> str:
    """Quote a remote path while preserving a leading home-directory token."""
    raw = str(value or "").strip()
    if raw == "~":
        return '"$HOME"'
    if raw.startswith("~/"):
        suffix = raw[2:]
        return '"$HOME"/' + shlex.quote(suffix)
    return shlex.quote(raw)


def _check_forbidden_keys(value: Any, prefix: str = "") -> None:
    if not isinstance(value, Mapping):
        return
    for key, child in value.items():
        dotted = f"{prefix}.{key}" if prefix else str(key)
        if str(key).lower() in _FORBIDDEN_KEYS:
            raise ConfigurationError(
                f"Secret-bearing key '{dotted}' is forbidden; use SSH agent/configuration."
            )
        _check_forbidden_keys(child, dotted)


def _as_list(value: Any, dotted: str) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ConfigurationError(f"{dotted} must be an array of strings.")
    return list(value)


@dataclass(frozen=True)
class RuntimeConfig:
    repository_root: Path
    runs_root: Path
    local_profile_name: str
    lrz_profile_name: str
    local: Mapping[str, Any]
    lrz: Mapping[str, Any]
    config_source: str

    @property
    def paths(self) -> AppPaths:
        def configured(name: str, default: str) -> Path:
            raw = str(self.local.get(name) or "").strip()
            return Path(raw).expanduser().resolve() if raw else Path(default).resolve()

        configured_data = str(self.local.get("data_dir") or "").strip()
        data = configured("data_dir", user_data_dir(
            "ExonDomainCompare", "ExonDomainCompare"))
        private_root = data if configured_data else None
        return AppPaths(
            repository=self.repository_root,
            data=data,
            config=configured(
                "config_dir",
                str(private_root / "config") if private_root else
                user_config_dir("ExonDomainCompare", "ExonDomainCompare"),
            ),
            cache=configured(
                "cache_dir",
                str(private_root / "cache") if private_root else
                user_cache_dir("ExonDomainCompare", "ExonDomainCompare"),
            ),
            logs=configured(
                "logs_dir",
                str(private_root / "logs") if private_root else
                user_log_dir("ExonDomainCompare", "ExonDomainCompare"),
            ),
            temp=configured("temp_dir", str(data / "tmp")),
            packages=configured("packages_dir", str(data / "packages")),
            datasets=configured("datasets_dir", str(data / "datasets")),
            registry=configured("registry_dir", str(data / "registry")),
            runs=self.runs_root,
            legacy_runs=self.repository_root / "runs",
            deleted=configured("deleted_dir", str(data / "deleted")),
            quarantine=configured("quarantine_dir", str(data / "quarantine")),
            migration_staging=configured(
                "migration_staging_dir", str(data / "tmp" / "migrations")),
        )

    @property
    def ssh_target(self) -> str:
        user = str(self.lrz.get("user", "")).strip()
        host = str(self.lrz.get("host", "")).strip()
        return f"{user}@{host}" if user and host else ""

    def executable_token(self, name: str) -> str:
        return str(self.local.get(name) or name)

    def executable(self, name: str) -> str | None:
        token = self.executable_token(name)
        if "/" in token:
            path = Path(token).expanduser()
            return str(path.resolve()) if path.is_file() and os.access(path, os.X_OK) else None
        return shutil.which(token)

    def local_python(self, *, current_executable: str | None = None) -> LocalPythonRuntime:
        """Resolve Python for local children without letting PATH replace the caller's venv."""
        current = str(Path(current_executable or sys.executable).expanduser().absolute())
        current_path = Path(current)
        if not current_path.is_file() or not os.access(current_path, os.X_OK):
            raise ConfigurationError(
                "The currently executing Python is not an executable file."
            )
        token = self.executable_token("python").strip() or "python"
        if token.lower() in _GENERIC_LOCAL_PYTHON_TOKENS:
            selected = current
            mode = "current_interpreter"
        elif "/" in token or "\\" in token:
            candidate = Path(token).expanduser()
            if not candidate.is_absolute():
                raise ConfigurationError(
                    "An explicit local Python path must be absolute; use 'python' or "
                    "'python3' to inherit the current interpreter."
                )
            if not candidate.is_file() or not os.access(candidate, os.X_OK):
                raise ConfigurationError(
                    "The explicitly configured local Python is not executable."
                )
            selected = str(candidate)
            mode = "explicit_path"
        else:
            resolved = shutil.which(token)
            if not resolved:
                raise ConfigurationError(
                    f"The configured local Python command {token!r} is unavailable."
                )
            selected = resolved
            mode = "specific_command"
        try:
            matches = os.path.samefile(selected, current)
        except OSError:
            matches = os.path.abspath(selected) == os.path.abspath(current)
        return LocalPythonRuntime(
            configured_token=token,
            selected=selected,
            current=current,
            selection_mode=mode,
            matches_current=matches,
        )

    def command(self, argv: Sequence[str]) -> str:
        return shlex.join([str(item) for item in argv])

    def remote_run_root(self, run_id: str) -> str:
        base = str(
            self.lrz.get("remote_work_dir") or self.lrz.get("remote_root") or ""
        ).rstrip("/")
        if not base:
            raise ConfigurationError(
                f"LRZ profile '{self.lrz_profile_name}' has no remote_root."
            )
        return f"{base}/{run_id}"

    def remote_directory_plan(self, run_id: str) -> dict[str, Any]:
        remote_run = self.remote_run_root(run_id)
        temp_base = str(self.lrz.get("remote_temp_dir") or "").rstrip("/")
        remote_temp = f"{temp_base}/{run_id}" if temp_base else f"{remote_run}/temp"
        directories = [
            remote_run, f"{remote_run}/input", f"{remote_run}/output",
            f"{remote_run}/slurm", remote_temp,
        ]
        return {
            "schema_version": SCHEMA_VERSION,
            "profile": self.lrz_profile_name,
            "run_id": run_id,
            "directories": [f"profile:{self.lrz_profile_name}:{index}"
                            for index, _ in enumerate(directories)],
            "command": "mkdir -p " + " ".join(
                remote_shell_path(path) for path in directories),
            "idempotent": True,
            "network_contacted": False,
        }

    def require_cluster(self) -> None:
        missing = [
            label for label, value in (
                ("user", self.lrz.get("user")),
                ("host", self.lrz.get("host")),
                ("remote_root", self.lrz.get("remote_root")),
                ("partition", self.lrz.get("partition")),
                ("InterProScan launcher", self.lrz.get("interproscan", {}).get("launcher")),
                ("pyTMHMM launcher", self.lrz.get("pytmhmm", {}).get("launcher")),
            ) if not str(value or "").strip()
        ]
        local_missing = [name for name in ("ssh", "scp") if not self.executable(name)]
        if missing or local_missing:
            parts = []
            if missing:
                parts.append("missing profile fields: " + ", ".join(missing))
            if local_missing:
                parts.append("missing local tools: " + ", ".join(local_missing))
            raise ConfigurationError("; ".join(parts))

    def require_cluster_connection(self) -> None:
        missing = [
            label for label, value in (
                ("user", self.lrz.get("user")),
                ("host", self.lrz.get("host")),
            ) if not str(value or "").strip()
        ]
        if not self.executable("ssh"):
            missing.append("local ssh executable")
        if missing:
            raise ConfigurationError(
                "cluster connection is incomplete: " + ", ".join(missing))

    def ssh_argv(self, remote_command: str) -> list[str]:
        self.require_cluster_connection()
        opts = _as_list(self.lrz.get("ssh_options"), "lrz.ssh_options")
        control_path = os.environ.get(CONTROL_PATH_ENV, "").strip()
        if control_path:
            opts += [
                "-o", "ControlMaster=auto", "-o", f"ControlPath={control_path}",
                "-o", "ControlPersist=30m", "-o", "ConnectTimeout=30",
            ]
        return [self.executable_token("ssh"), *opts, self.ssh_target, remote_command]

    def scp_argv(self, args: Sequence[str]) -> list[str]:
        self.require_cluster_connection()
        opts = _as_list(self.lrz.get("scp_options"), "lrz.scp_options")
        control_path = os.environ.get(CONTROL_PATH_ENV, "").strip()
        if control_path:
            opts += [
                "-o", "ControlMaster=auto", "-o", f"ControlPath={control_path}",
                "-o", "ControlPersist=30m", "-o", "ConnectTimeout=30",
            ]
        return [self.executable_token("scp"), *opts, *map(str, args)]

    def capabilities(self, *, redact_paths: bool = False) -> list[dict[str, Any]]:
        required_local = {"python", "node", "npm", "datasets", "mafft", "ssh", "scp"}
        runtime = self.local_python()
        python_configured = runtime.configured_token
        if redact_paths and ("/" in python_configured or "\\" in python_configured):
            python_configured = "<CONFIGURED_CAPABILITY>"
        rows = [{
            "capability": "python",
            "scope": "local",
            "required": True,
            "available": True,
            "configured": python_configured,
            "resolved": ("CURRENT_INTERPRETER" if runtime.selection_mode == "current_interpreter"
                         else f"CONFIGURED_INTERPRETER:{Path(runtime.selected).name}"),
        }]
        for name in (
            "node", "npm", "datasets", "mafft", "ssh", "scp", "rsync",
            "diamond", "blastp", "makeblastdb", "d2", "typst", "rsvg_convert",
        ):
            resolved = self.executable(name)
            configured = self.executable_token(name)
            if redact_paths and ("/" in configured or "\\" in configured):
                configured = "<CONFIGURED_CAPABILITY>"
            rows.append({
                "capability": name,
                "scope": "local",
                "required": name in required_local,
                "available": bool(resolved),
                "configured": configured,
                "resolved": f"PATH:{Path(resolved).name}" if resolved else "",
            })
        for name, configured in (
            ("slurm", bool(self.lrz.get("partition"))),
            ("interproscan", bool(self.lrz.get("interproscan", {}).get("launcher"))),
            ("pytmhmm", bool(self.lrz.get("pytmhmm", {}).get("launcher"))),
        ):
            rows.append({
                "capability": name,
                "scope": "remote",
                "required": True,
                "available": None,
                "configured": configured,
                "resolved": "not contacted (offline doctor)",
            })
        return rows

    def public_identity(self) -> dict[str, str]:
        return {
            "schema_version": SCHEMA_VERSION,
            "local_profile": self.local_profile_name,
            "lrz_profile": self.lrz_profile_name,
        }


def load_config(*, config_path: Path | str | None = None,
                repository_root: Path | str | None = None,
                local_profile: str | None = None,
                lrz_profile: str | None = None,
                explicit: Mapping[str, str] | None = None,
                env: Mapping[str, str] | None = None) -> RuntimeConfig:
    env = env or os.environ
    root = (Path(repository_root).resolve() if repository_root
            else discover_repository_root(env=env))
    default_path = root / "config" / "exondomain.default.toml"
    raw = _load_toml(default_path)
    source = str(default_path.relative_to(root))
    user_path = _user_config_path(root, config_path, env)
    if user_path:
        raw = _merge(raw, _load_toml(user_path))
        source = "explicit/user configuration"
    _check_forbidden_keys(raw)
    if str(raw.get("schema_version")) != SCHEMA_VERSION:
        raise ConfigurationError(
            f"Unsupported configuration schema_version {raw.get('schema_version')!r}; "
            f"expected {SCHEMA_VERSION!r}."
        )
    local_name = (local_profile or env.get(LOCAL_PROFILE_ENV)
                  or raw.get("active_local_profile"))
    lrz_name = (lrz_profile or env.get(LRZ_PROFILE_ENV)
                or raw.get("active_lrz_profile"))
    local_profiles = raw.get("local_profiles") or {}
    lrz_profiles = raw.get("lrz_profiles") or {}
    if local_name not in local_profiles:
        raise ConfigurationError(f"Unknown local profile: {local_name!r}.")
    if lrz_name not in lrz_profiles:
        raise ConfigurationError(f"Unknown LRZ profile: {lrz_name!r}.")
    local = copy.deepcopy(local_profiles[local_name])
    lrz = copy.deepcopy(lrz_profiles[lrz_name])
    for env_name, (section, key) in ENV_FIELDS.items():
        value = env.get(env_name)
        if value is None or value == "":
            continue
        if section == "local":
            local[key] = value
        elif section == "lrz":
            lrz[key] = value
        else:
            lrz.setdefault(section, {})[key] = value
    for dotted, value in (explicit or {}).items():
        section, _, key = dotted.partition(".")
        if section == "local" and key:
            local[key] = value
        elif section == "lrz" and key:
            lrz[key] = value
        else:
            raise ConfigurationError(f"Unsupported explicit override: {dotted}")
    runs_value = str(local.get("runs_root") or "").strip()
    if runs_value:
        runs_root = Path(runs_value).expanduser()
        if not runs_root.is_absolute():
            runs_root = root / runs_root
    else:
        data_value = str(local.get("data_dir") or "").strip()
        data_root = (
            Path(data_value).expanduser()
            if data_value else
            Path(user_data_dir("ExonDomainCompare", "ExonDomainCompare"))
        )
        runs_root = data_root / "runs"
    return RuntimeConfig(
        repository_root=root,
        runs_root=runs_root.resolve(),
        local_profile_name=str(local_name),
        lrz_profile_name=str(lrz_name),
        local=local,
        lrz=lrz,
        config_source=source,
    )
