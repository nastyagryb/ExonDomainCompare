#!/usr/bin/env python3
from __future__ import annotations

from typing import List

from exondomaincompare.config import RuntimeConfig, discover_repository_root, load_config

_CONFIG = load_config(repository_root=discover_repository_root(__file__))


def configure(config: RuntimeConfig) -> None:
    global _CONFIG
    _CONFIG = config


def ssh_target() -> str:
    return _CONFIG.ssh_target

# Substrings that identify LRZ MFA-banner / SLURM-noise lines (case-insensitive).
_NOISE_MARKERS = (
    "!!! mfa",
    "second factor for authentication",
    "doku.lrz.de",
    "note for using push tokens",
    "without pin directly press enter",
    "2fa prompt",
    "push token",
)


def clean_ssh_output(raw: str) -> str:
    if not raw:
        return ""
    kept: List[str] = []
    for line in raw.splitlines():
        s = line.strip()
        if not s:
            continue
        low = s.lower()
        if any(marker in low for marker in _NOISE_MARKERS):
            continue
        kept.append(s)
    return "\n".join(kept).strip()




def ssh_cmd(remote_command: str) -> List[str]:
    return _CONFIG.ssh_argv(remote_command)


def scp_cmd(extra: List[str]) -> List[str]:
    return _CONFIG.scp_argv(extra)
