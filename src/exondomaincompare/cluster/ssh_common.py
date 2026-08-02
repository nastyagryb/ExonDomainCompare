#!/usr/bin/env python3
"""Shared SSH helpers for the cluster round-trip scripts.

Two concerns are centralized here:

1. LRZ login nodes print an interactive MFA banner (and, for a finished job,
   ``squeue`` prints ``slurm_load_jobs error: Invalid job id specified``). Because
   the caller merges stderr into stdout, that noise otherwise gets parsed as if it
   were job data. ``clean_ssh_output`` strips those known noise lines so parsers
   only see real payload.

2. Optional SSH connection multiplexing (ControlMaster). ``run_cluster_roundtrip``
   exports ``FGFR2_SSH_CONTROL_PATH``; every ssh/scp call then reuses a single
   authenticated master connection, so the user is prompted for the LRZ password /
   MFA only once. If multiplexing is unavailable, ssh falls back to a normal
   connection automatically (ControlMaster=auto), so behavior never breaks.
"""
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
    """Return only meaningful lines from ssh stdout (drop MFA banner / blanks)."""
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
    """Build an ``ssh`` argv that reuses the multiplexed master when available."""
    return _CONFIG.ssh_argv(remote_command)


def scp_cmd(extra: List[str]) -> List[str]:
    """Build an ``scp`` argv (extra = source/dest args) reusing the master."""
    return _CONFIG.scp_argv(extra)
