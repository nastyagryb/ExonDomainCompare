"""Shared helpers for FGFR2-compatible index builders."""
from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

from exondomaincompare.config import discover_repository_root

REPO_ROOT = discover_repository_root(__file__)


def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_json(path: Path, default: Any = None) -> Any:
    if not path.is_file():
        return default
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return default


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
        fh.write("\n")


def to_int(v: Any) -> Optional[int]:
    if v is None or v == "":
        return None
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def shared_exon_group_id(genomic_start: Any, genomic_end: Any, strand: Any) -> str:
    """Genomic identity of one coding exon, shared by every transcript that carries it.

    Transcript-local exon numbers (E1, E2, …) are not identity: an alternative first
    exon renumbers everything behind it. Identity is therefore the genomic CDS interval
    plus the strand, which is invariant under upstream insertions and deletions — the
    property "Compare transcript models" relies on to avoid reporting unchanged
    downstream exons as different. Every builder must derive the id here so two paths
    over the same exon cannot disagree.
    """
    token = f"{genomic_start}:{genomic_end}:{strand or ''}".encode()
    return "SEG_" + hashlib.sha1(token).hexdigest()[:10]


def rel(run_dir: Path, path: Path) -> str:
    try:
        return str(path.relative_to(run_dir))
    except ValueError:
        return str(path)


def display_species(species_id: str) -> str:
    """``gallus_gallus`` -> ``Gallus gallus``.

    A binomial capitalises the genus only, so title case would misspell every
    species name the interface shows.
    """
    if not species_id:
        return ""
    text = species_id.replace("_", " ").strip()
    return text[:1].upper() + text[1:]


@dataclass
class SharedRunContext:
    run_dir: Path
    run_id: str
    gene_symbol: str = ""
    analysis_id: str = ""

    @property
    def core_dir(self) -> Path:
        return self.run_dir / "results" / "core_gene_analysis"

    @property
    def generic_dir(self) -> Path:
        return self.run_dir / "results" / "generic_gene_analysis"

    @property
    def website_indices(self) -> Path:
        return self.run_dir / "website_indices"

    @classmethod
    def from_run_dir(cls, run_dir: Path) -> "SharedRunContext":
        run_dir = run_dir.resolve()
        cfg = read_json(run_dir / "run_config.json", {})
        return cls(
            run_dir=run_dir,
            run_id=run_dir.name,
            gene_symbol=str(cfg.get("gene_symbol") or ""),
            analysis_id=str(cfg.get("analysis_id") or ""),
        )
