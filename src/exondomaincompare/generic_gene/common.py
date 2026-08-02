"""Shared helpers for the generic gene-analysis layer.

Gene-agnostic only. Reads a run's standardized model tables (produced by the
core gene-analysis runner) and writes the canonical generic output layer under
``runs/<run_id>/results/generic_gene_analysis/``.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List

from exondomaincompare.config import discover_repository_root

REPO_ROOT = discover_repository_root(__file__)
FREEZE_MARKERS = ("examples/FGFR2_final_pre_interpro_30species",
                  "results/final_30_until_interpro_prepare")

# Canonical generic output layer (relative to a run directory).
GENERIC_SUBDIR = Path("results/generic_gene_analysis")
CORE_SUBDIR = Path("results/core_gene_analysis")


def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: List[Dict[str, Any]], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
        fh.write("\n")


def read_json(path: Path, default: Any = None) -> Any:
    if not path.exists():
        return default
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return default


def read_fasta(path: Path) -> Dict[str, str]:
    """Minimal FASTA reader -> {header_id: sequence}. header_id = first token."""
    seqs: Dict[str, str] = {}
    if not path.exists():
        return seqs
    hdr = None
    buf: List[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if hdr is not None:
                seqs[hdr] = "".join(buf)
            hdr = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if hdr is not None:
        seqs[hdr] = "".join(buf)
    return seqs


@dataclass
class GenericContext:
    """Everything the generic builders need, resolved once per run."""
    run_dir: Path
    run_id: str
    analysis_id: str = ""
    gene_symbol: str = "gene"
    gene_display_name: str = ""
    reference_species: str = ""
    support_level: str = ""
    cluster_status: str = "pending"
    core_dir: Path = field(init=False)
    generic_dir: Path = field(init=False)
    figures_dir: Path = field(init=False)

    def __post_init__(self) -> None:
        self.core_dir = self.run_dir / CORE_SUBDIR
        self.generic_dir = self.run_dir / GENERIC_SUBDIR
        self.figures_dir = self.generic_dir / "figures"

    def core(self, name: str) -> Path:
        return self.core_dir / name

    def out(self, name: str) -> Path:
        return self.generic_dir / name

    def assert_not_freeze(self) -> None:
        rp = str(self.run_dir.resolve())
        for m in FREEZE_MARKERS:
            if m in rp:
                raise SystemExit(f"Refusing to write into the validated freeze: {self.run_dir}")


def load_context(run_id: str) -> GenericContext:
    # Honour an isolated test run root so a canary resolves its own run rather than
    # looking for it in — or worse, writing it into — the live registry.
    from exondomaincompare.config import load_config
    runs_root = load_config(repository_root=REPO_ROOT).runs_root
    run_dir = runs_root / run_id
    if not run_dir.exists():
        raise SystemExit(f"Run not found: {run_dir}")
    cfg = read_json(run_dir / "run_config.json", {}) or {}
    # gene_config.yaml carries analysis/gene metadata; read leniently.
    analysis_id = cfg.get("analysis_id") or cfg.get("analysis") or ""
    gene_symbol = cfg.get("gene_symbol") or cfg.get("gene") or "gene"
    ref_species = cfg.get("reference_species") or ""
    support = cfg.get("support_level") or ""
    display = cfg.get("gene_display_name") or ""

    # Prefer values recorded in the core gene_model_index if present.
    gmi = read_tsv(run_dir / CORE_SUBDIR / "gene_model_index.tsv")
    if gmi:
        analysis_id = analysis_id or gmi[0].get("analysis_id", "")
        gene_symbol = gmi[0].get("gene_symbol", gene_symbol) or gene_symbol

    ctx = GenericContext(
        run_dir=run_dir, run_id=run_id, analysis_id=analysis_id or run_id,
        gene_symbol=gene_symbol, gene_display_name=display, reference_species=ref_species,
        support_level=support,
    )
    ctx.cluster_status = _cluster_status(ctx)
    return ctx


def _cluster_status(ctx: GenericContext) -> str:
    """pending | complete | failed, from on-disk domain-annotation artifacts.

    Complete requires REAL InterProScan/pyTMHMM results (domain rows or output
    files), not merely the presence of empty stage sub-folders.
    """
    dom = ctx.core("domain_features.tsv")
    if dom.exists() and read_tsv(dom):
        return "complete"
    # Only actual result files under any InterProScan */output/ folder count.
    interpro_root = ctx.run_dir / "results" / "14_interproscan"
    if interpro_root.exists():
        for out_dir in interpro_root.glob("**/output"):
            if any(p.is_file() for p in out_dir.iterdir()):
                return "complete"
    st = read_json(ctx.run_dir / "status.json", {}) or {}
    if str(st.get("status", "")).endswith("failed"):
        return "failed"
    return "pending"


def display_species(sid: str) -> str:
    if not sid:
        return ""
    return sid.replace("_", " ").strip().capitalize()
