from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional

try:  # the gene → workflow decision this contract builds on
    from . import analysis_router as _ar  # type: ignore
except Exception:  # pragma: no cover - allow standalone import
    from exondomaincompare.framework import analysis_router as _ar  # type: ignore


FAMILY_FGFR2 = "fgfr2"
FAMILY_GENERIC = "generic"
SUPPORTED_FAMILIES = (FAMILY_FGFR2, FAMILY_GENERIC)

# Schema versions are part of a card's identity: a card rendered by an older
# renderer is not interchangeable with one rendered by the current renderer even
# when its output file still exists. Bump these when the contract changes shape.
GALLERY_SCHEMA_VERSION = 3
EXPLORER_SCHEMA_VERSION = 3
FIGURE_RENDERER_VERSION = 3

# The reference runs whose architecture each family must reproduce. These are
# documentation of intent, not routing inputs: no code may branch on them.
TEMPLATE_RUNS = {
    FAMILY_FGFR2: "FGFR2 IIIb/IIIc — 30 vertebrates (validated freeze)",
    FAMILY_GENERIC: "BCL2L1 Homo sapiens + Mus musculus",
}

# Gallery scope architecture, shared by both families: a dataset with two or more
# species gets a Comparative scope plus one scope per real species and opens on
# Comparative; a single-species dataset has no Comparative scope at all.
SCOPE_COMPARATIVE = "comparative"

_FGFR2_CATEGORIES = (
    "Comparative exon structure",
    "FGFR2 cassette evidence",
    "Comparative sequence analysis",
    "Comparative domain architecture",
    "FGFR2 IIIb/IIIc Boundary Consistency",
    "Comparative exon\u2013domain boundaries",
    "Comparative genomic context",
    "Framework and QC",
    "Supplements",
)

_GENERIC_CATEGORIES = (
    "Exon structure",
    "Isoform analysis",
    "Domain architecture",
    "Boundary",
    "Genomic context",
    "Comparative",
)

# Explorer components every generic run must offer, in display order. FGFR2 runs
# offer the same set plus its cassette/event views.
_GENERIC_EXPLORER_COMPONENTS = (
    "Summary",
    "Isoforms",
    "Data & Downloads",
    "Exon Map",
    "Domain Architecture",
    "Boundary",
    "MSA",
    "Synteny",
    "Exploratory Candidate Evidence",
)

_FGFR2_EXPLORER_COMPONENTS = _GENERIC_EXPLORER_COMPONENTS + (
    "Cassette Evidence",
    "IIIb/IIIc Boundary Consistency",
)


class UnknownAnalysisFamily(ValueError):
    pass


@dataclass
class ProductionContract:
    gene_symbol: str
    analysis_family: str
    gallery_schema_version: int
    explorer_schema_version: int
    figure_renderer_version: int
    template_run: str
    categories: List[str] = field(default_factory=list)
    explorer_components: List[str] = field(default_factory=list)
    routed_by: str = "canonical_gene_symbol"
    reason: str = ""

    @property
    def is_fgfr2(self) -> bool:
        return self.analysis_family == FAMILY_FGFR2

    def scopes(self, species_ids: Optional[List[str]]) -> List[str]:
        species = [s for s in (species_ids or []) if s]
        if len(species) >= 2:
            return [SCOPE_COMPARATIVE, *species]
        return list(species)

    def default_scope(self, species_ids: Optional[List[str]]) -> str:
        scopes = self.scopes(species_ids)
        return scopes[0] if scopes else ""

    def identity(self) -> Dict[str, Any]:
        return {
            "analysis_family": self.analysis_family,
            "gallery_schema_version": self.gallery_schema_version,
            "explorer_schema_version": self.explorer_schema_version,
            "figure_renderer_version": self.figure_renderer_version,
        }

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def normalize_gene_symbol(symbol: Optional[str]) -> str:
    return _ar.normalize_gene_symbol(symbol)


def analysis_family_for_gene(gene_symbol: Optional[str]) -> str:
    return FAMILY_FGFR2 if normalize_gene_symbol(gene_symbol) == "FGFR2" else FAMILY_GENERIC


def require_supported_family(family: Optional[str]) -> str:
    value = str(family or "").strip().lower()
    if value not in SUPPORTED_FAMILIES:
        raise UnknownAnalysisFamily(
            f"analysis_family {family!r} is not supported by this build "
            f"(supported: {', '.join(SUPPORTED_FAMILIES)}). Refusing to fall back to a "
            "legacy implementation."
        )
    return value


def resolve(gene_symbol: Optional[str]) -> ProductionContract:
    sym = normalize_gene_symbol(gene_symbol)
    family = analysis_family_for_gene(sym)
    fgfr2 = family == FAMILY_FGFR2
    return ProductionContract(
        gene_symbol=sym,
        analysis_family=family,
        gallery_schema_version=GALLERY_SCHEMA_VERSION,
        explorer_schema_version=EXPLORER_SCHEMA_VERSION,
        figure_renderer_version=FIGURE_RENDERER_VERSION,
        template_run=TEMPLATE_RUNS[family],
        categories=list(_FGFR2_CATEGORIES if fgfr2 else _GENERIC_CATEGORIES),
        explorer_components=list(_FGFR2_EXPLORER_COMPONENTS if fgfr2
                                 else _GENERIC_EXPLORER_COMPONENTS),
        reason=(
            f"Canonical gene symbol {sym} is the validated event gene; routed to the modern "
            "FGFR2 production architecture."
            if fgfr2 else
            f"Canonical gene symbol {sym or '(none)'} is not FGFR2; routed to the modern "
            "generic production architecture."
        ),
    )


def resolve_for_run(run_config: Optional[Dict[str, Any]]) -> ProductionContract:
    config = run_config or {}
    contract = resolve(config.get("gene_symbol"))
    recorded = config.get("analysis_family")
    if recorded and require_supported_family(recorded) != contract.analysis_family:
        raise UnknownAnalysisFamily(
            f"run records analysis_family={recorded!r} but gene "
            f"{contract.gene_symbol!r} routes to {contract.analysis_family!r}; "
            "the gene symbol is authoritative."
        )
    return contract


def stamp(run_config: Dict[str, Any]) -> Dict[str, Any]:
    contract = resolve_for_run(run_config)
    run_config.update(contract.identity())
    run_config["production_template"] = contract.template_run
    run_config["analysis_family_routed_by"] = contract.routed_by
    return run_config


def _main(argv: Optional[List[str]] = None) -> int:  # pragma: no cover - CLI helper
    import argparse
    import json

    ap = argparse.ArgumentParser(description="Resolve the production contract for a gene.")
    ap.add_argument("gene", help="Gene symbol, e.g. FGFR2 or BCL2L1.")
    ap.add_argument("--species", nargs="*", default=None)
    args = ap.parse_args(argv)
    contract = resolve(args.gene)
    out = contract.to_dict()
    if args.species is not None:
        out["scopes"] = contract.scopes(args.species)
        out["default_scope"] = contract.default_scope(args.species)
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
