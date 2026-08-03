"""Tests for dynamic, gene-agnostic run creation (arbitrary genes, no YAML).

These cover the core of "any valid gene symbol + species starts a generic
exploratory run without a pre-existing configs/genes/** file":

* a generic core-only gene config is synthesized in memory and serializes to
  valid, reloadable YAML;
* the gene locus is resolved from a whole-genome GFF3 by exact symbol, by
  synonym, or reported as not-found / ambiguous (all gene-agnostic);
* FGFR2 still routes to the validated workflow while every other gene routes to
  the shared exploratory workflow.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from exondomaincompare.framework import gene_config as gc  # noqa: E402
from exondomaincompare.framework import analysis_router as router  # noqa: E402
from exondomaincompare.framework.run_core_gene_analysis import resolve_gene_locus  # noqa: E402


# --------------------------------------------------------------------------- #
# generic config synthesis (no file required)
# --------------------------------------------------------------------------- #
def test_build_generic_gene_config_is_core_only_and_fileless():
    cfg = gc.build_generic_gene_config("FOXP1", generated_by="test")
    assert cfg.gene_symbol == "FOXP1"
    assert cfg.analysis_id == "FOXP1_core_only_pilot"
    assert cfg.source_path is None            # synthesized, not from a repo YAML
    assert cfg.has_event is False             # no validated event region
    assert cfg.is_core_only_pilot is True
    assert cfg.support_level == "core_only_pilot"
    assert cfg.experimental is True
    assert cfg.raw["provenance"]["generated"] is True


def test_generic_config_lowercases_to_canonical_symbol():
    cfg = gc.build_generic_gene_config("foxp1")
    assert cfg.gene_symbol == "FOXP1"


def test_empty_symbol_rejected():
    with pytest.raises(gc.GeneConfigError):
        gc.build_generic_gene_config("   ")


def test_generic_config_serializes_to_reloadable_yaml(tmp_path):
    cfg = gc.build_generic_gene_config("TP53", generated_by="test")
    text = gc.gene_config_to_yaml(cfg)
    assert text.lstrip().startswith("#")      # provenance header
    out = tmp_path / "gene_config.yaml"
    out.write_text(text, encoding="utf-8")
    reloaded = gc.load_gene_config_lenient(out)
    assert reloaded.gene_symbol == "TP53"
    assert reloaded.analysis_id == "TP53_core_only_pilot"
    assert reloaded.has_event is False


# --------------------------------------------------------------------------- #
# gene locus resolution (exact / synonym / none / ambiguous) — gene-agnostic
# --------------------------------------------------------------------------- #
def _write_gff(tmp_path: Path, gene_lines: list[str]) -> Path:
    gff = tmp_path / "genomic.gff"
    gff.write_text("##gff-version 3\n" + "\n".join(gene_lines) + "\n", encoding="utf-8")
    return gff


def _gene(seqid, start, end, name, synonym=""):
    attrs = f"ID=gene-{name};Name={name};gene={name};gene_biotype=protein_coding"
    if synonym:
        attrs += f";gene_synonym={synonym}"
    return f"{seqid}\tRefSeq\tgene\t{start}\t{end}\t.\t+\t.\t{attrs}"


def test_resolve_locus_exact_symbol(tmp_path):
    gff = _write_gff(tmp_path, [_gene("chr1", 100, 200, "TP53")])
    r = resolve_gene_locus(gff, "TP53")
    assert r["matched_by"] == "symbol"
    assert r["gene"]["gene_symbol"] == "TP53"


def test_resolve_locus_by_synonym(tmp_path):
    # zebrafish-style: the queried symbol lives only in gene_synonym.
    gff = _write_gff(tmp_path, [_gene("chr7", 10, 90, "foxp1b", synonym="foxp1,wu:fc83a06")])
    r = resolve_gene_locus(gff, "FOXP1")
    assert r["matched_by"] == "synonym"
    assert r["gene"]["gene_symbol"] == "foxp1b"


def test_resolve_locus_not_found(tmp_path):
    gff = _write_gff(tmp_path, [_gene("chr1", 100, 200, "TP53")])
    r = resolve_gene_locus(gff, "NOSUCHGENE")
    assert r["matched_by"] == "none"
    assert r["gene"] is None


def test_resolve_locus_ambiguous_synonym(tmp_path):
    gff = _write_gff(tmp_path, [
        _gene("chr1", 100, 200, "geneA", synonym="myquery"),
        _gene("chr2", 300, 400, "geneB", synonym="myquery"),
    ])
    r = resolve_gene_locus(gff, "MYQUERY")
    assert r["matched_by"] == "ambiguous"
    assert r["gene"] is None
    assert set(r["candidates"]) == {"geneA", "geneB"}


# --------------------------------------------------------------------------- #
# workflow routing: FGFR2 validated, everything else shared exploratory
# --------------------------------------------------------------------------- #
def test_fgfr2_routes_to_validated():
    wf = router.resolve_gene_workflow("FGFR2")
    assert wf.is_validated is True
    assert wf.workflow == router.WORKFLOW_VALIDATED
    assert wf.creator == "run_pre_interpro_for_run.py"


@pytest.mark.parametrize("sym", ["FOXP1", "TP53", "TPM1", "FGFR1", "SOMENEWGENE"])
def test_non_fgfr2_routes_to_shared_exploratory(sym):
    wf = router.resolve_gene_workflow(sym)
    assert wf.is_validated is False
    assert wf.workflow == router.WORKFLOW_SHARED
    assert wf.creator == "run_core_gene_analysis.py"
    assert wf.has_event is False
