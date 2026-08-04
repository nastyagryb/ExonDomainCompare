#!/usr/bin/env python3
from __future__ import annotations

import ast
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
SOURCE = ROOT / "src" / "exondomaincompare"
sys.path.insert(0, str(SCRIPTS))

from exondomaincompare.shared_gene_analysis import gene_locus_resolution as glr  # noqa: E402

REAL_RUN = ROOT / "runs" / "2026-07-29_1347_hba_panthera_leo"
FREEZE = ROOT / "results" / "final_30_until_interpro_prepare"


# --------------------------------------------------------------------------- #
# GFF3 fixtures
# --------------------------------------------------------------------------- #
def _gff(rows: list[str]) -> str:
    return "##gff-version 3\n" + "\n".join(rows) + "\n"


def _gene(seqid: str, start: int, end: int, gid: str, symbol: str, geneid: str,
          *, synonym: str = "", description: str = "", biotype: str = "protein_coding",
          ftype: str = "gene") -> str:
    attrs = [f"ID=gene-{gid}", f"Dbxref=GeneID:{geneid}", f"Name={symbol}",
             f"gene={symbol}", f"gene_biotype={biotype}"]
    if synonym:
        attrs.append(f"gene_synonym={synonym}")
    if description:
        attrs.append(f"description={description}")
    return f"{seqid}\tGnomon\t{ftype}\t{start}\t{end}\t.\t-\t.\t" + ";".join(attrs)


def _transcript(seqid: str, start: int, end: int, gid: str, tid: str, pid: str,
                product: str, symbol: str, *, n_cds: int = 3) -> list[str]:
    rows = [f"{seqid}\tGnomon\tmRNA\t{start}\t{end}\t.\t-\t.\t"
            f"ID=rna-{tid};Parent=gene-{gid};Dbxref=Genbank:{tid};Name={tid};"
            f"gbkey=mRNA;gene={symbol};product={product};transcript_id={tid}"]
    span = max(1, (end - start) // max(1, n_cds))
    phases = ["0", "1", "0", "2", "1"]
    for i in range(n_cds):
        s = start + i * span
        e = s + span - 10
        rows.append(f"{seqid}\tGnomon\texon\t{s}\t{e}\t.\t-\t.\t"
                    f"ID=exon-{tid}-{i + 1};Parent=rna-{tid};gene={symbol};"
                    f"product={product};transcript_id={tid}")
        rows.append(f"{seqid}\tGnomon\tCDS\t{s}\t{e}\t.\t-\t{phases[i % len(phases)]}\t"
                    f"ID=cds-{pid};Parent=rna-{tid};Name={pid};gbkey=CDS;"
                    f"gene={symbol};product={product};protein_id={pid}")
    return rows


@pytest.fixture()
def loc_annotation(tmp_path: Path) -> Path:
    rows = [
        _gene("NC_056694.1", 41223621, 41224454, "LOC122209634", "LOC122209634",
              "122209634"),
        *_transcript("NC_056694.1", 41223621, 41224454, "LOC122209634",
                     "XM_042921676.1", "XP_042777610.1",
                     "hemoglobin subunit alpha", "LOC122209634"),
        _gene("NC_056694.1", 41228112, 41228945, "LOC122209636", "LOC122209636",
              "122209636"),
        *_transcript("NC_056694.1", 41228112, 41228945, "LOC122209636",
                     "XM_042921681.1", "XP_042777615.1",
                     "hemoglobin subunit alpha-like", "LOC122209636"),
        _gene("NC_056694.1", 41300000, 41301000, "LOC122210362", "LOC122210362",
              "122210362"),
        *_transcript("NC_056694.1", 41300000, 41301000, "LOC122210362",
                     "XM_042922000.1", "XP_042778000.1",
                     "hemoglobin subunit zeta", "LOC122210362"),
        _gene("NC_056694.1", 41400000, 41401000, "HBM", "HBM", "122209919"),
        *_transcript("NC_056694.1", 41400000, 41401000, "HBM",
                     "XM_042922100.1", "XP_042778100.1",
                     "hemoglobin subunit mu", "HBM"),
    ]
    path = tmp_path / "loc_genomic.gff"
    path.write_text(_gff(rows), encoding="utf-8")
    return path


def _hba_record() -> glr.NcbiGeneRecord:
    return glr.NcbiGeneRecord(
        gene_id="122209636", official_symbol="LOC122209636",
        description="hemoglobin subunit alpha-like", aliases=["HBA"],
        other_designations=["hemoglobin subunit alpha-like", "alpha globin"],
        chromosome="E3", chr_accession="NC_056694.1", exon_count=3)


def _lookup(records, status: str = "ok"):
    return lambda sym, name, taxid: (list(records), status)


NO_HIT = _lookup([], "no_hit")


# --------------------------------------------------------------------------- #
# The annotation's own symbol (routes 1-3)
# --------------------------------------------------------------------------- #
def test_an_exact_symbol_in_the_annotation_resolves_without_asking_the_source(tmp_path):
    gff = tmp_path / "g.gff"
    gff.write_text(_gff([_gene("chr1", 100, 900, "TP53", "TP53", "7157"),
                         *_transcript("chr1", 100, 900, "TP53", "NM_000546.6",
                                      "NP_000537.3", "tumor protein p53", "TP53")]),
                   encoding="utf-8")

    def must_not_run(*_args):
        raise AssertionError("the source must not be consulted when the symbol matches")

    res = glr.resolve_gene_locus(gff, "TP53", gene_lookup=must_not_run)
    assert res.status == glr.RESOLVED
    assert res.locus.symbol == "TP53"
    assert res.identity.resolution_method == glr.ROUTE_EXACT_SYMBOL
    assert res.identity.as_dict()["symbol_differs_from_source"] is False


def test_a_case_normalised_symbol_resolves(tmp_path):
    gff = tmp_path / "g.gff"
    gff.write_text(_gff([_gene("chr1", 100, 900, "foxp1b", "foxp1b", "5555"),
                         *_transcript("chr1", 100, 900, "foxp1b", "NM_1.1", "NP_1.1",
                                      "forkhead box P1b", "foxp1b")]), encoding="utf-8")
    res = glr.resolve_gene_locus(gff, "FOXP1B", gene_lookup=NO_HIT)
    assert res.status == glr.RESOLVED
    assert res.identity.resolution_method == glr.ROUTE_NORMALIZED_SYMBOL
    assert res.locus.symbol == "foxp1b"


def test_an_annotation_provided_alias_resolves(tmp_path):
    gff = tmp_path / "g.gff"
    gff.write_text(_gff([_gene("chr1", 100, 900, "FGFR2", "FGFR2", "2263",
                               synonym="BEK,KGFR,CD332"),
                         *_transcript("chr1", 100, 900, "FGFR2", "NM_1.1", "NP_1.1",
                                      "fibroblast growth factor receptor 2", "FGFR2")]),
                   encoding="utf-8")

    def must_not_run(*_args):
        raise AssertionError("an annotation-local alias must not need the source")

    res = glr.resolve_gene_locus(gff, "BEK", gene_lookup=must_not_run)
    assert res.status == glr.RESOLVED
    assert res.identity.resolution_method == glr.ROUTE_ANNOTATION_ALIAS
    assert res.locus.symbol == "FGFR2"


# --------------------------------------------------------------------------- #
# NCBI Gene alias -> GeneID -> annotation (routes 4-6): the HBA case
# --------------------------------------------------------------------------- #
def test_hba_resolves_through_the_ncbi_alias_and_the_geneid_dbxref(loc_annotation):
    res = glr.resolve_gene_locus(loc_annotation, "HBA",
                                 scientific_name="Panthera leo", taxid="9689",
                                 gene_lookup=_lookup([_hba_record()]),
                                 assembly_accession="GCF_018350215.1")
    assert res.status == glr.RESOLVED
    assert res.identity.resolution_method == glr.ROUTE_NCBI_ALIAS_GENEID
    assert res.identity.resolved_gene_id == "122209636"
    assert res.identity.resolved_official_symbol == "LOC122209636"
    assert res.locus.gene_id == "gene-LOC122209636"


def test_the_alias_route_is_only_reached_after_the_annotation_routes(loc_annotation):
    res = glr.resolve_gene_locus(loc_annotation, "HBA",
                                 gene_lookup=_lookup([_hba_record()]))
    assert res.routes_attempted[:3] == [glr.ROUTE_EXACT_SYMBOL,
                                        glr.ROUTE_NORMALIZED_SYMBOL,
                                        glr.ROUTE_ANNOTATION_ALIAS]
    assert glr.ROUTE_NCBI_ALIAS_GENEID in res.routes_attempted


def test_a_returned_official_loc_symbol_maps_when_the_dbxref_is_absent(tmp_path):
    rows = [_gene("chr1", 100, 900, "LOC999", "LOC999", ""),
            *_transcript("chr1", 100, 900, "LOC999", "XM_9.1", "XP_9.1",
                         "widget-like", "LOC999")]
    gff = tmp_path / "g.gff"
    gff.write_text(_gff(rows).replace("Dbxref=GeneID:;", ""), encoding="utf-8")
    record = glr.NcbiGeneRecord(gene_id="999", official_symbol="LOC999",
                                description="widget-like", aliases=["WIDGET"])
    res = glr.resolve_gene_locus(gff, "WIDGET", gene_lookup=_lookup([record]))
    assert res.status == glr.RESOLVED
    assert res.identity.resolution_method == glr.ROUTE_NCBI_ALIAS_SYMBOL


def test_the_user_never_has_to_know_the_loc_identifier(loc_annotation):
    res = glr.resolve_gene_locus(loc_annotation, "HBA",
                                 gene_lookup=_lookup([_hba_record()]))
    assert res.resolved
    # And asking by the LOC symbol directly still works, through route 1.
    direct = glr.resolve_gene_locus(loc_annotation, "LOC122209636", gene_lookup=NO_HIT)
    assert direct.resolved
    assert direct.identity.resolution_method == glr.ROUTE_EXACT_SYMBOL
    assert direct.locus.gene_id == res.locus.gene_id


def test_a_predicted_xm_xp_only_model_is_accepted(loc_annotation):
    res = glr.resolve_gene_locus(loc_annotation, "HBA",
                                 gene_lookup=_lookup([_hba_record()]))
    assert res.locus.transcript_count == 1
    assert res.locus.protein_count == 1


# --------------------------------------------------------------------------- #
# Requested versus source identity (part 3)
# --------------------------------------------------------------------------- #
def test_the_display_symbol_stays_the_requested_one(loc_annotation):
    res = glr.resolve_gene_locus(loc_annotation, "HBA",
                                 gene_lookup=_lookup([_hba_record()]))
    d = res.identity.as_dict()
    assert d["requested_gene_symbol"] == "HBA"
    assert d["resolved_display_symbol"] == "HBA"
    assert d["resolved_official_symbol"] == "LOC122209636"
    assert d["symbol_differs_from_source"] is True
    assert d["source_description"] == "hemoglobin subunit alpha-like"


def test_the_source_symbol_is_not_concealed(loc_annotation):
    res = glr.resolve_gene_locus(loc_annotation, "HBA",
                                 gene_lookup=_lookup([_hba_record()]))
    payload = json.dumps(res.as_dict())
    assert "LOC122209636" in payload and "122209636" in payload


# --------------------------------------------------------------------------- #
# Paralog and gene-family safety (part 4)
# --------------------------------------------------------------------------- #
def test_a_description_match_alone_never_selects_a_locus(loc_annotation):
    res = glr.resolve_gene_locus(loc_annotation, "HBA", gene_lookup=NO_HIT)
    assert res.status == glr.GENE_NOT_FOUND
    assert res.locus is None


def test_the_description_route_is_declared_supporting_only():
    assert glr.ROUTE_DESCRIPTION in glr.SUPPORTING_ONLY
    assert glr.ROUTE_NCBI_ALIAS_GENEID not in glr.SUPPORTING_ONLY


def test_the_rejected_family_members_are_recorded_on_success(loc_annotation):
    res = glr.resolve_gene_locus(loc_annotation, "HBA",
                                 gene_lookup=_lookup([_hba_record()]))
    rejected = {c.locus.symbol: c for c in res.candidates if c.decision == "rejected"}
    assert "LOC122209634" in rejected, "the neighbouring alpha locus must be recorded"
    assert rejected["LOC122209634"].reason
    # Zeta and mu are a different product and are not confused with alpha.
    assert "LOC122210362" not in rejected
    assert "HBM" not in rejected


def test_two_current_records_naming_the_symbol_are_ambiguous_not_arbitrary(loc_annotation):
    second = glr.NcbiGeneRecord(gene_id="122209634", official_symbol="LOC122209634",
                                description="hemoglobin subunit alpha", aliases=["HBA"])
    res = glr.resolve_gene_locus(loc_annotation, "HBA",
                                 gene_lookup=_lookup([_hba_record(), second]))
    assert res.status == glr.AMBIGUOUS_FAMILY
    assert res.locus is None
    assert len(res.candidates) == 2
    assert "HBA" in res.message()


def test_two_annotation_loci_with_the_same_symbol_are_ambiguous(tmp_path):
    rows = [_gene("chr1", 100, 900, "HBA1a", "HBA", "1"),
            *_transcript("chr1", 100, 900, "HBA1a", "NM_1.1", "NP_1.1", "alpha", "HBA"),
            _gene("chr1", 2000, 2900, "HBA1b", "HBA", "2"),
            *_transcript("chr1", 2000, 2900, "HBA1b", "NM_2.1", "NP_2.1", "alpha", "HBA")]
    gff = tmp_path / "g.gff"
    gff.write_text(_gff(rows), encoding="utf-8")
    res = glr.resolve_gene_locus(gff, "HBA", gene_lookup=NO_HIT)
    assert res.status == glr.AMBIGUOUS_FAMILY
    assert res.locus is None


def test_a_pseudogene_family_member_is_rejected_as_a_pseudogene(tmp_path):
    rows = [_gene("chr1", 100, 900, "LOC1", "LOC1", "1"),
            *_transcript("chr1", 100, 900, "LOC1", "XM_1.1", "XP_1.1",
                         "hemoglobin subunit alpha-like", "LOC1"),
            _gene("chr1", 5000, 5900, "LOC2", "LOC2", "2",
                  biotype="pseudogene", ftype="pseudogene"),
            *_transcript("chr1", 5000, 5900, "LOC2", "XR_2.1", "", 
                         "hemoglobin subunit alpha-like", "LOC2")]
    gff = tmp_path / "g.gff"
    gff.write_text(_gff(rows), encoding="utf-8")
    record = glr.NcbiGeneRecord(gene_id="1", official_symbol="LOC1",
                                description="hemoglobin subunit alpha-like",
                                aliases=["HBA"])
    res = glr.resolve_gene_locus(gff, "HBA", gene_lookup=_lookup([record]))
    assert res.resolved and res.locus.gene_id == "gene-LOC1"
    pseudo = [c for c in res.candidates if c.locus.symbol == "LOC2"]
    assert pseudo and pseudo[0].reason == "pseudogene"


def test_a_replaced_record_is_sent_to_review_not_used(loc_annotation):
    stale = glr.NcbiGeneRecord(gene_id="111", official_symbol="LOC111",
                               description="hemoglobin subunit alpha-like",
                               aliases=["HBA"], status="replaced", current_id="122209636")
    res = glr.resolve_gene_locus(loc_annotation, "HBA", gene_lookup=_lookup([stale]))
    assert res.status == glr.REVIEW_REQUIRED
    assert res.locus is None


# --------------------------------------------------------------------------- #
# Distinguishable failure states (part 5)
# --------------------------------------------------------------------------- #
def test_gene_not_found_and_an_unmapped_alias_are_different_states(loc_annotation):
    absent = glr.resolve_gene_locus(loc_annotation, "SOMETHINGELSE",
                                    gene_lookup=NO_HIT)
    assert absent.status == glr.GENE_NOT_FOUND
    assert "No SOMETHINGELSE-related record exists" in absent.message()

    elsewhere = glr.NcbiGeneRecord(gene_id="999999", official_symbol="LOC999999",
                                   description="hemoglobin subunit alpha-like",
                                   aliases=["HBA"])
    unmapped = glr.resolve_gene_locus(loc_annotation, "HBA",
                                      gene_lookup=_lookup([elsewhere]),
                                      assembly_accession="GCF_018350215.1")
    assert unmapped.status == glr.ALIAS_MAPPING_FAILED
    assert "found in NCBI Gene as an alias" in unmapped.message()
    assert "GCF_018350215.1" in unmapped.detail
    assert absent.message() != unmapped.message()


def test_a_parser_failure_is_not_reported_as_gene_absence(tmp_path):
    broken = tmp_path / "broken.gff"
    broken.write_text("##gff-version 3\n"
                      "chr1\tGnomon\tgene\tNOT_A_NUMBER\t900\t.\t-\t.\tID=gene-X;gene=X\n",
                      encoding="utf-8")
    res = glr.resolve_gene_locus(broken, "X", gene_lookup=NO_HIT)
    assert res.status == glr.PARSER_FAILED
    assert res.status != glr.GENE_NOT_FOUND
    assert "could not be parsed" in res.message()


def test_a_missing_annotation_file_is_a_parser_failure(tmp_path):
    res = glr.resolve_gene_locus(tmp_path / "nope.gff", "X", gene_lookup=NO_HIT)
    assert res.status == glr.PARSER_FAILED


def test_an_unreachable_source_is_not_reported_as_gene_absence(loc_annotation):
    res = glr.resolve_gene_locus(loc_annotation, "HBA",
                                 gene_lookup=_lookup([], "failed:URLError"))
    assert res.status == glr.SOURCE_UNAVAILABLE
    assert res.status != glr.GENE_NOT_FOUND


def test_every_declared_status_has_a_distinct_message(loc_annotation):
    rendered = set()
    for status in glr.STATUSES:
        if status == glr.RESOLVED:
            continue
        res = glr.Resolution(status, glr.GeneIdentity("HBA"), detail="d")
        rendered.add(res.message())
    assert len(rendered) == len(glr.STATUSES) - 1


# --------------------------------------------------------------------------- #
# Annotation parsing: the coding-exon count comes from the CDS model
# --------------------------------------------------------------------------- #
def test_the_coding_exon_count_comes_from_the_cds_not_the_gene_page(tmp_path):
    from exondomaincompare.framework.run_core_gene_analysis import parse_gene_models

    rows = [_gene("chr1", 1000, 6000, "G", "G", "42")]
    rows.append("chr1\tGnomon\tmRNA\t1000\t6000\t.\t+\t.\t"
                "ID=rna-XM_1.1;Parent=gene-G;gene=G;product=widget;transcript_id=XM_1.1")
    # Five exons; the first and last are untranslated.
    for i, (s, e) in enumerate([(1000, 1100), (2000, 2100), (3000, 3100),
                                (4000, 4100), (5000, 5100)]):
        rows.append(f"chr1\tGnomon\texon\t{s}\t{e}\t.\t+\t.\t"
                    f"ID=exon-XM_1.1-{i + 1};Parent=rna-XM_1.1;gene=G")
    for phase, (s, e) in zip(["0", "2", "1"], [(2000, 2100), (3000, 3100), (4000, 4100)]):
        rows.append(f"chr1\tGnomon\tCDS\t{s}\t{e}\t.\t+\t{phase}\t"
                    f"ID=cds-XP_1.1;Parent=rna-XM_1.1;gene=G;protein_id=XP_1.1")
    gff = tmp_path / "g.gff"
    gff.write_text(_gff(rows), encoding="utf-8")

    models = parse_gene_models(gff, "G", gene_lookup=NO_HIT, allow_network=False)
    tx = models["transcripts"][0]
    assert len(tx["exons"]) == 5, "five annotated exons"
    assert len(tx["cds"]) == 3, "three coding exons, derived from the CDS model"
    assert tx["protein_id"] == "XP_1.1"
    assert [c["phase"] for c in tx["cds"]] == ["0", "2", "1"]
    # Internal coding-exon boundaries, which is what boundary views need.
    assert len(tx["cds"]) - 1 == 2


def test_transcript_exon_and_cds_parents_link_to_the_resolved_gene(loc_annotation):
    from exondomaincompare.framework.run_core_gene_analysis import parse_gene_models

    models = parse_gene_models(loc_annotation, "HBA",
                               gene_lookup=_lookup([_hba_record()]),
                               allow_network=False)
    assert models["gene"]["gene_id"] == "gene-LOC122209636"
    assert models["gene"]["strand"] == "-"
    assert len(models["transcripts"]) == 1, "only the resolved gene's transcripts"
    tx = models["transcripts"][0]
    assert tx["transcript_id"] == "XM_042921681.1"
    assert tx["protein_id"] == "XP_042777615.1"
    assert len(tx["exons"]) == 3 and len(tx["cds"]) == 3


# --------------------------------------------------------------------------- #
# The runner's adapter keeps its old contract
# --------------------------------------------------------------------------- #
def test_the_runner_adapter_reports_the_structured_status(loc_annotation):
    from exondomaincompare.framework.run_core_gene_analysis import resolve_gene_locus as adapter

    ok = adapter(loc_annotation, "HBA", gene_lookup=_lookup([_hba_record()]),
                 allow_network=False)
    assert ok["gene"]["gene_symbol"] == "LOC122209636"
    assert ok["matched_by"] == "ncbi_alias"
    assert ok["status"] == glr.RESOLVED
    assert ok["identity"]["resolved_display_symbol"] == "HBA"

    bad = adapter(loc_annotation, "NOPE", gene_lookup=NO_HIT, allow_network=False)
    assert bad["gene"] is None
    assert bad["matched_by"] == "none"
    assert bad["status"] == glr.GENE_NOT_FOUND
    assert bad["message"]


def test_the_old_catch_all_message_is_gone_from_the_runner():
    src = (SOURCE / "framework" / "run_core_gene_analysis.py").read_text(encoding="utf-8")
    assert "no locus with this symbol or synonym" not in src


# --------------------------------------------------------------------------- #
# No species- or gene-specific production branch (part 8)
# --------------------------------------------------------------------------- #
def _executable_source(path: Path) -> str:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)):
            body = node.body
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                node.body = body[1:] or [ast.Pass()]
    return ast.unparse(tree)


@pytest.mark.parametrize("module", [
    SOURCE / "shared_gene_analysis" / "gene_locus_resolution.py",
    SOURCE / "framework" / "run_core_gene_analysis.py",
])
def test_no_production_module_branches_on_the_species_or_the_gene(module):
    code = _executable_source(module).lower()
    for token in ("panthera", "122209636", "loc122209636", "gcf_018350215",
                  "hemoglobin", "xm_042921681", "xp_042777615"):
        assert token not in code, f"{module.name} hardcodes {token!r}"


def test_the_cascade_is_generic_over_the_gene_symbol(tmp_path):
    rows = [_gene("chrZ", 500, 1500, "LOC777", "LOC777", "777"),
            *_transcript("chrZ", 500, 1500, "LOC777", "XM_7.1", "XP_7.1",
                         "myoglobin-like", "LOC777")]
    gff = tmp_path / "g.gff"
    gff.write_text(_gff(rows), encoding="utf-8")
    record = glr.NcbiGeneRecord(gene_id="777", official_symbol="LOC777",
                                description="myoglobin-like", aliases=["MB"])
    res = glr.resolve_gene_locus(gff, "MB", scientific_name="Some species",
                                 gene_lookup=_lookup([record]))
    assert res.resolved
    assert res.identity.resolution_method == glr.ROUTE_NCBI_ALIAS_GENEID
    assert res.identity.resolved_official_symbol == "LOC777"


# --------------------------------------------------------------------------- #
# Retry in place (part 6 / part 9)
# --------------------------------------------------------------------------- #
def test_the_runner_offers_an_in_place_retry():
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    out = subprocess.run(
        [sys.executable, "-m", "exondomaincompare.framework.run_core_gene_analysis",
         "--help"], capture_output=True, text=True, cwd=str(ROOT), env=env)
    assert "--reuse-run-id" in out.stdout


def test_the_retry_endpoint_reuses_the_original_run_id():
    src = (ROOT / "webapp" / "backend" / "main.py").read_text(encoding="utf-8")
    start = src.index("def retry_precluster(")
    body = src[start:src.index("\n@app.", start)]
    assert "reuse_run_id=run_id" in body, "the retry must repair the run it was given"
    assert '"run_id": run_id' in body


def test_invalidation_keeps_what_the_retry_did_not_invalidate(tmp_path):
    from exondomaincompare.framework.run_core_gene_analysis import invalidate_derived_stages

    run = tmp_path / "run"
    registry = run / "results" / "01_species_registry" / "species_registry.tsv"
    cache = (run / "results" / "02_models" / "_ncbi_datasets_cache" / "x" / "genomic.gff")
    interpro = run / "results" / "14_interproscan" / "primary" / "out.tsv"
    tm = (run / "results" / "15_exon_domain_boundary_post_interpro" / "pytmhmm_primary"
          / "output" / "tm.tsv")
    derived = run / "results" / "02_models" / "gene_candidates_x.tsv"
    stage = run / "results" / "13_final_pre_interpro_closure" / "out.tsv"
    index = run / "website_indices" / "overview_index.json"
    for p in (registry, cache, interpro, tm, derived, stage, index):
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("x", encoding="utf-8")

    invalidate_derived_stages(run)
    assert registry.is_file(), "taxon resolution already succeeded; keep it"
    assert cache.is_file(), "a 450 MB annotation must not be re-downloaded"
    assert interpro.is_file(), "results fetched from the cluster are inputs, not derived"
    assert tm.is_file(), "results fetched from the cluster are inputs, not derived"
    assert not derived.exists()
    assert not stage.parent.exists()
    assert not index.parent.exists()


# --------------------------------------------------------------------------- #
# Isoform differences with one protein model (part 7)
# --------------------------------------------------------------------------- #
def test_one_protein_model_is_not_applicable_rather_than_a_failure(tmp_path):
    from exondomaincompare.shared_gene_analysis.indices.msa import build_msa_index

    wi = tmp_path / "website_indices"
    wi.mkdir(parents=True)
    (wi / "isoform_alignment_index.json").write_text(
        json.dumps({"status": "unavailable", "sequence_count": 1, "sequences": []}),
        encoding="utf-8")
    (wi / "primary_selection_index.json").write_text(
        json.dumps({"species_id": "panthera_leo"}), encoding="utf-8")
    (wi / "gene_explorer_index.json").write_text(json.dumps({"n_species": 1}),
                                                 encoding="utf-8")

    class Ctx:
        website_indices = wi
        run_dir = tmp_path

    index = build_msa_index(Ctx())
    assert index["available"] is False
    block = index["availability"]
    assert block["state"] == "not_applicable"
    assert "At least two distinct protein sequences" in block["reason"]
    assert block["reason_code"] == "single_unique_protein_sequence"
    assert "failed" not in block["reason"].lower()


def test_the_frontend_names_the_not_applicable_state():
    common = (ROOT / "webapp" / "frontend" / "src" / "pages" / "viewers"
              / "common.js").read_text(encoding="utf-8")
    assert "not_applicable" in common
    assert "Not applicable to this run" in common


# --------------------------------------------------------------------------- #
# Regression against the real repaired run
# --------------------------------------------------------------------------- #
requires_real_run = pytest.mark.skipif(not REAL_RUN.is_dir(),
                                       reason="the real HBA run is not present")


@requires_real_run
def test_the_real_run_resolved_hba_to_the_loc_locus():
    record = json.loads((REAL_RUN / "results" / "02_models"
                         / "gene_resolution_panthera_leo.json").read_text(encoding="utf-8"))
    assert record["status"] == glr.RESOLVED
    identity = record["identity"]
    assert identity["requested_gene_symbol"] == "HBA"
    assert identity["resolved_gene_id"] == "122209636"
    assert identity["resolved_official_symbol"] == "LOC122209636"
    assert identity["resolved_display_symbol"] == "HBA"
    assert identity["source_description"] == "hemoglobin subunit alpha-like"
    assert identity["resolution_method"] == glr.ROUTE_NCBI_ALIAS_GENEID
    assert record["locus"]["seqid"] == "NC_056694.1"


@requires_real_run
def test_the_real_run_recorded_the_rejected_neighbouring_alpha_locus():
    rows = (REAL_RUN / "results" / "02_models"
            / "gene_candidates_panthera_leo.tsv").read_text(encoding="utf-8").splitlines()
    header = rows[0].split("\t")
    entries = [dict(zip(header, r.split("\t"))) for r in rows[1:]]
    by_symbol = {e["symbol"]: e for e in entries}
    assert by_symbol["LOC122209636"]["decision"] == "accepted"
    assert by_symbol["LOC122209634"]["decision"] == "rejected"
    assert by_symbol["LOC122209634"]["reason"]


@requires_real_run
def test_the_real_run_produced_one_translated_primary_protein():
    faa = REAL_RUN / "results" / "core_gene_analysis" / "proteins_primary.faa"
    text = faa.read_text(encoding="utf-8")
    assert text.count(">") == 1
    assert "XP_042777615.1" in text
    sequence = "".join(l.strip() for l in text.splitlines() if not l.startswith(">"))
    assert len(sequence) == 142, "alpha globin is 142 aa"
    assert sequence.startswith("MVLS")
    assert "*" not in sequence


@requires_real_run
def test_the_real_run_keeps_its_id_gene_and_species():
    config = json.loads((REAL_RUN / "run_config.json").read_text(encoding="utf-8"))
    assert config["run_id"] == "2026-07-29_1347_hba_panthera_leo"
    assert config["gene_symbol"] == "HBA"
    assert config["species_ids"] == ["panthera_leo"]
    assert config["gene_identity"]["requested_gene_symbol"] == "HBA"
    assert (config["gene_identity_by_species"]["panthera_leo"]
            ["resolved_official_symbol"] == "LOC122209636")


@requires_real_run
def test_the_real_run_exposes_both_symbols_to_the_ui():
    overview = json.loads((REAL_RUN / "website_indices"
                           / "overview_index.json").read_text(encoding="utf-8"))
    identity = overview["gene_identity"]
    assert identity["requested_gene_symbol"] == "HBA"
    assert identity["source_symbols_by_species"]["panthera_leo"] == "LOC122209636"
    assert identity["any_symbol_differs_from_source"] is True
    assert overview["gene_symbol"] == "HBA", "the headline stays the requested symbol"


@requires_real_run
def test_the_real_run_passed_the_primary_fasta_gate():
    status = json.loads((REAL_RUN / "status.json").read_text(encoding="utf-8"))
    assert status["primary_fasta_status"] == "available"
    assert status["primary_fasta_count"] == 1
    assert status["pre_interpro_status"] == "complete"
    assert not status.get("failed_step")


@requires_real_run
def test_the_real_run_boundary_scope_follows_the_coding_exon_count():
    index = json.loads((REAL_RUN / "website_indices"
                        / "exon_domain_boundaries_index.json").read_text(encoding="utf-8"))
    assert index["scope"] == "internal_coding_exon_boundaries"
    assert index["n_boundaries"] <= 2


@requires_real_run
def test_the_real_run_shows_isoform_differences_as_not_applicable():
    index = json.loads((REAL_RUN / "website_indices"
                        / "msa_index.json").read_text(encoding="utf-8"))
    assert index["available"] is False
    assert index["availability"]["state"] == "not_applicable"


# --------------------------------------------------------------------------- #
# Existing runs and the validated freeze
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("gene", ["FGFR1", "TP53", "TPM1", "FGFR2"])
def test_a_curated_symbol_still_resolves_by_the_first_route(tmp_path, gene):
    gff = tmp_path / f"{gene}.gff"
    gff.write_text(_gff([_gene("chr1", 100, 900, gene, gene, "1"),
                         *_transcript("chr1", 100, 900, gene, "NM_1.1", "NP_1.1",
                                      f"{gene} product", gene)]), encoding="utf-8")

    def must_not_run(*_args):
        raise AssertionError(f"{gene} must resolve without a source lookup")

    res = glr.resolve_gene_locus(gff, gene, gene_lookup=must_not_run)
    assert res.status == glr.RESOLVED
    assert res.identity.resolution_method == glr.ROUTE_EXACT_SYMBOL


def test_the_cascade_never_writes_anything(loc_annotation):
    before = {p: p.stat().st_mtime_ns for p in loc_annotation.parent.rglob("*")
              if p.is_file()}
    glr.resolve_gene_locus(loc_annotation, "HBA", gene_lookup=_lookup([_hba_record()]))
    after = {p: p.stat().st_mtime_ns for p in loc_annotation.parent.rglob("*")
             if p.is_file()}
    assert before == after, "resolution is read-only"


@pytest.mark.skipif(not FREEZE.is_dir(), reason="the validated freeze is not present")
def test_the_validated_freeze_is_untouched():
    out = subprocess.run(["git", "status", "--porcelain", "--", str(FREEZE)],
                         capture_output=True, text=True, cwd=str(ROOT))
    assert out.stdout.strip() == "", f"the freeze changed:\n{out.stdout}"
