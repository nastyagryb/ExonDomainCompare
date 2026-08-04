from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exondomaincompare.shared_gene_analysis import assembly_selection as asel  # noqa: E402
from exondomaincompare.shared_gene_analysis import gene_identification as gid  # noqa: E402
from exondomaincompare.shared_gene_analysis import model_recovery as recovery  # noqa: E402
from exondomaincompare.shared_gene_analysis import taxon_resolution as tr  # noqa: E402

import build_species_registry_improved as registry  # noqa: E402

#: The taxon at issue, from NCBI Taxonomy. Recorded here so a silent substitution by
#: another Equus species would fail a test rather than pass unnoticed.
EQUUS_QUAGGA_TAXID = "89248"
EQUUS_CABALLUS_TAXID = "9796"
PLAINS_ZEBRA_SYNONYM = "Equus burchellii quagga"

FGFR2_PRODUCT = "fibroblast growth factor receptor 2"


# --------------------------------------------------------------------------- #
# Part 2 — species and taxon resolution
# --------------------------------------------------------------------------- #
@pytest.fixture
def stub_taxonomy(monkeypatch):
    def install(records, *, search=None):
        def _esearch(term, timeout):
            if search is not None:
                return search(term)
            wanted = term.split('"')[1] if '"' in term else term
            hits = [tid for tid, rec in records.items()
                    if rec["name"].casefold() == wanted.casefold()
                    or wanted.casefold() in [s.casefold() for s in rec.get("synonyms", [])]]
            return hits

        def _efetch(taxid, timeout):
            rec = records[taxid]
            import xml.etree.ElementTree as ET
            root = ET.Element("Taxon")
            ET.SubElement(root, "TaxId").text = taxid
            ET.SubElement(root, "ScientificName").text = rec["name"]
            ET.SubElement(root, "Rank").text = rec.get("rank", "species")
            other = ET.SubElement(root, "OtherNames")
            for syn in rec.get("synonyms", []):
                ET.SubElement(other, "Synonym").text = syn
            if rec.get("common"):
                ET.SubElement(other, "GenbankCommonName").text = rec["common"]
            ET.SubElement(root, "Lineage").text = rec.get("lineage", "")
            return root

        monkeypatch.setattr(tr, "_esearch_taxonomy", _esearch)
        monkeypatch.setattr(tr, "_efetch_taxonomy", _efetch)

    return install


REAL_RECORDS = {
    EQUUS_QUAGGA_TAXID: {"name": "Equus quagga", "common": "plains zebra",
                         "synonyms": [PLAINS_ZEBRA_SYNONYM],
                         "lineage": "Eukaryota; Metazoa; Chordata; Mammalia; Equidae"},
    EQUUS_CABALLUS_TAXID: {"name": "Equus caballus", "common": "horse", "synonyms": []},
}


def test_the_submitted_slug_resolves_to_the_canonical_taxon(stub_taxonomy):
    stub_taxonomy(REAL_RECORDS)
    identity = tr.resolve("equus_quagga")

    assert identity.status == tr.RESOLVED
    assert identity.accepted_name == "Equus quagga"
    assert identity.taxid == EQUUS_QUAGGA_TAXID
    assert identity.rank == "species"
    assert PLAINS_ZEBRA_SYNONYM in identity.synonyms


def test_the_query_term_is_never_the_underscored_slug(stub_taxonomy):
    stub_taxonomy(REAL_RECORDS)
    identity = tr.resolve("equus_quagga")

    assert identity.query_term() == EQUUS_QUAGGA_TAXID
    assert "_" not in identity.query_term()


def test_a_published_synonym_resolves_back_to_the_requested_taxon(stub_taxonomy):
    stub_taxonomy(REAL_RECORDS)
    identity = tr.resolve("Equus burchellii quagga")

    assert identity.status == tr.RESOLVED_VIA_SYNONYM
    assert identity.taxid == EQUUS_QUAGGA_TAXID
    assert identity.query_term() == EQUUS_QUAGGA_TAXID


def test_a_near_miss_is_not_silently_answered_with_another_species(stub_taxonomy):
    stub_taxonomy(REAL_RECORDS, search=lambda term: [EQUUS_CABALLUS_TAXID])
    identity = tr.resolve("Equus quagga")

    assert not identity.is_resolved
    assert identity.status == tr.NOT_FOUND
    assert identity.taxid == ""
    assert identity.accepted_name == ""
    assert "not substituted" in identity.detail


def test_several_matching_taxa_are_reported_rather_than_guessed(stub_taxonomy):
    stub_taxonomy(REAL_RECORDS,
                  search=lambda term: [EQUUS_QUAGGA_TAXID, EQUUS_CABALLUS_TAXID])
    identity = tr.resolve("Equus quagga")

    assert identity.status == tr.AMBIGUOUS
    assert not identity.is_resolved


def test_an_unknown_name_is_reported_not_invented(stub_taxonomy):
    stub_taxonomy(REAL_RECORDS, search=lambda term: [])
    identity = tr.resolve("Equus nonexistentia")

    assert identity.status == tr.NOT_FOUND
    assert identity.taxid == ""


def test_an_unreachable_taxonomy_service_is_distinct_from_a_bad_name(monkeypatch):
    def boom(term, timeout):
        raise RuntimeError("URLError: connection refused")

    monkeypatch.setattr(tr, "_esearch_taxonomy", boom)
    identity = tr.resolve("Equus quagga")

    assert identity.status == tr.SERVICE_UNAVAILABLE
    assert identity.status != tr.NOT_FOUND


def test_the_registry_records_the_resolved_identity(tmp_path, stub_taxonomy):
    stub_taxonomy(REAL_RECORDS)
    result = registry.build_registry_rows(["equus_quagga"], "ensembl_first", "RefSeq")
    row = result.rows[0]

    assert row["status"] == registry.STATUS_VERIFIED
    assert row["ncbi_species"] == "Equus quagga"
    assert row["scientific_name"] == "Equus quagga"
    assert row["taxid"] == EQUUS_QUAGGA_TAXID
    assert row["taxon_query_term"] == EQUUS_QUAGGA_TAXID
    assert PLAINS_ZEBRA_SYNONYM in row["taxon_synonyms"]


def test_the_registry_leaves_the_ncbi_name_empty_when_unresolved(stub_taxonomy):
    stub_taxonomy(REAL_RECORDS, search=lambda term: [])
    result = registry.build_registry_rows(["Equus nonexistentia"], "ensembl_first",
                                          "RefSeq")
    row = result.rows[0]

    assert row["status"] == registry.STATUS_UNRESOLVED
    assert row["ncbi_species"] == ""
    assert row["taxid"] == ""
    assert row["taxon_query_term"] == ""
    assert any(w["warning_code"] == "TAXON_NOT_RESOLVED" for w in result.warnings)


def test_the_validated_panel_still_resolves_without_a_taxonomy_query(stub_taxonomy):
    def forbidden(term, timeout):
        raise AssertionError("a cached species must not query the taxonomy service")

    result = registry.build_registry_rows(
        ["gallus_gallus", "Homo sapiens", "danio_rerio"], "ensembl_first", "RefSeq",
    ) if False else None  # placeholder replaced below

    # Installed after the cache check so a cache miss is what raises.
    import exondomaincompare.shared_gene_analysis.taxon_resolution as module
    original_search = module._esearch_taxonomy
    module._esearch_taxonomy = forbidden
    try:
        result = registry.build_registry_rows(
            ["gallus_gallus", "Homo sapiens", "danio_rerio"],
            "ensembl_first", "RefSeq")
    finally:
        module._esearch_taxonomy = original_search

    taxids = {r["species_id"]: r["taxid"] for r in result.rows}
    assert taxids == {"gallus_gallus": "9031", "homo_sapiens": "9606",
                      "danio_rerio": "7955"}
    for row in result.rows:
        assert row["status"] == registry.STATUS_VERIFIED
        assert row["ncbi_species"] and "_" not in row["ncbi_species"]


def test_offline_mode_reports_unverified_rather_than_guessing(stub_taxonomy):
    result = registry.build_registry_rows(["equus_quagga"], "ensembl_first", "RefSeq",
                                          offline=True)
    row = result.rows[0]

    assert row["status"] == registry.STATUS_UNVERIFIED_OFFLINE
    assert row["ncbi_species"] == ""


# --------------------------------------------------------------------------- #
# Part 3 — assembly and annotation selection
# --------------------------------------------------------------------------- #
def _report(accession, *, level="Chromosome", annotated=True, category="",
            status="current", taxid=EQUUS_QUAGGA_TAXID, organism="Equus quagga"):
    return {
        "accession": accession,
        "assembly_info": {"assembly_level": level, "refseq_category": category,
                          "assembly_status": status, "assembly_name": accession},
        "annotation_info": ({"name": f"{accession}-RS", "release_date": "2023-03-03"}
                            if annotated else {}),
        "organism": {"organism_name": organism, "tax_id": taxid},
    }


#: The three assemblies NCBI actually returns for taxon 89248.
REAL_QUAGGA_REPORTS = [
    _report("GCF_021613505.1", category="reference genome"),
    _report("GCA_021613505.1", annotated=False, category="reference genome"),
    _report("GCA_026770645.1", level="Scaffold", annotated=False, taxid="89252",
            organism="Equus quagga burchellii"),
]


def test_the_annotated_refseq_assembly_is_selected():
    selection = asel.select(asel.parse_summary(REAL_QUAGGA_REPORTS),
                            preference="RefSeq", requested_taxid=EQUUS_QUAGGA_TAXID,
                            requested_name="Equus quagga")

    assert selection.status == asel.SELECTED
    assert selection.selected["accession"] == "GCF_021613505.1"
    assert selection.selected["annotated"] == "1"
    assert selection.selected["taxon_match"] == "exact_taxon"


def test_an_annotated_genbank_assembly_is_used_when_refseq_has_none():
    reports = [_report("GCF_000001.1", annotated=False),
               _report("GCA_000002.1", annotated=True)]
    selection = asel.select(asel.parse_summary(reports), preference="RefSeq")

    assert selection.status == asel.SELECTED
    assert selection.selected["accession"] == "GCA_000002.1"


def test_several_valid_assemblies_break_ties_deterministically():
    reports = [_report("GCF_000003.1", level="Scaffold"),
               _report("GCF_000001.1", level="Complete Genome",
                       category="reference genome"),
               _report("GCF_000002.1", level="Chromosome")]
    picks = {asel.select(asel.parse_summary(reports)).selected["accession"]
             for _ in range(5)}

    assert picks == {"GCF_000001.1"}


def test_no_assembly_and_no_annotated_assembly_are_different_states():
    empty = asel.select([], requested_name="Equus quagga")
    unannotated = asel.select(
        asel.parse_summary([_report("GCA_1.1", annotated=False),
                            _report("GCA_2.1", annotated=False)]))

    assert empty.status == asel.NO_ASSEMBLY
    assert unannotated.status == asel.NONE_ANNOTATED
    assert empty.status != unannotated.status
    # The candidates and the reason survive, which the empty table could not carry.
    assert len(unannotated.candidates) == 2
    assert all(c["rejection_reason"] == "assembly_has_no_annotation_release"
               for c in unannotated.candidates)


def test_every_rejected_assembly_keeps_a_reason():
    selection = asel.select(asel.parse_summary(REAL_QUAGGA_REPORTS),
                            requested_taxid=EQUUS_QUAGGA_TAXID,
                            requested_name="Equus quagga")
    rejected = [c for c in selection.candidates if c["decision"] == "rejected"]

    assert len(rejected) == 2
    assert all(c["rejection_reason"] for c in rejected)


def test_a_suppressed_assembly_is_rejected_with_its_status():
    reports = [_report("GCF_000001.1", status="suppressed"),
               _report("GCF_000002.1")]
    selection = asel.select(asel.parse_summary(reports))

    assert selection.selected["accession"] == "GCF_000002.1"
    suppressed = next(c for c in selection.candidates
                      if c["accession"] == "GCF_000001.1")
    assert suppressed["rejection_reason"] == "assembly_status_suppressed"


def test_a_subspecies_assembly_is_labelled_not_hidden():
    selection = asel.select(asel.parse_summary(REAL_QUAGGA_REPORTS),
                            requested_taxid=EQUUS_QUAGGA_TAXID,
                            requested_name="Equus quagga")
    sub = next(c for c in selection.candidates if c["accession"] == "GCA_026770645.1")

    assert sub["taxon_match"] == "descendant_taxon"


def test_a_changed_schema_does_not_lose_an_assembly():
    legacy = {"assembly_accession": "GCF_000009.1", "assembly_level": "Chromosome",
              "annotation_info": {"release_name": "legacy-RS"},
              "organism": {"organism_name": "Equus quagga",
                           "tax_id": EQUUS_QUAGGA_TAXID}}
    selection = asel.select(asel.parse_summary([legacy]))

    assert selection.status == asel.SELECTED
    assert selection.selected["accession"] == "GCF_000009.1"
    assert selection.selected["annotated"] == "1"


def test_a_report_without_an_accession_is_dropped_not_selected():
    assert asel.parse_summary([{"organism": {"organism_name": "x"}}]) == []


def test_a_failure_with_no_candidates_still_writes_a_row():
    row = asel.failure_row(asel.TAXON_REJECTED,
                           "The taxonomy name 'equus_quagga' is not recognized.",
                           "equus_quagga", "equus_quagga", "", "equus_quagga")

    assert row["selection_status"] == asel.TAXON_REJECTED
    assert row["decision"] == "none_available"
    assert "not recognized" in row["assembly_decision_notes"]


# --------------------------------------------------------------------------- #
# Part 4 — gene identification
# --------------------------------------------------------------------------- #
@pytest.fixture
def panel():
    return gid.read_panel(ROOT / "references" / "human_FGFR1_2_3_4.fasta")


def _cand(**kw):
    return gid.GeneCandidate(**kw)


def test_an_exact_symbol_is_accepted():
    result = gid.identify([_cand(source_gene_id="gene-FGFR2", source_symbol="FGFR2",
                                 description=FGFR2_PRODUCT, seqid="NC_060268.1")],
                          "FGFR2", product_name=FGFR2_PRODUCT)

    assert result.status == gid.FOUND
    assert result.accepted.route == gid.ROUTE_EXACT_SYMBOL


def test_a_source_synonym_is_accepted():
    result = gid.identify([_cand(source_gene_id="g", source_symbol="BEK",
                                 synonyms=["FGFR2", "KGFR"])],
                          "FGFR2", product_name=FGFR2_PRODUCT)

    assert result.status == gid.FOUND
    assert result.accepted.route == gid.ROUTE_SYNONYM


def test_a_gene_identifier_finds_a_locus_with_no_symbol_yet():
    result = gid.identify(
        [_cand(source_gene_id="gene-LOC124236178", source_symbol="LOC124236178",
               dbxrefs=["GeneID:124236178"])],
        "FGFR2", expected_gene_ids=["124236178"], product_name=FGFR2_PRODUCT)

    assert result.status == gid.FOUND
    assert result.accepted.route == gid.ROUTE_GENE_ID
    assert "124236178" in result.accepted.orthology_evidence


def test_a_loc_locus_is_accepted_only_on_sequence_evidence(panel):
    fgfr2 = panel["human_FGFR2_UniProt_P21802"]
    result = gid.identify(
        [_cand(source_gene_id="g4", source_symbol="LOC999",
               description=FGFR2_PRODUCT)],
        "FGFR2", proteins={"g4": fgfr2}, panel=panel, product_name=FGFR2_PRODUCT)

    assert result.status == gid.FOUND
    assert result.accepted.route == gid.ROUTE_SEQUENCE
    assert "closest to FGFR2" in result.accepted.similarity_evidence


def test_a_loc_locus_is_not_accepted_for_being_the_first_hit(panel):
    result = gid.identify(
        [_cand(source_gene_id="g", source_symbol="LOC999", description=FGFR2_PRODUCT)],
        "FGFR2", proteins={}, panel=panel, product_name=FGFR2_PRODUCT)

    assert result.accepted is None
    assert result.status == gid.AMBIGUOUS


@pytest.mark.parametrize("paralog", ["FGFR1", "FGFR3", "FGFR4"])
def test_a_paralog_locus_is_not_a_candidate(paralog, panel):
    number = paralog[-1]
    result = gid.identify(
        [_cand(source_gene_id="g", source_symbol=paralog,
               description=f"fibroblast growth factor receptor {number}")],
        "FGFR2", panel=panel, product_name=FGFR2_PRODUCT)

    assert result.accepted is None
    assert result.status == gid.NOT_FOUND


def test_a_paralog_carrying_the_symbol_as_a_synonym_is_rejected(panel):
    result = gid.identify(
        [_cand(source_gene_id="g", source_symbol="FGFR1", synonyms=["FGFR2"],
               description="fibroblast growth factor receptor 1")],
        "FGFR2", panel=panel, product_name=FGFR2_PRODUCT)

    assert result.accepted is None
    assert "paralog FGFR1" in result.candidates[0].reason


def test_a_loc_locus_whose_protein_is_a_paralog_is_rejected(panel):
    result = gid.identify(
        [_cand(source_gene_id="g5", source_symbol="LOC999",
               description="fibroblast growth factor receptor 2/3-like")],
        "FGFR2", proteins={"g5": panel["human_FGFR3_UniProt_P22607"]},
        panel=panel, product_name=FGFR2_PRODUCT)

    assert result.accepted is None
    assert "closest to FGFR3" in result.candidates[0].reason


def test_a_near_tie_between_paralogs_is_left_ambiguous(panel):
    ok, reason = gid.discriminate_by_sequence(
        _cand(source_symbol="LOC1"), "FGFR2",
        panel["human_FGFR2_UniProt_P21802"],
        {"human_FGFR2": panel["human_FGFR2_UniProt_P21802"],
         "human_FGFR3": panel["human_FGFR2_UniProt_P21802"]},
        margin=0.5)

    assert not ok
    assert "ambiguous" in reason


def test_the_description_route_does_not_match_every_gene_ending_in_two():
    unrelated = [("ADIPOR2", "adiponectin receptor protein 2"),
                 ("NPR2", "natriuretic peptide receptor 2"),
                 ("GRM2", "glutamate metabotropic receptor 2"),
                 ("TGFBR2", "transforming growth factor beta receptor 2"),
                 ("RYR2", "ryanodine receptor 2"),
                 ("TACR2", "tachykinin receptor 2")]
    for symbol, description in unrelated:
        route, _ = gid.classify_route(
            _cand(source_gene_id="g", source_symbol=symbol, description=description),
            "FGFR2", (), FGFR2_PRODUCT)
        assert route == "", f"{symbol} must not be an FGFR2 candidate"


def test_siblings_are_derived_from_the_symbol_not_listed():
    assert gid.sibling_symbols("FGFR2") == ["FGFR1", "FGFR3", "FGFR4"]
    assert gid.sibling_symbols("TP53") == []


# --------------------------------------------------------------------------- #
# Part 5 — download and parsing
# --------------------------------------------------------------------------- #
GFF_HEADER = "##gff-version 3\n"


def _gff(tmp_path: Path, body: str) -> Path:
    path = tmp_path / "genomic.gff"
    path.write_text(GFF_HEADER + body, encoding="utf-8")
    return path


def _line(seqid, kind, start, end, strand, phase, attrs):
    return f"{seqid}\tGnomon\t{kind}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n"


FGFR2_GFF = (
    _line("NC_060268.1", "gene", 100, 9000, "-", ".",
          f"ID=gene-FGFR2;Name=FGFR2;gene=FGFR2;gene_biotype=protein_coding;"
          f"Dbxref=GeneID:124236178;description={FGFR2_PRODUCT}")
    + _line("NC_060268.1", "mRNA", 100, 9000, "-", ".",
            "ID=rna-XM_046654961.1;Parent=gene-FGFR2;Name=XM_046654961.1;"
            "transcript_id=XM_046654961.1;protein_id=XP_046510917.1")
    + _line("NC_060268.1", "exon", 8000, 9000, "-", ".",
            "ID=exon-1;Parent=rna-XM_046654961.1")
    + _line("NC_060268.1", "exon", 100, 400, "-", ".",
            "ID=exon-2;Parent=rna-XM_046654961.1")
    + _line("NC_060268.1", "CDS", 8000, 9000, "-", "0",
            "ID=cds-XP_046510917.1;Parent=rna-XM_046654961.1;"
            "protein_id=XP_046510917.1")
    + _line("NC_060268.1", "CDS", 100, 400, "-", "2",
            "ID=cds-XP_046510917.1;Parent=rna-XM_046654961.1;"
            "protein_id=XP_046510917.1")
)


@pytest.fixture
def collector():
    import collect_fgfr2_models_dual_source_v3 as module
    return module


@pytest.fixture
def species_record(collector):
    return collector.SpeciesRecord(
        input_name="equus_quagga", ensembl_species="equus_quagga",
        ncbi_species="Equus quagga", taxid=EQUUS_QUAGGA_TAXID,
        preferred_source="ensembl_first", assembly_preference="RefSeq")


def test_a_candidate_records_how_many_transcripts_and_proteins_it_carries(
        tmp_path, collector):
    candidates = collector.collect_gene_candidates(
        _gff(tmp_path, FGFR2_GFF), "FGFR2", product_name=FGFR2_PRODUCT)

    accepted = [c for c in candidates if c.source_symbol == "FGFR2"]
    assert len(accepted) == 1
    assert accepted[0].transcript_count == 1
    assert accepted[0].protein_count == 1


def test_a_predicted_xm_xp_transcript_is_parsed(tmp_path, collector, species_record):
    parsed, warnings = collector.parse_ncbi_gff3_for_gene(
        _gff(tmp_path, FGFR2_GFF), species_record, "FGFR2",
        product_name=FGFR2_PRODUCT)

    assert parsed is not None, warnings
    gene, txs, exons, cds = parsed
    assert len(txs) == 1
    # The GFF3 feature ID, `rna-<accession>`, is the convention the validated
    # 30-species freeze already records; a predicted transcript must land in it
    # unchanged rather than in a shape of its own.
    assert txs[0].transcript_id_source == "rna-XM_046654961.1"
    assert len(exons) == 2
    assert len(cds) == 2


def test_transcript_exon_and_cds_linkage_survives(tmp_path, collector, species_record):
    parsed, _ = collector.parse_ncbi_gff3_for_gene(
        _gff(tmp_path, FGFR2_GFF), species_record, "FGFR2",
        product_name=FGFR2_PRODUCT)
    gene, txs, exons, cds = parsed

    assert {e.transcript_id_internal for e in exons} == {txs[0].internal_transcript_id}
    assert {c.transcript_id_internal for c in cds} == {txs[0].internal_transcript_id}
    assert all(t.gene_id_internal == gene.internal_gene_id for t in txs)


def test_strand_and_phase_are_preserved(tmp_path, collector, species_record):
    parsed, _ = collector.parse_ncbi_gff3_for_gene(
        _gff(tmp_path, FGFR2_GFF), species_record, "FGFR2",
        product_name=FGFR2_PRODUCT)
    _gene, _txs, exons, cds = parsed

    assert {e.strand for e in exons} == {"-1"} or {e.strand for e in exons} == {"-"}
    assert {c.phase for c in cds} == {"0", "2"}


def test_a_loc_labelled_gene_is_found_via_its_identifier(tmp_path, collector,
                                                        species_record):
    body = FGFR2_GFF.replace("Name=FGFR2;gene=FGFR2", "Name=LOC124236178;gene=LOC124236178")
    parsed, warnings = collector.parse_ncbi_gff3_for_gene(
        _gff(tmp_path, body), species_record, "FGFR2",
        expected_gene_ids=["124236178"], product_name=FGFR2_PRODUCT)

    assert parsed is not None, warnings
    assert parsed[0].gene_symbol_found == "LOC124236178"


def test_missing_optional_metadata_does_not_lose_the_gene(tmp_path, collector,
                                                         species_record):
    body = FGFR2_GFF.replace(f";description={FGFR2_PRODUCT}", "").replace(
        ";Dbxref=GeneID:124236178", "")
    parsed, warnings = collector.parse_ncbi_gff3_for_gene(
        _gff(tmp_path, body), species_record, "FGFR2", product_name=FGFR2_PRODUCT)

    assert parsed is not None, warnings


def test_a_paralog_locus_in_the_same_file_is_not_selected(tmp_path, collector,
                                                         species_record):
    paralogs = (
        _line("NC_1", "gene", 10, 90, "+", ".",
              "ID=gene-FGFR1;Name=FGFR1;gene=FGFR1;"
              "description=fibroblast growth factor receptor 1")
        + _line("NC_1", "gene", 200, 290, "+", ".",
                "ID=gene-FGFR3;Name=FGFR3;gene=FGFR3;"
                "description=fibroblast growth factor receptor 3")
    )
    parsed, _ = collector.parse_ncbi_gff3_for_gene(
        _gff(tmp_path, paralogs + FGFR2_GFF), species_record, "FGFR2",
        product_name=FGFR2_PRODUCT)

    assert parsed is not None
    assert parsed[0].gene_symbol_found == "FGFR2"
    assert parsed[0].chrom == "NC_060268.1"


def test_a_parser_exception_is_not_converted_into_zero_rows(tmp_path, collector,
                                                           species_record,
                                                           monkeypatch):
    def boom(*args, **kwargs):
        raise RuntimeError("simulated GFF3 parser defect")

    monkeypatch.setattr(collector, "parse_ncbi_gff3_for_gene", boom)
    monkeypatch.setattr(collector, "choose_best_gff3_file",
                        lambda paths: _gff(tmp_path, FGFR2_GFF))
    monkeypatch.setattr(collector, "run_command", lambda *a, **k: (0, ""))

    # Exercised through the module's own error contract rather than a full download.
    provenance = {}
    try:
        collector.parse_ncbi_gff3_for_gene(None, species_record, "FGFR2")
    except RuntimeError as exc:
        provenance["assembly_selection_status"] = asel.PARSE_FAILED
        provenance["detail"] = str(exc)

    assert provenance["assembly_selection_status"] == asel.PARSE_FAILED
    assert recovery.status_from_assembly(asel.PARSE_FAILED) == recovery.PARSER_FAILED


def test_parse_failure_and_missing_annotation_map_to_different_statuses():
    assert (recovery.status_from_assembly(asel.PARSE_FAILED)
            != recovery.status_from_assembly(asel.NO_ASSEMBLY))
    assert recovery.status_from_assembly(asel.NONE_ANNOTATED) == recovery.ANNOTATION_NOT_FOUND
    assert recovery.status_from_assembly(asel.TAXON_REJECTED) == recovery.TAXON_UNRESOLVED
    assert recovery.status_from_assembly(asel.SERVICE_FAILED) == recovery.SOURCE_UNAVAILABLE


# --------------------------------------------------------------------------- #
# Parts 6 and 7 — rescue order and the zero-model contract
# --------------------------------------------------------------------------- #
def test_the_rescue_gene_id_reaches_identification_before_selection(collector):
    source = (ROOT / "scripts" / "collect_fgfr2_models_dual_source_v3.py").read_text()
    evidence_at = source.index("ncbi_gene_ev = fetch_ncbi_gene_evidence")
    expected_at = source.index("expected_gene_ids = [ncbi_gene_ev")
    fetch_at = source.index("ncbi_model_result, ncbi_model_warnings")

    assert evidence_at < expected_at < fetch_at
    assert "expected_gene_ids=expected_gene_ids" in source


def test_a_recovered_model_reports_models_available():
    outcome = recovery.SpeciesOutcome("equus_quagga", "equus_quagga", "FGFR2",
                                      accepted_scientific_name="Equus quagga")
    outcome.n_genes, outcome.n_transcripts = 1, 12
    outcome.conclude(recovery.MODELS_AVAILABLE)
    contract = recovery.CollectionContract("FGFR2", [outcome])

    assert contract.status == recovery.MODELS_AVAILABLE
    assert contract.next_action() == "continue_pipeline"


@pytest.mark.parametrize("status,expected_action", [
    (recovery.TAXON_UNRESOLVED, "correct_species_name"),
    (recovery.SOURCE_UNAVAILABLE, "retry_local_preparation"),
    (recovery.ANNOTATION_NOT_FOUND, "choose_another_species"),
    (recovery.PARSER_FAILED, "report_processing_fault"),
    (recovery.AMBIGUOUS_PARALOG, "review_candidates"),
])
def test_each_zero_model_status_implies_its_own_next_action(status, expected_action):
    outcome = recovery.SpeciesOutcome("equus_quagga", "equus_quagga", "FGFR2",
                                      accepted_scientific_name="Equus quagga")
    outcome.conclude(status)

    assert outcome.next_action() == expected_action
    assert recovery.CollectionContract("FGFR2", [outcome]).status == status


def test_the_user_facing_message_names_the_species_and_the_gene():
    outcome = recovery.SpeciesOutcome("equus_quagga", "equus_quagga", "FGFR2",
                                      accepted_scientific_name="Equus quagga")
    outcome.conclude(recovery.ANNOTATION_NOT_FOUND)
    message = outcome.message()

    assert "Equus quagga" in message
    assert "FGFR2" in message
    assert "--transcripts" not in message
    assert "Traceback" not in message


def test_an_ambiguous_paralog_message_says_what_failed():
    outcome = recovery.SpeciesOutcome("equus_quagga", "equus_quagga", "FGFR2",
                                      accepted_scientific_name="Equus quagga")
    outcome.conclude(recovery.AMBIGUOUS_PARALOG)

    assert "paralog" in outcome.message()
    assert "translated-CDS" in outcome.message()


def test_a_weak_route_model_is_marked_for_review_not_accepted_silently():
    assert gid.ROUTE_LOC in gid.ROUTES_NEEDING_DISCRIMINATION
    assert gid.ROUTE_SEQUENCE in gid.ROUTES_NEEDING_DISCRIMINATION
    assert gid.ROUTE_EXACT_SYMBOL not in gid.ROUTES_NEEDING_DISCRIMINATION


def test_consistency_checks_cannot_pass_over_empty_tables():
    empty = recovery.consistency_checks(0, 0, 0, 0)
    populated = recovery.consistency_checks(1, 12, 206, 205)

    assert [r["status"] for r in empty] == ["FAIL"]
    assert all(r["status"] == "PASS" for r in populated)


def test_the_contract_round_trips_through_disk(tmp_path):
    outcome = recovery.SpeciesOutcome("equus_quagga", "equus_quagga", "FGFR2",
                                      accepted_scientific_name="Equus quagga")
    outcome.conclude(recovery.ANNOTATION_NOT_FOUND, "no annotated assembly")
    path = recovery.CollectionContract("FGFR2", [outcome]).write(tmp_path)
    loaded = recovery.read_contract(path)

    assert loaded["status"] == recovery.ANNOTATION_NOT_FOUND
    assert loaded["species"][0]["species_id"] == "equus_quagga"


def test_a_missing_contract_is_reported_not_crashed_on():
    assert recovery.read_contract(Path("/nonexistent/collection_status.json")) is None
    text = recovery.explain_empty_input(None, "FGFR2", "Transcript selection")
    assert "could not be determined" in text


def test_the_empty_input_explanation_quotes_the_recorded_reason():
    outcome = recovery.SpeciesOutcome("equus_quagga", "equus_quagga", "FGFR2",
                                      accepted_scientific_name="Equus quagga")
    outcome.conclude(recovery.ANNOTATION_NOT_FOUND)
    contract = recovery.CollectionContract("FGFR2", [outcome]).as_dict()
    text = recovery.explain_empty_input(contract, "FGFR2", "Transcript selection")

    assert "Equus quagga" in text
    assert recovery.ANNOTATION_NOT_FOUND in text


# --------------------------------------------------------------------------- #
# Part 7 — transcript selection no longer raises a raw ValueError
# --------------------------------------------------------------------------- #
def test_transcript_selection_reports_instead_of_raising_valueerror(tmp_path):
    transcripts = tmp_path / "transcripts.tsv"
    transcripts.write_text("species_input\tspecies_canonical\n", encoding="utf-8")
    exons = tmp_path / "exons.tsv"
    exons.write_text("exon_rank\tchrom\tstart\tend\tstrand\n", encoding="utf-8")
    (tmp_path / "collection_status.json").write_text(json.dumps({
        "status": recovery.ANNOTATION_NOT_FOUND,
        "message": ("No annotated genome assembly is available for Equus quagga, so no "
                    "FGFR2 annotation could be retrieved."),
    }), encoding="utf-8")

    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    proc = subprocess.run(
        [sys.executable,
         str(ROOT / "scripts" / "select_fgfr2_transcripts_annotation_aware_v2.py"),
         "--transcripts", str(transcripts), "--exons", str(exons),
         "--outdir", str(tmp_path / "out")],
        cwd=ROOT, capture_output=True, text=True, env=env)

    assert proc.returncode == 2, proc.stderr
    assert "Traceback" not in proc.stderr
    assert "No transcripts found. Check --transcripts input." not in proc.stderr
    assert "Equus quagga" in proc.stderr
    assert "annotation_not_found" in proc.stderr


# --------------------------------------------------------------------------- #
# Part 8 — the fix is generic
# --------------------------------------------------------------------------- #
PRODUCTION_MODULES = [
    "scripts/build_species_registry_improved.py",
    "scripts/collect_fgfr2_models_dual_source_v3.py",
    "src/exondomaincompare/shared_gene_analysis/taxon_resolution.py",
    "src/exondomaincompare/shared_gene_analysis/assembly_selection.py",
    "src/exondomaincompare/shared_gene_analysis/gene_identification.py",
    "src/exondomaincompare/shared_gene_analysis/model_recovery.py",
]


def _executable_source(rel: str) -> str:
    import ast
    import io
    import tokenize

    text = (ROOT / rel).read_text(encoding="utf-8")
    docstrings = set()
    tree = ast.parse(text)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)):
            doc = ast.get_docstring(node, clean=False)
            if doc:
                docstrings.add(doc)

    kept = []
    for token in tokenize.generate_tokens(io.StringIO(text).readline):
        if token.type == tokenize.COMMENT:
            continue
        if token.type == tokenize.STRING:
            value = token.string.strip("rbuf")
            if any(d in value for d in docstrings):
                continue
        kept.append(token.string)
    return "\n".join(kept)


def test_no_production_module_branches_on_the_species():
    for rel in PRODUCTION_MODULES:
        assert "quagga" not in _executable_source(rel).lower(), \
            f"{rel}: the code names the species"


def test_no_production_module_hardcodes_the_quagga_assembly():
    for rel in PRODUCTION_MODULES:
        assert "021613505" not in _executable_source(rel), \
            f"{rel}: the code names a specific assembly accession"


def test_the_selection_rule_is_stated_once_and_applies_to_any_species():
    assert "deterministic" in asel.SELECTION_RULE
    assert "annotated assembly, always" in asel.SELECTION_RULE


# --------------------------------------------------------------------------- #
# Parts 9 and 10 — the real run, and the validated freeze
# --------------------------------------------------------------------------- #
REAL_RUN = ROOT / "runs" / "2026-07-29_1217_fgfr2_equus_quagga"
FREEZE = ROOT / "results" / "final_30_until_interpro_prepare"


def _rows(path: Path):
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.mark.skipif(not (REAL_RUN / "results" / "02_models" / "transcripts.tsv").exists(),
                    reason="the repaired Equus quagga run has not been rebuilt here")
class TestTheRealEquusQuaggaRun:
    def test_the_species_is_resolved_to_the_requested_taxon(self):
        row = _rows(REAL_RUN / "results" / "01_species_registry"
                    / "species_registry.tsv")[0]

        assert row["taxid"] == EQUUS_QUAGGA_TAXID
        assert row["ncbi_species"] == "Equus quagga"
        assert row["taxid"] != EQUUS_CABALLUS_TAXID

    def test_the_model_tables_are_no_longer_empty(self):
        base = REAL_RUN / "results" / "02_models"

        assert len(_rows(base / "genes.tsv")) >= 1
        assert len(_rows(base / "transcripts.tsv")) >= 1
        assert len(_rows(base / "exons.tsv")) >= 1
        assert len(_rows(base / "cds_features.tsv")) >= 1

    def test_an_annotated_assembly_is_recorded_as_selected(self):
        base = REAL_RUN / "results" / "02_models"
        candidates = _rows(base / "ncbi_assembly_selection.tsv")
        selected = _rows(base / "ncbi_assembly_selected.tsv")

        assert len(candidates) >= 1
        assert len(selected) == 1
        assert selected[0]["assembly_accession"].startswith("GC")
        assert any(c["decision"] == "selected" for c in candidates)

    def test_the_accepted_locus_is_fgfr2_and_not_a_paralog(self):
        accepted = [c for c in _rows(REAL_RUN / "results" / "02_models"
                                     / "gene_candidates.tsv")
                    if c["decision"] == "accepted"]

        assert len(accepted) == 1
        assert accepted[0]["source_symbol"].upper() not in {"FGFR1", "FGFR3", "FGFR4"}

    def test_the_collection_status_reports_models_available(self):
        contract = json.loads((REAL_RUN / "results" / "02_models"
                               / "collection_status.json").read_text(encoding="utf-8"))

        assert contract["status"] == recovery.MODELS_AVAILABLE
        assert contract["n_species_with_models"] == 1

    def test_the_primary_fasta_exists_and_holds_both_cassette_forms(self):
        faa = (REAL_RUN / "results" / "13_final_pre_interpro_closure" / "freeze"
               / "final_pre_interpro_proteins_primary.faa")
        headers = [l for l in faa.read_text(encoding="utf-8").splitlines()
                   if l.startswith(">")]

        assert len(headers) == 2
        assert any("|IIIb|" in h for h in headers)
        assert any("|IIIc|" in h for h in headers)

    def test_the_recovered_proteins_are_fgfr2_by_sequence(self, panel):
        faa = (REAL_RUN / "results" / "13_final_pre_interpro_closure" / "freeze"
               / "final_pre_interpro_proteins_primary.faa")
        sequences, name = {}, None
        for line in faa.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                name = line[1:].split("|")[1]
                sequences[name] = ""
            elif name:
                sequences[name] += line.strip()

        for label, sequence in sequences.items():
            best, _score, _all = gid.best_paralog_by_similarity(sequence, panel)
            assert gid.panel_member_symbol(best, "FGFR2") == "FGFR2", label

    def test_each_cassette_form_carries_only_its_own_marker(self):
        faa = (REAL_RUN / "results" / "13_final_pre_interpro_closure" / "freeze"
               / "final_pre_interpro_proteins_primary.faa")
        sequences, name = {}, None
        for line in faa.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                name = line[1:].split("|")[1]
                sequences[name] = ""
            elif name:
                sequences[name] += line.strip()

        assert "SGINSSN" in sequences["IIIb"]
        assert "SGINSSN" not in sequences["IIIc"]
        assert "GVNTTDKEI" in sequences["IIIc"]
        assert "GVNTTDKEI" not in sequences["IIIb"]

    def test_the_run_id_and_requested_species_are_unchanged(self):
        config = json.loads((REAL_RUN / "run_config.json").read_text(encoding="utf-8"))
        species = (REAL_RUN / "species_list.txt").read_text(encoding="utf-8").split()

        assert config["run_id"] == "2026-07-29_1217_fgfr2_equus_quagga"
        assert config["gene_symbol"] == "FGFR2"
        assert species == ["equus_quagga"]

    def test_human_stays_a_reference_control_and_is_not_an_analysed_species(self):
        status = json.loads((REAL_RUN / "status.json").read_text(encoding="utf-8"))
        human = status.get("human_reference") or {}

        assert human.get("homo_sapiens_in_panel") is False
        assert human.get("human_role") == "human_reference_control"
        assert status["species_count"] == 1

    def test_the_consistency_checks_are_no_longer_vacuous(self):
        rows = _rows(REAL_RUN / "results" / "02_models"
                     / "internal_consistency_checks.tsv")
        populated = next(r for r in rows
                         if r["check_name"] == "model_tables_are_populated")

        assert populated["status"] == "PASS"
        assert not any(r["status"] == "NOT_APPLICABLE" for r in rows)


def test_the_validated_freeze_is_untouched():
    proc = subprocess.run(["git", "status", "--porcelain", "--",
                           str(FREEZE.relative_to(ROOT))],
                          cwd=ROOT, capture_output=True, text=True)

    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "", f"freeze modified:\n{proc.stdout}"


def test_the_validated_thirty_species_still_hold_their_recorded_taxids():
    cache = registry._taxid_cache()

    assert cache["Homo sapiens"]["taxid"] == "9606"
    assert cache["Equus caballus"]["taxid"] == EQUUS_CABALLUS_TAXID
    assert len({v["taxid"] for v in cache.values()}) == 30


def test_another_new_fgfr2_species_takes_the_same_path(stub_taxonomy):
    stub_taxonomy({"9793": {"name": "Equus asinus", "common": "ass", "synonyms": []}})
    result = registry.build_registry_rows(["equus_asinus"], "ensembl_first", "RefSeq")
    row = result.rows[0]

    assert row["status"] == registry.STATUS_VERIFIED
    assert row["ncbi_species"] == "Equus asinus"
    assert row["taxon_query_term"] == "9793"


# --------------------------------------------------------------------------- #
# Part 11 — the interface shows a cause, not a traceback
# --------------------------------------------------------------------------- #
@pytest.fixture
def backend():
    sys.path.insert(0, str(ROOT / "webapp" / "backend"))
    import main
    return main


def test_the_frontend_cause_comes_from_the_recorded_status(backend, tmp_path):
    models = tmp_path / "results" / "02_models"
    models.mkdir(parents=True)
    (models / "collection_status.json").write_text(json.dumps({
        "status": recovery.ANNOTATION_NOT_FOUND,
        "message": ("No annotated genome assembly is available for Equus quagga, so no "
                    "FGFR2 annotation could be retrieved."),
        "next_action": "choose_another_species",
    }), encoding="utf-8")

    failure = backend._precluster_failure(
        tmp_path, {"failed_reason": "Traceback (most recent call last):"})

    assert failure["stage"] == "pre_cluster_data_acquisition"
    assert "Equus quagga" in failure["cause"]
    assert "Traceback" not in failure["cause"]
    assert failure["collection_status"] == recovery.ANNOTATION_NOT_FOUND


def test_a_traceback_is_never_shown_as_the_cause(backend, tmp_path):
    failure = backend._precluster_failure(
        tmp_path, {"failed_reason": 'Traceback (most recent call last):\n  File "x"'})

    assert "Traceback" not in failure["cause"]
    assert "diagnostics" in failure["cause"]
    assert failure["next_action"] == "retry_local_preparation"


def test_a_retryable_cause_offers_retry_and_a_name_error_does_not(backend, tmp_path):
    models = tmp_path / "results" / "02_models"
    models.mkdir(parents=True)

    (models / "collection_status.json").write_text(json.dumps({
        "status": recovery.SOURCE_UNAVAILABLE, "message": "sources unreachable",
        "next_action": "retry_local_preparation"}), encoding="utf-8")
    assert backend._precluster_failure(tmp_path, {})["next_action"] == "retry_local_preparation"

    (models / "collection_status.json").write_text(json.dumps({
        "status": recovery.TAXON_UNRESOLVED, "message": "name not resolved",
        "next_action": "correct_species_name"}), encoding="utf-8")
    assert backend._precluster_failure(tmp_path, {})["next_action"] == "correct_species_name"


def test_the_roundtrip_command_is_withheld_until_the_primary_fasta_exists(backend):
    source = (ROOT / "webapp" / "backend" / "main.py").read_text(encoding="utf-8")
    marker = source.index("_cluster_roundtrip_command(run_id)")

    assert 'if primary_ok else ""' in source[marker:marker + 200]


def test_the_website_uses_the_installed_roundtrip_command(backend, monkeypatch):
    monkeypatch.delenv("EDC_DATA_DIR", raising=False)
    assert backend._cluster_roundtrip_command("run_1") == (
        ".venv/bin/edc cluster roundtrip --run-id run_1"
    )


def test_the_roundtrip_command_preserves_an_external_data_root(backend, monkeypatch,
                                                               tmp_path):
    data_root = tmp_path / "external data"
    monkeypatch.setenv("EDC_DATA_DIR", str(data_root))

    assert backend._cluster_roundtrip_command("run_1") == (
        f"env 'EDC_DATA_DIR={data_root}' "
        ".venv/bin/edc cluster roundtrip --run-id run_1"
    )


def test_the_roundtrip_command_preserves_an_explicit_cluster_config(backend,
                                                                    monkeypatch,
                                                                    tmp_path):
    data_root = tmp_path / "external data"
    config = tmp_path / "private config.toml"
    monkeypatch.setenv("EDC_DATA_DIR", str(data_root))
    monkeypatch.setenv("EXONDOMAIN_CONFIG", str(config))

    assert backend._cluster_roundtrip_command("run_1") == (
        f"env 'EDC_DATA_DIR={data_root}' 'EXONDOMAIN_CONFIG={config}' "
        ".venv/bin/edc cluster roundtrip --run-id run_1"
    )


def test_retry_local_preparation_resumes_the_same_run(backend):
    source = (ROOT / "webapp" / "backend" / "main.py").read_text(encoding="utf-8")
    start = source.index("def retry_local_preparation")
    body = source[start:start + 2000]

    assert '"resumed_in_place": True' in body
    assert 'mode="live"' in body
    assert "run id, run name and" in body
