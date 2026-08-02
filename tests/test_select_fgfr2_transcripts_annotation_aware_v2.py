import csv
import json
from argparse import Namespace

from script_paths import load_script_module

sel = load_script_module("select_fgfr2_transcripts_annotation_aware_v2.py", "sel_v2")


def write_tsv(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, '') for k in fields})


def read_tsv(path):
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def test_normalize_biotype_and_appris_are_strict():
    assert sel.normalize_biotype('protein coding') == 'coding'
    assert sel.normalize_biotype('mRNA') == 'coding'
    assert sel.normalize_biotype('retained_intron') == 'decay_or_retained_intron'
    assert sel.is_appris_principal_annotation({'appris': 'principal1'}) is True
    assert sel.is_appris_principal_annotation({'appris': 'alternative1'}) is False


def test_refseq_select_does_not_use_generic_select_column():
    assert sel.is_refseq_select_annotation({'select': '1'}) is False
    assert sel.is_refseq_select_annotation({'refseq_select': 'select'}) is True


def test_flexible_exon_linking_supports_multiple_column_names():
    exons = [
        {'internal_transcript_id': 'txA', 'exon_rank': '1', 'chrom': '1', 'start': '10', 'end': '20', 'strand': '+'},
        {'parent_transcript_id': 'txB,txC', 'exon_rank': '1', 'chrom': '1', 'start': '30', 'end': '40', 'strand': '+'},
        {'transcript_id_source': 'NM_001.2', 'exon_rank': '1', 'chrom': '1', 'start': '50', 'end': '60', 'strand': '+'},
    ]
    idx = sel.build_exon_index(exons)
    assert len(sel.exons_for({'internal_transcript_id': 'txA', 'transcript_id_source': ''}, idx)) == 1
    assert len(sel.exons_for({'internal_transcript_id': 'txB', 'transcript_id_source': ''}, idx)) == 1
    assert len(sel.exons_for({'internal_transcript_id': '', 'transcript_id_source': 'NM_001'}, idx)) == 1


def test_duplicate_transcripts_are_flagged():
    warnings = []
    rows = [
        {'internal_transcript_id': 'a', 'transcript_id_source': 'NM_1', 'species_canonical': 'X'},
        {'internal_transcript_id': 'a', 'transcript_id_source': 'NM_2', 'species_canonical': 'X'},
        {'internal_transcript_id': 'b', 'transcript_id_source': 'NM_1', 'species_canonical': 'X'},
    ]
    flags = sel.detect_duplicate_transcripts(rows, warnings)
    assert 'duplicate_internal_transcript_id' in flags['a']
    assert any(w['warning_code'] == 'duplicate_source_transcript_id' for w in warnings)


def test_integration_creates_reports_plots_and_overlap(tmp_path):
    transcripts = tmp_path / 'transcripts.tsv'
    exons = tmp_path / 'exons.tsv'
    annotations = tmp_path / 'annotations.tsv'
    iso = tmp_path / 'iso.tsv'
    domains = tmp_path / 'domains.tsv'
    out = tmp_path / 'out'

    tx_fields = ['species_input','species_canonical','source_db','gene_id_internal','transcript_id_source','internal_transcript_id','transcript_name','transcript_biotype','translation_id_source','protein_length_aa','is_canonical_source','support_level','completeness_flags']
    ex_fields = ['transcript_id_internal','exon_rank','chrom','start','end','strand']
    ann_fields = ['internal_transcript_id','mane_select','refseq_select','appris','transcript_name']
    iso_fields = ['internal_transcript_id','isoform_class','evidence']
    dom_fields = ['internal_transcript_id','domain_name','interpro_id']

    write_tsv(transcripts, [
        {'species_input':'Human','species_canonical':'Homo sapiens','source_db':'Ensembl','gene_id_internal':'g1','transcript_id_source':'ENST_REF','internal_transcript_id':'tx_ref','transcript_name':'FGFR2-IIIb','transcript_biotype':'protein_coding','translation_id_source':'prot1','protein_length_aa':'821','is_canonical_source':'1','support_level':'TSL1','completeness_flags':''},
        {'species_input':'Human','species_canonical':'Homo sapiens','source_db':'Ensembl','gene_id_internal':'g1','transcript_id_source':'ENST_C','internal_transcript_id':'tx_c','transcript_name':'FGFR2-IIIc','transcript_biotype':'protein_coding','translation_id_source':'prot2','protein_length_aa':'820','is_canonical_source':'0','support_level':'TSL1','completeness_flags':''},
        {'species_input':'Human','species_canonical':'Homo sapiens','source_db':'Ensembl','gene_id_internal':'g1','transcript_id_source':'ENST_ALT','internal_transcript_id':'tx_alt','transcript_name':'FGFR2 alt','transcript_biotype':'protein_coding','translation_id_source':'prot3','protein_length_aa':'790','is_canonical_source':'0','support_level':'TSL2','completeness_flags':''},
    ], tx_fields)
    write_tsv(exons, [
        {'transcript_id_internal':'tx_ref','exon_rank':'1','chrom':'10','start':'1','end':'100','strand':'+'},
        {'transcript_id_internal':'tx_c','exon_rank':'1','chrom':'10','start':'101','end':'200','strand':'+'},
        {'transcript_id_internal':'tx_alt','exon_rank':'1','chrom':'10','start':'201','end':'300','strand':'+'},
    ], ex_fields)
    write_tsv(annotations, [
        {'internal_transcript_id':'tx_ref','mane_select':'1','refseq_select':'','appris':'principal1','transcript_name':'FGFR2-IIIb'},
        {'internal_transcript_id':'tx_c','mane_select':'','refseq_select':'','appris':'principal2','transcript_name':'FGFR2-IIIc'},
    ], ann_fields)
    write_tsv(iso, [
        {'internal_transcript_id':'tx_ref','isoform_class':'IIIb','evidence':'curated mutually exclusive exon'},
        {'internal_transcript_id':'tx_c','isoform_class':'IIIc','evidence':'curated mutually exclusive exon'},
    ], iso_fields)
    write_tsv(domains, [
        {'internal_transcript_id':'tx_ref','domain_name':'Ig-like domain III','interpro_id':'IPR003599'},
        {'internal_transcript_id':'tx_c','domain_name':'Ig-like domain III','interpro_id':'IPR003599'},
    ], dom_fields)

    args = Namespace(transcripts=transcripts, exons=exons, annotations=annotations, domains=domains, isoform_evidence=iso, outdir=out, min_reference_score=350, max_alternatives_per_species=1, gene_symbol='FGFR2', strict=False, no_plots=False, plot_top_n=5)
    meta = sel.run(args)

    selected = read_tsv(out / 'selected_transcripts.tsv')
    roles = {r['selection_role'] for r in selected}
    assert {'reference', 'FGFR2_IIIb_candidate', 'FGFR2_IIIc_candidate'} <= roles
    ref_rows = [r for r in selected if r['internal_transcript_id'] == 'tx_ref']
    assert any(r['overlapping_roles'] == '1' for r in ref_rows)
    assert (out / 'transcript_selection_report.md').exists()
    assert (out / 'transcript_selection_report.html').exists()
    assert (out / 'run_metadata.json').exists()
    assert list((out / 'plots').glob('top_scores_*.png'))
    metadata = json.loads((out / 'run_metadata.json').read_text())
    assert metadata['output_row_counts']['summary'] == 1


def test_medium_isoform_hint_is_provisional_not_candidate(tmp_path):
    transcripts = tmp_path / 'transcripts.tsv'
    exons = tmp_path / 'exons.tsv'
    out = tmp_path / 'out'
    tx_fields = ['species_input','species_canonical','source_db','gene_id_internal','transcript_id_source','internal_transcript_id','transcript_name','transcript_biotype','translation_id_source','protein_length_aa','is_canonical_source','support_level','completeness_flags']
    ex_fields = ['transcript_id_internal','exon_rank','chrom','start','end','strand']
    write_tsv(transcripts, [
        {'species_input':'Test','species_canonical':'Test species','source_db':'NCBI','gene_id_internal':'g','transcript_id_source':'NM_B','internal_transcript_id':'txb','transcript_name':'FGFR2 isoform b','transcript_biotype':'mRNA','translation_id_source':'p','protein_length_aa':'800','is_canonical_source':'0','support_level':'','completeness_flags':''},
    ], tx_fields)
    write_tsv(exons, [{'transcript_id_internal':'txb','exon_rank':'1','chrom':'1','start':'1','end':'100','strand':'+'}], ex_fields)
    args = Namespace(transcripts=transcripts, exons=exons, annotations=None, domains=None, isoform_evidence=None, outdir=out, min_reference_score=0, max_alternatives_per_species=1, gene_symbol='FGFR2', strict=False, no_plots=True, plot_top_n=5)
    sel.run(args)
    selected = read_tsv(out / 'selected_transcripts.tsv')
    roles = {r['selection_role'] for r in selected}
    assert 'FGFR2_IIIb_candidate' not in roles
    assert 'FGFR2_IIIb_provisional' in roles
