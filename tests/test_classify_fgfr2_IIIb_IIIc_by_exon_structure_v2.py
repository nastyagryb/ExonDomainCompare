import csv

import pytest

from script_paths import load_script_module

mod = load_script_module(
    "classify_fgfr2_IIIb_IIIc_by_exon_structure_v2_3_human_calibrated.py", "iso_v2"
)


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({f: r.get(f, '') for f in fields})


def read_tsv(path):
    with path.open(newline='', encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def tx_rows(strand='+'):
    return [
        dict(species_canonical='Homo sapiens', species_input='Homo sapiens', source_db='Ensembl', gene_id_internal='gene1', gene_symbol='FGFR2', internal_transcript_id='tx_b', transcript_id_source='T_B', transcript_biotype='protein_coding', protein_length_aa='821'),
        dict(species_canonical='Homo sapiens', species_input='Homo sapiens', source_db='Ensembl', gene_id_internal='gene1', gene_symbol='FGFR2', internal_transcript_id='tx_c', transcript_id_source='T_C', transcript_biotype='protein_coding', protein_length_aa='821'),
    ]


def make_exons_for_tx(tx, alt_start, strand='+', id_field='transcript_id_internal'):
    # 9 exons. exon 5 is the alternative exon, with shared exons 4 and 6 as anchors.
    starts = [100, 300, 500, 700, alt_start, 1300, 1500, 1700, 1900]
    if strand == '-':
        # Rank still transcript order; coordinates decrease along transcript order for minus strand.
        starts = [3000, 2800, 2600, 2400, alt_start, 1800, 1600, 1400, 1200]
    rows = []
    for i, s in enumerate(starts, start=1):
        length = 151 if i == 5 else 99
        rows.append({id_field: tx, 'exon_rank': str(i), 'chrom': 'chr10', 'start': str(s), 'end': str(s + length - 1), 'strand': strand})
    return rows


def test_plus_strand_classifies_lower_alt_as_iiib(tmp_path):
    transcripts = tx_rows('+')
    exons = make_exons_for_tx('tx_b', 900, '+') + make_exons_for_tx('tx_c', 1100, '+')
    exons_path = tmp_path / 'exons.tsv'
    tx_path = tmp_path / 'transcripts.tsv'
    out = tmp_path / 'out'
    write_tsv(tx_path, transcripts, list(transcripts[0].keys()))
    write_tsv(exons_path, exons, ['transcript_id_internal','exon_rank','chrom','start','end','strand'])
    mod.main(['--transcripts', str(tx_path), '--exons', str(exons_path), '--outdir', str(out)])
    rows = read_tsv(out / 'fgfr2_isoform_evidence.tsv')
    calls = {r['internal_transcript_id']: r['isoform_class'] for r in rows}
    assert calls['tx_b'] == 'IIIb'
    assert calls['tx_c'] == 'IIIc'
    assert rows[0]['assignment_rule'] == (
        'conservative_exon_order_rule_applied_only_to_exact_two_alt_exon_mutually_exclusive_slot;'
        'first=putative_FGFR2_exon8_IIIb;second=putative_FGFR2_exon9_IIIc'
    )
    # The lower-coordinate alternative exon is the first one in transcript order
    # on the plus strand and therefore the IIIb candidate.
    assert rows[0]['iiib_exon_sig'] == 'chr10:900-1050:+'
    assert rows[0]['iiic_exon_sig'] == 'chr10:1100-1250:+'
    legacy = {r['internal_transcript_id']: r['legacy_order_based_isoform_assignment'] for r in rows}
    assert legacy == {'tx_b': 'IIIb', 'tx_c': 'IIIc'}
    # No cds_features and no human cassette references are supplied here, so the
    # order rule must stay provisional instead of claiming sequence calibration.
    assert {r['direction_assignment_method'] for r in rows} == {'order_rule_provisional'}
    assert {r['direction_validation_status'] for r in rows} == {'direction_unresolved_no_sequence'}


def test_minus_strand_classifies_higher_coordinate_alt_as_iiib(tmp_path):
    transcripts = tx_rows('-')
    exons = make_exons_for_tx('tx_b', 2200, '-') + make_exons_for_tx('tx_c', 2000, '-')
    tx_path = tmp_path / 'transcripts.tsv'
    exons_path = tmp_path / 'exons.tsv'
    out = tmp_path / 'out'
    write_tsv(tx_path, transcripts, list(transcripts[0].keys()))
    write_tsv(exons_path, exons, ['transcript_id_internal','exon_rank','chrom','start','end','strand'])
    mod.main(['--transcripts', str(tx_path), '--exons', str(exons_path), '--outdir', str(out)])
    rows = read_tsv(out / 'fgfr2_isoform_evidence.tsv')
    calls = {r['internal_transcript_id']: r['isoform_class'] for r in rows}
    assert calls['tx_b'] == 'IIIb'
    assert calls['tx_c'] == 'IIIc'


def test_flexible_exon_id_linking_parent_transcript_id(tmp_path):
    transcripts = tx_rows('+')
    exons = make_exons_for_tx('tx_b', 900, '+', id_field='parent_transcript_id') + make_exons_for_tx('tx_c', 1100, '+', id_field='parent_transcript_id')
    tx_path = tmp_path / 'transcripts.tsv'
    exons_path = tmp_path / 'exons.tsv'
    out = tmp_path / 'out'
    write_tsv(tx_path, transcripts, list(transcripts[0].keys()))
    write_tsv(exons_path, exons, ['parent_transcript_id','exon_rank','chrom','start','end','strand'])
    mod.main(['--transcripts', str(tx_path), '--exons', str(exons_path), '--outdir', str(out)])
    rows = read_tsv(out / 'fgfr2_isoform_evidence.tsv')
    assert {r['isoform_class'] for r in rows} == {'IIIb', 'IIIc'}


def test_transcript_with_both_alternative_exons_is_ambiguous(tmp_path):
    transcripts = tx_rows('+') + [dict(tx_rows('+')[0], internal_transcript_id='tx_both', transcript_id_source='T_BOTH')]
    exons = make_exons_for_tx('tx_b', 900, '+') + make_exons_for_tx('tx_c', 1100, '+')
    both = make_exons_for_tx('tx_both', 900, '+')
    # insert second alternative exon into same slot between shared flanks; rank 6 and shift later ranks
    both = both[:5] + [{'transcript_id_internal':'tx_both','exon_rank':'6','chrom':'chr10','start':'1100','end':'1250','strand':'+'}] + [dict(r, exon_rank=str(int(r['exon_rank'])+1)) for r in both[5:]]
    exons += both
    tx_path = tmp_path / 'transcripts.tsv'
    exons_path = tmp_path / 'exons.tsv'
    out = tmp_path / 'out'
    write_tsv(tx_path, transcripts, list(transcripts[0].keys()))
    write_tsv(exons_path, exons, ['transcript_id_internal','exon_rank','chrom','start','end','strand'])
    mod.main(['--transcripts', str(tx_path), '--exons', str(exons_path), '--outdir', str(out)])
    rows = read_tsv(out / 'fgfr2_isoform_evidence.tsv')
    call = {r['internal_transcript_id']: r['isoform_class'] for r in rows}
    assert call['tx_both'] == 'ambiguous'


def test_more_than_two_alt_exons_warns(tmp_path):
    transcripts = tx_rows('+') + [dict(tx_rows('+')[0], internal_transcript_id='tx_x', transcript_id_source='T_X')]
    exons = make_exons_for_tx('tx_b', 900, '+') + make_exons_for_tx('tx_c', 1100, '+') + make_exons_for_tx('tx_x', 1000, '+')
    tx_path = tmp_path / 'transcripts.tsv'
    exons_path = tmp_path / 'exons.tsv'
    out = tmp_path / 'out'
    write_tsv(tx_path, transcripts, list(transcripts[0].keys()))
    write_tsv(exons_path, exons, ['transcript_id_internal','exon_rank','chrom','start','end','strand'])
    mod.main(['--transcripts', str(tx_path), '--exons', str(exons_path), '--outdir', str(out)])
    warnings = read_tsv(out / 'fgfr2_isoform_warnings.tsv')
    assert any(w['warning_type'] == 'slot_with_more_than_two_alternative_exons' for w in warnings)


def test_strict_mode_aborts_on_missing_schema(tmp_path):
    tx_path = tmp_path / 'bad_tx.tsv'
    exons_path = tmp_path / 'bad_exons.tsv'
    out = tmp_path / 'out'
    write_tsv(tx_path, [{'x':'1'}], ['x'])
    write_tsv(exons_path, [{'y':'1'}], ['y'])
    with pytest.raises(SystemExit):
        mod.main(['--transcripts', str(tx_path), '--exons', str(exons_path), '--outdir', str(out), '--strict'])
    warnings = read_tsv(out / 'fgfr2_isoform_warnings.tsv')
    assert any(w['warning_type'] == 'missing_required_column' for w in warnings)


def test_reports_metadata_and_plot_are_written(tmp_path):
    transcripts = tx_rows('+')
    exons = make_exons_for_tx('tx_b', 900, '+') + make_exons_for_tx('tx_c', 1100, '+')
    tx_path = tmp_path / 'transcripts.tsv'
    exons_path = tmp_path / 'exons.tsv'
    out = tmp_path / 'out'
    write_tsv(tx_path, transcripts, list(transcripts[0].keys()))
    write_tsv(exons_path, exons, ['transcript_id_internal','exon_rank','chrom','start','end','strand'])
    mod.main(['--transcripts', str(tx_path), '--exons', str(exons_path), '--outdir', str(out)])
    assert (out / 'run_metadata.json').exists()
    assert (out / 'fgfr2_isoform_report.md').exists()
    assert (out / 'fgfr2_isoform_report.html').exists()
    # matplotlib is a declared project dependency, so the summary plot is required.
    assert (out / 'plots' / 'isoform_summary_by_species.png').exists()
