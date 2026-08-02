from script_paths import load_script_module

mod = load_script_module("collect_fgfr2_models_dual_source_v3.py", "collector_v3")


def species(preferred_source='ensembl_first'):
    return mod.SpeciesRecord(
        input_name='Homo sapiens',
        ensembl_species='homo_sapiens',
        ncbi_species='Homo sapiens',
        taxid='9606',
        preferred_source=preferred_source,
        assembly_preference='RefSeq',
    )


def gene(source, confidence='high', biotype='protein_coding'):
    return mod.GeneModel(
        source_db=source,
        species_input='Homo sapiens',
        species_canonical='homo_sapiens',
        taxid='9606',
        assembly_accession='',
        assembly_name='',
        gene_symbol_requested='FGFR2',
        gene_symbol_found='FGFR2',
        gene_id_source=f'{source}_gene',
        gene_biotype=biotype,
        chrom='10',
        start='1000',
        end='5000',
        strand='1',
        model_status='model_found_high_confidence',
        model_confidence=confidence,
        lookup_method='test',
        internal_gene_id=f'homo_sapiens|{source}|FGFR2',
    )


def tx(source, gene_id, tx_id, length='821'):
    return mod.TranscriptModel(
        source_db=source,
        species_canonical='homo_sapiens',
        gene_id_internal=gene_id,
        transcript_id_source=tx_id,
        transcript_name=tx_id,
        transcript_biotype='protein_coding',
        translation_id_source=f'{tx_id}_prot',
        protein_length_aa=length,
        protein_length_source='test',
        is_canonical_source='',
        support_level='',
        completeness_flags='',
        transcript_model_confidence='high',
        chrom='10',
        start='1000',
        end='5000',
        strand='1',
        internal_transcript_id=f'{gene_id}|tx|{tx_id}',
    )


def exon(source, tx_id, rank):
    return mod.ExonModel(
        source_db=source,
        species_canonical='homo_sapiens',
        transcript_id_internal=tx_id,
        exon_id_source=f'ex{rank}',
        exon_rank=str(rank),
        chrom='10',
        start=str(1000 + rank * 100),
        end=str(1050 + rank * 100),
        strand='1',
        phase='',
        end_phase='',
        internal_exon_id=f'{tx_id}|exon|{rank}',
    )


def test_gff3_attributes_are_url_decoded_and_split():
    attrs = mod.parse_gff3_attributes('ID=gene-FGFR2;gene_synonym=BEK%2CKGFR%2CFGFR2;Name=FGFR2')
    assert attrs['gene_synonym'] == 'BEK,KGFR,FGFR2'
    assert mod.split_gff3_multi_value(attrs['gene_synonym']) == ['BEK', 'KGFR', 'FGFR2']
    assert mod.symbol_matches_gene(attrs, 'FGFR2') is True


def test_ncbi_gff3_parser_handles_synonyms_and_multi_parent(tmp_path):
    gff = tmp_path / 'genomic.gff3'
    gff.write_text(
        '##gff-version 3\n'
        'NC_000010.11\tRefSeq\tgene\t100\t1000\t.\t+\t.\tID=gene-1;Name=not_main;gene_synonym=BEK,FGFR2,KGFR;gbkey=Gene\n'
        'NC_000010.11\tRefSeq\tmRNA\t100\t1000\t.\t+\t.\tID=tx1;Parent=gene-1;Name=tx1;protein_id=XP_1\n'
        'NC_000010.11\tRefSeq\tmRNA\t100\t1000\t.\t+\t.\tID=tx2;Parent=gene-1;Name=tx2;protein_id=XP_2\n'
        'NC_000010.11\tRefSeq\texon\t100\t199\t.\t+\t.\tID=ex_shared;Parent=tx1,tx2\n'
        'NC_000010.11\tRefSeq\tCDS\t100\t199\t.\t+\t0\tID=cds_shared;Parent=tx1,tx2;protein_id=XP_shared\n',
        encoding='utf-8',
    )
    parsed, warnings = mod.parse_ncbi_gff3_for_gene(gff, species(), 'FGFR2')
    assert parsed is not None
    # The parser also returns the CDS features it collected.
    gene_model, txs, exons, _cds = parsed
    assert gene_model.gene_biotype == 'coding_gene'
    assert len(txs) == 2
    assert len(exons) == 2
    assert all(t.protein_length_source == 'NCBI_GFF3_CDS_span_estimate_stop_excluded'
               for t in txs)
    # 100 bp is not a whole number of codons, so this CDS is incomplete and carries
    # no terminator to discount: 33 codons, 33 residues. The terminator is only
    # subtracted from a complete CDS, where counting it as a residue is what let an
    # exon project one amino acid past the end of its protein.
    assert all(t.protein_length_aa == '33' for t in txs)
    assert not any('expected exactly 1' in w for w in warnings)


def test_ncbi_model_confidence_drops_when_parser_warnings_exist(tmp_path):
    gff = tmp_path / 'genomic.gff3'
    gff.write_text(
        '##gff-version 3\n'
        'ctg\tRefSeq\tgene\t100\t1000\t.\t+\t.\tID=gene-1;Name=FGFR2;gbkey=Gene\n'
        'ctg\tRefSeq\tCDS\t100\t399\t.\t+\t0\tID=cds1;Parent=synthetic_tx;protein_id=XP_1\n',
        encoding='utf-8',
    )
    parsed, warnings = mod.parse_ncbi_gff3_for_gene(gff, species(), 'FGFR2')
    assert parsed is not None
    # The parser also returns the CDS features it collected.
    gene_model, txs, exons, _cds = parsed
    assert any('synthesizing transcripts' in w for w in warnings)
    assert gene_model.model_confidence == 'medium'


def test_assembly_selection_honors_refseq_preference():
    candidates = [
        {'accession': 'GCA_000001', 'assembly_name': 'genbank', 'assembly_level': 'Chromosome', 'refseq_category': '', 'annotated': '1'},
        {'accession': 'GCF_000002', 'assembly_name': 'refseq', 'assembly_level': 'Scaffold', 'refseq_category': 'representative genome', 'annotated': '1'},
    ]
    best, scored = mod.pick_best_ncbi_assembly_from_summary(candidates, 'RefSeq')
    assert best['accession'] == 'GCF_000002'
    assert 'matches_species_registry_assembly_preference' in best['assembly_decision_notes']


def test_choose_best_gff3_file_prefers_genomic_and_larger(tmp_path):
    a = tmp_path / 'small.gff3'
    b = tmp_path / 'genomic.gff'
    a.write_text('x', encoding='utf-8')
    b.write_text('x' * 10, encoding='utf-8')
    assert mod.choose_best_gff3_file([a, b]) == b


def test_compare_models_applies_ncbi_first_preference():
    eg = gene('Ensembl')
    ng = gene('NCBI', biotype='Gene')
    et = [tx('Ensembl', eg.internal_gene_id, 'e1')]
    nt = [tx('NCBI', ng.internal_gene_id, 'n1')]
    ee = [exon('Ensembl', et[0].internal_transcript_id, 1)]
    ne = [exon('NCBI', nt[0].internal_transcript_id, 1)]
    selected, conf, reason, row, warnings = mod.compare_models(eg, et, ee, ng, nt, ne, preferred_source='ncbi_first')
    assert selected == 'NCBI'
    assert row['preferred_source_applied'] == 'ncbi_first'
    assert int(row['ncbi_source_score']) > int(row['ensembl_source_score'])
    assert row['same_biotype'] == '1'


def test_compare_models_applies_strict_ensembl_preference():
    eg = gene('Ensembl', confidence='medium')
    ng = gene('NCBI', confidence='high')
    et = [tx('Ensembl', eg.internal_gene_id, 'e1', length='760')]
    nt = [tx('NCBI', ng.internal_gene_id, 'n1', length='900'), tx('NCBI', ng.internal_gene_id, 'n2', length='880')]
    selected, conf, reason, row, warnings = mod.compare_models(eg, et, [], ng, nt, [], preferred_source='strict_ensembl')
    assert selected == 'Ensembl'
    assert reason == 'strict_ensembl_preference_applied'
