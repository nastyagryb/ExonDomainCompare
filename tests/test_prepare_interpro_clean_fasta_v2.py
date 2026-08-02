import csv
import json

from script_paths import load_script_module


def load_module():
    return load_script_module("prepare_interpro_clean_fasta_v2.py", "prep")


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_parse_v2_header():
    mod = load_module()
    header = "fgfr2prot_000001|species=homo_sapiens|source=Ensembl|role=reference|transcript=ENST1|protein=ENSP1|isoform=IIIb"
    parsed = mod.parse_header(header)
    assert parsed["original_fasta_id"] == "fgfr2prot_000001"
    assert parsed["species_canonical"] == "homo_sapiens"
    assert parsed["source_db"] == "Ensembl"
    assert parsed["selection_role"] == "reference"
    assert parsed["transcript_id"] == "ENST1"
    assert parsed["protein_id"] == "ENSP1"
    assert parsed["isoform"] == "IIIb"


def test_deduplicates_sequences_and_writes_mapping(tmp_path):
    mod = load_module()
    seq = "M" + "A" * 599
    fasta = tmp_path / "input.faa"
    fasta.write_text(
        f">fgfr2prot_000001|species=homo_sapiens|source=Ensembl|role=reference|transcript=T1|protein=P1|isoform=IIIb\n{seq}\n"
        f">fgfr2prot_000002|species=homo_sapiens|source=Ensembl|role=FGFR2_IIIb_candidate|transcript=T1|protein=P1|isoform=IIIb\n{seq}\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "out"
    assert mod.main(["--input", str(fasta), "--outdir", str(outdir), "--prefix", "FGFR2", "--split_size", "1"]) == 0
    unique_fasta = outdir / "fgfr2_interpro_clean_unique.fasta"
    assert unique_fasta.read_text(encoding="utf-8").count(">") == 1
    mapping = read_tsv(outdir / "fgfr2_interpro_id_mapping.tsv")
    assert len(mapping) == 2
    assert mapping[0]["unique_id"] == mapping[1]["unique_id"]
    assert mapping[0]["duplicate_group_size"] == "2"


def test_split_unique_fasta(tmp_path):
    mod = load_module()
    fasta = tmp_path / "input.faa"
    seqs = ["M" + aa * 599 for aa in ["A", "C", "D"]]
    fasta.write_text("".join(f">id{i}|species=s|source=NCBI|role=reference\n{seq}\n" for i, seq in enumerate(seqs, start=1)), encoding="utf-8")
    outdir = tmp_path / "out"
    mod.main(["--input", str(fasta), "--outdir", str(outdir), "--split_size", "2"])
    assert (outdir / "fgfr2_interpro_unique_part01.fasta").exists()
    assert (outdir / "fgfr2_interpro_unique_part02.fasta").exists()
    assert (outdir / "fgfr2_interpro_unique_part01.fasta").read_text(encoding="utf-8").count(">") == 2
    assert (outdir / "fgfr2_interpro_unique_part02.fasta").read_text(encoding="utf-8").count(">") == 1


def test_invalid_sequence_strict_aborts_and_writes_warnings(tmp_path):
    mod = load_module()
    fasta = tmp_path / "bad.faa"
    fasta.write_text(">bad|species=s\nMAAA*123\n", encoding="utf-8")
    outdir = tmp_path / "out"
    try:
        mod.main(["--input", str(fasta), "--outdir", str(outdir), "--strict"])
    except SystemExit as e:
        assert "Strict mode" in str(e)
    else:
        raise AssertionError("strict mode should abort")
    warnings = read_tsv(outdir / "fgfr2_interpro_prepare_warnings.tsv")
    types = {row["warning_type"] for row in warnings}
    assert "invalid_amino_acid_symbols" in types
    assert "stop_symbol_in_protein_sequence" in types


def test_ambiguous_symbols_are_warning_not_error(tmp_path):
    mod = load_module()
    fasta = tmp_path / "amb.faa"
    fasta.write_text(">x|species=s\nM" + "A" * 100 + "X" + "C" * 500 + "\n", encoding="utf-8")
    outdir = tmp_path / "out"
    mod.main(["--input", str(fasta), "--outdir", str(outdir), "--strict"])
    warnings = read_tsv(outdir / "fgfr2_interpro_prepare_warnings.tsv")
    assert any(row["warning_type"] == "ambiguous_amino_acid_symbols" for row in warnings)
    assert (outdir / "fgfr2_interpro_clean_unique.fasta").read_text(encoding="utf-8").count(">") == 1


def test_metadata_and_reports_created(tmp_path):
    mod = load_module()
    fasta = tmp_path / "input.faa"
    fasta.write_text(">fgfr2prot_000001|species=homo_sapiens|source=Ensembl|role=reference\nM" + "A" * 599 + "\n", encoding="utf-8")
    outdir = tmp_path / "out"
    mod.main(["--input", str(fasta), "--outdir", str(outdir)])
    meta = json.loads((outdir / "run_metadata.json").read_text(encoding="utf-8"))
    assert meta["script_name"] == "prepare_interpro_clean_fasta_v2.py"
    assert meta["input_records_accepted"] == 1
    assert (outdir / "fgfr2_interpro_prepare_report.md").exists()
    assert (outdir / "fgfr2_interpro_prepare_report.html").exists()
