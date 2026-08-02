import csv
import json
import subprocess
import sys

from script_paths import load_script_module, script_path

SCRIPT_NAME = "export_selected_fgfr2_proteins_complete_v2_1_region_qc.py"
SCRIPT = script_path(SCRIPT_NAME)
mod = load_script_module(SCRIPT_NAME, "export_v2")


def write_tsv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def make_seq(n=841):
    aas = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(aas[i % len(aas)] for i in range(n))


def make_ncbi_cache(tmp_path, protein_acc="XP_001.1", tx="XM_001.1", seq_len=841, product="fibroblast growth factor receptor 2 isoform X1"):
    data = tmp_path / "cache" / "ncbi_dataset" / "data" / "GCF_TEST"
    data.mkdir(parents=True)
    gff = data / "genomic.gff"
    gff.write_text(
        "##gff-version 3\n"
        f"chr\tRefSeq\tCDS\t1\t100\t.\t+\t0\tID=cds-{protein_acc};Parent=rna-{tx};protein_id={protein_acc};product={product}\n",
        encoding="utf-8",
    )
    faa = data / "protein.faa"
    faa.write_text(f">{protein_acc} {product} [Homo sapiens]\n{make_seq(seq_len)}\n", encoding="utf-8")
    return tmp_path / "cache"


def make_selected(tmp_path, rows):
    fields = [
        "species_input", "species_canonical", "source_db", "selection_role", "internal_transcript_id",
        "transcript_id_source", "translation_id_source", "protein_length_aa", "iii_isoform_assignment"
    ]
    path = tmp_path / "selected.tsv"
    write_tsv(path, rows, fields)
    return path


def run_export(tmp_path, selected, cache, *extra):
    out = tmp_path / "proteins.faa"
    report = tmp_path / "protein_export_report.tsv"
    cmd = [sys.executable, str(SCRIPT), "--selected", str(selected), "--cache", str(cache), "--out", str(out), "--report", str(report)] + list(extra)
    return subprocess.run(cmd, text=True, capture_output=True), out, report


def test_ncbi_exact_accession_export_writes_unique_id_and_high_confidence(tmp_path):
    cache = make_ncbi_cache(tmp_path)
    selected = make_selected(tmp_path, [{
        "species_input": "Homo sapiens", "species_canonical": "homo_sapiens", "source_db": "NCBI",
        "selection_role": "reference", "internal_transcript_id": "tx1", "transcript_id_source": "XM_001.1",
        "translation_id_source": "", "protein_length_aa": "841", "iii_isoform_assignment": "IIIb",
    }])
    res, out, report = run_export(tmp_path, selected, cache, "--strict")
    assert res.returncode == 0, res.stderr
    fasta = out.read_text()
    assert fasta.startswith(">fgfr2prot_000001|species=homo_sapiens|source=NCBI|role=reference")
    rows = read_tsv(report)
    assert rows[0]["match_method"] == "ncbi_exact_accession"
    assert rows[0]["match_confidence"] == "high"
    assert rows[0]["length_check_status"] == "exact"


def test_ncbi_product_rescue_is_medium_confidence_and_warned(tmp_path):
    cache = make_ncbi_cache(tmp_path, protein_acc="XP_999.1", tx="XM_OTHER.1", seq_len=840)
    selected = make_selected(tmp_path, [{
        "species_input": "Homo sapiens", "species_canonical": "homo_sapiens", "source_db": "NCBI",
        "selection_role": "FGFR2_IIIc_candidate", "internal_transcript_id": "tx2", "transcript_id_source": "XM_001.1",
        "translation_id_source": "", "protein_length_aa": "841", "iii_isoform_assignment": "IIIc",
    }])
    res, out, report = run_export(tmp_path, selected, cache)
    assert res.returncode == 0, res.stderr
    rows = read_tsv(report)
    assert rows[0]["match_method"] == "ncbi_product_species_length_rescue"
    assert rows[0]["match_confidence"] == "medium"
    warnings = read_tsv(tmp_path / "protein_export_warnings.tsv")
    assert any(w["warning_code"] == "ncbi_rescue_match_used" for w in warnings)


def test_ensembl_export_can_be_mocked_without_live_rest(tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    cache.mkdir()
    selected = make_selected(tmp_path, [{
        "species_input": "Homo sapiens", "species_canonical": "homo_sapiens", "source_db": "Ensembl",
        "selection_role": "reference", "internal_transcript_id": "tx3", "transcript_id_source": "ENST1",
        "translation_id_source": "ENSP1", "protein_length_aa": "600", "iii_isoform_assignment": "unclassified",
    }])
    class Args: pass
    args = Args()
    args.selected = selected; args.cache = cache; args.out = tmp_path / "out.faa"; args.report = tmp_path / "report.tsv"
    args.warnings = tmp_path / "warnings.tsv"; args.metadata = tmp_path / "meta.json"; args.md_report = tmp_path / "report.md"; args.html_report = tmp_path / "report.html"
    args.region_qc = tmp_path / "region_qc.tsv"; args.iii_region_start = 250; args.iii_region_end = 430
    args.roles = mod.DEFAULT_ROLES; args.strict = True; args.no_ensembl_rest = False; args.disable_ncbi_rescue = False
    args.ensembl_sleep = 0; args.ensembl_timeout = 1; args.min_protein_len = 500; args.max_protein_len = 1200
    monkeypatch.setattr(mod, "fetch_ensembl_protein", lambda translation_id, sleep=0, timeout=1: (make_seq(600), "mocked"))
    report_rows, warnings, metadata = mod.export_selected(args)
    assert report_rows[0]["match_method"] == "ensembl_rest_translation"
    assert report_rows[0]["match_confidence"] == "high"
    assert metadata["fasta_records_written"] == 1
    assert args.out.read_text().startswith(">fgfr2prot_000001")
    # v2.1 additionally writes the fixed-window III-region QC table.
    region_rows = read_tsv(args.region_qc)
    assert [r["output_id"] for r in region_rows] == ["fgfr2prot_000001"]
    assert region_rows[0]["iii_region_window_start_1based"] == "250"
    assert region_rows[0]["iii_region_window_end_1based"] == "430"
    assert metadata["region_qc_rows"] == 1


def test_major_length_difference_creates_warning(tmp_path):
    cache = make_ncbi_cache(tmp_path, seq_len=600)
    selected = make_selected(tmp_path, [{
        "species_input": "Homo sapiens", "species_canonical": "homo_sapiens", "source_db": "NCBI",
        "selection_role": "reference", "internal_transcript_id": "tx1", "transcript_id_source": "XM_001.1",
        "translation_id_source": "", "protein_length_aa": "841", "iii_isoform_assignment": "IIIb",
    }])
    res, out, report = run_export(tmp_path, selected, cache)
    assert res.returncode == 0
    rows = read_tsv(report)
    assert rows[0]["length_check_status"] == "major_difference"
    warnings = read_tsv(tmp_path / "protein_export_warnings.tsv")
    assert any(w["warning_code"] == "major_length_difference" for w in warnings)


def test_strict_mode_fails_on_missing_required_column(tmp_path):
    selected = tmp_path / "bad.tsv"
    write_tsv(selected, [{"species_canonical": "homo_sapiens"}], ["species_canonical"])
    cache = tmp_path / "cache"; cache.mkdir()
    res, out, report = run_export(tmp_path, selected, cache, "--strict")
    assert res.returncode != 0
    assert "Missing required columns" in res.stderr


def test_reports_metadata_and_warnings_files_are_written(tmp_path):
    cache = make_ncbi_cache(tmp_path)
    selected = make_selected(tmp_path, [{
        "species_input": "Homo sapiens", "species_canonical": "homo_sapiens", "source_db": "NCBI",
        "selection_role": "reference", "internal_transcript_id": "tx1", "transcript_id_source": "XM_001.1",
        "translation_id_source": "", "protein_length_aa": "841", "iii_isoform_assignment": "IIIb",
    }])
    res, out, report = run_export(tmp_path, selected, cache)
    assert res.returncode == 0, res.stderr
    assert (tmp_path / "protein_export_warnings.tsv").exists()
    assert (tmp_path / "run_metadata.json").exists()
    assert (tmp_path / "protein_export_report.md").exists()
    assert (tmp_path / "protein_export_report.html").exists()
    meta = json.loads((tmp_path / "run_metadata.json").read_text())
    assert meta["script_name"] == "export_selected_fgfr2_proteins_complete_v2.py"
    assert meta["fasta_records_written"] == 1
