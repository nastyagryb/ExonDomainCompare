import csv
import subprocess
import sys

from script_paths import load_script_module, script_path

SCRIPT_NAME = "build_species_registry_improved.py"
SCRIPT = script_path(SCRIPT_NAME)
mod = load_script_module(SCRIPT_NAME, "build_species_registry_improved")


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_slug_species_normalizes_spaces_and_punctuation():
    assert mod.slug_species("  Homo   sapiens ") == "homo_sapiens"
    assert mod.slug_species("Species-name test") == "species_name_test"


def test_known_species_row_contains_taxid_and_status():
    """A cached species needs no lookup, so this holds without network access."""
    result = mod.build_registry_rows(["Homo sapiens"], "ensembl_first", "RefSeq",
                                     offline=True)
    assert result.rows[0]["taxid"] == "9606"
    assert result.rows[0]["ensembl_species"] == "homo_sapiens"
    assert result.rows[0]["status"] == "taxon_verified"
    assert result.warnings == []


def test_unresolvable_species_is_kept_but_carries_no_invented_identity():
    """The registry must not echo the submitted slug back as a scientific name.

    Doing so is what sent ``equus_quagga`` to NCBI Datasets as if it were a taxon
    name. An unresolved species keeps its row, but with empty identity fields.
    """
    result = mod.build_registry_rows(["Example vertebrata"], "ensembl_first", "RefSeq",
                                     offline=True)
    row = result.rows[0]
    assert row["taxid"] == ""
    assert row["ncbi_species"] == ""
    assert row["status"] == "taxon_unverified_offline"
    assert result.warnings[0]["warning_code"] == "TAXON_UNVERIFIED_OFFLINE"


def test_duplicate_species_is_deduplicated_and_warned():
    result = mod.build_registry_rows(["Homo sapiens", "homo sapiens", "Homo sapiens"], "ensembl_first", "RefSeq", offline=True)
    assert len(result.rows) == 1
    warning_codes = [warning["warning_code"] for warning in result.warnings]
    assert "DUPLICATE_SPECIES" in warning_codes


def test_cli_writes_expected_files(tmp_path):
    species_list = tmp_path / "species.txt"
    species_list.write_text("# comment\nHomo sapiens\nUnknown animal\n", encoding="utf-8")
    outdir = tmp_path / "out"

    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--species_list", str(species_list),
         "--outdir", str(outdir), "--offline"],
        check=True,
        text=True,
        capture_output=True,
    )
    assert "[OK] species_registry.tsv rows: 2" in completed.stdout

    registry = read_tsv(outdir / "species_registry.tsv")
    warnings = read_tsv(outdir / "species_registry_warnings.tsv")
    assert [row["species_id"] for row in registry] == ["homo_sapiens", "unknown_animal"]
    assert warnings[0]["warning_code"] == "TAXON_UNVERIFIED_OFFLINE"


def test_strict_mode_fails_on_unknown_species(tmp_path):
    species_list = tmp_path / "species.txt"
    species_list.write_text("Unknown animal\n", encoding="utf-8")
    outdir = tmp_path / "out"

    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--species_list", str(species_list),
         "--outdir", str(outdir), "--strict", "--offline"],
        text=True,
        capture_output=True,
    )
    assert completed.returncode != 0
    assert "Strict mode failed" in completed.stderr
