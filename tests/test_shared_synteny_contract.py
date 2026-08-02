"""The shared local-synteny contract, checked against real run data.

The regression these tests exist for: the FGFR1 Mus musculus neighbourhood
displayed five upstream loci, the target and only four downstream loci. All ten
neighbours were present in ``synteny_neighbors.tsv`` and in the index; the tenth
was lost in the browser, where the loci were laid out as a horizontally
scrolling row of fixed-width buttons and the outermost one fell outside the
visible area. So the assertions here run in both directions: every real locus
reaches the contract, and no locus that is not in the source table appears in it.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

RUNS = {
    "fgfr1_two_species": "2026-07-26_2157_fgfr1_gallus_mus_core_pilot",
    "fgfr1_gallus": "2026-07-23_1100_fgfr1_gallus_core_pilot",
    "tp53_danio": "2026-07-21_1436_custom_run",
    "tpm1_two_species": "2026-07-16_1638_tpm1_human_mouse_twospecies",
}
FREEZE_CLOSURE = (ROOT / "results" / "final_30_until_interpro_prepare"
                  / "13_final_pre_interpro_closure")


def _run_dir(key: str) -> Path:
    path = ROOT / "runs" / RUNS[key]
    if not path.is_dir():
        pytest.skip(f"reference run {RUNS[key]} is not present")
    return path


def _index(key: str) -> dict:
    from shared_gene_analysis.common import SharedRunContext
    from shared_gene_analysis.indices.synteny_locus import build_synteny_locus_index
    return build_synteny_locus_index(SharedRunContext.from_run_dir(_run_dir(key)))


def _row(index: dict, species_id: str) -> dict:
    row = next((r for r in index["species"] if r["species_id"] == species_id), None)
    assert row is not None, f"{species_id} missing from the synteny index"
    return row


def _source_counts(run_dir: Path, species_id: str) -> Counter:
    path = run_dir / "results" / "core_gene_analysis" / "synteny_neighbors.tsv"
    with path.open(newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r.get("species_id") == species_id]
    return Counter(r["side"] for r in rows)


# --------------------------------------------------------------------------- #
# The reported case: FGFR1 Mus musculus
# --------------------------------------------------------------------------- #
def test_fgfr1_mus_keeps_all_five_downstream_neighbours():
    """The missing fifth downstream locus was a rendering loss, not missing data."""
    run_dir = _run_dir("fgfr1_two_species")
    source = _source_counts(run_dir, "mus_musculus")
    assert source["upstream"] == 5 and source["downstream"] == 5

    row = _row(_index("fgfr1_two_species"), "mus_musculus")
    assert row["displayed_upstream_count"] == 5
    assert row["displayed_downstream_count"] == 5
    assert row["truncation_status"] == "complete"
    assert row["counts_label"] == "10 flanking loci shown · 5 upstream · 5 downstream"
    assert [n["symbol"] for n in row["downstream"]] == [
        "Letm2", "Nsd3", "Plpp5", "Ddhd2", "Bag4"]


def test_fgfr1_mus_target_is_explicit_and_centred():
    row = _row(_index("fgfr1_two_species"), "mus_musculus")
    assert row["target_slot"] == 0
    assert row["target_symbol"] == "FGFR1"
    assert row["target_gene_id"] == "gene-Fgfr1"
    assert row["target"]["is_target"] is True
    # Exactly one central slot, and it sits in the middle of the display order.
    slots = [n["slot_x"] for n in row["loci"]]
    assert slots == sorted(slots)
    assert slots.count(0) == 1
    assert row["loci"][len(row["loci"]) // 2]["is_target"] is True


def test_fgfr1_gallus_matches_its_source_table():
    run_dir = _run_dir("fgfr1_two_species")
    source = _source_counts(run_dir, "gallus_gallus")
    row = _row(_index("fgfr1_two_species"), "gallus_gallus")
    assert row["displayed_upstream_count"] == source["upstream"]
    assert row["displayed_downstream_count"] == source["downstream"]
    assert [n["symbol"] for n in row["upstream"]] == [
        "ADAM9", "TM2D2", "PLEKHA2", "LOC121107413", "TACC1"]
    # An NCBI LOC identifier is a real locus with no curated symbol, not a
    # failed assignment.
    loc = next(n for n in row["upstream"] if n["symbol"] == "LOC121107413")
    assert loc["placeholder"] is True
    assert loc["orthology_class"] == "placeholder"


# --------------------------------------------------------------------------- #
# The same contract for other genes
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("key,species", [
    ("fgfr1_gallus", "gallus_gallus"),
    ("tp53_danio", "danio_rerio"),
    ("tpm1_two_species", "homo_sapiens"),
    ("tpm1_two_species", "mus_musculus"),
])
def test_rendered_locus_count_equals_source_table_count(key, species):
    run_dir = _run_dir(key)
    source = _source_counts(run_dir, species)
    row = _row(_index(key), species)
    assert row["displayed_upstream_count"] == min(source["upstream"], 5)
    assert row["displayed_downstream_count"] == min(source["downstream"], 5)
    assert row["upstream_count_available"] == source["upstream"]
    assert row["downstream_count_available"] == source["downstream"]


@pytest.mark.parametrize("key", list(RUNS))
def test_target_is_never_counted_as_a_neighbour(key):
    for row in _index(key)["species"]:
        flanking = row["upstream"] + row["downstream"]
        assert all(n["is_target"] is False for n in flanking)
        assert len(flanking) == row["displayed_flanking_count"]
        assert row["displayed_flanking_count"] == (
            row["displayed_upstream_count"] + row["displayed_downstream_count"])
        # The target is in `loci` exactly once, and only there.
        assert sum(1 for n in row["loci"] if n["is_target"]) == 1


@pytest.mark.parametrize("key", list(RUNS))
def test_no_fabricated_loci(key):
    """Every displayed symbol traces back to a row of the source table."""
    run_dir = _run_dir(key)
    path = run_dir / "results" / "core_gene_analysis" / "synteny_neighbors.tsv"
    with path.open(newline="", encoding="utf-8") as fh:
        source = {(r["species_id"], r["side"], int(r["order"]), r["neighbor_symbol"])
                  for r in csv.DictReader(fh, delimiter="\t")}
    for row in _index(key)["species"]:
        for locus in row["upstream"] + row["downstream"]:
            key_tuple = (row["species_id"], locus["side"], locus["rank"],
                         locus["source_symbol"])
            assert key_tuple in source, f"{key_tuple} is not in the source table"


@pytest.mark.parametrize("key", list(RUNS))
def test_complete_five_and_five_case_is_reported_as_complete(key):
    rows = _index(key)["species"]
    full = [r for r in rows if r["upstream_count_available"] == 5
            and r["downstream_count_available"] == 5]
    assert full, "no species in this run has the full five-and-five neighbourhood"
    for row in full:
        assert row["truncation_status"] == "complete"
        assert row["omission_reason"] == ""
        assert row["displayed_flanking_count"] == 10


def test_unequal_neighbour_counts_are_stated_rather_than_padded():
    """A genome with fewer real loci on one side must say so, not invent a gene."""
    from shared_gene_analysis import synteny_contract as sc
    neighbours = (
        [sc.neighbour_locus(side="upstream", rank=i, source_symbol=f"UP{i}",
                            resolved_symbol=f"UP{i}", strand="+") for i in range(1, 6)]
        + [sc.neighbour_locus(side="downstream", rank=i, source_symbol=f"DN{i}",
                              resolved_symbol=f"DN{i}", strand="-") for i in range(1, 5)]
    )
    row = sc.species_row("test_species", gene_symbol="GENEX",
                         target=sc.target_locus(gene_symbol="GENEX", strand="+"),
                         neighbours=neighbours)
    assert row["displayed_upstream_count"] == 5
    assert row["displayed_downstream_count"] == 4
    assert row["displayed_flanking_count"] == 9
    assert row["truncation_status"] == "fewer_available"
    assert row["counts_label"] == "9 flanking loci shown · 5 upstream · 4 downstream"
    assert "only 5 upstream and 4 downstream" in row["omission_reason"]
    assert len(row["loci"]) == 10  # nine real loci plus the target, nothing padded


def test_upstream_reads_outward_to_inward_so_the_display_order_is_genomic():
    from shared_gene_analysis import synteny_contract as sc
    neighbours = [
        sc.neighbour_locus(side="upstream", rank=i, source_symbol=f"UP{i}")
        for i in range(1, 4)
    ] + [
        sc.neighbour_locus(side="downstream", rank=i, source_symbol=f"DN{i}")
        for i in range(1, 4)
    ]
    row = sc.species_row("s", gene_symbol="G",
                         target=sc.target_locus(gene_symbol="G"),
                         neighbours=neighbours)
    assert [n["symbol"] for n in row["loci"]] == [
        "UP3", "UP2", "UP1", "G", "DN1", "DN2", "DN3"]


def test_status_and_orthology_strings_carry_readable_labels():
    """Internal state strings never reach the interface as the primary label."""
    from shared_gene_analysis import synteny_contract as sc
    for status in sc.STATUS_DISPLAY:
        label, definition = sc.status_display(status)
        assert "_" not in label
        assert definition
    for cls in sc.ORTHOLOGY_DISPLAY:
        label, definition = sc.orthology_display(cls)
        assert "_" not in label
        assert definition
    assert sc.orthology_display("placeholder")[1].startswith(
        "Placeholder locus label; curated gene symbol unavailable.")


# --------------------------------------------------------------------------- #
# FGFR2 uses the same contract
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not FREEZE_CLOSURE.is_dir(), reason="FGFR2 freeze not present")
def test_fgfr2_uses_the_same_contract_with_curated_symbols():
    import build_website_indices as bwi
    index = bwi.build_synteny_locus_index(FREEZE_CLOSURE)
    assert index["contract"] == "shared_synteny_v1"
    assert len(index["species"]) == 30
    for row in index["species"]:
        assert row["target_symbol"] == "FGFR2"
        assert row["target_slot"] == 0
        assert row["displayed_upstream_count"] <= 5
        assert row["displayed_downstream_count"] <= 5
        assert all(not n["is_target"] for n in row["upstream"] + row["downstream"])
    # The curated 5-neighbour panel must reach the displayed loci; before the
    # repair the view fell back to the raw supplement and showed bare identifiers.
    human = _row(index, "homo_sapiens")
    assert [n["symbol"] for n in human["downstream"]][:3] == [
        "WDR11", "PLPP4", "SEC23IP"]
    assert human["target"]["gene_id"] == "GeneID:2263"


# --------------------------------------------------------------------------- #
# The browser-side layout, exercised in Node against real indices
# --------------------------------------------------------------------------- #
HARNESS = Path(__file__).with_name("check_synteny_renderer.mjs")


@pytest.mark.parametrize("key", list(RUNS))
def test_renderer_centres_the_target_and_clips_nothing(key, tmp_path):
    """The regression lived in the browser, so the layout is checked in Node."""
    import shutil
    if shutil.which("node") is None:
        pytest.skip("node is required for the renderer checks")
    index_path = tmp_path / "synteny_locus_index.json"
    index_path.write_text(json.dumps(_index(key)), encoding="utf-8")
    proc = subprocess.run([shutil.which("node"), str(HARNESS), str(index_path)],
                          capture_output=True, text=True, cwd=str(ROOT))
    assert proc.returncode == 0, f"{proc.stdout}\n{proc.stderr}"
    assert "FAIL" not in proc.stdout


@pytest.mark.skipif(not FREEZE_CLOSURE.is_dir(), reason="FGFR2 freeze not present")
def test_renderer_handles_the_thirty_species_fgfr2_dataset(tmp_path):
    import shutil
    if shutil.which("node") is None:
        pytest.skip("node is required for the renderer checks")
    import build_website_indices as bwi
    index_path = tmp_path / "synteny_locus_index.json"
    index_path.write_text(json.dumps(bwi.build_synteny_locus_index(FREEZE_CLOSURE)),
                          encoding="utf-8")
    proc = subprocess.run([shutil.which("node"), str(HARNESS), str(index_path)],
                          capture_output=True, text=True, cwd=str(ROOT))
    assert proc.returncode == 0, f"{proc.stdout}\n{proc.stderr}"


def test_the_frontend_has_exactly_one_synteny_renderer():
    """A second drawing implementation is free to disagree with this one."""
    src = ROOT / "webapp" / "frontend" / "src"
    drawing = [p for p in src.rglob("*.jsx")
               if "st-genes" in p.read_text(encoding="utf-8")]
    assert drawing == [], f"a second synteny track renderer survives in {drawing}"
    viewer = (src / "pages" / "viewers" / "SyntenyViewer.jsx").read_text(encoding="utf-8")
    assert "SyntenyNeighbourhood" in viewer


@pytest.mark.skipif(not FREEZE_CLOSURE.is_dir(), reason="FGFR2 freeze not present")
def test_rebuild_command_leaves_the_freeze_untouched(tmp_path):
    out = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "fgfr2" / "rebuild_fgfr2_gallery.py"),
         "--dataset", "example", "--derived-root", str(tmp_path / "derived")],
        capture_output=True, text=True, cwd=ROOT,
        env={"PYTHONPATH": str(ROOT / "scripts"), "PATH": "/usr/bin:/bin"})
    assert out.returncode == 0, out.stderr
    payload = json.loads(out.stdout)[0]
    assert payload["freeze_unchanged"] is True
    assert payload["freeze_changed_files"] == []
    assert payload["output"] == "<DERIVED_ROOT>/example/website_indices"
