#!/usr/bin/env python3
"""One documented command must rebuild every species from the returned cluster outputs.

The user-facing contract is a single command:

    .venv/bin/edc cluster roundtrip --run-id <run_id>

whose post-cluster half is ``run_core_gene_analysis.py --post``. This test starts from a
run that has real InterProScan and pyTMHMM outputs but no derived tables, runs that one
command, and checks that everything a two-species dataset needs comes back — for *both*
species, plus the comparative layer.

It is an integration test and takes about a minute and a half, because it really does
regenerate the publication figures. That cost is the point: the defect it guards against
was a rebuild that silently produced results for one species only, and a mocked test
would not have caught it.
"""
from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

REFERENCE = ROOT / "runs" / "2026-07-26_2157_fgfr1_gallus_mus_core_pilot"
SANDBOX_ID = "test_roundtrip_integration"

# Derived artifacts that the post phase must reconstruct. Deleting them is what makes
# this a rebuild test rather than a check that the reference run happens to be complete.
DERIVED_TABLES = (
    "interpro_annotations.tsv",
    "domain_features.tsv",
    "tm_features.tsv",
    "exon_domain_boundary_distances.tsv",
)

pytestmark = pytest.mark.skipif(
    not (REFERENCE / "results/14_interproscan/primary/output/input.fasta.tsv").is_file(),
    reason="reference two-species run with returned cluster outputs not present")


def _tsv(path: Path):
    if not path.is_file():
        return []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def rebuilt(tmp_path_factory):
    """Run the documented post-cluster command over a stripped copy of the run.

    The copy lives under a disposable configured data root.  The repository reference
    is read-only input and neither the repository runs root nor a user's real data root
    is ever a write target.
    """
    data_root = tmp_path_factory.mktemp("edc-roundtrip-data")
    sandbox = data_root / "runs" / SANDBOX_ID
    sandbox.parent.mkdir(parents=True)
    shutil.copytree(REFERENCE, sandbox)

    # The Phase C adapter deliberately rejects a copied directory whose durable
    # identity still names the source run.  Keep the source fixture untouched and
    # make the identity/path projection self-consistent only inside this disposable
    # copy before invoking the production post-cluster entry point.
    source_id = REFERENCE.name
    for relative in ("run_config.json", "status.json"):
        path = sandbox / relative
        path.write_text(
            path.read_text(encoding="utf-8").replace(source_id, SANDBOX_ID),
            encoding="utf-8",
        )

    core = sandbox / "results" / "core_gene_analysis"
    for name in DERIVED_TABLES:
        (core / name).unlink(missing_ok=True)
    (sandbox / "website_indices" / "generic" / "protein_coordinate_model.json").unlink(
        missing_ok=True)

    env = os.environ.copy()
    env["EDC_DATA_DIR"] = str(data_root)
    env["MPLCONFIGDIR"] = str(data_root / "matplotlib")
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts/framework/run_core_gene_analysis.py"),
         "--post", "--run-id", SANDBOX_ID],
        cwd=str(ROOT), capture_output=True, text=True, env=env)
    try:
        yield sandbox, proc
    finally:
        shutil.rmtree(sandbox, ignore_errors=True)


def test_the_documented_command_succeeds(rebuilt):
    sandbox, proc = rebuilt
    assert proc.returncode == 0, (
        f"post phase exited {proc.returncode}\n"
        f"--- stdout tail ---\n{proc.stdout[-2500:]}\n"
        f"--- stderr tail ---\n{proc.stderr[-2500:]}")


def test_no_undocumented_second_command_is_needed(rebuilt):
    """Everything below is asserted after exactly one command has run."""
    sandbox, _ = rebuilt
    model = sandbox / "website_indices" / "generic" / "protein_coordinate_model.json"
    assert model.is_file(), (
        "the coordinate model must be rebuilt by the documented command, not by a "
        "separate manual step")


# --------------------------------------------------------------------------- #
# every species is rebuilt through every stage
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("table", DERIVED_TABLES)
def test_every_rebuilt_table_covers_both_species(rebuilt, table):
    sandbox, _ = rebuilt
    rows = _tsv(sandbox / "results" / "core_gene_analysis" / table)
    assert rows, f"{table} was not rebuilt"
    species = {r.get("species_id") for r in rows}
    assert species == {"gallus_gallus", "mus_musculus"}, (
        f"{table} covers {species}; a rebuild that produces one species is the defect "
        f"this test exists for")


def test_boundaries_are_rebuilt_for_each_species_separately(rebuilt):
    sandbox, _ = rebuilt
    rows = _tsv(sandbox / "results/core_gene_analysis/exon_domain_boundary_distances.tsv")
    per_species = {}
    for r in rows:
        per_species.setdefault(r.get("species_id"), 0)
        per_species[r["species_id"]] += 1
    assert set(per_species) == {"gallus_gallus", "mus_musculus"}
    # Counts may legitimately differ between species; only presence is contractual.
    assert all(n > 0 for n in per_species.values()), per_species


def test_no_species_overwrote_another(rebuilt):
    """Long-format tables must keep both species' rows, not one species twice."""
    sandbox, _ = rebuilt
    rows = _tsv(sandbox / "results/core_gene_analysis/domain_features.tsv")
    by_species = {}
    for r in rows:
        by_species.setdefault(r["species_id"], set()).add(r.get("protein_id"))
    assert len(by_species) == 2
    proteins = [p for s in by_species.values() for p in s]
    assert len(set(proteins)) == len(by_species), (
        f"each species must contribute its own protein; got {by_species}")


# --------------------------------------------------------------------------- #
# returned-sequence inventory
# --------------------------------------------------------------------------- #
def test_the_returned_sequence_inventory_reports_every_species(rebuilt):
    sandbox, _ = rebuilt
    qc = json.loads(
        (sandbox / "results/15_domain_architecture/post_cluster_qc.json").read_text())
    inv = qc.get("returned_sequence_inventory")
    assert inv, "the post phase must record which species' sequences came back"
    assert inv["n_species_submitted"] == 2
    assert inv["n_species_with_features"] == 2
    for row in inv["species"]:
        assert row["interproscan_status"] == "available", row
        assert row["protein_id"], row
        assert row["resolved_from"], "the inventory must record how the primary resolved"


def test_the_inventory_distinguishes_no_features_from_no_sequence():
    """A protein with no predicted TM region is not a missing analysis."""
    from framework.run_core_gene_analysis import _returned_sequence_inventory
    core = REFERENCE / "results" / "core_gene_analysis"
    inv = _returned_sequence_inventory(
        core, {"rows": [{"protein_id": "NP_990841.2"}]}, {"rows": []})
    statuses = {r["species_id"]: r for r in inv["species"]}
    assert statuses["gallus_gallus"]["interproscan_status"] == "available"
    assert statuses["gallus_gallus"]["pytmhmm_status"] == "no_transmembrane_predicted"
    # Mus had no returned rows in this synthetic call, but rows existed overall, so it
    # must read as "returned_no_features", never as a silently dropped species.
    assert statuses["mus_musculus"]["interproscan_status"] == "returned_no_features"


# --------------------------------------------------------------------------- #
# per-species coordinate models and comparative layer
# --------------------------------------------------------------------------- #
def test_both_species_get_a_complete_coordinate_model(rebuilt):
    sandbox, _ = rebuilt
    model = json.loads(
        (sandbox / "website_indices/generic/protein_coordinate_model.json").read_text())
    by_species = {m["species_id"]: m for m in model["models"]}
    assert set(by_species) == {"gallus_gallus", "mus_musculus"}
    for sid, m in by_species.items():
        assert m["status"] == "available", f"{sid}: {m['status']}"
        assert m["representative_domains"], f"{sid}: no domains"
        assert m["families_superfamilies"], f"{sid}: no family layer"
        assert m["member_signatures"], f"{sid}: no member signatures"
        assert m["tm_regions"], f"{sid}: no TM region"
        assert m["exon_boundaries"], f"{sid}: no boundaries"


def test_the_comparative_mapping_and_index_are_rebuilt(rebuilt):
    sandbox, _ = rebuilt
    model = json.loads(
        (sandbox / "website_indices/generic/protein_coordinate_model.json").read_text())

    assert model["msa_coordinate_map"]["available"] is True
    for sid, rep in model["msa_boundary_mapping"].items():
        assert rep["boundaries_mapped"] == rep["boundaries_total"] > 0, f"{sid}: {rep}"

    dash = model["boundary_dashboard"]
    assert dash["page_mode"] == "generic_multi_species_results_ready"
    ms = dash["multi_species"]
    assert ms["available"] is True
    assert ms["comparable_boundary_groups"], "no comparable groups after rebuild"
    assert len(ms["boundary_matrix"]) == 2
    assert len(ms["distance_statistics"]) == len(ms["comparable_boundary_groups"])


def test_publication_figures_are_generated_for_both_species(rebuilt):
    sandbox, _ = rebuilt
    figures = sandbox / "results/generic_gene_analysis/figures/main"
    if not figures.is_dir():
        pytest.skip("figure rendering unavailable in this environment")
    names = {p.name for p in figures.glob("*.svg")}
    for sid in ("gallus_gallus", "mus_musculus"):
        assert any(sid in n for n in names), (
            f"no main figure for {sid}; figures were generated for "
            f"{sorted({n.split('_')[1] for n in names})}")


# --------------------------------------------------------------------------- #
# status aggregation
# --------------------------------------------------------------------------- #
def test_the_run_status_is_ready_only_because_both_species_are_complete(rebuilt):
    sandbox, _ = rebuilt
    st = json.loads((sandbox / "status.json").read_text())
    assert st["status"] == "results_ready", st.get("failed_reason")
    species_status = st.get("species_status") or {}
    assert set(species_status) == {"gallus_gallus", "mus_musculus"}
    for sid, rec in species_status.items():
        assert rec["complete"] is True, f"{sid}: {rec}"
        assert rec["domain_architecture"] == "available"
        assert rec["boundary"] == "available"


def test_the_roundtrip_accepts_a_partial_run_as_an_end_state():
    """``post_cluster_partial`` is a scientific outcome, not a roundtrip failure."""
    src = (ROOT / "scripts/interpro_cluster/run_cluster_roundtrip.py").read_text()
    assert "post_cluster_partial" in src, (
        "finalize() must accept the partial status the per-species model produces, "
        "otherwise the roundtrip aborts on exactly the runs it is meant to report")


# --------------------------------------------------------------------------- #
# no silent first-species fallback anywhere in the rebuild
# --------------------------------------------------------------------------- #
def test_an_unattributable_result_row_stops_the_phase():
    src = (ROOT / "src/exondomaincompare/framework/run_core_gene_analysis.py").read_text()
    assert "default_species" not in src, (
        "attributing an unresolvable protein to the first species produces a "
        "complete-looking table in which one species carries another's domains")
    assert "Refusing to assign them to an" in src


def test_the_primary_protein_list_covers_every_species():
    from framework.run_core_gene_analysis import _primary_protein_ids
    ids = _primary_protein_ids(REFERENCE / "results" / "core_gene_analysis")
    assert len(ids) == 2, f"expected one primary per species, got {ids}"
    assert len(set(ids)) == 2
