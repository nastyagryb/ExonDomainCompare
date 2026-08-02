"""A newly created FGFR2 run must reach the same state and the same Gallery as the freeze.

Three defects met in one run — ``2026-07-29_1634_fgfr2_homo_sapiens_felis_catus`` — and each
one is pinned here:

1. Readiness was computed circularly. The post-cluster views were added to the required list
   only when the cluster outputs they check were *already present*, so a run whose round-trip
   had not happened was measured against the pre-cluster views alone and published as
   "Results ready" while its own page asked for the round-trip.

2. The cluster round-trip could never finish an FGFR2 run. ``post_interpro()`` called
   ``finalize()`` outside the core-only branch and returned, leaving the FGFR2 post-cluster
   builder below it unreachable: the annotation was fetched and then judged without ever
   being turned into domain and boundary layers.

3. The modern Gallery was wired into a migration script and the validated freeze's own paths,
   not into the path a new run takes, so every new FGFR2 run fell back to the legacy per-file
   catalogue.

The states a card may report are also pinned, because "pending cluster", "not applicable for
one species", "no comparable evidence" and "output missing" are four different findings and
collapsing them is what made the Gallery unreadable.
"""
from __future__ import annotations

import ast
import json
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "webapp" / "backend"))

from fgfr2 import coordinate_model as cm  # noqa: E402
from fgfr2 import gallery_catalogue as gc  # noqa: E402
from fgfr2 import run_gallery as rg  # noqa: E402
from shared_gene_analysis import run_availability as ra  # noqa: E402

REAL_RUN = ROOT / "runs" / "2026-07-29_1634_fgfr2_homo_sapiens_felis_catus"
FREEZE_INDICES = ROOT / "results" / "derived" / "example" / "website_indices"
FREEZE = ROOT / "results" / "final_30_until_interpro_prepare"

#: The two species the affected run was created for.
SPECIES = (("homo_sapiens", "IIIb", "ENSP1"), ("homo_sapiens", "IIIc", "ENSP2"),
           ("felis_catus", "IIIb", "XP1"), ("felis_catus", "IIIc", "XP2"))


# --------------------------------------------------------------------------- #
# fixture: a two-species FGFR2 run stopped exactly where the real one stopped
# --------------------------------------------------------------------------- #
def _truth_table() -> str:
    head = ("species\tisoform\tprotein_id\ttranscript_id\tfinal_isoform_label\t"
            "protein_length\tselection_role\n")
    rows = "".join(f"{sp}\t{iso}\t{pid}\ttx_{pid}\t{iso}\t800\tprimary\n"
                   for sp, iso, pid in SPECIES)
    return head + rows


CLOSURE_INPUTS = {
    "final_pre_interpro_truth_table.tsv": _truth_table(),
    "tables/figure3C_exon_to_protein_cassette_coordinate_map.tsv":
        "species\tisoform\tfinal_isoform_label\tprotein_length\tcassette_start_aa\t"
        "cassette_end_aa\tfeature_type\texon_or_cds_id\tfeature_label\tblock_start_aa\t"
        "block_end_aa\n"
        "homo_sapiens\tIIIb\tIIIb\t800\t383\t428\tCDS\tcds-1\texon 1\t1\t100\n"
        "felis_catus\tIIIb\tIIIb\t800\t383\t428\tCDS\tcds-1\texon 1\t1\t100\n",
    "tables/figure6B_species_resolved_IIIb_IIIc_cassette_residue_map.tsv":
        "species\tisoform\thuman_reference_residue_index\n",
    "tables/figure3B_IIIb_IIIc_cassette_amino_acid_motif_map.tsv":
        "human_reference_residue_index\thuman_IIIb_aa_one_letter\n1\tS\n",
    "MSA/final_fgfr2_full_length_protein_msa.aln.faa": "".join(
        f">{sp}|{iso}|{pid}|main_figure\nMDHTSV\n" for sp, iso, pid in SPECIES),
}

INDEX_PAYLOADS = {
    "run_index.json": {"run_id": "r", "gate_status": "pass", "kpi": {"species": 2}},
    "species_index.json": [{"species": "homo_sapiens"}, {"species": "felis_catus"}],
    "coordinate_track_index.json": {"available": True, "species": [{"species": "hs"}]},
    "cassette_residue_index.json": {"available": True, "species": [{"species": "hs"}],
                                    "evidence_level": "sequence_marker"},
    "msa_index.json": {"available": True,
                       "alignments": {"full_length": {"available": True, "rows": [{}]}}},
    "synteny_locus_index.json": {"available": True, "species": [{"species": "hs"}]},
    "figure_index.json": {"figures": [{"key": "f1"}]},
    "download_index.json": [{"group": "Tables", "label": "Final truth table"}],
}

#: Written only once the cluster annotation is back, because that is what produces them.
POST_CLUSTER_INDEX_PAYLOADS = {
    "species_domain_architecture.json": {"available": True, "species": [{}]},
    "domain_architecture_index.json": {"available": True, "species": [{}]},
    "boundary_consistency_summary.json": {"available": True, "total_primary_proteins": 4},
}


def _make_run(tmp_path: Path, *, with_cluster: bool = False) -> Path:
    """A two-species FGFR2 run with its pre-InterPro closure and primary FASTA."""
    run = tmp_path / "runs" / "2026-07-29_1634_fgfr2_homo_sapiens_felis_catus"
    closure = run / "results" / "13_final_pre_interpro_closure"
    for rel, text in CLOSURE_INPUTS.items():
        target = closure / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")

    freeze = closure / "freeze"
    freeze.mkdir(parents=True, exist_ok=True)
    (freeze / "final_pre_interpro_proteins_primary.faa").write_text(
        "".join(f">{sp}|{iso}|{iso}|{pid}|tx_{pid}|primary\nMDHTSV\n"
                for sp, iso, pid in SPECIES), encoding="utf-8")
    (run / "status.json").write_text(
        json.dumps({"species_count": 2, "status": "pre_interpro_complete"}),
        encoding="utf-8")
    (run / "run_config.json").write_text(
        json.dumps({"gene_symbol": "FGFR2", "has_event": True,
                    "species": ["Homo sapiens", "Felis catus"]}), encoding="utf-8")

    if with_cluster:
        ips = run / "results" / "14_interproscan" / "primary" / "output"
        ips.mkdir(parents=True, exist_ok=True)
        (ips / "input.fasta.tsv").write_text("ENSP1\tPF00069\n", encoding="utf-8")
        tm = (run / "results" / "15_exon_domain_boundary_post_interpro"
              / "pytmhmm_primary" / "output")
        tm.mkdir(parents=True, exist_ok=True)
        (tm / "pytmhmm_transmembrane_hits.tsv").write_text("ENSP1\t1\t20\n",
                                                           encoding="utf-8")

    time.sleep(0.02)
    wi = run / "website_indices"
    wi.mkdir(parents=True, exist_ok=True)
    payloads = dict(INDEX_PAYLOADS)
    if with_cluster:
        payloads.update(POST_CLUSTER_INDEX_PAYLOADS)
    for name, payload in payloads.items():
        (wi / name).write_text(json.dumps(payload), encoding="utf-8")
    return run


@pytest.fixture()
def pre_cluster_run(tmp_path: Path) -> Path:
    return _make_run(tmp_path)


@pytest.fixture()
def post_cluster_run(tmp_path: Path) -> Path:
    return _make_run(tmp_path, with_cluster=True)


# --------------------------------------------------------------------------- #
# 1. the readiness regression
# --------------------------------------------------------------------------- #
def test_post_cluster_views_are_required_whether_or_not_their_outputs_exist():
    """The root cause, pinned at its source.

    Making the requirement conditional on the outputs meant a run could satisfy the
    contract by not having done the work.
    """
    for view in ra.REQUIRED_FOR_READY_WITH_DOMAINS:
        assert view in ra.REQUIRED_FOR_READY + ra.REQUIRED_FOR_READY_WITH_DOMAINS

    source = (ROOT / "src" / "exondomaincompare" / "shared_gene_analysis"
              / "run_availability.py").read_text()
    tree = ast.parse(source)
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "readiness")
    guarded = [n for n in ast.walk(fn) if isinstance(n, ast.If)
               and "REQUIRED_FOR_READY_WITH_DOMAINS" in ast.dump(n)]
    assert not guarded, ("the domain requirements must not sit behind a condition on the "
                         "outputs they are there to check")


def test_a_pre_cluster_two_species_run_is_not_results_ready(pre_cluster_run: Path):
    verdict = ra.readiness(pre_cluster_run, n_species=2, has_event=True)
    assert verdict.ready is False
    blocking = {s.view for s in verdict.views if not s.available}
    assert blocking & {"domain_architecture", "exon_domain_boundaries"}, verdict.reason


def test_missing_cluster_annotation_is_pending_not_not_applicable(pre_cluster_run: Path):
    """A round-trip that has not run is work outstanding, not biology that does not apply."""
    states = {s.view: s for s in ra.view_states(pre_cluster_run, n_species=2)}
    for view in ("domain_architecture", "exon_domain_boundaries"):
        assert states[view].state == ra.PENDING, states[view].reason
        assert states[view].state != "not_applicable"
        assert states[view].available is False


def test_the_primary_fasta_is_available_before_the_cluster(pre_cluster_run: Path):
    fasta = (pre_cluster_run / "results" / "13_final_pre_interpro_closure" / "freeze"
             / "final_pre_interpro_proteins_primary.faa")
    headers = [ln for ln in fasta.read_text().splitlines() if ln.startswith(">")]
    assert len(headers) == 4
    assert len(set(headers)) == 4, "a duplicate header would mean an overwritten model"
    assert {h.split("|")[0].lstrip(">") for h in headers} == {"homo_sapiens", "felis_catus"}
    assert {h.split("|")[1] for h in headers} == {"IIIb", "IIIc"}


def test_a_run_with_real_cluster_annotation_becomes_ready(post_cluster_run: Path):
    verdict = ra.readiness(post_cluster_run, n_species=2, has_event=True)
    assert verdict.ready is True, verdict.reason


@pytest.fixture()
def status_model(pre_cluster_run: Path, monkeypatch: pytest.MonkeyPatch) -> dict:
    """The status model My Runs and the dataset selector both read."""
    import main as backend

    # The backend reports paths relative to the project root, so a run outside it has to
    # be given a matching root rather than a patched reporter.
    root = pre_cluster_run.parents[1]
    monkeypatch.setattr(backend, "PROJECT_ROOT", root)
    monkeypatch.setattr(backend, "LOCAL_RUNS_ROOT", root / "runs")
    return backend.derive_status_model(pre_cluster_run)


def test_my_runs_and_the_dataset_selector_read_one_cluster_required_verdict(
        status_model: dict):
    """Both surfaces project one model, so they cannot disagree by construction."""
    assert status_model["status"] == "cluster_required"
    assert status_model["status_label"] == "Cluster input ready"
    assert status_model["cluster_required"] is True


def test_the_roundtrip_command_is_offered_for_a_cluster_required_run(status_model: dict):
    assert status_model["cluster_command"] == (
        ".venv/bin/python scripts/edc.py cluster roundtrip "
        "--run-id 2026-07-29_1634_fgfr2_homo_sapiens_felis_catus")
    assert status_model["next_action"] == "run_cluster_roundtrip_command"


def test_the_command_is_withheld_when_the_primary_fasta_is_missing(
        pre_cluster_run: Path, monkeypatch: pytest.MonkeyPatch):
    """Part 4's gate: an invalid pre-cluster state must not offer cluster work."""
    import main as backend

    (pre_cluster_run / "results" / "13_final_pre_interpro_closure" / "freeze"
     / "final_pre_interpro_proteins_primary.faa").unlink()
    root = pre_cluster_run.parents[1]
    monkeypatch.setattr(backend, "PROJECT_ROOT", root)
    monkeypatch.setattr(backend, "LOCAL_RUNS_ROOT", root / "runs")

    model = backend.derive_status_model(pre_cluster_run)
    assert model["primary_fasta_status"] == "not_available"
    assert model["cluster_command"] == ""


def test_domain_architecture_reads_as_pending_before_the_roundtrip(
        pre_cluster_run: Path, monkeypatch: pytest.MonkeyPatch):
    """Part 12: the view must say the round-trip is outstanding, not that nothing was found."""
    import main as backend

    root = pre_cluster_run.parents[1]
    monkeypatch.setattr(backend, "PROJECT_ROOT", root)
    monkeypatch.setattr(backend, "LOCAL_RUNS_ROOT", root / "runs")

    state = backend.derive_status_model(pre_cluster_run)["view_availability"]
    assert state["domain_architecture"]["state"] == ra.PENDING
    assert state["exon_domain_boundaries"]["state"] == ra.PENDING


# --------------------------------------------------------------------------- #
# 2. the round-trip could not finish an FGFR2 run
# --------------------------------------------------------------------------- #
def test_the_roundtrip_reaches_the_fgfr2_post_cluster_builder():
    """No statement in ``post_interpro`` may sit after an unconditional return."""
    source = (ROOT / "scripts" / "interpro_cluster" / "run_cluster_roundtrip.py").read_text()
    fn = next(n for n in ast.walk(ast.parse(source))
              if isinstance(n, ast.FunctionDef) and n.name == "post_interpro")
    for index, node in enumerate(fn.body):
        if isinstance(node, ast.Return):
            assert index == len(fn.body) - 1, (
                "statements after this return are unreachable, which is how the FGFR2 "
                "post-cluster builder came to be skipped")
    assert "POST" in ast.dump(fn), "the FGFR2 post-cluster builder must be reachable"


# --------------------------------------------------------------------------- #
# 3. the Gallery a new run gets
# --------------------------------------------------------------------------- #
def _catalogue(run: Path, cluster_ready: bool) -> dict:
    index = cm.build_index(sources=cm.run_sources(run))
    return rg.catalogue_for(run, index, cluster_ready)


def test_a_new_two_species_run_gets_the_modern_gallery(pre_cluster_run: Path):
    catalogue = _catalogue(pre_cluster_run, cluster_ready=False)
    assert catalogue["default_scope"] == "comparative"
    assert set(catalogue["species_scopes"]) == {"homo_sapiens", "felis_catus"}
    assert catalogue["dataset"].startswith("run:")


def test_the_new_run_has_no_legacy_framework_catalogue(pre_cluster_run: Path):
    index = gc.flatten_for_gallery(_catalogue(pre_cluster_run, cluster_ready=False))
    categories = {f.get("category") for f in index["figures"]}
    assert "Framework" not in categories
    assert not any(str(f.get("title", "")).startswith("Figure ")
                   for f in index["figures"]), "numbered legacy cards are not a catalogue"
    assert categories & {"Comparative exon structure", "FGFR2 cassette evidence"}


def test_post_cluster_cards_are_pending_before_the_roundtrip(pre_cluster_run: Path):
    index = gc.flatten_for_gallery(_catalogue(pre_cluster_run, cluster_ready=False))
    pending = [f for f in index["figures"]
               if f["availability"]["status"] == gc.PENDING_CLUSTER]
    assert pending, "the domain and boundary cards must be visible and honest, not absent"
    for card in pending:
        assert not card.get("formats"), "a pending card must not offer a rendered file"
        assert "cluster" in card["availability"]["reason"].lower()


def test_no_domain_figure_claims_availability_before_the_cluster(pre_cluster_run: Path):
    index = gc.flatten_for_gallery(_catalogue(pre_cluster_run, cluster_ready=False))
    for card in index["figures"]:
        if card.get("category") in ("Domain architecture", "Comparative domain architecture"):
            assert card["availability"]["status"] != "available", card["figure_id"]


def test_both_species_get_their_own_scope(pre_cluster_run: Path):
    catalogue = _catalogue(pre_cluster_run, cluster_ready=False)
    scopes = catalogue["species_scopes"]
    assert set(scopes) == {"homo_sapiens", "felis_catus"}
    for sid, scope in scopes.items():
        assert scope["cards"], sid
        for card in scope["cards"]:
            assert card["species_id"] == sid, "a scope may not show another species' output"


def test_a_one_species_run_reports_comparison_as_not_applicable(tmp_path: Path):
    """One species is a scientific limit of the run, not a missing file."""
    run = _make_run(tmp_path, with_cluster=True)
    index = cm.build_index(sources=cm.run_sources(run))
    for model in index["models"]:
        model["species_id"] = "homo_sapiens"
    catalogue = rg.catalogue_for(run, index, True)
    comparative = [c for c in catalogue["comparative_cards"]
                   if str(c.get("category", "")).startswith("Comparative")
                   and not c.get("formats")]
    assert comparative
    for card in comparative:
        assert card["availability"]["status"] == "not_applicable"
        assert "one species" in card["availability"]["reason"]


# --------------------------------------------------------------------------- #
# 4. production wiring: no manual migration for a new run
# --------------------------------------------------------------------------- #
def test_the_production_index_builder_writes_the_modern_gallery(post_cluster_run: Path):
    import build_website_indices as bwi

    closure = post_cluster_run / "results" / "13_final_pre_interpro_closure"
    bwi.write_all(closure, post_cluster_run / "website_indices")

    index = json.loads((post_cluster_run / "website_indices" / "figure_index.json")
                       .read_text())
    assert index.get("default_scope") == "comparative"
    assert {f.get("species_id") for f in index["figures"]} >= {"homo_sapiens",
                                                               "felis_catus"}
    assert (post_cluster_run / "website_indices" / "figure_catalogue.json").is_file()


def test_the_modern_gallery_is_not_written_for_a_non_fgfr2_run(tmp_path: Path):
    run = tmp_path / "runs" / "2026-01-01_0000_hba_panthera_leo"
    (run / "results" / "13_final_pre_interpro_closure").mkdir(parents=True)
    assert rg.is_fgfr2_closure_run(run) is False


# --------------------------------------------------------------------------- #
# 5. regression: the validated dataset and the real run
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not FREEZE_INDICES.is_dir(), reason="validated dataset not present")
def test_the_validated_thirty_species_catalogue_is_unchanged():
    """Parameterising the builder may not move one card of the validated dataset."""
    model_index = json.loads(
        (FREEZE_INDICES / "protein_coordinate_model.json").read_text())
    derived = FREEZE_INDICES.parent
    built = gc.build_catalogue(model_index,
                               main_dir=derived / "figures" / "main",
                               comparative_dir=derived / "figures" / "comparative")
    stored = json.loads((FREEZE_INDICES / "figure_catalogue.json").read_text())

    def comparable(catalogue: dict) -> dict:
        copy = json.loads(json.dumps(catalogue))
        copy.pop("provenance", None)      # timestamp
        copy.pop("cluster_ready", None)   # added for runs that have a round-trip to wait for
        return copy

    assert comparable(built) == comparable(stored)


@pytest.mark.skipif(not FREEZE.is_dir(), reason="freeze not present")
def test_the_freeze_is_not_written_to():
    """The builders read the freeze; nothing here may write into it."""
    for module in (gc, cm, rg):
        source = Path(module.__file__).read_text()
        for verb in ("write_text(", "mkdir(", "open(\"w\""):
            for line in source.splitlines():
                if verb in line and "final_30_until_interpro_prepare" in line:
                    pytest.fail(f"{module.__name__} writes into the freeze: {line.strip()}")


@pytest.mark.skipif(not REAL_RUN.is_dir(), reason="the affected run is not present")
def test_the_affected_real_run_is_repaired():
    """The integration case: the run the regression was found on."""
    index = json.loads((REAL_RUN / "website_indices" / "figure_index.json").read_text())
    assert index["default_scope"] == "comparative"
    assert {f.get("species_id") for f in index["figures"] if f.get("species_id")} == {
        "homo_sapiens", "felis_catus"}
    assert "Framework" not in {f.get("category") for f in index["figures"]}

    status = json.loads((REAL_RUN / "status.json").read_text())
    assert "failed_reason" not in status, "a finished run may not keep an old failure note"

    verdict = ra.readiness(REAL_RUN, n_species=2, has_event=True)
    assert verdict.ready is True, verdict.reason
    # Its annotation is real, so every card either has a figure or names why it has none.
    for card in index["figures"]:
        state = card["availability"]["status"]
        assert state != gc.PENDING_CLUSTER, card["figure_id"]
        if state != "available":
            assert card["availability"]["reason"], card["figure_id"]
