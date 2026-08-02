"""Derived state must follow its inputs, and an empty view must name its own cause.

The *Equus quagga* run showed three empty pages over data that was complete on disk: its
website indices predated the closure directory they summarize, and nothing compared the
two. These tests pin the two properties that would have caught it — a derived artefact
older than its inputs is stale, and a view whose source table is missing says so instead
of implying the species lacks the biology.
"""
from __future__ import annotations

import ast
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "webapp" / "backend"))

from shared_gene_analysis import run_availability as ra  # noqa: E402
from exondomaincompare.contracts import write_freshness_contract  # noqa: E402

REAL_RUN = ROOT / "runs" / "2026-07-29_1217_fgfr2_equus_quagga"
FREEZE = ROOT / "results" / "final_30_until_interpro_prepare"


# --------------------------------------------------------------------------- #
# fixtures: a minimal run skeleton with the artefacts the views need
# --------------------------------------------------------------------------- #
CLOSURE_INPUTS = {
    "final_pre_interpro_truth_table.tsv":
        "species\tisoform\tprotein_id\tfinal_isoform_label\n"
        "equus_quagga\tIIIb\tXP_1\tIIIb\n",
    "tables/figure3C_exon_to_protein_cassette_coordinate_map.tsv":
        "species\tisoform\tfinal_isoform_label\tprotein_length\tcassette_start_aa\t"
        "cassette_end_aa\tfeature_type\texon_or_cds_id\tfeature_label\tblock_start_aa\t"
        "block_end_aa\n"
        "equus_quagga\tIIIb\tIIIb\t887\t383\t428\tCDS\tcds-1\texon 1\t1\t100\n",
    "tables/figure6B_species_resolved_IIIb_IIIc_cassette_residue_map.tsv":
        "species\tisoform\thuman_reference_residue_index\n",
    "tables/figure3B_IIIb_IIIc_cassette_amino_acid_motif_map.tsv":
        "human_reference_residue_index\thuman_IIIb_aa_one_letter\n1\tS\n",
    "MSA/final_fgfr2_full_length_protein_msa.aln.faa":
        ">equus_quagga|IIIb|XP_1|main_figure\nMDHTSV\n",
}

INDEX_PAYLOADS = {
    "run_index.json": {"run_id": "r", "gate_status": "pass", "kpi": {"species": 1}},
    "species_index.json": [{"species": "equus_quagga"}],
    "coordinate_track_index.json": {"available": True, "species": [{"species": "eq"}]},
    "cassette_residue_index.json": {"available": True, "species": [{"species": "eq"}],
                                    "evidence_level": "sequence_marker"},
    "msa_index.json": {"available": True,
                       "alignments": {"full_length": {"available": True, "rows": [{}]}}},
    "synteny_locus_index.json": {"available": True, "species": [{"species": "eq"}]},
    "figure_index.json": {"figures": [{"key": "f1"}]},
    "download_index.json": [{"group": "Tables", "label": "Final truth table"}],
    "species_domain_architecture.json": {"available": True, "species": [{}]},
    "domain_architecture_index.json": {"available": True, "species": [{}]},
    "boundary_consistency_summary.json": {"available": True, "total_primary_proteins": 2},
}


def _make_run(tmp_path: Path, *, with_cluster: bool = False,
              indices: bool = True) -> Path:
    run = tmp_path / "runs" / "2026-01-01_0000_fgfr2_test_species"
    closure = run / "results" / "13_final_pre_interpro_closure"
    for rel, text in CLOSURE_INPUTS.items():
        target = closure / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")
    freeze = closure / "freeze"
    freeze.mkdir(parents=True, exist_ok=True)
    (freeze / "final_pre_interpro_proteins_primary.faa").write_text(
        ">equus_quagga|IIIb|XP_1\nMDHTSV\n", encoding="utf-8")
    (run / "status.json").write_text(json.dumps({"species_count": 1}), encoding="utf-8")

    if with_cluster:
        ips = run / "results" / "14_interproscan" / "primary" / "output"
        ips.mkdir(parents=True, exist_ok=True)
        (ips / "input.fasta.tsv").write_text("XP_1\tPF00069\n", encoding="utf-8")
        tm = (run / "results" / "15_exon_domain_boundary_post_interpro"
              / "pytmhmm_primary" / "output")
        tm.mkdir(parents=True, exist_ok=True)
        (tm / "pytmhmm_transmembrane_hits.tsv").write_text("XP_1\t1\t20\n",
                                                           encoding="utf-8")

    if indices:
        time.sleep(0.02)
        wi = run / "website_indices"
        wi.mkdir(parents=True, exist_ok=True)
        for name, payload in INDEX_PAYLOADS.items():
            (wi / name).write_text(json.dumps(payload), encoding="utf-8")
    return run


@pytest.fixture()
def run_dir(tmp_path: Path) -> Path:
    return _make_run(tmp_path)


@pytest.fixture()
def complete_run(tmp_path: Path) -> Path:
    """A run that has been all the way through, cluster round-trip included.

    ``run_dir`` stops before the round-trip, so it is not a finished run: its domain and
    boundary views are still waiting for annotation. Treating it as complete is what let a
    pre-cluster run be published as "Results ready".
    """
    return _make_run(tmp_path, with_cluster=True)


# --------------------------------------------------------------------------- #
# staleness: the defect that emptied three views
# --------------------------------------------------------------------------- #
def test_indices_built_after_their_inputs_are_current(run_dir: Path):
    stale, reason = ra.indices_are_stale(run_dir)
    assert stale is False, reason


def test_indices_older_than_the_closure_are_stale(run_dir: Path):
    """Exactly the *Equus quagga* situation: the closure is rewritten, indices are not."""
    later = time.time() + 60
    target = (run_dir / "results" / "13_final_pre_interpro_closure"
              / "final_pre_interpro_truth_table.tsv")
    os.utime(target, (later, later))

    stale, reason = ra.indices_are_stale(run_dir)
    assert stale is True
    assert "13_final_pre_interpro_closure" in reason


def test_a_stale_index_is_not_reported_as_an_available_view(run_dir: Path):
    later = time.time() + 60
    os.utime(run_dir / "results" / "13_final_pre_interpro_closure"
             / "MSA" / "final_fgfr2_full_length_protein_msa.aln.faa", (later, later))

    states = {s.view: s for s in ra.view_states(run_dir, n_species=1)}
    assert states["msa"].state == ra.STALE
    assert states["msa"].available is False


def test_an_upstream_rebuild_that_changed_nothing_downstream_is_not_stale(run_dir: Path):
    """A registry rebuild alone must not condemn indices that still match the closure."""
    later = time.time() + 60
    registry = run_dir / "results" / "01_species_registry"
    registry.mkdir(parents=True, exist_ok=True)
    (registry / "species_registry.tsv").write_text("species\n", encoding="utf-8")
    os.utime(registry / "species_registry.tsv", (later, later))

    stale, _ = ra.indices_are_stale(run_dir)
    assert stale is False


def test_final_run_manifest_does_not_invalidate_current_indices(run_dir: Path):
    write_freshness_contract(
        run_dir, run_dir / "website_indices", generator="test")
    setup = run_dir / "results" / "00_run_setup"
    setup.mkdir(parents=True, exist_ok=True)
    (setup / "post_interpro_run_manifest.json").write_text(
        json.dumps({"status": "complete"}), encoding="utf-8")

    stale, reason = ra.indices_are_stale(run_dir)
    assert stale is False, reason


# --------------------------------------------------------------------------- #
# technically missing vs. scientifically unavailable
# --------------------------------------------------------------------------- #
def test_a_missing_source_table_is_technically_missing(run_dir: Path):
    (run_dir / "results" / "13_final_pre_interpro_closure" / "tables"
     / "figure3C_exon_to_protein_cassette_coordinate_map.tsv").unlink()

    states = {s.view: s for s in ra.view_states(run_dir, n_species=1)}
    assert states["exon_map"].state == ra.TECHNICALLY_MISSING
    assert "Retry local analysis" in states["exon_map"].reason
    assert states["exon_map"].missing_inputs


def test_a_header_only_source_table_is_scientifically_unavailable(run_dir: Path):
    """Written and empty is a finding about the data, not about the software."""
    (run_dir / "results" / "13_final_pre_interpro_closure" / "tables"
     / "figure3C_exon_to_protein_cassette_coordinate_map.tsv").write_text(
        "species\tisoform\n", encoding="utf-8")

    states = {s.view: s for s in ra.view_states(run_dir, n_species=1)}
    assert states["exon_map"].state == ra.SCIENTIFICALLY_UNAVAILABLE
    assert "Retry" not in states["exon_map"].reason


def test_the_two_states_are_never_conflated():
    assert ra.TECHNICALLY_MISSING != ra.SCIENTIFICALLY_UNAVAILABLE
    assert ra.TECHNICALLY_MISSING in ra.BLOCKING_STATES
    # An honest scientific absence must not block a run from being finished.
    assert ra.SCIENTIFICALLY_UNAVAILABLE not in ra.BLOCKING_STATES
    assert ra.NOT_APPLICABLE not in ra.BLOCKING_STATES


def test_a_missing_index_over_present_sources_is_technically_missing(run_dir: Path):
    (run_dir / "website_indices" / "msa_index.json").unlink()

    states = {s.view: s for s in ra.view_states(run_dir, n_species=1)}
    assert states["msa"].state == ra.TECHNICALLY_MISSING


def test_a_gene_without_an_event_region_has_no_cassette_view(run_dir: Path):
    states = {s.view: s for s in ra.view_states(run_dir, n_species=1, has_event=False)}
    assert states["event_region"].state == ra.NOT_APPLICABLE


def test_domain_views_are_pending_before_the_cluster_returns(run_dir: Path):
    states = {s.view: s for s in ra.view_states(run_dir, n_species=1)}
    assert states["domain_architecture"].state == ra.PENDING
    assert "cluster" in states["domain_architecture"].reason.lower()


# --------------------------------------------------------------------------- #
# readiness
# --------------------------------------------------------------------------- #
def test_a_complete_run_is_ready(complete_run: Path):
    verdict = ra.readiness(complete_run, n_species=1)
    assert verdict.ready is True, verdict.reason


def test_a_run_awaiting_its_cluster_roundtrip_is_not_ready(run_dir: Path):
    """The pre-cluster state is not a finished state, however complete its local stages."""
    verdict = ra.readiness(run_dir, n_species=1)
    assert verdict.ready is False
    assert "domain_architecture" in [v.view for v in verdict.blocking]


def test_results_ready_is_blocked_while_a_required_output_is_missing(run_dir: Path):
    (run_dir / "results" / "13_final_pre_interpro_closure" / "MSA"
     / "final_fgfr2_full_length_protein_msa.aln.faa").unlink()

    verdict = ra.readiness(run_dir, n_species=1)
    assert verdict.ready is False
    assert "MSA" in verdict.reason
    assert [v.view for v in verdict.blocking] == ["msa"]


def test_results_ready_is_blocked_while_the_indices_are_stale(run_dir: Path):
    later = time.time() + 60
    os.utime(run_dir / "results" / "13_final_pre_interpro_closure"
             / "final_pre_interpro_truth_table.tsv", (later, later))

    assert ra.readiness(run_dir, n_species=1).ready is False


def test_domain_views_are_required_once_the_cluster_output_is_in(tmp_path: Path):
    run = _make_run(tmp_path, with_cluster=True)
    (run / "website_indices" / "species_domain_architecture.json").write_text(
        json.dumps({"available": False, "species": []}), encoding="utf-8")

    verdict = ra.readiness(run, n_species=1)
    assert verdict.ready is False
    assert "domain_architecture" in [v.view for v in verdict.blocking]


def test_an_honest_scientific_absence_does_not_block_readiness(complete_run: Path):
    """A run whose synteny found nothing is still a finished run."""
    (complete_run / "website_indices" / "synteny_locus_index.json").write_text(
        json.dumps({"available": False, "species": []}), encoding="utf-8")

    verdict = ra.readiness(complete_run, n_species=1)
    assert verdict.ready is True, verdict.reason


# --------------------------------------------------------------------------- #
# dependency invalidation in the pipeline
# --------------------------------------------------------------------------- #
def test_the_pre_interpro_runner_invalidates_the_derived_layer():
    import run_pre_interpro_for_run as runner

    assert hasattr(runner, "invalidate_derived_layer")
    assert hasattr(runner, "refresh_derived_layer")
    source = (ROOT / "scripts" / "run_pre_interpro_for_run.py").read_text(encoding="utf-8")
    # A successful pre-InterPro run must not leave the run's own success flag from a
    # previous post-cluster pass in place.
    assert 'post_interpro_status="stale"' in source
    assert "refresh_derived_layer(run_dir, args.run_id)" in source


def test_invalidation_marks_the_derived_layer_stale(tmp_path: Path):
    import run_pre_interpro_for_run as runner

    run = _make_run(tmp_path)
    (run / "status.json").write_text(
        json.dumps({"post_interpro_status": "complete",
                    "website_indices_status": "complete"}), encoding="utf-8")

    runner.invalidate_derived_layer(run)

    st = json.loads((run / "status.json").read_text(encoding="utf-8"))
    assert st["post_interpro_status"] == "stale"
    assert st["website_indices_status"] == "stale"
    assert st["derived_layer_invalidated_at"]


def test_a_rebuild_reuses_returned_cluster_output_rather_than_resubmitting(tmp_path: Path):
    import run_pre_interpro_for_run as runner

    without = _make_run(tmp_path / "a")
    with_cluster = _make_run(tmp_path / "b", with_cluster=True)
    assert runner._cluster_outputs_valid(without) is False
    assert runner._cluster_outputs_valid(with_cluster) is True


def test_empty_cluster_output_files_do_not_count_as_valid(tmp_path: Path):
    import run_pre_interpro_for_run as runner

    run = _make_run(tmp_path, with_cluster=True)
    (run / "results" / "14_interproscan" / "primary" / "output"
     / "input.fasta.tsv").write_text("", encoding="utf-8")
    assert runner._cluster_outputs_valid(run) is False


# --------------------------------------------------------------------------- #
# roundtrip finalize ordering
# --------------------------------------------------------------------------- #
def test_finalize_writes_the_state_before_setting_the_phase():
    """The phase must not be set from a status the same method is about to write."""
    source = (ROOT / "scripts" / "interpro_cluster" / "run_cluster_roundtrip.py"
              ).read_text(encoding="utf-8")
    tree = ast.parse(source)
    finalize = next(n for n in ast.walk(tree)
                    if isinstance(n, ast.FunctionDef) and n.name == "finalize")
    calls = []
    for node in ast.walk(finalize):
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Attribute):
            calls.append(node.func.attr)
        elif isinstance(node.func, ast.Name):
            calls.append(node.func.id)
    assert "write_json" in calls and "set_phase" in calls
    # write_json before set_phase, so a recorded phase can never contradict the file.
    assert calls.index("write_json") < calls.index("set_phase")


def test_post_interpro_no_longer_sets_complete_before_finalize():
    source = (ROOT / "scripts" / "interpro_cluster" / "run_cluster_roundtrip.py"
              ).read_text(encoding="utf-8")
    tree = ast.parse(source)
    post = next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "post_interpro")
    for node in ast.walk(post):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "set_phase" and node.args
                and isinstance(node.args[0], ast.Constant)):
            assert node.args[0].value != "complete", (
                "post_interpro must leave the end state to finalize()")


def test_finalize_derives_the_end_state_from_the_artefacts(tmp_path: Path, monkeypatch):
    sys.path.insert(0, str(ROOT / "scripts" / "interpro_cluster"))
    import run_cluster_roundtrip as rt_mod

    run = _make_run(tmp_path, with_cluster=True)
    (run / "logs").mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(rt_mod, "REPO", ROOT)
    rt = rt_mod.Roundtrip.__new__(rt_mod.Roundtrip)
    rt.run_id = run.name
    rt.run_dir = run
    rt.status_path = run / "status.json"
    rt.log_path = run / "logs" / "cluster_roundtrip.log"
    rt.log = lambda *a, **k: None
    recorded = {}
    rt.set_phase = lambda phase, **kw: recorded.setdefault("phase", phase)

    # The post-cluster runner's own vocabulary, which the old check rejected outright.
    rt.status_path.write_text(json.dumps({"status": "complete", "species_count": 1}),
                              encoding="utf-8")
    rt.finalize()

    st = json.loads(rt.status_path.read_text(encoding="utf-8"))
    assert st["status"] == "results_ready"
    assert recorded["phase"] == "complete"


def test_shared_finalizer_clears_a_superseded_roundtrip_failure(
        complete_run: Path, monkeypatch: pytest.MonkeyPatch):
    from shared_gene_analysis import finalize_run_status as finalizer

    path = complete_run / "status.json"
    path.write_text(json.dumps({
        "species_count": 1,
        "status": "failed",
        "blocking_analyses": ["overview=stale"],
        "cluster_roundtrip": {
            "phase": "post_interpro_failed",
            "reason": "indices became stale",
        },
    }), encoding="utf-8")
    monkeypatch.setattr(finalizer, "evaluate_run", lambda _run: {
        "run_id": complete_run.name,
        "status": finalizer.RESULTS_READY,
        "reason": "Every required view has current data.",
        "blocking": [],
        "cluster_outputs": "current",
        "cluster_report": {},
        "decided": True,
    })

    finalizer.finalize(complete_run)
    status = json.loads(path.read_text(encoding="utf-8"))
    assert status["status"] == "results_ready"
    assert status["cluster_roundtrip"]["phase"] == "complete"
    assert status["cluster_roundtrip"]["reason"] == ""
    assert "blocking_analyses" not in status


def test_finalize_refuses_a_run_whose_required_views_are_missing(tmp_path: Path,
                                                                monkeypatch):
    sys.path.insert(0, str(ROOT / "scripts" / "interpro_cluster"))
    import run_cluster_roundtrip as rt_mod

    run = _make_run(tmp_path, with_cluster=True)
    (run / "logs").mkdir(parents=True, exist_ok=True)
    (run / "results" / "13_final_pre_interpro_closure" / "MSA"
     / "final_fgfr2_full_length_protein_msa.aln.faa").unlink()
    monkeypatch.setattr(rt_mod, "REPO", ROOT)
    rt = rt_mod.Roundtrip.__new__(rt_mod.Roundtrip)
    rt.run_id = run.name
    rt.run_dir = run
    rt.status_path = run / "status.json"
    rt.log = lambda *a, **k: None
    recorded = {}
    rt.set_phase = lambda phase, **kw: recorded.setdefault("phase", phase)
    rt.status_path.write_text(json.dumps({"status": "complete", "species_count": 1}),
                              encoding="utf-8")

    with pytest.raises(SystemExit):
        rt.finalize()
    assert recorded["phase"] == "post_interpro_failed"


# --------------------------------------------------------------------------- #
# indices carry their own reason; the frontend consumes it
# --------------------------------------------------------------------------- #
def test_the_three_repaired_indices_carry_an_availability_block():
    import build_website_indices as bwi

    closure = REAL_RUN / "results" / "13_final_pre_interpro_closure"
    if not closure.is_dir():
        pytest.skip("real Equus quagga run not present")
    for builder in (bwi.build_coordinate_track_index, bwi.build_msa_index,
                    bwi.build_cassette_residue_index):
        payload = builder(closure)
        block = payload.get("availability")
        assert block, f"{builder.__name__} carries no availability block"
        assert block["state"] in ra.STATES
        assert block["reason"]


def test_a_missing_source_table_makes_the_index_say_technically_missing(tmp_path: Path):
    import build_website_indices as bwi

    run = _make_run(tmp_path)
    closure = run / "results" / "13_final_pre_interpro_closure"
    (closure / "tables" / "figure3C_exon_to_protein_cassette_coordinate_map.tsv").unlink()

    block = bwi.build_coordinate_track_index(closure)["availability"]
    assert block["state"] == ra.TECHNICALLY_MISSING
    assert block["missing_inputs"]


def test_the_frontend_reads_the_canonical_reason_instead_of_hardcoding_one():
    common = (ROOT / "webapp" / "frontend" / "src" / "pages" / "viewers"
              / "common.js").read_text(encoding="utf-8")
    assert "export function unavailableState" in common
    for state in (ra.TECHNICALLY_MISSING, ra.SCIENTIFICALLY_UNAVAILABLE, ra.PENDING,
                  ra.STALE, ra.NOT_APPLICABLE):
        assert state in common, f"{state} is not handled in the frontend"

    viewers = ROOT / "webapp" / "frontend" / "src" / "pages" / "viewers"
    for name in ("CoordinateTrack.jsx", "CassetteExplorer.jsx", "MsaExplorer.jsx"):
        text = (viewers / name).read_text(encoding="utf-8")
        assert "unavailableState" in text, f"{name} still hardcodes its empty state"
        # The misleading message that told users to wait for a step that had finished.
        assert "becomes available once the pre-InterPro" not in text
        assert "No alignment files were found in this run." not in text
        assert "Figure 3C source table was not found." not in text


# --------------------------------------------------------------------------- #
# Data & Downloads inventory
# --------------------------------------------------------------------------- #
def test_downloads_expose_more_than_the_final_truth_table():
    import build_website_indices as bwi

    closure = REAL_RUN / "results" / "13_final_pre_interpro_closure"
    if not closure.is_dir():
        pytest.skip("real Equus quagga run not present")
    items = bwi.build_download_index(closure)
    groups = {i["group"] for i in items}
    for expected in ("Freeze", "Tables", "Alignment", "Gene models", "Selected models",
                     "Domain annotation", "Boundary analysis", "Figures"):
        assert expected in groups, f"Data & Downloads has no {expected} group"
    assert len(items) > 20, f"only {len(items)} downloadable files"


def test_downloads_reach_the_model_tables_and_cluster_results():
    import build_website_indices as bwi

    closure = REAL_RUN / "results" / "13_final_pre_interpro_closure"
    if not closure.is_dir():
        pytest.skip("real Equus quagga run not present")
    names = {Path(i["path"]).name for i in bwi.build_download_index(closure)}
    for expected in ("genes.tsv", "transcripts.tsv", "exons.tsv", "cds_features.tsv",
                     "final_fgfr2_full_length_protein_msa.aln.faa",
                     "figure3C_exon_to_protein_cassette_coordinate_map.tsv",
                     "final_pre_interpro_proteins_primary.faa",
                     "input.fasta.tsv", "pytmhmm_transmembrane_hits.tsv"):
        assert expected in names, f"{expected} is not downloadable"


def test_the_upstream_role_fasta_is_labelled_as_provenance():
    """Its header role is inverted relative to the final label; the label must warn."""
    import build_website_indices as bwi

    closure = REAL_RUN / "results" / "13_final_pre_interpro_closure"
    if not closure.is_dir():
        pytest.skip("real Equus quagga run not present")
    entry = next((i for i in bwi.build_download_index(closure)
                  if Path(i["path"]).name == "selected_fgfr2_proteins.faa"), None)
    assert entry is not None
    assert "provenance" in entry["label"].lower()
    assert "upstream" in entry["label"].lower()


# --------------------------------------------------------------------------- #
# final label vs. upstream role
# --------------------------------------------------------------------------- #
MARKER_TRUTH = {"XP_046510917.1": "IIIb", "XP_046510919.1": "IIIc"}


def _labelled_refs(node, out):
    if isinstance(node, dict):
        pid = next((node.get(k) for k in ("protein_id", "protein")
                    if isinstance(node.get(k), str)), None)
        if pid in MARKER_TRUTH:
            for key in ("isoform", "final_isoform_label", "panel", "label"):
                value = node.get(key)
                if isinstance(value, str) and value in ("IIIb", "IIIc"):
                    out.append((pid, key, value))
        for value in node.values():
            _labelled_refs(value, out)
    elif isinstance(node, list):
        for value in node:
            _labelled_refs(value, out)
    return out


@pytest.mark.skipif(not REAL_RUN.is_dir(), reason="real Equus quagga run not present")
def test_no_index_takes_the_isoform_label_from_the_inverted_upstream_role():
    checked = 0
    for path in sorted((REAL_RUN / "website_indices").glob("*.json")):
        for pid, key, value in _labelled_refs(json.loads(path.read_text(encoding="utf-8")),
                                              []):
            assert value == MARKER_TRUTH[pid], (
                f"{path.name} labels {pid} as {value} via '{key}'; the marker evidence "
                f"says {MARKER_TRUTH[pid]}")
            checked += 1
    assert checked > 0, "no isoform-labelled protein reference found to check"


@pytest.mark.skipif(not REAL_RUN.is_dir(), reason="real Equus quagga run not present")
def test_every_alignment_header_carries_the_marker_based_label():
    msa = REAL_RUN / "results" / "13_final_pre_interpro_closure" / "MSA"
    seen = 0
    for path in sorted(msa.glob("*.aln.faa")):
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.startswith(">"):
                continue
            parts = line[1:].split("|")
            pid = next((p for p in parts if p in MARKER_TRUTH), None)
            if pid:
                assert parts[1] == MARKER_TRUTH[pid], f"{path.name}: {line}"
                seen += 1
    assert seen >= 4


@pytest.mark.skipif(not REAL_RUN.is_dir(), reason="real Equus quagga run not present")
def test_the_markers_themselves_still_decide_the_label():
    primary = (REAL_RUN / "results" / "13_final_pre_interpro_closure" / "freeze"
               / "final_pre_interpro_proteins_primary.faa")
    seqs, header = {}, None
    for line in primary.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            header = next((p for p in line[1:].split("|") if p in MARKER_TRUTH), None)
            if header:
                seqs[header] = ""
        elif header:
            seqs[header] += line.strip()
    assert "SGINSSN" in seqs["XP_046510917.1"]
    assert "GVNTTDKEI" not in seqs["XP_046510917.1"]
    assert "GVNTTDKEI" in seqs["XP_046510919.1"]
    assert "SGINSSN" not in seqs["XP_046510919.1"]


# --------------------------------------------------------------------------- #
# the real run, end to end
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not REAL_RUN.is_dir(), reason="real Equus quagga run not present")
def test_the_repaired_run_has_no_stale_indices_and_is_ready():
    stale, reason = ra.indices_are_stale(REAL_RUN)
    assert stale is False, reason
    verdict = ra.readiness(REAL_RUN, n_species=1, has_event=True)
    assert verdict.ready is True, verdict.reason


@pytest.mark.skipif(not REAL_RUN.is_dir(), reason="real Equus quagga run not present")
def test_the_three_previously_empty_views_carry_real_content():
    wi = REAL_RUN / "website_indices"

    coord = json.loads((wi / "coordinate_track_index.json").read_text(encoding="utf-8"))
    assert coord["available"] is True
    panels = coord["species"][0]["panels"]
    assert set(panels) == {"IIIb", "IIIc"}
    for iso, panel in panels.items():
        assert panel["final_isoform_label"] == iso
        assert panel["blocks"], f"{iso} panel has no exon blocks"
        assert panel["cassette_start_aa"] and panel["cassette_end_aa"]

    cassette = json.loads((wi / "cassette_residue_index.json").read_text(encoding="utf-8"))
    assert cassette["available"] is True
    assert cassette["evidence_level"] != "none"
    assert cassette["discriminating"]

    msa = json.loads((wi / "msa_index.json").read_text(encoding="utf-8"))
    assert msa["available"] is True
    full = msa["alignments"]["full_length"]
    assert full["available"] is True
    assert {r["protein_id"] for r in full["rows"]} == set(MARKER_TRUTH)


@pytest.mark.skipif(not REAL_RUN.is_dir(), reason="real Equus quagga run not present")
def test_the_persisted_roundtrip_status_agrees_with_the_derived_status():
    st = json.loads((REAL_RUN / "status.json").read_text(encoding="utf-8"))
    assert st["status"] == "results_ready"
    assert (st.get("cluster_roundtrip") or {}).get("phase") == "complete"

    import main as backend
    model = backend.derive_status_model(REAL_RUN.resolve())
    assert model["status"] == "results_ready"
    assert model["readiness"]["ready"] is True


@pytest.mark.skipif(not REAL_RUN.is_dir(), reason="real Equus quagga run not present")
def test_the_api_reports_every_view_as_available_for_the_repaired_run():
    import main as backend

    model = backend.derive_status_model(REAL_RUN.resolve())
    for view, node in (model["view_availability"] or {}).items():
        assert node["state"] == ra.AVAILABLE, f"{view}: {node['state']} — {node['reason']}"


@pytest.mark.skipif(not REAL_RUN.is_dir(), reason="real Equus quagga run not present")
def test_the_recovered_model_counts_are_unchanged_by_this_repair():
    models = REAL_RUN / "results" / "02_models"
    counts = {}
    for name in ("genes.tsv", "transcripts.tsv", "exons.tsv", "cds_features.tsv"):
        rows = (models / name).read_text(encoding="utf-8").splitlines()
        counts[name] = max(len(rows) - 1, 0)
    assert counts == {"genes.tsv": 1, "transcripts.tsv": 12, "exons.tsv": 206,
                      "cds_features.tsv": 205}


# --------------------------------------------------------------------------- #
# regressions elsewhere
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not FREEZE.is_dir(), reason="validated freeze not present")
def test_the_validated_freeze_is_untouched():
    proc = subprocess.run(["git", "status", "--porcelain", "--", str(FREEZE)],
                          cwd=str(ROOT), capture_output=True, text=True, check=False)
    assert proc.stdout.strip() == "", f"freeze modified:\n{proc.stdout}"


@pytest.mark.skipif(not FREEZE.is_dir(), reason="validated freeze not present")
def test_the_validated_dataset_still_reports_every_view_as_available():
    """The extended download index and availability blocks must not regress the freeze."""
    import build_website_indices as bwi

    closure = next((p for p in FREEZE.rglob("*final_pre_interpro_closure")
                    if p.is_dir()), None)
    if closure is None:
        pytest.skip("freeze holds no closure directory")
    for builder in (bwi.build_coordinate_track_index, bwi.build_msa_index,
                    bwi.build_cassette_residue_index):
        payload = builder(closure)
        assert payload["available"] is True, builder.__name__
        assert payload["availability"]["state"] == ra.AVAILABLE

    items = bwi.build_download_index(closure)
    assert len(items) >= 13
