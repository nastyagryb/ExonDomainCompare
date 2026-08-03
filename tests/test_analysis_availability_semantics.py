#!/usr/bin/env python3
"""Analysis applicability, and how it differs from run completion.

Chicken MC1R is the case that exposed the confusion. Its protein is encoded by a single
coding exon, so it has zero internal coding-exon boundaries and no exon–domain boundary
analysis can exist for it. The pipeline handled that correctly — it produced a boundary
table with no rows because no boundary exists — and the completion check, which demanded at
least one row, concluded the run was ``post_cluster_partial``. The gene's exon structure was
reported as a pipeline defect.

The tests below fix the distinction in place:

* a *run* is complete when every **applicable** stage finished;
* an *analysis* is applicable when its biological prerequisites exist.

``not_applicable`` is a successfully resolved answer and must never block completion.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from shutil import which

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "webapp" / "backend"))

from exondomaincompare.shared_gene_analysis import analysis_availability as aa  # noqa: E402

RUNS = ROOT / "runs"
MC1R = RUNS / "2026-07-29_1306_mc1r_gallus_gallus"
HBA = RUNS / "2026-07-29_1347_hba_panthera_leo"
NKD2 = RUNS / "2026-07-29_1502_nkd2_panthera_onca"
AKT1 = RUNS / "2026-07-29_1526_akt1_mus_musculus"

def _needs(path: Path):
    return pytest.mark.skipif(not path.is_dir(), reason=f"{path.name} is not present")


# --------------------------------------------------------------------------- #
# The arithmetic of internal boundaries
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("coding_exons,expected", [
    (0, 0), (1, 0), (2, 1), (3, 2), (10, 9), (17, 16),
])
def test_internal_boundaries_are_one_fewer_than_the_coding_exons(coding_exons, expected):
    assert aa.internal_boundary_count(coding_exons) == expected


def test_one_coding_exon_produces_zero_internal_boundaries():
    # The single fact the whole MC1R correction rests on.
    assert aa.internal_boundary_count(1) == 0


def test_a_negative_or_absent_count_never_produces_a_negative_boundary_count():
    assert aa.internal_boundary_count(-3) == 0


# --------------------------------------------------------------------------- #
# Canonical states
# --------------------------------------------------------------------------- #
def test_the_seven_canonical_states_stay_distinct():
    assert len(set(aa.STATES)) == 7
    assert "unavailable" not in aa.STATES, "the states must not collapse into one"


def test_a_resolved_not_applicable_analysis_does_not_block_a_run():
    assert aa.NOT_APPLICABLE not in aa.BLOCKING_STATES
    assert aa.SCIENTIFICALLY_UNAVAILABLE not in aa.BLOCKING_STATES


def test_a_genuine_gap_still_blocks_a_run():
    for state in (aa.TECHNICALLY_MISSING, aa.PENDING, aa.STALE, aa.FAILED):
        assert state in aa.BLOCKING_STATES


# --------------------------------------------------------------------------- #
# Synthetic runs: applicability from counts alone
# --------------------------------------------------------------------------- #
def _write_run(tmp_path: Path, *, exons: int, proteins: dict, transcripts=None,
              boundary_rows: int = 0, domain_rows: int = 1,
              cluster: bool = True) -> Path:
    """A minimal run directory with only the canonical tables the module reads."""
    run = tmp_path / "run"
    core = run / aa.CORE
    core.mkdir(parents=True)

    pid = next(iter(proteins))
    rows = ["\t".join(["species_id", "protein_id", "transcript_id", "exon_id",
                       "exon_number"])]
    for i in range(1, exons + 1):
        rows.append(f"sp\t{pid}\tT1\t{pid}:cds{i}\t{i}")
    (core / "exon_protein_map.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    iso = ["\t".join(["species_id", "protein_id", "transcript_id", "protein_length",
                      "primary_status"])]
    for n, (protein, seq) in enumerate(proteins.items()):
        iso.append(f"sp\t{protein}\tT{n + 1}\t{len(seq)}\t"
                   f"{'primary' if n == 0 else 'alternative'}")
    (core / "protein_isoform_index.tsv").write_text("\n".join(iso) + "\n",
                                                    encoding="utf-8")

    gm = ["\t".join(["species_id", "gene_id", "transcript_id", "protein_id"])]
    for transcript, protein in (transcripts or
                                [(f"T{n + 1}", p) for n, p in enumerate(proteins)]):
        gm.append(f"sp\tg1\t{transcript}\t{protein}")
    (core / "gene_model_index.tsv").write_text("\n".join(gm) + "\n", encoding="utf-8")

    (core / "proteins_all_isoforms.faa").write_text(
        "".join(f">{p}\n{s}\n" for p, s in proteins.items()), encoding="utf-8")

    bnd = ["\t".join(["species_id", "protein_id", "exon_boundary_id"])]
    bnd += [f"sp\t{pid}\tb{i}" for i in range(1, boundary_rows + 1)]
    (core / "exon_domain_boundary_distances.tsv").write_text("\n".join(bnd) + "\n",
                                                            encoding="utf-8")
    dom = ["\t".join(["species_id", "protein_id", "interpro_accession"])]
    dom += [f"sp\t{pid}\tIPR{i:06d}" for i in range(1, domain_rows + 1)]
    (core / "domain_features.tsv").write_text("\n".join(dom) + "\n", encoding="utf-8")

    if cluster:
        out = run / "results" / "14_interproscan" / "primary" / "output"
        out.mkdir(parents=True)
        (out / "primary.tsv").write_text("hit\n", encoding="utf-8")
    return run


def test_a_single_exon_protein_makes_the_boundary_analysis_not_applicable(tmp_path):
    run = _write_run(tmp_path, exons=1, proteins={"P1": "MKV"})
    state = aa.build_manifest(run).by_name()["boundary_analysis"]
    assert state.status == aa.NOT_APPLICABLE
    assert state.reason_code == aa.SINGLE_CODING_EXON
    assert state.prerequisite_name == "internal_coding_exon_boundary_count"
    assert state.prerequisite_count == 0


def test_the_single_exon_message_and_badge_are_the_agreed_wording(tmp_path):
    run = _write_run(tmp_path, exons=1, proteins={"P1": "MKV"})
    state = aa.build_manifest(run).by_name()["boundary_analysis"]
    assert state.user_message == (
        "The selected protein is encoded by one coding exon and therefore has no "
        "internal coding-exon boundaries to analyse.")
    assert state.badge == "No internal boundaries"


def test_a_single_exon_run_can_still_be_ready(tmp_path):
    manifest = aa.build_manifest(_write_run(tmp_path, exons=1, proteins={"P1": "MKV"}))
    assert manifest.ready
    assert manifest.blocking == []


def test_no_synthetic_boundary_is_invented_for_a_single_exon_protein(tmp_path):
    run = _write_run(tmp_path, exons=1, proteins={"P1": "MKV"}, boundary_rows=0)
    aa.write_manifest(run)
    table = run / aa.CORE / "exon_domain_boundary_distances.tsv"
    body = [ln for ln in table.read_text(encoding="utf-8").splitlines()[1:] if ln.strip()]
    assert body == [], "a boundary row was fabricated to make the run look complete"


def test_writing_the_manifest_creates_a_status_record_not_a_scientific_table(tmp_path):
    run = _write_run(tmp_path, exons=1, proteins={"P1": "MKV"})
    written = aa.write_manifest(run)
    assert written.name == "analysis_availability.json"
    payload = json.loads(written.read_text(encoding="utf-8"))
    assert payload["ready"] is True
    assert "boundary_analysis" in payload["not_applicable"]


def test_two_coding_exons_make_the_boundary_analysis_applicable(tmp_path):
    run = _write_run(tmp_path, exons=2, proteins={"P1": "MKV"}, boundary_rows=1)
    state = aa.build_manifest(run).by_name()["boundary_analysis"]
    assert state.status == aa.AVAILABLE
    assert state.prerequisite_count == 1


# --------------------------------------------------------------------------- #
# One protein sequence
# --------------------------------------------------------------------------- #
def test_one_unique_protein_sequence_makes_isoform_comparison_not_applicable(tmp_path):
    run = _write_run(tmp_path, exons=3, proteins={"P1": "MKVLA"}, boundary_rows=2)
    state = aa.build_manifest(run).by_name()["protein_isoform_comparison"]
    assert state.status == aa.NOT_APPLICABLE
    assert state.user_message == (
        "Only one distinct translated protein sequence is available. At least two "
        "distinct protein sequences are required to detect protein-level isoform "
        "differences.")


def test_one_unique_protein_sequence_makes_candidate_analysis_not_applicable(tmp_path):
    run = _write_run(tmp_path, exons=3, proteins={"P1": "MKVLA"}, boundary_rows=2)
    manifest = aa.build_manifest(run)
    state = manifest.by_name()["protein_difference_candidate_analysis"]
    assert state.status == aa.NOT_APPLICABLE
    assert manifest.prerequisites.protein_difference_candidate_count == 0


def test_several_transcripts_encoding_one_protein_are_one_isoform(tmp_path):
    # Two transcript models, one amino-acid sequence: a transcript-level difference, not a
    # protein-level one. Counting transcripts as isoforms would invent a difference.
    run = _write_run(tmp_path, exons=3, proteins={"P1": "MKVLA"},
                     transcripts=[("T1", "P1"), ("T2", "P1")], boundary_rows=2)
    pre = aa.prerequisites(run)
    assert pre.transcript_model_count == 2
    assert pre.unique_protein_sequence_count == 1
    assert pre.protein_difference_candidate_count == 0


def test_transcript_only_variation_stays_reportable(tmp_path):
    run = _write_run(tmp_path, exons=3, proteins={"P1": "MKVLA"},
                     transcripts=[("T1", "P1"), ("T2", "P1"), ("T3", "P1")],
                     boundary_rows=2)
    pre = aa.prerequisites(run)
    assert pre.transcript_only_variant_count == 2
    assert pre.transcripts_by_protein_sequence["P1"] == ["T1", "T2", "T3"]


def test_two_distinct_sequences_make_isoform_comparison_applicable(tmp_path):
    run = _write_run(tmp_path, exons=3, proteins={"P1": "MKVLA", "P2": "MKVLAQQ"},
                     boundary_rows=2)
    names = aa.build_manifest(run).by_name()
    assert names["protein_isoform_comparison"].status == aa.AVAILABLE
    assert names["protein_difference_candidate_analysis"].status == aa.AVAILABLE


def test_two_transcripts_are_not_two_protein_isoforms(tmp_path):
    """The distinction the counts exist to preserve."""
    same = _write_run(tmp_path / "a", exons=3, proteins={"P1": "MKVLA"},
                      transcripts=[("T1", "P1"), ("T2", "P1")], boundary_rows=2)
    differ = _write_run(tmp_path / "b", exons=3, proteins={"P1": "MKVLA", "P2": "MKVLAQ"},
                        boundary_rows=2)
    assert aa.prerequisites(same).transcript_model_count == 2
    assert aa.prerequisites(same).unique_protein_sequence_count == 1
    assert aa.prerequisites(differ).unique_protein_sequence_count == 2


# --------------------------------------------------------------------------- #
# What still blocks a ready run
# --------------------------------------------------------------------------- #
def test_a_missing_domain_table_blocks_readiness(tmp_path):
    run = _write_run(tmp_path, exons=3, proteins={"P1": "MKVLA"}, boundary_rows=2)
    (run / aa.CORE / "domain_features.tsv").unlink()
    manifest = aa.build_manifest(run)
    assert manifest.by_name()["domain_architecture"].status == aa.TECHNICALLY_MISSING
    assert not manifest.ready


def test_a_missing_boundary_table_blocks_readiness_when_boundaries_exist(tmp_path):
    run = _write_run(tmp_path, exons=3, proteins={"P1": "MKVLA"}, boundary_rows=2)
    (run / aa.CORE / "exon_domain_boundary_distances.tsv").unlink()
    manifest = aa.build_manifest(run)
    assert manifest.by_name()["boundary_analysis"].status == aa.TECHNICALLY_MISSING
    assert not manifest.ready


def test_an_unfinished_cluster_leaves_an_applicable_analysis_pending(tmp_path):
    run = _write_run(tmp_path, exons=3, proteins={"P1": "MKVLA"}, boundary_rows=2,
                     cluster=False)
    manifest = aa.build_manifest(run)
    assert manifest.by_name()["boundary_analysis"].status == aa.PENDING
    assert not manifest.ready


def test_an_output_older_than_the_models_is_stale_and_blocks_readiness(tmp_path):
    import os
    import time

    run = _write_run(tmp_path, exons=3, proteins={"P1": "MKVLA"}, boundary_rows=2)
    # Re-collecting the models rewrites the primary FASTA; the analyses that described the
    # previous proteins are no longer current results.
    later = time.time() + 60
    faa = run / aa.CORE / "proteins_primary.faa"
    faa.write_text(">P1\nMKVLA\n", encoding="utf-8")
    os.utime(faa, (later, later))

    manifest = aa.build_manifest(run)
    names = manifest.by_name()
    assert names["boundary_analysis"].status == aa.STALE
    assert names["domain_architecture"].status == aa.STALE
    assert not manifest.ready
    assert "older than the data it summarises" in names["boundary_analysis"].user_message


def test_a_not_applicable_analysis_cannot_become_stale(tmp_path):
    import os
    import time

    # There is nothing to recompute for a single-exon protein, so a newer FASTA does not
    # make the answer out of date.
    run = _write_run(tmp_path, exons=1, proteins={"P1": "MKV"})
    later = time.time() + 60
    faa = run / aa.CORE / "proteins_primary.faa"
    faa.write_text(">P1\nMKV\n", encoding="utf-8")
    os.utime(faa, (later, later))
    assert aa.build_manifest(run).by_name()["boundary_analysis"].status == aa.NOT_APPLICABLE


def test_an_unfinished_cluster_does_not_make_a_single_exon_gene_pending(tmp_path):
    # Applicability is settled before the cluster is consulted, otherwise a single-exon
    # gene would wait for an annotation that could never change the answer.
    run = _write_run(tmp_path, exons=1, proteins={"P1": "MKV"}, cluster=False)
    state = aa.build_manifest(run).by_name()["boundary_analysis"]
    assert state.status == aa.NOT_APPLICABLE


def test_a_run_without_the_core_layout_keeps_its_own_contract(tmp_path):
    # FGFR2 event-pipeline runs store models elsewhere. Reading the core tables there gives
    # zero of everything, and zero coding exons must not be reported as "not applicable"
    # for a run whose analyses are complete.
    empty = tmp_path / "fgfr2_style"
    (empty / "results").mkdir(parents=True)
    assert not aa.has_core_tables(empty)
    assert aa.annotate_dataset_model({"figures": []}, empty) == {"figures": []}


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def test_a_figure_for_an_inapplicable_analysis_is_not_registered(tmp_path):
    run = _write_run(tmp_path, exons=1, proteins={"P1": "MKV"})
    model = {"figures": [
        {"figure_id": "main_sp_boundary_class_summary", "status": "unavailable"},
        {"figure_id": "exploratory_candidate_ranking", "status": "unavailable"},
        {"figure_id": "transcript_exon_structure", "status": "available"},
    ]}
    ids = [f["figure_id"] for f in aa.annotate_dataset_model(model, run)["figures"]]
    assert ids == ["transcript_exon_structure"]


def test_a_genuinely_pending_figure_stays_visible(tmp_path):
    # A gap must remain reported; only inapplicable figures are dropped.
    run = _write_run(tmp_path, exons=3, proteins={"P1": "MKVLA"}, boundary_rows=2,
                     cluster=False)
    model = {"figures": [{"figure_id": "main_sp_boundary_class_summary",
                          "status": "pending_cluster"}]}
    ids = [f["figure_id"] for f in aa.annotate_dataset_model(model, run)["figures"]]
    assert ids == ["main_sp_boundary_class_summary"]


def test_an_available_figure_is_never_dropped(tmp_path):
    run = _write_run(tmp_path, exons=1, proteins={"P1": "MKV"})
    model = {"figures": [{"figure_id": "candidate_thing", "status": "available"}]}
    assert len(aa.annotate_dataset_model(model, run)["figures"]) == 1


# --------------------------------------------------------------------------- #
# Identity carried with every model
# --------------------------------------------------------------------------- #
def test_the_index_version_changes_when_the_indices_change(tmp_path):
    run = _write_run(tmp_path, exons=1, proteins={"P1": "MKV"})
    indices = run / "website_indices"
    indices.mkdir(parents=True, exist_ok=True)
    (indices / "a.json").write_text("{}", encoding="utf-8")
    before = aa.index_version(run)
    (indices / "b.json").write_text('{"x":1}', encoding="utf-8")
    assert aa.index_version(run) != before


# --------------------------------------------------------------------------- #
# Real runs
# --------------------------------------------------------------------------- #
@_needs(MC1R)
def test_mc1r_gallus_gallus_has_one_coding_exon_and_no_internal_boundary():
    pre = aa.prerequisites(MC1R)
    assert pre.coding_exon_count == 1
    assert pre.internal_coding_exon_boundary_count == 0
    assert pre.unique_protein_sequence_count == 1


@_needs(MC1R)
def test_mc1r_boundary_is_not_applicable_and_the_run_is_ready():
    manifest = aa.build_manifest(MC1R)
    names = manifest.by_name()
    assert names["exon_map"].status == aa.AVAILABLE
    assert names["domain_architecture"].status == aa.AVAILABLE
    assert names["boundary_analysis"].status == aa.NOT_APPLICABLE
    assert names["protein_isoform_comparison"].status == aa.NOT_APPLICABLE
    assert names["protein_difference_candidate_analysis"].status == aa.NOT_APPLICABLE
    assert manifest.ready


@_needs(MC1R)
def test_mc1r_is_results_ready_and_not_post_interpro_incomplete():
    from exondomaincompare.framework.core_run_milestones import evaluate_core_run
    report = evaluate_core_run(MC1R)
    assert report["inferred_status"] == "results_ready"
    persisted = json.loads((MC1R / "status.json").read_text(encoding="utf-8"))
    assert persisted["status"] == "results_ready"
    assert not persisted.get("failed_reason")


@_needs(MC1R)
def test_mc1r_boundary_is_not_blocking_at_the_species_level():
    from exondomaincompare.framework.species_completion import build_species_completion
    record = build_species_completion(MC1R)["gallus_gallus"]
    assert record["boundary_analysis"] == "not_applicable"
    assert record["blocking_analyses"] == []
    assert record["complete"]


@_needs(HBA)
def test_hba_panthera_leo_keeps_boundary_available_but_candidates_not_applicable():
    pre = aa.prerequisites(HBA)
    assert pre.coding_exon_count == 3
    assert pre.internal_coding_exon_boundary_count == 2
    assert pre.unique_protein_sequence_count == 1

    names = aa.build_manifest(HBA).by_name()
    assert names["exon_map"].status == aa.AVAILABLE
    assert names["domain_architecture"].status == aa.AVAILABLE
    assert names["boundary_analysis"].status == aa.AVAILABLE
    assert names["protein_isoform_comparison"].status == aa.NOT_APPLICABLE
    assert names["protein_difference_candidate_analysis"].status == aa.NOT_APPLICABLE
    assert aa.build_manifest(HBA).ready


@_needs(HBA)
def test_hba_absent_candidates_do_not_make_the_run_incomplete():
    persisted = json.loads((HBA / "status.json").read_text(encoding="utf-8"))
    assert persisted["status"] == "results_ready"


@_needs(HBA)
def test_hba_keeps_its_validated_protein_and_transcript():
    rows = (HBA / aa.CORE / "protein_isoform_index.tsv").read_text(encoding="utf-8")
    assert "XP_042777615.1" in rows
    assert "XM_042921681.1" in rows
    assert "\t142\t" in rows


@_needs(NKD2)
def test_nkd2_panthera_onca_has_nine_internal_boundaries():
    pre = aa.prerequisites(NKD2)
    assert pre.coding_exon_count == 10
    assert pre.internal_coding_exon_boundary_count == 9
    assert aa.build_manifest(NKD2).by_name()["boundary_analysis"].status == aa.AVAILABLE


@_needs(AKT1)
def test_a_run_with_several_distinct_sequences_keeps_isoform_comparison_available():
    # The counter-case: eleven transcripts, three distinct proteins. Isoform comparison and
    # candidate analysis are genuinely applicable here.
    pre = aa.prerequisites(AKT1)
    assert pre.transcript_model_count > pre.unique_protein_sequence_count
    assert pre.unique_protein_sequence_count >= 2
    names = aa.build_manifest(AKT1).by_name()
    assert names["protein_isoform_comparison"].status == aa.AVAILABLE
    assert names["protein_difference_candidate_analysis"].status == aa.AVAILABLE


# --------------------------------------------------------------------------- #
# The served model
# --------------------------------------------------------------------------- #
@_needs(MC1R)
def test_the_served_model_carries_identity_and_the_manifest():
    import main

    model = main.current_dataset_model(dataset=f"run:{MC1R.name}")
    assert model["dataset_id"] == f"run:{MC1R.name}"
    assert model["run_id"] == MC1R.name
    assert model["index_version"], "the model must name the index version it came from"

    manifest = model["analysis_availability"]
    assert manifest["ready"] is True
    states = {a["analysis_name"]: a for a in manifest["analyses"]}
    assert states["boundary_analysis"]["status"] == aa.NOT_APPLICABLE
    assert states["boundary_analysis"]["reason_code"] == aa.SINGLE_CODING_EXON
    assert states["boundary_analysis"]["blocks_results_ready"] is False


@_needs(MC1R)
def test_the_manifest_row_carries_every_field_the_downloads_table_shows():
    import main

    model = main.current_dataset_model(dataset=f"run:{MC1R.name}")
    row = {a["analysis_name"]: a
           for a in model["analysis_availability"]["analyses"]}["boundary_analysis"]
    for field in ("analysis_name", "status", "prerequisite_name", "prerequisite_count",
                  "reason_code", "user_message"):
        assert field in row, field
    assert row["prerequisite_name"] == "internal_coding_exon_boundary_count"
    assert row["prerequisite_count"] == 0


@_needs(MC1R)
def test_the_served_model_registers_no_empty_boundary_or_candidate_figure():
    import main

    model = main.current_dataset_model(dataset=f"run:{MC1R.name}")
    figures = model.get("figures") or []
    figures = figures if isinstance(figures, list) else figures.get("figures") or []
    assert figures, "the gallery must still hold the figures that were produced"
    for figure in figures:
        assert figure.get("status") == "available", figure.get("figure_id")
        fid = figure.get("figure_id", "").lower()
        assert "boundary" not in fid and "candidate" not in fid, fid


@_needs(MC1R)
def test_the_single_exon_run_is_results_ready_over_the_status_endpoint():
    import main

    status = main.dataset_status(f"run:{MC1R.name}")
    assert status["status"] == "results_ready"
    assert status["available_views"]["exon_domain_boundaries"] is False
    assert status["available_views"]["domain_architecture"] is True


@_needs(HBA)
def test_a_gene_with_internal_boundaries_keeps_the_boundary_view():
    import main

    status = main.dataset_status(f"run:{HBA.name}")
    assert status["status"] == "results_ready"
    assert status["available_views"]["exon_domain_boundaries"] is True


# --------------------------------------------------------------------------- #
# Availability is stamped as an index is served, so runs need no rebuild
# --------------------------------------------------------------------------- #
@_needs(MC1R)
def test_the_boundary_index_reports_not_applicable_without_being_rebuilt():
    import main

    index = main.current_shared_index("exon-domain-boundaries", dataset=f"run:{MC1R.name}")
    assert index["status"] == aa.NOT_APPLICABLE
    assert index["available"] is False
    block = index["availability"]
    assert block["state"] == aa.NOT_APPLICABLE
    assert block["reason_code"] == aa.SINGLE_CODING_EXON
    assert block["badge"] == "No internal boundaries"
    assert block["reason"] == aa.MESSAGES[aa.SINGLE_CODING_EXON]
    assert block["prerequisite_count"] == 0
    # The file on disk still carries the pipeline's own wording: the correction is in what
    # is served, not a rewrite of a produced index.
    on_disk = json.loads(
        (MC1R / "website_indices" / "exon_domain_boundary_index.json")
        .read_text(encoding="utf-8"))
    assert on_disk.get("availability") is None


@_needs(HBA)
def test_a_produced_boundary_index_keeps_its_own_verdict():
    import main

    index = main.current_shared_index("exon-domain-boundaries", dataset=f"run:{HBA.name}")
    assert index["available"] is True
    assert index["status"] == "available"
    assert index["n_boundaries"] == 2


@_needs(MC1R)
def test_the_isoform_alignment_index_uses_the_agreed_one_protein_wording():
    import main

    index = main.current_msa(dataset=f"run:{MC1R.name}")
    block = index["availability"]
    assert block["state"] == aa.NOT_APPLICABLE
    assert block["label"] == "Protein isoform comparison"
    assert block["reason"] == aa.MESSAGES[aa.SINGLE_PROTEIN_SEQUENCE]
    assert block["badge"] == "One protein sequence"


@_needs(AKT1)
def test_a_real_isoform_comparison_is_not_marked_not_applicable():
    import main

    index = main.current_msa(dataset=f"run:{AKT1.name}")
    assert index["available"] is True
    assert index.get("status") != aa.NOT_APPLICABLE


def test_serving_an_index_for_a_non_core_run_leaves_it_untouched(tmp_path):
    import main

    payload = {"available": False, "status": "pending_cluster"}
    ds = {"run_base": tmp_path}          # no core tables at all
    assert main._with_availability(payload, "msa_index.json", ds) == payload


# --------------------------------------------------------------------------- #
# Regressions for the validated biology
# --------------------------------------------------------------------------- #
def test_the_validated_fgfr2_freeze_keeps_its_own_availability_contract():
    import main

    model = main.current_dataset_model(dataset="example")
    assert model["dataset_id"] == "example"
    # The freeze is read-only and is not re-judged by the shared manifest.
    assert "analysis_availability" not in model
    figures = model.get("figures") or []
    figures = figures if isinstance(figures, list) else figures.get("figures") or []
    assert len(figures) > 100, "the validated gallery must be untouched"


@pytest.mark.parametrize("run_name,gene", [
    ("2026-07-26_2157_fgfr1_gallus_mus_core_pilot", "FGFR1"),
    ("2026-07-16_1642_tp53_human_core_pilot", "TP53"),
])
def test_multi_exon_multi_isoform_runs_keep_every_analysis_available(run_name, gene):
    run = RUNS / run_name
    if not run.is_dir():
        pytest.skip(f"{run_name} is not present")
    pre = aa.prerequisites(run)
    assert pre.internal_coding_exon_boundary_count > 0, gene
    assert pre.unique_protein_sequence_count >= 2, gene
    for state in aa.build_manifest(run).analyses:
        assert state.status == aa.AVAILABLE, f"{gene}: {state.analysis_name}"


def test_the_fgfr2_validated_freeze_is_untouched():
    freeze = ROOT / "results" / "final_30_until_interpro_prepare"
    if not freeze.is_dir():
        pytest.skip("no freeze directory in this checkout")
    proc = subprocess.run(["git", "status", "--porcelain", "--", str(freeze)],
                          capture_output=True, text=True, cwd=ROOT)
    assert proc.stdout.strip() == "", f"the validated freeze changed:\n{proc.stdout}"


# --------------------------------------------------------------------------- #
# Frontend freshness, driven through node against the real modules
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(which("node") is None, reason="node is required")
def test_the_frontend_data_layer_keeps_dataset_identity():
    proc = subprocess.run([which("node"), str(ROOT / "tests"
                                              / "check_frontend_freshness.mjs")],
                          capture_output=True, text=True, cwd=ROOT)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "FAIL" not in proc.stdout, proc.stdout


# --------------------------------------------------------------------------- #
# Wiring that has to stay in place
# --------------------------------------------------------------------------- #
APP = ROOT / "webapp" / "frontend" / "src" / "App.jsx"
WORKFLOW = ROOT / "webapp" / "frontend" / "src" / "pages" / "RunWorkflowPage.jsx"


def test_switching_datasets_clears_the_previous_datasets_data_first():
    body = APP.read_text(encoding="utf-8").split("const selectDataset")[1] \
        .split("const refreshAll")[0]
    clear = body.index("setDatasetModel(null)")
    load = body.index("loadInto(target, epoch)")
    assert clear < load, "the previous model must be cleared before the new one is loaded"
    assert "setDatasetInfo(null)" in body
    assert "setGeneTarget(null)" in body


def test_a_dataset_load_is_applied_only_when_it_is_still_current():
    body = APP.read_text(encoding="utf-8").split("const loadInto")[1] \
        .split("const loadDatasets")[0]
    assert "epoch !== epochRef.current" in body
    assert "payloadMatchesDataset(info, target)" in body
    assert "payloadMatchesDataset(model, target)" in body


def test_a_transient_model_load_is_retried_before_the_gallery_is_shown():
    body = APP.read_text(encoding="utf-8").split("const loadInto")[1] \
        .split("const loadDatasets")[0]
    assert "firstModel" in body
    assert body.count("client.datasetModel()") == 2
    assert "signal?.aborted" in body


def test_the_gallery_rejects_a_scope_from_another_dataset():
    gallery = (ROOT / "webapp" / "frontend" / "src" / "pages"
               / "FigureGallery.jsx").read_text(encoding="utf-8")
    assert 'speciesList.some((s) => s.id === scope)' in gallery
    assert '? scope : "comparative"' in gallery


def test_a_superseded_dataset_load_is_aborted():
    body = APP.read_text(encoding="utf-8")
    assert "abortRef.current?.abort()" in body
    assert "new AbortController()" in body


def test_polling_stops_only_at_a_stable_state():
    body = APP.read_text(encoding="utf-8")
    assert "TERMINAL_RUN_STATES.has(datasetInfo?.status)" in body
    # A state change must refresh the registry too, so My Runs and the selector move with
    # the page rather than needing a browser refresh.
    poll = body.split("const timer = window.setInterval")[1].split("POLL_INTERVAL_MS")[0]
    assert "loadDatasets()" in poll
    assert "setRefreshNonce" in poll


def test_a_new_run_is_shown_before_the_registry_is_asked_again():
    body = WORKFLOW.read_text(encoding="utf-8").split("async function handleStart")[1] \
        .split("\n  }")[0]
    insert = body.index("insertRun(res)")
    reload = body.index("loadRuns({ preserveIds })")
    assert insert < reload, "the created run must appear without waiting for the registry"


def test_a_provisioning_run_is_discovered_without_a_manual_refresh():
    body = WORKFLOW.read_text(encoding="utf-8").split("async function handleStart")[1] \
        .split("\n  }")[0]
    assert "RUN_DISCOVERY_ATTEMPTS" in body
    assert "knownIds.has(run.run_id)" in body
    assert "setTimeout(resolve, RUN_DISCOVERY_POLL_MS)" in body
    assert "preserveIds" in body


def test_the_selected_run_is_polled_until_it_is_stable():
    body = WORKFLOW.read_text(encoding="utf-8")
    assert "isActiveRunState(status)" in body, \
        "polling must cover every non-terminal state, not only the local running ones"


def test_a_detail_refresh_updates_the_matching_run_card():
    body = WORKFLOW.read_text(encoding="utf-8").split("const loadDetail")[1] \
        .split("useEffect(() => { loadRuns();")[0]
    assert "ref.summary" in body
    assert "run.run_id === runId" in body
    assert "runsRef.current = next" in body


def test_a_resolved_analysis_loses_its_pending_navigation_badge():
    body = APP.read_text(encoding="utf-8")
    assert 'isResolvedAnalysis(datasetModel, "boundary_analysis")' in body
