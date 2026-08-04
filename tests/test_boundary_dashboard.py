from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SGA = ROOT / "scripts" / "shared_gene_analysis"
if str(SGA) not in sys.path:
    sys.path.insert(0, str(SGA))

from exondomaincompare.shared_gene_analysis import boundary_dashboard as bd  # noqa: E402
from exondomaincompare.shared_gene_analysis import protein_coordinate_model as pcm  # noqa: E402

# Both dashboards are built from committed fixtures rather than from the live run
# registry: the Core tables are verbatim copies of the runs they were taken from
# (see each fixture's `fixture_provenance.json`), so the scientific assertions below
# are unchanged, but deleting or re-running a production run can no longer decide
# whether this file passes. `tests/fixtures/build_boundary_fixtures.py` regenerates
# them. Paths are repository-relative in the served dashboard because the project
# root is passed explicitly.
FIXTURES = Path(__file__).resolve().parent / "fixtures"
TP53_DANIO_RUN = FIXTURES / "tp53_danio_boundary_run"
FGFR1_RUN = FIXTURES / "fgfr1_gallus_boundary_run"


@pytest.fixture(scope="module")
def tp53_index():
    return pcm.build_models_for_run(TP53_DANIO_RUN, project_root=ROOT)


@pytest.fixture(scope="module")
def fgfr1_index():
    # Post repair-and-simplification pass the FGFR1 run is a completed cluster
    # round-trip (results-ready). Its pre-cluster/pending behaviour is exercised
    # separately by `pending_index` so the honest-pending code paths stay covered.
    return pcm.build_models_for_run(FGFR1_RUN, project_root=ROOT)


def _synthetic_pending_index():
    bnds = [{
        "id": f"b{i}", "protein_position": 40 * i + 10,
        "left_exon_label": f"E{i}", "right_exon_label": f"E{i + 1}",
        "label": f"E{i} → E{i + 1}",
        "boundary_class": None, "signed_distance": None,
        "nearest_domain_id": None, "mapping_status": "pending_cluster",
    } for i in range(1, 17)]
    model = {
        "species_id": "gallus_gallus", "scientific_name": "Gallus gallus",
        "protein_id": "NP_990841.2", "transcript_id": "NM_205510.2",
        "protein_length": 817, "near_edge_threshold_aa": 5,
        "coordinate_system": "protein_1_based_inclusive",
        "status": "pending_cluster",
        "exons": [{"id": f"e{i}"} for i in range(1, 18)],
        "representative_domains": [], "candidate_regions": [],
        "exon_boundaries": bnds,
    }
    idx = {"models": [model], "gene_symbol": "FGFR1",
           "coordinate_system": "protein_1_based_inclusive"}
    idx["boundary_dashboard"] = bd.build_boundary_dashboard(idx)
    return idx


@pytest.fixture(scope="module")
def pending_index():
    return _synthetic_pending_index()


# --------------------------------------------------------------------------- #
# Part 1 — page-mode resolution (data-driven, no gene-symbol branching)
# --------------------------------------------------------------------------- #
def test_page_mode_single_species(tp53_index):
    assert bd.resolve_page_mode(tp53_index) == bd.PAGE_SINGLE


def test_page_mode_pending(pending_index):
    assert bd.resolve_page_mode(pending_index) == bd.PAGE_PENDING


def test_page_mode_fgfr1_results_ready_after_repair(fgfr1_index):
    # Regression guard for the post-cluster discovery repair (Part A): the real
    # FGFR1 run is a completed round-trip and must route to the single-species
    # results-ready page, never to the pending page.
    assert bd.resolve_page_mode(fgfr1_index) == bd.PAGE_SINGLE


def test_page_mode_unavailable_empty():
    assert bd.resolve_page_mode({"models": []}) == bd.PAGE_UNAVAILABLE
    assert bd.resolve_page_mode({}) == bd.PAGE_UNAVAILABLE


def test_page_mode_validated_event_takes_priority(tp53_index):
    # FGFR2 freeze protection: a validated event never routes to the generic page.
    assert bd.resolve_page_mode(tp53_index, event_layer_type="validated") == bd.PAGE_VALIDATED_EVENT


def test_page_mode_multi_species_requires_two_available():
    idx = {"models": [
        {"species_id": "a", "status": "available"},
        {"species_id": "b", "status": "available"},
    ]}
    assert bd.resolve_page_mode(idx) == bd.PAGE_MULTI
    idx2 = {"models": [
        {"species_id": "a", "status": "available"},
        {"species_id": "b", "status": "pending_cluster"},
    ]}
    assert bd.resolve_page_mode(idx2) == bd.PAGE_SINGLE


# --------------------------------------------------------------------------- #
# Part 4 / 8 — TP53 single-species dashboard data
# --------------------------------------------------------------------------- #
def test_tp53_dashboard_header_and_counts(tp53_index):
    dash = tp53_index["boundary_dashboard"]
    assert dash["page_mode"] == bd.PAGE_SINGLE
    ss = dash["single_species"]
    h = ss["header"]
    assert h["gene"] == "TP53"
    assert h["scientific_name"] == "Danio rerio"
    assert h["protein_id"] == "NP_001258749.1"
    assert h["transcript_id"] == "NM_001271820.1"
    assert h["protein_length"] == 374
    assert h["n_coding_exons"] == 10
    assert h["n_internal_boundaries"] == 9
    assert h["near_edge_threshold_aa"] == 5
    assert h["status_badge"] == "results_ready"


def test_tp53_class_summary_5_inside_4_outside(tp53_index):
    s = tp53_index["boundary_dashboard"]["single_species"]["summary"]
    assert s["total"] == 9
    assert s["inside_domain"] == 5
    assert s["outside_annotated_domains"] == 4
    # No fabricated exact/near values.
    assert s["exact_domain_edge"] == 0
    assert s["near_domain_edge"] == 0
    assert s["unavailable_or_uncertain"] == 0


def test_tp53_inspection_cases_are_real_and_selectable(tp53_index):
    dash = tp53_index["boundary_dashboard"]
    cases = dash["single_species"]["inspection_cases"]
    assert cases, "TP53 should surface at least one real inspection case"
    # every inspection case must reference an existing boundary (globally selectable)
    model = tp53_index["models"][0]
    bids = {b["id"] for b in model["exon_boundaries"]}
    for c in cases:
        assert c["boundary_id"] in bids
        assert c["kind"] in {
            "large_domain_edge_distance", "mapping_requires_inspection",
            "representative_annotation_unavailable", "incomplete_evidence",
            "candidate_associated",
        }
    # TP53 has candidate-overlapping and large-distance cases, no error kinds
    kinds = {c["kind"] for c in cases}
    assert "candidate_associated" in kinds or "large_domain_edge_distance" in kinds


# --------------------------------------------------------------------------- #
# Part 11 — FGFR1 honest pending state (no fabricated classifications)
# --------------------------------------------------------------------------- #
def test_pending_dashboard(pending_index):
    dash = pending_index["boundary_dashboard"]
    assert dash["page_mode"] == bd.PAGE_PENDING
    ss = dash["single_species"]
    h = ss["header"]
    assert h["protein_id"] == "NP_990841.2"
    assert h["protein_length"] == 817
    assert h["n_coding_exons"] == 17
    assert h["n_internal_boundaries"] == 16
    assert h["status_badge"] == "pending_cluster"
    assert "pending" in h["representative_domain_source"].lower()
    # No classifications, no inspection cases before the cluster round-trip.
    s = ss["summary"]
    assert s["total"] == 16
    assert s["unavailable_or_uncertain"] == 16
    assert s["inside_domain"] == 0 and s["outside_annotated_domains"] == 0
    assert ss["inspection_cases"] == []


def test_pending_no_fake_classifications_on_boundaries(pending_index):
    model = pending_index["models"][0]
    for b in model["exon_boundaries"]:
        assert b["boundary_class"] is None
        assert b["signed_distance"] is None
        assert b["nearest_domain_id"] is None


def test_fgfr1_results_ready_dashboard_after_repair(fgfr1_index):
    # Post-cluster discovery repair (Part A): the real FGFR1 run must now expose
    # real normalized domains and rebuilt boundary classifications — never fake,
    # never pending.
    dash = fgfr1_index["boundary_dashboard"]
    assert dash["page_mode"] == bd.PAGE_SINGLE
    model = fgfr1_index["models"][0]
    assert model["status"] == "available"
    assert model["representative_domains"], "repaired run must carry real domains"
    ss = dash["single_species"]
    assert ss["header"]["status_badge"] in ("results_ready", "partial")
    # at least one boundary is classified against a real domain (no fabrication,
    # but also no all-pending state)
    classified = [b for b in model["exon_boundaries"] if b["boundary_class"] is not None]
    assert classified, "boundary classifications must be rebuilt from real domains"


# --------------------------------------------------------------------------- #
# Parts 12–13 — multi-species contract + comparable-boundary matching
# --------------------------------------------------------------------------- #
def test_multi_species_contract_schema(tp53_index):
    ms = tp53_index["boundary_dashboard"]["multi_species"]
    for key in ("species_rows", "comparable_boundary_groups", "mapping_methods",
                "boundary_matrix", "distance_statistics", "inspection_cases"):
        assert key in ms
    assert ms["mapping_methods"] == bd.MAPPING_METHOD_PRIORITY
    # single species → no comparative results (honest empty)
    assert ms["available"] is False
    assert ms["comparable_boundary_groups"] == []
    assert len(ms["species_rows"]) == 1
    row = ms["species_rows"][0]
    for key in ("species_id", "scientific_name", "taxonomic_group", "primary_protein",
                "transcript", "protein_length", "analysis_status"):
        assert key in row


def test_mapping_method_priority_order():
    assert bd.MAPPING_METHOD_PRIORITY[0] == "shared_exon_group"
    assert bd.MAPPING_METHOD_PRIORITY[1] == "msa_aligned_position"
    assert bd.MAPPING_METHOD_PRIORITY[-1] == "exon_rank"  # secondary only


def test_match_single_species_returns_empty(tp53_index):
    assert bd.match_comparable_boundaries(tp53_index["models"]) == []


def _fake_two_species_shared_group():
    def _mk(sp, pid, group):
        return {
            "species_id": sp, "protein_id": pid, "status": "available",
            "scientific_name": sp, "near_edge_threshold_aa": 5,
            "exons": [
                {"id": f"{pid}:e1", "tooltip": {"shared_exon_group": group}},
                {"id": f"{pid}:e2", "tooltip": {"shared_exon_group": group}},
            ],
            "exon_boundaries": [{
                "id": f"{pid}:b1", "left_exon_id": f"{pid}:e1", "right_exon_id": f"{pid}:e2",
                "protein_position": 100, "signed_distance": 2,
                "boundary_class": "near_domain_edge", "msa_column": None,
            }],
        }
    return [_mk("sp_a", "PA", "SEG7"), _mk("sp_b", "PB", "SEG7")]


def test_match_shared_exon_group_high_confidence():
    groups = bd.match_comparable_boundaries(_fake_two_species_shared_group())
    assert len(groups) == 1
    g = groups[0]
    assert g["mapping_method"] == "shared_exon_group"
    assert g["mapping_status"] == "high_confidence_comparable"
    assert g["mapping_status"] in bd.COMPARABLE_STATES
    assert len(g["per_species_native_positions"]) == 2


def test_match_never_groups_by_exon_rank_alone():
    # Two species, same ordinal boundary name but NO shared exon group / MSA column.
    def _mk(sp, pid):
        return {
            "species_id": sp, "protein_id": pid, "status": "available",
            "exons": [{"id": f"{pid}:e1", "tooltip": {}}, {"id": f"{pid}:e2", "tooltip": {}}],
            "exon_boundaries": [{
                "id": f"{pid}:b1", "left_exon_id": f"{pid}:e1", "right_exon_id": f"{pid}:e2",
                "protein_position": 100, "signed_distance": 2, "msa_column": None,
            }],
        }
    assert bd.match_comparable_boundaries([_mk("a", "PA"), _mk("b", "PB")]) == []


def test_boundary_matrix_and_consistency_present_for_real_groups():
    idx = {"models": _fake_two_species_shared_group(),
           "gene_symbol": "X", "coordinate_system": "protein_1_based_inclusive"}
    ms = bd.build_multi_species_contract(idx)
    assert ms["available"] is True
    assert len(ms["comparable_boundary_groups"]) == 1
    assert len(ms["boundary_matrix"]) == 2
    assert len(ms["distance_statistics"]) == 1
    stat = ms["distance_statistics"][0]
    assert stat["metric_label"] == "Boundary-position consistency"
    assert stat["species_with_mapped_boundary"] == 2


# --------------------------------------------------------------------------- #
# uncertain / unavailable mapping states (Part 9)
# --------------------------------------------------------------------------- #
def test_uncertain_and_unavailable_inspection_states():
    model = {
        "status": "available", "protein_length": 100,
        "exons": [], "representative_domains": [], "candidate_regions": [],
        "exon_boundaries": [
            {"id": "b1", "protein_position": 10, "mapping_status": "unmapped",
             "boundary_class": "outside_annotated_domains", "signed_distance": 30,
             "nearest_edge_position": 40, "label": "E1 → E2"},
            {"id": "b2", "protein_position": 20, "mapping_status": "unavailable",
             "boundary_class": "unavailable_or_uncertain", "signed_distance": None,
             "label": "E2 → E3"},
        ],
    }
    kinds = {c["kind"] for c in bd.build_inspection_cases(model)}
    assert "mapping_requires_inspection" in kinds
    assert "representative_annotation_unavailable" in kinds
    assert bd._status_badge(model) == "uncertain_mapping"


# --------------------------------------------------------------------------- #
# Part 20 — caption generation
# --------------------------------------------------------------------------- #
def test_caption_generation_available(tp53_index):
    model = tp53_index["models"][0]
    cap = bd.generate_caption(model, "TP53")
    assert "Danio rerio" in cap["text"] and "TP53" in cap["text"]
    assert "NP_001258749.1" in cap["text"]
    assert "5 amino acids" in cap["text"]
    f = cap["fields"]
    assert f["gene"] == "TP53"
    assert f["near_edge_threshold_aa"] == 5
    assert f["coordinate_system"] == "protein_1_based_inclusive"
    assert f["analysis_status"] == "available"


def test_caption_generation_pending(pending_index):
    model = pending_index["models"][0]
    cap = bd.generate_caption(model, "FGFR1")
    assert "pending" in cap["text"].lower() or "will be calculated" in cap["text"].lower()
    assert cap["fields"]["analysis_status"] == "pending_cluster"


# --------------------------------------------------------------------------- #
# export / no personal paths in the served dashboard
# --------------------------------------------------------------------------- #
def test_no_personal_paths_in_dashboard(tp53_index, fgfr1_index):
    for idx in (tp53_index, fgfr1_index):
        blob = json.dumps(idx["boundary_dashboard"])
        assert "/Users/" not in blob
        assert str(ROOT) not in blob
