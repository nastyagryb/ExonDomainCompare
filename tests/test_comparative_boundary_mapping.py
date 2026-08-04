#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

TWO_SPECIES_RUN = ROOT / "runs" / "2026-07-26_2157_fgfr1_gallus_mus_core_pilot"
ALIGNMENT = (TWO_SPECIES_RUN
             / "results/generic_gene_analysis/msa/primaries_msa.aln.faa")

pytestmark = pytest.mark.skipif(
    not TWO_SPECIES_RUN.is_dir(),
    reason=f"reference two-species run not present: {TWO_SPECIES_RUN.name}")


@pytest.fixture(scope="module")
def mapped():
    from exondomaincompare.shared_gene_analysis.msa_coordinates import (
        annotate_boundaries_with_columns, build_msa_coordinate_map,
    )
    from exondomaincompare.shared_gene_analysis.protein_coordinate_model import build_models_for_run
    idx = build_models_for_run(TWO_SPECIES_RUN)
    coord_map = build_msa_coordinate_map(ALIGNMENT)
    report = annotate_boundaries_with_columns(idx["models"], coord_map)
    return idx, coord_map, report


# --------------------------------------------------------------------------- #
# 11. MSA coordinate mapping
# --------------------------------------------------------------------------- #
def test_the_cross_species_alignment_holds_one_primary_per_species(mapped):
    _, coord_map, _ = mapped
    assert coord_map["available"], coord_map.get("reason")
    species = {s["species_id"] for s in coord_map["species"]}
    assert species == {"gallus_gallus", "mus_musculus"}
    assert len({s["protein_id"] for s in coord_map["species"]}) == 2


def test_the_mapping_reproduces_each_protein_length(mapped):
    idx, coord_map, _ = mapped
    by_species = {m["species_id"]: m for m in idx["models"]}
    for s in coord_map["species"]:
        assert s["protein_length"] == by_species[s["species_id"]]["protein_length"], (
            f"{s['species_id']}: the alignment contains {s['protein_length']} residues "
            f"but the model reports {by_species[s['species_id']]['protein_length']} aa; "
            f"a mismatch means a different protein was aligned")


def test_every_boundary_of_every_species_receives_a_column(mapped):
    _, _, report = mapped
    for sid, r in report.items():
        assert r["protein_matches_alignment"], (
            f"{sid}: model protein {r['protein_id']} is not the aligned protein "
            f"{r['aligned_protein_id']}; mapping through a different protein would "
            f"place boundaries at plausible but wrong columns")
        assert r["boundaries_mapped"] == r["boundaries_total"] > 0, r


def test_native_positions_and_columns_differ_where_the_alignment_has_gaps(mapped):
    idx, _, _ = mapped
    gallus = next(m for m in idx["models"] if m["species_id"] == "gallus_gallus")
    shifted = [b for b in gallus["exon_boundaries"]
               if b["msa_column"] is not None
               and b["msa_column"] != b.get("protein_position")]
    assert shifted, (
        "Gallus is shorter than Mus and carries alignment gaps, so at least some "
        "boundaries must sit at a column different from their native position")


def test_a_single_species_run_reports_no_mapping_instead_of_failing():
    from exondomaincompare.shared_gene_analysis.msa_coordinates import build_msa_coordinate_map
    cm = build_msa_coordinate_map(ROOT / "does" / "not" / "exist.faa")
    assert cm["available"] is False
    assert "single-species" in cm["reason"]
    assert "/Users/" not in cm["reason"], "must not leak a machine-local path"


# --------------------------------------------------------------------------- #
# 14, 17. comparable-boundary grouping: evidence priority, never exon rank
# --------------------------------------------------------------------------- #
def test_comparable_groups_exist_and_name_their_evidence(mapped):
    import exondomaincompare.shared_gene_analysis.boundary_dashboard as bd
    idx, _, _ = mapped
    groups = bd.match_comparable_boundaries(idx["models"])
    assert groups, "two fully processed species must yield comparable boundary groups"
    for g in groups:
        assert g["mapping_method"] in bd.MAPPING_METHOD_PRIORITY
        assert g["mapping_method"] != "exon_rank", (
            "exon rank is descriptive only and must never be a grouping key")
        assert g["mapping_status"] in bd.COMPARABLE_STATES, g["mapping_status"]
        species = {p["species_id"] for p in g["per_species_native_positions"]}
        assert len(species) >= 2, f"{g['comparable_boundary_group_id']} is single-species"


def test_grouping_is_not_by_exon_rank(mapped):
    import exondomaincompare.shared_gene_analysis.boundary_dashboard as bd
    idx, _, _ = mapped
    groups = bd.match_comparable_boundaries(idx["models"])
    for g in groups:
        members = g["per_species_native_positions"]
        cols = {p["msa_column"] for p in members if p["msa_column"] is not None}
        if len(cols) > 1:
            assert g["mapping_status"] == "tentative", (
                f"{g['comparable_boundary_group_id']} merges columns {sorted(cols)} and "
                f"must therefore be tentative, not {g['mapping_status']}")
            assert max(cols) - min(cols) <= bd.NEAR_COLUMN_TOLERANCE


def test_a_genomic_exon_group_is_not_used_as_cross_species_evidence():
    import exondomaincompare.shared_gene_analysis.boundary_dashboard as bd
    # Two species whose exon groups are distinct genomic hashes, but whose boundaries
    # share an alignment column. Grouping must succeed via the column, not the hash.
    def model(sid, group, col):
        return {
            "species_id": sid, "protein_id": f"P_{sid}", "status": "available",
            "exons": [{"id": f"{sid}:e1", "tooltip": {"shared_exon_group": f"{group}_a"}},
                      {"id": f"{sid}:e2", "tooltip": {"shared_exon_group": f"{group}_b"}}],
            "exon_boundaries": [{
                "id": f"{sid}:b1", "left_exon_id": f"{sid}:e1",
                "right_exon_id": f"{sid}:e2", "protein_position": 100,
                "boundary_position_aa": 100, "signed_distance": -2,
                "boundary_class": "near_domain_edge", "msa_column": col}],
        }
    groups = bd.match_comparable_boundaries(
        [model("gallus_gallus", "SEG_gal", 120), model("mus_musculus", "SEG_mus", 120)])
    assert len(groups) == 1, (
        "distinct genomic exon-group hashes must not prevent grouping by alignment "
        "column; treating them as the first-choice cross-species key discarded "
        "every group")
    assert groups[0]["mapping_method"] == "msa_aligned_position"


def test_a_shared_exon_group_is_still_used_when_it_really_spans_species():
    import exondomaincompare.shared_gene_analysis.boundary_dashboard as bd

    def model(sid):
        return {
            "species_id": sid, "protein_id": f"P_{sid}", "status": "available",
            "exons": [{"id": f"{sid}:e1", "tooltip": {"shared_exon_group": "SEG_x"}},
                      {"id": f"{sid}:e2", "tooltip": {"shared_exon_group": "SEG_y"}}],
            "exon_boundaries": [{
                "id": f"{sid}:b1", "left_exon_id": f"{sid}:e1",
                "right_exon_id": f"{sid}:e2", "protein_position": 100,
                "boundary_position_aa": 100, "signed_distance": 0,
                "boundary_class": "exact_domain_edge", "msa_column": None}],
        }
    groups = bd.match_comparable_boundaries([model("a"), model("b")])
    assert len(groups) == 1
    assert groups[0]["mapping_method"] == "shared_exon_group"


def test_boundaries_without_any_evidence_are_not_grouped():
    import exondomaincompare.shared_gene_analysis.boundary_dashboard as bd

    def model(sid, pos):
        return {
            "species_id": sid, "protein_id": f"P_{sid}", "status": "available",
            "exons": [{"id": f"{sid}:e1", "tooltip": {}},
                      {"id": f"{sid}:e2", "tooltip": {}}],
            "exon_boundaries": [{
                "id": f"{sid}:b1", "left_exon_id": f"{sid}:e1",
                "right_exon_id": f"{sid}:e2", "protein_position": pos,
                "boundary_position_aa": pos, "signed_distance": 0,
                "boundary_class": "inside_domain", "msa_column": None}],
        }
    assert bd.match_comparable_boundaries([model("a", 100), model("b", 100)]) == [], (
        "identical native positions are not evidence of comparability")


# --------------------------------------------------------------------------- #
# 18-19. matrix and signed distances
# --------------------------------------------------------------------------- #
def test_the_matrix_has_one_row_per_species_and_one_cell_per_group(mapped):
    import exondomaincompare.shared_gene_analysis.boundary_dashboard as bd
    idx, _, _ = mapped
    groups = bd.match_comparable_boundaries(idx["models"])
    matrix = bd.build_boundary_matrix(idx["models"], groups)
    assert len(matrix) == 2
    for row in matrix:
        assert len(row["cells"]) == len(groups)


def test_matrix_signed_distances_match_the_underlying_boundaries(mapped):
    import exondomaincompare.shared_gene_analysis.boundary_dashboard as bd
    idx, _, _ = mapped
    groups = bd.match_comparable_boundaries(idx["models"])
    by_species = {m["species_id"]: m for m in idx["models"]}
    truth = {sid: {b["id"]: b.get("signed_distance") for b in m["exon_boundaries"]}
             for sid, m in by_species.items()}

    for g in groups:
        for member in g["per_species_native_positions"]:
            sid, bid = member["species_id"], member["boundary_id"]
            assert member["signed_distance"] == truth[sid][bid], (
                f"{g['comparable_boundary_group_id']}/{sid}: comparative view reports "
                f"{member['signed_distance']} but the boundary record says "
                f"{truth[sid][bid]}")


def test_consistency_statistics_cover_every_group_without_inventing_values(mapped):
    import exondomaincompare.shared_gene_analysis.boundary_dashboard as bd
    idx, _, _ = mapped
    groups = bd.match_comparable_boundaries(idx["models"])
    stats = bd.boundary_position_consistency(idx["models"], groups)
    assert len(stats) == len(groups)
    ids = {s["comparable_boundary_group_id"] for s in stats}
    assert ids == {g["comparable_boundary_group_id"] for g in groups}
    for s in stats:
        assert 0 < s["mapping_coverage"] <= 1
        assert s["species_with_mapped_boundary"] >= 2


def test_the_contract_is_published_in_the_run_index():
    import json
    model = json.loads(
        (TWO_SPECIES_RUN / "website_indices/generic/protein_coordinate_model.json")
        .read_text())
    ms = (model.get("boundary_dashboard") or {}).get("multi_species") or {}
    assert ms.get("available") is True, "the published index still reports no comparison"
    assert ms.get("comparable_boundary_groups"), "no comparable groups were published"
    assert ms.get("boundary_matrix"), "no matrix was published"
    assert model.get("msa_coordinate_map", {}).get("available") is True
    assert model["boundary_dashboard"]["page_mode"] == "generic_multi_species_results_ready"
