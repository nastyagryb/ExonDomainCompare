"""Tests for the shared protein-coordinate contract, validator and boundary classifier.

Covers milestone Part 1 (coordinate model + validation), Part 5 (boundary-class
priority, signed-distance, exact/near/inside/outside/unavailable) and Part 16
(coordinate ranges, single-species schema, TP53 post-cluster, no absolute paths).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SGA = ROOT / "scripts" / "shared_gene_analysis"
if str(SGA) not in sys.path:
    sys.path.insert(0, str(SGA))

from exondomaincompare.shared_gene_analysis import boundary_classification as bc  # noqa: E402
from exondomaincompare.shared_gene_analysis import protein_coordinate_model as pcm  # noqa: E402
from exondomaincompare.shared_gene_analysis import validate_protein_coordinate_model as vc  # noqa: E402

TP53_DANIO_RUN = ROOT / "runs" / "2026-07-21_1436_custom_run"

DOMAINS = [
    {"id": "D1", "label": "p53_DNA-bd", "start": 69, "end": 258},
    {"id": "D2", "label": "p53_tetramer", "start": 294, "end": 334},
]


# --------------------------------------------------------------------------- #
# Part 5 — mutually-exclusive boundary classification priority
# --------------------------------------------------------------------------- #
def test_exact_domain_edge():
    r = bc.classify_boundary(69, DOMAINS, threshold=5)
    assert r["class"] == bc.EXACT
    assert r["absolute_distance"] == 0
    assert r["signed_distance"] == 0
    assert r["nearest_edge_type"] == "start"


def test_near_domain_edge_threshold():
    assert bc.classify_boundary(72, DOMAINS, threshold=5)["class"] == bc.NEAR   # 3 aa
    assert bc.classify_boundary(74, DOMAINS, threshold=5)["class"] == bc.NEAR   # exactly 5 aa
    # 6 aa is beyond threshold and inside the domain -> inside, not near
    assert bc.classify_boundary(75, DOMAINS, threshold=5)["class"] == bc.INSIDE


def test_inside_domain():
    r = bc.classify_boundary(150, DOMAINS, threshold=5)
    assert r["class"] == bc.INSIDE
    assert r["absolute_distance"] > 5


def test_outside_annotated_domains():
    r = bc.classify_boundary(18, DOMAINS, threshold=5)
    assert r["class"] == bc.OUTSIDE
    assert r["signed_distance"] == 18 - 69  # negative, N-terminal of first domain


def test_unavailable_when_no_domain():
    r = bc.classify_boundary(100, [], threshold=5)
    assert r["class"] == bc.UNAVAILABLE
    assert r["nearest_domain_id"] is None


def test_signed_distance_sign_convention():
    left = bc.classify_boundary(60, DOMAINS)   # left of start 69
    right = bc.classify_boundary(270, DOMAINS)  # right of end 258
    assert left["signed_distance"] < 0
    assert right["signed_distance"] > 0


def test_legacy_class_mapping_is_total():
    for legacy in ("exact_edge", "near_edge", "inside_domain", "outside_domain", "unknown", ""):
        assert bc.canonical_class(legacy) in bc.CANONICAL_CLASSES


# --------------------------------------------------------------------------- #
# Part 1 — validator catches inconsistent coordinates loudly
# --------------------------------------------------------------------------- #
def _minimal_model(**over):
    m = {
        "schema_version": 1,
        "species_id": "danio_rerio",
        # A model states its own identity and role; nothing downstream may infer
        # which protein it is from list position or file name.
        "model_id": "gene:danio_rerio:primary",
        "model_role": "primary_reference",
        "is_primary_reference": True,
        "protein_id": "NP_x",
        "protein_length": 100,
        "coordinate_system": "protein_1_based_inclusive",
        "status": "available",
        "exons": [
            {"id": "e1", "label": "E1", "start": 1, "end": 50, "source": "gff",
             "source_file": "results/core_gene_analysis/exon_protein_map.tsv", "status": "ok"},
            {"id": "e2", "label": "E2", "start": 50, "end": 100, "source": "gff",
             "source_file": "results/core_gene_analysis/exon_protein_map.tsv", "status": "ok"},
        ],
        "exon_boundaries": [
            {"id": "b1", "label": "b1", "start": 50, "end": 50, "source": "core",
             "source_file": "x.tsv", "status": "inside_domain"},
        ],
        "representative_domains": [], "families_superfamilies": [], "member_signatures": [],
        "functional_sites": [], "disorder_regions": [], "tm_regions": [], "candidate_regions": [],
    }
    m.update(over)
    return m


def test_validator_accepts_consistent_model():
    assert validate_ok(_minimal_model())


def test_validator_rejects_out_of_range():
    m = _minimal_model(exons=[{"id": "e", "label": "E", "start": 1, "end": 200,
                               "source": "gff", "source_file": "x.tsv", "status": "ok"}])
    assert not validate_ok(m)


def test_validator_rejects_start_after_end():
    m = _minimal_model(representative_domains=[{"id": "d", "label": "D", "start": 60, "end": 40,
                                                "source": "InterProScan", "source_file": "x.tsv",
                                                "status": "representative_domain"}])
    assert not validate_ok(m)


def test_validator_rejects_boundary_off_exon_edge():
    m = _minimal_model(exon_boundaries=[{"id": "b", "label": "b", "start": 37, "end": 37,
                                         "source": "core", "source_file": "x.tsv", "status": "x"}])
    assert not validate_ok(m)


def test_validator_rejects_two_primaries_for_one_species():
    idx = {"schema_version": 1, "models": [_minimal_model(), _minimal_model()]}
    assert vc.validate_index(idx)  # non-empty -> violations


def validate_ok(model) -> bool:
    idx = {"schema_version": 1, "models": [model]}
    return not vc.validate_index(idx)


# --------------------------------------------------------------------------- #
# Part 16 — TP53 Danio post-cluster real data
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not TP53_DANIO_RUN.is_dir(), reason="TP53 Danio run not present")
def test_tp53_danio_model_matches_reference():
    idx = pcm.build_models_for_run(TP53_DANIO_RUN)
    assert idx["schema_version"] == 1
    assert idx["n_models"] == 1
    m = idx["models"][0]
    assert m["species_id"] == "danio_rerio"
    assert m["protein_id"] == "NP_001258749.1"
    assert m["transcript_id"] == "NM_001271820.1"
    assert m["protein_length"] == 374
    assert len(m["exons"]) == 10
    assert len(m["representative_domains"]) == 2
    assert len(m["tm_regions"]) == 0
    assert m["tm_analysis"]["performed"] is True
    assert m["tm_analysis"]["tm_region_count"] == 0
    assert len([b for b in m["exon_boundaries"]]) == 9
    # every boundary carries a canonical class
    for b in m["exon_boundaries"]:
        assert b["class"] in bc.CANONICAL_CLASSES


@pytest.mark.skipif(not TP53_DANIO_RUN.is_dir(), reason="TP53 Danio run not present")
def test_tp53_danio_model_validates():
    idx = pcm.build_models_for_run(TP53_DANIO_RUN)
    core = TP53_DANIO_RUN / "results" / "core_gene_analysis"
    assert vc.validate_index(idx, core_dir=core) == []


@pytest.mark.skipif(not TP53_DANIO_RUN.is_dir(), reason="TP53 Danio run not present")
def test_no_absolute_personal_paths_in_model():
    idx = pcm.build_models_for_run(TP53_DANIO_RUN)
    blob = str(idx)
    assert "/Users/" not in blob
    assert str(ROOT) not in blob


# --------------------------------------------------------------------------- #
# Part 11 / Part 16 — signed-distance + boundary-on-architecture plots
# --------------------------------------------------------------------------- #
def _load_plots():
    import os
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcache")
    sys.path.insert(0, str(ROOT / "scripts" / "plotting"))
    from exondomaincompare.presentation import shared_gene_plots as sp  # noqa
    sp.apply_style()
    return sp


def test_signed_distance_plot_writes_all_formats(tmp_path):
    sp = _load_plots()
    bounds = [
        {"exon_boundary_id": "b1", "boundary_position_aa": 69, "signed_distance_aa": 0,
         "class": "exact_domain_edge", "nearest_edge": "start"},
        {"exon_boundary_id": "b2", "boundary_position_aa": 18, "signed_distance_aa": -51,
         "class": "outside_annotated_domains", "nearest_edge": "start"},
    ]
    assert sp.plot_signed_boundary_distances(tmp_path, "sd", gene_symbol="TP53",
                                             species_name="Danio rerio", protein_id="NP_x",
                                             boundaries=bounds, threshold=5)
    for ext in ("svg", "pdf", "png"):
        assert (tmp_path / f"sd.{ext}").is_file()


def test_boundary_on_architecture_plot_writes_all_formats(tmp_path):
    sp = _load_plots()
    domains = [{"label": "p53_DNA-bd", "start": 69, "end": 258}]
    bounds = [{"boundary_position_aa": 18, "class": "outside_annotated_domains",
               "nearest_edge_position": 69}]
    assert sp.plot_boundary_on_architecture(tmp_path, "arch", gene_symbol="TP53",
                                            species_name="Danio rerio", protein_id="NP_x",
                                            protein_length=374, domains=domains,
                                            exon_boundaries=bounds,
                                            exon_blocks=[{"start": 1, "end": 18}], threshold=5)
    for ext in ("svg", "pdf", "png"):
        assert (tmp_path / f"arch.{ext}").is_file()


# --------------------------------------------------------------------------- #
# Part 2 — Exon Map static publication figures + model transcript-model layer
# --------------------------------------------------------------------------- #
_TMODELS = [
    {"protein_id": "P1", "transcript_id": "T1", "protein_length": 200, "is_primary": True,
     "curation_status": "curated", "exon_count": 3,
     "blocks": [{"id": "P1:e1", "label": "E1", "start": 1, "end": 60, "exon_number": 1,
                 "shared_exon_group_id": "g1"},
                {"id": "P1:e2", "label": "E2", "start": 61, "end": 130, "exon_number": 2,
                 "shared_exon_group_id": "g2"},
                {"id": "P1:e3", "label": "E3", "start": 131, "end": 200, "exon_number": 3,
                 "shared_exon_group_id": "g3"}]},
    {"protein_id": "P2", "transcript_id": "T2", "protein_length": 170, "is_primary": False,
     "curation_status": "predicted", "exon_count": 2,
     "blocks": [{"id": "P2:e1", "label": "E1", "start": 1, "end": 60, "exon_number": 1,
                 "shared_exon_group_id": "g1"},
                {"id": "P2:e3", "label": "E2", "start": 61, "end": 170, "exon_number": 2,
                 "shared_exon_group_id": "g3"}]},
]
_CAND = {"id": "C1", "start": 61, "end": 130, "candidate_type": "deletion"}


def test_transcript_model_comparison_all_writes_all_formats(tmp_path):
    sp = _load_plots()
    assert sp.plot_transcript_model_comparison(tmp_path, "cmp_all", gene_symbol="TP53",
                                               models=_TMODELS, candidate=_CAND, diff_only=False)
    for ext in ("svg", "pdf", "png"):
        assert (tmp_path / f"cmp_all.{ext}").is_file()


def test_transcript_model_comparison_diff_writes_all_formats(tmp_path):
    sp = _load_plots()
    assert sp.plot_transcript_model_comparison(tmp_path, "cmp_diff", gene_symbol="TP53",
                                               models=_TMODELS, candidate=_CAND, diff_only=True)
    for ext in ("svg", "pdf", "png"):
        assert (tmp_path / f"cmp_diff.{ext}").is_file()


def test_selected_candidate_detail_writes_all_formats(tmp_path):
    sp = _load_plots()
    blocks = [{"protein_start_aa": 1, "protein_end_aa": 60, "exon_number": 1},
              {"protein_start_aa": 61, "protein_end_aa": 130, "exon_number": 2}]
    assert sp.plot_selected_candidate_detail(tmp_path, "cand", gene_symbol="TP53",
                                             primary_id="P1", exon_blocks=blocks, candidate=_CAND)
    for ext in ("svg", "pdf", "png"):
        assert (tmp_path / f"cand.{ext}").is_file()


def test_domain_feature_stack_writes_all_formats(tmp_path):
    sp = _load_plots()
    exons = [{"start": 1, "end": 68, "label": "E1"}, {"start": 69, "end": 258, "label": "E2"}]
    cands = [{"id": "C1", "start": 61, "end": 130}]
    assert sp.plot_domain_feature_stack(
        tmp_path, "stack", gene_symbol="TP53", species_name="Danio rerio", protein_id="P1",
        protein_length=374, domains=DOMAINS, families=[{"start": 55, "end": 347, "label": "p53_fam"}],
        tms=[], exons=exons, boundaries=[69, 259], candidates=cands,
        lanes=("domains", "families", "tms", "exons", "boundaries", "candidates"))
    for ext in ("svg", "pdf", "png"):
        assert (tmp_path / f"stack.{ext}").is_file()


def test_member_signature_supplement_writes_all_formats(tmp_path):
    sp = _load_plots()
    sigs = [{"start": 69, "end": 258, "source": "PFAM", "interpro_accession": "IPR011615"},
            {"start": 55, "end": 347, "source": "PRINTS", "interpro_accession": "IPR002117"}]
    assert sp.plot_member_signature_supplement(
        tmp_path, "sig", gene_symbol="TP53", species_name="Danio rerio", protein_id="P1",
        protein_length=374, signatures=sigs)
    for ext in ("svg", "pdf", "png"):
        assert (tmp_path / f"sig.{ext}").is_file()


@pytest.mark.skipif(not TP53_DANIO_RUN.exists(), reason="TP53 Danio run not present")
def test_tp53_feature_layers_are_classified_from_real_interpro():
    model = pcm.build_models_for_run(TP53_DANIO_RUN)["models"][0]
    assert model["status"] == "available"
    # representative structural domains stay on their own layer
    labels = {d["label"] for d in model["representative_domains"]}
    assert any("DNA-bd" in x for x in labels)
    assert any("tetramer" in x.lower() for x in labels)
    # families/superfamilies, member signatures, sites and disorder are all real
    assert model["families_superfamilies"], "families/superfamilies must be populated"
    assert model["member_signatures"], "member signatures must be populated"
    assert model["disorder_regions"], "MOBIDB_LITE disorder must be populated"
    # real pyTMHMM result: zero TM, reported explicitly (never a blank absence)
    assert model["tm_analysis"]["performed"] is True
    assert model["tm_analysis"]["tm_region_count"] == 0
    assert "No transmembrane region predicted by pyTMHMM" in model["tm_analysis"]["message"]
    assert model["pending_info"] is None


FGFR1_RUN = ROOT / "runs" / "2026-07-23_1100_fgfr1_gallus_core_pilot"

# Required Part-1 boundary contract fields on every boundary object.
_BOUNDARY_FIELDS = (
    "id", "boundary_id", "label", "protein_position", "left_exon_id", "left_exon_label",
    "right_exon_id", "right_exon_label", "nearest_domain_id", "nearest_domain_label",
    "nearest_domain_start", "nearest_domain_end", "nearest_edge_type", "nearest_edge_position",
    "signed_distance", "absolute_distance", "boundary_class", "near_threshold", "mapping_status",
    "source", "source_file",
)


# --------------------------------------------------------------------------- #
# Boundary milestone — Part 1 data contract + Part 8 TP53 real result
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not TP53_DANIO_RUN.is_dir(), reason="TP53 Danio run not present")
def test_tp53_boundary_contract_and_class_summary():
    m = pcm.build_models_for_run(TP53_DANIO_RUN)["models"][0]
    bounds = m["exon_boundaries"]
    assert len(bounds) == 9
    length = m["protein_length"]
    for b in bounds:
        for f in _BOUNDARY_FIELDS:
            assert f in b, f"boundary missing field {f}"
        # every position is inside the protein coordinate system (Part 1 validation)
        assert 1 <= b["protein_position"] <= length
        assert b["near_threshold"] == 5
        assert b["mapping_status"] == "mapped"
        assert b["boundary_class"] in bc.CANONICAL_CLASSES
    summary = {c: 0 for c in bc.CANONICAL_CLASSES}
    for b in bounds:
        summary[b["boundary_class"]] += 1
    # current validated model: 5 inside_domain, 4 outside_annotated_domains
    assert summary[bc.INSIDE] == 5
    assert summary[bc.OUTSIDE] == 4
    assert summary[bc.EXACT] == 0 and summary[bc.NEAR] == 0 and summary[bc.UNAVAILABLE] == 0


@pytest.mark.skipif(not TP53_DANIO_RUN.is_dir(), reason="TP53 Danio run not present")
def test_tp53_every_boundary_matches_shared_classifier():
    """Check each boundary against the normalized feature coordinates via the
    shared classifier (Part 8: verify every boundary manually)."""
    m = pcm.build_models_for_run(TP53_DANIO_RUN)["models"][0]
    domains = [{"id": d["id"], "interpro_accession": d.get("interpro_accession"),
                "label": d["label"], "start": d["start"], "end": d["end"]}
               for d in m["representative_domains"]]
    for b in m["exon_boundaries"]:
        ref = bc.classify_boundary(b["protein_position"], domains, threshold=5)
        assert b["signed_distance"] == ref["signed_distance"]
        assert b["absolute_distance"] == ref["absolute_distance"]
        assert b["boundary_class"] == ref["class"]
        assert b["nearest_edge_type"] == ref["nearest_edge_type"]


@pytest.mark.skipif(not TP53_DANIO_RUN.is_dir(), reason="TP53 Danio run not present")
def test_tp53_boundary_left_right_exon_and_edge_position():
    m = pcm.build_models_for_run(TP53_DANIO_RUN)["models"][0]
    by_pos = {b["protein_position"]: b for b in m["exon_boundaries"]}
    # E1→E2 junction at aa 18 is outside the DNA-binding domain, N-terminal (signed<0)
    e1e2 = by_pos[18]
    assert e1e2["left_exon_label"] == "E1" and e1e2["right_exon_label"] == "E2"
    assert e1e2["signed_distance"] < 0
    assert e1e2["nearest_edge_type"] == "start"
    assert e1e2["nearest_edge_position"] == e1e2["nearest_domain_start"]
    # E3→E4 junction at aa 94 lies inside the DNA-binding domain
    e3e4 = by_pos[94]
    assert e3e4["boundary_class"] == bc.INSIDE
    assert e3e4["left_exon_label"] == "E3" and e3e4["right_exon_label"] == "E4"


@pytest.mark.skipif(not TP53_DANIO_RUN.is_dir(), reason="TP53 Danio run not present")
def test_tp53_boundary_selection_and_filtered_export_contract():
    """The fields that drive linked selection (nearest domain / adjacent exons)
    and the filtered-TSV export are all present and consistent."""
    m = pcm.build_models_for_run(TP53_DANIO_RUN)["models"][0]
    dom_ids = {d["id"] for d in m["representative_domains"]}
    exon_ids = {e["id"] for e in m["exons"]}
    for b in m["exon_boundaries"]:
        assert b["nearest_domain_id"] in dom_ids              # domain highlight target exists
        assert b["left_exon_id"] in exon_ids                  # adjacent-exon highlight targets exist
        assert b["right_exon_id"] in exon_ids
        # the exporter reads these — they must all be serialisable scalars
        assert isinstance(b["signed_distance"], int)
        assert isinstance(b["absolute_distance"], int)


# --------------------------------------------------------------------------- #
# Boundary milestone — Part 9 FGFR1 honest pending state
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not FGFR1_RUN.exists(), reason="FGFR1 Gallus pilot run not present")
def test_fgfr1_boundary_results_ready_positions():
    # Post-cluster discovery repair (Part A): the FGFR1 round-trip is complete, so
    # the coordinate model now carries real, classified boundaries — not the old
    # pre-cluster pending placeholders.
    m = pcm.build_models_for_run(FGFR1_RUN)["models"][0]
    bounds = m["exon_boundaries"]
    assert m["status"] == "available"
    assert len(bounds) == 16                                   # 17 coding exons → 16 internal
    length = m["protein_length"]
    for b in bounds:
        assert 1 <= b["protein_position"] <= length            # honest positions on the axis
        assert b["left_exon_label"] and b["right_exon_label"]  # still a real exon junction
    # boundaries are now classified against the real normalized domains
    assert any(b["boundary_class"] is not None for b in bounds)
    assert any(b["signed_distance"] is not None for b in bounds)


# --------------------------------------------------------------------------- #
# Boundary milestone — Part 10 static figures + Part 13 no-personal-paths
# --------------------------------------------------------------------------- #
def test_boundary_class_summary_plot_writes_all_formats(tmp_path):
    sp = _load_plots()
    bounds = [{"class": "inside_domain"}, {"class": "outside_annotated_domains"},
              {"class": "inside_domain"}]
    assert sp.plot_boundary_class_summary(tmp_path, "sum", gene_symbol="TP53",
                                          species_name="Danio rerio", protein_id="NP_x",
                                          boundaries=bounds, threshold=5)
    for ext in ("svg", "pdf", "png"):
        assert (tmp_path / f"sum.{ext}").is_file()


def test_selected_boundary_detail_plot_writes_all_formats(tmp_path):
    sp = _load_plots()
    b = {"label": "E1 → E2", "protein_position": 18, "signed_distance": -51,
         "absolute_distance": 51, "class": "outside_annotated_domains", "nearest_edge_position": 69}
    assert sp.plot_selected_boundary_detail(tmp_path, "det", gene_symbol="TP53",
                                            species_name="Danio rerio", protein_id="NP_x",
                                            boundary=b, domains=DOMAINS,
                                            exon_blocks=[{"start": 1, "end": 18, "label": "E1"},
                                                         {"start": 18, "end": 35, "label": "E2"}])
    for ext in ("svg", "pdf", "png"):
        assert (tmp_path / f"det.{ext}").is_file()


def test_boundary_evidence_supplement_plot_writes_all_formats(tmp_path):
    sp = _load_plots()
    bounds = [{"label": "E1 → E2", "protein_position": 18, "signed_distance": -51,
               "class": "outside_annotated_domains", "nearest_domain_label": "p53_DNA-bd",
               "nearest_edge_type": "start"}]
    assert sp.plot_boundary_evidence_supplement(tmp_path, "supp", gene_symbol="TP53",
                                                species_name="Danio rerio", protein_id="NP_x",
                                                boundaries=bounds, threshold=5)
    for ext in ("svg", "pdf", "png"):
        assert (tmp_path / f"supp.{ext}").is_file()


@pytest.mark.skipif(not TP53_DANIO_RUN.is_dir(), reason="TP53 Danio run not present")
def test_no_personal_paths_in_boundary_objects():
    m = pcm.build_models_for_run(TP53_DANIO_RUN)["models"][0]
    blob = str(m["exon_boundaries"])
    assert "/Users/" not in blob and str(ROOT) not in blob
    for b in m["exon_boundaries"]:
        assert b["source_file"].startswith("results/") or "core_gene_analysis" in b["source_file"]


def test_fgfr2_canonical_vocabulary_is_stable():
    """FGFR2 regression guard: the generic canonical vocabulary the Boundary tab
    uses must stay exactly these five mutually-exclusive classes (the frozen FGFR2
    Boundary Consistency vocabulary lives elsewhere and is never imported here)."""
    assert bc.CANONICAL_CLASSES == (
        "exact_domain_edge", "near_domain_edge", "inside_domain",
        "outside_annotated_domains", "unavailable_or_uncertain")


@pytest.mark.skipif(not (ROOT / "runs" / "2026-07-23_1100_fgfr1_gallus_core_pilot").exists(),
                    reason="FGFR1 Gallus pilot run not present")
def test_fgfr1_model_carries_transcript_models_and_candidate_c1():
    index = pcm.build_models_for_run(ROOT / "runs" / "2026-07-23_1100_fgfr1_gallus_core_pilot")
    model = index["models"][0]
    assert model["protein_id"] == "NP_990841.2"
    assert model["protein_length"] == 817
    assert model["n_transcript_models"] == 8
    assert any(t["is_primary"] for t in model["transcript_models"])
    c1 = next(c for c in model["candidate_regions"] if c["id"] == "C1")
    assert (c1["start"], c1["end"]) == (31, 118)
    assert c1["candidate_type"] == "exon_aligned_insertion"
    assert c1["affected_proteins"]  # real affected isoforms, not fabricated
    # Post-cluster discovery repair (Part A): the round-trip is complete, so the
    # domain/TM layers are now real results — no longer pending, never fabricated.
    assert model["status"] == "available"
    assert model["representative_domains"], "repaired run carries real domains"
    assert model["tm_analysis"].get("pending") is not True
