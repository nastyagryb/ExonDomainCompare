#!/usr/bin/env python3
"""The Comparative Exon–Domain Boundary Explorer, checked against the real two-species run.

Two layers are covered here:

* the backend contract in ``src/exondomaincompare/shared_gene_analysis/boundary_dashboard.py``, which
  is the single place that decides which boundaries are comparable, and
* the browser-side behaviour (filtering, matrix cells, pair connections, exports), which
  is exercised by ``tests/check_comparative_explorer.mjs`` in Node and reported here.

The division matters: the frontend must not contain a second comparable-boundary
algorithm, so these tests also assert that it does not grow one.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

MODEL_REL = Path("website_indices") / "generic" / "protein_coordinate_model.json"
FGFR1_MULTI = ROOT / "runs" / "2026-07-26_2157_fgfr1_gallus_mus_core_pilot"
FGFR1_SINGLE = ROOT / "runs" / "2026-07-23_1100_fgfr1_gallus_core_pilot"
VIEWERS = ROOT / "webapp" / "frontend" / "src" / "pages" / "viewers"
HARNESS = Path(__file__).with_name("check_comparative_explorer.mjs")

pytestmark = pytest.mark.skipif(not (FGFR1_MULTI / MODEL_REL).is_file(),
                                reason="two-species coordinate model not present")


@pytest.fixture(scope="module")
def multi():
    index = json.loads((FGFR1_MULTI / MODEL_REL).read_text(encoding="utf-8"))
    return index["boundary_dashboard"]["multi_species"]


@pytest.fixture(scope="module")
def index():
    return json.loads((FGFR1_MULTI / MODEL_REL).read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- #
# the browser-side behaviour
# --------------------------------------------------------------------------- #
def test_the_explorer_behaviour_holds_on_the_real_dataset():
    """Filters, matrix cells, pair connections and exports, exercised in Node."""
    if shutil.which("node") is None:
        pytest.skip("node is required to exercise the explorer's JavaScript")
    proc = subprocess.run(["node", str(HARNESS), str(FGFR1_MULTI / MODEL_REL)],
                          capture_output=True, text=True, cwd=str(ROOT))
    assert proc.returncode == 0, f"{proc.stdout}\n{proc.stderr}"
    assert "comparative explorer behaviour verified" in proc.stdout
    # A harness that silently checked nothing would also exit zero.
    assert proc.stdout.count("  ok    ") >= 45, proc.stdout


# --------------------------------------------------------------------------- #
# the comparative index is the single source of the comparability decision
# --------------------------------------------------------------------------- #
def test_the_frontend_holds_no_second_comparability_algorithm():
    """The browser may filter and draw; it may not decide what is comparable."""
    sources = [p for p in VIEWERS.glob("*.js*") if p.name.startswith("comparative")
               or p.name.startswith("Comparative")]
    assert sources, "the comparative explorer sources are missing"
    joined = "\n".join(p.read_text(encoding="utf-8") for p in sources)
    # Matching boundaries across species means bucketing them by an evidence key. If any
    # of that appears in the browser, the backend has stopped being the single authority.
    for forbidden in ("match_comparable", "matchComparable", "buildComparableGroups",
                      "NEAR_COLUMN_TOLERANCE"):
        assert forbidden not in joined, (
            f"{forbidden} in the frontend duplicates the backend's comparability rule")
    for source in sources:
        text = source.read_text(encoding="utf-8")
        assert "exon_number" not in text and "exon_rank" not in text, (
            f"{source.name} appears to compare boundaries by exon rank")


def test_comparable_groups_never_rest_on_exon_rank(multi):
    methods = {g["mapping_method"] for g in multi["comparable_boundary_groups"]}
    assert methods, "no comparable groups to check"
    assert methods <= {"msa_aligned_position", "shared_exon_group"}, methods
    for g in multi["comparable_boundary_groups"]:
        evidence = g.get("supporting_evidence") or {}
        assert evidence, f"{g['comparable_boundary_group_id']} carries no evidence"
        assert "exon_number" not in evidence and "exon_rank" not in evidence


def test_every_group_spans_at_least_two_species(multi):
    for g in multi["comparable_boundary_groups"]:
        species = {o["species_id"] for o in g["per_species_native_positions"]}
        assert len(species) >= 2, (
            f"{g['comparable_boundary_group_id']} is a single-species group and must not "
            f"be published as comparable")


# --------------------------------------------------------------------------- #
# the observation contract the matrix hover and the detail rows depend on
# --------------------------------------------------------------------------- #
REQUIRED_OBSERVATION_FIELDS = (
    "species_id", "scientific_name", "taxonomic_group", "protein_id", "boundary_id",
    "exon_transition", "native_position", "msa_column", "nearest_domain_instance_id",
    "nearest_domain_label", "nearest_domain_start", "nearest_domain_end",
    "nearest_edge", "nearest_edge_position", "signed_distance", "absolute_distance",
    "boundary_class", "domain_annotation_available", "mapping_method",
    "mapping_status", "mapping_confidence",
)


def test_each_observation_carries_every_field_the_ui_shows(multi):
    for g in multi["comparable_boundary_groups"]:
        for o in g["per_species_native_positions"]:
            missing = [f for f in REQUIRED_OBSERVATION_FIELDS if f not in o]
            assert not missing, (
                f"{g['comparable_boundary_group_id']}/{o['species_id']} is missing "
                f"{missing}; the hover would have to invent them")


def test_the_taxonomic_group_is_real(multi):
    """The taxonomic-group filter needs a value, and it must come from the registry."""
    groups = multi["filter_options"]["taxonomic_groups"]
    assert groups, "no taxonomic groups published for the filter"
    assert "Analysed species" not in groups, (
        "a species fell back to the placeholder group, so the registry lookup failed")


def test_matrix_cells_mirror_the_group_detail(multi):
    groups = {g["comparable_boundary_group_id"]: g
              for g in multi["comparable_boundary_groups"]}
    for row in multi["boundary_matrix"]:
        for cell in row["cells"]:
            g = groups[cell["comparable_boundary_group_id"]]
            obs = next((o for o in g["per_species_native_positions"]
                        if o["species_id"] == row["species_id"]), None)
            if obs is None:
                assert not cell["observed"]
                assert cell["signed_distance"] is None
            else:
                assert cell["observed"]
                assert cell["signed_distance"] == obs["signed_distance"]
                assert cell["observation"] == obs, (
                    "the cell and the detail row must be the same record")


def test_a_cell_without_an_observation_carries_no_distance(index):
    """The absent case, forced.

    In the real two-species run every comparable group is mapped in both species, so no
    unobserved cell exists to check. That is exactly why this test constructs one: a
    missing observation rendered as ``0`` would read as "sits exactly on a domain edge",
    which is a real boundary class and the opposite of missing data. The failure mode is
    invisible on a fully mapped dataset and would surface only on a run with a gap.
    """
    from exondomaincompare.shared_gene_analysis import boundary_dashboard as bd

    models = index["models"]
    groups = bd.match_comparable_boundaries(models)
    assert groups
    # Drop one species' observation from the first group.
    orphaned = json.loads(json.dumps(groups))
    removed = orphaned[0]["per_species_native_positions"].pop()
    matrix = bd.build_boundary_matrix(models, orphaned)

    row = next(r for r in matrix if r["species_id"] == removed["species_id"])
    cell = row["cells"][0]
    assert cell["observed"] is False
    assert cell["state"] == "boundary_absent_or_unmapped"
    for field in ("native_position", "signed_distance", "absolute_distance",
                  "observation"):
        assert cell[field] is None, (
            f"an unobserved cell must not carry {field}={cell[field]!r}")


# --------------------------------------------------------------------------- #
# consistency statistics stay honest at n = 2
# --------------------------------------------------------------------------- #
def test_two_species_statistics_lead_with_the_raw_pair(multi):
    for s in multi["distance_statistics"]:
        raw = s["raw_signed_distances"]
        if len(raw) != 2:
            continue
        assert s["primary_statistic"] == "raw_pair", (
            "with two observations the median is just their mean; the raw pair must be "
            "the primary statistic")
        assert s["cross_species_difference"] == abs(
            raw[0]["signed_distance"] - raw[1]["signed_distance"])


def test_the_statistics_avoid_conservation_claims(multi):
    for s in multi["distance_statistics"]:
        assert s["metric_label"] == "Boundary-position consistency"
        assert "conserv" not in json.dumps(s).lower(), (
            "a consistent boundary position is not evidence of conservation")


def test_annotation_gaps_are_reported_not_hidden(multi):
    for s in multi["distance_statistics"]:
        expected = all(r["domain_annotation_available"]
                       for r in s["raw_signed_distances"])
        assert s["domain_annotation_available_in_all"] == expected


# --------------------------------------------------------------------------- #
# inspection cases
# --------------------------------------------------------------------------- #
def test_inspection_cases_cover_the_real_discrepancies(multi):
    cases = multi["inspection_cases"]
    assert cases, "the run has differing classes and a tentative mapping to report"
    types = {c["case_type"] for c in cases}
    assert "different_boundary_classes" in types
    assert "tentative_mapping" in types
    for c in cases:
        assert c["comparable_boundary_group_id"]
        assert c["detail"] and c["label"]
        assert c["severity"] in {"review", "caution"}


def test_a_tentative_case_explains_why_it_is_tentative(multi):
    tentative = [c for c in multi["inspection_cases"]
                 if c["case_type"] == "tentative_mapping"]
    assert tentative
    for c in tentative:
        assert "not established" in c["detail"] or "close but not identical" in c["detail"]


def test_inspection_cases_do_not_call_biology_a_technical_error(multi):
    banned = ("error", "bug", "broken", "wrong", "invalid")
    for c in multi["inspection_cases"]:
        text = f"{c['label']} {c['detail']}".lower()
        for word in banned:
            assert word not in text, f"{c['case_id']} calls a discrepancy a {word}"


# --------------------------------------------------------------------------- #
# filter vocabularies are published, not hardcoded in the UI
# --------------------------------------------------------------------------- #
def test_the_filter_vocabularies_come_from_the_data(multi):
    opts = multi["filter_options"]
    for key in ("species", "taxonomic_groups", "boundary_classes",
                "representative_domain_groups", "mapping_statuses", "edges",
                "inspection_case_types"):
        assert key in opts, f"filter option {key} is not published"
    observed_classes = {o["boundary_class"] for g in multi["comparable_boundary_groups"]
                        for o in g["per_species_native_positions"]}
    assert set(opts["boundary_classes"]) == observed_classes, (
        "the class filter must offer exactly the classes that occur")
    assert set(opts["edges"]) <= {"start", "end"}


# --------------------------------------------------------------------------- #
# rendered figures
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def rendered(tmp_path_factory):
    if shutil.which("node") is None:
        pytest.skip("node is required to render the comparative figures")
    out = tmp_path_factory.mktemp("comparative_figures")
    proc = subprocess.run(
        ["node", str(ROOT / "scripts/plotting/render_comparative_figures.mjs"),
         str(FGFR1_MULTI / MODEL_REL), str(out)],
        capture_output=True, text=True, cwd=str(ROOT))
    assert proc.returncode == 0, f"{proc.stdout}\n{proc.stderr}"
    return out


@pytest.mark.parametrize("stem", [
    "comparative_boundary_matrix",
    "comparative_paired_signed_distance",
    "comparative_boundary_consistency",
])
def test_each_comparative_figure_is_a_real_vector_file(rendered, stem):
    svg = rendered / f"{stem}.svg"
    pdf = rendered / f"{stem}.pdf"
    assert svg.is_file() and pdf.is_file()
    text = svg.read_text(encoding="utf-8")
    assert "<svg" in text
    # A screenshot-based export would contain a raster payload instead of shapes.
    assert "data:image/png" not in text and "<image" not in text
    assert 'class="' not in text, "a CSS-dependent SVG renders black outside the browser"
    head = pdf.read_bytes()[:8]
    assert head.startswith(b"%PDF-")
    assert pdf.stat().st_size > 1000


def test_the_renderer_reports_no_layout_warnings(rendered):
    summary = json.loads((rendered / "comparative_render_summary.json").read_text())
    assert summary
    for fig in summary:
        assert not fig["warnings"], f"{fig['stem']}: {fig['warnings']}"
        assert fig["marks"] > 20, f"{fig['stem']} drew only {fig['marks']} marks"


def test_the_comparative_figures_ship_their_source_tables(rendered):
    for stem in ("comparative_boundary_matrix", "comparative_paired_signed_distance",
                 "comparative_boundary_consistency"):
        tsv = rendered / f"{stem}.tsv"
        assert tsv.is_file(), f"{stem} has no source table"
        lines = tsv.read_text(encoding="utf-8").strip().split("\n")
        assert len(lines) > 1, f"{stem}.tsv has no rows"


def test_a_single_species_run_produces_no_comparative_figures(tmp_path):
    """Refusing to render is correct: an empty comparative panel is a claim."""
    if shutil.which("node") is None:
        pytest.skip("node is required")
    if not (FGFR1_SINGLE / MODEL_REL).is_file():
        pytest.skip("single-species reference run not present")
    proc = subprocess.run(
        ["node", str(ROOT / "scripts/plotting/render_comparative_figures.mjs"),
         str(FGFR1_SINGLE / MODEL_REL), str(tmp_path)],
        capture_output=True, text=True, cwd=str(ROOT))
    assert proc.returncode != 0
    assert "nothing to render" in proc.stderr


# --------------------------------------------------------------------------- #
# the validated FGFR2 explorer and the single-species views are untouched
# --------------------------------------------------------------------------- #
def test_the_validated_fgfr2_boundary_page_is_not_rerouted():
    page = (ROOT / "webapp/frontend/src/pages/BoundaryPage.jsx").read_text(encoding="utf-8")
    assert 'eventType !== "validated"' in page, (
        "the generic dashboard must stay behind the non-validated branch so the "
        "validated FGFR2 Boundary Consistency Explorer keeps its own page")
    assert "ValidatedOrPendingBoundaryPage" in page


def test_the_comparative_explorer_does_not_touch_the_frozen_vocabulary():
    for name in ("ComparativeBoundaryDashboard.jsx", "ComparativeBoundaryMatrix.jsx",
                 "ComparativePairedPlot.jsx", "comparativeFigures.js",
                 "comparativeFilters.js"):
        text = (VIEWERS / name).read_text(encoding="utf-8")
        assert "fgfr2Styles" not in text and "boundary.js" not in text, (
            f"{name} reaches into the frozen FGFR2 vocabulary")


def test_the_single_species_dashboard_still_renders_its_own_mode():
    dash = (ROOT / "webapp/frontend/src/pages/GlobalBoundaryDashboard.jsx").read_text(
        encoding="utf-8")
    assert "generic_multi_species_results_ready" in dash
    assert "SingleSpeciesOverview" in dash, (
        "the single-species overview must survive alongside the comparative explorer")
