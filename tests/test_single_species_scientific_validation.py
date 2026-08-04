from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pdf_probe import probe_pdf, probe_png, probe_svg  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FGFR1_RUN = ROOT / "runs" / "2026-07-23_1100_fgfr1_gallus_core_pilot"
TP53_RUN = ROOT / "runs" / "2026-07-21_1436_custom_run"
RENDERER = ROOT / "scripts" / "plotting" / "render_main_figures.mjs"
FIGURE_DIR = Path("results") / "generic_gene_analysis" / "figures" / "main"

SHARED_MAIN_FIGURES = [
    "primary_exon_projection",
    "integrated_domain_architecture",
    "boundary_on_architecture",
    "signed_boundary_distances",
    "boundary_class_summary",
]

# Cards this phase removed as scientifically redundant or misleading. Their content
# is either present in an integrated main figure or is an on-demand export.
REMOVED_FIGURE_IDS = [
    "generic_domain_architecture",
    "generic_exon_domain_boundary_distribution",
    "primary_protein_exon_projection",
    "domain_arch_gallus_gallus_representative_architecture",
    "domain_arch_gallus_gallus_domain_exon_projection",
    "domain_arch_gallus_gallus_domain_boundary_overlay",
    "domain_arch_gallus_gallus_domain_candidate_overlay",
    "boundary_gallus_gallus_boundary_on_architecture",
    "boundary_gallus_gallus_signed_boundary_distances",
    "boundary_gallus_gallus_boundary_class_summary",
    "boundary_gallus_gallus_selected_boundary_detail",
    "boundary_gallus_gallus_boundary_evidence_supplement",
]


def _model(run_dir: Path) -> Path:
    return run_dir / "website_indices" / "generic" / "protein_coordinate_model.json"


def _gallery(run_dir: Path) -> dict:
    return json.loads((run_dir / "website_indices" / "figures_index.json").read_text())


@pytest.fixture(scope="module")
def fgfr1() -> Path:
    if not FGFR1_RUN.exists():
        pytest.skip(f"reference run missing: {FGFR1_RUN}")
    return FGFR1_RUN


@pytest.fixture(scope="module")
def gallery(fgfr1) -> dict:
    return _gallery(fgfr1)


@pytest.fixture(scope="module")
def cards(gallery) -> list[dict]:
    return gallery.get("figures") or []


# --------------------------------------------------------------------------- #
# The Gene Explorer figure and the Gallery figure are one figure
# --------------------------------------------------------------------------- #

def test_shipped_gallery_figures_are_reproduced_byte_for_byte(fgfr1, tmp_path):
    if shutil.which("node") is None:
        pytest.skip("node is required to re-render the figures")
    result = subprocess.run(
        ["node", str(RENDERER), str(_model(fgfr1)), str(tmp_path)],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert result.returncode == 0, result.stderr

    shipped_dir = fgfr1 / FIGURE_DIR
    compared = 0
    for svg in sorted(tmp_path.glob("main_*.svg")):
        shipped = shipped_dir / svg.name
        assert shipped.exists(), f"{svg.name} is not shipped in the run directory"
        assert hashlib.sha256(shipped.read_bytes()).hexdigest() == \
            hashlib.sha256(svg.read_bytes()).hexdigest(), \
            f"{svg.name} differs from a fresh render of the coordinate model"
        compared += 1
    assert compared == len(SHARED_MAIN_FIGURES)


def test_every_shared_main_figure_ships_all_formats_and_a_source_table(fgfr1):
    fig_dir = fgfr1 / FIGURE_DIR
    for kind in SHARED_MAIN_FIGURES:
        stem = f"main_gallus_gallus_{kind}"
        for ext in ("svg", "pdf", "png", "tsv"):
            assert (fig_dir / f"{stem}.{ext}").exists(), f"missing {stem}.{ext}"


def test_every_main_card_links_its_source_table(cards):
    for card in cards:
        if card["figure_id"].startswith("main_"):
            assert card.get("table_url"), f"{card['figure_id']} offers no source table"


def test_shipped_pdfs_are_vector_with_publication_page_sizes(fgfr1):
    for pdf in sorted((fgfr1 / FIGURE_DIR).glob("*.pdf")):
        info = probe_pdf(pdf)
        assert info.ok, f"{pdf.name}: {info.error}"
        assert info.images == [], f"{pdf.name} embeds raster data"
        assert info.has_vector_text, f"{pdf.name} has no selectable text"
        assert 2.0 <= info.width_in <= 14.0, f"{pdf.name} is {info.width_in:.1f}in wide"
        assert 1.0 <= info.height_in <= 14.0, f"{pdf.name} is {info.height_in:.1f}in tall"


def test_shipped_svgs_are_standalone(fgfr1):
    for svg in sorted((fgfr1 / FIGURE_DIR).glob("*.svg")):
        probe = probe_svg(svg)
        assert not probe["has_css_var"], f"{svg.name} contains var(--…)"
        assert not probe["has_class_attr"], f"{svg.name} relies on CSS classes"
        assert not probe["has_foreign_object"]
        assert probe["unpainted"] == [], f"{svg.name} has unpainted marks"


def test_shipped_pngs_are_three_hundred_dpi(fgfr1):
    for png in sorted((fgfr1 / FIGURE_DIR).glob("*.png")):
        probe = probe_png(png)
        assert probe["dpi"] == (300, 300), f"{png.name} reports dpi {probe['dpi']}"
        assert not probe["is_blank"]
        assert probe["frac_light"] > 0.5, f"{png.name} is mostly dark"


# --------------------------------------------------------------------------- #
# Gallery information architecture
# --------------------------------------------------------------------------- #

def test_gallery_has_one_card_per_figure_and_no_format_duplicates(cards):
    ids = [c["figure_id"] for c in cards]
    assert len(ids) == len(set(ids)), "duplicate figure_id in the Gallery index"
    for card in cards:
        # A card offers formats; a format must never become its own card.
        assert not re.search(r"\b(svg|pdf|png)\b", card["figure_id"].lower()), \
            f"format-specific card: {card['figure_id']}"
        assert not re.search(r"\.(svg|pdf|png|tsv)$", card.get("title", "").lower()), \
            f"filename used as a title: {card.get('title')!r}"


def test_removed_figures_are_absent_from_the_gallery(cards):
    present = {c["figure_id"] for c in cards}
    still_there = [f for f in REMOVED_FIGURE_IDS if f in present]
    assert not still_there, f"removed figures reappeared: {still_there}"


def test_the_shared_main_figures_are_registered_as_gallery_cards(cards):
    present = {c["figure_id"] for c in cards}
    for kind in SHARED_MAIN_FIGURES:
        assert f"main_gallus_gallus_{kind}" in present, f"main figure {kind} is not a card"


def test_main_cards_are_filed_under_the_expected_categories(cards):
    by_id = {c["figure_id"]: c for c in cards}
    expected = {
        "primary_exon_projection": "Exon structure",
        "integrated_domain_architecture": "Domain architecture",
        "boundary_on_architecture": "Exon–domain boundaries",
        "signed_boundary_distances": "Exon–domain boundaries",
        "boundary_class_summary": "Exon–domain boundaries",
    }
    for kind, category in expected.items():
        card = by_id[f"main_gallus_gallus_{kind}"]
        assert card.get("category") == category, \
            f"{kind} is filed under {card.get('category')!r}"


def test_every_main_card_is_scientifically_self_describing(cards):
    for card in cards:
        if not card["figure_id"].startswith("main_"):
            continue
        fid = card["figure_id"]
        assert card.get("title"), f"{fid} has no title"
        question = card.get("scientific_question") or ""
        assert question.endswith("?"), f"{fid} states no scientific question"
        assert len(card.get("interpretation") or "") > 40, \
            f"{fid} has no cautious interpretation"
        assert card.get("gene_symbol"), f"{fid} names no gene"
        assert card.get("species"), f"{fid} names no species"
        assert card.get("protein_id"), f"{fid} names no protein"
        assert card.get("stage") in ("pre_cluster", "post_cluster"), \
            f"{fid} declares no analysis stage"
        for key in ("svg_url", "pdf_url", "png_url"):
            assert card.get(key), f"{fid} offers no {key}"


def test_every_main_card_carries_a_caption_derived_from_its_own_data(cards, primary):
    by_id = {c["figure_id"]: c for c in cards}
    length = str(primary["protein_length"])
    protein = primary["protein_id"]
    n_boundaries = str(len(primary["exon_boundaries"]))

    for kind in SHARED_MAIN_FIGURES:
        caption = by_id[f"main_gallus_gallus_{kind}"].get("caption") or ""
        assert len(caption) > 60, f"{kind} has no usable caption"
        assert "FGFR1" in caption and "Gallus gallus" in caption
        assert protein in caption, f"{kind} caption names no protein"

    assert length in by_id["main_gallus_gallus_primary_exon_projection"]["caption"]
    assert "17 coding exons" in \
        by_id["main_gallus_gallus_primary_exon_projection"]["caption"]
    for kind in ("boundary_on_architecture", "signed_boundary_distances",
                 "boundary_class_summary"):
        caption = by_id[f"main_gallus_gallus_{kind}"]["caption"]
        assert n_boundaries in caption, f"{kind} caption states no boundary count"
    summary = by_id["main_gallus_gallus_boundary_class_summary"]["caption"]
    for expected in ("8 inside a domain", "6 near a domain edge",
                     "2 outside annotated domains"):
        assert expected in summary, f"class summary caption lacks: {expected}"


def test_candidate_related_cards_use_cautious_wording(cards):
    for card in cards:
        blob = " ".join(str(card.get(k, "")) for k in
                        ("title", "scientific_question", "interpretation", "caption"))
        if "candidate" not in blob.lower():
            continue
        low = blob.lower()
        assert "not validated" in low or "exploratory" in low, \
            f"{card['figure_id']} presents candidates without a caution"


def test_no_card_leaks_implementation_detail(cards):
    forbidden = ("schema_version", "json", "dataframe", "matplotlib", "html2canvas",
                 "jspdf", "figureSpec", "todo", "fixme")
    for card in cards:
        blob = " ".join(str(card.get(k, "")) for k in
                        ("title", "scientific_question", "interpretation")).lower()
        for word in forbidden:
            assert word.lower() not in blob, \
                f"{card['figure_id']} exposes implementation wording: {word}"


def test_no_absolute_personal_paths_leak_into_the_indices(fgfr1):
    for name in ("figures_index.json", "generic/figures_index.json",
                 "generic/protein_coordinate_model.json"):
        fp = fgfr1 / "website_indices" / name
        if fp.exists():
            assert "/Users/" not in fp.read_text(), f"{name} contains an absolute path"


# --------------------------------------------------------------------------- #
# Domain-instance identity and boundary bookkeeping
# --------------------------------------------------------------------------- #

def test_every_figure_names_a_candidate_the_same_way(fgfr1, primary):
    canonical = {(int(c["start"]), int(c["end"])): c["id"]
                 for c in primary.get("candidate_regions") or []}
    assert canonical, "the coordinate model lists no candidate regions"

    fig_dir = fgfr1 / FIGURE_DIR
    labelled = re.compile(r"(C\d+)\s*·\s*aa\s*(\d+)[–-](\d+)")
    checked = 0
    for svg in sorted(fig_dir.glob("main_*.svg")):
        text = svg.read_text()
        for label, start, end in labelled.findall(text):
            expected = canonical.get((int(start), int(end)))
            assert expected, \
                f"{svg.name} draws candidate aa {start}–{end}, which the " \
                f"coordinate model does not list"
            assert label == expected, (
                f"{svg.name} calls aa {start}–{end} {label!r}, but the coordinate "
                f"model calls it {expected!r}")
            checked += 1
    assert checked, "no figure labels a candidate, so nothing was verified"


@pytest.fixture(scope="module")
def primary(fgfr1) -> dict:
    doc = json.loads(_model(fgfr1).read_text())
    return next((m for m in doc["models"] if m.get("role") == "primary"), doc["models"][0])


def test_every_domain_instance_has_a_unique_stable_id(primary):
    domains = primary["representative_domains"]
    ids = [d["domain_instance_id"] for d in domains]
    assert len(ids) == len(set(ids))
    for d in domains:
        assert d["domain_instance_id"] == \
            f"{d['interpro_accession']}:{d['start']}-{d['end']}"


def test_repeated_accession_instances_are_numbered_in_coordinate_order(primary):
    ig = [d for d in primary["representative_domains"]
          if d["interpro_accession"] == "IPR007110"]
    assert len(ig) == 3
    assert [d["instance_number"] for d in ig] == [1, 2, 3]
    assert [d["start"] for d in ig] == sorted(d["start"] for d in ig)
    # Each instance is separately addressable and separately labelled.
    assert len({d["short_label"] for d in ig}) == 3


def test_boundary_records_store_the_instance_actually_used(primary):
    by_id = {d["domain_instance_id"]: d for d in primary["representative_domains"]}
    classified = [b for b in primary["exon_boundaries"]
                  if b.get("signed_distance") is not None]
    assert classified
    for b in classified:
        inst = by_id.get(b["nearest_domain_instance_id"])
        assert inst is not None, \
            f"{b['label']} references unknown instance {b.get('nearest_domain_instance_id')}"
        # The stored coordinates must be the coordinates of that instance, not of
        # some other instance sharing the accession.
        assert b["nearest_domain_start"] == inst["start"]
        assert b["nearest_domain_end"] == inst["end"]


def test_signed_distance_matches_the_stored_edge_coordinate(primary):
    for b in primary["exon_boundaries"]:
        if b.get("signed_distance") is None:
            continue
        pos = b["boundary_position_aa"]
        edge = b["nearest_edge_type"] if "nearest_edge_type" in b else b.get("nearest_edge")
        reference = b["nearest_domain_start"] if edge == "start" else b["nearest_domain_end"]
        assert b["signed_distance"] == pos - reference, \
            f"{b['label']}: signed distance {b['signed_distance']} does not match " \
            f"position {pos} minus {edge} edge {reference}"
        assert b["absolute_distance"] == abs(b["signed_distance"])


def test_inside_domain_boundaries_really_lie_inside_the_stored_domain(primary):
    inside = [b for b in primary["exon_boundaries"]
              if (b.get("boundary_class") or b.get("classification")) == "inside_domain"]
    assert inside
    for b in inside:
        assert b["nearest_domain_start"] <= b["boundary_position_aa"] <= b["nearest_domain_end"]


def test_near_edge_boundaries_respect_the_declared_threshold(primary):
    threshold = primary["near_edge_threshold_aa"]
    near = [b for b in primary["exon_boundaries"]
            if (b.get("boundary_class") or b.get("classification")) == "near_domain_edge"]
    assert near
    for b in near:
        assert abs(b["signed_distance"]) <= threshold


def test_boundary_labels_use_exon_transition_form(primary):
    for b in primary["exon_boundaries"]:
        assert re.fullmatch(r"E\d+ → E\d+", b["label"]), \
            f"unexpected boundary label: {b['label']!r}"


# --------------------------------------------------------------------------- #
# Regression: TP53 and the immutable FGFR2 reference
# --------------------------------------------------------------------------- #

def test_tp53_gallery_carries_the_same_main_figure_set():
    if not TP53_RUN.exists():
        pytest.skip("TP53 regression run missing")
    ids = {c["figure_id"] for c in (_gallery(TP53_RUN).get("figures") or [])}
    for kind in SHARED_MAIN_FIGURES:
        assert any(i.endswith(kind) and i.startswith("main_") for i in ids), \
            f"TP53 Gallery lacks a {kind} card"


def test_tp53_domain_instances_and_boundaries_are_internally_consistent():
    if not TP53_RUN.exists():
        pytest.skip("TP53 regression run missing")
    doc = json.loads(_model(TP53_RUN).read_text())
    for model in doc["models"]:
        by_id = {d["domain_instance_id"]: d for d in model.get("representative_domains") or []}
        assert len(by_id) == len(model.get("representative_domains") or [])
        for b in model.get("exon_boundaries") or []:
            if b.get("signed_distance") is None:
                continue
            inst = by_id.get(b["nearest_domain_instance_id"])
            assert inst is not None
            assert b["nearest_domain_start"] == inst["start"]
            assert b["nearest_domain_end"] == inst["end"]


def test_the_fgfr2_figure_index_is_untouched_by_this_phase():
    result = subprocess.run(
        ["git", "status", "--porcelain", "results/final_30_until_interpro_prepare"],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert result.returncode == 0
    assert result.stdout.strip() == "", \
        f"the FGFR2 freeze has local modifications:\n{result.stdout}"


def test_no_run_figure_directory_contains_a_full_page_raster_pdf():
    offenders = []
    for run in (FGFR1_RUN, TP53_RUN):
        if not run.exists():
            continue
        for pdf in (run / FIGURE_DIR).glob("*.pdf"):
            info = probe_pdf(pdf)
            if info.is_single_raster_page or info.width_pt > 2000:
                offenders.append(f"{run.name}/{pdf.name} "
                                 f"({info.width_pt:.0f}pt, images={info.images})")
    assert not offenders, f"raster or oversized PDFs found: {offenders}"
