"""Artefact validation for the redesigned single-species Gallery figures.

These tests inspect the exported files — a standalone SVG, a vector PDF, a 300 dpi
PNG and the source table — rather than the code that draws them, because the
defects being guarded against were invisible in the source: labels that overlapped
or were clipped, a genomic axis annotated in amino acids, "10/10 neighbours
resolved" on a single-species locus map, and every downstream exon coloured as
altered because an upstream deletion had shifted the protein coordinates.

Every figure is regenerated into a throw-away copy of the real run, so the test
also proves the generators still run end to end for both reference datasets:
FGFR1 / Gallus gallus (post-cluster) and TP53 / Danio rerio.
"""

from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pdf_probe import probe_pdf, probe_png, probe_svg  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
VENV_PY = ROOT / "venv" / "bin" / "python"
FGFR1_RUN = ROOT / "runs" / "2026-07-23_1100_fgfr1_gallus_core_pilot"
TP53_RUN = ROOT / "runs" / "2026-07-21_1436_custom_run"

# The generators that own the figures under test, in pipeline order.
GENERATORS = (
    "scripts/plotting/generate_exon_map_figures.py",
    "scripts/plotting/generate_domain_figures.py",
    "scripts/plotting/generate_boundary_figures.py",
)

ANALYSIS = "results/generic_gene_analysis"
EXON_MAP_DIR = f"{ANALYSIS}/figures/exon_map"
DOMAIN_DIR = f"{ANALYSIS}/figures/domain_architecture"

# figure_id, output directory, file stem, Gallery category, card kind.
# One entry per figure — one Gallery card per figure, never one card per format.
OWNED_FIGURES = [
    ("transcript_exon_structure", EXON_MAP_DIR,
     "exon_map_{sp}_transcript_and_protein_structure", "Exon structure", "main"),
    ("transcript_model_comparison", EXON_MAP_DIR,
     "exon_map_{sp}_model_comparison_all", "Exon structure", "main"),
    ("transcript_model_comparison_differences", EXON_MAP_DIR,
     "exon_map_{sp}_model_comparison_differences", "Exon structure", "main"),
    ("local_gene_neighbourhood", EXON_MAP_DIR,
     "exon_map_{sp}_local_gene_neighbourhood", "Genomic context", "main"),
    ("exploratory_candidate_ranking", EXON_MAP_DIR,
     "exon_map_{sp}_exploratory_candidate_ranking", "Exploratory candidates", "main"),
    ("generic_candidate_domain_context", EXON_MAP_DIR,
     "exon_map_{sp}_candidate_domain_context", "Exploratory candidates", "main"),
    ("domain_arch_{sp}_member_signature_supplement", DOMAIN_DIR,
     "domain_arch_{sp}_member_signature_supplement", "Domain architecture",
     "supplement"),
]
FIGURE_IDS = [f[0] for f in OWNED_FIGURES]

# Cards these generators retire. None of them may come back, in either index.
RETIRED_FIGURE_IDS = [
    "exon_map_{sp}_primary_projection",
    "exon_map_{sp}_selected_candidate_detail",
    "generic_exploratory_event_candidates",
    "generic_synteny_neighbourhood",
    "domain_arch_{sp}_representative_architecture",
    "domain_arch_{sp}_domain_exon_projection",
    "domain_arch_{sp}_domain_boundary_overlay",
    "domain_arch_{sp}_domain_candidate_overlay",
    "generic_domain_architecture",
    "boundary_{sp}_boundary_on_architecture",
    "boundary_{sp}_signed_boundary_distances",
    "boundary_{sp}_boundary_class_summary",
    "boundary_{sp}_selected_boundary_detail",
    "boundary_{sp}_boundary_evidence_supplement",
    "generic_exon_domain_boundary_distribution",
]

# Publication page geometry: narrower than a journal column or wider than a
# landscape page is a layout defect, not a figure.
MIN_PAGE_IN = 2.0
MAX_PAGE_IN = 14.0


# --------------------------------------------------------------------------- #
# Regenerating the artefacts
# --------------------------------------------------------------------------- #

def _sandbox_run(source: Path, dest: Path) -> Path:
    """A minimal copy of a real run: the inputs the figure generators read.

    The generators write figures, rewrite the export manifest inside the served
    coordinate model and register Gallery cards, so they are pointed at a copy
    rather than at the reference run itself.

    The copy keeps the run's directory name. Card registration rejects a card whose
    output paths name a different run, so renaming the copy would make every card
    copied in from the reference index look like it came from somewhere else.
    """
    dest = dest / source.name
    (dest / ANALYSIS).mkdir(parents=True, exist_ok=True)
    for tsv in (source / ANALYSIS).glob("*.tsv"):
        shutil.copyfile(tsv, dest / ANALYSIS / tsv.name)
    for name in ("figures_index.json", "generic/figures_index.json",
                 "generic/protein_coordinate_model.json"):
        src = source / "website_indices" / name
        if src.is_file():
            out = dest / "website_indices" / name
            out.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, out)
    return dest


def _generate(run_dir: Path) -> None:
    env = dict(os.environ)
    env["PYTHONPYCACHEPREFIX"] = str(ROOT / "tmp" / "pycache")
    env["MPLCONFIGDIR"] = str(ROOT / "tmp" / "mpl")
    python = str(VENV_PY) if VENV_PY.exists() else sys.executable
    for script in GENERATORS:
        result = subprocess.run([python, script, str(run_dir)],
                                capture_output=True, text=True, cwd=ROOT, env=env)
        if result.returncode != 0:
            pytest.fail(f"{script} failed for {run_dir.name}:\n"
                        f"{result.stdout}\n{result.stderr}")


def _species_id(run_dir: Path) -> str:
    doc = json.loads((run_dir / "website_indices" / "generic"
                      / "protein_coordinate_model.json").read_text())
    models = doc["models"]
    primary = next((m for m in models if m.get("role") == "primary"), models[0])
    return primary.get("species_id") or "sp"


class Figures:
    """The regenerated artefacts of one run, addressed by figure id."""

    def __init__(self, run_dir: Path, cards_before: set | None = None):
        self.run_dir = run_dir
        self.cards_before = cards_before or set()
        self.species_id = _species_id(run_dir)
        self.model = json.loads((run_dir / "website_indices" / "generic"
                                 / "protein_coordinate_model.json").read_text())

    def stem(self, figure_id: str) -> Path:
        for fid, directory, stem, _cat, _kind in OWNED_FIGURES:
            if fid.format(sp=self.species_id) == figure_id.format(sp=self.species_id):
                return (self.run_dir / directory
                        / stem.format(sp=self.species_id))
        raise KeyError(figure_id)

    def path(self, figure_id: str, ext: str) -> Path:
        return self.stem(figure_id).with_suffix(f".{ext}")

    def card_id(self, figure_id: str) -> str:
        return figure_id.format(sp=self.species_id)

    def index(self, name: str = "figures_index.json") -> dict:
        return json.loads((self.run_dir / "website_indices" / name).read_text())

    def cards(self, name: str = "figures_index.json") -> list:
        return self.index(name).get("figures") or []

    def card(self, figure_id: str) -> dict:
        wanted = self.card_id(figure_id)
        found = [c for c in self.cards() if c.get("figure_id") == wanted]
        assert found, f"no Gallery card for {wanted}"
        return found[0]

    def table(self, figure_id: str) -> list:
        with self.path(figure_id, "tsv").open(encoding="utf-8") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))


def _card_ids(run_dir: Path) -> set:
    path = run_dir / "website_indices" / "figures_index.json"
    if not path.is_file():
        return set()
    doc = json.loads(path.read_text())
    return {f.get("figure_id") for f in (doc.get("figures") or [])}


@pytest.fixture(scope="module")
def fgfr1(tmp_path_factory) -> Figures:
    if not FGFR1_RUN.exists():
        pytest.skip(f"reference run missing: {FGFR1_RUN}")
    run = _sandbox_run(FGFR1_RUN, tmp_path_factory.mktemp("fgfr1_run"))
    before = _card_ids(run)
    _generate(run)
    return Figures(run, before)


@pytest.fixture(scope="module")
def tp53(tmp_path_factory) -> Figures:
    if not TP53_RUN.exists():
        pytest.skip(f"regression run missing: {TP53_RUN}")
    run = _sandbox_run(TP53_RUN, tmp_path_factory.mktemp("tp53_run"))
    before = _card_ids(run)
    _generate(run)
    return Figures(run, before)


# --------------------------------------------------------------------------- #
# Text helpers: matplotlib writes PDF text glyph by glyph
# --------------------------------------------------------------------------- #

def _squeeze(text: str) -> str:
    """Comparable text: matplotlib writes PDF strings as two-byte glyph codes."""
    return re.sub(r"[\s\x00]+", "", text).lower()


def pdf_text(path: Path) -> str:
    """The PDF text layer with whitespace removed, for substring assertions."""
    return _squeeze(probe_pdf(path).text)


def svg_text(path: Path) -> str:
    return _squeeze("\n".join(t for t in probe_svg(path)["texts"] if t))


def has(haystack: str, needle: str) -> bool:
    return _squeeze(needle) in haystack


# --------------------------------------------------------------------------- #
# The reference datasets really are the datasets we claim to plot
# --------------------------------------------------------------------------- #

def test_reference_run_is_the_real_fgfr1_gallus_dataset(fgfr1):
    primary = fgfr1.model["models"][0]
    assert fgfr1.model["gene_symbol"] == "FGFR1"
    assert primary["species_id"] == "gallus_gallus"
    assert primary["protein_id"] == "NP_990841.2"
    assert primary["transcript_id"] == "NM_205510.2"
    assert primary["protein_length"] == 817
    assert len(primary["exons"]) == 17
    assert len(primary["transcript_models"]) == 8
    assert [d["domain_instance_id"] for d in primary["representative_domains"]] == [
        "IPR007110:33-118", "IPR007110:145-244", "IPR007110:253-355",
        "IPR001245:476-750",
    ]


def test_regression_run_is_the_real_tp53_danio_dataset(tp53):
    primary = tp53.model["models"][0]
    assert tp53.model["gene_symbol"] == "TP53"
    assert primary["species_id"] == "danio_rerio"


# --------------------------------------------------------------------------- #
# Every figure ships every format, from one card
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_figure_exports_svg_pdf_png_and_a_source_table(fgfr1, figure_id):
    for ext in ("svg", "pdf", "png", "tsv"):
        path = fgfr1.path(figure_id, ext)
        assert path.is_file(), f"{path.name} was not exported"
        assert path.stat().st_size > 0, f"{path.name} is empty"
    assert len(fgfr1.table(figure_id)) >= 1, f"{figure_id} source table has no rows"


@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_figure_exports_every_format_for_the_tp53_run(tp53, figure_id):
    for ext in ("svg", "pdf", "png", "tsv"):
        path = tp53.path(figure_id, ext)
        assert path.is_file(), f"{path.name} was not exported for TP53"
        assert path.stat().st_size > 0, f"{path.name} is empty for TP53"


# --------------------------------------------------------------------------- #
# PDF: true vector at a sensible physical size
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_pdf_is_one_well_formed_page(fgfr1, figure_id):
    info = probe_pdf(fgfr1.path(figure_id, "pdf"))
    assert info.ok, info.error
    assert info.header.startswith("%PDF-")
    assert info.n_pages == 1


@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_pdf_embeds_no_raster_image(fgfr1, figure_id):
    info = probe_pdf(fgfr1.path(figure_id, "pdf"))
    assert info.images == [], f"{figure_id} PDF embeds raster data: {info.images}"
    assert not info.is_single_raster_page


@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_pdf_has_selectable_text_in_a_referenced_font(fgfr1, figure_id):
    info = probe_pdf(fgfr1.path(figure_id, "pdf"))
    assert info.fonts, f"{figure_id} PDF references no font"
    assert info.has_vector_text
    assert info.n_text_ops >= 10


@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_pdf_page_has_publication_dimensions(fgfr1, figure_id):
    info = probe_pdf(fgfr1.path(figure_id, "pdf"))
    assert MIN_PAGE_IN <= info.width_in <= MAX_PAGE_IN, \
        f"{figure_id} PDF is {info.width_in:.1f}in wide"
    assert MIN_PAGE_IN <= info.height_in <= MAX_PAGE_IN, \
        f"{figure_id} PDF is {info.height_in:.1f}in tall"


@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_pdf_draws_vector_geometry(fgfr1, figure_id):
    info = probe_pdf(fgfr1.path(figure_id, "pdf"))
    assert info.n_path_ops >= 10, f"{figure_id} PDF has almost no vector geometry"
    assert info.n_fill_ops >= 5


@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_tp53_pdf_is_vector_too(tp53, figure_id):
    info = probe_pdf(tp53.path(figure_id, "pdf"))
    assert info.ok, info.error
    assert info.images == []
    assert info.fonts and info.n_text_ops >= 10
    assert MIN_PAGE_IN <= info.width_in <= MAX_PAGE_IN
    assert MIN_PAGE_IN <= info.height_in <= MAX_PAGE_IN


# --------------------------------------------------------------------------- #
# SVG: valid, self-contained, explicitly white
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_svg_is_valid_xml_with_real_marks_and_text(fgfr1, figure_id):
    info = probe_svg(fgfr1.path(figure_id, "svg"))  # raises on invalid XML
    assert info["viewBox"], f"{figure_id} SVG has no viewBox"
    assert info["n_marks"] >= 10
    assert info["n_text"] >= 8


@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_svg_is_self_contained_and_white_backed(fgfr1, figure_id):
    info = probe_svg(fgfr1.path(figure_id, "svg"))
    assert not info["has_css_var"], f"{figure_id} SVG carries unresolved var(--…)"
    assert not info["has_external_css"]
    assert not info["has_foreign_object"]
    # The figure patch is painted pure white, so the SVG does not inherit a dark
    # page background when it is opened outside the app.
    assert "#ffffff" in info["raw"].lower(), \
        f"{figure_id} SVG has no explicit white background"


# --------------------------------------------------------------------------- #
# PNG: 300 dpi, light, readable
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_png_is_300_dpi_and_readable(fgfr1, figure_id):
    info = probe_png(fgfr1.path(figure_id, "png"))
    assert info["dpi"] == (300, 300), f"{figure_id} PNG dpi is {info['dpi']}"
    assert not info["is_blank"]
    assert info["frac_light"] > 0.6, \
        f"{figure_id} PNG is not a light figure ({info['frac_light']:.2f} light)"
    assert info["frac_dark"] < 0.25
    assert info["contrast"] > 0.5
    assert info["width"] >= 600 and info["height"] >= 300


# --------------------------------------------------------------------------- #
# A. two panels, two coordinate systems
# --------------------------------------------------------------------------- #

def test_transcript_figure_has_a_nucleotide_and_an_amino_acid_axis(fgfr1):
    text = pdf_text(fgfr1.path("transcript_exon_structure", "pdf"))
    assert has(text, "Genomic coordinate on the assembly"), \
        "panel A does not label a nucleotide axis"
    assert has(text, "Amino-acid position in each translated protein model"), \
        "panel B does not label an amino-acid axis"
    assert has(text, "panel A in nucleotides"), "the two panels are not distinguished"
    assert has(text, "amino acids")


def test_transcript_figure_never_labels_the_nucleotide_axis_in_amino_acids(fgfr1):
    """The original defect: one axis, labelled as both at once."""
    text = pdf_text(fgfr1.path("transcript_exon_structure", "pdf"))
    assert not has(text, "Amino-acid coordinate on the assembly")
    assert not has(text, "Genomic position (aa)")
    # The nucleotide axis is stated in kb, the protein axis in aa, and neither
    # unit appears on the other axis label.
    assert has(text, "(kb, NCBI RefSeq annotation)")


def test_transcript_figure_names_the_transcript_and_the_protein(fgfr1):
    text = pdf_text(fgfr1.path("transcript_exon_structure", "pdf"))
    assert has(text, "NM_205510.2")
    assert has(text, "NP_990841.2")
    assert has(text, "primary")
    assert has(text, "817 aa")
    assert has(text, "minus strand"), "the transcription direction is not stated"


def test_transcript_figure_suppresses_labels_it_cannot_place(fgfr1):
    """Narrow exons drop their label instead of overprinting the neighbour."""
    labels = [t.strip() for t
              in probe_svg(fgfr1.path("transcript_exon_structure", "svg"))["texts"]
              if re.fullmatch(r"E\d+", (t or "").strip())]
    assert labels, "no exon labels at all"
    counts = Counter(labels)
    # 17 coding exons in each panel, so a label placed in both panels appears
    # twice. On the nucleotide axis the narrowest exons are too small to label.
    assert set(counts) == {f"E{i}" for i in range(1, 18)}
    assert max(counts.values()) == 2
    suppressed = [label for label, n in counts.items() if n == 1]
    assert 1 <= len(suppressed) <= 6, \
        f"expected a few narrow exons to stay unlabelled, got {suppressed}"


# --------------------------------------------------------------------------- #
# B. exon identity is compared on genomic coordinates
# --------------------------------------------------------------------------- #

def _identity_counts(figures: Figures, figure_id: str) -> Counter:
    return Counter(r["exon_identity"] for r in figures.table(figure_id))


def test_fgfr1_exon_identity_counts_are_the_real_ones(fgfr1):
    """The real FGFR1 numbers, per model, from the exported source table.

    140 exon rows: 17 primary exons plus 123 rows for the 7 alternative models.
    Only 14 of those 123 differ from the primary — 4 alternative exons, 4 shifted
    boundaries and 6 absent exons. The pre-redesign figure marked whole tails of
    exons as altered because it compared protein offsets.
    """
    rows = fgfr1.table("transcript_model_comparison")
    assert len(rows) == 140
    per_model: dict[str, Counter] = {}
    for r in rows:
        per_model.setdefault(r["protein_id"], Counter())[r["exon_identity"]] += 1

    assert per_model["NP_990841.2"] == Counter({"primary": 17})
    assert per_model["XP_015152847.2"] == Counter(
        {"shared": 16, "alternative": 1, "missing": 1})
    assert per_model["XP_015152849.2"] == Counter(
        {"shared": 15, "alternative": 1, "shifted": 1, "missing": 1})
    assert per_model["XP_015152852.2"] == Counter(
        {"shared": 14, "alternative": 1, "shifted": 1, "missing": 2})
    assert per_model["XP_024998538.1"] == Counter({"shared": 16, "missing": 1})
    assert per_model["XP_040507404.1"] == Counter({"shared": 16, "shifted": 1})
    assert per_model["XP_040507405.1"] == Counter({"shared": 17})
    assert per_model["XP_046759308.1"] == Counter(
        {"shared": 15, "alternative": 1, "shifted": 1, "missing": 1})

    counts = _identity_counts(fgfr1, "transcript_model_comparison")
    assert counts["primary"] == 17
    assert counts["shared"] == 109
    assert counts["alternative"] == 4
    assert counts["shifted"] == 4
    assert counts["missing"] == 6


def test_transcript_comparison_does_not_flag_every_downstream_exon(fgfr1):
    """A protein-coordinate shift is not an exon difference."""
    counts = _identity_counts(fgfr1, "transcript_model_comparison")
    altered = counts["alternative"] + counts["shifted"] + counts["missing"]
    non_primary = sum(counts.values()) - counts["primary"]
    assert non_primary == 123
    assert altered == 14, "the real FGFR1 figure flags 14 of 123 exon rows"
    # The broken comparison marked well over half of the alternative-model exons.
    assert altered / non_primary < 0.2
    assert counts["shared"] > altered * 5


def test_transcript_comparison_states_lengths_and_a_shared_amino_acid_axis(fgfr1):
    text = pdf_text(fgfr1.path("transcript_model_comparison", "pdf"))
    assert has(text, "Amino-acid position in each protein model")
    assert has(text, "common scale")  # one axis for every model

    assert has(text, "817 aa")  # the primary
    assert has(text, "821 aa")  # a longer alternative model
    assert has(text, "NP_990841.2") and has(text, "primary")
    assert has(text, "reference model")


def test_transcript_comparison_legend_names_every_identity_class(fgfr1):
    text = pdf_text(fgfr1.path("transcript_model_comparison", "pdf"))
    for entry in ("shared genomic exon", "shifted exon boundary",
                  "missing protein region", "coding exon of the primary model"):
        assert has(text, entry), f"legend entry {entry!r} missing"


def test_differences_only_figure_drops_the_identical_models(fgfr1):
    """XP_040507405.1 shares all 17 exons, so it belongs only in the all-models view."""
    all_text = pdf_text(fgfr1.path("transcript_model_comparison", "pdf"))
    diff_text = pdf_text(fgfr1.path("transcript_model_comparison_differences", "pdf"))
    assert has(all_text, "XP_040507405.1")
    assert not has(diff_text, "XP_040507405.1")
    # the primary stays, as the reference the differences are measured against
    assert has(diff_text, "NP_990841.2")


# --------------------------------------------------------------------------- #
# C. member-database signature supplement
# --------------------------------------------------------------------------- #

def test_signature_supplement_is_a_supplement_grouped_by_member_database(fgfr1):
    card = fgfr1.card("domain_arch_{sp}_member_signature_supplement")
    assert card["kind"] != "main", "the signature supplement is not a main figure"
    assert card["category"] == "Domain architecture"
    text = pdf_text(fgfr1.path("domain_arch_{sp}_member_signature_supplement", "pdf"))
    databases = {r["member_database"]
                 for r in fgfr1.table("domain_arch_{sp}_member_signature_supplement")}
    assert len(databases) >= 3
    for db in databases:
        assert has(text, db), f"member database {db!r} has no track header"


def test_signature_supplement_separates_integrated_from_unintegrated(fgfr1):
    rows = fgfr1.table("domain_arch_{sp}_member_signature_supplement")
    statuses = {r["integration_status"] for r in rows}
    assert "integrated" in statuses
    text = pdf_text(fgfr1.path("domain_arch_{sp}_member_signature_supplement", "pdf"))
    assert has(text, "integrated into an InterPro entry")
    if "unintegrated" in statuses:
        assert has(text, "not integrated")


def test_signature_supplement_rows_are_grouped_and_ordered(fgfr1):
    rows = fgfr1.table("domain_arch_{sp}_member_signature_supplement")
    blocks: list[str] = []
    for r in rows:
        if not blocks or blocks[-1] != r["member_database"]:
            blocks.append(r["member_database"])
    # Every member database is one contiguous block of tracks, never interleaved.
    assert len(blocks) == len(set(blocks)), f"interleaved databases: {blocks}"
    for db in set(blocks):
        starts = [int(r["start_aa"]) for r in rows if r["member_database"] == db]
        assert starts == sorted(starts), f"{db} tracks are not in positional order"


# --------------------------------------------------------------------------- #
# D. local genomic neighbourhood
# --------------------------------------------------------------------------- #

def test_neighbourhood_figure_titles_and_counts_come_from_the_data(fgfr1):
    text = pdf_text(fgfr1.path("local_gene_neighbourhood", "pdf"))
    rows = fgfr1.table("local_gene_neighbourhood")
    up = sum(1 for r in rows if r["side"] == "upstream")
    down = sum(1 for r in rows if r["side"] == "downstream")
    assert up == 5 and down == 5, "the real FGFR1 locus has 5 flanking loci per side"
    card = fgfr1.card("local_gene_neighbourhood")
    assert card["title"] == "FGFR1 · Local genomic neighbourhood"
    assert card["category"] == "Genomic context"
    assert has(text, f"{up + down} flanking loci shown")
    assert has(text, f"{up} upstream") and has(text, f"{down} downstream")


def test_neighbourhood_figure_makes_no_orthology_claim(fgfr1):
    """A single-species locus map cannot carry orthology confidence."""
    pdf = pdf_text(fgfr1.path("local_gene_neighbourhood", "pdf"))
    svg = svg_text(fgfr1.path("local_gene_neighbourhood", "svg"))
    for forbidden in ("neighbours resolved", "neighbors resolved", "ortholog",
                      "confidence", "high confidence", "unresolved neighbour",
                      "10/10"):
        assert not has(pdf, forbidden), f"{forbidden!r} appears in the PDF"
        assert not has(svg, forbidden), f"{forbidden!r} appears in the SVG"


def test_neighbourhood_figure_names_the_placeholder_locus(fgfr1):
    text = pdf_text(fgfr1.path("local_gene_neighbourhood", "pdf"))
    rows = fgfr1.table("local_gene_neighbourhood")
    placeholders = [r for r in rows if r["locus_kind"] == "placeholder_locus"]
    assert placeholders, "the real FGFR1 neighbourhood contains LOC121107413"
    assert any(has(text, r["locus_symbol"]) for r in placeholders)
    assert has(text, "LOC121107413")
    assert has(text, "placeholder locus")


def test_neighbourhood_legend_names_the_four_things_it_draws(fgfr1):
    text = pdf_text(fgfr1.path("local_gene_neighbourhood", "pdf"))
    for entry in ("target gene", "annotated neighbouring gene",
                  "placeholder locus", "transcription direction"):
        assert has(text, entry), f"legend entry {entry!r} missing"


def test_neighbourhood_directions_match_the_annotation(fgfr1):
    rows = fgfr1.table("local_gene_neighbourhood")
    target = [r for r in rows if r["locus_kind"] == "target"]
    assert len(target) == 1
    # FGFR1 is annotated on the minus strand in this assembly.
    assert target[0]["transcription_direction"] == "-"
    assert {r["transcription_direction"] for r in rows} <= {"+", "-"}


# --------------------------------------------------------------------------- #
# E. exploratory candidates — a score is not a validation
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("figure_id", ["exploratory_candidate_ranking",
                                       "generic_candidate_domain_context"])
def test_candidate_figure_is_labelled_exploratory_and_not_validated(fgfr1, figure_id):
    text = pdf_text(fgfr1.path(figure_id, "pdf"))
    assert has(text, "not validated"), f"{figure_id} does not state its status"
    assert has(text, "exploratory candidate")
    card = fgfr1.card(figure_id)
    assert card["category"] == "Exploratory candidates"
    assert "not validated" in card["interpretation"].lower()


@pytest.mark.parametrize("figure_id", ["exploratory_candidate_ranking",
                                       "generic_candidate_domain_context"])
def test_candidate_figure_never_calls_a_score_a_validation(fgfr1, figure_id):
    text = pdf_text(fgfr1.path(figure_id, "pdf"))
    for forbidden in ("validated event", "confirmed", "validated candidate",
                      "validation score", "biological validation: high",
                      "experimentally"):
        assert not has(text, forbidden), f"{forbidden!r} appears in {figure_id}"


def test_candidate_ranking_shows_every_ranking_column(fgfr1):
    text = pdf_text(fgfr1.path("exploratory_candidate_ranking", "pdf"))
    for header in ("Rank", "Candidate", "Interval (aa)", "Length (aa)",
                   "Affected isoforms", "Supporting comparisons", "Evidence score",
                   "Evidence strength", "Biological validation"):
        assert has(text, header), f"column {header!r} missing"
    rows = fgfr1.table("exploratory_candidate_ranking")
    assert [int(r["rank"]) for r in rows] == list(range(1, len(rows) + 1))
    assert all(r["biological_validation"] == "not validated" for r in rows)
    top = rows[0]
    # C1 is the real top-ranked candidate: aa 31–118, 88 aa long.
    assert top["candidate_label"] == "C1"
    assert (int(top["aa_start"]), int(top["aa_end"])) == (31, 118)
    assert int(top["length_aa"]) == 88
    assert has(text, "31") and has(text, "118")


def test_candidate_context_resolves_the_real_domain_instance(fgfr1):
    """C1 overlaps Ig-like domain 1 (aa 33–118), not a collapsed "IPR007110"."""
    rows = fgfr1.table("generic_candidate_domain_context")
    c1 = next(r for r in rows if r["candidate_label"] == "C1")
    # The instance id carries the real coordinates of that one instance.
    assert c1["overlapping_domain_instance_ids"] == "IPR007110:33-118"
    assert c1["overlapping_domain_labels"] == "Ig-like domain 1"
    assert c1["nearest_domain_instance_id"] == "IPR007110:33-118"
    text = pdf_text(fgfr1.path("generic_candidate_domain_context", "pdf"))
    assert has(text, "Ig-like domain 1")
    # all three Ig-like instances stay distinct in the figure
    for label in ("Ig-like domain 1", "Ig-like domain 2", "Ig-like domain 3"):
        assert has(text, label)


def test_candidate_context_shows_exons_boundaries_and_the_nearest_domain_edge(fgfr1):
    text = pdf_text(fgfr1.path("generic_candidate_domain_context", "pdf"))
    assert has(text, "Coding exons and internal boundaries")
    assert has(text, "nearest edge")
    assert has(text, "Amino-acid position on NP_990841.2")
    rows = fgfr1.table("generic_candidate_domain_context")
    c1 = next(r for r in rows if r["candidate_label"] == "C1")
    assert c1["nearest_edge_type"] in {"start", "end"}
    assert int(c1["distance_to_nearest_edge_aa"]) >= 0
    assert c1["coincident_exon_boundary"], \
        "C1 coincides with a coding-exon boundary in the real data"


# --------------------------------------------------------------------------- #
# The Gallery contract
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("figure_id,_dir,_stem,category,kind", OWNED_FIGURES)
def test_each_figure_has_exactly_one_card_in_the_right_category(
        fgfr1, figure_id, _dir, _stem, category, kind):
    wanted = fgfr1.card_id(figure_id)
    matches = [c for c in fgfr1.cards() if c.get("figure_id") == wanted]
    assert len(matches) == 1, f"{wanted} has {len(matches)} cards, expected 1"
    card = matches[0]
    assert card["category"] == category
    assert card["kind"] == kind
    assert card["status"] == "available"


def test_no_card_is_a_per_format_duplicate(fgfr1):
    ids = [c["figure_id"] for c in fgfr1.cards()]
    assert len(ids) == len(set(ids)), "duplicate Gallery cards"
    for suffix in ("_svg", "_pdf", "_png", "_tsv"):
        assert not [i for i in ids if i.endswith(suffix)], \
            "a file format was registered as its own card"


@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_every_card_carries_a_question_and_an_interpretation(fgfr1, figure_id):
    card = fgfr1.card(figure_id)
    assert card["title"] and not card["title"].endswith((".svg", ".pdf", ".png"))
    assert len(card["scientific_question"]) > 20
    assert card["scientific_question"].endswith("?")
    assert len(card["interpretation"]) > 40
    assert card["caption"]
    assert card["svg_url"] and card["pdf_url"] and card["png_url"]
    assert card["table_url"]


@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_card_urls_point_at_files_that_exist(fgfr1, figure_id):
    card = fgfr1.card(figure_id)
    for key in ("svg_url", "pdf_url", "png_url", "table_url"):
        rel = card[key].split("path=", 1)[1].split("&", 1)[0]
        assert (fgfr1.run_dir / rel).is_file(), f"{key} points at a missing file"


def test_retired_cards_do_not_reappear(fgfr1):
    sp = fgfr1.species_id
    retired = {t.format(sp=sp) for t in RETIRED_FIGURE_IDS}
    for index_name in ("figures_index.json", "generic/figures_index.json"):
        try:
            cards = fgfr1.cards(index_name)
        except FileNotFoundError:
            continue
        present = {c.get("figure_id") for c in cards} & retired
        assert not present, f"retired cards back in {index_name}: {sorted(present)}"


def test_cards_owned_by_other_stages_survive(fgfr1):
    """The alignment and integrated main figures belong to other stages.

    A generator may replace its own cards and retire the ones it supersedes, but it
    must never delete a card it does not own — an earlier version of this stage
    wiped whole Gallery sections.
    """
    sp = fgfr1.species_id
    retired = {t.format(sp=sp) for t in RETIRED_FIGURE_IDS}
    mine = {fgfr1.card_id(f) for f in FIGURE_IDS}
    foreign = {i for i in fgfr1.cards_before if i} - retired - mine
    assert foreign, "the reference index carried no foreign cards to preserve"
    ids = {c.get("figure_id") for c in fgfr1.cards()}
    assert foreign <= ids, \
        f"cards owned by other stages were dropped: {sorted(foreign - ids)}"
    assert any("alignment" in i for i in ids), "the alignment figures were dropped"


def test_tp53_registers_the_same_cards(tp53):
    ids = {c.get("figure_id") for c in tp53.cards()}
    for figure_id in FIGURE_IDS:
        assert tp53.card_id(figure_id) in ids, \
            f"{figure_id} has no card in the TP53 run"


def test_tp53_neighbourhood_and_candidate_wording_hold_too(tp53):
    nb = pdf_text(tp53.path("local_gene_neighbourhood", "pdf"))
    assert not has(nb, "neighbours resolved")
    assert not has(nb, "ortholog")
    assert has(nb, "flanking loci shown")
    for figure_id in ("exploratory_candidate_ranking",
                      "generic_candidate_domain_context"):
        assert has(pdf_text(tp53.path(figure_id, "pdf")), "not validated")
