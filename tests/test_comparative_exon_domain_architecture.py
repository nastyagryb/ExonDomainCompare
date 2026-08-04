from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from shutil import which

import pytest

ROOT = Path(__file__).resolve().parents[1]
FRONTEND = ROOT / "webapp" / "frontend" / "src"
BUILDER = FRONTEND / "pages" / "viewers" / "comparativeGalleryFigures.js"
GENERATOR = ROOT / "scripts" / "plotting" / "generate_comparative_gallery_figures.py"
RENDERER = ROOT / "scripts" / "plotting" / "render_comparative_gallery_figures.mjs"

CARD_IDS = ("cmp_exon_domain_architecture_native", "cmp_exon_domain_architecture_msa")


def test_one_builder_serves_both_coordinate_systems():
    text = BUILDER.read_text(encoding="utf-8")
    assert text.count("export function comparativeExonDomainArchitectureFigureSpec") == 1
    body = text.split("export function comparativeExonDomainArchitectureFigureSpec", 1)[1]
    body = body.split("\nexport function ", 1)[0]
    assert 'const useMsa = mode === "msa"' in body
    # Nothing branches on a gene symbol: the gene is data, not a code path.
    for gene in ("FGFR2", "FGFR1", "TP53", "TPM1"):
        assert gene not in body, f"the builder special-cases {gene}"


def test_the_two_tracks_of_a_species_share_one_axis():
    body = BUILDER.read_text(encoding="utf-8").split(
        "export function comparativeExonDomainArchitectureFigureSpec", 1)[1]
    body = body.split("\nexport function ", 1)[0]
    # One scale for the whole figure, and a connector spanning both tracks.
    assert body.count("const scale = (v) =>") == 1
    assert "fig.line(x, domainY - 9, x, exonY + 6," in body


def test_the_boundary_class_comes_from_the_boundary_analysis():
    body = BUILDER.read_text(encoding="utf-8").split(
        "export function comparativeExonDomainArchitectureFigureSpec", 1)[1]
    body = body.split("\nexport function ", 1)[0]
    assert "canonClass(b.boundary_class || b.class)" in body
    assert "boundaryClassColour(cls)" in body
    assert "boundaryClassLabel(cls)" in body


def test_species_are_grouped_in_the_canonical_order():
    body = BUILDER.read_text(encoding="utf-8").split(
        "export function comparativeExonDomainArchitectureFigureSpec", 1)[1]
    body = body.split("\nexport function ", 1)[0]
    assert "speciesOrder(models)" in body


@pytest.mark.skipif(which("node") is None, reason="node is required")
def test_the_rendered_geometry_holds():
    proc = subprocess.run([which("node"), str(ROOT / "tests"
                                              / "check_exon_domain_architecture.mjs")],
                          capture_output=True, text=True, cwd=ROOT)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "FAIL" not in proc.stdout, proc.stdout


def test_both_cards_are_registered_with_a_question_and_an_interpretation():
    text = GENERATOR.read_text(encoding="utf-8")
    for card in CARD_IDS:
        assert f'"{card}"' in text, f"{card} is not in the comparative inventory"
    meta = text.split("FIGURE_META", 1)[1]
    for card in CARD_IDS:
        block = meta.split(f'"{card}": {{', 1)[1].split("},", 1)[0]
        for field in ("title", "category", "kind", "question", "interpretation"):
            assert f'"{field}"' in block, f"{card} has no {field}"
        assert "Comparative exon–domain architecture" in block


def test_both_cards_declare_their_source_table():
    text = GENERATOR.read_text(encoding="utf-8")
    sources = text.split("SOURCE_ARTEFACT", 1)[1]
    for card in CARD_IDS:
        assert f'"{card}": "msa_aligned_domains"' in sources


def test_the_renderer_emits_both_cards():
    text = RENDERER.read_text(encoding="utf-8")
    assert "comparativeExonDomainArchitectureFigureSpec" in text
    assert 'emit("cmp_exon_domain_architecture_native"' in text
    assert 'emit("cmp_exon_domain_architecture_msa"' in text
    # The MSA panel is only emitted when the alignment actually has columns.
    msa_block = text.split('emit("cmp_exon_domain_architecture_msa"', 1)[0]
    assert msa_block.rstrip().endswith("if (nColumns) {")


def _multi_species_runs() -> list[Path]:
    out = []
    for card in sorted((ROOT / "runs").glob("*/website_indices/figures_index.json")):
        try:
            doc = json.loads(card.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        ids = {f.get("figure_id") for f in (doc.get("figures") or [])}
        if "cmp_domain_architecture_native" in ids:
            out.append(card.parents[1])
    return out


def test_every_multi_species_dataset_carries_both_cards():
    runs = _multi_species_runs()
    if not runs:
        pytest.skip("no multi-species dataset with a comparative gallery")
    for run_dir in runs:
        doc = json.loads((run_dir / "website_indices"
                          / "figures_index.json").read_text(encoding="utf-8"))
        ids = {f.get("figure_id") for f in (doc.get("figures") or [])}
        missing = [c for c in CARD_IDS if c not in ids]
        assert not missing, f"{run_dir.name} is missing {missing}"


def test_the_card_files_exist_in_every_published_format():
    runs = _multi_species_runs()
    if not runs:
        pytest.skip("no multi-species dataset with a comparative gallery")
    for run_dir in runs:
        doc = json.loads((run_dir / "website_indices"
                          / "figures_index.json").read_text(encoding="utf-8"))
        by_id = {f.get("figure_id"): f for f in (doc.get("figures") or [])}
        for card in CARD_IDS:
            entry = by_id.get(card)
            assert entry, f"{run_dir.name}: {card} missing"
            for url_key in ("png_url", "svg_url", "pdf_url", "table_url"):
                assert entry.get(url_key), f"{run_dir.name}: {card} has no {url_key}"


def test_a_single_species_dataset_gets_no_comparative_card():
    generator = GENERATOR.read_text(encoding="utf-8")
    assert 'return {"figures": 0, "cards": 0, "skipped": "single_species"}' in generator
    renderer = RENDERER.read_text(encoding="utf-8")
    assert "comparative gallery requires at least two species" in renderer
