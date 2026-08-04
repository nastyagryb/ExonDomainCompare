from __future__ import annotations

import importlib
import json
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

RUNNER = ROOT / "src" / "exondomaincompare" / "framework" / "run_core_gene_analysis.py"
SEQUENCE = ROOT / "scripts" / "plotting" / "figure_sequence.py"
FGFR1_RUN = ROOT / "runs" / "2026-07-23_1100_fgfr1_gallus_core_pilot"
ANALYSIS = Path("results") / "generic_gene_analysis"

# Every figure stage the pipeline is expected to run, in invocation order. The
# shared main figures must come first: they retire the cards they supersede, so the
# supplements are registered into an already-cleaned set. The comparative stage
# skips cleanly on a single-species run. Registration is last and is the only owner
# of the availability record, so no producing stage can leave it describing a card
# set the figures have moved on from.
FIGURE_STAGES = (
    "plotting.generate_shared_main_figures",
    "plotting.generate_exon_map_figures",
    "plotting.generate_domain_figures",
    "plotting.generate_boundary_figures",
    "plotting.generate_alignment_figures",
    "plotting.generate_comparative_gallery_figures",
    "plotting.figure_registration",
)

# The final single-species set, keyed by the analysis each card answers.
EXPECTED_CARDS = {
    "Exon structure": 4,
    "Isoform analysis": 3,
    "Domain architecture": 2,
    "Exon–domain boundaries": 3,
    "Genomic context": 1,
    "Exploratory candidates": 2,
}


def test_the_pipeline_names_every_figure_stage():
    from plotting.figure_sequence import FIGURE_STAGES as declared
    assert tuple(module for _label, module in declared) == FIGURE_STAGES
    assert "run_figure_stages" in RUNNER.read_text(), \
        "the pipeline no longer runs the shared figure sequence"


def test_the_shared_main_figures_run_before_the_supplements():
    text = SEQUENCE.read_text()
    positions = [text.index(s) for s in FIGURE_STAGES if s in text]
    assert positions == sorted(positions), \
        "the figure stages are not invoked in the documented order"


def _phase_body(name: str) -> str:
    text = RUNNER.read_text()
    assert f"def {name}(" in text, f"{RUNNER.name} has no {name}()"
    return text.split(f"def {name}(", 1)[1].split("\ndef ", 1)[0]


@pytest.mark.parametrize("phase", ["phase_create", "phase_post"])
def test_the_figure_stages_run_after_the_generic_pipeline(phase):
    body = _phase_body(phase)
    generic = body.find("_run_generic_pipeline(")
    figures = body.find("_build_coordinate_model_and_figures(")
    assert generic != -1, f"{phase} no longer runs the generic pipeline"
    assert figures != -1, f"{phase} no longer builds the publication figures"
    assert figures > generic, (
        f"in {phase} the publication figure stages run before the generic pipeline, "
        f"which rewrites figures_index.json and would discard their cards")


def test_the_generic_orchestrator_keeps_cards_it_does_not_own(tmp_path):
    sys.path.insert(0, str(ROOT / "scripts"))
    from generic_gene.run_generic_gene_analysis import _keep_publication_cards

    index = tmp_path / "figures_index.json"
    index.write_text(json.dumps({"figures": [
        {"figure_id": "main_sp_primary_exon_projection", "category": "Exon structure",
         "title": "Projection", "status": "available", "svg_url": "u"},
        {"figure_id": "primary_protein_exon_projection", "title": "old pre-cluster card"},
    ]}))

    rebuilt = {
        "figures": [{"figure_id": "primary_protein_exon_projection",
                     "title": "old pre-cluster card", "status": "available"}],
        "available": [{"id": "primary_protein_exon_projection"}],
    }
    _keep_publication_cards(index, rebuilt)

    ids = {f["figure_id"] for f in rebuilt["figures"]}
    assert "main_sp_primary_exon_projection" in ids, \
        "a card registered by a figure stage was thrown away"
    assert "main_sp_primary_exon_projection" in {a["id"] for a in rebuilt["available"]}
    # The manifest-derived card is left for the figure stages to retire.
    assert "primary_protein_exon_projection" in ids


@pytest.fixture(scope="module")
def fresh_gallery(tmp_path_factory) -> list[dict]:
    if not FGFR1_RUN.exists():
        pytest.skip(f"reference run missing: {FGFR1_RUN}")
    if shutil.which("node") is None:
        pytest.skip("node is required by the shared figure renderers")

    dest = tmp_path_factory.mktemp("fresh_run")
    (dest / ANALYSIS).mkdir(parents=True, exist_ok=True)
    for tsv in (FGFR1_RUN / ANALYSIS).glob("*.tsv"):
        shutil.copyfile(tsv, dest / ANALYSIS / tsv.name)
    for name in ("figures_index.json", "generic/figures_index.json",
                 "generic/protein_coordinate_model.json",
                 "isoform_alignment_index.json"):
        src = FGFR1_RUN / "website_indices" / name
        if src.is_file():
            out = dest / "website_indices" / name
            out.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, out)

    # Empty the Gallery, so what remains afterwards was created by the stages.
    for name in ("figures_index.json", "generic/figures_index.json"):
        path = dest / "website_indices" / name
        if path.is_file():
            doc = json.loads(path.read_text())
            doc["figures"], doc["available"] = [], []
            path.write_text(json.dumps(doc, indent=2))

    model = dest / "website_indices" / "generic" / "protein_coordinate_model.json"
    for stage in FIGURE_STAGES:
        importlib.import_module(stage).generate(dest, model)

    doc = json.loads((dest / "website_indices" / "figures_index.json").read_text())
    return doc.get("figures") or []


def test_a_fresh_run_ships_the_complete_figure_set(fresh_gallery):
    assert len(fresh_gallery) == sum(EXPECTED_CARDS.values()), \
        f"expected {sum(EXPECTED_CARDS.values())} cards, got " \
        f"{[c['figure_id'] for c in fresh_gallery]}"


def test_a_fresh_run_fills_every_gallery_category(fresh_gallery):
    counts: dict[str, int] = {}
    for card in fresh_gallery:
        counts[card.get("category", "?")] = counts.get(card.get("category", "?"), 0) + 1
    assert counts == EXPECTED_CARDS


def test_every_card_of_a_fresh_run_is_self_describing(fresh_gallery):
    for card in fresh_gallery:
        fid = card.get("figure_id")
        assert card.get("category"), f"{fid} declares no category"
        assert card.get("caption"), f"{fid} carries no caption"
        assert card.get("scientific_question"), f"{fid} states no scientific question"
        for key in ("svg_url", "pdf_url", "png_url"):
            assert card.get(key), f"{fid} offers no {key}"


def test_a_fresh_run_registers_the_main_figures_the_gene_explorer_exports(fresh_gallery):
    main = [c for c in fresh_gallery if c["figure_id"].startswith("main_")]
    assert len(main) == 8, f"expected 8 shared main figures, got {len(main)}"
    for card in main:
        assert card.get("renderer") == "shared_figure_specification", \
            f"{card['figure_id']} was not drawn by the shared renderer"
