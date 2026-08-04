from __future__ import annotations

import ast
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import plotting  # noqa: E402


FIGURE_FUNCTIONS = {
    "plot_gene_model_overview",
    "plot_transcript_exon_structure",
    "plot_isoform_alignment_overview",
    "plot_protein_exon_architecture",
    "plot_synteny_neighbourhood",
    "plot_evidence_regions_on_protein",
    "plot_domain_architecture",
    "plot_exon_domain_boundary_distribution",
    "plot_candidate_domain_context",
}
BUILDER_FIGURE_FUNCTIONS = FIGURE_FUNCTIONS - {"plot_gene_model_overview"}


def _tree(relative_path: str) -> ast.AST:
    return ast.parse((ROOT / relative_path).read_text(encoding="utf-8"))


def _attribute_calls(tree: ast.AST, owner: str) -> set[str]:
    return {
        node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id == owner
    }


def test_package_exposes_versioned_canonical_api():
    assert plotting.API_VERSION == "1"
    expected = FIGURE_FUNCTIONS | {
        "apply_style",
        "figure_title",
        "shared_legend",
        "legend_patch",
        "save_figure_all_formats",
    }
    assert expected <= set(plotting.__all__)
    assert all(callable(getattr(plotting, name)) for name in expected)


def test_generic_builder_uses_only_canonical_plotting_api():
    tree = _tree("scripts/generic_gene/build_generic_precluster_figures.py")
    imports = [
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.Import, ast.ImportFrom))
    ]
    assert any(
        isinstance(node, ast.Import)
        and any(alias.name == "plotting" and alias.asname == "plots" for alias in node.names)
        for node in imports
    )
    assert not any(
        isinstance(node, ast.ImportFrom)
        and (
            node.module == "plotting.shared_gene_plots"
            or any(alias.name == "shared_gene_plots" for alias in node.names)
        )
        for node in imports
    )
    assert not any(
        isinstance(node, ast.Import)
        and any(alias.name in {"matplotlib", "matplotlib.pyplot", "fgfr2_plot_style"}
                for alias in node.names)
        for node in imports
    )
    assert _attribute_calls(tree, "plots") == BUILDER_FIGURE_FUNCTIONS | {"apply_style"}


def test_shared_figures_use_common_helpers_and_fgfr2_primitives():
    tree = _tree("src/exondomaincompare/presentation/shared_gene_plots.py")
    calls = _attribute_calls(tree, "ps")
    assert {
        "apply_rcparams",
        "compact_legend",
        "gene_arrow",
        "legend_patch",
        "savefig",
        "title",
    } <= calls

    functions = {
        node.name: node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    for name in FIGURE_FUNCTIONS:
        local_calls = {
            node.func.id
            for node in ast.walk(functions[name])
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }
        assert {"figure_title", "shared_legend", "save_figure_all_formats"} <= local_calls


def test_validated_event_layer_scripts_stay_on_same_base_primitives():
    event_scripts = [
        "scripts/make_fgfr2_post_interpro_exon_domain_figures.py",
        "scripts/make_fgfr2_synteny_figures_paper.py",
        "scripts/make_fgfr2_final_closure_figures.py",
        "scripts/make_fgfr2_final_framework_figure.py",
    ]
    for path in event_scripts:
        tree = _tree(path)
        assert any(
            isinstance(node, ast.ImportFrom)
            and node.module == "exondomaincompare.presentation"
            and any(alias.name == "fgfr2_plot_style" for alias in node.names)
            for node in ast.walk(tree)
        ), path


def test_canonical_api_writes_all_formats(tmp_path):
    plotting.apply_style()
    written = plotting.plot_gene_model_overview(
        tmp_path,
        "smoke",
        gene_symbol="GENE",
        isoforms=[
            {"protein_id": "P1", "protein_length": 100, "primary_status": "primary"},
            {"protein_id": "P2", "protein_length": 80, "primary_status": "alternative"},
        ],
    )
    assert written is True
    assert all((tmp_path / f"smoke.{ext}").is_file() for ext in ("svg", "pdf", "png"))
