"""Canonical, gene-agnostic plotting API.

Import plotting functions from this package, not from the implementation module
``shared_gene_plots``.  The implementation is built on the validated visual
primitives in :mod:`fgfr2_plot_style`; keeping this boundary explicit lets
generic figures share the established style without changing FGFR2 event-layer
figure scripts.
"""

from .shared_gene_plots import (
    apply_style,
    figure_title,
    legend_patch,
    plot_evidence_regions_on_protein,
    plot_candidate_domain_context,
    plot_domain_architecture,
    plot_exon_domain_boundary_distribution,
    plot_gene_model_overview,
    plot_isoform_alignment_overview,
    plot_protein_exon_architecture,
    plot_synteny_neighbourhood,
    plot_transcript_exon_structure,
    save_figure_all_formats,
    shared_legend,
)

API_VERSION = "1"

__all__ = [
    "API_VERSION",
    "apply_style",
    "figure_title",
    "legend_patch",
    "plot_evidence_regions_on_protein",
    "plot_candidate_domain_context",
    "plot_domain_architecture",
    "plot_exon_domain_boundary_distribution",
    "plot_gene_model_overview",
    "plot_isoform_alignment_overview",
    "plot_protein_exon_architecture",
    "plot_synteny_neighbourhood",
    "plot_transcript_exon_structure",
    "save_figure_all_formats",
    "shared_legend",
]
