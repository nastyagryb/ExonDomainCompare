# Canonical plotting API

New gene-agnostic figures must import `scripts/plotting` and call its public
functions. `scripts/plotting/shared_gene_plots.py` is the implementation module;
`src/exondomaincompare/presentation/fgfr2_plot_style.py` remains the single source for palette, typography,
line weights, gene arrows, legends, titles, and multi-format export.

Public API version 1:

- `apply_style`
- `plot_gene_model_overview`
- `plot_protein_exon_architecture`
- `plot_synteny_neighbourhood`
- `plot_evidence_regions_on_protein`
- `figure_title`, `shared_legend`, `legend_patch`
- `save_figure_all_formats`

The generic pre-cluster builder imports only this package API. It contains data
adaptation and manifest handling, but no Matplotlib drawing.

## Optional validated FGFR2 event-layer figures

The following existing scripts are validated FGFR2-specific event-layer
figures, not part of the gene-agnostic API:

- `scripts/make_fgfr2_post_interpro_exon_domain_figures.py`
- `scripts/make_fgfr2_synteny_figures_paper.py`
- `scripts/make_fgfr2_final_closure_figures.py`
- `scripts/make_fgfr2_final_framework_figure.py`

They remain unchanged to preserve validated output. They import
`src/exondomaincompare/presentation/fgfr2_plot_style.py` directly and therefore run on the same base visual
primitives as the canonical generic API. Migrating them through the package API
is optional and requires a separate pixel-identity check against validated
reference artifacts.

Callers must never direct either plotting layer to a freeze directory. The
generic builder enforces this through `GenericContext.assert_not_freeze()`
before creating any figure directory.
