# Detailed architecture appendix: figure and code traceability

This file is the verification bridge between the detailed appendix figures and
the implementation. Paths are relative to the repository root. The figures are
descriptive views of the implementation; the referenced source modules and tests
remain the authoritative executable evidence.

## Recommended thesis placement

| Figure | Recommended placement | Purpose in the thesis |
| --- | --- | --- |
| 02 Detailed End-to-End Analysis Pipeline | Appendix A; referenced from Section 3.2 | Complete execution sequence from run creation to canonical outputs |
| 03 Detailed Source-Code Component and Dependency Map | Appendix A; referenced from Section 3.2 | Mapping from architectural responsibilities to source packages |
| 05 Detailed C4 Container and Deployment View | Appendix A; referenced from Section 3.2 | Local application, storage, external resources and explicit LRZ boundary |
| 08 Detailed Post-InterPro Evidence and Boundary Pipeline | Appendix B; referenced from Section 3.3 | Scientifically critical annotation roles, boundary calculation and uncertainty handling |
| 10 Detailed Run, Storage and State Contracts | Appendix A; referenced from Sections 3.2 and 3.3 | Portable layout, registry, legacy boundary, lifecycle and availability semantics |
| 11 Detailed Web, API, Selection and Export Integration | Appendix C; referenced from Section 3.5 | Dataset-scoped API, linked selection, views and reproducible exports |

The four compact thesis figures remain in the main text. The detailed figures
should normally occupy landscape appendix pages and be cited from the relevant
main-text section with a short sentence such as: "A module-level representation
and the corresponding implementation traceability are provided in Appendix A."

## Figure 02 - End-to-end pipeline

| Diagram area | Primary implementation evidence | Representative tests |
| --- | --- | --- |
| Run creation and routing | `src/exondomaincompare/framework/analysis_router.py`; `production_contract.py`; `gene_config.py`; `src/exondomaincompare/runs/layout.py`; `registry.py` | `tests/test_generic_gene_run_creation.py`; `test_fgfr2_future_run_production_wiring.py`; `test_run_management.py` |
| Source models and coordinate projection | `framework/run_core_gene_analysis.py`; `shared_gene_analysis/gene_locus_resolution.py`; `model_recovery.py`; `protein_coordinate_model.py`; `strand.py` | `test_gene_locus_resolution.py`; `test_protein_coordinate_model.py`; `test_species_taxon_and_model_recovery.py` |
| Pre-/post-cluster coordinate audit | `framework/coordinate_evidence_register.py`; `framework/run_core_gene_analysis.py`; `adapters/fgfr2_core_analysis_adapter.py` | `test_coordinate_evidence_register.py` |
| Pre-cluster comparative analysis | `framework/scan_isoform_event_candidates.py`; `shared_gene_analysis/comparative_dataset.py`; `indices/msa.py`; `synteny_contract.py` | `test_multispecies_pipeline_repair.py`; `test_shared_synteny_contract.py`; `test_isoform_alignment_figures.py` |
| Cluster round trip | `src/exondomaincompare/cluster/ssh_common.py`; `src/exondomaincompare/cli.py`; compatibility scripts under `scripts/interpro_cluster/` | `test_multispecies_roundtrip_integration.py`; `test_no_machine_local_paths.py` |
| Post-cluster products and publication utilities | `framework/interpro_annotations.py`; `build_core_gene_indices.py`; `shared_gene_analysis/boundary_observations.py`; `comparative_dataset.py`; `shared_gene_analysis/package_builder.py`; `runs/outputs.py` | `test_post_interpro_generic.py`; `test_interpro_annotation_layers.py`; `test_comparative_boundary_mapping.py`; `test_comparative_gallery_and_packages.py` |

## Figure 03 - Source-code components

| Component group | Primary implementation evidence | Representative tests |
| --- | --- | --- |
| Configuration and entry points | `src/exondomaincompare/config.py`; `config_cli.py`; `cli.py`; `cluster_setup.py` | `test_no_machine_local_paths.py`; `test_release_e3_setup.py` |
| Run persistence | `src/exondomaincompare/runs/layout.py`; `registry.py`; `legacy.py`; `migration.py`; `outputs.py` | `test_run_management.py`; `test_cross_pipeline_release_repair.py`; `test_no_machine_local_paths.py` |
| Framework and scientific modules | `src/exondomaincompare/framework/`; `shared_gene_analysis/`; `generic_gene/`; `scientific/`; `framework/coordinate_evidence_register.py` | `test_analysis_availability_semantics.py`; `test_boundary_dashboard.py`; `test_single_species_scientific_validation.py`; `test_coordinate_evidence_register.py` |
| FGFR2 adapters | `src/exondomaincompare/adapters/`; FGFR2 validation scripts under `scripts/` | `test_classify_fgfr2_IIIb_IIIc_by_exon_structure_v2.py`; `test_select_fgfr2_transcripts_annotation_aware_v2.py`; `test_fgfr2_boundary_catalogue_correction.py` |
| Application and presentation | `webapp/backend/`; `webapp/frontend/src/`; `src/exondomaincompare/presentation/` | `test_backend_api_v1.py`; `test_final_ui_cleanup.py`; `test_publication_plot_system.py` |

## Figure 05 - Containers and deployment

| Container or boundary | Primary implementation evidence | Representative tests |
| --- | --- | --- |
| React frontend and FastAPI backend | `webapp/frontend/src/App.jsx`; `api.js`; `webapp/backend/main.py` | `test_backend_api_v1.py`; `test_final_ui_cleanup.py` |
| Local scientific computation | `src/exondomaincompare/framework/`; `shared_gene_analysis/`; `presentation/` | `test_post_interpro_generic.py`; `test_plotting_api.py` |
| Run-local coordinate audit artifacts | `framework/coordinate_evidence_register.py`; `results/core_gene_analysis/evidence_register/` | `test_coordinate_evidence_register.py` |
| Portable storage and configuration | `src/exondomaincompare/config.py`; `runs/` | `test_no_machine_local_paths.py`; `test_run_management.py` |
| LRZ cluster boundary | `src/exondomaincompare/cluster/ssh_common.py`; `cluster_setup.py`; `scripts/interpro_cluster/` | `test_multispecies_roundtrip_integration.py`; `test_no_machine_local_paths.py` |

## Figure 08 - Post-InterPro evidence

| Evidence stage | Primary implementation evidence | Representative tests |
| --- | --- | --- |
| InterProScan ingestion and annotation roles | `framework/interpro_annotations.py` | `test_interpro_annotation_layers.py`; `test_post_interpro_generic.py` |
| Boundary identity, distance and classification | `shared_gene_analysis/boundary_classification.py`; `boundary_observations.py`; `boundary_dashboard.py` | `test_domain_instance_boundary_repair.py`; `test_boundary_dashboard.py` |
| Additive coordinate evidence register | `framework/coordinate_evidence_register.py`; `framework/run_core_gene_analysis.py`; `adapters/fgfr2_core_analysis_adapter.py` | `test_coordinate_evidence_register.py` |
| Comparative mapping | `shared_gene_analysis/comparative_dataset.py`; `msa_coordinates.py`; `species_order.py` | `test_comparative_boundary_mapping.py`; `test_comparative_boundary_explorer.py`; `test_species_order.py` |
| Availability and freshness | `shared_gene_analysis/analysis_availability.py`; `cluster_output_freshness.py`; `run_availability.py` | `test_analysis_availability_semantics.py`; `test_run_availability_and_dependency_invalidation.py` |
| Figure and table parity | `webapp/frontend/src/pages/viewers/figureSpec.js`; `figureExport.js`; `semanticStyles.js`; `src/exondomaincompare/presentation/` | `test_figure_parity_contract.py`; `test_figure_format_parity.py`; `test_publication_plot_system.py` |

## Figure 10 - Run, storage and state contracts

| Contract | Primary implementation evidence | Representative tests |
| --- | --- | --- |
| Runtime configuration and application paths | `src/exondomaincompare/config.py` | `test_no_machine_local_paths.py`; `test_release_e3_setup.py` |
| Registry and collision handling | `src/exondomaincompare/runs/registry.py` | `test_run_management.py`; `test_external_run_status_regression.py` |
| Canonical layout and legacy adapter | `src/exondomaincompare/runs/layout.py`; `legacy.py`; `migration.py` | `test_run_management.py`; `test_cross_pipeline_release_repair.py`; `test_no_machine_local_paths.py` |
| Lifecycle inference | `framework/core_run_milestones.py`; `webapp/frontend/src/runStates.js`; `pages/runworkflow/runStatus.js` | `test_external_run_status_regression.py`; `test_run_management.py` |
| Analysis availability | `shared_gene_analysis/analysis_availability.py`; `run_availability.py` | `test_analysis_availability_semantics.py`; `test_run_availability_and_dependency_invalidation.py` |
| Coordinate-record QC and provenance | `framework/coordinate_evidence_register.py`; `configs/framework/core_gene_analysis_contract.yaml` | `test_coordinate_evidence_register.py` |
| Output and package identity | `src/exondomaincompare/runs/outputs.py`; `shared_gene_analysis/package_builder.py` | `test_comparative_gallery_and_packages.py`; `test_gallery_downloads_correction.py` |

## Figure 11 - Web, API, selection and exports

| Integration area | Primary implementation evidence | Representative tests |
| --- | --- | --- |
| Dataset loading and navigation | `webapp/frontend/src/App.jsx`; `datasetStatus.js`; `runStates.js` | `test_final_ui_cleanup.py`; `test_external_run_status_regression.py` |
| Dataset-scoped API client | `webapp/frontend/src/api.js` | `test_backend_api_v1.py`; `test_final_ui_cleanup.py` |
| FastAPI endpoint families | `webapp/backend/main.py`; `canonical_dataset.py` | `test_backend_api_v1.py`; `test_dataset_model_loading.py` |
| Coordinate-register provenance catalogue | `webapp/backend/canonical_dataset.py`; `framework/coordinate_evidence_register.py` | `test_coordinate_evidence_register.py`; `test_dataset_model_loading.py` |
| Linked scientific selection | `webapp/frontend/src/components/ScientificSelectionContext.jsx` | `test_comparative_boundary_explorer.py`; `test_comparative_exon_domain_architecture.py` |
| Scientific viewers | `webapp/frontend/src/pages/`; `pages/viewers/` | `test_comparative_boundary_explorer.py`; `test_gallery_figure_redesign.py`; `test_isoform_alignment_figures.py` |
| Rendering, downloads and manifests | `pages/viewers/figureSpec.js`; `figureExport.js`; `plotExport.js`; `src/exondomaincompare/runs/outputs.py`; `shared_gene_analysis/package_builder.py` | `test_figure_parity_contract.py`; `test_figure_format_parity.py`; `test_comparative_gallery_and_packages.py` |

## Reproducibility note

The D2 source files under `docs/architecture/d2/` are the editable source of
truth for the diagrams. SVG and PDF files are rebuilt with the D2 and
`rsvg-convert` commands documented in `docs/architecture/README.md`. When the
implementation changes, the traceability tables and the corresponding D2 source
must be reviewed together.

The atlas was audited against ExonDomainCompare `main` at commit `62a29a9`; the
coordinate-evidence additions originate from commit `afda7ae` (`Add run-level
coordinate evidence register`). The register
is an additive audit projection: source-table values and classifications remain
authoritative, while the register supplies model linkage, explicit record-level QC,
source-row pointers, run reports and checksums.
