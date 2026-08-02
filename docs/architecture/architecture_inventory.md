# Architecture Inventory & Traceability

This document maps every box and major pipeline step in the Architecture Atlas
diagrams (`docs/architecture/*.pdf`) to the relevant repository locations. It is
the traceability bridge between the visual documentation and the implementation.

All paths are relative to the repository root. Where no single file corresponds
to a conceptual component, the small set of relevant files is listed. Items
marked _planned_, _optional_ or _pending cluster_ are not yet completed.

---

## Cross-cutting components

| Component | Repository location(s) |
| --- | --- |
| Analysis router (FGFR2 → validated, else shared exploratory) | `scripts/framework/analysis_router.py` (`resolve_gene_workflow`); backend `webapp/backend/main.py` (`_resolve_workflow`) |
| Gene / event configuration | `configs/genes/FGFR2_IIIb_IIIc.yaml`; `configs/genes/drafts/{FGFR1,TPM1,TP53}_core_only_pilot.yaml`; loader `scripts/framework/gene_config.py`; contracts `configs/framework/*.yaml` |
| Run lifecycle / milestones | `scripts/create_new_run.py`; `scripts/framework/core_run_milestones.py`; per-run `runs/<RUN_ID>/status.json`, `run_config.json` |
| Canonical dataset model (read-only adapter) | `webapp/backend/canonical_dataset.py` (`build_canonical_dataset_model`, `adapt_fgfr2_legacy`, `adapt_shared_run`) |

---

## 01 — Overall System Architecture

| Diagram box | Repository location(s) |
| --- | --- |
| Researcher / Web application (React frontend) | `webapp/frontend/src/` (`App.jsx`, `main.jsx`, `pages/`, `api.js`) |
| FastAPI backend | `webapp/backend/main.py` |
| Run configuration / analysis router | `scripts/framework/analysis_router.py`; `configs/genes/` |
| Generic local pipeline | `scripts/framework/run_core_gene_analysis.py`; `scripts/generic_gene/run_generic_gene_analysis.py` |
| Pre-cluster package (primary FASTA freeze) | `runs/<RUN_ID>/results/13_final_pre_interpro_closure/freeze/final_pre_interpro_proteins_primary.faa` |
| LRZ cluster (InterProScan, pyTMHMM) | `scripts/interpro_cluster/run_cluster_roundtrip.py` (+ `submit/check/fetch`) |
| Post-InterPro normalization + boundary | `scripts/framework/interpro_annotations.py`; `run_core_gene_analysis.py` `phase_post`; `build_core_gene_indices.py` |
| Canonical scientific dataset | `runs/<RUN_ID>/website_indices/` (+ `generic/`); `canonical_dataset.py` |
| Interactive exploration / figures / downloads | `webapp/frontend/src/pages/`; `scripts/plotting/`; `scripts/generic_gene/build_generic_precluster_figures.py` |
| External: NCBI Datasets, MAFFT, InterProScan, pyTMHMM | cached NCBI GFF/FASTA via `run_core_gene_analysis.find_cached_annotation`; MAFFT via `scripts/run_fgfr2_mafft_alignments.py` (`run_mafft`); InterProScan/pyTMHMM on the cluster |

---

## 02 — End-to-End Analysis Pipeline

| Stage | Repository location(s) |
| --- | --- |
| 1 Gene + species input | `run_core_gene_analysis.py` (`phase_create`, `resolve_gene_config`); `scripts/create_new_run.py` |
| 2 Gene / transcript resolution | `run_core_gene_analysis.py` (`parse_gene_models`) |
| 3 Transcript + protein retrieval (cached NCBI) | `run_core_gene_analysis.py` (`find_cached_annotation`, `load_fasta`) |
| 4 Protein-coding filtering | `run_core_gene_analysis.py` (`_collect_species_rows`) |
| 5 Primary selection | `run_core_gene_analysis.py` (`_select_primary_pids`); `scripts/framework/primary_selection.py` |
| 6 Exon / CDS extraction | `run_core_gene_analysis.py`; `scripts/collect_fgfr2_models_dual_source_v3.py` (`build_cds_features_from_parts`) |
| 7 Exon-to-protein mapping | `run_core_gene_analysis.py` (`_exon_protein_map_for_transcript`); `scripts/generic_gene/build_exon_protein_architecture.py` |
| 8 Within-species isoform comparison | `scripts/framework/scan_isoform_event_candidates.py` |
| 9 Exploratory candidate detection _(optional)_ | `scan_isoform_event_candidates.py`; `build_event_region_evidence.py`; `cluster_event_region_evidence.py` |
| 10 Within-species isoform alignment (MAFFT) | `scripts/generic_gene/build_generic_msa_index.py` |
| 11 Cross-species primary MSA _(multi-species only)_ | `scripts/generic_gene/build_generic_msa_index.py` (`primaries_msa`) |
| 12 Local synteny | `run_core_gene_analysis.py` (`extract_synteny_neighbors`); `scripts/generic_gene/build_synteny_neighbourhood.py` |
| 13 Pre-cluster freeze + package | `run_core_gene_analysis.py` (`build_core_contract`); `scripts/framework/build_core_gene_indices.py` |
| 14 Submit | `scripts/interpro_cluster/submit_cluster_analysis.py` |
| 15 InterProScan | cluster sbatch via `submit_cluster_analysis.py` |
| 16 pyTMHMM | cluster sbatch via `submit_cluster_analysis.py` |
| 17 Poll + fetch | `check_cluster_analysis.py`; `fetch_cluster_analysis.py` |
| 18 Annotation normalization | `scripts/framework/interpro_annotations.py`; `run_core_gene_analysis.py` (`_parse_interproscan`, `_protein_species_map`) |
| 19 Representative domain layer | `interpro_annotations.py` (`representative_domains`) |
| 20 Boundary calculation | `run_core_gene_analysis.py` (`_classify_boundary`); `build_core_gene_indices.py` (`NEAR_EDGE_THRESHOLD_AA`) |
| 21 Canonical dataset / index (re)generation | `build_core_gene_indices.py`; `build_generic_website_indices.py` |
| 22 Frontend presentation | `webapp/frontend/src/pages/GeneExplorer.jsx` + viewers |
| 23 Figure / report generation | `scripts/plotting/`; `build_generic_precluster_figures.py` |
| FGFR2 specialization points | `scripts/select_fgfr2_transcripts_annotation_aware_v2.py`; `scripts/classify_fgfr2_IIIb_IIIc_by_exon_structure_v2_3_human_calibrated.py`; `scripts/analyze_exon_domain_boundary_consistency.py` |

---

## 03 — Software Component Diagram

| Group / component | Repository location(s) |
| --- | --- |
| Frontend: shell, selector, Gene Explorer, viewers, pages | `webapp/frontend/src/App.jsx`, `pages/DatasetSwitcher.jsx`, `pages/GeneExplorer.jsx`, `pages/viewers/`, `pages/{Overview,BoundaryPage,FigureGallery,RunWorkflowPage}.jsx` |
| Frontend: shared selection state | `webapp/frontend/src/components/ScientificSelectionContext.jsx` |
| Frontend: API client | `webapp/frontend/src/api.js` |
| Backend: app, discovery, resolver, serialization, files | `webapp/backend/main.py` |
| Backend: canonical adapter | `webapp/backend/canonical_dataset.py` |
| Shared analysis: retrieval / parsing / selection / mapping | `scripts/framework/run_core_gene_analysis.py`; `scripts/framework/primary_selection.py`; `scripts/generic_gene/build_exon_protein_architecture.py` |
| Shared analysis: candidates | `scripts/framework/scan_isoform_event_candidates.py` + `build/cluster_event_region_evidence.py` |
| Shared analysis: alignment / synteny | `scripts/generic_gene/build_generic_msa_index.py`; `scripts/generic_gene/build_synteny_neighbourhood.py` |
| Shared analysis: InterPro / pyTMHMM normalization | `scripts/framework/interpro_annotations.py`; `run_core_gene_analysis.py` `phase_post` |
| Shared analysis: boundary + indices | `build_core_gene_indices.py`; `scripts/generic_gene/build_generic_website_indices.py`; `scripts/shared_gene_analysis/` |
| Shared analysis: figures | `scripts/plotting/`; `build_generic_precluster_figures.py` |
| Cluster round trip | `scripts/interpro_cluster/` |
| Specialized (FGFR2) layer | `scripts/collect_fgfr2_models_dual_source_v3.py`, `scripts/select_fgfr2_transcripts_annotation_aware_v2.py`, `scripts/classify_fgfr2_IIIb_IIIc_by_exon_structure_v2_3_human_calibrated.py`, `scripts/export_selected_fgfr2_proteins_complete_v2_1_region_qc.py`, `scripts/build_website_indices.py`, `scripts/analyze_exon_domain_boundary_consistency.py` |
| Persistence / artifacts | `runs/<RUN_ID>/`; frozen example `results/final_30_until_interpro_prepare/13_final_pre_interpro_closure/` |

---

## 04 — Scientific Data Flow

| Data object (backing file) | Producer |
| --- | --- |
| Gene + species config (`gene_config.yaml`, `species_list.txt`) | `phase_create` |
| Gene / transcript / CDS models (`gene_model_index.tsv`, `02_models/`) | `parse_gene_models` |
| Protein sequences (`proteins_primary.faa`, `proteins_all_isoforms.faa`) | `load_fasta`, `build_core_contract` |
| Isoforms (`protein_isoform_index.tsv`, `primary_selection_evidence.tsv`) | `_select_primary_pids`, `primary_selection.py` |
| Exon-protein map (`exon_protein_map.tsv`) | `_exon_protein_map_for_transcript` |
| Candidate regions (`event_candidate_regions.tsv`, `event_region_evidence.tsv`) | `scan_isoform_event_candidates.py` + evidence builders |
| Alignments (`isoform_msa__<species>.aln.faa`, `primaries_msa.aln.faa`) | `build_generic_msa_index.py` |
| InterProScan raw (`input.fasta.{tsv,json,gff3}` → `interpro_annotations.tsv`) | cluster + `interpro_annotations.py` |
| Normalized domain layers (`domain_features.tsv`) | `interpro_annotations.py`, `phase_post` |
| pyTMHMM topology (`tm_features.tsv`) | `phase_post` |
| Boundary rows (`exon_domain_boundary_distances.tsv`, `_summary.tsv`) | `_classify_boundary`, `phase_post` |
| Synteny neighbours (`synteny_neighbors.tsv`) | `extract_synteny_neighbors` |
| Canonical dataset (`website_indices/*.json`) | `build_core_gene_indices.py` + `canonical_dataset.py` |

---

## 05 — C4 Container Diagram

| Container | Repository location(s) |
| --- | --- |
| Web frontend | `webapp/frontend/` |
| Backend API | `webapp/backend/main.py` |
| Local pipeline / orchestrator | `scripts/framework/run_core_gene_analysis.py`, `scripts/generic_gene/` |
| Annotation + boundary layer | `scripts/framework/interpro_annotations.py`, `phase_post`, `build_core_gene_indices.py` |
| Cluster workflow | `scripts/interpro_cluster/` |
| Figure / report generator | `scripts/plotting/`, `build_generic_precluster_figures.py` |
| Run-directory storage | `runs/<RUN_ID>/` |
| External systems | NCBI Datasets; Ensembl (FGFR2 via `collect_fgfr2_models_dual_source_v3.py`); LRZ Slurm (`ssh_common.py`); InterProScan; pyTMHMM; MAFFT |

---

## 06 — Generic Framework vs FGFR2 Specialization

| Element | Repository location(s) |
| --- | --- |
| Shared foundation | `scripts/framework/`, `scripts/generic_gene/`, `scripts/shared_gene_analysis/`, `scripts/plotting/`, `webapp/` |
| Generic exploratory mode | `scripts/generic_gene/`, `scripts/framework/scan_isoform_event_candidates.py`, generic branches in `webapp/frontend/src/pages/viewers/DomainArchitecture.jsx` (`GenericDomainArchitecture`) |
| FGFR2 validated use case | `configs/genes/FGFR2_IIIb_IIIc.yaml`; `scripts/*fgfr2*`; FGFR2 branches in `DomainArchitecture.jsx`, `CoordinateTrack.jsx`, `MsaExplorer.jsx`; `scripts/build_website_indices.py` |
| Routing decision | `scripts/framework/analysis_router.py` |

---

## 07 — Researcher Interaction Workflow

| Element | Repository location(s) |
| --- | --- |
| Journey (run → gene/species → views → export) | `webapp/frontend/src/App.jsx`; `pages/RunWorkflowPage.jsx`; `pages/GeneExplorer.jsx` |
| Shared selection (species / protein / candidate) | `webapp/frontend/src/components/ScientificSelectionContext.jsx` |
| Coordinated views | `pages/viewers/{CoordinateTrack,DomainArchitecture,MsaExplorer,SyntenyViewer}.jsx`; `pages/BoundaryPage.jsx` |
| Multi-species MSA vs within-species alignment | `scripts/generic_gene/build_generic_msa_index.py`; `pages/viewers/MsaExplorer.jsx` |

---

## 08 — Post-InterPro Annotation & Boundary Pipeline

| Step | Repository location(s) |
| --- | --- |
| Primary-protein FASTA | `13_final_pre_interpro_closure/freeze/final_pre_interpro_proteins_primary.faa` |
| InterProScan raw / pyTMHMM results | fetched to `results/14_interproscan/primary/output/`, `results/15_exon_domain_boundary_post_interpro/pytmhmm_primary/output/` |
| Protein-ID normalization + ownership | `run_core_gene_analysis.py` (`_parse_interproscan`, `_protein_species_map`) |
| Layer separation (domain/family/feature/raw/topology) | `scripts/framework/interpro_annotations.py` (`layer_for`, `representative_domains`, `family_annotations`, `feature_annotations`) |
| Nearest-domain-edge + classification | `run_core_gene_analysis.py` (`_classify_boundary`); `build_core_gene_indices.py` (`NEAR_EDGE_THRESHOLD_AA = 5`) |
| Consumers: Domain Architecture / Boundary summary / global page | `webapp/frontend/src/pages/viewers/DomainArchitecture.jsx`; `pages/GeneExplorer.jsx`; `pages/BoundaryPage.jsx`; indices `domain_architecture_index.json`, `exon_domain_boundary_index.json` |

---

## 09 — Overview Poster

Composite of 01–08; see the mappings above. Source
`docs/architecture/sources/09_exondomaincompare_architecture_poster.d2`.

---

## Tests

| Area | Test file |
| --- | --- |
| InterPro annotation layers, representative domains, boundary scope, FGFR2 preservation | `tests/test_interpro_annotation_layers.py` |
| Generic post-InterPro flow, species assignment, boundary/TM, FGFR2 cassette preservation | `tests/test_post_interpro_generic.py` |
| Canonical dataset model (FGFR2 legacy + shared run) | `tests/test_canonical_dataset.py` |
| Plotting API contract | `tests/test_plotting_api.py` |
| Backend API smoke | `tests/test_backend_api_v1.py` |

---

## Build

All diagrams and the guide are produced from the sources in
`docs/architecture/sources/` and `docs/architecture/guide/` by:

```bash
python scripts/docs/build_architecture_atlas.py        # or: make architecture
```

Toolchain: `d2`, `rsvg-convert` (librsvg), `typst` — installed via
`brew install d2 librsvg typst`.
