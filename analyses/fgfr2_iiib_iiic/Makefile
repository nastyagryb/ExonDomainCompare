PYTHON ?= python
IQTREE ?= iqtree3
OUT ?= results/reproduced
ALIGNMENT = data/framework_snapshot/alignments/final_fgfr2_IIIb_IIIc_combined_cassette_msa.aln.faa
TRUTH = data/framework_snapshot/cohort/fgfr2_post_rescue_final_truth_table.tsv
CALLS = data/framework_snapshot/protein_annotations/interpro_ensemble_coordinate_support_calls.tsv

.PHONY: verify core iqtree cross-annotation loco manifest

verify:
	$(PYTHON) scripts/verify_release.py

core:
	$(PYTHON) scripts/analysis/01_domain_caller_bias.py --calls $(CALLS) --outdir $(OUT)/01_domain_caller_bias --bootstrap 5000
	$(PYTHON) scripts/analysis/02_robustness_surface.py --calls $(CALLS) --outdir $(OUT)/02_robustness_surface
	$(PYTHON) scripts/analysis/03_weighted_jsd_permutation.py --alignment $(ALIGNMENT) --truth-table $(TRUTH) --outdir $(OUT)/03_weighted_jsd --permutations 10000
	$(PYTHON) scripts/analysis/03b_select_structure_barcode.py --jsd $(OUT)/03_weighted_jsd/weighted_jsd_positions.tsv --alignment-evidence data/framework_snapshot/alignments/final_isoform_discriminating_residues.tsv --structure-map data/external/structures/structure_mapping_17_discriminating_positions.tsv --outdir $(OUT)/03b_structure_barcode_selection
	$(PYTHON) scripts/analysis/04_structure_interface_mapping.py --positions $(OUT)/03b_structure_barcode_selection/structure_mapping_17_discriminating_positions.tsv --config config/structure_mapping.yaml --pdb-dir data/external/structures --no-download --outdir $(OUT)/04_structure_mapping --permutations 10000

iqtree:
	$(PYTHON) scripts/analysis/05_iqtree_topology_tests.py --alignment $(ALIGNMENT) --truth-table $(TRUTH) --outdir $(OUT)/05_iqtree --iqtree $(IQTREE) --run

cross-annotation:
	$(PYTHON) scripts/analysis/06_cross_annotation_replication.py --manifest config/cross_annotation_species.tsv --references data/external/curated_human_FGFR2_IIIb_IIIc_cassette_reference.faa --outdir $(OUT)/06_cross_annotation

loco:
	$(PYTHON) scripts/analysis/09_loco_validation.py --outdir $(OUT)/08_loco_validation

manifest:
	$(PYTHON) scripts/build_manifest.py
