#!/usr/bin/env bash
set -euo pipefail
iqtree3 -s cassette_complete_pairs.safe_ids.faa -st AA -m MFP -B 1000 --alrt 1000 -T AUTO -seed 107252 --prefix unconstrained
BEST_MODEL=$(awk '/Best-fit model according to BIC:/ {print $NF}' unconstrained.iqtree)
iqtree3 -s cassette_complete_pairs.safe_ids.faa -st AA -m "$BEST_MODEL" -g constraint_isoform_monophyly.nwk -T AUTO -seed 995853 --prefix isoform_constraint
iqtree3 -s cassette_complete_pairs.safe_ids.faa -st AA -m "$BEST_MODEL" -g constraint_species_pairs.nwk -T AUTO -seed 460166 --prefix species_pair_constraint
cat unconstrained.treefile isoform_constraint.treefile species_pair_constraint.treefile > candidate_topologies.trees
iqtree3 -s cassette_complete_pairs.safe_ids.faa -st AA -m "$BEST_MODEL" -n 0 -z candidate_topologies.trees -zb 10000 -au -seed 160107 --prefix topology_AU_test
