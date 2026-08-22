# Quarantined ancestral reconstruction

This directory preserves an earlier combined LOCO/ancestral-sequence-reconstruction implementation and its outputs for audit purposes. The ancestral results are invalid and must not be cited as scientific results.

The reconstruction attempted to name an ancestral IIIb/IIIc split from an unrooted gene-tree topology. That topology does not identify the required directional ancestral node, so the labels assigned to the reconstructed nodes are not biologically justified. Re-running the same calculation would reproduce numbers without correcting the identifiability problem.

The valid leave-one-clade-out component has been separated, re-run independently, and released under `results/08_loco_validation`. A future ancestral analysis would require an explicit, justified rooting or species-tree reconciliation strategy, a documented target-node definition, and sensitivity analyses to alternative roots and topologies.

Files in this directory are excluded from release verification and from all active result summaries.
