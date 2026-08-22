# Software versions and execution provenance

- ExonDomainCompare release: 0.2.1
- FGFR2 case-study package: 1.0.1
- Python used for the frozen downstream analysis: 3.13.1
- IQ-TREE: 3.1.3, as recorded in the released `.iqtree` reports
- IQ-TREE model for the topology comparison: LG+G4
- IQ-TREE seeds: 107252, 995853, 460166 and 160107 for the unconstrained, isoform-constrained, species-pair-constrained and AU-test stages
- Python analysis dependencies: exact versions in `requirements.txt`
- pyTMHMM version installed by the current application setup: 1.3.6

The frozen framework snapshot contains normalized pyTMHMM output but no independent runtime-version record from its historical cluster job. Version 1.3.6 is therefore the reproducible pin for new installations, not a retroactive claim about an unrecorded historical environment.

No archival DOI is assigned in this repository. A DOI must be added to the citation metadata only after an actual archival deposit. Until then, the version tag and commit hash identify the released code state.
