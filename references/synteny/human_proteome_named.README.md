# Shared human reference proteome (synteny broad-homology naming)

`human_proteome_named.faa` is the whole **human** RefSeq proteome
(assembly GCF_000001405.40 / GRCh38.p14, taxid 9606) with `SYMBOL|protein_id`
FASTA headers. It is used ONLY as a curated HUMAN reference/control layer to name
LOC/uncharacterized synteny neighbors of custom-run species by their best human
homolog — exactly as the Example dataset does when homo_sapiens is an analysed species.

- This is HUMAN reference data (the only reusable curated reference layer).
- It is NOT a non-human Example output and carries no non-human results.
- Custom runs that do not include homo_sapiens as an analysed species reuse this file
  (see `build_human_proteome_db` fallback in
  `scripts/validate_fgfr2_local_synteny_neighborhood.py`).
- Reproducible via: `datasets download genome accession GCF_000001405.40 --include gff3,protein`
  then relabel protein.faa headers as `SYMBOL|protein_id` from genomic.gff.
