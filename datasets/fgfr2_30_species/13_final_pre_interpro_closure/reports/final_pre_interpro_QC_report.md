# Final pre-InterProScan QC report

Run `20260623T102528Z`. Cross-table consistency gate: **PASS**  
(`gates/final_pre_interpro_cross_table_consistency_gate.tsv`).

Hard gates: post-rescue cross-table consistency, synteny validation, and final closure cross-table gate.

Primary figures include only `primary_claim_supported` or `primary_claim_supported_with_minor_flags`
rows per `final_pre_interpro_truth_table.tsv`. Rescued-and-validated cases (Gorilla, Canis IIIb,
Pongo IIIc) appear as accepted primary cases in main figures; rescue provenance is in tables/reports
and Figure 4, not repeated as problem styling in architecture/MSA/synteny figures.

Review-case explanation (`tables/final_review_case_explanation.tsv`): **Pongo abelii IIIb** and
**Canis lupus familiaris IIIc** remain supplement/review. Their locus, orthology, synteny, MSA,
coordinates and protein integrity all pass, but no source-compatible externally validated
isoform-specific candidate was recovered (sequence support only), so the isoform-specific claim is
retained as review with provenance rather than asserted as primary. These are genuine unresolved
cases, **not** rescued ones (no row is simultaneously rescued and unresolved). Their isoform partners
**Pongo IIIc** and **Canis IIIb** are confirmed/rescued and appear as accepted primary rows.
