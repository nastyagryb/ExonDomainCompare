# Final pre-InterProScan results summary

**Run ID:** `20260623T102528Z`  
**Pipeline closure version:** 1.0  
**Generated:** 2026-06-23T11:18:24Z  
**Cross-table consistency gate:** PASS  
**full_clean_run_completed:** true  
**used_cached_v3_outputs:** false  
**used_cached_msa_outputs:** false

> Run mode: **FULL CLEAN END-TO-END RUN** — Steps 1-11, MSA/rescue/synteny and closure all executed (no cached stages).

## Dataset freeze statement

The pre-InterProScan dataset is frozen after label reconciliation, maximal suspicious-case rescue,
coordinate validation, MSA robustness analysis, and local synteny / gene-neighborhood validation.
Domain evidence is **not** inferred before InterProScan. Successfully rescued and validated cases are
included as final accepted primary cases; rescue provenance is retained in tables and reports.
Unresolved or non-recoverable cases are excluded from primary claims and retained only as
review/supplement where appropriate.

## Counts

| metric | value |
|---|---|
| species/isoform rows | 60 |
| primary claim (supported) | 0 |
| primary with minor flags | 58 |
| supplement/review | 2 |
| excluded | 0 |
| InterPro primary-ready | 58 |
| supplement-review-only | 2 |
| rescued candidates | 2 |
| upstream label swaps corrected | 55 |

## InterProScan submission

Submit **`freeze/final_pre_interpro_proteins_primary.faa`** first (58 sequences).  
Use **`freeze/final_pre_interpro_proteins_all_review_included.faa`** for exploratory/review runs.  
See **`freeze/final_pre_interpro_sequence_manifest.tsv`** for per-sequence MD5 checksums and readiness.

## Key outputs

- `final_pre_interpro_truth_table.tsv` — single source of truth
- `gates/final_pre_interpro_cross_table_consistency_gate.tsv` — must PASS
- `MSA/final_*` — frozen MSA snapshot
- `figures/` — final paper figure set (SVG/PDF/PNG)
- `archive/FGFR2_final_pre_interpro_freeze_20260623T102528Z.zip`

