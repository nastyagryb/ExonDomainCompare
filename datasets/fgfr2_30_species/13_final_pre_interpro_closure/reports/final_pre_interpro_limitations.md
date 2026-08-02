# Final pre-InterProScan limitations

- No InterProScan / domain annotation yet — domain-boundary claims remain downstream.
- Synteny validates locus context only; it does not assign IIIb/IIIc.
- MSA is a robustness layer; it does not relabel isoforms.
- Some distant teleost neighbors show partial synteny (not over-penalized).
- Uncharacterized LOC neighbors may receive loose human-homology names (with % identity); these are
  not curated orthology claims.
- MCScanX block-level synteny was intentionally omitted from this build.
