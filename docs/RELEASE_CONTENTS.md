# Release contents

## Included

- Python source and compatibility entry points
- frontend and backend source
- tests and GitHub CI
- configuration templates
- installation, LRZ and troubleshooting documentation
- MIT software license and a separate research-data notice
- scientific reference files
- the checksummed, read-only FGFR2 30-species dataset
- the checksummed BCL2L1 human/mouse demonstration

## Not included

- personal configuration
- passwords, keys or MFA data
- local Python and Node.js environments
- mutable runs and generated results
- caches, logs, temporary files and screenshots
- downloaded NCBI, InterProScan or pyTMHMM working data
- internal stabilization artifacts

The `.gitignore` file enforces these exclusions. Before release, the exact Git
file list, file sizes, metadata and dataset checksums are checked with:

```bash
make release-check
```
