# Gene analysis configurations

Most users do not need to edit this folder. New exploratory runs are created in
the website by entering a gene symbol and one or more species.

The configuration layer distinguishes two kinds of analysis:

- Generic core analysis works with supported protein-coding genes and reports
  gene models, protein isoforms, exploratory Candidate evidence, domains,
  synteny and exon-domain boundaries when the required input is available.
- Event-specific analysis requires a separately configured and validated
  biological event. FGFR2 IIIb/IIIc is the validated event-specific case.

Exploratory Candidates are protein-isoform difference evidence. They are not
validated splicing events.

## Files

| File | Purpose |
| --- | --- |
| `FGFR2_IIIb_IIIc.yaml` | Active configuration for the validated FGFR2 IIIb/IIIc analysis |
| `TEMPLATE_gene_event.yaml` | Schema example for a future event-specific analysis |
| `drafts/*.yaml` | Development examples and pilots; they are not validated event analyses |

Adding a YAML file does not validate a biological event and does not make it an
active event-specific analysis. A detector, adapter coverage and tests are also
required.

## Validate a configuration

From the repository root:

```bash
.venv/bin/python -m exondomaincompare.framework.gene_config \
  --config configs/genes/FGFR2_IIIb_IIIc.yaml \
  --validate
```

## Run workflow

Create a run in the website. When `My Runs` shows that cluster annotation is
required, run the displayed command from the repository root:

```bash
.venv/bin/edc cluster roundtrip --run-id RUN_ID
```

LRZ must be configured once before the first cluster run. See
[`docs/LRZ.md`](../../docs/LRZ.md).

## Source locations

- Configuration loading and validation:
  `src/exondomaincompare/framework/gene_config.py`
- Generic core analysis:
  `src/exondomaincompare/framework/run_core_gene_analysis.py`
- Exploratory Candidate scan:
  `src/exondomaincompare/framework/scan_isoform_event_candidates.py`
- Generic index builder:
  `src/exondomaincompare/framework/build_core_gene_indices.py`
- FGFR2 projection adapters: `src/exondomaincompare/adapters/`
- Event-detector registry: `configs/framework/event_detectors.yaml`

The file contracts are documented in
[`docs/core_gene_analysis_contract.md`](../../docs/core_gene_analysis_contract.md)
and
[`docs/event_detector_contract.md`](../../docs/event_detector_contract.md).

## Biological scope

The bundled FGFR2 30-species dataset is read-only. Custom runs write only under
the configured user data directory. Generic confidence scores rank exploratory
evidence; they do not establish biological validation.
