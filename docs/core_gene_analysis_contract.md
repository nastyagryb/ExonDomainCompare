# Core gene-analysis output contract (v1)

This document defines the **core, gene/event-agnostic** analysis contract. It is
deliberately **separate** from the event-detector contract
(`docs/event_detector_contract.md`).

## Two layers

| Layer | Applies to | Requires an event region? |
| --- | --- | --- |
| **Core gene analysis** | any supported protein-coding gene with collectable gene models/proteins | **No** |
| **Event-specific analysis** | genes with a configured/detected alternative or event region | **Yes** |

> **Key principle:** *"No event region"* is **not** a failure. A gene without a
> known cassette/event can still produce useful gene-model, domain-architecture,
> synteny and **all-exon boundary** results through the core contract. FGFR2
> IIIb/IIIc remains the validated *event-specific* case study layered on top of
> core outputs.

The machine-readable version is `configs/framework/core_gene_analysis_contract.yaml`.

## Output location

- Custom run: `runs/<run_id>/results/core_gene_analysis/`
- Example preview: a **safe artifact** folder only
  (e.g. `artifacts/core_gene_analysis/example/`) — never the freeze.

## Required outputs

### 1. `gene_model_index.tsv`
`analysis_id`, `gene_symbol`, `species_id`, `gene_id`, `transcript_id`,
`protein_id`, `source`, `protein_length`, `model_status`, `notes`.

### 2. `protein_isoform_index.tsv`
`species_id`, `protein_id`, `transcript_id`, `isoform_label` *(optional — only if
known; core analysis does not require labels)*, `protein_length`, `sequence_path`,
`primary_status`, `notes`.

### 3. `exon_protein_map.tsv`
`species_id`, `protein_id`, `transcript_id`, `exon_id`, `exon_number`, `cds_start`,
`cds_end`, `protein_start_aa`, `protein_end_aa`, `phase`, `confidence`, `source`.

### 4. `domain_features.tsv`
`species_id`, `protein_id`, `domain_source`, `domain_id`, `domain_name`,
`start_aa`, `end_aa`, `score`.

### 5. `tm_features.tsv`
`species_id`, `protein_id`, `start_aa`, `end_aa`, `source`.

### 6. `synteny_neighbors.tsv`
`species_id`, `gene_symbol`, `neighbor_symbol`, `side`, `order`, `orientation`,
`classification`, `source`, `status`.

### 7. `exon_domain_boundary_distances.tsv`  *(generic boundary analysis)*
`analysis_id`, `gene_symbol`, `species_id`, `protein_id`, `transcript_id`,
`exon_boundary_id`, `boundary_position_aa`, `nearest_domain_id`,
`nearest_domain_instance_id`, `nearest_domain_start_aa`, `nearest_domain_end_aa`,
`nearest_domain_name`, `nearest_domain_boundary_type`, `distance_aa`, `category`,
`source`.

`nearest_domain_instance_id` is `<interpro_accession>:<start_aa>-<end_aa>` and
identifies the **one domain feature instance** the distance was measured against.
An InterPro accession alone is not an identity: a repeated entry (e.g. the three
`IPR007110` Ig-like domains of FGFR1) has several instances at different
coordinates, so `nearest_domain_id` / `nearest_domain_accession` must never be used
to look a feature up.

`category` is one of:

- `exactly_aligned` — distance 0 aa to a domain boundary
- `near_boundary` — within `near_boundary_threshold_aa` (default 5 aa)
- `inside_domain` — the exon boundary falls strictly inside a domain
- `outside_annotated_domain` — no annotated domain nearby

> This is the **gene-agnostic** boundary analysis over **all** protein-coding
> exon boundaries. It does **not** need a configured event region. The FGFR2
> event-specific "Boundary Consistency" (cassette-vs-domain) view is an additional
> event-layer view and is unchanged.

### 8. `core_gene_report.json`
`analysis_id`, `gene_symbol`, `dataset_id`, `summary`, `inputs`, `outputs`,
`warnings`, `failures`, `available_views`.

### 9. `core_model_collection_report.json` *(collection provenance)*
Written by the gene-agnostic collector: `input_mode`, `source` (assembly / GFF /
FASTA / taxid), `gene_locus`, `selection_method_requested`,
`selection_rule_applied`, `n_transcripts`, `n_protein_coding`,
`n_primary_proteins`, `primary_protein_ids`, `exon_map_available` +
`exon_map_reason`, `synteny_available` + `synteny_reason`, `fasta_reason`,
`warnings`. When a required output is missing, the matching `*_reason` field says
exactly why (e.g. `exon_map_unavailable`, `synteny_unavailable`,
`no_protein_coding_transcript`).

### 10. `proteins_primary.faa` / `proteins_all_isoforms.faa`
`proteins_primary.faa` is the primary protein(s) sent to the cluster. The
`proteins_all_isoforms.faa` file holds every coding isoform with an available
sequence and feeds the optional exploratory event-candidate scan (it is **not**
sent to the cluster).

## Real live runner (experimental)

`src/exondomaincompare/framework/run_core_gene_analysis.py` is a real, gene-agnostic core
runner (not only a projection of FGFR2 outputs):

```bash
# create: collect gene models/proteins, build the core contract + primary FASTA
python -m exondomaincompare.framework.run_core_gene_analysis \
    --gene-config configs/genes/drafts/FGFR1_core_only_pilot.yaml \
    --species "Gallus gallus" --selection-method canonical_if_available \
    --input-mode auto

# after the cluster round-trip fetched InterProScan/pyTMHMM outputs:
.venv/bin/edc cluster roundtrip --run-id <run_id>
# (the round-trip calls run_core_gene_analysis.py --post automatically)
```

**Primary-protein selection rules** (documented, gene-agnostic):

| `--selection-method` | Rule |
| --- | --- |
| `canonical_if_available` *(default)* | MANE Select / RefSeq Select / curated (`NM_`/`NR_`) transcript, else the longest protein |
| `longest_protein` | the longest protein by sequence length |
| `all_protein_coding_transcripts` | every coding transcript (records a warning; no single canonical isoform) |

**Input modes** (`--input-mode`):

| Mode | Meaning | Status |
| --- | --- | --- |
| `auto` *(default)* | use a local NCBI/Ensembl cache if present, otherwise (future) download the required source annotation | cache path implemented; auto-download planned |
| `local_cache` | use only an existing local annotation/protein cache | implemented |
| `user_files` | user supplies `--gff` + `--protein-faa` (later: genome FASTA) | schema + CLI wired; full UI upload planned |

The chosen mode is stored in `run_config.json` and
`core_model_collection_report.json` so `user_files` can be added later without a
schema change.

## Required milestones and the validator

Every core run is classified by **file-based milestones** (single source of truth:
`src/exondomaincompare/framework/core_run_milestones.py`, used by both the CLI validator and the
webapp):

1. **run setup** — `run_config.json`, `gene_config.yaml`, `species_list.txt`, `status.json`
2. **model collection** — `gene_model_index.tsv`, `protein_isoform_index.tsv`
3. **primary FASTA** — `proteins_primary.faa` (or a clear reason)
4. **exon map** — `exon_protein_map.tsv` (or `exon_map_unavailable`)
5. **synteny** — `synteny_neighbors.tsv` (or `synteny_unavailable`)
6. **cluster input** — `cluster_input_fasta` in `run_config.json` + one round-trip command
7. **post-domain** — `domain_features.tsv`, `tm_features.tsv`,
   `exon_domain_boundary_distances.tsv`, `core_gene_report.json`, generic indices

```bash
python -m exondomaincompare.framework.validate_core_gene_run --run-id <run_id>
```

**Honest status semantics** (a run never looks analysis-ready when required core
outputs are missing):

| Inferred status | Meaning |
| --- | --- |
| `created_not_started` | folder scaffolded, no core collection yet |
| `core_model_collection_failed` | collection ran but produced no gene models |
| `incomplete` | required core outputs (e.g. primary FASTA) missing |
| `cluster_required` | core collection complete, cluster annotation pending |
| `cluster_running` | cluster jobs submitted / running |
| `post_interpro_incomplete` | cluster outputs missing when `--post` ran (nothing fabricated) |
| `results_ready` | domain architecture + generic indices built |

## Optional, exploratory event-candidate scan

`src/exondomaincompare/framework/scan_isoform_event_candidates.py` compares a gene's protein
isoforms within a species and records candidate isoform-specific regions to
`event_candidate_regions.tsv`. These are **candidates only — not validated
events**, they never enable event-specific views, and the run succeeds with or
without the scan.

## What core analysis can show (no event required)

`overview`, `gene_models` / `isoforms`, `msa` (if alignment exists),
`domain_architecture` (after domains), `synteny` (if locus annotation available),
`exon_domain_boundaries` (after domains + exon map).

## What requires an event region

`event_region`, `event_specific_comparison`, `event_discriminating_columns`,
`event_domain_relation` (a.k.a. event-vs-domain boundary relation). See
`docs/event_detector_contract.md`.

## FGFR2 mapping

The existing FGFR2 outputs can be projected into these core files by
`src/exondomaincompare/adapters/fgfr2_core_analysis_adapter.py` (projection only — no biology
recomputed). FGFR2 additionally has the event layer (IIIb/IIIc) on top.

## Support levels

| Support level | Meaning |
| --- | --- |
| `validated_event_analysis` | FGFR2 IIIb/IIIc — runnable, validated, event-specific (core outputs **plus** the event layer) |
| `core_only_pilot` | Core gene/domain/synteny/exon-boundary analysis, **no** validated event (experimental until validated) |
| `exploratory_event_candidates` | candidate isoform-specific regions only (from the optional scan); **not** validated events |
| `draft_not_runnable` | a config exists but no working runner yet |

User-defined event regions (`event.status: user_defined` with
`reference_protein` + `region_start_aa`/`region_end_aa`) are **planned**: the
schema and feasibility probe already accept them (the probe then reports
`event_user_defined_possible`), but the interactive definition UI is not built
yet. See `configs/genes/drafts/TEMPLATE_user_defined_event.yaml`.

## Building generic indices from core outputs

`src/exondomaincompare/framework/build_core_gene_indices.py` turns the core contract TSVs into
the generic website indices, **without** any event outputs:

```bash
# from a run's core outputs
python -m exondomaincompare.framework.build_core_gene_indices --run-id <run_id>

# from an explicit core dir (e.g. the synthetic mock), no event region
python -m exondomaincompare.framework.build_core_gene_indices \
    --core-dir artifacts/core_gene_analysis/mock \
    --config configs/genes/drafts/TPM1_core_only_pilot.yaml \
    --dataset-id mock:tpm1_core_only \
    --out artifacts/generic_indices/mock_core_only
```

For a no-event config it writes `event_region_index.json` with
`available: false, reason: "no_event_configured"`, disables the event-specific
boundary view, and enables the generic Exon–Domain Boundaries view. A synthetic
no-event **mock** dataset lives under `artifacts/generic_indices/mock_core_only/`
for UI smoke testing only (not a biological result).

## What core-only mode does NOT do

Core-only mode is **experimental until validated**. It does **not**:

- claim any biological event/cassette for a non-FGFR2 gene,
- enable event-specific views (event region, event-specific boundary
  consistency, event-discriminating columns, human event comparison) unless an
  event is configured / user-defined,
- treat exploratory candidate regions as validated events, or
- modify or depend on the FGFR2 validated freeze.

FGFR2 IIIb/IIIc remains the only validated *event-specific* case study. A
core-only run for another gene is a genuine gene-level analysis (models, domains,
synteny, all-exon boundaries) but its biological interpretation is not validated.
