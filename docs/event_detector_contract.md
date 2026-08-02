# Event-detector output contract (v1)

This document defines the **gene/event-agnostic contract** that any event
detector must satisfy. It is the boundary between *biological event detection*
(gene-specific) and the *generic framework* (gene-agnostic indices, views, API).

> **This is the OPTIONAL event layer.** The framework separates two analysis
> layers:
>
> - **Core gene analysis** — generic, runs for any protein-coding gene, does
>   **not** require an event region. See `docs/core_gene_analysis_contract.md`.
> - **Event-specific analysis** (this contract) — only applies when a gene has a
>   *configured or detected* event region (e.g. FGFR2 IIIb/IIIc).
>
> A gene with **no** event region is **not** a failure: it still gets core
> gene-model / domain / synteny / all-exon-boundary results. FGFR2 IIIb/IIIc is
> the validated event-specific case study layered on top of core outputs.
>
> Event regions can arrive three ways: **configured** (a supported detector, e.g.
> FGFR2), **user_defined** (a user supplies `reference_protein` +
> `region_start_aa`/`region_end_aa`; planned — the feasibility probe already
> reports `event_user_defined_possible`), or **exploratory candidates** from the
> optional `scan_isoform_event_candidates.py` scan. Only *configured* /
> *user_defined* regions enable event-specific views; exploratory candidates are
> shown as candidates only and never as validated events.

> **Why this matters.** Today the biological event detection is FGFR2-specific.
> By defining a stable output contract, the existing FGFR2 IIIb/IIIc logic can be
> *wrapped* (projected) into a generic shape, and a future gene becomes runnable
> by implementing a detector that produces the *same* outputs — without touching
> the generic layer. The machine-readable version lives at
> `configs/framework/event_detector_contract.yaml`.

## Principles

- A detector **produces files**; it does not need to expose internal logic.
- The generic index builder and all downstream views read **only** these files.
- FGFR2 IIIb/IIIc is the first detector, provided as a **projection adapter** over
  the existing validated outputs (`scripts/adapters/fgfr2_event_detector_adapter.py`).
  It recomputes **nothing** and changes **no** biology.
- Nothing here writes into the example freeze
  (`results/final_30_until_interpro_prepare/`).

## Output location

- For a custom run: `runs/<run_id>/results/generic_event_detector/`
- For the example dataset preview: a **safe artifact** folder only
  (e.g. `artifacts/generic_event_detector/example/`) — never the freeze.

## Required outputs

### 1. `event_isoform_candidates.tsv`

One row per `(species, event label)` candidate protein/isoform.

| column | required | meaning |
| --- | --- | --- |
| `analysis_id` | yes | e.g. `FGFR2_IIIb_IIIc` |
| `gene_symbol` | yes | e.g. `FGFR2` |
| `species_id` | yes | canonical species id, e.g. `gallus_gallus` |
| `event_id` | yes | e.g. `FGFR2_IIIb_IIIc_cassette` |
| `event_label` | yes | e.g. `IIIb` / `IIIc` (generic label id) |
| `protein_id` | yes | protein accession |
| `transcript_id` | no | transcript accession |
| `gene_id` | no | gene id if known |
| `sequence_id` | no | stable sequence hash/id |
| `protein_length` | no | integer aa length |
| `candidate_status` | yes | detector's status (e.g. `primary_claim_supported…`) |
| `evidence_level` | no | coarse class (e.g. `primary`, `minor`, `review`) |
| `source` | yes | detector/source label |
| `notes` | no | free text |

### 2. `event_region_coordinates.tsv`

Native protein-coordinate span of the event region per candidate.

| column | required | meaning |
| --- | --- | --- |
| `analysis_id`, `gene_symbol`, `species_id`, `event_id`, `event_label`, `protein_id` | yes | keys |
| `region_start_aa` | no | 1-based start (aa) |
| `region_end_aa` | no | 1-based end (aa) |
| `region_length_aa` | no | length (aa) |
| `coordinate_basis` | yes | how coordinates were derived, e.g. `native_protein_aa` |
| `coordinate_confidence` | no | e.g. `high`/`medium`/`low` |
| `source_table` | no | provenance file |
| `notes` | no | free text |

### 3. `event_label_evidence.tsv`

Per `(species, event label)` marker / label evidence and validation.

| column | required | meaning |
| --- | --- | --- |
| `analysis_id`, `gene_symbol`, `species_id`, `event_id`, `event_label`, `protein_id` | yes | keys |
| `marker_sequence` | no | discriminating marker peptide from the gene config |
| `marker_status` | no | e.g. `present`, `reference`, `not_evaluated` |
| `sequence_evidence` | no | detector's sequence-level evidence summary |
| `reference_similarity` | no | similarity to the curated reference (e.g. discriminating-column agreement) |
| `validated_label` | yes | final label after the detector's reconciliation |
| `confidence` | no | coarse confidence class |
| `source` | yes | detector/source label |
| `notes` | no | free text |

### 4. `event_detector_report.json`

Run-level provenance and counts:

```json
{
  "detector_name": "fgfr2_iiib_iiic",
  "detector_version": "1.0",
  "contract_version": 1,
  "analysis_id": "FGFR2_IIIb_IIIc",
  "gene_config": "configs/genes/FGFR2_IIIb_IIIc.yaml",
  "dataset_id": "run:<run_id>",
  "input_files": ["…"],
  "output_files": ["…"],
  "n_species": 1,
  "n_candidate_proteins": 2,
  "n_accepted_primary_proteins": 2,
  "warnings": [],
  "failures": [],
  "reference_control_used": true
}
```

## What these outputs enable

The generic layer can build, from these files alone:

- **event region / cassette view** (`event_region_index.json`)
- **domain architecture** (protein ids + InterPro/TM annotation)
- **boundary relation** (event region vs domain boundaries)
- **synteny context** (gene locus neighbourhood)
- **generic website indices** (`scripts/adapters/fgfr2_to_generic_indices.py`)

## Mapping FGFR2 IIIb/IIIc onto the contract

| generic field | FGFR2 source |
| --- | --- |
| `event_label` | `final_isoform_label` (`IIIb` / `IIIc`) |
| `candidate_status` | `final_claim_status_after_rescue` |
| `evidence_level` | `readiness_class` / `interpro_included` |
| `region_*_aa` | `native_cassette_start_aa` / `native_cassette_end_aa` (MSA boundary projection) |
| `marker_sequence` | gene config `event.labels[].marker_sequence` (`SGINSSN` / `GVNTTDKEI`) |
| `reference_similarity` | discriminating-column agreement vs human (cassette residue index) |
| `validated_label` | `final_isoform_label` |

## Implementing a detector for a future gene

1. Produce the four outputs above for your gene/event.
2. Register the detector in `configs/framework/event_detectors.yaml` with
   `status: supported`.
3. Ensure the generic index adapter supports your event type.
4. Validate on the example + a custom run.
5. Only then enable the analysis in the UI.

A gene with a config but **no supported detector** stays `not runnable`.
