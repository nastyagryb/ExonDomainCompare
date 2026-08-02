# Gene / event analysis configs

This folder holds **gene/event analysis configurations** for the generalization
layer that sits *above* the validated FGFR2 pipeline.

> **Status:** `FGFR2_IIIb_IIIc` is the **only active, validated analysis**.
> The generalization layer is an *architecture preparation step*. Adding a config
> here does **not** make an arbitrary gene work automatically — a second gene
> needs its own pipeline coverage and adapter validation first.

## Two analysis layers: core vs event

The framework distinguishes two layers so that future genes do **not** require a
known FGFR2-like cassette/event region:

1. **Core gene analysis** *(generic, always possible)* — runs for any supported
   protein-coding gene where gene models/proteins can be collected. Produces gene
   models, isoforms, exon→protein maps, domain/TM features, synteny neighbours and
   the **generic all-exon exon–domain boundary** analysis. Contract:
   `docs/core_gene_analysis_contract.md`.
2. **Event-specific analysis** *(optional)* — runs **only** when a gene has a
   configured or detected event region (e.g. FGFR2 IIIb/IIIc). Adds the event
   region view, discriminating columns, reference comparison and the
   event-vs-domain boundary relation. Contract:
   `docs/event_detector_contract.md`.

> **"No event" is not a failure.** A gene without a known event region is still
> analysable through the core layer. The config expresses this with
> `analysis_modes` and `event.status`:
>
> ```yaml
> analysis_modes:
>   core_gene_analysis: true
>   event_analysis: configured   # configured | auto | user_defined | disabled | optional
> event:
>   status: configured           # configured | user_defined | not_configured
> ```
>
> - FGFR2: `event_analysis: configured`, `event.status: configured`.
> - Future gene, no known event: `event_analysis: optional`, `event.status: not_configured`.
> - User-supplied region: `event.status: user_defined` with `reference_protein` +
>   `region_start_aa` / `region_end_aa`.
>
> Event fields (`event.id`, `event.type`, labels/markers) are required **only**
> when an event region is actually configured; core-only configs validate without them.

## Files

| File | Purpose |
| --- | --- |
| `FGFR2_IIIb_IIIc.yaml` | Canonical, validated FGFR2 IIIb/IIIc analysis config (active). |
| `TEMPLATE_gene_event.yaml` | Documentation placeholder for a future gene/event config. **Not wired into the UI, not runnable.** |
| `drafts/FGFR1_draft.yaml` | Draft second-gene candidate. **`status: draft_not_runnable`, `runnable: false`.** |
| `drafts/TPM1_draft.yaml` | Second illustrative draft. **Not runnable.** |
| `drafts/TPM1_core_only_pilot.yaml` | **Core-only pilot (experimental).** Proves the core layer for a gene with **no** configured event region. `status: core_only_pilot`, `runnable: experimental`. |
| `drafts/FGFR1_core_only_pilot.yaml` | **Core-only pilot (experimental), wired to the real runner.** Used by `scripts/framework/run_core_gene_analysis.py` for a live FGFR1 core run. |
| `drafts/TEMPLATE_user_defined_event.yaml` | Design preview of a **user-defined** event region (`event.status: user_defined`). Documents the accepted schema; not runnable yet. |

## Support levels

Every analysis has one of three support levels (`GeneConfig.support_level`, also
exposed by `/api/analysis-capabilities`):

1. **`validated_event_analysis`** — runnable, validated, event-specific.
   Example: **FGFR2 IIIb/IIIc** (the only one). Requires an active config, a
   configured event region, and a supported event detector.
2. **`core_only_pilot`** — *experimental* core-only proof of concept. Core
   gene/domain/synteny/exon-boundary analysis **without** an event region. Not a
   normal runnable choice; only exposed behind an explicit "experimental" label.
   Example: `drafts/TPM1_core_only_pilot.yaml`.
3. **`draft_not_runnable`** — a config exists but there is no working runner yet.
   Examples: `drafts/FGFR1_draft.yaml`, `drafts/TPM1_draft.yaml`.

> Core-only pilots demonstrate that "no event region" is a first-class,
> non-failing case. They do **not** imply the gene is biologically validated.

## Current supported analysis

- **`FGFR2_IIIb_IIIc`** — the only runnable, validated analysis.

Draft configs under `drafts/` (e.g. `FGFR1_draft`, `TPM1_draft`) are **not runnable**
and are never offered as run choices in the UI. They only document how a future
gene would be configured.

## What a config describes

A config is a *description layer only*. It never changes biological logic and
never modifies the example freeze at `results/final_30_until_interpro_prepare/`.
It captures, in a gene/event-agnostic way:

- `analysis` — id, display name, description
- `gene` — symbol, display name, reference species
- `event` — id, type, display name, and the mutually-exclusive/alternative
  `labels` (e.g. `IIIb` / `IIIc`) with optional discriminating marker peptides
- `reference_control` — the curated human reference/control layer (reference only,
  never counted as an analysed species)
- `views` — which explore views this analysis enables
- `ui_labels` — the visible wording (e.g. FGFR2 shows "Cassette" / "Boundary Consistency")
- `canonical_outputs` — where the run's canonical artifacts live

## How the layer is consumed today

- **Loader:** `scripts/framework/gene_config.py` loads + validates a config and
  exposes a stable view (`analysis_id`, `gene_symbol`, `event_id`, `event_type`,
  `event_labels`, `reference_species`, `enabled_views`, `ui_labels`).
  Validate with:
  ```bash
  python scripts/framework/gene_config.py --config configs/genes/FGFR2_IIIb_IIIc.yaml --validate
  ```
- **Run metadata:** new runs record `analysis_id` / `event_id` / `event_type` /
  `gene_config` in `run_config.json`, and a copy of the config is written to
  `runs/<run_id>/gene_config.yaml`. Older runs without this metadata default to
  `FGFR2_IIIb_IIIc` for backward compatibility.
- **Generic indices:** `scripts/adapters/fgfr2_to_generic_indices.py` reads the
  existing FGFR2 website indices and writes a *parallel* set of gene/event-agnostic
  indices (`dataset_summary.json`, `gene_event_index.json`, `event_region_index.json`,
  `domain_architecture_index.json`, `synteny_index.json`, `boundary_relation_index.json`,
  `available_views.json`). It does **not** replace the canonical indices.
- **Frontend labels:** `webapp/frontend/src/labels.js` (`getDatasetLabels`) drives
  visible labels from the active dataset's `ui_labels`, falling back to the current
  FGFR2 wording — so the FGFR2 UI is unchanged.

## What makes a gene runnable?

There are two bars, matching the two layers:

**Core gene analysis (generic):** possible for any supported protein-coding gene
with `analysis_modes.core_gene_analysis: true`, once gene models can be collected
and the core adapter / views are wired. This provides overview, gene models,
domain architecture, synteny and all-exon boundary distances — **no event region
needed**.

**Full event-specific analysis:** a config alone is **not** enough. The
event-specific layer is runnable only when **all** of the following hold:

1. **A valid `gene_config` exists** (loads + validates via `scripts/framework/gene_config.py`).
2. **A configured event region** (`analysis_modes.event_analysis: configured` /
   `event.status: configured`) with a **supported event detector** for its
   `analysis_id`, registered as `status: supported` in
   `configs/framework/event_detectors.yaml`.
3. **The detector produces the contract outputs** defined in
   `docs/event_detector_contract.md` (`event_isoform_candidates.tsv`,
   `event_region_coordinates.tsv`, `event_label_evidence.tsv`,
   `event_detector_report.json`).
4. **The adapter / index builder supports the event type**
   (`scripts/adapters/fgfr2_to_generic_indices.py` can build the generic indices
   from the detector/core outputs).
5. **UI views are enabled and tested** (config promoted to `runnable: true`).

`scripts/framework/gene_config.py:discover_analyses()` treats an analysis as
runnable only when it has an active config **and** a supported detector; the
feasibility probe reports the same.

### Why FGFR1 is not runnable yet

`drafts/FGFR1_draft.yaml` (and `drafts/TPM1_draft.yaml`) satisfy requirement 1
only: they are valid draft configs. In the core-vs-event model their **core**
gene analysis is *potentially supported* (`core_gene_analysis: true`), but their
**event** analysis is `not_configured`: no curated event region, markers, or
`supported` detector exist. The feasibility probe reports
`core_gene_analysis: potentially_supported` **and**
`event_analysis: not_configured` (overall `requires_gene_specific_event_detector`),
and they stay `runnable: false`.

So the honest statement is: FGFR1 *could* in principle produce core gene-level
results (gene models, domains, synteny, all-exon boundaries) once wired end-to-end,
but its **event-specific** IIIb/IIIc-style views require a gene-specific detector
first — not just a config edit. FGFR2 remains the only validated, wired analysis.

## Event detectors

- **Contract:** `docs/event_detector_contract.md` (+ machine-readable
  `configs/framework/event_detector_contract.yaml`) defines the gene-agnostic
  output files every detector must produce.
- **Registry:** `configs/framework/event_detectors.yaml` maps detectors to
  analyses and records their status. `fgfr2_iiib_iiic` is `supported`
  (`legacy_fgfr2_adapter`); generic detectors are `planned`.
- **FGFR2 detector (first supported):** `scripts/adapters/fgfr2_event_detector_adapter.py`
  *projects* the existing validated FGFR2 outputs into the contract shape — it
  recomputes nothing and changes no biology.

## Core analysis (generic, event-independent)

- **Contract:** `docs/core_gene_analysis_contract.md` (+ machine-readable
  `configs/framework/core_gene_analysis_contract.yaml`) defines the gene-agnostic
  core outputs (`gene_model_index.tsv`, `protein_isoform_index.tsv`,
  `exon_protein_map.tsv`, `domain_features.tsv`, `tm_features.tsv`,
  `synteny_neighbors.tsv`, `exon_domain_boundary_distances.tsv`,
  `core_gene_report.json`).
- **FGFR2 core adapter:** `scripts/adapters/fgfr2_core_analysis_adapter.py`
  projects existing FGFR2 outputs into the core contract, including the **generic
  all-exon exon–domain boundary** distances (`exactly_aligned` / `near_boundary` /
  `inside_domain` / `outside_annotated_domain`). Only boundary geometry is
  computed; no biology is recomputed.
- The generic all-exon boundary view is available whenever domain architecture +
  the exon map exist — independent of any event region. FGFR2 additionally keeps
  its event-specific "Boundary Consistency" (cassette-vs-domain) view.

## Core-only pilot & generic index builder

- **Core index builder:** `scripts/framework/build_core_gene_indices.py` reads
  the core contract TSVs (from a run's `results/core_gene_analysis/` or an explicit
  `--core-dir`) and writes the generic website indices **without** any event
  outputs. For a no-event config it emits `event_region_index.json` with
  `available: false, reason: "no_event_configured"`, disables the event-specific
  boundary view, and enables the generic Exon–Domain Boundaries view.
- **Mock no-event dataset (UI test):** `artifacts/generic_indices/mock_core_only/`
  (built from `artifacts/core_gene_analysis/mock/`) is a *synthetic* one-gene /
  one-species / one-protein / two-domain dataset with no event region. It exists
  **only** to smoke-test the UI's `has_event = false` path and is not a biological
  result.

## Live core-only runner (experimental)

There is now a real, gene-agnostic core runner:
`scripts/framework/run_core_gene_analysis.py`. It collects gene models/proteins
for a gene + species (input modes `auto` / `local_cache` / `user_files`), builds
the core contract + primary FASTA, extracts synteny, and — after the cluster
round-trip — parses **real** InterProScan/pyTMHMM outputs into domain/TM features
and all-exon boundary distances. See `docs/core_gene_analysis_contract.md`.

```bash
python scripts/framework/run_core_gene_analysis.py \
    --gene-config configs/genes/drafts/FGFR1_core_only_pilot.yaml \
    --species "Gallus gallus"
python scripts/framework/validate_core_gene_run.py --run-id <run_id>
python scripts/interpro_cluster/run_cluster_roundtrip.py --run-id <run_id>
```

What is still needed to move a core-only gene beyond *experimental*:

1. **Broader multi-species collection + auto-download** (the `auto` mode currently
   uses a local NCBI cache; live download is planned).
2. **Validation** on multiple species before exposing anything beyond experimental.
3. *(Only if event-specific views are ever wanted)* a gene-specific event detector
   + curated markers, or a user-defined event region — otherwise the gene stays
   core-only.

## Future analysis workflow

Adding a future gene is a staged process. Each stage must pass before the next:

1. **Create a gene config** (start as a draft in `drafts/`, `runnable: false`).
2. **Run the feasibility probe:**
   ```bash
   python scripts/framework/probe_gene_event_feasibility.py \
       --gene <SYMBOL> --reference-species homo_sapiens \
       --outdir artifacts/gene_feasibility/<SYMBOL>
   ```
   This reports `supported` / `partially_supported` /
   `requires_gene_specific_event_detector` / `not_supported_yet`, plus which
   modules are reusable vs gene-specific (from `configs/framework/module_capabilities.yaml`).
3. **Implement the event detector if needed.** The biological event detection
   (finding the event region + discriminating markers) is the real work and is
   currently FGFR2-specific.
4. **Add adapter mapping** for the new event type in
   `scripts/adapters/` so the generic indices can be produced.
5. **Validate** with the example dataset preview and at least one custom run.
6. **Enable in the UI only after tests pass** (promote from draft to an active,
   `runnable: true` config).

## What is generic now vs still FGFR2-specific

**Generic now (reusable for future genes):**

- The **config layer** (`configs/genes/*.yaml`, `scripts/framework/gene_config.py`).
- **Run management** (folders, status model, stop/delete/refresh) and
  **species normalization** — independent of the gene.
- **External annotation tools** (InterProScan, pyTMHMM) and the
  **domain-architecture** roll-up — driven by protein sequence, not by the gene.
- The **generic website-index schema** and the webapp shell, whose labels are
  **config-driven** (`getDatasetLabels`, backend `ui_labels`).
- **Feasibility probing** and the **module capability inventory**.

**Still FGFR2-specific (must be implemented per gene):**

- **Biological event detection** — the FGFR2 IIIb/IIIc marker detection
  (`SGINSSN` / `GVNTTDKEI`) and IIIb/IIIc label reconciliation are not reusable.
- **Event region definition** — each gene needs its own event region + detector.
- The **FGFR2→generic adapter** reads FGFR2 field names; a new event type needs
  its own adapter mapping.
- The concrete **human reference/control residues** are gene-specific (the
  reference/control *pattern* is reusable).

> **Bottom line:** the webapp / config / index layer is becoming gene-agnostic,
> but the science (event detection) is still FGFR2-specific. A new gene is not
> runnable until it has its own event definition/detector, adapter coverage, and
> passes validation. Configs and drafts alone do **not** make a gene work.
