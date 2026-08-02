# ExonDomainCompare — Architecture Atlas

A publication-quality set of architecture diagrams and a multi-page guide that
explain the ExonDomainCompare framework: one generic, annotation-aware
exon–protein analysis framework, with FGFR2 as a validated specialization built
on top of it. All artifacts are derived from the current codebase and are fully
reproducible from editable sources.

## Contents

| File | What it is |
| --- | --- |
| `01_overall_system_architecture.pdf` | The complete system at a glance (researcher → results). |
| `02_analysis_pipeline.pdf` | Ordered run stages: local, cluster and post-cluster phases. |
| `03_component_diagram.pdf` | Reusable software components by architectural layer. |
| `04_data_flow_diagram.pdf` | How biological data objects are transformed. |
| `05_c4_container_diagram.pdf` | C4 container view (scientific-computing adaptation). |
| `06_generic_vs_fgfr2.pdf` | Shared foundation + generic and FGFR2 modes. |
| `07_scientific_user_workflow.pdf` | Researcher journey and selection synchronization. |
| `08_post_interpro_pipeline.pdf` | Post-cluster annotation layering + boundary classification. |
| `09_exondomaincompare_architecture_poster.pdf` | Landscape overview poster (presentations / appendix). |
| `ExonDomainCompare_Architecture_Guide.pdf` | ~20-page guide with every diagram, captions and repository references. |
| `architecture_inventory.md` | Traceability: each diagram component → repository location. |

Editable sources and intermediate renders:

```
docs/architecture/
├── README.md
├── architecture_inventory.md
├── sources/            # D2 diagram sources (+ _theme.d2 shared identity)
├── svg/                # rendered SVGs (intermediate)
├── guide/              # Typst source for the Architecture Guide
├── 0*.pdf              # distribution diagram PDFs (vector)
└── ExonDomainCompare_Architecture_Guide.pdf
```

## Building

One-time toolchain (open-source, via Homebrew):

```bash
brew install d2 librsvg typst
```

- **d2** — compiles the diagram sources to SVG.
- **rsvg-convert** (librsvg) — converts SVG to vector PDF (renders D2's text faithfully).
- **typst** — typesets the Architecture Guide and embeds the vector diagram PDFs.

Build everything (diagrams → SVG + PDF, then the guide):

```bash
python scripts/docs/build_architecture_atlas.py
# or:
make architecture
```

Subsets and checks:

```bash
python scripts/docs/build_architecture_atlas.py --diagrams   # diagrams only
python scripts/docs/build_architecture_atlas.py --guide      # guide only
python scripts/docs/build_architecture_atlas.py --check      # verify toolchain
```

## Visual design system

Every diagram shares one semantic palette with stable meaning (defined in
`sources/_theme.d2`):

- **User / frontend** (indigo) · **Local compute** (teal) · **LRZ cluster** (orange)
- **Scientific data** (violet) · **Annotation / boundary** (cyan) · **Generic mode** (green)
- **Outputs** (amber) · **FGFR2 specialization** (rose accent)

Shapes: rounded rectangles = software components; documents = files/datasets;
cylinders = persistent stores; clouds = external systems; person = researcher.
Arrows: solid = primary data flow; dashed = optional/conditional; dotted =
configuration/reference; rose = FGFR2 specialization.

## Accuracy notes

The diagrams reflect the current implementation. In particular:

- Species are always a collection; single- and multi-species use the same pipeline.
- Only **primary proteins** are submitted to the cluster for InterProScan / pyTMHMM.
- Domain Architecture and Boundary views require completed cluster annotation
  (they show a calm _pending_ state beforehand).
- Only `FGFR2_IIIb_IIIc` is a validated event; FGFR1 / TPM1 / TP53 are core-only pilots.
- Generic paths never use IIIb/IIIc or cassette terminology.

See `architecture_inventory.md` for the component-to-code mapping.
