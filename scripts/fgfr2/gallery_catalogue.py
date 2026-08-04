from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from fgfr2 import coordinate_model as cm  # noqa: E402
from exondomaincompare.shared_gene_analysis import model_roles as mr  # noqa: E402
from exondomaincompare.shared_gene_analysis import species_order as so  # noqa: E402

SCHEMA_VERSION = 1
GENE_SYMBOL = "FGFR2"

CLOSURE_FIGURES = cm.CLOSURE / "figures"
ARCH_PER_SPECIES = cm.ARCHITECTURE / "figures" / "per_species"
_BOUNDARY_ROOT = cm.FREEZE / "16_final_thesis_analyses" \
    / "exon_domain_boundary_consistency"
BOUNDARY_FIGURES = _BOUNDARY_ROOT / "figures"
BOUNDARY_TABLES = _BOUNDARY_ROOT / "tables"


def _asset_reference(path: Path) -> str:
    resolved = Path(path).resolve()
    for parent in resolved.parents:
        if (parent / "run_config.json").is_file():
            try:
                return resolved.relative_to(parent).as_posix()
            except ValueError:
                break
    return cm._rel(resolved)


PENDING_CLUSTER = "pending_cluster"

PENDING_CLUSTER_REASON = (
    "InterProScan and pyTMHMM annotation has not been returned from the cluster yet. "
    "This figure is produced by the cluster round-trip.")


PENDING_ANNOTATION = "pending_cluster_annotation"


COMPARATIVE_CATEGORY_ORDER = [
    "Comparative exon structure",
    "FGFR2 cassette evidence",
    "Comparative sequence analysis",
    "Comparative domain architecture",
    "FGFR2 IIIb/IIIc Boundary Consistency",
    "Comparative exon–domain boundaries",
    "Comparative genomic context",
]

GENERIC_BOUNDARY_SCOPE_NOTE = (
    "This analysis evaluates supported internal coding-exon boundaries across the "
    "complete protein architecture. It is separate from the validated FGFR2 "
    "IIIb/IIIc cassette-boundary analysis.")

#: Reading order inside one species scope.
SPECIES_CATEGORY_ORDER = [
    "Exon structure",
    "Isoform and cassette analysis",
    "Domain architecture",
    "Exon–domain boundaries",
    "Genomic context",
    "FGFR2 event evidence",
]

SUPPLEMENTS = "Supplements"

POST_CLUSTER_CATEGORIES = frozenset({
    "Comparative exon structure",
    "Comparative domain architecture",
    "FGFR2 IIIb/IIIc Boundary Consistency",
    "Comparative exon–domain boundaries",
    "Exon structure",
    "Isoform and cassette analysis",
    "Domain architecture",
    "Exon–domain boundaries",
    "FGFR2 event evidence",
})


@dataclass(frozen=True)
class Dataset:
    dataset_id: str
    closure: Path
    architecture: Path
    boundary_root: Path
    #: The tables a card cites as its source data, by the keys the card definitions use.
    tables: Dict[str, str]
    cluster_ready: bool = True
    coordinate_model_path: Optional[Path] = None

    @property
    def closure_figures(self) -> Path:
        return self.closure / "figures"

    @property
    def arch_overview(self) -> Path:
        return self.architecture / "figures" / "overview"

    @property
    def arch_per_species(self) -> Path:
        return self.architecture / "figures" / "per_species"

    @property
    def boundary_figures(self) -> Path:
        return self.boundary_root / "figures"

    @property
    def boundary_tables(self) -> Path:
        return self.boundary_root / "tables"

    @property
    def website_indices(self) -> Path:
        if self.dataset_id.startswith("run:"):
            return self.closure.parent.parent / "website_indices"
        return self.closure / "website_indices"

    def pending(self, category: str) -> str:
        if self.cluster_ready or category not in POST_CLUSTER_CATEGORIES:
            return ""
        return PENDING_CLUSTER_REASON


def freeze_dataset() -> Dataset:
    return Dataset(
        dataset_id="example",
        closure=cm.CLOSURE,
        architecture=cm.ARCHITECTURE,
        boundary_root=_BOUNDARY_ROOT,
        tables={
            "truth": cm._rel(cm.TRUTH_TSV),
            "architecture": cm._rel(cm.FEATURES_TSV),
            "interpro": cm._rel(cm.INTERPRO_TSV),
            "msa": cm._rel(cm.FULL_MSA),
        },
        cluster_ready=True,
        coordinate_model_path=(ROOT / "results" / "derived" / "example"
                               / "website_indices" / "protein_coordinate_model.json"),
    )


def run_dataset(run_dir: Path, *, cluster_ready: Optional[bool] = None) -> Dataset:
    run_dir = Path(run_dir)
    results = run_dir / "results"
    closure = results / "13_final_pre_interpro_closure"
    architecture = results / "15_exon_domain_boundary_post_interpro"
    boundary_root = results / "16_final_thesis_analyses" / "exon_domain_boundary_consistency"
    if cluster_ready is None:
        cluster_ready = (architecture / "tables"
                         / "exon_domain_architecture_features.tsv").is_file()
    return Dataset(
        dataset_id=f"run:{run_dir.name}",
        closure=closure,
        architecture=architecture,
        boundary_root=boundary_root,
        tables={
            "truth": _asset_reference(closure / "final_pre_interpro_truth_table.tsv"),
            "architecture": _asset_reference(
                architecture / "tables" / "exon_domain_architecture_features.tsv"),
            "interpro": _asset_reference(
                architecture / "tables" / "interpro_domain_features_normalized.tsv"),
            "msa": _asset_reference(
                closure / "MSA" / "final_fgfr2_full_length_protein_msa.aln.faa"),
        },
        cluster_ready=bool(cluster_ready),
        coordinate_model_path=(run_dir / "website_indices"
                               / "protein_coordinate_model.json"),
    )


EXPORT_FORMATS = ("png", "svg", "pdf", "tsv")


def _formats(directory: Path, stem: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for ext in EXPORT_FORMATS:
        path = directory / f"{stem}.{ext}"
        if path.is_file():
            out[ext] = _asset_reference(path)
    return out


def _mode(stem: str, label: str, directory: Optional[Path], *, description: str = "",
          default: bool = False) -> Optional[Dict[str, Any]]:
    if directory is None:
        return None
    formats = _formats(directory, stem)
    if not formats:
        return None
    return {
        "mode_id": stem,
        "label": label,
        "description": description,
        "is_default": default,
        "formats": formats,
        "thumbnail": formats.get("png", ""),
    }


def _card(figure_id: str, *, title: str, category: str, scope: str,
          figure_type: str, question: str, interpretation: str,
          renderer: str, source_data: Sequence[str],
          modes: Sequence[Optional[Dict[str, Any]]] = (),
          kind: str = "main", species_id: str = "",
          model_selection: Optional[Dict[str, Any]] = None,
          supersedes: Sequence[str] = (),
          availability: Optional[Dict[str, Any]] = None) -> Optional[Dict[str, Any]]:
    present = [m for m in modes if m]
    stated_absence = (availability or {}).get("unavailable_models")
    pending = PENDING_CLUSTER_REASON if category in POST_CLUSTER_CATEGORIES else ""
    if not present and not stated_absence and not pending:
        return None
    if present and not any(m["is_default"] for m in present):
        present[0]["is_default"] = True
    default = next((m for m in present if m["is_default"]), None)
    return {
        "figure_id": figure_id,
        "title": title,
        "category": category,
        "kind": kind,
        "scope": scope,
        "species_id": species_id,
        "figure_type": figure_type,
        "scientific_question": question,
        "interpretation": interpretation,
        "renderer": renderer,
        "source_data": list(source_data),
        "model_selection": model_selection,
        "modes": present,
        "n_modes": len(present),
        "export_formats": sorted({f for m in present for f in m["formats"]}),
        "formats": dict(default["formats"]) if default else {},
        "thumbnail": default["thumbnail"] if default else "",
        "supersedes": list(supersedes),
        # Preserve explicit unsupported verdicts; otherwise an unrendered card is pending.
        "availability": ({**(availability or {}), "status": PENDING_CLUSTER,
                          "reason": pending}
                         if pending and not present and not stated_absence else
                         availability or {"status": "available"}),
    }


FREEZE_TABLES = freeze_dataset().tables

_FREEZE_RENDERER = "scripts/make_all_figures.py (validated FGFR2 freeze)"
_SHARED_RENDERER = "scripts/plotting/render_comparative_gallery_figures.mjs"


def _comparative_cards(comparative_dir: Optional[Path],
                       ds: Dataset) -> List[Dict[str, Any]]:
    F, A = ds.closure_figures, ds.arch_overview
    FREEZE_TABLES = ds.tables
    BOUNDARY_TABLES, BOUNDARY_FIGURES = ds.boundary_tables, ds.boundary_figures
    cards: List[Optional[Dict[str, Any]]] = [

        _card(
            "fgfr2_cmp_all_species_exon_domain_architecture",
            title="All species · exon–domain architecture",
            category="Comparative exon structure",
            scope="comparative",
            figure_type="all_species_exon_domain_architecture",
            question="How do coding-exon structure, domain architecture and the "
                     "IIIb/IIIc cassette slot line up across the species panel?",
            interpretation="Rows are one protein per species and isoform, ordered "
                           "taxonomically. Ig1/Ig2/Ig3 are numbered by position, not "
                           "by an FGFR-specific signature, and the transmembrane "
                           "helix is a pyTMHMM prediction.",
            renderer=_FREEZE_RENDERER,
            source_data=[FREEZE_TABLES["architecture"]],
            modes=[
                _mode("Figure_10_all_species_FGFR2_exon_domain_architecture_primary",
                      "All species", A, default=True,
                      description="Every protein in the panel."),
                _mode("Figure_10A_IIIb_exon_domain_architecture_primary",
                      "IIIb only", A, description="Cassette filter: IIIb proteins."),
                _mode("Figure_10B_IIIc_exon_domain_architecture_primary",
                      "IIIc only", A, description="Cassette filter: IIIc proteins."),
                _mode("Figure_10C_mammals_exon_domain_architecture_primary",
                      "Mammals", A, description="Taxon filter: mammals."),
                _mode("Figure_10D_nonmammals_exon_domain_architecture_primary",
                      "Non-mammals", A, description="Taxon filter: non-mammals."),
            ],
            supersedes=["Figure_10A_IIIb_exon_domain_architecture_primary",
                        "Figure_10B_IIIc_exon_domain_architecture_primary",
                        "Figure_10C_mammals_exon_domain_architecture_primary",
                        "Figure_10D_nonmammals_exon_domain_architecture_primary"],
        ),

        _card(
            "fgfr2_cmp_cassette_zoom",
            title="IIIb/IIIc cassette · zoom and exon-to-protein coordinates",
            category="FGFR2 cassette evidence",
            scope="comparative",
            figure_type="cassette_zoom",
            question="Where does the mutually exclusive IIIb/IIIc cassette sit in "
                     "each protein, and how does its coding exon project onto the "
                     "protein sequence?",
            interpretation="The cassette interval is the validated projection from "
                           "the freeze. Codon phase at the cassette boundary is "
                           "unknown for part of the panel, which the per-protein "
                           "flags record.",
            renderer=_FREEZE_RENDERER,
            source_data=[FREEZE_TABLES["architecture"], FREEZE_TABLES["truth"]],
            modes=[
                _mode("Figure_3_final_IIIb_IIIc_cassette_zoom_pre_interpro",
                      "Cassette zoom", F, default=True,
                      description="The cassette region across the panel."),
                _mode("Figure_3C_exon_to_protein_cassette_coordinate_map",
                      "Exon-to-protein coordinate map", F,
                      description="The same cassette as a coordinate mapping."),
            ],
            supersedes=["Figure_3C_exon_to_protein_cassette_coordinate_map"],
        ),

        _card(
            "fgfr2_cmp_cassette_motif_map",
            title="IIIb/IIIc cassette · amino-acid motif map",
            category="FGFR2 cassette evidence",
            scope="comparative",
            figure_type="cassette_motif_map",
            question="Which amino-acid positions inside the cassette distinguish "
                     "the IIIb form from the IIIc form across species?",
            interpretation="Positions are read from the cassette alignment. A "
                           "conserved position is a conserved alignment column, "
                           "which is not by itself evidence of a functional role.",
            renderer=_FREEZE_RENDERER,
            source_data=[FREEZE_TABLES["msa"]],
            modes=[_mode("Figure_3B_IIIb_IIIc_cassette_amino_acid_motif_map",
                         "Motif map", F, default=True)],
        ),

        _card(
            "fgfr2_cmp_msa_integrity",
            title="Full-length protein alignment · integrity",
            category="Comparative sequence analysis",
            scope="comparative",
            figure_type="msa_integrity",
            question="Is the full-length FGFR2 protein alignment sound enough to "
                     "carry cross-species comparisons?",
            interpretation="Gap fraction and outlier status per sequence. A "
                           "sequence passing this check is comparable in alignment "
                           "columns; it says nothing about its biological quality.",
            renderer=_FREEZE_RENDERER,
            source_data=[FREEZE_TABLES["msa"], FREEZE_TABLES["truth"]],
            modes=[_mode("Figure_5_full_length_FGFR2_MSA_integrity_paper",
                         "Alignment integrity", F, default=True)],
        ),

        _card(
            "fgfr2_cmp_residue_agreement",
            title="IIIb/IIIc residue agreement",
            category="Comparative sequence analysis",
            scope="comparative",
            figure_type="residue_agreement",
            question="Do the cassette residues of each species agree with the "
                     "human reference assignment of IIIb and IIIc?",
            interpretation="Agreement is measured against the human reference "
                           "residues. A disagreement is a reason to inspect a "
                           "sequence, not a conclusion about the isoform.",
            renderer=_FREEZE_RENDERER,
            source_data=[FREEZE_TABLES["msa"], FREEZE_TABLES["truth"]],
            modes=[
                _mode("Figure_6_human_referenced_IIIb_IIIc_residue_agreement",
                      "Human-referenced", F, default=True,
                      description="Agreement against the human reference residues."),
                _mode("Figure_6B_species_resolved_IIIb_IIIc_cassette_residue_map",
                      "Species-resolved", F,
                      description="The same residues, resolved per species."),
            ],
            supersedes=["Figure_6B_species_resolved_IIIb_IIIc_cassette_residue_map"],
        ),

        _card(
            "fgfr2_cmp_isoform_discriminating_residues",
            title="Isoform-discriminating residues",
            category="Comparative sequence analysis",
            scope="comparative",
            figure_type="isoform_discriminating_residues",
            question="Which residues separate the IIIb and IIIc forms consistently "
                     "across the panel?",
            interpretation="Discriminating positions are those that differ between "
                           "the two forms in the aligned panel. Consistency across "
                           "species is positional evidence only.",
            renderer=_FREEZE_RENDERER,
            source_data=[FREEZE_TABLES["msa"]],
            modes=[_mode("Figure_7_isoform_discriminating_residues",
                         "Discriminating residues", F, default=True)],
        ),

        _card(
            "fgfr2_cmp_boundary_consistency_matrix",
            title="IIIb/IIIc Boundary Consistency Matrix",
            category="FGFR2 IIIb/IIIc Boundary Consistency",
            scope="comparative",
            figure_type="boundary_consistency_matrix",
            question="How does each species' cassette start and end sit relative to "
                     "the nearest protein-domain boundary?",
            interpretation="Classes are the frozen FGFR2 Boundary Consistency "
                           "vocabulary, which is separate from the generic "
                           "exon-boundary classes used in the species scopes. "
                           "Co-location is a positional observation.",
            renderer=_FREEZE_RENDERER,
            source_data=[_asset_reference(
                BOUNDARY_TABLES / "exon_domain_boundary_distances.tsv")],
            modes=[_mode("Figure_11_exon_domain_boundary_consistency_heatmap",
                         "Consistency matrix", BOUNDARY_FIGURES, default=True)],
        ),

        _card(
            "fgfr2_cmp_boundary_distance_distribution",
            title="IIIb/IIIc Boundary-distance distribution",
            category="FGFR2 IIIb/IIIc Boundary Consistency",
            scope="comparative",
            figure_type="boundary_distance_distribution",
            question="How far are cassette boundaries from the nearest protein-domain "
                     "boundary, across the panel?",
            interpretation="The distribution of absolute distances by isoform. The "
                           "shape depends on the near-edge threshold stated in the "
                           "figure.",
            renderer=_FREEZE_RENDERER,
            source_data=[_asset_reference(
                BOUNDARY_TABLES / "exon_domain_boundary_distances.tsv")],
            modes=[_mode("Figure_12_boundary_distance_distribution",
                         "Distance distribution", BOUNDARY_FIGURES, default=True)],
        ),

        _card(
            "fgfr2_cmp_synteny_neighbourhood",
            title="Local synteny neighbourhood",
            category="Comparative genomic context",
            scope="comparative",
            figure_type="synteny_neighbourhood",
            question="Is the FGFR2 locus in a conserved gene neighbourhood across "
                     "the panel?",
            interpretation="Neighbour loci as annotated in each assembly. A missing "
                           "neighbour is an annotation gap in that assembly, not "
                           "evidence that the gene is absent.",
            renderer=_FREEZE_RENDERER,
            source_data=[_asset_reference(ds.website_indices
                                          / "synteny_locus_index.json")],
            modes=[
                _mode("Figure_9A_FGFR2_local_synteny_5neighbor_paper",
                      "5 neighbours", F, default=True,
                      description="Five loci either side of the target."),
                _mode("Supplement_Figure_FGFR2_local_synteny_10neighbor_all_species",
                      "10 neighbours", F,
                      description="Ten loci either side, all species."),
            ],
            # Byte-identical to 9A; it was a second entry point to one figure.
            supersedes=["Figure_9_FGFR2_local_synteny_neighborhood"],
        ),

    ]

    G = "Comparative exon–domain boundaries"
    boundary_source = ([_asset_reference(ds.coordinate_model_path)]
                       if ds.coordinate_model_path else [])
    cards += [
        _card(
            "fgfr2_cmp_all_exon_boundary_matrix",
            title="All coding-exon Boundary matrix",
            category=G,
            scope="comparative",
            figure_type="cmp_boundary_matrix",
            question="For each comparable coding-exon boundary, how does every "
                     "species place it relative to the nearest representative "
                     "domain edge?",
            interpretation=GENERIC_BOUNDARY_SCOPE_NOTE + " One cell per species "
                           "and comparable-boundary group, showing that species' "
                           "own class: exact, near, inside, outside, uncertain or "
                           "unmapped. An empty cell means the species has no "
                           "supported observation of that group, not agreement.",
            renderer=_SHARED_RENDERER,
            source_data=boundary_source,
            modes=[_mode("cmp_boundary_matrix", "Boundary matrix",
                         comparative_dir, default=True)],
        ),
        _card(
            "fgfr2_cmp_all_exon_signed_distance",
            title="All coding-exon signed-distance comparison",
            category=G,
            scope="comparative",
            figure_type="cmp_paired_signed_distance",
            question="Do the species that share a coding-exon boundary place it "
                     "at the same signed distance from the same domain edge?",
            interpretation=GENERIC_BOUNDARY_SCOPE_NOTE + " Distances are signed "
                           "against the specific domain instance each species "
                           "measured, so direction and domain identity are both "
                           "preserved. Agreement across species is a positional "
                           "observation and is not evidence of conservation.",
            renderer=_SHARED_RENDERER,
            source_data=boundary_source,
            modes=[_mode("cmp_paired_signed_distance", "Signed distances",
                         comparative_dir, default=True)],
        ),
        _card(
            "fgfr2_cmp_all_exon_boundary_consistency",
            title="All coding-exon Boundary-position consistency",
            category=G,
            scope="comparative",
            figure_type="cmp_boundary_position_consistency",
            question="Which comparable coding-exon boundaries sit in a similar "
                     "position across species, and which vary?",
            interpretation=GENERIC_BOUNDARY_SCOPE_NOTE + " Per group: how many "
                           "species observed it, how the observations spread, and "
                           "how the group was mapped. A group mapped in every "
                           "species can still be tentative, so the mapping status "
                           "is shown rather than folded into a score.",
            renderer=_SHARED_RENDERER,
            source_data=boundary_source,
            modes=[_mode("cmp_boundary_position_consistency",
                         "Position consistency", comparative_dir, default=True)],
        ),
    ]

    cards.append(_card(
        "fgfr2_cmp_shared_exon_domain_architecture",
        title="Comparative exon–domain architecture (shared renderer)",
        category="Comparative domain architecture",
        scope="comparative",
        figure_type="cmp_exon_domain_architecture",
        question="Do exon boundaries fall at the same domain edges in every "
                 "species, in native coordinates and in alignment columns?",
        interpretation="Drawn by the same renderer every other gene uses, from "
                       "the shared coordinate model, with one row per species' "
                       "primary reference model. Native and alignment coordinates "
                       "are two views of the same boundaries.",
        renderer=_SHARED_RENDERER,
        source_data=boundary_source,
        modes=[
            _mode("cmp_exon_domain_architecture_native", "Native coordinates",
                  comparative_dir, default=True,
                  description="Amino-acid positions in each protein."),
            _mode("cmp_exon_domain_architecture_msa", "Alignment columns",
                  comparative_dir,
                  description="The same features in shared alignment columns."),
        ],
    ))

    return [c for c in cards if c]



WITHDRAWN_CARDS = [
    {
        "figure_id": "fgfr2_cmp_framework_evidence_stack",
        "title": "Framework evidence stack",
        "source_figures": ["Figure_8_final_framework_evidence_stack",
                           "Figure_Final_Framework_Evidence_Stack"],
        "retained_source_data": ["truth"],
        "reason": ("A per-protein QC view of which evidence layer set each final "
                   "isoform label. It documents how the analysis reached its "
                   "conclusions rather than answering a biological question, so it "
                   "belongs to the validation record and not to the figure "
                   "catalogue."),
    },
    {
        "figure_id": "fgfr2_cmp_synteny_conservation_matrix",
        "title": "Synteny neighbour conservation",
        "source_figures": ["Figure_9B_FGFR2_5neighbor_conservation_matrix_paper"],
        "retained_source_data": ["synteny_locus_index"],
        "reason": ("Per-neighbour conservation cells over the same loci the main "
                   "comparative Synteny neighbourhood figure already shows, where a "
                   "blank cell is an annotation gap in one assembly rather than "
                   "absence of the gene. The neighbourhood figure remains the "
                   "catalogue's genomic-context card."),
    },
]


_SUPPLEMENT_DEFS = [
    ("Figure_1_framework_overview", "Framework overview",
     "Which analysis steps produced the final isoform assignments?",
     "A methods diagram of the pipeline, not a result."),
    ("Figure_2_final_exon_to_protein_architecture_pre_interpro",
     "Exon-to-protein architecture (pre-InterPro)",
     "How did coding exons project onto the proteins before the domain layer "
     "existed?",
     "Superseded by the all-species exon–domain architecture, which adds the "
     "domain layer to the same projection. Kept because the thesis cites it."),
    ("Figure_4_label_reconciliation_and_rescue_summary",
     "Label reconciliation and rescue summary",
     "How many proteins needed their upstream label corrected or rescued?",
     "Counts of reconciliation outcomes. A rescued label is one supported by an "
     "external validated candidate, recorded per protein in the truth table."),
    ("Figure_9C_FGFR2_synteny_review_cases_paper", "Synteny review cases",
     "Which loci did the synteny analysis flag for manual review?",
     "The flagged cases only. A flag marks a case worth inspecting, not an error."),
    ("Supplement_all_species_cassette_zoom", "All species · cassette zoom",
     "How does the cassette region look in every protein of the panel?",
     "The per-protein form of the cassette zoom."),
    ("Supplement_all_species_exon_protein_architecture",
     "All species · exon-to-protein architecture",
     "How do coding exons project onto every protein of the panel?",
     "The per-protein form of the pre-InterPro projection."),
    ("Supplement_full_length_MSA_QC_histograms",
     "Full-length alignment · QC distributions",
     "How are gap fraction and outlier scores distributed across the alignment?",
     "Distributions behind the alignment integrity figure."),
    ("Supplement_Figure_6_residue_agreement_review_cases",
     "Residue agreement · review cases",
     "Which sequences disagreed with the human reference residues?",
     "The disagreeing cases only, for inspection."),
    ("Supplement_review_cases_pre_interpro", "Review cases",
     "Which proteins carried a review flag before the InterPro step?",
     "Flagged proteins with the flag that was raised."),
    ("Supplement_Figure_6B_review_rows_cassette_residue_map",
     "Review rows · cassette residue map",
     "What do the cassette residues look like for the flagged rows?",
     "The residue map restricted to flagged rows."),
    ("Supplement_review_unresolved_case_panels", "Unresolved case panels",
     "Which cases stayed unresolved, and what was seen for each?",
     "The unresolved cases with their recorded reason. An unresolved case is "
     "reported as unresolved and is not counted as a validated claim."),
]

#: Byte-identical to another supplement in the freeze; listing both made one
#: figure look like two independent pieces of evidence.
_SUPPLEMENT_DUPLICATES = {
    "Supplement_full_length_MSA_outliers": "Supplement_full_length_MSA_QC_histograms",
}


def _supplement_cards(ds: Dataset) -> List[Dict[str, Any]]:
    FREEZE_TABLES, CLOSURE_FIGURES = ds.tables, ds.closure_figures
    out = []
    for stem, title, question, interpretation in _SUPPLEMENT_DEFS:
        superseded = [k for k, v in _SUPPLEMENT_DUPLICATES.items() if v == stem]
        card = _card(
            f"fgfr2_supp_{stem}",
            title=title,
            category=SUPPLEMENTS,
            scope="comparative",
            figure_type=stem,
            question=question,
            interpretation=interpretation,
            renderer=_FREEZE_RENDERER,
            source_data=[FREEZE_TABLES["truth"]],
            modes=[_mode(stem, title, CLOSURE_FIGURES, default=True)],
            kind="supplement",
            supersedes=superseded,
        )
        if card:
            out.append(card)
    return out


_SPECIES_CARD_DEFS = [
    ("primary_exon_projection", "Coding exon projection", "Exon structure",
     "Which coding exons produce which regions of this protein?",
     "Coding exons projected onto the protein sequence. The cassette exon is "
     "marked as the validated event it is, not as an exploratory candidate."),
    ("integrated_domain_architecture", "Integrated domain architecture",
     "Domain architecture",
     "How do the Ig-like domains, the kinase domain, membrane topology and the "
     "coding exons align along this protein?",
     "Representative domain instances are resolved by coordinate, so the three "
     "Ig-like domains stay distinct. Ig numbering is positional; the "
     "transmembrane helix is a pyTMHMM prediction."),
    ("boundary_on_architecture", "Exon boundaries on the architecture",
     "Exon–domain boundaries",
     "Where do this protein's internal coding-exon boundaries fall relative to "
     "its domain edges?",
     "Generic exon-boundary classes, not the frozen FGFR2 Boundary Consistency "
     "vocabulary. Co-location of a boundary with a domain edge is positional."),
    ("signed_boundary_distances", "Signed boundary distances",
     "Exon–domain boundaries",
     "How far, and on which side, does each boundary sit from the nearest domain "
     "edge?",
     "Distances are signed against the specific domain instance measured, so "
     "direction and domain identity are both preserved."),
    ("boundary_class_summary", "Boundary-class summary", "Exon–domain boundaries",
     "How are this protein's exon boundaries distributed across the "
     "domain-relation classes?",
     "Counts of mutually exclusive classes for one protein. Class membership "
     "depends on the near-edge threshold stated in the figure."),
]

_CARD_REQUIRES = {
    "primary_exon_projection": "exon_structure",
    "integrated_domain_architecture": "domain_architecture",
    "boundary_on_architecture": "exon_domain_boundaries",
    "signed_boundary_distances": "exon_domain_boundaries",
    "boundary_class_summary": "exon_domain_boundaries",
}


def _model_supports(model: Dict[str, Any], layer: str) -> bool:
    if layer in (model.get("unavailable_layers") or []):
        return False
    if layer == "domain_architecture":
        return bool(model.get("representative_domains"))
    if layer == "exon_structure":
        return any(not e.get("is_cassette_exon") for e in model.get("exons") or [])
    if layer == "exon_domain_boundaries":
        return any(b.get("signed_distance") is not None
                   for b in model.get("exon_boundaries") or [])
    return True


def _model_stem(model: Dict[str, Any], figure_type: str) -> str:
    species = model["species_id"]
    key = species if model.get("is_primary_reference") else \
        f"{species}_{model.get('final_isoform_label') or model.get('isoform')}"
    return f"main_{key}_{figure_type}"


def _species_scope(species_id: str, models: Sequence[Dict[str, Any]],
                   main_dir: Optional[Path],
                   availability: Dict[str, Any],
                   ds: Dataset) -> Dict[str, Any]:
    FREEZE_TABLES, ARCH_PER_SPECIES = ds.tables, ds.arch_per_species
    per_species = (availability.get("per_species") or {}).get(species_id) or {}
    unavailable = [u for u in availability.get("unavailable_combinations") or []
                   if u["species_id"] == species_id]

    model_entries = [{
        "model_id": m["model_id"],
        "model_role": m["model_role"],
        "role_label": mr.role_label(m["model_role"]),
        "isoform": m.get("final_isoform_label") or m.get("isoform"),
        "protein_id": m["protein_id"],
        "transcript_id": m.get("transcript_id") or "",
        "protein_length": m.get("protein_length"),
        "is_primary_reference": m["is_primary_reference"],
        "availability_status": m["availability_status"],
        "unavailable_layers": m.get("unavailable_layers") or [],
        "unavailable_reason": m.get("unavailable_reason") or "",
        "reconstruction_status": m.get("reconstruction_status") or "",
        "review_status": m.get("review_status") or "",
    } for m in models]

    for u in unavailable:
        model_entries.append({
            "model_id": "",
            "model_role": "",
            "role_label": f"{u['isoform']} model unavailable",
            "isoform": u["isoform"],
            "protein_id": u.get("protein_id") or "",
            "transcript_id": "",
            "protein_length": None,
            "is_primary_reference": False,
            "availability_status": u["availability_status"],
            "unavailable_layers": ["exon_structure", "domain_architecture",
                                   "exon_domain_boundaries"],
            "unavailable_reason": u["omission_reason"],
            "reconstruction_status": "",
            "review_status": u.get("review_status") or "",
        })

    all_pending = bool(models) and not unavailable and all(
        m.get("availability_status") == PENDING_ANNOTATION for m in models)

    cards: List[Dict[str, Any]] = []
    for figure_type, title, category, question, interpretation in _SPECIES_CARD_DEFS:
        layer = _CARD_REQUIRES[figure_type]
        supporting = [m for m in models if _model_supports(m, layer)]
        modes = []
        for m in supporting:
            label = m.get("final_isoform_label") or m.get("isoform") or "primary"
            modes.append(_mode(
                _model_stem(m, figure_type), label, main_dir,
                description=mr.role_label(m["model_role"]),
                default=m["is_primary_reference"]) if main_dir else None)
        withheld = [{
            "isoform": m.get("final_isoform_label") or m.get("isoform"),
            "model_id": m["model_id"],
            "reason": m.get("unavailable_reason")
                      or f"this model has no {layer.replace('_', ' ')} layer",
        } for m in models if m not in supporting]
        withheld += [{
            "isoform": u["isoform"], "model_id": "", "reason": u["omission_reason"],
        } for u in unavailable]
        card = _card(
            f"fgfr2_{species_id}_{figure_type}",
            title=title,
            category=category,
            scope="species",
            species_id=species_id,
            figure_type=figure_type,
            question=question,
            interpretation=interpretation,
            renderer="scripts/plotting/render_main_figures.mjs",
            source_data=[FREEZE_TABLES["architecture"]],
            modes=modes,
            model_selection={
                "kind": "isoform_model",
                "options": [{
                    "model_id": m["model_id"],
                    "label": m.get("final_isoform_label") or m.get("isoform"),
                    "model_role": m["model_role"],
                    "is_default": m["is_primary_reference"],
                } for m in supporting],
                "unavailable": withheld,
            },
            availability=({"status": "available", "unavailable_models": withheld}
                          if modes and any(modes) else
                          {"status": PENDING_CLUSTER, "reason": PENDING_CLUSTER_REASON,
                           "unavailable_models": withheld} if all_pending else
                          {"status": "unavailable", "unavailable_models": withheld}),
        )
        if card:
            cards.append(card)

    combined = [m for m in models if _model_supports(m, "exon_structure")]
    if len(combined) > 1 and main_dir:
        modes = [_mode(_model_stem(m, "primary_exon_projection"),
                       f"{m.get('final_isoform_label') or m.get('isoform')} exon series",
                       main_dir, default=m["is_primary_reference"],
                       description=mr.role_label(m["model_role"]))
                 for m in combined]
        card = _card(
            f"fgfr2_{species_id}_isoform_model_comparison",
            title="IIIb and IIIc · isoform model comparison",
            category="Isoform and cassette analysis",
            scope="species",
            species_id=species_id,
            figure_type="isoform_model_comparison",
            question="How do this species' IIIb and IIIc proteins differ in their "
                     "coding-exon series and in the position of the cassette?",
            interpretation="Two separate validated proteins, shown as two models. "
                           "Their coordinates are not merged: the cassette sits at a "
                           "different exon index in each transcript, so a single "
                           "combined numbering would describe neither protein.",
            renderer="scripts/plotting/render_main_figures.mjs",
            source_data=[FREEZE_TABLES["architecture"], FREEZE_TABLES["truth"]],
            modes=modes,
            model_selection={
                "kind": "isoform_model",
                "options": [{
                    "model_id": m["model_id"],
                    "label": m.get("final_isoform_label") or m.get("isoform"),
                    "model_role": m["model_role"],
                    "is_default": m["is_primary_reference"],
                } for m in combined],
                "unavailable": [],
            },
        )
        if card:
            cards.append(card)

    validated_modes = [
        _mode(f"{species_id}_{m.get('final_isoform_label') or m.get('isoform')}"
              f"_exon_domain_architecture",
              m.get("final_isoform_label") or m.get("isoform"), ARCH_PER_SPECIES,
              description=mr.role_label(m["model_role"]),
              default=m["is_primary_reference"])
        for m in models]
    validated = _card(
        f"fgfr2_{species_id}_validated_architecture",
        title="Validated exon–domain architecture (freeze figure)",
        category="FGFR2 event evidence",
        scope="species",
        species_id=species_id,
        figure_type="validated_exon_domain_architecture",
        question="What did the validated FGFR2 analysis publish for this species' "
                 "proteins?",
        interpretation="The figure as the freeze rendered it, unchanged. It is the "
                       "citable form; the interactive views above are drawn by the "
                       "shared renderer from the same validated coordinates.",
        renderer=_FREEZE_RENDERER,
        source_data=[FREEZE_TABLES["architecture"]],
        modes=validated_modes,
        model_selection={
            "kind": "isoform_model",
            "options": [{
                "model_id": m["model_id"],
                "label": m.get("final_isoform_label") or m.get("isoform"),
                "model_role": m["model_role"],
                "is_default": m["is_primary_reference"],
            } for m in models],
            "unavailable": [{
                "isoform": u["isoform"], "model_id": "",
                "reason": u["omission_reason"],
            } for u in unavailable],
        },
    )
    if validated:
        cards.append(validated)

    return {
        "species_id": species_id,
        "scientific_name": so.scientific_name(species_id),
        "taxon_group": so.taxon_group(species_id),
        "clade": so.clade_of(species_id),
        "models": model_entries,
        "n_models": len(models),
        "isoforms_available": per_species.get("isoforms_available") or [],
        "isoforms_unavailable": per_species.get("isoforms_unavailable") or [],
        "cards": cards,
        "n_cards": len(cards),
        "categories": [c for c in SPECIES_CATEGORY_ORDER
                       if any(card["category"] == c for card in cards)],
    }


CROSS_SPECIES_BOUNDARY_CATEGORY = "Comparative exon–domain boundaries"

NO_COMPARABLE_BOUNDARIES_REASON = (
    "No exon–domain boundary could be matched across species: no two species share an "
    "exon group or an aligned protein position at a boundary. There is nothing to "
    "compare, so this figure has no data rather than a missing file.")

SINGLE_SPECIES_REASON = (
    "This run analyses one species. A cross-species comparison needs at least two, so "
    "this figure does not apply to this run rather than being absent from it.")


def _refine_pending(cards: Sequence[Dict[str, Any]], ds: Dataset,
                    model_index: Optional[Dict[str, Any]] = None) -> None:
    if not ds.cluster_ready:
        return
    comparable = (((model_index or {}).get("boundary_dashboard") or {})
                  .get("multi_species") or {}).get("comparable_boundary_groups")
    no_comparable = model_index is not None and not comparable
    one_species = len({m.get("species_id") for m in
                       ((model_index or {}).get("models") or [])}) < 2
    for card in cards:
        if card["availability"].get("status") != PENDING_CLUSTER:
            continue
        comparative = str(card.get("category") or "").startswith("Comparative")
        if one_species and comparative:
            card["availability"] = {"status": "not_applicable",
                                    "reason": SINGLE_SPECIES_REASON}
            continue
        if no_comparable and card.get("category") == CROSS_SPECIES_BOUNDARY_CATEGORY:
            card["availability"] = {"status": "scientifically_unavailable",
                                    "reason": NO_COMPARABLE_BOUNDARIES_REASON}
            continue
        card["availability"] = {
            "status": "technically_missing",
            "reason": ("Cluster annotation was returned but this figure was not derived "
                       "from it. Rebuild the run's figures."),
        }


def build_catalogue(model_index: Dict[str, Any], *,
                    main_dir: Optional[Path] = None,
                    comparative_dir: Optional[Path] = None,
                    dataset: Optional[Dataset] = None) -> Dict[str, Any]:
    ds = dataset or freeze_dataset()
    availability = model_index.get("availability") or {}
    species_ids = model_index.get("species_scope") or []
    by_species = cm.models_by_species(model_index)

    comparative = _comparative_cards(comparative_dir, ds)
    supplements = _supplement_cards(ds)
    scopes = {sid: _species_scope(sid, by_species.get(sid) or [], main_dir,
                                  availability, ds)
              for sid in species_ids}
    _refine_pending(comparative + supplements
                    + [c for s in scopes.values() for c in s["cards"]], ds, model_index)

    retired = sorted({s for c in comparative + supplements for s in c["supersedes"]})
    species_card_counts = sorted({v["n_cards"] for v in scopes.values()})

    return {
        "schema_version": SCHEMA_VERSION,
        "gene_symbol": GENE_SYMBOL,
        "dataset": ds.dataset_id,
        "cluster_ready": ds.cluster_ready,
        "multi_species": len(species_ids) > 1,
        # A multi-species dataset opens on the comparative view: the panel is the
        # result, and one of thirty species is not a sensible starting point.
        "default_scope": "comparative" if len(species_ids) > 1 else (
            species_ids[0] if species_ids else "comparative"),
        "comparative_cards": comparative,
        "species_scopes": scopes,
        "supplements": supplements,
        "filters": {
            "species": [{
                "species_id": sid,
                "scientific_name": so.scientific_name(sid),
                "taxon_group": so.taxon_group(sid),
                "clade": so.clade_of(sid),
                "n_models": len(by_species.get(sid) or []),
            } for sid in species_ids],
            "taxon_groups": list(dict.fromkeys(
                so.taxon_group(sid) for sid in species_ids)),
            "comparative_categories": [c for c in COMPARATIVE_CATEGORY_ORDER
                                       if any(x["category"] == c for x in comparative)],
            "species_categories": SPECIES_CATEGORY_ORDER,
            "supplements_hidden_by_default": True,
            "isoform_models": list(cm.VALIDATED_ISOFORMS),
        },
        "availability": availability,
        "counts": {
            "n_species": len(species_ids),
            "n_comparative_cards": len(comparative),
            "n_comparative_main": sum(1 for c in comparative if c["kind"] == "main"),
            "n_supplement_cards": len(supplements),
            "n_species_scopes": len(scopes),
            "species_scope_card_counts": species_card_counts,
            "n_species_scope_cards_total": sum(v["n_cards"] for v in scopes.values()),
            "n_retired_or_merged_source_figures": len(retired),
        },
        "retired_source_figures": retired,
        "withdrawn_cards": WITHDRAWN_CARDS,
        "retained_source_data": {
            key: ds.tables[key] for card in WITHDRAWN_CARDS
            for key in card["retained_source_data"] if key in ds.tables
        },
        "provenance": {
            "generated_by": "scripts/fgfr2/gallery_catalogue.py",
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "freeze_is_read_only": True,
            "note": ("Rendered files are derived assets; a card is a scientific "
                     "question. One card can own several rendered files as modes and "
                     "export formats, and a file with no question behind it is not a "
                     "card."),
        },
    }


def withdrawn_figure_stems() -> set:
    return {stem for card in WITHDRAWN_CARDS for stem in card["source_figures"]}


def flatten_for_gallery(catalogue: Dict[str, Any]) -> Dict[str, Any]:
    figures: List[Dict[str, Any]] = []
    figures += catalogue["comparative_cards"]
    for scope in catalogue["species_scopes"].values():
        figures += scope["cards"]
    figures += catalogue["supplements"]

    # Expose card status directly so the frontend needs no FGFR2-specific rules.
    scientific = {s["species_id"]: s["scientific_name"]
                  for s in catalogue["filters"]["species"]}
    for card in figures:
        card["status"] = card["availability"]["status"]
        card["species"] = scientific.get(card["species_id"], "") \
            if card["species_id"] else "Comparative"
        card["source_table"] = (card["source_data"] or [""])[0]
        card["caption"] = card["interpretation"]
        card["section"] = card["category"]

    return {
        "schema_version": SCHEMA_VERSION,
        "gene_symbol": catalogue["gene_symbol"],
        "multi_species": catalogue["multi_species"],
        "default_scope": catalogue["default_scope"],
        "categories": (catalogue["filters"]["comparative_categories"]
                       + catalogue["filters"]["species_categories"] + [SUPPLEMENTS]),
        "figures": figures,
        "species": catalogue["filters"]["species"],
        "supplements_hidden_by_default": True,
        "counts": catalogue["counts"],
        "availability": catalogue["availability"],
        "provenance": catalogue["provenance"],
    }


def write_catalogue(outdir: Path, *, main_dir: Optional[Path] = None,
                    comparative_dir: Optional[Path] = None,
                    dataset: Optional[Dataset] = None,
                    model_index: Optional[Dict[str, Any]] = None) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if model_index is None:
        model_index = json.loads(
            (outdir / "protein_coordinate_model.json").read_text(encoding="utf-8"))
    catalogue = build_catalogue(model_index, main_dir=main_dir,
                                comparative_dir=comparative_dir, dataset=dataset)
    path = outdir / "figure_catalogue.json"
    path.write_text(json.dumps(catalogue, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    (outdir / "figure_index.json").write_text(
        json.dumps(flatten_for_gallery(catalogue), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8")
    return path


def main(argv: Optional[List[str]] = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description='The curated FGFR2 Figure Gallery catalogue.',
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    derived = ROOT / "results" / "derived" / "example"
    ap.add_argument("--indices", default=str(derived / "website_indices"))
    ap.add_argument("--main-figures", default=str(derived / "figures" / "main"))
    ap.add_argument("--comparative-figures",
                    default=str(derived / "figures" / "comparative"))
    args = ap.parse_args(argv)

    main_dir = Path(args.main_figures)
    cmp_dir = Path(args.comparative_figures)
    path = write_catalogue(
        Path(args.indices),
        main_dir=main_dir if main_dir.is_dir() else None,
        comparative_dir=cmp_dir if cmp_dir.is_dir() else None)
    catalogue = json.loads(path.read_text(encoding="utf-8"))
    print(json.dumps({"written": cm._rel(path), **catalogue["counts"],
                      "default_scope": catalogue["default_scope"],
                      "retired": catalogue["retired_source_figures"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
