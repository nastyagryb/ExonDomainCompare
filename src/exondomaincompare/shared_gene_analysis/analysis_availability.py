#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional
from exondomaincompare.contracts import file_sha256

# --------------------------------------------------------------------------- #
# Canonical states. Deliberately not collapsible into one "unavailable".
# --------------------------------------------------------------------------- #
#: Prerequisites exist and valid scientific output was produced.
AVAILABLE = "available"
#: The biological or mathematical prerequisites do not exist, determined successfully.
NOT_APPLICABLE = "not_applicable"
#: The evidence could not be supported from the available annotation or models.
SCIENTIFICALLY_UNAVAILABLE = "scientifically_unavailable"
#: The required stage has not finished.
PENDING = "pending"
#: An expected output of a supposedly completed applicable stage is absent.
TECHNICALLY_MISSING = "technically_missing"
#: The output exists but no longer matches the current canonical inputs.
STALE = "stale"
#: An applicable analysis attempted to run and failed.
FAILED = "failed"

STATES = (AVAILABLE, NOT_APPLICABLE, SCIENTIFICALLY_UNAVAILABLE, PENDING,
          TECHNICALLY_MISSING, STALE, FAILED)

#: States that mean the run is not finished. ``not_applicable`` and
#: ``scientifically_unavailable`` are settled answers and are absent on purpose.
BLOCKING_STATES = (TECHNICALLY_MISSING, PENDING, STALE, FAILED)

#: States that are settled, whether or not they produced a plot.

# --------------------------------------------------------------------------- #
# Reason codes — stable identifiers for the machine, messages for the reader.
# --------------------------------------------------------------------------- #
SINGLE_CODING_EXON = "single_coding_exon"
NO_INTERNAL_BOUNDARIES = "no_internal_coding_exon_boundaries"
SINGLE_PROTEIN_SEQUENCE = "single_unique_protein_sequence"
SINGLE_SPECIES = "single_species_dataset"
CLUSTER_PENDING = "cluster_annotation_pending"
OUTPUT_MISSING = "expected_output_missing"
OUTPUT_STALE = "output_older_than_inputs"
PRODUCED = "valid_output_produced"

MESSAGES = {
    SINGLE_CODING_EXON: ("The selected protein is encoded by one coding exon and therefore "
                         "has no internal coding-exon boundaries to analyse."),
    NO_INTERNAL_BOUNDARIES: ("The selected protein has no internal coding-exon boundaries "
                            "to analyse."),
    SINGLE_PROTEIN_SEQUENCE: ("Only one distinct translated protein sequence is available. "
                              "At least two distinct protein sequences are required to "
                              "detect protein-level isoform differences."),
    SINGLE_SPECIES: ("This analysis compares species and needs at least two species in the "
                     "dataset."),
    CLUSTER_PENDING: "The cluster annotation for this run has not finished yet.",
    OUTPUT_MISSING: ("Expected outputs were not generated even though the stage is marked "
                     "complete. Retry local analysis."),
    OUTPUT_STALE: ("This result is older than the data it summarises. Rebuild the run's "
                   "indices."),
    PRODUCED: "",
}

# Short labels for a compact badge next to a resolved not-applicable analysis.
BADGES = {
    SINGLE_CODING_EXON: "No internal boundaries",
    NO_INTERNAL_BOUNDARIES: "No internal boundaries",
    SINGLE_PROTEIN_SEQUENCE: "One protein sequence",
    SINGLE_SPECIES: "Single species",
}

# --------------------------------------------------------------------------- #
# Prerequisite counts, read from the canonical tables
# --------------------------------------------------------------------------- #
CORE = "results/core_gene_analysis"


def _rows(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        return []
    with open(path, encoding="utf-8", newline="") as handle:
        return [r for r in csv.DictReader(handle, delimiter="\t") if any(r.values())]


def _fasta(path: Path) -> Dict[str, str]:
    if not path.is_file():
        return {}
    out: Dict[str, str] = {}
    key, buf = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if key:
                out[key] = "".join(buf)
            key, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if key:
        out[key] = "".join(buf)
    return out


@dataclass
class Prerequisites:
    species_count: int = 0
    #: Distinct annotated transcript models, whether or not they are coding.
    transcript_model_count: int = 0
    #: Transcript models carrying a CDS.
    coding_transcript_count: int = 0
    #: Distinct protein accessions.
    protein_product_count: int = 0
    #: Distinct amino-acid sequences — the count that governs isoform comparison.
    unique_protein_sequence_count: int = 0
    #: Transcripts that differ from another transcript of the same protein sequence.
    transcript_only_variant_count: int = 0
    #: Coding exons of the selected primary model, from its CDS parts.
    coding_exon_count: int = 0
    #: ``max(coding_exon_count - 1, 0)``. Never inferred from a total exon count.
    internal_coding_exon_boundary_count: int = 0
    #: Protein-level differences actually supported by two or more sequences.
    protein_difference_candidate_count: int = 0
    #: Protein accession → the transcripts that encode that exact sequence.
    transcripts_by_protein_sequence: Dict[str, List[str]] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "species_count": self.species_count,
            "transcript_model_count": self.transcript_model_count,
            "coding_transcript_count": self.coding_transcript_count,
            "protein_product_count": self.protein_product_count,
            "unique_protein_sequence_count": self.unique_protein_sequence_count,
            "transcript_only_variant_count": self.transcript_only_variant_count,
            "coding_exon_count": self.coding_exon_count,
            "internal_coding_exon_boundary_count":
                self.internal_coding_exon_boundary_count,
            "protein_difference_candidate_count": self.protein_difference_candidate_count,
            "transcripts_by_protein_sequence":
                {k: list(v) for k, v in self.transcripts_by_protein_sequence.items()},
        }


def has_core_tables(run_dir: Path) -> bool:
    core = Path(run_dir) / CORE
    return ((core / "exon_protein_map.tsv").is_file()
            and (core / "protein_isoform_index.tsv").is_file())


def internal_boundary_count(coding_exon_count: int) -> int:
    return max(int(coding_exon_count) - 1, 0)


def prerequisites(run_dir: Path, species_id: str = "") -> Prerequisites:
    run_dir = Path(run_dir)
    core = run_dir / CORE
    exon_map = _rows(core / "exon_protein_map.tsv")
    isoforms = _rows(core / "protein_isoform_index.tsv")
    gene_models = _rows(core / "gene_model_index.tsv")
    sequences = _fasta(core / "proteins_all_isoforms.faa") \
        or _fasta(core / "proteins_primary.faa")

    species = sorted({r.get("species_id", "") for r in exon_map
                      or isoforms or gene_models if r.get("species_id")})
    want = species_id or (species[0] if species else "")

    scoped_iso = [r for r in isoforms if not want or r.get("species_id") == want]
    scoped_models = [r for r in gene_models if not want or r.get("species_id") == want]

    transcripts = {r.get("transcript_id", "") for r in scoped_models
                   if r.get("transcript_id")}
    coding = {r.get("transcript_id", "") for r in scoped_models
              if r.get("transcript_id") and r.get("protein_id")}
    proteins = {r.get("protein_id", "") for r in scoped_iso if r.get("protein_id")} \
        or {r.get("protein_id", "") for r in scoped_models if r.get("protein_id")}

    # Group transcripts by the amino-acid sequence they encode, so identical products are
    # one isoform with several transcript models rather than several isoforms.
    by_sequence: Dict[str, List[str]] = {}
    representative: Dict[str, str] = {}
    for row in (scoped_models or scoped_iso):
        pid, tid = row.get("protein_id", ""), row.get("transcript_id", "")
        if not pid:
            continue
        seq = sequences.get(pid, "") or f"__unknown__{pid}"
        by_sequence.setdefault(seq, [])
        if tid and tid not in by_sequence[seq]:
            by_sequence[seq].append(tid)
        representative.setdefault(seq, pid)

    unique_sequences = len(by_sequence) if by_sequence else len(proteins)
    transcript_only = sum(max(len(t) - 1, 0) for t in by_sequence.values())

    # The primary model's coding exons. The exon map holds one row per coding exon, so its
    # row count for that protein *is* the coding-exon count.
    primary = ""
    for row in scoped_iso:
        if (row.get("primary_status") or "").lower() == "primary":
            primary = row.get("protein_id", "")
            break
    if not primary and scoped_iso:
        primary = scoped_iso[0].get("protein_id", "")
    if not primary and scoped_models:
        primary = scoped_models[0].get("protein_id", "")
    coding_exons = len({r.get("exon_id") or r.get("exon_number")
                        for r in exon_map
                        if (not want or r.get("species_id") == want)
                        and (not primary or r.get("protein_id") == primary)})

    return Prerequisites(
        species_count=len(species),
        transcript_model_count=len(transcripts),
        coding_transcript_count=len(coding),
        protein_product_count=len(proteins),
        unique_protein_sequence_count=unique_sequences,
        transcript_only_variant_count=transcript_only,
        coding_exon_count=coding_exons,
        internal_coding_exon_boundary_count=internal_boundary_count(coding_exons),
        # A protein-level difference needs two sequences to differ. With one sequence the
        # count is zero by definition, not by failure to look.
        protein_difference_candidate_count=0 if unique_sequences < 2 else -1,
        transcripts_by_protein_sequence={representative[s]: t
                                         for s, t in by_sequence.items()},
    )


# --------------------------------------------------------------------------- #
# Analysis records
# --------------------------------------------------------------------------- #
@dataclass
class AnalysisState:
    analysis_name: str
    label: str
    status: str
    prerequisite_name: str = ""
    prerequisite_count: Optional[int] = None
    reason_code: str = ""
    user_message: str = ""
    badge: str = ""
    #: False for analyses a run may finish without, e.g. boundary or isoform comparison.
    required: bool = False

    @property
    def blocks_ready(self) -> bool:
        return self.status in BLOCKING_STATES

    @property
    def applicable(self) -> bool:
        return self.status != NOT_APPLICABLE

    def as_dict(self) -> Dict[str, Any]:
        return {
            "analysis_name": self.analysis_name,
            "label": self.label,
            "status": self.status,
            "prerequisite_name": self.prerequisite_name,
            "prerequisite_count": self.prerequisite_count,
            "reason_code": self.reason_code,
            "user_message": self.user_message,
            "badge": self.badge,
            "required": self.required,
            "blocks_results_ready": self.blocks_ready,
        }


def _state(name: str, label: str, status: str, *, reason: str = "",
           prerequisite: str = "", count: Optional[int] = None,
           required: bool = False, message: str = "") -> AnalysisState:
    return AnalysisState(
        analysis_name=name, label=label, status=status,
        prerequisite_name=prerequisite, prerequisite_count=count,
        reason_code=reason,
        user_message=message or MESSAGES.get(reason, ""),
        badge=BADGES.get(reason, "") if status == NOT_APPLICABLE else "",
        required=required)


def _cluster_outputs_present(run_dir: Path) -> bool:
    out = run_dir / "results" / "14_interproscan" / "primary" / "output"
    return out.is_dir() and any(p.stat().st_size > 0 for p in out.glob("*.tsv"))


#: The models every downstream analysis describes. An analysis written before these were
#: last collected no longer describes them.
CANONICAL_INPUT = "proteins_primary.faa"


def _is_stale(run_dir: Path, output: str) -> bool:
    core = Path(run_dir) / CORE
    models, produced = core / CANONICAL_INPUT, core / output
    if not (models.is_file() and produced.is_file()):
        return False
    try:
        return produced.stat().st_mtime < models.stat().st_mtime
    except OSError:
        return False


def boundary_analysis(run_dir: Path, pre: Prerequisites,
                      cluster_ready: Optional[bool] = None) -> AnalysisState:
    name, label = "boundary_analysis", "Exon–domain boundary analysis"
    count = pre.internal_coding_exon_boundary_count
    if count == 0:
        reason = (SINGLE_CODING_EXON if pre.coding_exon_count == 1
                  else NO_INTERNAL_BOUNDARIES)
        return _state(name, label, NOT_APPLICABLE, reason=reason,
                      prerequisite="internal_coding_exon_boundary_count", count=count)

    if cluster_ready is None:
        cluster_ready = _cluster_outputs_present(run_dir)
    if not cluster_ready:
        return _state(name, label, PENDING, reason=CLUSTER_PENDING,
                      prerequisite="internal_coding_exon_boundary_count", count=count)

    table = run_dir / CORE / "exon_domain_boundary_distances.tsv"
    if not table.is_file():
        return _state(name, label, TECHNICALLY_MISSING, reason=OUTPUT_MISSING,
                      prerequisite="internal_coding_exon_boundary_count", count=count)
    if not _rows(table):
        # Boundaries exist but none was classified: the domain layer had nothing to relate
        # them to. A scientific gap, not a missing file.
        return _state(name, label, SCIENTIFICALLY_UNAVAILABLE,
                      prerequisite="internal_coding_exon_boundary_count", count=count,
                      message="No representative domain could be related to this "
                              "protein's coding-exon boundaries.")
    if _is_stale(run_dir, "exon_domain_boundary_distances.tsv"):
        return _state(name, label, STALE, reason=OUTPUT_STALE,
                      prerequisite="internal_coding_exon_boundary_count", count=count)
    return _state(name, label, AVAILABLE, reason=PRODUCED,
                  prerequisite="internal_coding_exon_boundary_count", count=count)


def protein_isoform_comparison(pre: Prerequisites) -> AnalysisState:
    name = "protein_isoform_comparison"
    label = "Protein isoform comparison"
    count = pre.unique_protein_sequence_count
    if count < 2:
        return _state(name, label, NOT_APPLICABLE, reason=SINGLE_PROTEIN_SEQUENCE,
                      prerequisite="unique_protein_sequence_count", count=count)
    return _state(name, label, AVAILABLE, reason=PRODUCED,
                  prerequisite="unique_protein_sequence_count", count=count)


def candidate_analysis(pre: Prerequisites) -> AnalysisState:
    name = "protein_difference_candidate_analysis"
    label = "Protein-difference candidate analysis"
    count = pre.unique_protein_sequence_count
    if count < 2:
        return _state(name, label, NOT_APPLICABLE, reason=SINGLE_PROTEIN_SEQUENCE,
                      prerequisite="unique_protein_sequence_count", count=count)
    return _state(name, label, AVAILABLE, reason=PRODUCED,
                  prerequisite="unique_protein_sequence_count", count=count)


def exon_map_analysis(run_dir: Path, pre: Prerequisites) -> AnalysisState:
    name, label = "exon_map", "Exon map"
    count = pre.coding_exon_count
    if count == 0:
        return _state(name, label, SCIENTIFICALLY_UNAVAILABLE,
                      prerequisite="coding_exon_count", count=0,
                      message="No coding exon was recovered for the selected model.")
    if not (Path(run_dir) / CORE / "exon_protein_map.tsv").is_file():
        return _state(name, label, TECHNICALLY_MISSING, reason=OUTPUT_MISSING,
                      prerequisite="coding_exon_count", count=count, required=True)
    return _state(name, label, AVAILABLE, reason=PRODUCED,
                  prerequisite="coding_exon_count", count=count, required=True)


def domain_architecture_analysis(run_dir: Path, pre: Prerequisites,
                                 cluster_ready: Optional[bool] = None) -> AnalysisState:
    name, label = "domain_architecture", "Domain architecture"
    run_dir = Path(run_dir)
    if cluster_ready is None:
        cluster_ready = _cluster_outputs_present(run_dir)
    if not cluster_ready:
        return _state(name, label, PENDING, reason=CLUSTER_PENDING,
                      prerequisite="unique_protein_sequence_count",
                      count=pre.unique_protein_sequence_count)
    table = run_dir / CORE / "domain_features.tsv"
    if not table.is_file():
        return _state(name, label, TECHNICALLY_MISSING, reason=OUTPUT_MISSING,
                      prerequisite="unique_protein_sequence_count",
                      count=pre.unique_protein_sequence_count, required=True)
    if not _rows(table):
        return _state(name, label, SCIENTIFICALLY_UNAVAILABLE,
                      prerequisite="unique_protein_sequence_count",
                      count=pre.unique_protein_sequence_count,
                      message="No domain signature was matched on this protein.")
    if _is_stale(run_dir, "domain_features.tsv"):
        return _state(name, label, STALE, reason=OUTPUT_STALE,
                      prerequisite="unique_protein_sequence_count",
                      count=pre.unique_protein_sequence_count, required=True)
    return _state(name, label, AVAILABLE, reason=PRODUCED,
                  prerequisite="unique_protein_sequence_count",
                  count=pre.unique_protein_sequence_count, required=True)




# --------------------------------------------------------------------------- #
# Manifest
# --------------------------------------------------------------------------- #
@dataclass
class Manifest:
    run_id: str
    species_id: str
    prerequisites: Prerequisites
    analyses: List[AnalysisState]

    @property
    def blocking(self) -> List[AnalysisState]:
        return [a for a in self.analyses if a.blocks_ready]

    @property
    def ready(self) -> bool:
        return not self.blocking

    def reason(self) -> str:
        if self.ready:
            resolved = [a for a in self.analyses if a.status == NOT_APPLICABLE]
            if resolved:
                names = ", ".join(a.label.lower() for a in resolved)
                return (f"every applicable analysis is complete; not applicable to this "
                        f"model: {names}")
            return "every applicable analysis is complete"
        first = self.blocking[0]
        return f"{first.label}: {first.user_message or first.status}"

    def by_name(self) -> Dict[str, AnalysisState]:
        return {a.analysis_name: a for a in self.analyses}

    def as_dict(self) -> Dict[str, Any]:
        return {
            "schema_version": 1,
            "run_id": self.run_id,
            "species_id": self.species_id,
            "ready": self.ready,
            "reason": self.reason(),
            "prerequisites": self.prerequisites.as_dict(),
            "analyses": [a.as_dict() for a in self.analyses],
            "blocking": [a.analysis_name for a in self.blocking],
            "not_applicable": [a.analysis_name for a in self.analyses
                               if a.status == NOT_APPLICABLE],
        }


def build_manifest(run_dir: Path, species_id: str = "",
                   cluster_ready: Optional[bool] = None) -> Manifest:
    run_dir = Path(run_dir)
    pre = prerequisites(run_dir, species_id)
    if cluster_ready is None:
        cluster_ready = _cluster_outputs_present(run_dir)

    analyses = [
        exon_map_analysis(run_dir, pre),
        domain_architecture_analysis(run_dir, pre, cluster_ready),
        boundary_analysis(run_dir, pre, cluster_ready),
        protein_isoform_comparison(pre),
        candidate_analysis(pre),
    ]
    return Manifest(run_id=run_dir.name, species_id=species_id, prerequisites=pre,
                    analyses=analyses)


def write_manifest(run_dir: Path, species_id: str = "") -> Path:
    manifest = build_manifest(run_dir, species_id)
    out = Path(run_dir) / "website_indices" / "analysis_availability.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest.as_dict(), indent=2) + "\n", encoding="utf-8")
    return out


def index_version(run_dir: Path) -> str:
    import hashlib

    indices = Path(run_dir) / "website_indices"
    digest = hashlib.sha256()
    if indices.is_dir():
        for path in sorted(indices.rglob("*.json")):
            try:
                checksum = file_sha256(path)
            except OSError:
                continue
            digest.update(str(path.relative_to(indices)).encode("utf-8"))
            digest.update(checksum.encode("ascii"))
    status = Path(run_dir) / "status.json"
    if status.is_file():
        try:
            digest.update(file_sha256(status).encode("ascii"))
        except OSError:
            pass
    return digest.hexdigest()[:16]


#: Model section → the analysis whose availability governs it.
MODEL_SECTIONS = {
    "exon_domain_boundaries": "boundary_analysis",
    "isoform_alignment": "protein_isoform_comparison",
    "candidate_evidence": "protein_difference_candidate_analysis",
    "domain_architecture": "domain_architecture",
}

#: Substrings identifying the figures an analysis would have produced. Matched against the
#: figure id so a run's own naming does not have to be enumerated exhaustively.
FIGURE_MARKERS = {
    "boundary_analysis": ("boundary",),
    "protein_difference_candidate_analysis": ("candidate",),
    "protein_isoform_comparison": ("isoform_alignment", "msa_isoform"),
}


def _drop_inapplicable_figures(model: Dict[str, Any],
                               states: Dict[str, AnalysisState]) -> None:
    inapplicable = [name for name, state in states.items()
                    if state.status == NOT_APPLICABLE and name in FIGURE_MARKERS]
    if not inapplicable:
        return
    markers = tuple(m for name in inapplicable for m in FIGURE_MARKERS[name])

    def keep(figure: Any) -> bool:
        if not isinstance(figure, dict):
            return True
        if figure.get("status") == AVAILABLE:
            return True
        fid = str(figure.get("figure_id") or "").lower()
        return not any(marker in fid for marker in markers)

    for key in ("figures", "figure_gallery"):
        block = model.get(key)
        if isinstance(block, list):
            model[key] = [f for f in block if keep(f)]
        elif isinstance(block, dict) and isinstance(block.get("figures"), list):
            block["figures"] = [f for f in block["figures"] if keep(f)]


def annotate_dataset_model(model: Dict[str, Any], run_dir: Path) -> Dict[str, Any]:
    run_dir = Path(run_dir)
    if not has_core_tables(run_dir):
        return model
    try:
        manifest = build_manifest(run_dir)
    except Exception:
        return model

    model["analysis_availability"] = manifest.as_dict()
    states = manifest.by_name()

    for section, analysis in MODEL_SECTIONS.items():
        block = model.get(section)
        state = states.get(analysis)
        if not isinstance(block, dict) or state is None:
            continue
        produced = block.get("available") is True
        if state.status == NOT_APPLICABLE:
            block["available"] = False
            block["status"] = NOT_APPLICABLE
        # A section that genuinely produced output keeps its own verdict; only its
        # availability block is filled in so the frontend reads one uniform shape.
        block["availability"] = availability_block(
            state if not (produced and state.status != NOT_APPLICABLE)
            else _state(state.analysis_name, state.label, AVAILABLE, reason=PRODUCED,
                        prerequisite=state.prerequisite_name,
                        count=state.prerequisite_count))

    _drop_inapplicable_figures(model, states)
    return model


def availability_block(state: AnalysisState) -> Dict[str, Any]:
    return {
        "state": state.status,
        "label": state.label,
        "reason": state.user_message,
        "reason_code": state.reason_code,
        "badge": state.badge,
        "prerequisite_name": state.prerequisite_name,
        "prerequisite_count": state.prerequisite_count,
    }


__all__ = [
    "AVAILABLE", "NOT_APPLICABLE", "SCIENTIFICALLY_UNAVAILABLE", "PENDING",
    "TECHNICALLY_MISSING", "STALE", "FAILED", "STATES", "BLOCKING_STATES",
    "RESOLVED_STATES", "SINGLE_CODING_EXON", "NO_INTERNAL_BOUNDARIES",
    "SINGLE_PROTEIN_SEQUENCE", "SINGLE_SPECIES", "MESSAGES", "BADGES",
    "Prerequisites", "AnalysisState", "Manifest",
    "internal_boundary_count", "prerequisites", "boundary_analysis",
    "protein_isoform_comparison", "candidate_analysis", "exon_map_analysis",
    "domain_architecture_analysis", "cross_species_analysis",
    "build_manifest", "write_manifest", "availability_block",
    "index_version", "annotate_dataset_model", "MODEL_SECTIONS",
    "has_core_tables", "FIGURE_MARKERS",
]
