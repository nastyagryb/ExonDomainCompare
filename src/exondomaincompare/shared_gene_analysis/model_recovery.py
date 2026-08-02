"""model_recovery.py — the contract between model collection and transcript selection.

Step 2 used to hand Step 3 four TSV files and nothing else. When it recovered no model
it wrote those files with headers and no rows, printed ``[OK] genes.tsv rows: 0`` and
exited 0. Step 3 then read an empty transcript table and did the only honest thing
available to it, in the worst possible form:

    raise ValueError("No transcripts found. Check --transcripts input.")

That traceback became the user's explanation for the Equus quagga run. It named the
wrong file, pointed at the wrong stage, and said nothing about the zebra, the assembly
or the misspelled taxonomy query four stages upstream.

So Step 2 now states an outcome, and the outcome is the thing Step 3 reads first. Each
status carries a sentence a user can act on and a next action the interface can offer.
An empty transcript table stops being a mystery to be diagnosed from a traceback and
becomes a reported result with a reason.

The statuses are deliberately more numerous than "worked" and "failed", because they
imply different next steps: correct a name, wait for a service, accept that a genome
has no annotation, review an ambiguous paralog, or investigate a parser.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

#: A usable model was built.
MODELS_AVAILABLE = "models_available"
#: Nothing yet, but an unexhausted route remains. Only ever an intermediate state.
RESCUE_REQUIRED = "rescue_required"
#: The submitted name is not a taxon the source knows.
TAXON_UNRESOLVED = "taxon_unresolved"
#: The source could not be reached or refused the request for transient reasons.
SOURCE_UNAVAILABLE = "source_unavailable"
#: The taxon exists but has no annotated assembly to search.
ANNOTATION_NOT_FOUND = "annotation_not_found"
#: An archive arrived and could not be read. A defect, never an empty result.
PARSER_FAILED = "parser_failed"
#: Candidates found, none assignable to the requested gene rather than a paralog.
AMBIGUOUS_PARALOG = "ambiguous_paralog"
#: The gene was found but yielded no protein-coding transcript.
NO_VALID_TRANSLATED_CDS = "no_valid_translated_cds"
#: Something was recovered but not on evidence strong enough to analyse unreviewed.
REVIEW_REQUIRED = "review_required"


#: One sentence per status, addressed to the person who started the run. The gene and
#: species are filled in, because "no transcripts found" without them is what made the
#: original failure unreadable.
_MESSAGES = {
    MODELS_AVAILABLE:
        "{gene} models were recovered for {species}.",
    TAXON_UNRESOLVED:
        "{species} could not be resolved to a species in NCBI Taxonomy, so no "
        "annotation could be searched. Check the spelling of the scientific name.",
    SOURCE_UNAVAILABLE:
        "The annotation sources could not be reached while looking up {gene} for "
        "{species}. No data were retrieved; the run can be retried.",
    ANNOTATION_NOT_FOUND:
        "No annotated genome assembly is available for {species}, so no {gene} "
        "annotation could be retrieved.",
    PARSER_FAILED:
        "The annotation for {species} was downloaded but could not be read, so no "
        "{gene} model was built. This is a processing fault, not missing data.",
    AMBIGUOUS_PARALOG:
        "{gene} candidates were found for {species}, but none passed paralog and "
        "translated-CDS validation.",
    NO_VALID_TRANSLATED_CDS:
        "A {gene} locus was found for {species}, but it yielded no protein-coding "
        "transcript that could be translated.",
    REVIEW_REQUIRED:
        "A {gene} model for {species} was recovered only from weaker evidence and "
        "needs review before it is analysed.",
    RESCUE_REQUIRED:
        "No {gene} model for {species} from the preferred source; alternative routes "
        "are still being tried.",
}

#: What the interface should offer next. ``retry_local_preparation`` is for states a
#: retry can plausibly change; the others would only repeat themselves.
_NEXT_ACTIONS = {
    MODELS_AVAILABLE: "continue_pipeline",
    TAXON_UNRESOLVED: "correct_species_name",
    SOURCE_UNAVAILABLE: "retry_local_preparation",
    ANNOTATION_NOT_FOUND: "choose_another_species",
    PARSER_FAILED: "report_processing_fault",
    AMBIGUOUS_PARALOG: "review_candidates",
    NO_VALID_TRANSLATED_CDS: "review_candidates",
    REVIEW_REQUIRED: "review_candidates",
    RESCUE_REQUIRED: "retry_local_preparation",
}


@dataclass
class RouteAttempt:
    """One rung of the cascade: what was tried, what came back."""

    route: str
    status: str
    detail: str = ""
    n_transcripts: int = 0
    n_proteins: int = 0

    def as_row(self) -> Dict[str, str]:
        return {
            "route": self.route,
            "status": self.status,
            "detail": self.detail,
            "n_transcripts": str(self.n_transcripts),
            "n_proteins": str(self.n_proteins),
        }


@dataclass
class SpeciesOutcome:
    """The recovery story for one species, in the order it happened."""

    species_id: str
    species_input: str
    gene_symbol: str
    status: str = RESCUE_REQUIRED
    detail: str = ""
    attempts: List[RouteAttempt] = field(default_factory=list)
    accepted_route: str = ""
    assembly_accession: str = ""
    assembly_status: str = ""
    taxid: str = ""
    accepted_scientific_name: str = ""
    n_genes: int = 0
    n_transcripts: int = 0
    n_exons: int = 0
    n_cds_features: int = 0
    n_translated_proteins: int = 0

    def record(self, route: str, status: str, detail: str = "",
               n_transcripts: int = 0, n_proteins: int = 0) -> None:
        self.attempts.append(RouteAttempt(route, status, detail,
                                          n_transcripts, n_proteins))

    def conclude(self, status: str, detail: str = "") -> None:
        self.status = status
        self.detail = detail or self.detail

    @property
    def species_label(self) -> str:
        return (self.accepted_scientific_name or self.species_input
                or self.species_id).replace("_", " ")

    def message(self) -> str:
        template = _MESSAGES.get(self.status, "{gene} collection for {species} ended "
                                              "in an unrecognised state.")
        text = template.format(gene=self.gene_symbol, species=self.species_label)
        return f"{text} {self.detail}".strip() if self.detail else text

    def next_action(self) -> str:
        return _NEXT_ACTIONS.get(self.status, "review_candidates")

    def as_dict(self) -> Dict[str, Any]:
        return {
            "species_id": self.species_id,
            "species_input": self.species_input,
            "accepted_scientific_name": self.accepted_scientific_name,
            "taxid": self.taxid,
            "gene_symbol": self.gene_symbol,
            "status": self.status,
            "message": self.message(),
            "next_action": self.next_action(),
            "detail": self.detail,
            "accepted_route": self.accepted_route,
            "assembly_accession": self.assembly_accession,
            "assembly_status": self.assembly_status,
            "attempts": [a.as_row() for a in self.attempts],
            "counts": {
                "genes": self.n_genes,
                "transcripts": self.n_transcripts,
                "exons": self.n_exons,
                "cds_features": self.n_cds_features,
                "translated_proteins": self.n_translated_proteins,
            },
        }


def status_from_assembly(selection_status: str) -> str:
    """Map an assembly-selection outcome onto a collection status.

    Kept as an explicit mapping so that "the name was rejected" and "the species has no
    annotated genome" stay different answers all the way to the interface.
    """
    from exondomaincompare.shared_gene_analysis import assembly_selection as asel

    return {
        asel.NO_QUERY_TERM: TAXON_UNRESOLVED,
        asel.TAXON_REJECTED: TAXON_UNRESOLVED,
        asel.SERVICE_FAILED: SOURCE_UNAVAILABLE,
        asel.NO_ASSEMBLY: ANNOTATION_NOT_FOUND,
        asel.NONE_ANNOTATED: ANNOTATION_NOT_FOUND,
        asel.ALL_REJECTED: ANNOTATION_NOT_FOUND,
        asel.DOWNLOAD_FAILED: SOURCE_UNAVAILABLE,
        asel.PARSE_FAILED: PARSER_FAILED,
    }.get(selection_status, RESCUE_REQUIRED)


def status_from_identification(identification_status: str) -> str:
    from exondomaincompare.shared_gene_analysis import gene_identification as gid

    return {
        gid.FOUND: MODELS_AVAILABLE,
        gid.NOT_FOUND: ANNOTATION_NOT_FOUND,
        gid.AMBIGUOUS: AMBIGUOUS_PARALOG,
        gid.REJECTED_PARALOG: AMBIGUOUS_PARALOG,
    }.get(identification_status, RESCUE_REQUIRED)


@dataclass
class CollectionContract:
    """What Step 2 tells the rest of the pipeline, and the interface."""

    gene_symbol: str
    outcomes: List[SpeciesOutcome] = field(default_factory=list)

    @property
    def usable(self) -> List[SpeciesOutcome]:
        return [o for o in self.outcomes if o.status == MODELS_AVAILABLE]

    @property
    def status(self) -> str:
        if self.usable:
            return MODELS_AVAILABLE
        if not self.outcomes:
            return ANNOTATION_NOT_FOUND
        # With several species and no usable model, the most specific shared reason is
        # more useful than a generic failure. A single-species run simply reports its
        # own outcome.
        statuses = [o.status for o in self.outcomes]
        for candidate in (PARSER_FAILED, AMBIGUOUS_PARALOG, NO_VALID_TRANSLATED_CDS,
                          TAXON_UNRESOLVED, ANNOTATION_NOT_FOUND, SOURCE_UNAVAILABLE,
                          REVIEW_REQUIRED):
            if candidate in statuses:
                return candidate
        return statuses[0]

    def message(self) -> str:
        if self.usable and len(self.usable) == len(self.outcomes):
            return (f"{self.gene_symbol} models were recovered for all "
                    f"{len(self.outcomes)} species.")
        blocked = [o for o in self.outcomes if o.status != MODELS_AVAILABLE]
        if self.usable:
            return (f"{self.gene_symbol} models were recovered for {len(self.usable)} "
                    f"of {len(self.outcomes)} species. "
                    + " ".join(o.message() for o in blocked[:3]))
        return " ".join(o.message() for o in blocked[:3]) or "No models were recovered."

    def next_action(self) -> str:
        return _NEXT_ACTIONS.get(self.status, "review_candidates")

    def as_dict(self) -> Dict[str, Any]:
        return {
            "contract_version": 1,
            "gene_symbol": self.gene_symbol,
            "status": self.status,
            "message": self.message(),
            "next_action": self.next_action(),
            "n_species": len(self.outcomes),
            "n_species_with_models": len(self.usable),
            "species": [o.as_dict() for o in self.outcomes],
            "generated_by": "scripts/shared_gene_analysis/model_recovery.py",
        }

    def write(self, outdir: Path) -> Path:
        path = Path(outdir) / "collection_status.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.as_dict(), indent=2, ensure_ascii=False) + "\n",
                        encoding="utf-8")
        return path


def read_contract(path: Path) -> Optional[Dict[str, Any]]:
    """Read the contract, or ``None`` when a step ran before it existed."""
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def explain_empty_input(contract: Optional[Dict[str, Any]], gene_symbol: str,
                        stage: str) -> str:
    """The message a downstream step shows instead of raising on empty input.

    A step that finds no transcripts is not the place where the run went wrong; it is
    the place where the earlier failure became visible. Its job is to quote the recorded
    reason, not to invent one from the shape of its own input.
    """
    if not contract:
        return (f"{stage} received no {gene_symbol} transcripts and no collection "
                "status was recorded, so the reason could not be determined. See the "
                "run logs for the model collection step.")
    message = str(contract.get("message") or "").strip()
    status = str(contract.get("status") or "unknown")
    if message:
        return f"{message} ({status})"
    return (f"{stage} received no {gene_symbol} transcripts; model collection reported "
            f"status {status}.")


def consistency_checks(n_genes: int, n_transcripts: int, n_exons: int,
                       n_cds: int) -> List[Dict[str, str]]:
    """Checks that cannot pass vacuously.

    The failed run's ``internal_consistency_checks.tsv`` reported six of six PASS over
    four empty tables — ``orphan_transcripts=0``, ``genes_without_transcripts=0`` — so
    the one artefact whose job is to catch an inconsistent result actively reassured
    the reader that an empty run was sound. A check over no data is reported as
    NOT_APPLICABLE.
    """
    if n_genes == 0 and n_transcripts == 0:
        return [{
            "check_name": "model_tables_are_populated",
            "status": "FAIL",
            "affected_species": "",
            "details": ("no gene or transcript rows were written; downstream checks "
                        "over these tables are not applicable rather than passing"),
        }]
    rows = [{
        "check_name": "model_tables_are_populated",
        "status": "PASS",
        "affected_species": "",
        "details": (f"genes={n_genes} transcripts={n_transcripts} exons={n_exons} "
                    f"cds_features={n_cds}"),
    }]
    if n_cds == 0:
        rows.append({
            "check_name": "transcripts_have_cds_features",
            "status": "FAIL",
            "affected_species": "",
            "details": "transcripts were recovered but no CDS feature was parsed",
        })
    return rows


def summarise(outcomes: Sequence[SpeciesOutcome]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for outcome in outcomes:
        counts[outcome.status] = counts.get(outcome.status, 0) + 1
    return counts
