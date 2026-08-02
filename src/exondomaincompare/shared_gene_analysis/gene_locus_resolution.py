#!/usr/bin/env python3
"""Which locus in a genome annotation is the gene the user asked for.

The generic resolver had two routes: the exact symbol on the annotation's gene row, and a
``gene_synonym`` attribute on that same row. Both fail for a large class of assemblies, and
they fail silently in a way that reads like biology.

*Panthera leo* is the example. The user asks for HBA. The RefSeq assembly annotates the
locus as ``LOC122209636`` and — like every one of the 32 109 gene rows in that annotation —
carries neither ``description=`` nor ``gene_synonym=``. The product name lives on the child
mRNA rows, and the alias HBA exists only in NCBI Gene. So both routes miss, and the run
reported "no locus with this symbol or synonym", which is true of the GFF3 attribute and
false about the species: NCBI Gene resolves HBA to GeneID 122209636, and that GeneID sits in
the annotation's own ``Dbxref``.

Matching on the description would not have saved it, and would have been worse than failing.
The alpha-globin cluster on that chromosome holds ``LOC122209634`` "hemoglobin subunit
alpha", ``LOC122209636`` "hemoglobin subunit alpha-like", two zeta loci, ``HBM`` and
``HBQ1``. A description route looking for "hemoglobin subunit alpha" prefers the locus that
is *not* the one NCBI calls HBA. Only the GeneID is decisive; the description is at most
corroboration.
"""
from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

# --------------------------------------------------------------------------- #
# Resolution routes, in the order they are consulted
# --------------------------------------------------------------------------- #
ROUTE_EXACT_SYMBOL = "annotation_exact_symbol"
ROUTE_NORMALIZED_SYMBOL = "annotation_case_normalized_symbol"
ROUTE_ANNOTATION_ALIAS = "annotation_provided_alias"
ROUTE_NCBI_SYMBOL = "ncbi_gene_official_symbol"
ROUTE_NCBI_ALIAS_GENEID = "ncbi_alias_to_geneid_to_annotation"
ROUTE_NCBI_ALIAS_SYMBOL = "ncbi_alias_to_official_symbol"
ROUTE_DESCRIPTION = "source_description_supporting_only"
ROUTE_ORTHOLOGY = "orthology_aware_rescue"

#: Documentation of the cascade. Position in this tuple is the priority.

#: Routes whose evidence is a name similarity rather than an identifier. They may
#: corroborate a candidate but must never select one on their own — in the lion
#: alpha-globin cluster the closest description belongs to the wrong locus.
SUPPORTING_ONLY: Tuple[str, ...] = (ROUTE_DESCRIPTION,)

# --------------------------------------------------------------------------- #
# Outcome statuses
# --------------------------------------------------------------------------- #
RESOLVED = "gene_resolved"
#: No record for this gene in this species anywhere — annotation or NCBI Gene.
GENE_NOT_FOUND = "gene_not_found"
#: NCBI Gene knows the alias, but its GeneID and official symbol are both absent from the
#: annotation. A defect or an assembly mismatch, not an absent gene.
ALIAS_MAPPING_FAILED = "alias_resolved_but_annotation_mapping_failed"
#: Several loci are equally plausible. Better to say so than to pick one.
AMBIGUOUS_FAMILY = "ambiguous_gene_family"
TRANSCRIPT_NOT_FOUND = "transcript_not_found"
NO_VALID_TRANSLATED_CDS = "no_valid_translated_cds"
PARSER_FAILED = "parser_failed"
SOURCE_UNAVAILABLE = "source_unavailable"
REVIEW_REQUIRED = "review_required"

STATUSES = (RESOLVED, GENE_NOT_FOUND, ALIAS_MAPPING_FAILED, AMBIGUOUS_FAMILY,
            TRANSCRIPT_NOT_FOUND, NO_VALID_TRANSLATED_CDS, PARSER_FAILED,
            SOURCE_UNAVAILABLE, REVIEW_REQUIRED)

#: Human-readable cause per status. The first two are the pair the previous resolver
#: collapsed into one message.
STATUS_MESSAGE = {
    GENE_NOT_FOUND: ("No {gene}-related record exists for this species in the selected "
                     "annotation or in NCBI Gene."),
    ALIAS_MAPPING_FAILED: ("{gene} was found in NCBI Gene as an alias of {official}, but "
                           "the corresponding assembly locus could not be mapped in the "
                           "selected annotation."),
    AMBIGUOUS_FAMILY: ("{gene} matches several equally plausible loci in this genome "
                       "({candidates}). Enter the specific symbol you mean."),
    TRANSCRIPT_NOT_FOUND: "The {gene} locus was resolved but carries no transcript.",
    NO_VALID_TRANSLATED_CDS: ("The {gene} locus was resolved but no transcript carries a "
                              "translatable CDS."),
    PARSER_FAILED: "The genome annotation for this species could not be parsed: {detail}",
    SOURCE_UNAVAILABLE: "The gene-identity source could not be reached: {detail}",
    REVIEW_REQUIRED: ("{gene} resolved to a locus that needs manual review before it can "
                      "be analysed: {detail}"),
}

_SYMBOL_FIELDS = ("gene", "Name", "gene_name", "locus_tag")
_ALIAS_FIELDS = ("gene_synonym", "synonym", "Alias")
_LOC_RE = re.compile(r"^LOC\d+$", re.IGNORECASE)
_GENE_FEATURES = frozenset({"gene", "pseudogene"})
_TRANSCRIPT_FEATURES = frozenset({"mRNA", "transcript"})

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


# --------------------------------------------------------------------------- #
# Annotation model
# --------------------------------------------------------------------------- #
@dataclass
class AnnotationLocus:
    """One gene row of the annotation plus what its children say about it."""

    gene_id: str                      # GFF3 ID, e.g. gene-LOC122209636
    symbol: str
    source_gene_id: str = ""          # numeric GeneID from Dbxref
    dbxref: str = ""
    aliases: List[str] = field(default_factory=list)
    description: str = ""
    biotype: str = ""
    feature_type: str = "gene"
    seqid: str = ""
    start: int = 0
    end: int = 0
    strand: str = ""
    #: Product names of child transcripts. On Gnomon annotations this is the only place a
    #: description exists at all, so it is collected — as evidence, not as a key.
    products: List[str] = field(default_factory=list)
    transcript_count: int = 0
    protein_count: int = 0

    @property
    def is_loc_labelled(self) -> bool:
        return bool(_LOC_RE.match(self.symbol or ""))

    @property
    def is_pseudogene(self) -> bool:
        return (self.feature_type == "pseudogene"
                or "pseudogene" in (self.biotype or "").lower())

    def as_dict(self) -> Dict[str, Any]:
        return {
            "gene_id": self.gene_id, "symbol": self.symbol,
            "source_gene_id": self.source_gene_id, "dbxref": self.dbxref,
            "aliases": list(self.aliases), "description": self.description,
            "biotype": self.biotype, "feature_type": self.feature_type,
            "seqid": self.seqid, "start": self.start, "end": self.end,
            "strand": self.strand, "products": list(self.products),
            "transcript_count": self.transcript_count,
            "protein_count": self.protein_count,
            "is_loc_labelled": self.is_loc_labelled,
            "is_pseudogene": self.is_pseudogene,
        }


def _attrs(field9: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in field9.strip().rstrip(";").split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            out[key.strip()] = urllib.parse.unquote(value.strip())
    return out


def _geneid_from_dbxref(dbxref: str) -> str:
    for token in (dbxref or "").split(","):
        token = token.strip()
        if token.upper().startswith("GENEID:"):
            return token.split(":", 1)[1].strip()
    return ""


def _split_aliases(attr: Dict[str, str]) -> List[str]:
    out: List[str] = []
    for key in _ALIAS_FIELDS:
        for token in re.split(r"[,;|]", attr.get(key, "") or ""):
            token = token.strip()
            if token and token not in out:
                out.append(token)
    return out


def scan_annotation(gff_path: Path) -> Dict[str, AnnotationLocus]:
    """Every gene locus of the annotation, keyed by GFF3 ID.

    One pass. Child transcripts contribute their product names and their counts, because on
    a Gnomon annotation the gene row itself carries no description at all — which is
    precisely why a resolver that only reads gene rows cannot see what the locus is.
    """
    loci: Dict[str, AnnotationLocus] = {}
    products: Dict[str, List[str]] = {}
    transcripts: Dict[str, int] = {}
    proteins: Dict[str, set] = {}
    tx_parent: Dict[str, str] = {}

    with open(gff_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or "\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            if ftype in _GENE_FEATURES:
                attr = _attrs(f[8])
                gene_id = attr.get("ID", "")
                if not gene_id:
                    continue
                symbol = next((attr[k] for k in _SYMBOL_FIELDS if attr.get(k)), "")
                dbxref = attr.get("Dbxref", "")
                loci[gene_id] = AnnotationLocus(
                    gene_id=gene_id, symbol=symbol,
                    source_gene_id=_geneid_from_dbxref(dbxref), dbxref=dbxref,
                    aliases=_split_aliases(attr),
                    description=attr.get("description", "") or attr.get("Note", ""),
                    biotype=attr.get("gene_biotype", ""), feature_type=ftype,
                    seqid=f[0], start=int(f[3]), end=int(f[4]), strand=f[6])
            elif ftype in _TRANSCRIPT_FEATURES:
                attr = _attrs(f[8])
                parent = attr.get("Parent", "")
                if not parent:
                    continue
                transcripts[parent] = transcripts.get(parent, 0) + 1
                if attr.get("ID"):
                    tx_parent[attr["ID"]] = parent
                product = attr.get("product", "")
                if product:
                    products.setdefault(parent, [])
                    if product not in products[parent]:
                        products[parent].append(product)
            elif ftype == "CDS":
                attr = _attrs(f[8])
                parent = tx_parent.get(attr.get("Parent", ""), "")
                pid = attr.get("protein_id", "")
                if parent and pid:
                    proteins.setdefault(parent, set()).add(pid)

    for gene_id, locus in loci.items():
        locus.products = products.get(gene_id, [])
        locus.transcript_count = transcripts.get(gene_id, 0)
        locus.protein_count = len(proteins.get(gene_id, ()))
    return loci


# --------------------------------------------------------------------------- #
# NCBI Gene lookup
# --------------------------------------------------------------------------- #
@dataclass
class NcbiGeneRecord:
    gene_id: str
    official_symbol: str
    description: str = ""
    aliases: List[str] = field(default_factory=list)
    other_designations: List[str] = field(default_factory=list)
    chromosome: str = ""
    chr_accession: str = ""
    exon_count: Optional[int] = None
    status: str = ""
    current_id: str = ""

    @property
    def is_live(self) -> bool:
        """A replaced or discontinued record must not be treated as current."""
        return not (self.status or "").strip() and not (self.current_id or "").strip()

    def mentions(self, symbol: str) -> bool:
        want = (symbol or "").strip().upper()
        if not want:
            return False
        pool = [self.official_symbol, *self.aliases]
        return any((v or "").strip().upper() == want for v in pool)

    def as_dict(self) -> Dict[str, Any]:
        return {"gene_id": self.gene_id, "official_symbol": self.official_symbol,
                "description": self.description, "aliases": list(self.aliases),
                "other_designations": list(self.other_designations),
                "chromosome": self.chromosome, "chr_accession": self.chr_accession,
                "exon_count": self.exon_count, "status": self.status,
                "current_id": self.current_id, "is_live": self.is_live}


def _fetch_json(url: str, timeout: float, sleep_s: float) -> Any:
    if sleep_s:
        time.sleep(sleep_s)
    request = urllib.request.Request(url, headers={"User-Agent": "ExonDomainCompare/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8", errors="replace"))


def query_ncbi_gene(gene_symbol: str, scientific_name: str = "", taxid: str = "",
                    *, timeout: float = 30.0, sleep_s: float = 0.34,
                    retmax: int = 25) -> Tuple[List[NcbiGeneRecord], str]:
    """Species-qualified NCBI Gene lookup for a symbol, and how the lookup went.

    The query is qualified by organism so a symbol cannot resolve against another species.
    Both the official symbol and the source-provided aliases are returned, because the
    symbol the user knows is frequently only an alias — that is the whole point.
    """
    organism = (scientific_name or "").strip()
    if not organism and taxid:
        organism = f"txid{taxid}"
    if not organism:
        return [], "no_species_qualifier"
    term = f'"{organism}"[Organism] AND {gene_symbol}[All Fields]'
    try:
        search = _fetch_json(
            EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(
                {"db": "gene", "term": term, "retmode": "json", "retmax": str(retmax)}),
            timeout, 0.0)
        ids = ((search.get("esearchresult") or {}).get("idlist") or [])
        if not ids:
            return [], "no_hit"
        summary = _fetch_json(
            EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(
                {"db": "gene", "id": ",".join(ids), "retmode": "json"}),
            timeout, sleep_s)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError,
            ValueError) as exc:
        return [], f"failed:{type(exc).__name__}"

    result = summary.get("result") or {}
    records: List[NcbiGeneRecord] = []
    for uid in ids:
        row = result.get(uid)
        if not isinstance(row, dict):
            continue
        genomic = (row.get("genomicinfo") or [{}])
        first = genomic[0] if genomic else {}
        exon_count = first.get("exoncount")
        records.append(NcbiGeneRecord(
            gene_id=str(row.get("uid") or uid),
            official_symbol=str(row.get("name") or ""),
            description=str(row.get("description") or ""),
            aliases=[a.strip() for a in str(row.get("otheraliases") or "").split(",")
                     if a.strip()],
            other_designations=[d.strip() for d in
                                str(row.get("otherdesignations") or "").split("|")
                                if d.strip()],
            chromosome=str(row.get("chromosome") or ""),
            chr_accession=str(first.get("chraccver") or ""),
            exon_count=int(exon_count) if str(exon_count or "").isdigit() else None,
            status=str(row.get("status") or ""),
            current_id=str(row.get("currentid") or "")))
    return records, "ok"


#: Injection point. Tests and offline runs replace this with a fixture, so the cascade is
#: exercised without a network and the production algorithm stays identical.
GeneLookup = Callable[[str, str, str], Tuple[List[NcbiGeneRecord], str]]


# --------------------------------------------------------------------------- #
# Candidate inventory
# --------------------------------------------------------------------------- #
@dataclass
class Candidate:
    """One plausible locus and every reason for and against it."""

    locus: AnnotationLocus
    routes: List[str] = field(default_factory=list)
    ncbi_gene_id: str = ""
    ncbi_official_symbol: str = ""
    ncbi_aliases: List[str] = field(default_factory=list)
    ncbi_description: str = ""
    orthology_evidence: str = ""
    family_evidence: str = ""
    decision: str = ""
    reason: str = ""

    @property
    def decisive_routes(self) -> List[str]:
        return [r for r in self.routes if r not in SUPPORTING_ONLY]

    def as_dict(self) -> Dict[str, Any]:
        return {
            **self.locus.as_dict(),
            "routes": list(self.routes),
            "decisive_routes": self.decisive_routes,
            "ncbi_gene_id": self.ncbi_gene_id,
            "ncbi_official_symbol": self.ncbi_official_symbol,
            "ncbi_aliases": list(self.ncbi_aliases),
            "ncbi_description": self.ncbi_description,
            "orthology_evidence": self.orthology_evidence,
            "family_evidence": self.family_evidence,
            "decision": self.decision, "reason": self.reason,
        }


@dataclass
class GeneIdentity:
    """What the user asked for and what the source calls it — kept apart on purpose.

    Overwriting the requested symbol with ``LOC122209636`` would hide the question the user
    asked; hiding the source symbol would make the result unverifiable against NCBI. Both
    are recorded, and the display symbol stays the user's.
    """

    requested_gene_symbol: str
    resolved_gene_id: str = ""
    resolved_official_symbol: str = ""
    resolved_display_symbol: str = ""
    source_description: str = ""
    resolution_method: str = ""
    resolution_confidence: str = ""
    source_provenance: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "requested_gene_symbol": self.requested_gene_symbol,
            "resolved_gene_id": self.resolved_gene_id,
            "resolved_official_symbol": self.resolved_official_symbol,
            "resolved_display_symbol": (self.resolved_display_symbol
                                        or self.requested_gene_symbol),
            "source_description": self.source_description,
            "resolution_method": self.resolution_method,
            "resolution_confidence": self.resolution_confidence,
            "source_provenance": dict(self.source_provenance),
            # True when the annotation's own symbol differs from what the user typed, so
            # the UI knows to show both rather than silently swapping one for the other.
            "symbol_differs_from_source": bool(
                self.resolved_official_symbol
                and self.resolved_official_symbol.upper()
                != (self.requested_gene_symbol or "").upper()),
        }


@dataclass
class Resolution:
    status: str
    identity: GeneIdentity
    locus: Optional[AnnotationLocus] = None
    candidates: List[Candidate] = field(default_factory=list)
    ncbi_records: List[NcbiGeneRecord] = field(default_factory=list)
    ncbi_lookup_status: str = "not_attempted"
    detail: str = ""
    routes_attempted: List[str] = field(default_factory=list)

    @property
    def resolved(self) -> bool:
        return self.status == RESOLVED and self.locus is not None

    def message(self) -> str:
        template = STATUS_MESSAGE.get(self.status, "")
        if not template:
            return self.detail
        return template.format(
            gene=self.identity.requested_gene_symbol,
            official=(self.identity.resolved_official_symbol or "an assembly locus"),
            candidates=", ".join(c.locus.symbol for c in self.candidates
                                 if c.decision == "candidate") or "several loci",
            detail=self.detail or "no further detail recorded")

    def as_dict(self) -> Dict[str, Any]:
        return {
            "status": self.status,
            "message": self.message(),
            "detail": self.detail,
            "identity": self.identity.as_dict(),
            "locus": self.locus.as_dict() if self.locus else None,
            "candidates": [c.as_dict() for c in self.candidates],
            "ncbi_records": [r.as_dict() for r in self.ncbi_records],
            "ncbi_lookup_status": self.ncbi_lookup_status,
            "routes_attempted": list(self.routes_attempted),
        }


# --------------------------------------------------------------------------- #
# The cascade
# --------------------------------------------------------------------------- #
def _describes(text: str, needle: str) -> bool:
    return bool(needle) and needle.lower() in (text or "").lower()


def _family_stem(description: str) -> str:
    """The family part of a product name, with a trailing '-like' or variant dropped."""
    text = re.sub(r"\s*,?\s*transcript variant.*$", "", description or "",
                  flags=re.IGNORECASE)
    text = re.sub(r"[\s-]*like$", "", text.strip(), flags=re.IGNORECASE)
    return text.strip().lower()


def resolve_gene_locus(gff_path: Path, gene_symbol: str, *,
                       scientific_name: str = "", taxid: str = "",
                       gene_lookup: Optional[GeneLookup] = None,
                       allow_network: bool = True,
                       assembly_accession: str = "") -> Resolution:
    """Resolve ``gene_symbol`` to one locus of ``gff_path``, or say why not.

    The cascade is ordered by how much the evidence can be trusted: an identifier beats a
    symbol, a symbol beats a name similarity, and a name similarity never decides alone.
    """
    requested = (gene_symbol or "").strip()
    identity = GeneIdentity(requested_gene_symbol=requested,
                            resolved_display_symbol=requested)
    identity.source_provenance = {"assembly_accession": assembly_accession,
                                  "annotation_file": str(gff_path),
                                  "scientific_name": scientific_name, "taxid": taxid}
    attempted: List[str] = []

    # A parser failure has to stay a parser failure. Turning it into "gene not found" is
    # what let a broken annotation read as an absent gene.
    try:
        loci = scan_annotation(Path(gff_path))
    except OSError as exc:
        return Resolution(PARSER_FAILED, identity, detail=f"{type(exc).__name__}: {exc}",
                          routes_attempted=attempted)
    except Exception as exc:  # malformed coordinates, broken attributes
        return Resolution(PARSER_FAILED, identity,
                          detail=f"{type(exc).__name__}: {exc}", routes_attempted=attempted)
    if not loci:
        return Resolution(PARSER_FAILED, identity,
                          detail="the annotation contains no gene feature",
                          routes_attempted=attempted)

    want = requested.upper()
    by_symbol: Dict[str, List[AnnotationLocus]] = {}
    by_geneid: Dict[str, AnnotationLocus] = {}
    for locus in loci.values():
        by_symbol.setdefault((locus.symbol or "").upper(), []).append(locus)
        if locus.source_gene_id:
            by_geneid.setdefault(locus.source_gene_id, locus)

    def accept(locus: AnnotationLocus, route: str, confidence: str,
               record: Optional[NcbiGeneRecord] = None,
               extra_candidates: Sequence[Candidate] = (),
               ncbi_status: str = "not_attempted",
               records: Sequence[NcbiGeneRecord] = ()) -> Resolution:
        identity.resolved_official_symbol = locus.symbol
        identity.resolution_method = route
        identity.resolution_confidence = confidence
        identity.resolved_gene_id = (record.gene_id if record
                                     else locus.source_gene_id)
        identity.source_description = (
            (record.description if record else "")
            or locus.description
            or (locus.products[0] if locus.products else ""))
        chosen = Candidate(locus=locus, routes=[route], decision="accepted",
                           reason=f"selected through {route}")
        if record:
            chosen.ncbi_gene_id = record.gene_id
            chosen.ncbi_official_symbol = record.official_symbol
            chosen.ncbi_aliases = list(record.aliases)
            chosen.ncbi_description = record.description
        return Resolution(RESOLVED, identity, locus=locus,
                          candidates=[chosen, *extra_candidates],
                          ncbi_records=list(records), ncbi_lookup_status=ncbi_status,
                          routes_attempted=attempted)

    # 1 + 2. The annotation's own symbol. Exact case first, then case-normalized: these
    # differ under teleost and plant conventions (foxp1b vs FOXP1B).
    #
    # Duplicates are ambiguity, not a race. Returning whichever row came first in the file
    # would make the answer depend on annotation order.
    attempted.append(ROUTE_EXACT_SYMBOL)
    case_exact = [l for l in loci.values() if l.symbol == requested]
    if len(case_exact) == 1:
        return accept(case_exact[0], ROUTE_EXACT_SYMBOL, "exact_symbol")
    attempted.append(ROUTE_NORMALIZED_SYMBOL)
    exact = case_exact or by_symbol.get(want, [])
    route = ROUTE_EXACT_SYMBOL if case_exact else ROUTE_NORMALIZED_SYMBOL
    if len(exact) == 1:
        return accept(exact[0], route, "exact_symbol")
    if len(exact) > 1:
        cands = [Candidate(locus=l, routes=[route], decision="candidate",
                           reason="several loci carry this exact symbol") for l in exact]
        return Resolution(AMBIGUOUS_FAMILY, identity, candidates=cands,
                          detail=f"{len(exact)} loci carry the symbol {requested}",
                          routes_attempted=attempted)

    # 3. An alias on the annotation's gene row. Present on curated annotations; entirely
    # absent from Gnomon-annotated assemblies, which is why the cascade cannot stop here.
    attempted.append(ROUTE_ANNOTATION_ALIAS)
    alias_hits = [l for l in loci.values()
                  if any(a.upper() == want for a in l.aliases)]
    unique_alias = {l.symbol.upper(): l for l in alias_hits}
    if len(unique_alias) == 1:
        return accept(next(iter(unique_alias.values())), ROUTE_ANNOTATION_ALIAS,
                      "annotation_alias")
    if len(unique_alias) > 1:
        cands = [Candidate(locus=l, routes=[ROUTE_ANNOTATION_ALIAS],
                           decision="candidate",
                           reason="several loci list this symbol as an alias")
                 for l in unique_alias.values()]
        return Resolution(AMBIGUOUS_FAMILY, identity, candidates=cands,
                          detail=f"{len(unique_alias)} loci list {requested} as an alias",
                          routes_attempted=attempted)

    # 4-6. NCBI Gene. The symbol the user knows may exist only here, as an alias, and the
    # link back into the assembly is the GeneID.
    lookup = gene_lookup
    if lookup is None and allow_network:
        lookup = lambda sym, name, tid: query_ncbi_gene(sym, name, tid)  # noqa: E731
    records: List[NcbiGeneRecord] = []
    lookup_status = "not_attempted"
    if lookup is not None:
        attempted.append(ROUTE_NCBI_SYMBOL)
        records, lookup_status = lookup(requested, scientific_name, taxid)

    if lookup_status.startswith("failed"):
        # Not knowing is not the same as knowing there is nothing.
        return Resolution(SOURCE_UNAVAILABLE, identity,
                          ncbi_lookup_status=lookup_status,
                          detail=f"NCBI Gene lookup {lookup_status}",
                          routes_attempted=attempted)

    live = [r for r in records if r.is_live and r.mentions(requested)]
    if live:
        attempted.append(ROUTE_NCBI_ALIAS_GENEID)
        # Only records that actually name the requested symbol count, and only those whose
        # GeneID or official symbol is present in this annotation are usable.
        mapped: List[Tuple[NcbiGeneRecord, AnnotationLocus, str]] = []
        for record in live:
            locus = by_geneid.get(record.gene_id)
            if locus is not None:
                mapped.append((record, locus, ROUTE_NCBI_ALIAS_GENEID))
                continue
            by_official = by_symbol.get((record.official_symbol or "").upper(), [])
            if len(by_official) == 1:
                mapped.append((record, by_official[0], ROUTE_NCBI_ALIAS_SYMBOL))
        if ROUTE_NCBI_ALIAS_SYMBOL not in attempted:
            attempted.append(ROUTE_NCBI_ALIAS_SYMBOL)

        distinct = {locus.gene_id: (record, locus, route)
                    for record, locus, route in mapped}
        if len(distinct) == 1:
            record, locus, route = next(iter(distinct.values()))
            # The description corroborates but did not decide; record that it agreed.
            supporting = []
            if _describes(" ".join(locus.products), _family_stem(record.description)):
                supporting.append(ROUTE_DESCRIPTION)
            resolution = accept(locus, route,
                                "ncbi_geneid_mapped" if route == ROUTE_NCBI_ALIAS_GENEID
                                else "ncbi_official_symbol_mapped",
                                record=record,
                                extra_candidates=_family_inventory(
                                    loci, record, exclude=locus.gene_id),
                                ncbi_status=lookup_status, records=records)
            resolution.candidates[0].routes.extend(supporting)
            resolution.candidates[0].family_evidence = record.description
            return resolution
        if len(distinct) > 1:
            cands = [Candidate(locus=locus, routes=[route], decision="candidate",
                               ncbi_gene_id=record.gene_id,
                               ncbi_official_symbol=record.official_symbol,
                               ncbi_aliases=list(record.aliases),
                               ncbi_description=record.description,
                               reason="several current NCBI Gene records name this symbol")
                     for record, locus, route in distinct.values()]
            return Resolution(AMBIGUOUS_FAMILY, identity, candidates=cands,
                              ncbi_records=records, ncbi_lookup_status=lookup_status,
                              detail=(f"{len(distinct)} current NCBI Gene records name "
                                      f"{requested} and all map into this annotation"),
                              routes_attempted=attempted)

        # NCBI knows the symbol; the annotation does not contain the locus. That is a
        # mapping failure, and saying "gene not found" here would blame the species.
        identity.resolved_gene_id = live[0].gene_id
        identity.resolved_official_symbol = live[0].official_symbol
        identity.source_description = live[0].description
        identity.resolution_method = ROUTE_NCBI_ALIAS_GENEID
        identity.resolution_confidence = "unmapped"
        return Resolution(ALIAS_MAPPING_FAILED, identity,
                          candidates=_family_inventory(loci, live[0]),
                          ncbi_records=records, ncbi_lookup_status=lookup_status,
                          detail=(f"NCBI Gene {live[0].gene_id} "
                                  f"({live[0].official_symbol}) is absent from "
                                  f"{assembly_accession or 'the selected assembly'}"),
                          routes_attempted=attempted)

    # 7. Replaced or discontinued records only — say so rather than using them.
    superseded = [r for r in records if not r.is_live]
    if superseded and not live:
        return Resolution(REVIEW_REQUIRED, identity, ncbi_records=records,
                          ncbi_lookup_status=lookup_status,
                          detail=("the only NCBI Gene records for this symbol are "
                                  "replaced or discontinued: "
                                  + ", ".join(f"{r.gene_id}({r.status or 'replaced'})"
                                              for r in superseded)),
                          routes_attempted=attempted)

    return Resolution(GENE_NOT_FOUND, identity, ncbi_records=records,
                      ncbi_lookup_status=lookup_status,
                      detail=(f"no locus in the annotation and no current NCBI Gene "
                              f"record names {requested}"),
                      routes_attempted=attempted)


def _family_inventory(loci: Dict[str, AnnotationLocus], record: NcbiGeneRecord,
                      exclude: str = "") -> List[Candidate]:
    """The related loci that were considered and why each was not chosen.

    Recorded even on success: the reason the lion HBA locus is defensible is that the four
    other alpha-like loci on that chromosome were seen and rejected on identifier grounds,
    not overlooked.
    """
    stem = _family_stem(record.description)
    if not stem:
        return []
    out: List[Candidate] = []
    for locus in loci.values():
        if locus.gene_id == exclude:
            continue
        text = " ".join([locus.description, *locus.products])
        if not _describes(text, stem):
            continue
        candidate = Candidate(locus=locus, routes=[ROUTE_DESCRIPTION],
                              family_evidence=text[:200])
        if locus.is_pseudogene:
            candidate.decision = "rejected"
            candidate.reason = "pseudogene"
        elif locus.protein_count == 0:
            candidate.decision = "rejected"
            candidate.reason = "no translated CDS"
        else:
            candidate.decision = "rejected"
            candidate.reason = ("related family member; description similarity is not "
                                "decisive and its GeneID is not the resolved one")
        out.append(candidate)
    return sorted(out, key=lambda c: c.locus.symbol)[:25]


__all__ = [
    "ROUTE_EXACT_SYMBOL", "ROUTE_NORMALIZED_SYMBOL", "ROUTE_ANNOTATION_ALIAS",
    "ROUTE_NCBI_SYMBOL", "ROUTE_NCBI_ALIAS_GENEID", "ROUTE_NCBI_ALIAS_SYMBOL",
    "ROUTE_DESCRIPTION", "ROUTE_ORTHOLOGY", "ROUTE_ORDER", "SUPPORTING_ONLY",
    "RESOLVED", "GENE_NOT_FOUND", "ALIAS_MAPPING_FAILED", "AMBIGUOUS_FAMILY",
    "TRANSCRIPT_NOT_FOUND", "NO_VALID_TRANSLATED_CDS", "PARSER_FAILED",
    "SOURCE_UNAVAILABLE", "REVIEW_REQUIRED", "STATUSES", "STATUS_MESSAGE",
    "AnnotationLocus", "NcbiGeneRecord", "Candidate", "GeneIdentity", "Resolution",
    "scan_annotation", "query_ncbi_gene", "resolve_gene_locus",
]
