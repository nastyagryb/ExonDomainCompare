"""gene_identification.py — find the requested gene in an annotation, cautiously.

The old matcher accepted a gene locus when any of ``gene``, ``Name``,
``gene_synonym``, ``gene_name``, ``product`` or ``description`` held a token equal to
``FGFR2``. That is too strict in one direction and too loose in the other.

Too strict: a genome annotated before the symbol was assigned calls the locus
``LOC124236178`` and puts "fibroblast growth factor receptor 2" only in the
description. No token equals FGFR2, so the gene is invisible — even though it is the
right gene, with transcripts and proteins.

Too loose: it will accept a locus whose ``gene_synonym`` mentions FGFR2 while the
locus itself is FGFR1, and it takes whatever matched first. FGFR1, FGFR3 and FGFR4 are
close paralogs with the same domain architecture; confusing one for FGFR2 would not
look like an error downstream. It would look like a result.

So candidates are gathered through ordered routes, weakest evidence last, and each one
carries the evidence that admitted it and the reason it was accepted or rejected. A
candidate reached by a weak route must clear paralog discrimination before it can be
used; a LOC-labelled locus is never accepted merely for being the first hit.

Nothing here knows about any particular species.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

#: Routes, strongest first. The route that admitted a candidate is recorded, so a
#: reviewer can see whether a model rests on an official symbol or on a description.
ROUTE_EXACT_SYMBOL = "exact_gene_symbol"
ROUTE_GENE_ID = "official_gene_identifier"
ROUTE_SYNONYM = "source_provided_synonym"
ROUTE_DESCRIPTION = "orthology_aware_description"
ROUTE_LOC = "loc_labelled_candidate"
ROUTE_SEQUENCE = "sequence_similarity_rescue"

ROUTE_ORDER = (ROUTE_EXACT_SYMBOL, ROUTE_GENE_ID, ROUTE_SYNONYM,
               ROUTE_DESCRIPTION, ROUTE_LOC, ROUTE_SEQUENCE)

#: Routes whose evidence does not by itself identify which paralog was found. A
#: candidate admitted through one of these needs positive discrimination before use.
ROUTES_NEEDING_DISCRIMINATION = (ROUTE_DESCRIPTION, ROUTE_LOC, ROUTE_SEQUENCE)

#: Outcomes for the identification as a whole.
FOUND = "gene_identified"
NOT_FOUND = "no_gene_candidate_in_annotation"
AMBIGUOUS = "ambiguous_paralog_candidates"
REJECTED_PARALOG = "candidates_rejected_as_paralogs"

_SYMBOL_FIELDS = ("gene", "Name", "gene_name", "locus_tag")
_SYNONYM_FIELDS = ("gene_synonym", "synonym", "Alias")
_TEXT_FIELDS = ("product", "description", "Note", "gene_desc")

_LOC_RE = re.compile(r"^LOC\d+$", re.IGNORECASE)


def _tokens(value: Any) -> List[str]:
    return [t.strip() for t in re.split(r"[,|;]", str(value or "")) if t.strip()]


def _norm(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip()).lower()


def paralog_family(symbol: str) -> Tuple[str, str]:
    """Split a symbol like ``FGFR3`` into family and member number.

    Used only to tell members of one family apart. It deliberately does not try to
    interpret arbitrary gene names.
    """
    match = re.fullmatch(r"([A-Za-z]+?)(\d+)([A-Za-z0-9\-]*)", str(symbol or "").strip())
    if not match:
        return (str(symbol or "").strip().upper(), "")
    return (match.group(1).upper(), match.group(2))


def sibling_symbols(gene_symbol: str, members: Sequence[str] = ("1", "2", "3", "4")
                    ) -> List[str]:
    """The paralogs a candidate must be told apart from.

    For FGFR2 this is FGFR1, FGFR3 and FGFR4 — derived from the requested symbol
    rather than listed, so the same logic protects any numbered family.

    A symbol whose number is not one of the enumerated members has no such family:
    TP53 is not the fifty-third of a series, and inventing TP1…TP4 would have the
    discriminator reject loci for resembling genes that do not exist.
    """
    family, number = paralog_family(gene_symbol)
    if not number or number not in members:
        return []
    return [f"{family}{m}" for m in members if m != number]


def _product_variant(product_name: str, number: str) -> str:
    """A sibling's product name: the same wording with the trailing number swapped.

    "fibroblast growth factor receptor 2" → "fibroblast growth factor receptor 3".
    Used to recognise a paralog's description, which is otherwise identical to this
    gene's down to the last word.
    """
    return re.sub(r"\d+\s*$", number, _norm(product_name)).strip()


def _describes_product(text: str, gene_symbol: str,
                       product_name: str = "") -> bool:
    """Whether free text names this gene's product rather than a sibling's.

    Two ways to qualify: the symbol itself appears (``FGFR2``), or the source's own
    spelled-out product name appears (``fibroblast growth factor receptor 2``).

    There is deliberately no "…receptor 2" shortcut. That rule looked harmless and
    admitted sixty-six loci from one zebra annotation — ADIPOR2, NPR2, GRM2, TGFBR2,
    RYR2, TACR2 — every gene whose description happens to end in the number two. A
    route that matches most of the genome is not evidence of anything.
    """
    family, number = paralog_family(gene_symbol)
    if not number:
        return False
    lowered = _norm(text)
    if not lowered:
        return False
    if re.search(rf"\b{re.escape(family.lower())}\s*-?\s*{number}\b", lowered):
        return True
    if product_name and _product_variant(product_name, number) in lowered:
        return True
    return False


def _mentions_sibling(text: str, siblings: Sequence[str],
                      product_name: str = "") -> Optional[str]:
    lowered = _norm(text)
    for sib in siblings:
        family, number = paralog_family(sib)
        if re.search(rf"\b{re.escape(family.lower())}\s*-?\s*{number}\b", lowered):
            return sib
        if product_name and _product_variant(product_name, number) in lowered:
            return sib
    return None


@dataclass
class GeneCandidate:
    """One annotated locus considered as the requested gene."""

    source_gene_id: str = ""
    source_symbol: str = ""
    description: str = ""
    seqid: str = ""
    strand: str = ""
    start: str = ""
    end: str = ""
    biotype: str = ""
    dbxrefs: List[str] = field(default_factory=list)
    synonyms: List[str] = field(default_factory=list)
    transcript_count: int = 0
    protein_count: int = 0
    route: str = ""
    orthology_evidence: str = ""
    similarity_evidence: str = ""
    paralog_evidence: str = ""
    decision: str = ""
    reason: str = ""

    def as_row(self) -> Dict[str, str]:
        return {
            "source_gene_id": self.source_gene_id,
            "source_symbol": self.source_symbol,
            "description": self.description,
            "seqid": self.seqid,
            "strand": self.strand,
            "start": self.start,
            "end": self.end,
            "biotype": self.biotype,
            "dbxrefs": ";".join(self.dbxrefs),
            "synonyms": ";".join(self.synonyms),
            "transcript_count": str(self.transcript_count),
            "protein_count": str(self.protein_count),
            "identification_route": self.route,
            "orthology_evidence": self.orthology_evidence,
            "similarity_evidence": self.similarity_evidence,
            "paralog_discrimination_evidence": self.paralog_evidence,
            "decision": self.decision,
            "reason": self.reason,
        }


def candidate_from_attributes(attr: Dict[str, str], *, seqid: str = "", strand: str = "",
                              start: str = "", end: str = "") -> GeneCandidate:
    symbol = ""
    for key in _SYMBOL_FIELDS:
        if attr.get(key):
            symbol = _tokens(attr[key])[0] if _tokens(attr[key]) else ""
            break
    synonyms: List[str] = []
    for key in _SYNONYM_FIELDS:
        synonyms += _tokens(attr.get(key, ""))
    description = ""
    for key in _TEXT_FIELDS:
        if attr.get(key):
            description = str(attr[key])
            break
    return GeneCandidate(
        source_gene_id=attr.get("ID", "") or attr.get("Dbxref", ""),
        source_symbol=symbol,
        description=description,
        seqid=seqid, strand=strand, start=start, end=end,
        biotype=attr.get("gene_biotype", "") or attr.get("gbkey", ""),
        dbxrefs=_tokens(attr.get("Dbxref", "")),
        synonyms=sorted(dict.fromkeys(synonyms)),
    )


def classify_route(candidate: GeneCandidate, gene_symbol: str,
                   expected_gene_ids: Sequence[str] = (),
                   product_name: str = "") -> Tuple[str, str]:
    """The strongest route that admits this candidate, and the evidence for it.

    Returns ``("", "")`` when nothing about the locus points at the requested gene.
    """
    target = _norm(gene_symbol)

    if _norm(candidate.source_symbol) == target:
        return ROUTE_EXACT_SYMBOL, f"gene symbol {candidate.source_symbol}"

    wanted = {str(g).strip() for g in expected_gene_ids if str(g).strip()}
    if wanted:
        for xref in candidate.dbxrefs:
            # A GeneID cross-reference is an identifier, not a name, so it survives a
            # locus that has no symbol yet — the LOC case.
            if ":" in xref and xref.split(":", 1)[1].strip() in wanted:
                return ROUTE_GENE_ID, f"Dbxref {xref}"
        if candidate.source_gene_id.replace("gene-", "").strip() in wanted:
            return ROUTE_GENE_ID, f"gene id {candidate.source_gene_id}"

    if any(_norm(s) == target for s in candidate.synonyms):
        return ROUTE_SYNONYM, f"source synonym list {';'.join(candidate.synonyms)}"

    if _describes_product(candidate.description, gene_symbol, product_name):
        return ROUTE_DESCRIPTION, f"description {candidate.description!r}"

    if _LOC_RE.match(candidate.source_symbol or "") and _describes_product(
            f"{candidate.description} {' '.join(candidate.synonyms)}", gene_symbol,
            product_name):
        return ROUTE_LOC, (f"LOC-labelled locus {candidate.source_symbol} described as "
                           f"{candidate.description!r}")
    return "", ""


def discriminate(candidate: GeneCandidate, gene_symbol: str,
                 product_name: str = "") -> Tuple[bool, str]:
    """Whether this candidate is the requested paralog rather than a sibling.

    Annotation-level discrimination only: text that names a sibling is disqualifying,
    an exact symbol or an official identifier is conclusive. Anything weaker is left
    undecided here and handed to the sequence step, because guessing between FGFR2 and
    FGFR3 from a description is exactly the error that would not show up later.
    """
    siblings = sibling_symbols(gene_symbol)
    haystack = " ".join([candidate.source_symbol, candidate.description,
                         *candidate.synonyms])

    if _norm(candidate.source_symbol) in {_norm(s) for s in siblings}:
        sib = candidate.source_symbol
        return False, f"locus is annotated as the paralog {sib}"

    if candidate.route in (ROUTE_EXACT_SYMBOL, ROUTE_GENE_ID, ROUTE_SYNONYM):
        return True, f"identified by {candidate.route}"

    sibling = _mentions_sibling(haystack, siblings, product_name)
    if sibling and not _describes_product(haystack, gene_symbol, product_name):
        return False, f"annotation text names the paralog {sibling}"
    if sibling:
        return False, (f"annotation text names both {gene_symbol} and the paralog "
                       f"{sibling}; the locus cannot be assigned from annotation alone")
    return False, "weak route: sequence discrimination required"


def best_paralog_by_similarity(protein: str, panel: Dict[str, str]
                               ) -> Tuple[str, float, Dict[str, float]]:
    """Which panel member a translated protein most resembles.

    A k-mer containment score, not an alignment. It is used for one decision only —
    which of a handful of paralogs a sequence is closest to — and for that a cheap
    symmetric measure is enough and needs no external aligner. The absolute value is
    not interpretable as identity, so it is never reported as one.
    """
    def kmers(seq: str, k: int = 5) -> set:
        seq = re.sub(r"[^A-Za-z]", "", seq or "").upper()
        return {seq[i:i + k] for i in range(max(0, len(seq) - k + 1))}

    query = kmers(protein)
    scores: Dict[str, float] = {}
    if not query:
        return "", 0.0, scores
    for name, seq in panel.items():
        ref = kmers(seq)
        if not ref:
            continue
        scores[name] = len(query & ref) / len(query | ref)
    if not scores:
        return "", 0.0, scores
    best = max(scores, key=scores.get)
    return best, scores[best], scores


def read_panel(path) -> Dict[str, str]:
    """Read a FASTA paralog panel into ``{header: sequence}``."""
    panel: Dict[str, str] = {}
    name = ""
    chunks: List[str] = []
    try:
        text = open(path, "r", encoding="utf-8").read()
    except OSError:
        return panel
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                panel[name] = "".join(chunks)
            name = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())
    if name:
        panel[name] = "".join(chunks)
    return panel


def panel_member_symbol(header: str, gene_symbol: str) -> str:
    """The symbol a panel header stands for, e.g. ``human_FGFR3_UniProt_...`` → FGFR3."""
    family, _ = paralog_family(gene_symbol)
    match = re.search(rf"{re.escape(family)}\s*-?\s*(\d+)", header, re.IGNORECASE)
    return f"{family}{match.group(1)}" if match else ""


def discriminate_by_sequence(candidate: GeneCandidate, gene_symbol: str,
                             protein: str, panel: Dict[str, str],
                             margin: float = 0.02) -> Tuple[bool, str]:
    """Accept a weak-route candidate only if its protein is closest to this paralog.

    The margin keeps a near-tie from being read as a decision. Two paralogs that score
    within it are reported as ambiguous rather than resolved, because the run is better
    off saying it could not tell than naming the wrong receptor.
    """
    if not protein:
        return False, "no translated protein available for sequence discrimination"
    if not panel:
        return False, "no paralog reference panel available"

    best, _score, scores = best_paralog_by_similarity(protein, panel)
    best_symbol = panel_member_symbol(best, gene_symbol)
    ranked = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
    evidence = ", ".join(f"{panel_member_symbol(h, gene_symbol) or h}={s:.3f}"
                         for h, s in ranked[:4])

    if best_symbol.upper() != gene_symbol.upper():
        return False, (f"translated protein is closest to {best_symbol or best} "
                       f"({evidence}); rejected as a paralog")
    if len(ranked) > 1 and (ranked[0][1] - ranked[1][1]) < margin:
        runner = panel_member_symbol(ranked[1][0], gene_symbol) or ranked[1][0]
        return False, (f"translated protein scores {ranked[0][1]:.3f} for {gene_symbol} "
                       f"and {ranked[1][1]:.3f} for {runner}, within the {margin} "
                       f"decision margin ({evidence}); left ambiguous")
    return True, (f"translated protein is closest to {gene_symbol} by k-mer "
                  f"containment ({evidence})")


@dataclass
class Identification:
    status: str
    accepted: Optional[GeneCandidate] = None
    candidates: List[GeneCandidate] = field(default_factory=list)
    detail: str = ""

    def rows(self) -> List[Dict[str, str]]:
        return [c.as_row() for c in self.candidates]


def identify(candidates: Sequence[GeneCandidate], gene_symbol: str, *,
             expected_gene_ids: Sequence[str] = (),
             proteins: Optional[Dict[str, str]] = None,
             panel: Optional[Dict[str, str]] = None,
             product_name: str = "") -> Identification:
    """Choose the locus that is the requested gene, recording every decision.

    ``proteins`` maps a candidate's ``source_gene_id`` to a translated sequence, used
    only for candidates whose annotation evidence is too weak to decide.

    ``product_name`` is the source's spelled-out name for this gene's product, e.g.
    "fibroblast growth factor receptor 2". It is what makes a description route
    specific: without it, only the symbol itself counts as description evidence.
    """
    considered: List[GeneCandidate] = []
    for cand in candidates:
        route, evidence = classify_route(cand, gene_symbol, expected_gene_ids,
                                         product_name)
        if not route:
            continue
        cand.route = route
        cand.orthology_evidence = evidence
        considered.append(cand)

    if not considered:
        return Identification(
            status=NOT_FOUND,
            detail=(f"no locus in the annotation carries the symbol {gene_symbol}, an "
                    "official identifier for it, a source synonym, or a description of "
                    "its product"))

    considered.sort(key=lambda c: ROUTE_ORDER.index(c.route))

    accepted: Optional[GeneCandidate] = None
    ambiguous: List[GeneCandidate] = []
    for cand in considered:
        ok, reason = discriminate(cand, gene_symbol, product_name)
        cand.paralog_evidence = reason
        if ok:
            if accepted is None:
                accepted, cand.decision, cand.reason = cand, "accepted", reason
            else:
                cand.decision = "rejected"
                cand.reason = (f"a stronger candidate was already accepted "
                               f"({accepted.source_symbol or accepted.source_gene_id} "
                               f"via {accepted.route})")
            continue

        if reason.startswith("weak route"):
            sequence = (proteins or {}).get(cand.source_gene_id, "")
            ok_seq, seq_reason = discriminate_by_sequence(
                cand, gene_symbol, sequence, panel or {})
            cand.similarity_evidence = seq_reason
            cand.route = ROUTE_SEQUENCE if ok_seq else cand.route
            if ok_seq and accepted is None:
                accepted, cand.decision, cand.reason = cand, "accepted", seq_reason
            elif ok_seq:
                cand.decision = "rejected"
                cand.reason = "a stronger candidate was already accepted"
            else:
                cand.decision = "rejected"
                cand.reason = seq_reason
                if "ambiguous" in seq_reason or "no translated protein" in seq_reason:
                    ambiguous.append(cand)
            continue

        cand.decision = "rejected"
        cand.reason = reason

    if accepted is not None:
        return Identification(
            status=FOUND, accepted=accepted, candidates=considered,
            detail=(f"{accepted.source_symbol or accepted.source_gene_id} on "
                    f"{accepted.seqid} accepted via {accepted.route}: {accepted.reason}"))

    if ambiguous:
        return Identification(
            status=AMBIGUOUS, candidates=considered,
            detail=("candidates were found but none could be assigned to "
                    f"{gene_symbol} rather than a paralog: "
                    + "; ".join(f"{c.source_symbol or c.source_gene_id}: {c.reason}"
                                for c in ambiguous)))

    return Identification(
        status=REJECTED_PARALOG, candidates=considered,
        detail=(f"every candidate was rejected: "
                + "; ".join(f"{c.source_symbol or c.source_gene_id}: {c.reason}"
                            for c in considered)))
