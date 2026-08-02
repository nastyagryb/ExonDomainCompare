"""taxon_resolution.py — resolve a submitted species name to a real NCBI taxon.

The species registry used to answer this question from a literal table of thirty
entries: the validated FGFR2 panel. For anything else it copied the string it was
given into the field later used as an NCBI query term. Since the web form submits
Ensembl-style slugs, that field became ``equus_quagga``, and

    datasets summary genome taxon equus_quagga
    Error: The taxonomy name 'equus_quagga' is not recognized.

which is where the Equus quagga run died — four stages before the empty transcript
table that was reported as its cause. The failure was certain for every species
outside the table and invisible for every species inside it.

So resolution is asked of the service that owns the answer. What comes back is the
accepted scientific name, the numeric taxid, the rank and the synonyms the source
itself publishes. The numeric taxid is the important one: it is unambiguous, it needs
no spelling, and every downstream query prefers it.

Two guarantees this module owes the caller:

*Never substitute a different species.* A name that cannot be resolved is returned
unresolved, with a reason. Silently answering an *Equus quagga* request with *Equus
caballus* — same genus, a published genome, a plausible FGFR2 — would produce a run
that looks successful and is about the wrong animal.

*A synonym is only an alias if it resolves back.* Querying under a synonym is allowed
when the source maps that synonym to the requested taxon; it is not allowed as a way
of reaching a neighbouring taxon that happens to have better data.
"""
from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

#: Resolution outcomes. A caller must be able to tell "this is not a species" from
#: "the service was unreachable"; the first is a data problem the user can fix by
#: correcting the name, the second is a transient one that a retry may fix.
RESOLVED = "resolved"
RESOLVED_VIA_SYNONYM = "resolved_via_synonym"
NOT_FOUND = "taxon_not_found"
AMBIGUOUS = "taxon_ambiguous"
SERVICE_UNAVAILABLE = "taxonomy_service_unavailable"
OFFLINE = "taxonomy_lookup_skipped_offline"


@dataclass
class TaxonIdentity:
    """One resolved (or unresolved) species identity, as recorded in the registry."""

    submitted_name: str
    species_id: str = ""
    status: str = NOT_FOUND
    accepted_name: str = ""
    taxid: str = ""
    rank: str = ""
    common_name: str = ""
    synonyms: List[str] = field(default_factory=list)
    query_used: str = ""
    lineage: str = ""
    detail: str = ""

    @property
    def is_resolved(self) -> bool:
        return self.status in (RESOLVED, RESOLVED_VIA_SYNONYM)

    def query_term(self) -> str:
        """What to send to a source service.

        The numeric taxid where it is known, because it cannot be misspelled and
        cannot drift when a name is revised. Otherwise the accepted name — never the
        underscored slug, which is what broke the original run.
        """
        return self.taxid or self.accepted_name or self.submitted_name

    def as_row(self) -> Dict[str, str]:
        return {
            "submitted_name": self.submitted_name,
            "species_id": self.species_id,
            "resolution_status": self.status,
            "accepted_scientific_name": self.accepted_name,
            "taxid": self.taxid,
            "rank": self.rank,
            "common_name": self.common_name,
            "synonyms": "; ".join(self.synonyms),
            "query_used": self.query_used,
            "lineage": self.lineage,
            "detail": self.detail,
        }


def normalise_name(name: str) -> str:
    """A submitted name as a scientific name.

    The web form and the species lists use ``equus_quagga``; NCBI wants ``Equus
    quagga``. Underscores become spaces and the genus is capitalised. Anything that
    already looks like a scientific name is left alone, so a caller cannot damage a
    correctly spelled subspecies trinomial by passing it through here.
    """
    cleaned = re.sub(r"[_\s]+", " ", (name or "").strip())
    if not cleaned:
        return ""
    if cleaned == cleaned.lower():
        parts = cleaned.split(" ")
        parts[0] = parts[0].capitalize()
        cleaned = " ".join(parts)
    return cleaned


def species_id(name: str) -> str:
    """The pipeline's internal identifier: lowercase, underscored, filesystem-safe."""
    return re.sub(r"[^a-z0-9]+", "_", (name or "").strip().lower()).strip("_")


def _get(url: str, timeout: float, retries: int = 2, sleep_s: float = 0.34) -> bytes:
    last: Optional[Exception] = None
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as handle:
                return handle.read()
        except (urllib.error.URLError, OSError) as exc:
            last = exc
            if attempt < retries:
                time.sleep(sleep_s * (attempt + 1))
    raise RuntimeError(f"{type(last).__name__}: {last}")


def _esearch_taxonomy(term: str, timeout: float) -> List[str]:
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(
        {"db": "taxonomy", "term": term, "retmode": "json", "retmax": "20"})
    payload = json.loads(_get(url, timeout).decode("utf-8", "replace"))
    return list((payload.get("esearchresult") or {}).get("idlist") or [])


def _efetch_taxonomy(taxid: str, timeout: float) -> Optional[ET.Element]:
    url = EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(
        {"db": "taxonomy", "id": taxid, "retmode": "xml"})
    root = ET.fromstring(_get(url, timeout).decode("utf-8", "replace"))
    return root.find("Taxon")


def _names_from(node: ET.Element) -> Dict[str, Any]:
    other = node.find("OtherNames")
    synonyms: List[str] = []
    common = ""
    if other is not None:
        for tag in ("Synonym", "EquivalentName", "Includes", "Anamorph"):
            synonyms += [e.text.strip() for e in other.findall(tag)
                         if e.text and e.text.strip()]
        for tag in ("GenbankCommonName", "CommonName"):
            hit = other.find(tag)
            if hit is not None and hit.text and not common:
                common = hit.text.strip()
    lineage = node.findtext("Lineage") or ""
    return {
        "accepted_name": (node.findtext("ScientificName") or "").strip(),
        "taxid": (node.findtext("TaxId") or "").strip(),
        "rank": (node.findtext("Rank") or "").strip(),
        "common_name": common,
        "synonyms": sorted(dict.fromkeys(synonyms)),
        "lineage": lineage,
    }


def resolve(name: str, *, timeout: float = 20.0, offline: bool = False,
            known: Optional[Dict[str, Dict[str, str]]] = None) -> TaxonIdentity:
    """Resolve one submitted species name.

    ``known`` is consulted first purely as a cache: the validated panel's taxids are
    already recorded and re-querying them on every run would be thirty needless
    requests against a service that rate-limits. It is a shortcut to the same answer,
    never a substitute for resolution — a name the cache does not hold still gets
    resolved rather than being echoed back as though it were a scientific name.
    """
    submitted = (name or "").strip()
    identity = TaxonIdentity(submitted_name=submitted, species_id=species_id(submitted))
    if not submitted:
        identity.detail = "empty species name"
        return identity

    canonical = normalise_name(submitted)

    if known:
        hit = _from_cache(canonical, submitted, known)
        if hit is not None:
            hit.species_id = species_id(hit.accepted_name) or identity.species_id
            return hit

    if offline:
        identity.status = OFFLINE
        identity.detail = ("offline mode: no taxonomy lookup was attempted, so the "
                           "species identity is unverified")
        return identity

    try:
        ids = _esearch_taxonomy(f'"{canonical}"[Scientific Name]', timeout)
        query = f'"{canonical}"[Scientific Name]'
        if not ids:
            # A name may be a synonym, a historical name or a subspecies rendering.
            # An unfielded search finds those; the check below is what stops it from
            # wandering to a different animal.
            ids = _esearch_taxonomy(canonical, timeout)
            query = canonical
    except RuntimeError as exc:
        identity.status = SERVICE_UNAVAILABLE
        identity.detail = str(exc)
        return identity

    identity.query_used = query
    if not ids:
        identity.status = NOT_FOUND
        identity.detail = (f"NCBI Taxonomy returned no taxon for {canonical!r}; "
                           "check the spelling of the scientific name")
        return identity

    try:
        node = _efetch_taxonomy(ids[0], timeout)
    except RuntimeError as exc:
        identity.status = SERVICE_UNAVAILABLE
        identity.detail = str(exc)
        return identity
    if node is None:
        identity.status = SERVICE_UNAVAILABLE
        identity.detail = f"taxonomy record {ids[0]} could not be parsed"
        return identity

    record = _names_from(node)
    identity.accepted_name = record["accepted_name"]
    identity.taxid = record["taxid"]
    identity.rank = record["rank"]
    identity.common_name = record["common_name"]
    identity.synonyms = record["synonyms"]
    identity.lineage = record["lineage"]
    identity.species_id = species_id(identity.accepted_name) or identity.species_id

    if len(ids) > 1:
        # Several taxa match. That is reportable, not resolvable by guessing: picking
        # the first hit is how a request for one animal becomes an analysis of
        # another. The first hit is still described so a reader can see the choice.
        identity.status = AMBIGUOUS
        identity.detail = (f"{len(ids)} taxa match {canonical!r} "
                           f"(first: {identity.accepted_name}, taxid {identity.taxid}); "
                           "submit the exact scientific name")
        return identity

    if _same_taxon(canonical, identity):
        identity.status = RESOLVED
        return identity

    if _is_published_synonym(canonical, identity):
        identity.status = RESOLVED_VIA_SYNONYM
        identity.detail = (f"{canonical!r} is a synonym NCBI Taxonomy maps to "
                           f"{identity.accepted_name} (taxid {identity.taxid})")
        return identity

    # The service answered with a taxon whose name is neither the submitted one nor
    # one of its published synonyms. That is a near-miss, and accepting it is exactly
    # the silent species substitution this module exists to prevent.
    identity.status = NOT_FOUND
    identity.detail = (
        f"NCBI Taxonomy answered {canonical!r} with {identity.accepted_name!r} "
        f"(taxid {identity.taxid}), which is not the requested taxon or a published "
        "synonym of it; the requested species was not substituted")
    identity.accepted_name = ""
    identity.taxid = ""
    return identity


def _from_cache(canonical: str, submitted: str,
                known: Dict[str, Dict[str, str]]) -> Optional[TaxonIdentity]:
    for key in (canonical, submitted, species_id(submitted), species_id(canonical)):
        entry = known.get(key) or known.get(key.casefold())
        if not entry or not entry.get("taxid"):
            continue
        return TaxonIdentity(
            submitted_name=submitted,
            status=RESOLVED,
            accepted_name=entry.get("ncbi_species") or canonical,
            taxid=str(entry["taxid"]),
            rank="species",
            common_name=entry.get("common_name", ""),
            query_used="internal_verified_cache",
            detail="taxid from the project's verified species table",
        )
    return None


def _norm(value: str) -> str:
    return re.sub(r"\s+", " ", (value or "").strip()).casefold()


def _same_taxon(canonical: str, identity: TaxonIdentity) -> bool:
    return _norm(canonical) == _norm(identity.accepted_name)


def _is_published_synonym(canonical: str, identity: TaxonIdentity) -> bool:
    """Whether the source itself lists the submitted name for this taxon.

    The test is on the synonyms NCBI publishes, not on string similarity. *Equus
    quagga* and *Equus caballus* share a genus and eight characters; only one of them
    is a zebra.
    """
    target = _norm(canonical)
    return any(_norm(s) == target for s in identity.synonyms)


