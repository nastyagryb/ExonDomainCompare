#!/usr/bin/env python3
"""
build_species_registry.py

Create a validated species_registry.tsv from a plain-text species list.

The registry is the controlled entry point for downstream FGFR2 collection,
transcript selection, and comparative exon-domain mapping. Its job is to turn
whatever a user typed into identifiers the source services actually accept.

It used to do that from ``KNOWN_SPECIES`` alone — thirty entries, the validated FGFR2
panel — and for anything else it copied the submitted string into every name field and
left the taxid empty. Since runs are created with Ensembl-style slugs, the field later
used as an NCBI taxonomy query term became ``equus_quagga``, which no taxonomy service
recognises. That is what killed the Equus quagga run, four stages upstream of the
empty transcript table reported as its cause, and it would have killed every run for a
species outside the table.

Unknown species are now resolved against NCBI Taxonomy, which owns the answer. The
built-in table is kept, but only as a cache of already-verified taxids so a
thirty-species run does not make thirty needless requests. A species that resolves is
recorded with its accepted name, numeric taxid and published synonyms. A species that
does not resolve is recorded as unresolved with a reason — never as a different animal
that happens to have a better genome.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

ALLOWED_PREFERRED_SOURCES = {"ensembl_first", "ncbi_first", "ensembl_only", "ncbi_only"}
ALLOWED_ASSEMBLY_PREFERENCES = {"RefSeq", "GenBank", "Ensembl", "best_available"}

KNOWN_SPECIES: Dict[str, Dict[str, str]] = {
    "Homo sapiens": {"ensembl_species": "homo_sapiens", "ncbi_species": "Homo sapiens", "taxid": "9606", "common_name": "human", "clade": "mammal"},
    "Pan troglodytes": {"ensembl_species": "pan_troglodytes", "ncbi_species": "Pan troglodytes", "taxid": "9598", "common_name": "chimpanzee", "clade": "mammal"},
    "Gorilla gorilla gorilla": {"ensembl_species": "gorilla_gorilla_gorilla", "ncbi_species": "Gorilla gorilla gorilla", "taxid": "9595", "common_name": "gorilla", "clade": "mammal"},
    "Pongo abelii": {"ensembl_species": "pongo_abelii", "ncbi_species": "Pongo abelii", "taxid": "9601", "common_name": "orangutan", "clade": "mammal"},
    "Callithrix jacchus": {"ensembl_species": "callithrix_jacchus", "ncbi_species": "Callithrix jacchus", "taxid": "9483", "common_name": "marmoset", "clade": "mammal"},
    "Macaca mulatta": {"ensembl_species": "macaca_mulatta", "ncbi_species": "Macaca mulatta", "taxid": "9544", "common_name": "rhesus macaque", "clade": "mammal"},
    "Mus musculus": {"ensembl_species": "mus_musculus", "ncbi_species": "Mus musculus", "taxid": "10090", "common_name": "mouse", "clade": "mammal"},
    "Rattus norvegicus": {"ensembl_species": "rattus_norvegicus", "ncbi_species": "Rattus norvegicus", "taxid": "10116", "common_name": "rat", "clade": "mammal"},
    "Oryctolagus cuniculus": {"ensembl_species": "oryctolagus_cuniculus", "ncbi_species": "Oryctolagus cuniculus", "taxid": "9986", "common_name": "rabbit", "clade": "mammal"},
    "Canis lupus familiaris": {"ensembl_species": "canis_lupus_familiaris", "ncbi_species": "Canis lupus familiaris", "taxid": "9615", "common_name": "dog", "clade": "mammal"},
    "Felis catus": {"ensembl_species": "felis_catus", "ncbi_species": "Felis catus", "taxid": "9685", "common_name": "cat", "clade": "mammal"},
    "Bos taurus": {"ensembl_species": "bos_taurus", "ncbi_species": "Bos taurus", "taxid": "9913", "common_name": "cow", "clade": "mammal"},
    "Sus scrofa": {"ensembl_species": "sus_scrofa", "ncbi_species": "Sus scrofa", "taxid": "9823", "common_name": "pig", "clade": "mammal"},
    "Equus caballus": {"ensembl_species": "equus_caballus", "ncbi_species": "Equus caballus", "taxid": "9796", "common_name": "horse", "clade": "mammal"},
    "Ovis aries": {"ensembl_species": "ovis_aries", "ncbi_species": "Ovis aries", "taxid": "9940", "common_name": "sheep", "clade": "mammal"},
    "Monodelphis domestica": {"ensembl_species": "monodelphis_domestica", "ncbi_species": "Monodelphis domestica", "taxid": "13616", "common_name": "opossum", "clade": "mammal"},
    "Ornithorhynchus anatinus": {"ensembl_species": "ornithorhynchus_anatinus", "ncbi_species": "Ornithorhynchus anatinus", "taxid": "9258", "common_name": "platypus", "clade": "mammal"},
    "Gallus gallus": {"ensembl_species": "gallus_gallus", "ncbi_species": "Gallus gallus", "taxid": "9031", "common_name": "chicken", "clade": "bird"},
    "Meleagris gallopavo": {"ensembl_species": "meleagris_gallopavo", "ncbi_species": "Meleagris gallopavo", "taxid": "9103", "common_name": "turkey", "clade": "bird"},
    "Taeniopygia guttata": {"ensembl_species": "taeniopygia_guttata", "ncbi_species": "Taeniopygia guttata", "taxid": "59729", "common_name": "zebra finch", "clade": "bird"},
    "Anolis carolinensis": {"ensembl_species": "anolis_carolinensis", "ncbi_species": "Anolis carolinensis", "taxid": "28377", "common_name": "green anole", "clade": "reptile"},
    "Alligator mississippiensis": {"ensembl_species": "alligator_mississippiensis", "ncbi_species": "Alligator mississippiensis", "taxid": "8496", "common_name": "alligator", "clade": "reptile"},
    "Chrysemys picta bellii": {"ensembl_species": "chrysemys_picta_bellii", "ncbi_species": "Chrysemys picta bellii", "taxid": "8479", "common_name": "painted turtle", "clade": "reptile"},
    "Xenopus tropicalis": {"ensembl_species": "xenopus_tropicalis", "ncbi_species": "Xenopus tropicalis", "taxid": "8364", "common_name": "frog", "clade": "amphibian"},
    "Ambystoma mexicanum": {"ensembl_species": "ambystoma_mexicanum", "ncbi_species": "Ambystoma mexicanum", "taxid": "8296", "common_name": "axolotl", "clade": "amphibian"},
    "Danio rerio": {"ensembl_species": "danio_rerio", "ncbi_species": "Danio rerio", "taxid": "7955", "common_name": "zebrafish", "clade": "fish"},
    "Oryzias latipes": {"ensembl_species": "oryzias_latipes", "ncbi_species": "Oryzias latipes", "taxid": "8090", "common_name": "medaka", "clade": "fish"},
    "Gasterosteus aculeatus": {"ensembl_species": "gasterosteus_aculeatus", "ncbi_species": "Gasterosteus aculeatus", "taxid": "69293", "common_name": "stickleback", "clade": "fish"},
    "Takifugu rubripes": {"ensembl_species": "takifugu_rubripes", "ncbi_species": "Takifugu rubripes", "taxid": "31033", "common_name": "fugu", "clade": "fish"},
    "Oreochromis niloticus": {"ensembl_species": "oreochromis_niloticus", "ncbi_species": "Oreochromis niloticus", "taxid": "8128", "common_name": "tilapia", "clade": "fish"},
}

REGISTRY_FIELDS = [
    "species_id",
    "input_name",
    "scientific_name",
    "ensembl_species",
    "ncbi_species",
    "taxid",
    "common_name",
    "clade",
    "preferred_source",
    "assembly_preference",
    "status",
    "notes",
    # How the identity was established, so a downstream reader (or a reviewer of a
    # finished run) can see whether a name was verified against the source or merely
    # echoed back. The original run had no way to express the difference.
    "taxon_resolution_status",
    "accepted_scientific_name",
    "taxon_synonyms",
    "taxon_rank",
    "taxon_query_term",
]
WARNING_FIELDS = ["input_name", "warning_code", "warning"]

#: Registry statuses. ``unresolved_taxon`` is the honest outcome for a name the
#: taxonomy service does not know: the run should stop and say so rather than proceed
#: with a query term that cannot work.
STATUS_VERIFIED = "taxon_verified"
STATUS_VERIFIED_SYNONYM = "taxon_verified_via_synonym"
STATUS_UNRESOLVED = "unresolved_taxon"
STATUS_UNVERIFIED_OFFLINE = "taxon_unverified_offline"

@dataclass(frozen=True)
class RegistryBuildResult:
    rows: List[Dict[str, str]]
    warnings: List[Dict[str, str]]


def normalize_species_name(name: str) -> str:
    return re.sub(r"\s+", " ", name.strip())


def slug_species(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", normalize_species_name(name).lower()).strip("_")


def _build_known_index() -> Dict[str, Tuple[str, Dict[str, str]]]:
    """Map several normalized key forms -> (scientific_name, entry).

    Users / the run workflow pass Ensembl-style slugs (e.g. ``gallus_gallus``), while the
    built-in table is keyed by scientific name (``Gallus gallus``). Without this index a
    cached species misses its cache entry and has to be resolved against NCBI Taxonomy
    again on every run. Indexing by the scientific name, the Ensembl slug and a slugified
    form makes the lookup format-robust.
    """
    idx: Dict[str, Tuple[str, Dict[str, str]]] = {}
    for sci, entry in KNOWN_SPECIES.items():
        keys = {
            sci.casefold(),
            slug_species(sci),
            entry["ensembl_species"].casefold(),
            slug_species(entry["ensembl_species"]),
            entry["ncbi_species"].casefold(),
            slug_species(entry["ncbi_species"]),
        }
        for key in keys:
            idx.setdefault(key, (sci, entry))
    return idx


_KNOWN_INDEX = _build_known_index()


def lookup_known_species(name: str) -> Tuple[str, Dict[str, str]] | Tuple[None, None]:
    """Resolve an input species name (scientific or Ensembl slug) to a known entry."""
    for key in (name.casefold(), slug_species(name)):
        hit = _KNOWN_INDEX.get(key)
        if hit:
            return hit
    return None, None


def read_species_list(path: Path) -> List[str]:
    species: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            value = normalize_species_name(line.split("#", 1)[0])
            if value:
                species.append(value)
    return species


def validate_unique_species(species_names: Sequence[str]) -> List[Dict[str, str]]:
    seen = set()
    warnings: List[Dict[str, str]] = []
    for name in species_names:
        key = name.casefold()
        if key in seen:
            warnings.append({
                "input_name": name,
                "warning_code": "DUPLICATE_SPECIES",
                "warning": "Duplicate species entry detected; only the first occurrence was kept.",
            })
        seen.add(key)
    return warnings


def deduplicate_preserve_order(species_names: Sequence[str]) -> List[str]:
    seen = set()
    unique: List[str] = []
    for name in species_names:
        key = name.casefold()
        if key not in seen:
            unique.append(name)
            seen.add(key)
    return unique


def _resolve_taxon(name: str, offline: bool):
    """Resolve one name against NCBI Taxonomy, or report why it could not be."""
    from exondomaincompare.shared_gene_analysis import taxon_resolution as tr
    return tr.resolve(name, offline=offline, known=_taxid_cache())


def _taxid_cache() -> Dict[str, Dict[str, str]]:
    """The built-in table as a lookup cache, keyed the several ways callers spell names.

    These thirty taxids are already verified and belong to the validated panel, so
    re-querying them on every run would be thirty requests for answers the project
    already has. This is a shortcut to the same result, not a substitute for
    resolution: a name that is not here is resolved against the service.
    """
    cache: Dict[str, Dict[str, str]] = {}
    for sci, entry in KNOWN_SPECIES.items():
        for key in (sci, sci.casefold(), slug_species(sci),
                    entry["ensembl_species"], entry["ncbi_species"].casefold()):
            cache.setdefault(key, entry)
    return cache


def build_registry_rows(
    species_names: Sequence[str],
    default_preferred_source: str,
    default_assembly_preference: str,
    offline: bool = False,
) -> RegistryBuildResult:
    if default_preferred_source not in ALLOWED_PREFERRED_SOURCES:
        raise ValueError(f"Invalid preferred source: {default_preferred_source}")
    if default_assembly_preference not in ALLOWED_ASSEMBLY_PREFERENCES:
        raise ValueError(f"Invalid assembly preference: {default_assembly_preference}")

    from exondomaincompare.shared_gene_analysis import taxon_resolution as tr

    warnings = validate_unique_species(species_names)
    rows: List[Dict[str, str]] = []

    for name in deduplicate_preserve_order(species_names):
        identity = _resolve_taxon(name, offline)
        _sci, entry = lookup_known_species(name)

        # The Ensembl key stays the slug of the submitted name. Ensembl uses its own
        # keys and does not carry every NCBI taxon, so a species can legitimately
        # resolve at NCBI and be absent from Ensembl — as Equus quagga is.
        ensembl_key = entry["ensembl_species"] if entry else slug_species(name)
        # The NCBI-facing name is the accepted one, or empty when unresolved. Writing
        # the submitted slug here is precisely the bug that produced
        # "The taxonomy name 'equus_quagga' is not recognized"; an empty field makes a
        # downstream step stop with a clear reason instead of issuing a doomed query.
        ncbi_name = identity.accepted_name

        if identity.is_resolved:
            status = (STATUS_VERIFIED if identity.status == tr.RESOLVED
                      else STATUS_VERIFIED_SYNONYM)
            notes = (f"resolved against NCBI Taxonomy as {identity.accepted_name} "
                     f"(taxid {identity.taxid})")
        elif identity.status == tr.OFFLINE:
            status = STATUS_UNVERIFIED_OFFLINE
            # Offline, the built-in table is all there is. A known species keeps its
            # verified name; an unknown one gets no NCBI name, because guessing one is
            # what caused the original failure.
            ncbi_name = entry["ncbi_species"] if entry else ""
            notes = identity.detail
        else:
            status = STATUS_UNRESOLVED
            notes = identity.detail

        rows.append({
            "species_id": (entry["ensembl_species"] if entry
                           else (identity.species_id or slug_species(name))),
            "input_name": name,
            "scientific_name": (ncbi_name or tr.normalise_name(name)),
            "ensembl_species": ensembl_key,
            "ncbi_species": ncbi_name,
            "taxid": identity.taxid or (entry["taxid"] if entry else ""),
            "common_name": identity.common_name or (entry.get("common_name", "") if entry else ""),
            "clade": entry.get("clade", "") if entry else "",
            "preferred_source": default_preferred_source,
            "assembly_preference": default_assembly_preference,
            "status": status,
            "notes": notes,
            "taxon_resolution_status": identity.status,
            "accepted_scientific_name": identity.accepted_name,
            "taxon_synonyms": "; ".join(identity.synonyms),
            "taxon_rank": identity.rank,
            "taxon_query_term": identity.query_term() if identity.is_resolved else "",
        })

        if status == STATUS_UNRESOLVED:
            warnings.append({
                "input_name": name,
                "warning_code": "TAXON_NOT_RESOLVED",
                "warning": (f"{identity.detail} No source query was issued under an "
                            "unverified name."),
            })
        elif status == STATUS_UNVERIFIED_OFFLINE:
            warnings.append({
                "input_name": name,
                "warning_code": "TAXON_UNVERIFIED_OFFLINE",
                "warning": identity.detail,
            })
        elif status == STATUS_VERIFIED_SYNONYM:
            warnings.append({
                "input_name": name,
                "warning_code": "TAXON_RESOLVED_VIA_SYNONYM",
                "warning": identity.detail,
            })
    return RegistryBuildResult(rows=rows, warnings=warnings)


def validate_registry_schema(rows: Iterable[Dict[str, str]]) -> None:
    for i, row in enumerate(rows, start=1):
        missing = [field for field in REGISTRY_FIELDS if field not in row]
        if missing:
            raise ValueError(f"Row {i} is missing fields: {', '.join(missing)}")
        if not row["species_id"]:
            raise ValueError(f"Row {i} has empty species_id")
        if not re.fullmatch(r"[a-z0-9_]+", row["species_id"]):
            raise ValueError(f"Row {i} has invalid species_id: {row['species_id']}")
        if row["taxid"] and not row["taxid"].isdigit():
            raise ValueError(f"Row {i} has non-numeric taxid: {row['taxid']}")


def write_tsv(rows: Sequence[Dict[str, str]], path: Path, fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a validated species_registry.tsv from a species list.")
    parser.add_argument("--species_list", required=True, type=Path, help="Text file with one species name per line.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    parser.add_argument("--default_preferred_source", default="ensembl_first", choices=sorted(ALLOWED_PREFERRED_SOURCES))
    parser.add_argument("--default_assembly_preference", default="RefSeq", choices=sorted(ALLOWED_ASSEMBLY_PREFERENCES))
    parser.add_argument("--strict", action="store_true", help="Fail if unknown or duplicate species are detected.")
    parser.add_argument("--offline", action="store_true",
                        help="Skip the NCBI Taxonomy lookup. Species outside the "
                             "built-in table are then recorded as unverified rather "
                             "than resolved.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    species_names = read_species_list(args.species_list)
    if not species_names:
        raise ValueError("Species list is empty after removing blank lines and comments.")

    result = build_registry_rows(
        species_names,
        default_preferred_source=args.default_preferred_source,
        default_assembly_preference=args.default_assembly_preference,
        offline=args.offline,
    )
    validate_registry_schema(result.rows)

    if args.strict and result.warnings:
        for warning in result.warnings:
            print(f"[WARNING] {warning['input_name']}: {warning['warning']}", file=sys.stderr)
        raise ValueError("Strict mode failed because warnings were generated.")

    outdir = args.outdir.resolve()
    write_tsv(result.rows, outdir / "species_registry.tsv", REGISTRY_FIELDS)
    write_tsv(result.warnings, outdir / "species_registry_warnings.tsv", WARNING_FIELDS)

    print(f"[OK] species_registry.tsv rows: {len(result.rows)}")
    print(f"[OK] species_registry_warnings.tsv rows: {len(result.warnings)}")
    print(f"[OK] Wrote: {outdir / 'species_registry.tsv'}")
    print(f"[OK] Wrote: {outdir / 'species_registry_warnings.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
