"""human_reference_control.py — the canonical human FGFR2 IIIb/IIIc reference control.

Human-referenced FGFR2 figures compare an analysed species against the validated
human IIIb/IIIc cassette. Those figures must therefore work when *Homo sapiens* is
not part of the analysed dataset — a rat + rabbit run still has a human reference,
it simply has it as an external control rather than as an analysed row.

This module builds one immutable, versioned reference-control object from the
validated freeze and validates it before anyone may render against it.

Why it is built here and not read from the old cache
----------------------------------------------------
The previous control object was copied out of the example run's derived
``cassette_residue_index.json``. That index numbers IIIb and IIIc on a *single*
shared ``human_reference_residue_index`` axis, deduplicated by one shared set. IIIb
is 46 cassette positions and IIIc is 48, so the shared axis silently truncated IIIc
to 46 and dropped residues from IIIb — destroying the ``GVNTTDKEI`` IIIc marker and
turning IIIb positions 9–10 into gaps. A reference read through a lossy derived
index is not a reference, so this module reads the validated per-isoform source
table directly and keeps IIIb and IIIc on their own axes.

The freeze itself is never written; only this derived object outside it is.
"""
from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

SCHEMA_VERSION = 2
CONTROL_FILENAME = "human_reference_control.json"

# Source of truth: the validated 30-species freeze, read-only.
FREEZE_SOURCE_TABLE = (
    "results/final_30_until_interpro_prepare/13_final_pre_interpro_closure/"
    "tables/figure6B_species_resolved_IIIb_IIIc_cassette_residue_map.tsv"
)
FREEZE_TRUTH_TABLE = (
    "results/final_30_until_interpro_prepare/13_final_pre_interpro_closure/"
    "final_pre_interpro_truth_table.tsv"
)
BUNDLED_DATASET_ROOT = "datasets/fgfr2_30_species"

REFERENCE_SPECIES = "homo_sapiens"
REFERENCE_DISPLAY_NAME = "Homo sapiens"
REFERENCE_TAXON_ID = "9606"
REFERENCE_GENE = "FGFR2"

PANELS = ("IIIb", "IIIc")

# Regression markers. A reference that has lost these is not the validated human
# cassette and must never be rendered as one.
MARKERS = {"IIIb": "SGINSSN", "IIIc": "GVNTTDKEI"}

# The validated cassette lengths. IIIb and IIIc genuinely differ; a build that
# reports equal lengths has collapsed them onto one axis.
EXPECTED_LENGTHS = {"IIIb": 46, "IIIc": 48}

_AA_PROPERTY = {
    "A": "hydrophobic", "V": "hydrophobic", "L": "hydrophobic", "I": "hydrophobic",
    "M": "hydrophobic", "P": "hydrophobic",
    "F": "aromatic", "W": "aromatic", "Y": "aromatic",
    "S": "polar", "T": "polar", "N": "polar", "Q": "polar", "C": "polar",
    "D": "negative", "E": "negative",
    "K": "positive", "R": "positive", "H": "positive",
    "G": "special_case",
}


class ReferenceControlError(ValueError):
    """Raised when the human reference control fails its integrity contract."""


@dataclass
class PanelReference:
    """One validated human cassette panel (IIIb or IIIc) on its own axis."""

    panel: str
    reference_species: str
    taxon_id: str
    gene: str
    final_isoform_label: str
    transcript_accession: str
    protein_accession: str
    cassette_sequence: str
    cassette_length: int
    validated_marker: str
    source_identity: str
    source_table_checksum: str
    sequence_checksum: str
    schema_version: int
    residues: List[Dict[str, Any]]
    discriminating_positions: List[int]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _file_checksum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def aa_property(aa: Optional[str]) -> str:
    letter = str(aa or "").strip().upper()
    if not letter or letter in ("-", "."):
        return "gap"
    return _AA_PROPERTY.get(letter, "other")


# A combined-alignment column separates IIIb from IIIc when it is classified as one of
# these, or when a motif-map row flags it directly.
DISCRIMINATING_POSITION_CLASSES = (
    "isoform_discriminating_conserved",
    "IIIb_specific_conserved",
    "IIIc_specific_conserved",
)


def _int_or_none(value: Any) -> Optional[int]:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def is_discriminating_column(row: Dict[str, Any]) -> bool:
    """Whether one combined-alignment column distinguishes IIIb from IIIc."""
    flag = str(row.get("is_isoform_discriminating", "")).strip().lower()
    if flag in ("true", "1", "yes"):
        return True
    return str(row.get("position_class", "")).strip() in DISCRIMINATING_POSITION_CLASSES


def discriminating_positions_by_panel(rows: List[Dict[str, Any]]) -> Dict[str, set]:
    """Discriminating cassette positions **per panel**, never one shared set.

    The IIIb/IIIc comparison lives on the combined cassette alignment, the one axis
    both panels share. A panel's own residue numbering does not: IIIc carries a
    two-residue insertion, so from that point on the same alignment column is a
    different residue number in IIIb than in IIIc, and two columns exist only in IIIc.

    Collapsing both onto one numeric set — as the Figure 6B overlay and the cassette
    residue index both did — therefore painted the gold discriminating marker on the
    wrong residues: IIIb was marked at 16 and 17, which are IIIc-specific columns it
    does not have at all, and every IIIc annotation behind the insertion was shifted.
    Each combined column is mapped through *its own* panel index instead, which is why
    this returns one set per panel and callers must ask for the panel they are drawing.

    ``rows`` are the discriminating/motif rows carrying ``human_IIIb_reference_index``
    and ``human_IIIc_reference_index``; a row with no index for a panel (an
    isoform-specific column) contributes nothing to that panel.
    """
    out: Dict[str, set] = {panel: set() for panel in PANELS}
    fields = {"IIIb": "human_IIIb_reference_index", "IIIc": "human_IIIc_reference_index"}
    for row in rows:
        if not is_discriminating_column(row):
            continue
        for panel in PANELS:
            index = _int_or_none(row.get(fields[panel]))
            if index:
                out[panel].add(index)
    return out


def has_panel_indices(rows: List[Dict[str, Any]]) -> bool:
    """Whether the rows already carry both per-panel residue-index columns."""
    return any(_int_or_none(r.get("human_IIIb_reference_index")) is not None
               or _int_or_none(r.get("human_IIIc_reference_index")) is not None
               for r in rows)


def panel_indices_from_combined_alignment(
        rows: List[Dict[str, Any]],
        *,
        column_fields: Sequence[str] = ("combined_alignment_col", "MSA_column", "alignment_col"),
        aa_fields: Dict[str, Sequence[str]] = None) -> List[Dict[str, Any]]:
    """Add the per-panel residue indices to rows that only carry combined columns.

    The validated freeze predates the two ``human_III?_reference_index`` columns, and
    the freeze is read-only. The indices are recoverable without it: walking the
    combined alignment columns in order and counting the non-gap residues of each panel
    reproduces exactly the numbering the newer analysis writes out (verified against a
    run that carries both). Returns new dicts; the input rows are not modified.
    """
    fields = aa_fields or {
        "IIIb": ("human_IIIb_aa", "human_IIIb_aa_one_letter"),
        "IIIc": ("human_IIIc_aa", "human_IIIc_aa_one_letter"),
    }

    def column_of(row: Dict[str, Any]) -> int:
        for name in column_fields:
            value = _int_or_none(row.get(name))
            if value is not None:
                return value
        return 0

    def aa_of(row: Dict[str, Any], panel: str) -> str:
        for name in fields[panel]:
            if name in row:
                return str(row.get(name) or "").strip()
        return ""

    ordered = sorted(rows, key=column_of)
    seen = {panel: 0 for panel in PANELS}
    out: List[Dict[str, Any]] = []
    for row in ordered:
        enriched = dict(row)
        for panel in PANELS:
            aa = aa_of(row, panel)
            key = f"human_{panel}_reference_index"
            if aa and aa not in ("-", "."):
                seen[panel] += 1
                enriched[key] = seen[panel]
            else:
                enriched.setdefault(key, "")
        out.append(enriched)
    return out


def build(repo_root: Optional[Path] = None) -> Dict[str, Any]:
    """Build the reference-control object from the validated freeze.

    Reads the freeze read-only and returns the derived object; the caller decides
    where to persist it. Raises :class:`ReferenceControlError` when the freeze does
    not yield a usable reference rather than returning a degraded one.
    """
    root = Path(repo_root or _repo_root())
    def validated_path(relative: str) -> Path:
        legacy = root / relative
        if legacy.is_file():
            return legacy
        closure_relative = relative.removeprefix(
            "results/final_30_until_interpro_prepare/")
        return root / BUNDLED_DATASET_ROOT / closure_relative

    source = validated_path(FREEZE_SOURCE_TABLE)
    truth = validated_path(FREEZE_TRUTH_TABLE)
    if not source.is_file():
        raise ReferenceControlError(f"validated source table missing: {FREEZE_SOURCE_TABLE}")
    if not truth.is_file():
        raise ReferenceControlError(f"validated truth table missing: {FREEZE_TRUTH_TABLE}")

    source_checksum = _file_checksum(source)
    rows = _read_tsv(source)
    truth_rows = {
        r.get("isoform"): r for r in _read_tsv(truth)
        if r.get("species") == REFERENCE_SPECIES
    }
    # The freeze's Figure 6B table was written with the collapsed discriminating rule
    # and the freeze is read-only, so its per-row flag would give both panels the same
    # positions. The freeze's own discriminating analysis holds the truth; its
    # per-panel indices are recovered from the combined alignment it also stores.
    freeze_disc = validated_path(FREEZE_SOURCE_TABLE.replace(
        "tables/figure6B_species_resolved_IIIb_IIIc_cassette_residue_map.tsv",
        "MSA/final_isoform_discriminating_residues.tsv"))
    discriminating_by_panel: Dict[str, set] = {panel: set() for panel in PANELS}
    if freeze_disc.is_file():
        disc_rows = _read_tsv(freeze_disc)
        if not has_panel_indices(disc_rows):
            disc_rows = panel_indices_from_combined_alignment(disc_rows)
        discriminating_by_panel = discriminating_positions_by_panel(disc_rows)

    panels: Dict[str, Dict[str, Any]] = {}
    for panel in PANELS:
        human = [r for r in rows
                 if r.get("species") == REFERENCE_SPECIES and r.get("isoform") == panel]
        if not human:
            raise ReferenceControlError(f"freeze has no human {panel} cassette rows")
        # Each panel keeps its own residue axis; never merge the two indices.
        human.sort(key=lambda r: int(r.get("human_reference_residue_index") or 0))
        residues: List[Dict[str, Any]] = []
        for r in human:
            i = int(r.get("human_reference_residue_index") or 0)
            aa = (r.get("human_reference_aa") or "").strip().upper()
            residues.append({"i": i, "aa": aa, "property": aa_property(aa)})
        discriminating = sorted(discriminating_by_panel[panel]) or [
            int(r.get("human_reference_residue_index") or 0) for r in human
            if str(r.get("is_discriminating_position", "")).lower() == "true"]
        sequence = "".join(res["aa"] or "-" for res in residues)
        meta = truth_rows.get(panel, {})
        panels[panel] = PanelReference(
            panel=panel,
            reference_species=REFERENCE_SPECIES,
            taxon_id=REFERENCE_TAXON_ID,
            gene=REFERENCE_GENE,
            final_isoform_label=(meta.get("final_isoform_label") or panel),
            transcript_accession=(meta.get("transcript_id") or ""),
            protein_accession=(meta.get("protein_id") or ""),
            cassette_sequence=sequence,
            cassette_length=len(residues),
            validated_marker=MARKERS[panel],
            source_identity="FGFR2 IIIb/IIIc — 30 vertebrates (validated freeze)",
            source_table_checksum=source_checksum,
            sequence_checksum=hashlib.sha256(sequence.encode("utf-8")).hexdigest(),
            schema_version=SCHEMA_VERSION,
            residues=residues,
            discriminating_positions=discriminating,
        ).to_dict()

    data = {
        "schema_version": SCHEMA_VERSION,
        "contract": "fgfr2_human_reference_control_v2",
        "role": "human_reference_control",
        "validation_status": "validated",
        "reference_species": REFERENCE_SPECIES,
        "display_species_name": REFERENCE_DISPLAY_NAME,
        "taxon_id": REFERENCE_TAXON_ID,
        "gene": REFERENCE_GENE,
        "source_identity": "FGFR2 IIIb/IIIc — 30 vertebrates (validated freeze)",
        "source_table": FREEZE_SOURCE_TABLE,
        "source_table_checksum": source_checksum,
        "panels": list(PANELS),
        "reference": panels,
        "note": (
            "Validated human FGFR2 IIIb/IIIc cassette, reused as an external "
            "reference/control. Never counted as an analysed species and never a "
            "substitute for an analysed human model."
        ),
    }
    validate(data)
    return data


def validate(data: Optional[Dict[str, Any]]) -> List[str]:
    """Raise unless ``data`` is a complete, non-degenerate reference control.

    Returns the (empty) problem list on success so callers can also use it as a
    predicate via :func:`problems`.
    """
    issues = problems(data)
    if issues:
        raise ReferenceControlError("; ".join(issues))
    return issues


def problems(data: Optional[Dict[str, Any]]) -> List[str]:
    """Every reason ``data`` may not be rendered as the human reference."""
    issues: List[str] = []
    if not isinstance(data, dict) or not data:
        return ["reference control is empty"]
    if data.get("schema_version") != SCHEMA_VERSION:
        issues.append(f"schema_version {data.get('schema_version')} != {SCHEMA_VERSION}")
    reference = data.get("reference") or {}
    for panel in PANELS:
        entry = reference.get(panel) or {}
        tag = f"{panel}"
        if not entry:
            issues.append(f"{tag}: missing panel")
            continue
        sequence = str(entry.get("cassette_sequence") or "")
        if not sequence:
            issues.append(f"{tag}: empty reference sequence")
            continue
        if set(sequence) <= {"-", ".", " "}:
            issues.append(f"{tag}: reference is all gaps")
        if MARKERS[panel] not in sequence:
            issues.append(f"{tag}: validated marker {MARKERS[panel]} absent from reference")
        if len(sequence) != EXPECTED_LENGTHS[panel]:
            issues.append(f"{tag}: cassette length {len(sequence)} != "
                          f"expected {EXPECTED_LENGTHS[panel]}")
        if entry.get("cassette_length") != len(sequence):
            issues.append(f"{tag}: recorded cassette_length {entry.get('cassette_length')} "
                          f"disagrees with sequence ({len(sequence)})")
        digest = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
        if entry.get("sequence_checksum") != digest:
            issues.append(f"{tag}: sequence checksum mismatch")
        if not entry.get("protein_accession") or not entry.get("transcript_accession"):
            issues.append(f"{tag}: missing protein/transcript accession")

    # A swap is the one corruption both panels can pass individually: each panel
    # must carry its own marker and not the other panel's.
    for panel in PANELS:
        other = "IIIc" if panel == "IIIb" else "IIIb"
        sequence = str((reference.get(panel) or {}).get("cassette_sequence") or "")
        if sequence and MARKERS[other] in sequence and MARKERS[panel] not in sequence:
            issues.append(f"{panel}: carries the {other} marker (IIIb/IIIc reference swapped)")
    return issues


def is_available(data: Optional[Dict[str, Any]]) -> bool:
    """Whether a human-referenced figure may be registered as available."""
    return not problems(data)


def control_path(repo_root: Optional[Path] = None) -> Path:
    root = Path(repo_root or _repo_root())
    return root / "results" / "web_state" / CONTROL_FILENAME


def load(repo_root: Optional[Path] = None, *, rebuild: bool = False) -> Dict[str, Any]:
    """Return the validated reference control, rebuilding it when unusable.

    A cached object that fails validation is never returned: the whole point of
    the control is that a figure rendered against it is trustworthy. Normal reads
    build in memory so a release checkout remains unchanged.
    """
    path = control_path(repo_root)
    if not rebuild and path.is_file():
        try:
            cached = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            cached = None
        if cached and is_available(cached):
            return cached
    data = build(repo_root)
    if rebuild:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return data


def panel_sequence(data: Dict[str, Any], panel: str) -> str:
    return str(((data.get("reference") or {}).get(panel) or {}).get("cassette_sequence") or "")


def panel_residues(data: Dict[str, Any], panel: str) -> List[Dict[str, Any]]:
    return list(((data.get("reference") or {}).get(panel) or {}).get("residues") or [])


def _main(argv: Optional[List[str]] = None) -> int:  # pragma: no cover - CLI helper
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rebuild", action="store_true", help="rebuild from the freeze")
    ap.add_argument("--check", action="store_true", help="validate only, write nothing")
    args = ap.parse_args(argv)

    data = build() if args.check else load(rebuild=args.rebuild)
    for panel in PANELS:
        entry = data["reference"][panel]
        print(f"{panel}: len={entry['cassette_length']} marker={entry['validated_marker']} "
              f"protein={entry['protein_accession']}")
        print(f"     {entry['cassette_sequence']}")
    print(f"OK — reference control validated (schema v{data['schema_version']}).")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
