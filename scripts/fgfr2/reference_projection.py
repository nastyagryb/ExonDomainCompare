from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
import human_reference_control as HRC

GAP_CHARS = {"-", ".", " ", "*"}

SOURCE_ANALYSED = "analysed_species"
SOURCE_CONTROL = "canonical_reference_control"
SOURCE_UNAVAILABLE = "unavailable"

METHOD_ALIGNED_ROW = "analysed_human_alignment_row"
METHOD_IDENTITY = "canonical_control_column_identity"
METHOD_GLOBAL = "canonical_control_global_alignment"


@dataclass
class ReferenceProjection:

    panel: str
    source: str
    method: str
    reference_sequence: str = ""
    protein_accession: str = ""
    transcript_accession: str = ""
    by_column: Dict[int, Tuple[Optional[int], str]] = field(default_factory=dict)
    control_agreement: str = ""
    note: str = ""

    @property
    def available(self) -> bool:
        return bool(self.by_column) and self.source != SOURCE_UNAVAILABLE

    @property
    def n_reference_positions(self) -> int:
        return sum(1 for idx, _ in self.by_column.values() if idx is not None)

    def at(self, column: int) -> Tuple[Optional[int], str]:
        return self.by_column.get(column, (None, ""))


def _is_gap(ch: str) -> bool:
    return ch in GAP_CHARS


def _map_from_row(sequence: str) -> Dict[int, Tuple[Optional[int], str]]:
    out: Dict[int, Tuple[Optional[int], str]] = {}
    index = 0
    for column, ch in enumerate(sequence):
        if _is_gap(ch):
            out[column] = (None, "")
        else:
            index += 1
            out[column] = (index, ch.upper())
    return out


def _consensus(rows: Sequence[str], width: int) -> str:
    letters: List[str] = []
    for column in range(width):
        counts: Dict[str, int] = {}
        for row in rows:
            if column < len(row) and not _is_gap(row[column]):
                ch = row[column].upper()
                counts[ch] = counts.get(ch, 0) + 1
        letters.append(max(counts, key=counts.get) if counts else "-")
    return "".join(letters)


def _global_align(reference: str, consensus: str) -> Optional[Tuple[str, str]]:
    try:
        from Bio import Align
        from Bio.Align import substitution_matrices
    except Exception:
        return None
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    try:
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    except Exception:
        aligner.match_score, aligner.mismatch_score = 1.0, -1.0
    aligner.open_gap_score, aligner.extend_gap_score = -11.0, -1.0
    try:
        alignment = next(iter(aligner.align(reference, consensus.replace("-", "X"))))
        return str(alignment[0]), str(alignment[1])
    except Exception:
        return None


def _project_control(reference: str, rows: Sequence[str],
                     width: int) -> Tuple[Dict[int, Tuple[Optional[int], str]], str]:
    if width == len(reference):
        return _map_from_row(reference), METHOD_IDENTITY

    aligned = _global_align(reference, _consensus(rows, width))
    if aligned is None:
        return {}, METHOD_GLOBAL
    ref_row, cons_row = aligned
    out: Dict[int, Tuple[Optional[int], str]] = {}
    ref_index, column = 0, 0
    for ref_ch, cons_ch in zip(ref_row, cons_row):
        if not _is_gap(ref_ch):
            ref_index += 1
        if _is_gap(cons_ch):
            continue
        if column < width:
            out[column] = ((ref_index, ref_ch.upper()) if not _is_gap(ref_ch) else (None, ""))
        column += 1
    for missing in range(column, width):
        out[missing] = (None, "")
    return out, METHOD_GLOBAL


def resolve(panel: str, items: Sequence[Tuple[str, str]], *,
            control: Optional[Dict[str, object]] = None,
            repo_root: Optional[Path] = None) -> ReferenceProjection:
    rows = [seq for _, seq in items]
    width = max((len(seq) for seq in rows), default=0)

    human_row = next((seq for header, seq in items
                      if header.lower().startswith(f"{HRC.REFERENCE_SPECIES}|{panel.lower()}")
                      or header.lower().startswith(f"{HRC.REFERENCE_SPECIES}|{panel}")), None)

    if control is None:
        try:
            control = HRC.load(repo_root)
        except HRC.ReferenceControlError:
            control = None
    control_ok = control is not None and HRC.is_available(control)
    entry = ((control or {}).get("reference") or {}).get(panel) or {} if control_ok else {}
    control_sequence = str(entry.get("cassette_sequence") or "")

    if human_row is not None:
        by_column = _map_from_row(human_row)
        analysed = "".join(ch for ch in human_row if not _is_gap(ch)).upper()
        if not control_ok:
            agreement = "control_unavailable"
        elif analysed == control_sequence:
            agreement = "analysed_human_matches_canonical_control"
        else:
            agreement = "analysed_human_differs_from_canonical_control"
        return ReferenceProjection(
            panel=panel, source=SOURCE_ANALYSED, method=METHOD_ALIGNED_ROW,
            reference_sequence=analysed, by_column=by_column,
            protein_accession=str(entry.get("protein_accession") or ""),
            transcript_accession=str(entry.get("transcript_accession") or ""),
            control_agreement=agreement,
            note=("Homo sapiens is an analysed species; its own model is the reference. "
                  "The canonical control is reported for comparison only."),
        )

    if not control_ok or not control_sequence:
        return ReferenceProjection(
            panel=panel, source=SOURCE_UNAVAILABLE, method="",
            note=("Homo sapiens is not analysed and the canonical human reference control "
                  "failed validation; human-referenced figures stay unavailable."),
        )

    by_column, method = _project_control(control_sequence, rows, width)
    if not by_column:
        return ReferenceProjection(
            panel=panel, source=SOURCE_UNAVAILABLE, method=method,
            reference_sequence=control_sequence,
            note=("The canonical human reference could not be projected onto this run's "
                  f"{panel} cassette columns ({width} columns vs {len(control_sequence)} "
                  "reference residues)."),
        )
    return ReferenceProjection(
        panel=panel, source=SOURCE_CONTROL, method=method,
        reference_sequence=control_sequence, by_column=by_column,
        protein_accession=str(entry.get("protein_accession") or ""),
        transcript_accession=str(entry.get("transcript_accession") or ""),
        control_agreement="external_control",
        note=("Homo sapiens is not part of this dataset; the validated human FGFR2 "
              f"{panel} cassette is shown as an external reference control."),
    )

