#!/usr/bin/env python3
"""One canonical reading of a genomic strand, for every source this project fetches.

The sources disagree on how to spell the same fact. Ensembl's REST API returns the strand
as the integer ``-1``; NCBI Datasets and RefSeq/Gnomon GFF3 write ``-``; some cached tables
carry the string ``"reverse"``. A production check written as ``strand == "-"`` is therefore
true for one source and false for another *describing the same gene on the same strand*.

That is not a cosmetic difference. When such a check decides whether a transcript's CDS
parts are read 5'→3', a minus-strand gene from the source with the unrecognised spelling is
assembled in reverse: its exons are projected onto the protein from the C-terminus back, so
its exon-to-protein map, its internal coding-exon boundaries and every cross-species
comparison built on them are wrong — while the gene from the other source is right. Two
orthologues then have no comparable boundary at all, which is how a two-species comparison
came to produce zero comparable groups.

So the spelling is normalised once, here, and nowhere else decides what a strand means.

    >>> normalize_strand("-1") == normalize_strand("-") == MINUS
    True
    >>> is_reverse(-1) and is_reverse("reverse")
    True
    >>> normalize_strand("") is None
    True
"""
from __future__ import annotations

from typing import Any, Optional

#: The canonical values. Integers, so arithmetic on an orientation stays possible.
PLUS = 1
MINUS = -1

#: Every spelling of the forward strand this repository has met, lowercased.
_PLUS_TOKENS = frozenset({"+", "+1", "1", "plus", "forward", "fwd", "f", "sense", "true"})

#: Every spelling of the reverse strand this repository has met, lowercased.
_MINUS_TOKENS = frozenset({"-", "-1", "minus", "reverse", "rev", "r", "antisense"})

#: Spellings that state the strand is not known. Distinguished from a missing value only in
#: that the source said so explicitly; both normalise to None.
_UNKNOWN_TOKENS = frozenset({"", ".", "?", "na", "n/a", "none", "null", "unknown", "0"})


def normalize_strand(value: Any) -> Optional[int]:
    """The strand as ``PLUS``, ``MINUS``, or ``None`` when the source does not say.

    Accepts the integers, floats and strings the adapters actually produce, in any case and
    with surrounding whitespace. ``None`` means unknown, and callers must decide what an
    unknown strand implies for them rather than defaulting it to forward here: silently
    treating "no information" as "forward" is how an unordered transcript passes for an
    ordered one.
    """
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        if value > 0:
            return PLUS
        if value < 0:
            return MINUS
        return None
    token = str(value).strip().lower()
    if token in _PLUS_TOKENS:
        return PLUS
    if token in _MINUS_TOKENS:
        return MINUS
    if token in _UNKNOWN_TOKENS:
        return None
    # A numeric string the token sets do not cover, e.g. "+2" from a malformed table.
    try:
        number = float(token)
    except ValueError:
        return None
    return PLUS if number > 0 else (MINUS if number < 0 else None)


def is_reverse(value: Any) -> bool:
    """Whether the feature is on the reverse strand. Unknown counts as not reverse.

    Use this only where "unknown" and "forward" genuinely lead to the same handling, such as
    drawing an arrow. Where the answer changes the biology — the order in which CDS parts are
    concatenated — call ``normalize_strand`` and handle ``None`` explicitly.
    """
    return normalize_strand(value) == MINUS


def is_forward(value: Any) -> bool:
    """Whether the feature is on the forward strand. Unknown counts as not forward."""
    return normalize_strand(value) == PLUS


def strand_symbol(value: Any) -> str:
    """``"+"``, ``"-"`` or ``""``, for display and for tables that store the symbol."""
    normalized = normalize_strand(value)
    return "" if normalized is None else ("+" if normalized == PLUS else "-")


def strand_sign(value: Any, default: int = PLUS) -> int:
    """``+1`` or ``-1`` for use as a multiplier, falling back to ``default`` when unknown."""
    normalized = normalize_strand(value)
    return default if normalized is None else normalized


def same_strand(*values: Any) -> bool:
    """Whether all given values describe one strand, across source spellings.

    ``"-"`` and ``-1`` are the same strand; an unknown value matches nothing, including
    another unknown, because two absences of information are not an agreement.
    """
    normalized = [normalize_strand(v) for v in values]
    if not normalized or any(n is None for n in normalized):
        return False
    return len(set(normalized)) == 1


__all__ = ["PLUS", "MINUS", "normalize_strand", "is_reverse", "is_forward",
           "strand_symbol", "strand_sign", "same_strand"]
