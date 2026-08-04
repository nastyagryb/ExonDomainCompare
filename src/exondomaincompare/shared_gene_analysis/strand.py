#!/usr/bin/env python3
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
    return normalize_strand(value) == MINUS


def is_forward(value: Any) -> bool:
    return normalize_strand(value) == PLUS


def strand_symbol(value: Any) -> str:
    normalized = normalize_strand(value)
    return "" if normalized is None else ("+" if normalized == PLUS else "-")


def strand_sign(value: Any, default: int = PLUS) -> int:
    normalized = normalize_strand(value)
    return default if normalized is None else normalized


def same_strand(*values: Any) -> bool:
    normalized = [normalize_strand(v) for v in values]
    if not normalized or any(n is None for n in normalized):
        return False
    return len(set(normalized)) == 1


__all__ = ["PLUS", "MINUS", "normalize_strand", "is_reverse", "is_forward",
           "strand_symbol", "strand_sign", "same_strand"]
