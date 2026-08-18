#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""Find fuzzy occurrences of a sequence using :mod:`regex`."""

from typing import TypedDict

import regex


class FindResult(TypedDict):
    """Description of one fuzzy match."""

    fuzzy_counts: tuple[int, int, int]
    start: int
    stop: int
    cigar: str
    alignment: str


def _alignment(
    needle: str,
    haystack: str,
    start: int,
    stop: int,
    changes: tuple[list[int], list[int], list[int]],
) -> tuple[str, str]:
    """Return an extended CIGAR string and a display alignment."""
    substitutions, insertions, deletions = map(set, changes)
    query: list[str] = []
    target: list[str] = []
    indicators: list[str] = []
    operations: list[str] = []
    query_position = 0
    target_position = start

    def append(query_base: str, indicator: str, target_base: str, op: str) -> None:
        query.append(query_base)
        indicators.append(indicator)
        target.append(target_base)
        operations.append(op)

    while query_position < len(needle) or target_position < stop:
        # regex reports a deletion at the target coordinate where the missing
        # query character would have occurred.
        if target_position in deletions and query_position < len(needle):
            append(needle[query_position], " ", "-", "I")
            query_position += 1
            deletions.remove(target_position)
            continue

        # An insertion is reported immediately after the extra target base.
        if target_position + 1 in insertions and target_position < stop:
            append("-", " ", haystack[target_position], "D")
            target_position += 1
            insertions.remove(target_position)
            continue

        if query_position >= len(needle) or target_position >= stop:
            # These branches handle edits at either end of a match.
            if query_position < len(needle):
                append(needle[query_position], " ", "-", "I")
                query_position += 1
            else:
                append("-", " ", haystack[target_position], "D")
                target_position += 1
            continue

        query_base = needle[query_position]
        target_base = haystack[target_position]
        if target_position in substitutions:
            append(query_base, ".", target_base, "X")
        else:
            append(query_base, "|", target_base, "=")
        query_position += 1
        target_position += 1

    cigar_parts: list[str] = []
    for operation in operations:
        if cigar_parts and cigar_parts[-1].endswith(operation):
            count = int(cigar_parts[-1][:-1]) + 1
            cigar_parts[-1] = f"{count}{operation}"
        else:
            cigar_parts.append(f"1{operation}")

    return "".join(cigar_parts), "\n".join(
        ("".join(query), "".join(indicators), "".join(target))
    )


def fuzzymatch(
    needle: str,
    haystack: str,
    limit: int = 15,
    substitutions: int = 0,
    insertions: int = 0,
    deletions: int = 0,
    flags: int = regex.ENHANCEMATCH | regex.IGNORECASE,
) -> list[FindResult]:
    """Find up to ``limit`` overlapping fuzzy matches of ``needle``.

    The three edit limits are independent. Coordinates are zero-based and
    ``stop`` is exclusive. In the CIGAR string, ``I`` consumes the query and
    ``D`` consumes the target, following Edlib's convention.

    Each :class:`regex.Match` exposes ``fuzzy_changes`` as three lists of
    absolute haystack coordinates: substitutions, insertions, and deletions.
    The coordinates are walked together with the characters of ``needle`` and
    the matched haystack interval. A substitution consumes one character from
    each sequence and produces ``X``. A regex deletion consumes only a query
    character, places a gap in the target, and produces ``I``. A regex
    insertion consumes only a target character, places a gap in the query,
    and produces ``D``; regex reports its coordinate immediately after that
    target character. Every other aligned pair produces ``=``. Consecutive
    identical operations are run-length encoded to form the extended CIGAR
    string, for example ``4=1X3=``.

    The human-readable alignment is constructed during the same walk. Its
    first row is the gapped query, its third row is the gapped target, and its
    middle row contains ``|`` for matches, ``.`` for substitutions, and a
    space for gaps. Thus the CIGAR and alignment describe exactly the same
    edit path.
    """
    if not isinstance(needle, str) or not isinstance(haystack, str):
        raise TypeError("needle and haystack must be strings")
    if not needle:
        raise ValueError("needle must not be empty")
    if not isinstance(limit, int):
        raise TypeError("limit must be an integer")
    if limit < 0:
        raise ValueError("limit must not be negative")

    edit_limits = (substitutions, insertions, deletions)
    if any(not isinstance(value, int) for value in edit_limits):
        raise TypeError("edit limits must be integers")
    if any(value < 0 for value in edit_limits):
        raise ValueError("edit limits must not be negative")
    if limit == 0:
        return []

    expression = (
        f"(?:{regex.escape(needle)})"
        f"{{s<={substitutions},i<={insertions},d<={deletions}}}"
    )
    results: list[FindResult] = []

    for match in regex.finditer(expression, haystack, flags=flags, overlapped=True):
        start, stop = match.span()
        cigar, alignment = _alignment(
            needle, haystack, start, stop, match.fuzzy_changes
        )
        results.append(
            {
                "fuzzy_counts": match.fuzzy_counts,
                "start": start,
                "stop": stop,
                "cigar": cigar,
                "alignment": alignment,
            }
        )
        if len(results) == limit:
            break

    return results
