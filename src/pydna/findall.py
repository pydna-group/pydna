#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2013-2026 Björn Johansson
# SPDX-FileCopyrightText: 2023-2026 The Project Contributors
# SPDX-License-Identifier: BSD-3-Clause
"""
Find approximate occurrences of a subsequence within a linear or circular sequence.

This module uses `edlib <https://github.com/Martinsos/edlib>`_ to identify every
interval whose Levenshtein distance from a query sequence does not exceed a
specified threshold. Matches may contain substitutions, insertions, and
deletions. Results include zero-based coordinates, an extended
`CIGAR string <https://en.wikipedia.org/wiki/Sequence_alignment#CIGAR_Format>`_,
and a human-readable alignment.
"""

from typing import TypedDict

import edlib


IUPAC_DNA: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "U": frozenset("T"),
    "R": frozenset("GA"),
    "Y": frozenset("TC"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "B": frozenset("GTC"),
    "D": frozenset("GAT"),
    "H": frozenset("ACT"),
    "V": frozenset("GCA"),
    "N": frozenset("AGCT"),
}


IUPAC_EQUALITIES: list[tuple[str, str]] = [
    (left, right)
    for left, left_bases in IUPAC_DNA.items()
    for right, right_bases in IUPAC_DNA.items()
    if left < right and left_bases & right_bases
]


class FindResult(TypedDict):
    distance: int
    start: int
    stop: int
    cigar: str
    alignment: str


def findall(
    needle: str,
    haystack: str,
    max_edits: int = 0,
    circular: bool = False,
) -> list[FindResult]:
    """Find all approximate occurrences of a sequence.

    Every distinct interval in ``haystack`` whose Levenshtein distance from
    ``needle`` is less than or equal to ``max_edits`` is returned. Edit
    operations may be substitutions, insertions, or deletions. Comparisons use
    the extended IUPAC DNA alphabet, so ambiguous symbols match when their
    possible concrete bases overlap.

    Parameters
    ----------
    needle : str
        Sequence to search for.
    haystack : str
        Sequence to search.
    max_edits : int, default=0
        Maximum permitted Levenshtein distance.
    circular : bool, default=False
        Whether to allow matches to cross the boundary between the end and
        beginning of ``haystack``.

    Returns
    -------
    list of FindResult
        Matching intervals sorted by start coordinate, edit distance, and stop
        coordinate. Each result contains:

        ``distance``
            Levenshtein distance between ``needle`` and the matched interval.
        ``start``
            Zero-based inclusive start coordinate in ``haystack``.
        ``stop``
            Zero-based inclusive stop coordinate in ``haystack``. For a
            circular match crossing the origin, ``stop`` is smaller than
            ``start``.
        ``cigar``
            Extended CIGAR string produced by Edlib, using ``=``, ``X``, ``I``,
            and ``D`` operations.
        ``alignment``
            Three-line alignment containing the aligned query, match
            indicators, and aligned target.

    Raises
    ------
    TypeError
        If ``needle`` or ``haystack`` is not a string, or if ``max_edits`` is
        not an integer.
    ValueError
        If ``needle`` is empty or ``max_edits`` is negative.

    Notes
    -----
    The search is case-sensitive; input sequences are compared exactly as
    provided. IUPAC ambiguity is only applied to uppercase IUPAC symbols.
    Circular matches are restricted to at most one traversal of ``haystack``.

    Examples
    --------
    >>> result, *_ = findall(
    ...     "ACGTACGT", "TTTACGTTACGTAAA", max_edits=1, circular=False
    ... )
    >>> result["distance"]
    1
    >>> result["start"], result["stop"]
    (3, 12)
    >>> result["cigar"]
    '4=1D4='

    A circular match may cross the origin:

    >>> result = findall("TAAC", "ACGTTA", circular=True)[0]
    >>> result["start"], result["stop"]
    (4, 2)
    """
    if not isinstance(needle, str) or not isinstance(haystack, str):
        raise TypeError("needle and haystack must be strings")

    if not needle:
        raise ValueError("needle must not be empty")

    if not haystack:
        raise ValueError("haystack must not be empty")

    if not isinstance(max_edits, int):
        raise TypeError("max_edits must be an integer")

    if max_edits < 0:
        raise ValueError("max_edits must be a positive integer")

    needle_length = len(needle)
    haystack_length = len(haystack)

    # An interval differing by at most k edits can differ in length
    # from the needle by at most k characters.
    minimum_length = max(1, needle_length - max_edits)
    maximum_length = needle_length + max_edits

    if circular:
        # Do not allow a match to traverse the circular haystack more than once.
        maximum_length = min(maximum_length, haystack_length)

        if minimum_length > haystack_length:
            return []

        repetitions = (maximum_length - 1 + haystack_length - 1) // haystack_length
        searchable = haystack + haystack * repetitions
        starts = range(haystack_length)
    else:
        maximum_length = min(maximum_length, haystack_length)

        if minimum_length > haystack_length:
            return []

        searchable = haystack
        starts = range(haystack_length)

    matches: list[FindResult] = []

    for start in starts:
        for length in range(minimum_length, maximum_length + 1):
            end = start + length

            if not circular and end > haystack_length:
                break

            target = searchable[start:end]

            result = edlib.align(
                needle,
                target,
                mode="NW",
                task="path",
                k=max_edits,
                additionalEqualities=IUPAC_EQUALITIES,
            )

            distance = result["editDistance"]

            # Edlib returns -1 when the distance is greater than k.
            if distance < 0:
                continue

            nice = edlib.getNiceAlignment(
                result,
                needle,
                target,
            )

            alignment = "\n".join(
                (
                    nice["query_aligned"],
                    nice["matched_aligned"],
                    nice["target_aligned"],
                )
            )

            stop = end

            if circular:
                stop %= haystack_length

            matches.append(
                {
                    "distance": distance,
                    "start": start,
                    "stop": stop,
                    "cigar": result["cigar"],
                    "alignment": alignment,
                }
            )

    return sorted(
        matches,
        key=lambda match: (
            match["start"],
            match["distance"],
            match["stop"],
        ),
    )
