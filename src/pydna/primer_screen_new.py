#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fast primer screening
---------------------

This module provides fast primer screening using the Aho–Corasick string-search
algorithm. It is useful for PCR diagnostic purposes when given a list of primers
and a single sequence or list of sequences to analyze.

The primer list can consist of :class:`pydna.primer.Primer` objects returned by
:func:`pydna.parsers.parse_primers` or any objects with a ``seq`` attribute, such
as :class:`pydna.seqrecord.SeqRecord` or :class:`Bio.SeqRecord.SeqRecord`.

The Aho–Corasick algorithm efficiently finds all occurrences of a set of sequences
within a larger text. If the same primer list is used repeatedly, creating an
:class:`ahocorasick.Automaton` greatly speeds up repeated searches. See
:func:`make_automaton` for information on creating, saving, and loading such
automata.

Functions
---------

- :func:`forward_primers`
- :func:`reverse_primers`
- :func:`primer_pairs`
- :func:`flanking_primer_pairs`
- :func:`diff_primer_pairs`
- :func:`diff_primer_triplets`

References
----------

Aho–Corasick algorithm:
    https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_algorithm

This module uses `pyahocorasick`:
    Documentation: https://pyahocorasick.readthedocs.io/en/latest
    GitHub: https://github.com/WojciechMula/pyahocorasick
    PyPI: https://pypi.python.org/pypi/pyahocorasick
"""

# TODO: circular templates

from itertools import product
from itertools import combinations
from itertools import pairwise
from collections import defaultdict
from collections import Counter
from collections import namedtuple
from collections.abc import Callable
from collections.abc import Sequence

from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer

import ahocorasick
import warnings

warnings.warn(
    "The primer_screen module is experimental "
    "and not yet extensively tested. "
    "api may change in future versions.",
    category=FutureWarning,
)

amplicon_tuple = namedtuple(
    typename="amplicon_tuple", field_names="fp, rp, fposition, rposition, size"
)
primer_tuple = namedtuple(typename="primer_tuple", field_names="seq, fp, rp, size")


def closest_diff(nums: list[int]) -> int:
    """
    Smallest difference between two consecutive integers in a sorted list.

    Given a list of integers ex. 1, 5, 7, 11, 19, return the smallest absolute
    difference, in this case 7-5 = 2.

    >>> closest_diff([1, 5, 7, 11, 19])
    2

    Parameters
    ----------
    nums : list[int]
        List of integers.

    Raises
    ------
    ValueError
        At least two numbers are required.

    Returns
    -------
    int
        Diff, always >= 0.
    """
    if len(nums) < 2:
        raise ValueError("Need at least two numbers")

    nums = sorted(nums)
    min_diff = float("inf")

    for a, b in zip(nums, nums[1:]):
        diff = abs(a - b)
        if diff < min_diff:
            min_diff = diff
            x, y = a, b

    return abs(x - y)


def expand_iupac_to_dna(seq: str) -> list[str]:
    """
    Expand an extended IUPAC DNA string to unambiguous IUPAC nucleotide alphabet.

    Expands a string containing extended IUPAC code (ACGTURYSWKMBDHVN) including
    U for uracil into all possible DNA strings using only AGCT. Uses
    :func:`itertools.product` to compute the Cartesian expansion.

    Returns a list of strings.

    Example:

    >>> expand_iupac_to_dna("ATNG")
    ['ATAG', 'ATCG', 'ATGG', 'ATTG']

    Parameters
    ----------
    seq : str
        String containing extended IUPAC DNA.

    Returns
    -------
    list[str]
        List of strings in unambiguous IUPAC nucleotide alphabet.
    """
    _IUPAC_TO_DNA = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "U": "T",
        "R": "AG",  # puRine
        "Y": "CT",  # pYrimidine
        "S": "GC",  # Strong (3 H-bonds)
        "W": "AT",  # Weak (2 H-bonds)
        "K": "GT",  # Keto
        "M": "AC",  # aMino
        "B": "CGT",  # not A
        "D": "AGT",  # not C
        "H": "ACT",  # not G
        "V": "ACG",  # not T
        "N": "ACGT",  # any
        "X": "",
    }

    choices_per_pos = [_IUPAC_TO_DNA[ch] for ch in seq.upper()]
    return ["".join(t) for t in product(*choices_per_pos)]


def make_automaton(
    primer_list: Sequence[Primer | None], limit: int = 16
) -> ahocorasick.Automaton:
    """
    Create an :class:`ahocorasick.Automaton` for fast primer screening.

    This automaton can be reused across calls to :func:`forward_primers`,
    :func:`reverse_primers`, :func:`primer_pairs`, :func:`flanking_primer_pairs`,
    :func:`diff_primer_pairs`, and :func:`diff_primer_triplets`.

    Parameters
    ----------
    primer_list : :class:`collections.abc.Sequence` of :class:`pydna.primer.Primer`
        A list or tuple of primer objects.
    limit : int, optional
        Length of the 3' suffix used for matching.

    Returns
    -------
    :class:`ahocorasick.Automaton`
        Aho–Corasick automaton.
    """
    automaton = ahocorasick.Automaton()
    suffix_dict = defaultdict(list)

    for i, s in enumerate(primer_list):
        if not s or len(s) < limit:
            continue
        for footprint in expand_iupac_to_dna(str(s.seq)[-limit:].upper()):
            suffix_dict[footprint].append(i)

    for footprint, indices in suffix_dict.items():
        automaton.add_word(footprint, tuple(indices))

    automaton.make_automaton()
    return automaton


def callback(a: int, b: int) -> bool:
    """
    PCR product size quality control.

    This function accepts two integers and returns whether the difference meets a
    20% threshold. Used by :func:`diff_primer_pairs` and
    :func:`diff_primer_triplets`.

    Parameters
    ----------
    a, b : int
        Product sizes.

    Returns
    -------
    bool
        True if distinguishable, False otherwise.
    """
    return abs(a - b) >= 0.2 * max((a, b))


# --- Continued: Updated functions with enhanced Sphinx references ---


def forward_primers(
    seq: Dseqrecord,
    primer_list: Sequence[Primer | None],
    limit: int = 16,
    automaton: ahocorasick.Automaton = None,
) -> dict[int, list[int]]:
    """
    Identify forward primers annealing to a template.

    Uses :func:`make_automaton` to find all primer suffix matches.

    Parameters
    ----------
    seq : :class:`pydna.dseqrecord.Dseqrecord`
        Target DNA sequence.
    primer_list : :class:`collections.abc.Sequence` of :class:`pydna.primer.Primer`
        Primers to search.
    limit : int, optional
        3' annealing suffix length.
    automaton : :class:`ahocorasick.Automaton`, optional
        Automaton built by :func:`make_automaton`.

    Returns
    -------
    dict[int, list[int]]
        Primer-index → list of annealing positions.
    """

    automaton = automaton or make_automaton(primer_list, limit=limit)
    limit = automaton.get_stats()["longest_word"]

    fps = defaultdict(list)

    for end_index, ids in automaton.iter(str(seq.seq).upper()):
        for i in ids:
            fps[i].append(end_index + 1)

    return dict(fps)


def reverse_primers(
    seq: Dseqrecord,
    primer_list: Sequence[Primer],
    limit: int = 16,
    automaton: ahocorasick.Automaton = None,
) -> dict[int, list[int]]:
    """
    Identify reverse primers annealing to a template.

    Parameters
    ----------
    seq : :class:`pydna.dseqrecord.Dseqrecord`
        Target sequence.
    primer_list : :class:`collections.abc.Sequence` of :class:`pydna.primer.Primer`
        Primers to search.
    limit : int, optional
        3' annealing suffix length.
    automaton : :class:`ahocorasick.Automaton`, optional
        Automaton built by :func:`make_automaton`.

    Returns
    -------
    dict[int, list[int]]
        Primer-index → list of annealing positions.
    """

    automaton = automaton or make_automaton(primer_list, limit=limit)
    limit = automaton.get_stats()["longest_word"]

    rps = defaultdict(list)
    ln = len(seq)

    for end_index, ids in automaton.iter(str(seq.seq.reverse_complement()).upper()):
        for i in ids:
            rps[i].append(ln - (end_index + 1))

    return dict(rps)


def primer_pairs(
    seq: Dseqrecord,
    primer_list: Sequence[Primer],
    short: int = 500,
    long: int = 2000,
    limit: int = 16,
    automaton: ahocorasick.Automaton = None,
) -> list[amplicon_tuple]:
    """
    Identify valid forward/reverse primer pairs.

    Notes
    -----
    Uses :func:`forward_primers` and :func:`reverse_primers` internally.

    Parameters
    ----------
    seq : :class:`pydna.dseqrecord.Dseqrecord`
    primer_list : :class:`collections.abc.Sequence` of :class:`pydna.primer.Primer`
    short, long : int
        Minimum and maximum product sizes.
    limit : int
        3' matching suffix size.
    automaton : :class:`ahocorasick.Automaton`, optional

    Returns
    -------
    list[:class:`amplicon_tuple`]
        A list of PCR product descriptions.
    """

    automaton = automaton or make_automaton(primer_list, limit=limit)
    limit = automaton.get_stats()["longest_word"]

    fps = {
        fp: pos[0]
        for fp, pos in forward_primers(
            seq, primer_list, limit=limit, automaton=automaton
        ).items()
        if len(pos) == 1
    }

    rps = {
        rp: pos[0]
        for rp, pos in reverse_primers(
            seq, primer_list, limit=limit, automaton=automaton
        ).items()
        if len(pos) == 1
    }

    products = []
    for fp, fposition in fps.items():
        for rp, rposition in rps.items():
            size = len(primer_list[fp]) + rposition - fposition + len(primer_list[rp])
            if short <= size <= long and fposition <= rposition:
                products.append(amplicon_tuple(fp, rp, fposition, rposition, size))

    return products


def flanking_primer_pairs(
    seq: Dseqrecord,
    primer_list: Sequence[Primer],
    target: tuple[int, int],
    limit: int = 16,
    automaton: ahocorasick.Automaton = None,
) -> list[amplicon_tuple]:
    """
    Identify primer pairs flanking a region.

    Parameters
    ----------
    seq : :class:`pydna.dseqrecord.Dseqrecord`
    primer_list : :class:`collections.abc.Sequence` of :class:`pydna.primer.Primer`
    target : tuple[int, int]
        Start and end of target region.
    limit : int
    automaton : :class:`ahocorasick.Automaton`, optional

    Returns
    -------
    list[:class:`amplicon_tuple`]
    """

    automaton = automaton or make_automaton(primer_list, limit=limit)
    limit = automaton.get_stats()["longest_word"]

    begin, end = target
    assert begin < end

    amplicons = primer_pairs(
        seq,
        primer_list,
        short=end - begin,
        long=len(seq),
        limit=limit,
        automaton=automaton,
    )

    products = [a for a in amplicons if a.fposition >= begin and a.rposition >= end]
    return products[::-1]


def diff_primer_pairs(
    sequences: Sequence[Dseqrecord],
    primer_list: Sequence[Primer],
    short: int = 500,
    long: int = 1500,
    limit: int = 16,
    automaton: ahocorasick.Automaton = None,
    callback: Callable[[list], bool] = callback,
) -> list[tuple[primer_tuple]]:
    """
    Diagnostic primer pairs.

    Uses :func:`primer_pairs` and :func:`callback` to identify
    sequence-specific PCR product sizes.

    Parameters
    ----------
    sequences : Sequence[:class:`pydna.dseqrecord.Dseqrecord`]
    primer_list : Sequence[:class:`pydna.primer.Primer`]
    callback : :class:`collections.abc.Callable`

    Returns
    -------
    list[tuple[:class:`primer_tuple`]]
    """

    automaton = automaton or make_automaton(primer_list, limit=limit)
    limit = automaton.get_stats()["longest_word"]

    primer_pair_dict = defaultdict(dict)
    number_of_sequences = len(sequences)

    for seq in sequences:
        for fp, rp, *_, size in primer_pairs(
            seq, primer_list, short=short, long=long, limit=limit, automaton=automaton
        ):
            primer_pair_dict[frozenset((fp, rp))][size] = fp, rp, seq

    primer_pair_dict = {
        k: v for k, v in primer_pair_dict.items() if len(v) == number_of_sequences
    }
    primer_pair_dict = {
        k: v
        for k, v in primer_pair_dict.items()
        if all(callback(a, b) for a, b in pairwise(v.keys()))
    }

    result = []
    for primer_pair, seqd in primer_pair_dict.items():
        result.append(
            (
                closest_diff(seqd.keys()),
                tuple(
                    primer_tuple(s, fp, rp, size) for size, (fp, rp, s) in seqd.items()
                ),
            )
        )

    result.sort(reverse=True)
    return [b for a, b in result]


def diff_primer_triplets(
    sequences: Sequence[Dseqrecord],
    primer_list: Sequence[Primer],
    limit: int = 16,
    short: int = 500,
    long: int = 1500,
    automaton: ahocorasick.Automaton = None,
    callback: Callable[[list], bool] = callback,
) -> list[tuple[primer_tuple]]:
    """
        Diagnostic primer triplets.

        Uses :func:`primer_pairs`, :func:`callback`, and
        :func:`itertools.combinations`.

        Returns
        -------
        list[tuple[:class:`primer

    # (Full updated content inserted below.)

    # --- Final continuation of updated module ---

        …
    """

    automaton = automaton or make_automaton(primer_list, limit=limit)
    limit = automaton.get_stats()["longest_word"]
    number_of_sequences = len(sequences)

    # Collect all primer pairs for each sequence
    pp = {}
    for seq in sequences:
        pp[seq] = primer_pairs(
            seq,
            primer_list,
            short=short,
            long=long,
            limit=limit,
            automaton=automaton,
        )

    # Count occurrences of primer pairs across sequences
    pair_counter = Counter()
    for seq, tups in pp.items():
        for t in tups:
            pair = frozenset(t[:2])  # unordered FP/RP
            pair_counter[pair] += 1

    # Remove non-unique pairs
    pairs_to_remove = {pair for pair, count in pair_counter.items() if count > 1}
    for seq in pp:
        pp[seq] = [t for t in pp[seq] if frozenset(t[:2]) not in pairs_to_remove]

    primertrios = defaultdict(dict)

    # Combine primer pairs into triplets
    for seq1, seq2 in combinations(sequences, 2):
        for fp1, rp1, *_, size1 in pp[seq1]:
            for fp2, rp2, *_, size2 in pp[seq2]:
                primertrio = frozenset((fp1, rp1, fp2, rp2))
                if len(primertrio) == 3 and callback(size1, size2):
                    # If the same trio already exists, remove it (must be unique)
                    if primertrios[primertrio]:
                        del primertrios[primertrio]
                    else:
                        primertrios[primertrio][size1] = (fp1, rp1, seq1)
                        primertrios[primertrio][size2] = (fp2, rp2, seq2)

    # Construct final result list
    result = []
    for primertrio, seqd in primertrios.items():
        if len(seqd) == number_of_sequences and set(sequences) == set(
            s for *_, s in seqd.values()
        ):
            result.append(
                (
                    closest_diff(seqd.keys()),
                    tuple(
                        primer_tuple(s, fp, rp, size)
                        for size, (fp, rp, s) in seqd.items()
                    ),
                )
            )

    # Sort highest separation score first
    result.sort(key=lambda item: item[0], reverse=True)

    # Return just the primer_tuple collections
    return [b for a, b in result]


# --- END OF UPDATED MODULE ---
