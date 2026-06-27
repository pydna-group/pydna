#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2013-2026 Björn Johansson
# SPDX-FileCopyrightText: 2023-2026 The Project Contributors
# SPDX-License-Identifier: BSD-3-Clause

"""
Fast primer screening
---------------------

This module provides fast primer screening using the Aho-Corasick string-search
algorithm. It is useful for PCR diagnostic purposes when given a list of primers
and a single sequence or list of sequences to analyze.

The primer list can consist of `Primer` objects returned by :func:`pydna.parsers.parse_primers`
or any objects with a ``seq`` attribute, such as :class:`pydna.seqrecord.SeqRecord`
or :class:`Bio.SeqRecord.SeqRecord`.

The Aho-Corasick algorithm efficiently finds all occurrences of a set of sequences
within a larger text. If the same primer list is used repeatedly, creating an
automaton greatly speeds up repeated searches. See :func:`make_automaton` for
information on creating, saving, and loading such automata.

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

Aho-Corasick algorithm:
    https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_algorithm

This module uses `pyahocorasick`:
    Documentation: https://pyahocorasick.readthedocs.io/en/latest
    GitHub: https://github.com/WojciechMula/pyahocorasick
    PyPI: https://pypi.python.org/pypi/pyahocorasick
"""

from operator import attrgetter
from itertools import product
from itertools import combinations
from collections import defaultdict
from collections import namedtuple
from collections.abc import Callable
from collections.abc import Sequence

from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer
from pydna.utils import _directed_interval_overlap

import ahocorasick

import warnings

from Bio.Data.IUPACData import ambiguous_dna_values

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


def contained(a: int, b: int, x: int, y: int, L: int, circular=True) -> bool:
    """
    Test whether interval (a, b) is contained in interval (x, y)
    in a sequence of length L. The sequence may be linear or circular.

    Coordinates must satisfy 0 <= coordinate <= L.

    If circular=False, intervals are treated as ordinary linear intervals.

    If circular=True, intervals are treated as directed circular intervals:
    - (4, 8) means 4 -> 8
    - (9, 1) means 9 -> L -> 0 -> 1
    """

    assert L >= 1, "L must be a positive integer"

    def valid_coordinate(n: int) -> bool:
        return 0 <= n <= L

    if not all(valid_coordinate(n) for n in (a, b, x, y)):
        return False

    if not circular:
        assert a <= b, "In a linear interval, a <= b"
        assert x <= y, "In a linear interval, x <= y"
        return x <= a and b <= y

    # The directed circular interval (a, b) is contained in (x, y) when arc 1
    # is fully contained in arc 2 (see pydna.utils._directed_interval_overlap).
    return bool(_directed_interval_overlap(a, b, x, y, L).contained_1_in_2)


def closest_pair_and_diff(nums: list[int]) -> int:
    """
    Smallest difference between two consecutive integers in a sorted list.

    Given a list of integers eg. 1, 5, 7, 11, 19, return the consequtive pair
    with the smallest difference and the difference. If more than one pair
    gives the same difference, return the larger pair.

    >>> nums = [1, 5, 7, 11, 19, 21]
    >>> closest_pair_and_diff(nums)
    ((19, 21), 2)


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
    assert len(nums) > 1, "Need at least two numbers"
    nums.sort()
    diff_dict = {b - a: (a, b) for a, b in zip(nums, nums[1:])}
    diff = min(diff_dict.keys())
    return diff_dict[diff], diff


def expand_iupac_to_dna(seq: str) -> list[str]:
    """
    Expand an extended IUPAC DNA string to unambiguous IUPAC nucleotide alphabet.

    Expands a string containing extended IUPAC code (ACGTURYSWKMBDHVN) including
    U for uracil into all possible DNA strings using only AGCT.

    Returns a list of strings.

    Example:

    >>> expand_iupac_to_dna("ATNG")
    ['ATGG', 'ATAG', 'ATTG', 'ATCG']
    >>> x = expand_iupac_to_dna("ACGTURYSWKMBDHVN")
    >>> len(x)
    20736


    Parameters
    ----------
    seq : str
        String containing extended IUPAC DNA.

    Returns
    -------
    list[str]
        List of strings in unambiguous IUPAC nucleotide alphabet.

    """
    custom_dict = {**ambiguous_dna_values}
    # Include RNA
    custom_dict["U"] = "T"
    choices_per_pos = [custom_dict[ch] for ch in seq.upper()]
    # Cartesian product of all position choices
    return ["".join(tup) for tup in product(*choices_per_pos)]


def make_automaton(
    primer_list: Sequence[Primer | None], limit: str = 16
) -> ahocorasick.Automaton:
    """
    Aho-Corasick automaton for a list of primers.

    An automaton `here <https://github.com/WojciechMula/pyahocorasick>`__ can
    be made prior to primer screening for a list of Primer
    objects for faster primer search.


    This automaton can be reused as an optional argument across calls to :func:`forward_primers`,
    :func:`reverse_primers`, :func:`primer_pairs`, :func:`flanking_primer_pairs`,
    :func:`diff_primer_pairs`, and :func:`diff_primer_triplets`.

    The primer list can contain None, this can be used to remove primers
    from the primer_list for the automaton, while keeping the original index
    for each primer.

    The limit is the part of the primer used to find annealing positions.
    The automaton processes the uppercase 3' part of each primer up to `limit`.
    It has to be rebuilt if a different limit is needed.

    The primers can contain ambiguous bases from the extended IUPAC DNA alphabet.

    The automaton can be saved and loaded like this (from the pyahocorasick docs):

    ::

        import pickle
        from pydna import primer_screen

        # build automaton
        atm = make_automaton(pl, limit = 16)

        # save automaton
        atm.save("atm.automaton", pickle.dumps)

        # load automaton
        import ahocorasick
        atm = ahocorasick.load(path, pickle.loads)

        # use automaton
        fps = forward_primers(template, primer_list, automaton=atm)


    Parameters
    ----------
    primer_list : list[Primer] | tuple[Primer]
        This is a list of pydna.primer.Primer objects or
        any object with a seq property such as Bio.SeqRecord.SeqRecord.
    limit : str, optional
        This is the primer part in  the 3'-end that has to
        anneal. The default is 16.

    Returns
    -------
    ahocorasick.Automaton
        pyahocorasick automaton made for the list of Primer objects.

    """
    automaton = ahocorasick.Automaton()

    suffix_dict = defaultdict(list)

    for i, s in enumerate(primer_list):
        # filter for primers that evaluate to False such as None
        # or primers that are too short.
        if not s or (len(s) < limit):
            continue
        # Primers may share suffix, so primer indices pertaining to a
        # certain suffix are collected together.
        for footprint in expand_iupac_to_dna(str(s.seq)[-limit:].upper()):
            suffix_dict[footprint].append(i)

    for footprint, indices in suffix_dict.items():
        automaton.add_word(footprint, tuple(indices))

    automaton.make_automaton()

    return automaton


def callback(a: int, b: int) -> bool:
    """
    PCR product sizes quality control.

    This function accepts two integers representing PCR product sizes
    and returns True or False indicating the ease with which the size
    differences can be distinguished on a typical agarose gel.

    Parameters
    ----------
    a : int
        One size.
    b : int
        Another size.

    Returns
    -------
    bool
        True if successful, False otherwise.

    """
    # The length difference has to be 20%
    # of the size of the larger fragment
    return abs(a - b) >= 0.2 * max((a, b))


def forward_primers(
    seq: Dseqrecord,
    primer_list: Sequence[Primer | None],
    limit: int = 16,
    automaton: ahocorasick.Automaton = None,
) -> dict[int, list[int]]:
    """
    Forward primers from `primer_list` annealing to `seq` with at least `limit`
    base pairs.

    The optional automaton can speed up the primer search if the same primer
    list is often used, see :func:`make_automaton` for more information.

    The resulting dict has the form:

    ::

        { primer_A_index : [location1, location2, ...]
          primer_B_index : [location1, location2, ...] }

    Where a key such as primer_A_index (integer) is the index for a primer
    in `primer_list` and the value is a list of locations (integers) where
    the primer binds.

    The concept of location is the same as used in :mod:`pydna.primer`.
    The forward primer in the figure below anneals at position 14 on the
    template.

    ::

         5-gtcatgatctagtcgatgtta-3
          |||||||||||||||||||||

                 5'-tagtcg-3' = forward primer, location = 14
                    ||||||
          |||||||||||||||||||||
         3-cagtactagatcagctacaat-5
                         |
           012345678911111111112 position
                     01234567890



    Parameters
    ----------
    seq : Dseqrecord
        Target sequence to find primer annealing positions.
    primer_list : list[Primer] | tuple[Primer]
        This is a list of pydna.primer.Primer objects or any object
        with a seq property such as Bio.SeqRecord.SeqRecord.
    limit : str, optional
        This is the part at the 3'-end of each primer that has to
        anneal. The default is 16.
    automaton : ahocorasick.Automaton, optional
        Automaton made with the :func:`make_automaton`. The default is None.

    Returns
    -------
    dict[int, list[int]]
        Dict of lists where keys are primer indices in primer_list and
        values are lists with primer locations.

    """
    assert primer_list, "primer_list must not be empty."

    # if no automaton is given, we make one.
    automaton = automaton or make_automaton(primer_list, limit=limit)

    # The limit is taken from automaton stats.
    # If the automaton is given, the limit argument will be ignored.
    limit = automaton.get_stats()["longest_word"]

    # A defaultdict of lists is used to collect primer locations since
    # different primers can anneal in the same place.
    fps = defaultdict(list)
    # if seq is circular, we have to look across the origin of the sequence
    seq = seq[:] + seq[: limit - 1] if seq.circular else seq

    for end_index, ids in automaton.iter(str(seq.seq).upper()):
        for i in ids:
            fps[i].append(end_index + 1)

    return dict(fps)


def reverse_primers(
    seq: Dseqrecord,
    primer_list: Sequence[Primer | None],
    limit: int = 16,
    automaton: ahocorasick.Automaton = None,
) -> dict[int, list[int]]:
    """
    Primers from `primer_list` annealing in reverse to `seq` with at least
    `limit` base pairs.

    The optional automaton can speed up the primer search if the same primer
    list is often used, see :func:`make_automaton` for more information.

    The resulting dict has the form:

    ::

        { primer_A_index : [location1, location2, ...]
          primer_B_index : [location1, location2, ...] }

    Where a key such as primer_A_index (integer) is the index for a primer
    in `primer_list` and the value is a list of locations (integers) where
    the primer binds.

    The concept of location is the same as used in :mod:`pydna.primer`.
    The reverse primer below anneals at position 9.

    ::

        5-gtcatgatctagtcgatgtta-3
          |||||||||||||||||||||
                   ||||||
                 3-atcagc-5 = reverse primer, location = 9

          |||||||||||||||||||||
        3-cagtactagatcagctacaat-5
                   |
          012345678911111111112 position
                    01234567890


    Parameters
    ----------
    seq : Dseqrecord
        Target sequence to find primer annealing positions.
    primer_list : list[Primer] | tuple[Primer]
        This is a list of pydna.primer.Primer objects or any object
        with a seq property such as Bio.SeqRecord.SeqRecord.
    limit : str, optional
        This is the part in the 3'-end of each primer that has to
        anneal. The default is 16.
    automaton : ahocorasick.Automaton, optional
        Automaton made with the :func:`make_automaton`. The default is None.

    Returns
    -------
    dict[int, list[int]]
        Dict of lists where keys are primer indices in primer_list and
        values are lists with primer locations.

    """
    assert primer_list, "primer_list must not be empty."

    # if no automaton is given, we make one.
    automaton = automaton or make_automaton(primer_list, limit=limit)

    # The limit is taken from automaton stats.
    # If the automaton is given, the limit argument will be ignored.
    limit = automaton.get_stats()["longest_word"]

    # A defaultdict of lists is used to collect primer locations since
    # different primers can anneal in the same place.
    rps = defaultdict(list)
    # if seq is circular, we have to look across the origin of the sequence
    seq = seq[:] + seq[: limit - 1] if seq.circular else seq
    ln = len(seq)

    # We use the reverse complement of the sequence instead of taking the
    # reverse complement of each primer.
    for end_index, ids in automaton.iter(str(seq.seq.reverse_complement()).upper()):
        for i in ids:
            rps[i].append(ln - (end_index + 1))

    return dict(rps)


def primer_pairs(
    seq: Dseqrecord,
    primer_list: Sequence[Primer | None],
    short: int = 500,
    long: int = 2000,
    limit: int = 16,
    automaton: ahocorasick.Automaton = None,
) -> list[amplicon_tuple[int, int, int, int, int]]:
    """
    Primer pairs that form PCR products larger than `short` and smaller
    than `long`.

    The PCR product size includes the PCR primers. Only unique primer pairs
    are returned. This means that the forward and reverse primers can only
    bind in one position on the template each.

    If you suspect that primers bind on multiple locations, use the
    :func:`forward_primers` and :func:`reverse_primers` functions.

    The function returns a list of flat 5-namedtuples of integers and
    integers with this form:

    ::

        [
         ((index_fp1, index_rp1, position_fp1, position_rp1, size1),
         ((index_fp2, index_rp2, position_fp2, position_rp2, size2),
          ]


    The indices are the `primer_list` indices and positions are the positions of
    the primers as described in :func:`forward_primers` and :func:`reverse_primers`
    functions.
    The size includes the length of each primer, so it is the true total length
    of the PCR product.

    Parameters
    ----------
    seq : Dseqrecord
        Target sequence to find primer annealing positions.
    primer_list : list[Primer] | tuple[Primer]
        This is a list of pydna.primer.Primer objects or any object
        with a seq property such as Bio.SeqRecord.SeqRecord.
    limit : str, optional
        This is the part in the 3'-end of each primer that has to
        anneal. The default is 16.
    short : int, optional
        Lower limit for the size of the PCR products. The default is 500.
    long : int, optional
        Upper limit for the size of the PCR products. The default is 1500.
    automaton : ahocorasick.Automaton, optional
        Automaton made with the :func:`make_automaton`. The default is None.

    Returns
    -------
    list[tuple(int, int, int, int, int)]
        List of tuples (index_fp, position_fp, index_rp, position_rp, size)

    """
    assert primer_list, "primer_list must not be empty."

    # if no automaton is given, we make one.
    automaton = automaton or make_automaton(primer_list, limit=limit)

    # The limit is taken from automaton stats.
    # If the automaton is given, the limit argument will be ignored.
    limit = automaton.get_stats()["longest_word"]

    # Unique forward primers are collected
    fps = {
        fp: pos[0]
        for fp, pos in forward_primers(seq, primer_list, automaton=automaton).items()
        if len(pos) == 1
    }

    # Unique reverse primers are collected
    rps = {
        rp: pos[0]
        for rp, pos in reverse_primers(seq, primer_list, automaton=automaton).items()
        if len(pos) == 1
    }
    products = []

    for fp, fposition in fps.items():
        for rp, rposition in rps.items():
            # We calculate the size of a potential PCR product
            size = len(primer_list[fp]) + rposition - fposition + len(primer_list[rp])
            # If the size falls within long and short, the data is kept.
            if short <= size <= long and fposition <= rposition:
                products.append(amplicon_tuple(fp, rp, fposition, rposition, size))
    # if sequence is circular we also look at forward primers that sits after the
    # reverse primer and amplify across the origin
    if seq.circular:
        ln = len(seq)
        for fp, fposition in fps.items():
            for rp, rposition in rps.items():
                if fposition > rposition:
                    size = (
                        len(primer_list[fp])
                        + rposition
                        + (ln - fposition)
                        + len(primer_list[rp])
                    )
                    if short <= size <= long:
                        products.append(
                            amplicon_tuple(
                                fp, rp, fposition % len(seq), rposition % len(seq), size
                            )
                        )
    return products


def flanking_primer_pairs(
    seq: Dseqrecord,
    primer_list: Sequence[Primer | None],
    target: tuple[int, int] = (None, None),
    limit: int = 16,
    automaton: ahocorasick.Automaton = None,
) -> list[amplicon_tuple[int, int, int, int, int]]:
    """
    Primer pairs that flank a target position (begin..end). This means that
    forward primers have to bind before or at the begin position and reverse primers
    have to bind at or after the end position.

    The function returns a list of the same flat 5-namedtuples of integers returned
    from the :func:`primer_pairs` function.

    ::

        [
         (index_fp1, position_fp1, index_rp1, position_rp1, size1),
         (index_fp2, position_fp2, index_rp2, position_rp2, size2),
         ]


    Parameters
    ----------
    seq : Dseqrecord
        Target sequence to find primer annealing positions.
    primer_list : list[Primer] | tuple[Primer]
        This is a list of pydna.primer.Primer objects or any object
        with a seq property such as Bio.SeqRecord.SeqRecord.
    target : tuple[int, int]
        Start and stop position for target sequence.
    limit : str, optional
        This is the part in the 3'-end of each primer that has to
        anneal. The default is 16.
    automaton : ahocorasick.Automaton, optional
        Automaton made with the :func:`make_automaton`. The default is None.


    Returns
    -------
    list[tuple[int, int, int, int, int]]
        List of tuples (index_fp, position_fp, index_rp, position_rp, size).

    """

    begin, end = target

    assert begin is not None, "begin has to be set."
    assert end is not None, "end has to be set."

    amplicons = primer_pairs(
        seq,
        primer_list,
        short=end - begin,
        long=len(seq),
        limit=limit,
        automaton=automaton,
    )
    products = []

    for amplicon in amplicons:
        if contained(
            begin,
            end,
            amplicon.fposition,
            amplicon.rposition,
            len(seq),
            circular=seq.circular,
        ):
            products.append(amplicon)

    return sorted(products, key=attrgetter("size"))  # results sorted by size


def diff_primer_pairs(
    sequences: list[Dseqrecord] | tuple[Dseqrecord],
    primer_list: Sequence[Primer | None],
    short: int = 500,
    long: int = 1500,
    limit: int = 16,
    automaton: ahocorasick.Automaton = None,
    callback: Callable[[list], bool] = callback,
) -> tuple[tuple[Dseqrecord, int, int, int]]:
    """
    Primer pairs for diagnostic PCR.

    Given an iterable of sequences and a primer list, primers are selected that result in
    unique product sizes from each of the input sequences.

    Primers 1 and 2 both form PCR products from sequenceA and B below, but of
    different sizes. Primers 1 and 2 could be used to verify genetic modifications such
    as cloning an insert into a plasmid vector.

    ::

         1>              <2
        -------NNNNNNNNN----  sequenceA


         1>           <2
        -------XXXXX--------  sequenceB


    The callback function is used to return true or false for the PCR products. This score is
    meant to filter for PCR products that are likely to migrate to
    sufficiently distinct locations to be distinguishable on a typical agarose gel.

    Only products larger than `short` and smaller than `long` are returned.

    An example of the output for two sequences (Dseqrecord(-3308), Dseqrecord(-3613)).
    Primers 501 and 1806 would yield a 933 bp product with the 3308 bp sequence and the same
    primer pair would give 1212 bp with the 3613 bp sequence.

    A list of named 4-tuples is returned (Sequence, forward_primer, reverse_primer, size_bp),
    where each tuple has one entry for each sequence in the input argument.

    ::

        [
            ((Dseqrecord(-3308), 501, 1806, 933), (Dseqrecord(-3613), 501, 1806, 1212)),
        ]


    Parameters
    ----------
    sequences : list[Dseqrecord] | tuple[Dseqrecord]
        Target sequence to find primer annealing positions.
    primer_list : list[Primer] | tuple[Primer]
        This is a list of pydna.primer.Primer objects or any object
        with a seq property such as Bio.SeqRecord.SeqRecord.
    limit : str, optional
        This is the part in the 3'-end of each primer that has to
        anneal. The default is 16.
    short : int, optional
        Lower limit for the size of the PCR products. The default is 500.
    long : int, optional
        Upper limit for the size of the PCR products. The default is 1500.
    automaton : ahocorasick.Automaton, optional
        Automaton made with the :func:`make_automaton`. The default is None.
    callback : callable[[list], bool], optional
        A function accepting a list of integers and returning True or False.
        The default is callback.

    Returns
    -------
    list[tuple[Dseqrecord, int, int, int]]
        (Sequence, forward_primer, reverse_primer, size_bp)

    """

    automaton = automaton or make_automaton(primer_list, limit=limit)
    limit = automaton.get_stats()["longest_word"]
    primer_pair_dict = defaultdict(list)

    sequences = list(dict.fromkeys(sequences))
    number_of_sequences = len(sequences)

    # primer pairs from all sequences are collected in a defaultdict.
    #
    for seq in sequences:
        for pp in primer_pairs(
            seq, primer_list, short=short, long=long, limit=limit, automaton=automaton
        ):
            primer_pair_dict[frozenset((pp.fp, pp.rp))].append((pp, seq))

    result = []
    for k, v in primer_pair_dict.items():
        # verify one pcr product per sequence
        if len(v) == number_of_sequences:
            # we need the closest pair to see if the bands are likely to resolve
            closest_pair, diff = closest_pair_and_diff([pp.size for pp, _ in v])
            if callback(*closest_pair):
                # we add the diff here so we can sort the results by this value
                # we normally want a large difference in size
                result.append(
                    (
                        diff,
                        tuple(
                            primer_tuple(seq, pp.fp, pp.rp, pp.size) for pp, seq in v
                        ),
                    )
                )
    result.sort(reverse=True)
    return [pts for diff, pts in result]


def diff_primer_triplets(
    sequences: list[Dseqrecord] | tuple[Dseqrecord],
    primer_list: Sequence[Primer | None],
    limit: int = 16,
    short: int = 500,
    long: int = 1500,
    automaton: ahocorasick.Automaton = None,
    callback: Callable[[list], bool] = callback,
) -> tuple[tuple[tuple[Dseqrecord, int, int, int]]]:
    """
    Primer triplets for diagnostic PCR.

    Given a list of sequences and a primer list, primer triplets are selected that result in
    PCR products of different sizes from each of the input sequences.

    Primers 1, 2 and 3 form PCR products from sequenceA and B below, but of
    different sizes. Primer 1 binds both sequences while primers 2 and 3 bind one
    sequence each. This primer triplet could be used to verify genetic
    modifications.

    ::

         1>        <2
        -------AAAAAAAAA----  sequenceA

         1>     <3
        -------BBBBB--------  sequenceB

        ...

         1>      <3
        -------NNNNN--------  sequenceN-1

         2>      <3
        -AA----NNNNN--------  sequenceN

    The callback function is used to give a score for the PCR products. This score can
    be used to decide if a collection of PCR products are likely to migrate to distinct
    locations on a typical agarose gel.

    Only products larger than `short` and smaller than `long` are returned.

    An example of the output for two sequences = [Dseqrecord(-7664), Dseqrecord(-3613)].
    Primer pair 701, 700 would produce a 724 bp product with the 7664 bp sequence while
    the primer pair 701, 1564 would give a 1450 bp product with the 3613 bp sequence.

    ::

        [
            ((Dseqrecord(-7664), 701, 700, 724), (Dseqrecord(-3613), 701, 1564, 1450)),
         ]

    Parameters
    ----------
    sequences : list[Dseqrecord] | tuple[Dseqrecord]
        Target sequence to find primer annealing positions.
    primer_list : list[Primer] | tuple[Primer]
        This is a list of pydna.primer.Primer objects or any object
        with a seq property such as Bio.SeqRecord.SeqRecord.
    limit : str, optional
        This is the part in the 3'-end of each primer that has to
        anneal. The default is 16.
    short : int, optional
        Lower limit for the size of the PCR products. The default is 500.
    long : int, optional
        Upper limit for the size of the PCR products. The default is 2000.
    automaton : ahocorasick.Automaton, optional
        Automaton made with the :func:`make_automaton`. The default is None.
    callback : callable[[list], bool], optional
        A function accepting a list of integers and returning True or False.
        The default is callback.

    Returns
    -------
    list[tuple[Dseqrecord, int, int, int]]
        (Sequence, forward_primer, reverse_primer, size_bp)

    """

    automaton = automaton or make_automaton(primer_list, limit=limit)
    limit = automaton.get_stats()["longest_word"]
    # number_of_sequences = len(sequences)
    sequences = list(dict.fromkeys(sequences))
    sequence_set = set(sequences)
    # primer pairs for each sequence are collected together with the
    # primer pairs in a defaultdict(list) where the key of a frozenset
    # of the primer pair numbers.

    pairs = defaultdict(dict)

    for seq in sequences:
        primerpairs_for_one_sequence = primer_pairs(
            seq, primer_list, short=short, long=long, limit=limit, automaton=automaton
        )
        for pp in primerpairs_for_one_sequence:
            pairs[seq][pp.fp, pp.rp] = pp

    trios = defaultdict(list)

    # All collected primer pairs are compared and collected if they share a
    # primer (a trio). The trios are collected in another defaultdict(list)

    for one, two in combinations(pairs, 2):
        for p1, pt1 in pairs[one].items():
            for p2, pt2 in pairs[two].items():
                if len(trio := frozenset(p1 + p2)) == 3:
                    trios[trio].extend(((one, pt1), (two, pt2)))

    # We loop through all values in trios and filter:
    # at least one primer pair for each sequence
    trios = {
        trio: values
        for trio, values in trios.items()
        # filter dict for trios that have a PCR product from all sequences
        if set(s for s, _ in values) == sequence_set
    }

    # Each key, value pair in trios look like this:
    #
    #   trios[(1,2,3)] = [(seq1, pair1), (seq2, pair2), ... (seqN, pairN)]
    #
    # (1,2,3) is a triplet of primer numbers (integers)
    # seq1 .. seqN are sequences fed to the function
    # pair1 .. pairN is an amplicon_tuple(fp, rp, fposition, rposition, size)
    #
    # - visible differences between band sizes = callback true for all pairs if product sizes
    #
    # collect in intermediate_result.

    intermediate_result = []

    for seq_primer_tuples in trios.values():
        closest_pair, diff = closest_pair_and_diff(
            [pair.size for seq, pair in seq_primer_tuples]
        )
        if callback(*closest_pair):
            intermediate_result.append((diff, seq_primer_tuples))

    # The closest difference between any two product sizes is also collected
    # so we can sort results for primer combos that produce the largest
    # product differences:

    # sorting in-place for closest difference in reverse
    intermediate_result.sort(key=lambda item: item[0], reverse=True)

    # remove the closest difference as it is no lionger needed.
    intermediate_result = [b for a, b in intermediate_result]

    # Express result as a list of tuples of primer tuples

    primer_tuple_list = []

    for tuple_list in intermediate_result:
        primer_tuple_list.append(
            tuple(primer_tuple(s, at.fp, at.rp, at.size) for s, at in tuple_list)
        )

    return primer_tuple_list
