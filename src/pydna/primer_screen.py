#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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


def closest_diff(nums: list[int]) -> int:
    """
    Smallest difference between two consecutive integers in a sorted list.

    Given a list of integers eg. 1, 5, 7, 11, 19, return the smallest
    absolute difference, in this case 7-5 = 2.

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

    # if no automaton is given, we make one.
    automaton = automaton or make_automaton(primer_list, limit=limit)

    # The limit is taken from automaton stats.
    limit = automaton.get_stats()["longest_word"]

    # A defaultdict of lists is used to collect primer locations since
    # different primers can anneal in the same place.
    fps = defaultdict(list)

    for end_index, ids in automaton.iter(str(seq.seq).upper()):
        for i in ids:
            fps[i].append(end_index + 1)

    return dict(fps)


def reverse_primers(
    seq: Dseqrecord,
    primer_list: list[Primer] | tuple[Primer],
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
    # if no automaton is given, we make one.
    automaton = automaton or make_automaton(primer_list, limit=limit)

    # The limit is taken from automaton stats.
    # If the automaton is given, the limit argument will be ignored.
    limit = automaton.get_stats()["longest_word"]

    # A defaultdict of lists is used to collect primer locations since
    # different primers can anneal in the same place.
    rps = defaultdict(list)
    ln = len(seq)

    # We use the reverse complement of the sequence instead of taking the
    # reverse complement of each primer.
    for end_index, ids in automaton.iter(str(seq.seq.reverse_complement()).upper()):
        for i in ids:
            rps[i].append(ln - (end_index + 1))

    return dict(rps)


def primer_pairs(
    seq: Dseqrecord,
    primer_list: list[Primer] | tuple[Primer],
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
    automaton = automaton or make_automaton(primer_list, limit=limit)
    limit = automaton.get_stats()["longest_word"]

    # Unique forward primers are collected
    fps = {
        fp: pos[0]
        for fp, pos in forward_primers(
            seq, primer_list, limit=limit, automaton=automaton
        ).items()
        if len(pos) == 1
    }

    # Unique reverse primers are collected
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
            # We calculate the size of a potential PCR product
            size = len(primer_list[fp]) + rposition - fposition + len(primer_list[rp])
            # If the size falls within long and short, the data is kept.
            if short <= size <= long and fposition <= rposition:
                products.append(amplicon_tuple(fp, rp, fposition, rposition, size))
    return products


def flanking_primer_pairs(
    seq: Dseqrecord,
    primer_list: list[Primer] | tuple[Primer],
    target: tuple[int, int],
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

    automaton = automaton or make_automaton(primer_list, limit=limit)
    limit = automaton.get_stats()["longest_word"]

    begin, end = target

    assert begin < end, "begin has to be smaller than end."

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
        if amplicon.fposition >= begin and end <= amplicon.rposition:
            products.append(amplicon)

    return products[::-1]


def diff_primer_pairs(
    sequences: list[Dseqrecord] | tuple[Dseqrecord],
    primer_list: list[Primer] | tuple[Primer],
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
    sequences: list[Dseqrecord] | tuple[Dseqrecord],
    primer_list: list[Primer] | tuple[Primer],
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
        -------NNNNNNNNN----  sequenceA

         1>     <3
        -------XXXXX--------  sequenceB



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
    number_of_sequences = len(sequences)
    pp = {}
    # pp = { seq1: [(a,b,c,d,e), ...], seq2: [(i,j,k,l,m), ... ]}

    # All primer pairs for each sequence are collected.
    for seq in sequences:
        pp[seq] = primer_pairs(
            seq, primer_list, short=short, long=long, limit=limit, automaton=automaton
        )

    # We count all the times a specific pair occurs
    pair_counter = Counter()

    for seq, tuples in pp.items():
        for t in tuples:
            pair = frozenset(t[:2])  # first two integers, unordered
            pair_counter[pair] += 1

    # Pick pairs that appear more than once.
    pairs_to_remove = {pair for pair, count in pair_counter.items() if count > 1}

    # Remove pairs that appear more than once.
    for seq in pp:
        pp[seq] = [t for t in pp[seq] if frozenset(t[:2]) not in pairs_to_remove]

    primertrios = defaultdict(dict)

    for seq1, seq2 in combinations(sequences, 2):
        for fp1, rp1, *_, size1 in pp[seq1]:
            for fp2, rp2, *_, size2 in pp[seq2]:
                primertrio = frozenset((fp1, rp1, fp2, rp2))
                if len(primertrio) == 3 and callback(size1, size2):
                    if primertrios[primertrio]:
                        del primertrios[primertrio]
                    else:
                        primertrios[primertrio][size1] = (fp1, rp1, seq1)
                        primertrios[primertrio][size2] = (fp2, rp2, seq2)

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

    result.sort(key=lambda item: item[0], reverse=True)
    return [b for a, b in result]
