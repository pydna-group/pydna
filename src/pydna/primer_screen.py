#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Aho–Corasick algorithm
https://github.com/WojciechMula/pyahocorasick
https://pypi.org/project/pyahocorasick/
"""

from itertools import product
from itertools import permutations
from collections import defaultdict
from collections.abc import Callable

from pydna.readers import read
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer

import ahocorasick

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


def expand_iupac_to_dna(seq: str) -> list[str]:
    """
    Expand an extended-IUPAC DNA string (ACGTURYSWKMBDHVN) into all possible
    DNA strings using only A/C/G/T. Returns a list of strings.

    Example:
        expand_iupac_to_dna("ATNG") -> ["ATAG","ATCG","ATGG","ATTG"]
    """

    choices_per_pos = [_IUPAC_TO_DNA[ch] for ch in seq.upper()]
    # Cartesian product of all position choices
    return ["".join(tup) for tup in product(*choices_per_pos)]


def make_automaton(primer_list, limit=16):
    automaton = ahocorasick.Automaton()
    tempdict = defaultdict(list)
    for idx, key in enumerate(primer_list):
        footprint = str(key.seq)[-limit:].upper()
        if len(footprint) >= limit:
            tempdict[footprint].append(idx)
            for anchor in expand_iupac_to_dna(footprint):
                automaton.add_word(anchor, (*tempdict[footprint],))
    automaton.make_automaton()
    return automaton


def callback(s: list[int]):
    return max(s) - min(s) >= 0.2 * max(s)


def primer_screen(
    seq: Dseqrecord, primer_list: list[Primer] | tuple[Primer], limit: int = 16
):

    automaton = make_automaton(primer_list, limit=limit)
    forward_primers = {}
    reverse_primers = {}

    fps = defaultdict(list)

    # fps = {primer1: [psition1, position2, ...],
    #        primer2: [psition1, position2, ...], ... }
    # fps is a dict of lists where the key is the index of the primer
    # in the list used to make the automaton.

    for end_index, ids in automaton.iter(str(seq.seq).upper()):
        for i in ids:
            fps[i].append(end_index + 1)

    # automaton.iter return positions sorted in ascending order
    #
    # forward_primers is a dict with each input sequence as key
    # and the fps dict (see above as value)

    forward_primers[seq] = fps

    # The rps is a dict for reverse primers, similar to the fps

    rps = defaultdict(list)
    ln = len(seq)

    for end_index, ids in automaton.iter(str(seq.seq.rc()).upper()):
        for i in ids:
            rps[i].append(ln - (end_index + 1))

    reverse_primers[seq] = rps

    return dict(fps), dict(rps)


def primer_pairs(
    seq: Dseqrecord,
    primer_list: list[Primer] | tuple[Primer],
    limit: int = 16,
    low: int = 500,
    high: int = 2000,
):

    fps, rps = primer_screen(seq, primer_list)
    products = []

    for fp, fpositions in fps.items():
        for fposition in fpositions:
            for rp, rpositions in rps.items():
                for rposition in rpositions:
                    size = (
                        len(primer_list[fp])
                        + rposition
                        - fposition
                        + len(primer_list[rp])
                    )
                    if low <= size <= high:
                        products.append(((fp, fposition), (rp, rposition), size))
    return products


def diff_primer_pairs(
    sequences: list[Dseqrecord] | tuple[Dseqrecord],
    primer_list: list[Primer] | tuple[Primer],
    limit: int = 16,
    low: int = 500,
    high: int = 2000,
    callback: Callable[[list], bool] = callback,
):

    product_dict = {}
    number_of_sequences = len(sequences)

    for seq in sequences:

        products = [
            (fp, rp, size)
            for (fp, _), (rp, _), size in primer_pairs(
                seq, primer_list, limit, low, high
            )
        ]

        product_dict[seq] = defaultdict(list)

        for fp, rp, size in products:
            product_dict[seq][frozenset((fp, rp))].append((fp, rp, size))

    # all primer pairs that are common between the sequences
    common = set.intersection(*(set(inner.keys()) for inner in product_dict.values()))
    result = []

    # The length of each pcr product is collected across sequences.
    # The ones with different lengths as one distinct length per sequence
    # are kept.
    for primer_pair_set in common:
        sizes = defaultdict(list)
        for seq, d in product_dict.items():
            for fp, rp, size in d[primer_pair_set]:
                sizes[size].append(seq)
        # Callback function returns True or False and meant to
        # filter for size differences.
        if len(sizes) >= number_of_sequences and callback(sizes.keys()):
            result.append((max(sizes.keys()) - min(sizes.keys()), fp, rp, dict(sizes)))

    result.sort(reverse=True)

    return tuple((fp, rp, dct) for (s, fp, rp, dct) in result)


def diff_primer_trios(
    sequences: list[Dseqrecord] | tuple[Dseqrecord],
    primer_list: list[Primer] | tuple[Primer],
    limit: int = 16,
    low: int = 300,
    high: int = 2500,
    callback: Callable[[list], bool] = callback,
):

    product_dict = {}
    fprimers = {}
    rprimers = {}

    for seq in sequences:

        products = primer_pairs(seq, primer_list, limit, low, high)

        product_dict[seq] = defaultdict(list)

        for (fp, x), (rp, y), size in products:
            product_dict[seq][frozenset((fp, rp))].append((fp, rp, size, x, y))

        fprimers[seq], rprimers[seq] = primer_screen(seq, pl)

    # all primer pairs that are common between the sequences
    common = set.intersection(*(set(inner.keys()) for inner in product_dict.values()))
    result = []

    # The length of each pcr product is collected across sequences.
    # The ones with identical lengths are kept.

    for primer_pair_set in common:
        sizes = defaultdict(list)
        for seq, d in product_dict.items():
            for fp, rp, size, x, y in d[primer_pair_set]:
                sizes[size].append(seq)
        if len(sizes) == 1:
            result.append((fp, rp, size, x, y))

    for seq1, seq2 in permutations(sequences, 2):
        fps = fprimers[seq1].keys() - fprimers[seq2].keys()
        rps = rprimers[seq1].keys() - rprimers[seq2].keys()
        fps, rps
        # for fp, rp, size, x, y in results:
        #     for

    # breakpoint()


if __name__ == "__main__":

    from pydna_utils.myprimers import PrimerList

    # from pydna.amplify import pcr

    pl = PrimerList()

    s = read("/home/bjorn/Desktop/pydna/ahocorasick/FAS2_S288C_wild-type_locus.gb")

    primer_screen(s, pl)
    primer_pairs(s, pl)

    fs = read("/home/bjorn/Desktop/pydna/ahocorasick/fas2::NatMX4_locus.gb")
    ss = read("/home/bjorn/Desktop/pydna/ahocorasick/fas2::KanMX4_locus.gb")

    # result = diff_primer_pairs((fs, ss), pl)

    result = diff_primer_trios((s, fs), pl)

    # print(fp, rp, tuple((v, k) for k,v in sizes.items()))
    # [((1775, 1415), {1446: [File(.)(-3308)], 1751: [File(.)(-3613)]}, 305),
    #  ((701, 1421), {1858: [File(.)(-3308)], 2163: [File(.)(-3613)]}, 305),
    #  ((1477, 1447), {2216: [File(.)(-3308)], 2521: [File(.)(-3613)]}, 305),
    # fp = pl[701]
    # rp = pl[1421]
    # rp = pl[473]
    # rp = pl[189]
    # >1775_fw_ERG10_KanMX_del
    # GCGCCATATCATATATATTTATACAGATTAGACGTACTCAAAATGcagctgaagcttcgtacgc
    # >1415_ScFas2.3_loxP_rv
    # Aaacgatagaaaaataacaaagtaattactattatgtctgataaa
    # pcr(fp,rp,fs, limit=16).figure()
    # pcr(fp,rp,ss, limit=16).figure()
