#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Aho–Corasick algorithm
# https://github.com/WojciechMula/pyahocorasick
# https://pypi.org/project/pyahocorasick/

from itertools import product
from collections import defaultdict

import ahocorasick
from pydna.readers import read


# from Bio.SeqIO.FastaIO import SimpleFastaParser
# with open(p) as handle:
#     for values in SimpleFastaParser(handle):
#         print(values)

# primer_bind     1..35
#                 /PCR_conditions="primer
#                 sequence:GATCGGCCGGATCCAAATGACTGAATTCAAGGCCG"
#                 /locus_tag="1_5CYC1clone"
#                 /label="1_5CYC1clone"
#                 /ApEinfo_label="1_5CYC1clone"

# Map IUPAC -> sets over {A,C,G,T}. U is treated as T.

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


def iter_iupac_expansions(seq: str):
    """
    Yield DNA strings one-by-one instead of building a full list.
    """
    for tup in product(*(_IUPAC_TO_DNA[ch] for ch in seq.upper())):
        yield "".join(tup)


def primer_screen(
    first_sequence, second_sequence, primer_list, limit=16, low=300, high=2500
):

    automaton = ahocorasick.Automaton()

    for idx, key in enumerate(primer_list):
        footprint = str(key.seq)[-limit:].upper()
        if len(footprint) >= limit:
            for anchor in expand_iupac_to_dna(footprint):
                if anchor:
                    automaton.add_word(anchor, (idx, anchor))

    automaton.make_automaton()

    forward_primers = {}
    reverse_primers = {}
    products = {}
    product_dict = {}

    for seq in (first_sequence, second_sequence):

        haystack = str(seq.seq).upper()
        fps = dict()

        for end_index, (insert_order, original_value) in automaton.iter(haystack):

            start_index = end_index - limit + 1
            fps[insert_order] = start_index + limit

        forward_primers[seq] = fps  # sorted by position

        haystack_rc = str(seq.seq.rc()).upper()
        rps = dict()
        ln = len(haystack)

        for end_index, (insert_order, original_value) in automaton.iter(haystack_rc):

            start_index = end_index - limit + 1
            rps[insert_order] = ln - (start_index + limit)

        reverse_primers[seq] = rps  # reverse sorted by position

        products[seq] = [
            ((f_key, f_val), (r_key, r_val), r_val - f_val)
            for f_key, f_val in fps.items()
            for r_key, r_val in rps.items()
            if low <= r_val - f_val <= high
        ]

        product_dict[seq] = defaultdict(list)

        for (fprm, fpos), (rpro, rpos), diff in products[seq]:
            product_dict[seq][fprm, rpro].append(((fprm, fpos), (rpro, rpos), diff))

    # for key in forward_primers.keys():

    #     fps = {k:data[key][k] for k in common_fps}

    #     fps = dict(sorted(fps.items(), key=lambda item: item[1]))

    # for key in reverse_primers.keys():

    #     rps = {k:data[key][k] for k in common_rps}

    #     rps = dict(sorted(rps.items(), key=lambda item: item[1]))

    # data = forward_primers

    # common_fps = set.intersection(*(set(inner.keys()) for inner in data.values()))

    # data = reverse_primers

    # common_rps = set.intersection(*(set(inner.keys()) for inner in data.values()))


if __name__ == "__main__":

    from pydna_utils.myprimers import PrimerList

    pl = PrimerList()

    fs = read("/home/bjorn/Desktop/pydna/ahocorasick/FAS2_S288C_wild-type_locus.gb")
    ss = read("/home/bjorn/Desktop/pydna/ahocorasick/fas2::KanMX4_locus.gb")

    primer_screen(fs, ss, pl)
