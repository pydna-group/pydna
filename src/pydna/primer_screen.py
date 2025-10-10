# -*- coding: utf-8 -*-
# Aho–Corasick algorithm
# https://github.com/WojciechMula/pyahocorasick
# https://pypi.org/project/pyahocorasick/

import ahocorasick
from pydna.parsers import parse_primers
from pydna.readers import read
from itertools import product

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


def primer_screen(first_sequence, second_sequence, primer_list):

    p = "/home/bjorn/myvault/PRIMERS.md"

    lst = parse_primers(p)[::-1]

    del lst[582]

    limit = 16

    automaton = ahocorasick.Automaton()

    for idx, key in enumerate(lst):
        footprint = str(key.seq)[-limit:].upper()
        if len(footprint) >= limit:
            for anchor in expand_iupac_to_dna(footprint):
                if anchor:
                    automaton.add_word(anchor, (idx, anchor))

    automaton.make_automaton()

    paths = [
        "/home/bjorn/Desktop/pydna/ahocorasick/fas2::KanMX4_locus.gb",
        "/home/bjorn/Desktop/pydna/ahocorasick/fas2::NatMX4_locus.gb",
        "/home/bjorn/Desktop/pydna/ahocorasick/FAS2_S288C_wild-type_locus.gb",
    ]

    sequences = []

    for path in paths:
        sequences.append(read(path))

    del sequences[1]

    forward_primers = {}
    reverse_primers = {}

    for seq in sequences:

        haystack = str(seq.seq).upper()

        ln = len(haystack)

        fps = dict()

        for end_index, (insert_order, original_value) in automaton.iter(haystack):

            start_index = end_index - limit + 1

            fps[insert_order] = start_index + limit

        forward_primers[seq] = fps

        haystack_rc = str(seq.seq.rc()).upper()

        rps = dict()

        for end_index, (insert_order, original_value) in automaton.iter(haystack_rc):

            start_index = end_index - limit + 1

            rps[insert_order] = ln - (start_index + limit)

            reverse_primers[seq] = rps

    data = forward_primers

    common_fps = set.intersection(*(set(inner.keys()) for inner in data.values()))

    data = reverse_primers

    common_rps = set.intersection(*(set(inner.keys()) for inner in data.values()))

    print(common_fps)
    print(common_rps)


if __name__ == "__main__":

    from pydna_utils.myprimers import PrimerList

    pl = PrimerList()

    fs = read("/home/bjorn/Desktop/pydna/ahocorasick/FAS2_S288C_wild-type_locus.gb")
    ss = read("/home/bjorn/Desktop/pydna/ahocorasick/fas2::KanMX4_locus.gb")

    primer_screen(fs, ss, pl)
