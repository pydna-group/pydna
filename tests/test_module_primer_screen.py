#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import pickle
import ahocorasick
from pydna.readers import read
from pydna_utils.myprimers import PrimerList
from pydna.amplify import pcr

from pydna.primer_screen import make_automaton
from pydna.primer_screen import forward_primers
from pydna.primer_screen import reverse_primers
from pydna.primer_screen import primer_pairs
from pydna.primer_screen import flanking_primer_pairs
from pydna.primer_screen import diff_primer_pairs
from pydna.primer_screen import diff_primer_triplets


from pydna.parsers import parse_primers


primers = parse_primers ("""
>82_MSW_fwd this primer sits inside the loxP of pUG6
TCCTTGACAGTCTTGACG

>149_MX4rev (25-mer) Primer in the Ashbya gossypi TEF terminator in the reverse direction. NEW (2021-01-05)
ACAAATGACAAGTTCTTGAAAACAA

>700_sc_fas2-B1: (25-mer)
ATTTCTCTATGTAAAGACAGAGCAG

>701_sc_fas2-A1: (25-mer)
CTATATTTCTATTCTATCCGAACTC

>1564_KANMX_rev
CACTCGCATCAACCAAACC
""")

pl = [None for i in range(1565)]

for primer in primers:
    number, rest = primer.id.split("_", maxsplit=1)
    pl[int(number)] = primer

wt = read("FAS2_S288C_wild-type_locus.gb")
nat = read("fas2::NatMX4_locus.gb")
kan = read("fas2::KanMX4_locus.gb")

atm = None

def test_automaton():

    atm = make_automaton(pl)

    atm.save("atm.automaton", pickle.dumps)

    atm = ahocorasick.load("atm.automaton", pickle.loads)


def test_diff_primer_pairs():

    results1 = diff_primer_pairs((nat, kan), pl, automaton=atm)

    assert results1 == (
        ((nat, 82, 149, 944), (kan, 82, 149, 1181)),
        )


def test_diff_primer_triplets():
    results2 = diff_primer_triplets((wt, kan), pl, automaton=atm)

    assert results2 == (
        ((wt, 701, 700, 724), (kan, 701, 1564, 1450)),
        )


# result = flanking_primer_pairs(s, pl, (1001, 5661), automaton=atm)

# PCR products with primers 700 701 1564 sizes 724 bp (wt) 1450 (deletion locus)
# from pydna.amplify import pcr

# for k, v in result.items():
#     primers = [pl[p] for p in k]
#     product1 = pcr(primers, s, limit=16)
#     product2 = pcr(primers, ss, limit=16)
#     assert {s: len(product1), ss: len(product2)} == v
#     print(len(product1), len(product2))


# print(forward_primers(s, pl, automaton=atm))
# print(reverse_primers(s, pl, automaton=atm))
# print(primer_pairs(s, pl, automaton=atm, short=79))

# pl = [pl[700],pl[701], pl[1564]]

# common_forward_primers([wt, nat, kan], pl, automaton=atm)
# common_reverse_primers([wt, nat, kan], pl, automaton=atm)

# unique_forward_primers([wt, nat, kan], pl, automaton=atm)
# unique_reverse_primers([wt, nat, kan], pl, automaton=atm)



# s = set()

# for p1, p2, (size1, size2) in results:
#     product1 = pcr(pl[p1], pl[p2], nat, limit=16)
#     product2 = pcr(pl[p1], pl[p2], kan, limit=16)
#     assert len(product1)==size1 and len(product2)==size2
#     s.add((p1, p2))

# assert len(s) == len(results)


    # with pytest.raises(ValueError):
    #     assembly_fragments(frags, 20)
