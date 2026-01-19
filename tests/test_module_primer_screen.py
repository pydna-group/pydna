#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest
import tempfile
import pickle
import pathlib
import ahocorasick

from pydna.readers import read
from pydna.amplify import pcr
from pydna.parsers import parse_primers

from pydna.primer_screen import make_automaton
from pydna.primer_screen import forward_primers
from pydna.primer_screen import reverse_primers
from pydna.primer_screen import primer_pairs
from pydna.primer_screen import flanking_primer_pairs
from pydna.primer_screen import diff_primer_pairs
from pydna.primer_screen import diff_primer_triplets
from pydna.primer_screen import primer_tuple
from pydna.primer_screen import amplicon_tuple
from pydna.primer_screen import expand_iupac_to_dna


test_files = pathlib.Path(os.path.join(os.path.dirname(__file__)))


primers = parse_primers(
    """
>51_TefTermFwd A.gos 20-mer
CAGATGCGAAGTTAAGTGCG

>82_MSW_fwd this primer sits inside the loxP of pUG6
TCCTTGACAGTCTTGACG

>149_MX4rev (25-mer) Primer in the Ashbya gossypi TEF terminator in the reverse direction. NEW (2021-01-05)
ACAAATGACAAGTTCTTGAAAACAA

>255_kanC (22-mer)
TGATTTTGATGACGAGCGTAAT

>594_pRS306_rev 24-mer
tcgctttcttcccttcctttctcg

>595_pRS306_fwd
taaagggagcccccgatttagagc

>607_XDHcasrv (28-mer)
aaacagctatgaccatgattacgccaag

>700_sc_fas2-B1: (25-mer)
ATTTCTCTATGTAAAGACAGAGCAG

>701_sc_fas2-A1: (25-mer)
CTATATTTCTATTCTATCCGAACTC

>1215_URA3r 27-mer
CGGTTTCTTTGAAATTTTTTTGATTCG

>1564_KANMX_rev
CACTCGCATCAACCAAACC
"""
)

pl = [None for i in range(1565)]

for primer in primers:
    number, rest = primer.id.split("_", maxsplit=1)
    pl[int(number)] = primer

wt = read(test_files / "FAS2_S288C_wild-type_locus.gb")
nat = read(test_files / "fas2__NatMX4_locus.gb")
kan = read(test_files / "fas2__KanMX4_locus.gb")

pIL68 = read(test_files / "pIL68.gb")
pIL75 = read(test_files / "pIL75.gb")


atm = None


def test_automaton():

    atm = make_automaton(pl)
    # The NamedTemporaryFile(delete=False) may be necessary for the tests
    # to pass on windows.
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        # Save automaton to temporary file
        atm.save(tmp.name, pickle.dumps)
        # tmp.flush()

        # Load it back from the same file
        atm2 = ahocorasick.load(tmp.name, pickle.loads)

    # Verify that loading worked
    assert atm2 is not None
    assert [x for x in atm2.keys()] == [x for x in atm.keys()]
    assert [x for x in atm2.values()] == [x for x in atm.values()]


def test_forward_primers():
    result = forward_primers(kan, pl, automaton=atm)
    assert result == {701: [534], 82: [1168], 255: [1979], 51: [2434]}


def test_reverse_primers():
    result = reverse_primers(wt, pl, automaton=atm)
    assert result == {700: [1208]}


def test_primer_pairs():
    result = primer_pairs(kan, pl, automaton=atm, short=0)
    answer = [
        amplicon_tuple(fp=701, rp=149, fposition=534, rposition=2306, size=1822),
        amplicon_tuple(fp=701, rp=1564, fposition=534, rposition=1940, size=1450),
        amplicon_tuple(fp=82, rp=149, fposition=1168, rposition=2306, size=1181),
        amplicon_tuple(fp=82, rp=1564, fposition=1168, rposition=1940, size=809),
        amplicon_tuple(fp=255, rp=149, fposition=1979, rposition=2306, size=374),
    ]
    assert result == answer


def test_flanking_primer_pairs():
    result = flanking_primer_pairs(kan, pl, target=(550, 1200), automaton=atm)

    answer = [
        amplicon_tuple(fp=82, rp=1564, fposition=1168, rposition=1940, size=809),
        amplicon_tuple(fp=82, rp=149, fposition=1168, rposition=2306, size=1181),
    ]
    assert result == answer


def test_diff_primer_pairs():

    results = diff_primer_pairs((nat, kan), pl, automaton=atm)

    assert results == [
        (
            primer_tuple(seq=nat, fp=82, rp=149, size=944),
            primer_tuple(seq=kan, fp=82, rp=149, size=1181),
        ),
    ]


def test_diff_primer_triplets_1():

    results = diff_primer_triplets((wt, kan), pl, automaton=atm)

    assert results == [
        (
            primer_tuple(seq=wt, fp=701, rp=700, size=724),
            primer_tuple(seq=kan, fp=701, rp=1564, size=1450),
        ),
    ]


def test_diff_primer_triplets_2():

    triplets = diff_primer_triplets([pIL68, pIL75], pl)

    assert len(triplets) == 2

    answer = [
        (
            primer_tuple(seq=pIL68, fp=1215, rp=594, size=1474),
            primer_tuple(seq=pIL75, fp=51, rp=594, size=548),
        ),
        (
            primer_tuple(seq=pIL68, fp=1215, rp=594, size=1474),
            primer_tuple(seq=pIL75, fp=255, rp=594, size=1005),
        ),
    ]

    assert triplets == answer

    assert len(pcr(pl[1215], pl[594], pIL68)) == 1474
    assert len(pcr(pl[51], pl[594], pIL75)) == 548

    with pytest.raises(ValueError, match="No PCR product!"):
        pcr(pl[51], pl[594], pIL68)

    with pytest.raises(ValueError, match="No PCR product!"):
        pcr(pl[1215], pl[594], pIL75)
