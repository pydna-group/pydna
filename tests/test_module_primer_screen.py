#!/usr/bin/env python

import os
import pytest
import tempfile
import pickle
import pathlib
import ahocorasick

from pydna.primer import Primer
from pydna.readers import read
from pydna.amplify import pcr
from pydna.parsers import parse_primers
from pydna.dseqrecord import Dseqrecord
from pydna.primer_screen import contained
from pydna.primer_screen import closest_pair_and_diff
from pydna.primer_screen import make_automaton
from pydna.primer_screen import forward_primers
from pydna.primer_screen import reverse_primers
from pydna.primer_screen import primer_pairs
from pydna.primer_screen import flanking_primer_pairs
from pydna.primer_screen import diff_primer_pairs
from pydna.primer_screen import diff_primer_triplets
from pydna.primer_screen import primer_tuple
from pydna.primer_screen import amplicon_tuple

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


def test_contained():
    assert contained(4, 7, 4, 8, 10, circular=True) is True
    assert contained(5, 8, 4, 8, 10, circular=True) is True
    assert contained(4, 8, 4, 8, 10, circular=True) is True

    assert contained(4, 7, 4, 8, 10, circular=False) is True
    assert contained(5, 8, 4, 8, 10, circular=False) is True
    assert contained(4, 8, 4, 8, 10, circular=False) is True

    assert contained(9, 10, 9, 1, 10, circular=True) is True
    assert contained(9, 1, 9, 1, 10, circular=True) is True
    assert contained(0, 1, 9, 1, 10, circular=True) is True

    assert contained(1, 5, 6, 38, 39, circular=True) is False
    assert contained(1, 5, 38, 6, 39, circular=True) is True

    assert contained(1, 5, 6, 38, 1, circular=True) is False

    # Overlaps but not contained
    assert contained(1, 5, 2, 8, 39, circular=True) is False
    assert contained(1, 5, 8, 2, 39, circular=True) is False


def test_closest_pair_and_diff():

    assert closest_pair_and_diff([1, 5, 7, 11, 19]) == ((5, 7), 2)
    assert closest_pair_and_diff([1, 5, 7, 17, 19]) == ((17, 19), 2)

    with pytest.raises(AssertionError):
        closest_pair_and_diff([5])


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

    s = str(primers[0].seq + "aaa" + primers[1].rc().seq)

    assert s == "CAGATGCGAAGTTAAGTGCGaaaCGTCAAGACTGTCAAGGA"

    # >51_TefTermFwd A.gos 20-mer
    # CAGATGCGAAGTTAAGTGCG
    # >82_MSW_fwd this primer sits inside the loxP of pUG6
    # TCCTTGACAGTCTTGACG

    # CAGATGCGAAGTTAAGTGCG
    # ||||||||||||||||||||
    # CAGATGCGAAGTTAAGTGCGaaaCGTCAAGACTGTCAAGGA
    # GTCTACGCTTCAATTCACGCtttGCAGTTCTGACAGTTCCT
    #                        ||||||||||||||||||
    #                        GCAGTTCTGACAGTTCCT

    assert primer_pairs(Dseqrecord(s), pl, short=0) == [
        amplicon_tuple(fp=51, rp=82, fposition=20, rposition=23, size=41)
    ]

    assert primer_pairs(Dseqrecord(s, circular=True), pl, short=0) == [
        amplicon_tuple(fp=51, rp=82, fposition=20, rposition=23, size=41)
    ]

    #      CAGATGCGAAGTTAAGTGCG>
    #      ||||||||||||||||||||
    # AAGGACAGATGCGAAGTTAAGTGCGaaaCGTCAAGACTGTC                       Circular template
    # TTCCTGTCTACGCTTCAATTCACGCtttGCAGTTCTGACAG
    #                                          AAGGACAGATGCGAAGTTAAGTGCGaaaCGTCAAGACTGTC
    #                                          TTCCTGTCTACGCTTCAATTCACGCtttGCAGTTCTGACAG
    #                             ||||||||||||||||||
    #                            <GCAGTTCTGACAGTTCCT

    assert primer_pairs(Dseqrecord(s[-5:] + s[:-5], circular=True), pl, short=0) == [
        amplicon_tuple(fp=51, rp=82, fposition=25, rposition=28, size=41)
    ]

    #                                     CAGATGCGAAGTTAAGTGCG>
    #                                     ||||||||||||||||||||
    # GCGAAGTTAAGTGCGaaaCGTCAAGACTGTCAAGGACAGAT                        Circular template
    # CGCTTCAATTCACGCtttGCAGTTCTGACAGTTCCTGTCTA
    #                                          GCGAAGTTAAGTGCGaaaCGTCAAGACTGTCAAGGACAGAT
    #                                          CGCTTCAATTCACGCtttGCAGTTCTGACAGTTCCTGTCTA
    #                   ||||||||||||||||||
    #                  <GCAGTTCTGACAGTTCCT

    assert primer_pairs(Dseqrecord(s[5:] + s[:5], circular=True), pl, short=0) == [
        amplicon_tuple(fp=51, rp=82, fposition=15, rposition=18, size=41)
    ]

    # amplicon_tuple(fp=51, rp=82, fposition=15, rposition=18, size=41) !=
    # amplicon_tuple(fp=51, rp=82, fposition=56, rposition=18, size=41)

    # CTCACTTGAAGTAATG
    # ||||||||||||||||
    # CTCACTTGAAGTAATGTATCGTGCACCTACCAAACCTCT
    # GAGTGAACTTCATTACATAGCACGTGGATGGTTTGGAGA
    #                        ||||||||||||||||
    #                        GTGGATGGTTTGGAGA

    f = Primer("CTCACTTGAAGTAATG", name="1")
    r = Primer("AGAGGTTTGGTAGGTG", name="2")
    t = Dseqrecord("CTCACTTGAAGTAATGTATCGTGCACCTACCAAACCTCT")
    automaton2 = make_automaton([f, r])

    assert primer_pairs(t, [f, r], automaton=automaton2, short=0) == [
        amplicon_tuple(fp=0, rp=1, fposition=16, rposition=23, size=39)
    ]

    #                        CTCACTTGAAGTAATG>
    #                        ||||||||||||||||
    # TATCGTGCACCTACCAAACCTCTCTCACTTGAAGTAATG   Linear template
    # ATAGCACGTGGATGGTTTGGAGAGAGTGAACTTCATTAC   No product
    #        ||||||||||||||||
    #        GTGGATGGTTTGGAGA

    t = Dseqrecord("TATCGTGCACCTACCAAACCTCTCTCACTTGAAGTAATG")

    assert primer_pairs(t, [f, r], automaton=automaton2, short=0) == []

    #                        CTCACTTGAAGTAATG>
    #                        ||||||||||||||||
    # TATCGTGCACCTACCAAACCTCTCTCACTTGAAGTAATG   Circular template
    # ATAGCACGTGGATGGTTTGGAGAGAGTGAACTTCATTAC   On product across the ori
    #        ||||||||||||||||
    #       <GTGGATGGTTTGGAGA

    t = Dseqrecord("TATCGTGCACCTACCAAACCTCTCTCACTTGAAGTAATG", circular=True)

    assert primer_pairs(t, [f, r], automaton=automaton2, short=0) == [
        amplicon_tuple(fp=0, rp=1, fposition=0, rposition=7, size=39)
    ]

    #                           CTCACTTGAAGTAATG>
    #                           ||||||||||||||||           Circular template
    # ATGTATCGTGCACCTACCAAACCTCTCTCACTTGAAGTA              On product across the ori
    # TACATAGCACGTGGATGGTTTGGAGAGAGTGAACTTCAT
    #                                        ATGTATCGTGCACCTACCAAACCTCTCTCACTTGAAGTA
    #                                        TACATAGCACGTGGATGGTTTGGAGAGAGTGAACTTCAT
    #           ||||||||||||||||
    #          <GTGGATGGTTTGGAGA

    t = Dseqrecord("ATGTATCGTGCACCTACCAAACCTCTCTCACTTGAAGTA", circular=True)

    assert primer_pairs(t, [f, r], automaton=automaton2, short=0) == [
        amplicon_tuple(fp=0, rp=1, fposition=3, rposition=10, size=39)
    ]


def test_flanking_primer_pairs():
    result = flanking_primer_pairs(kan, pl, target=(550, 1200), automaton=atm)

    answer = [
        amplicon_tuple(fp=701, rp=1564, fposition=534, rposition=1940, size=1450),
        amplicon_tuple(fp=701, rp=149, fposition=534, rposition=2306, size=1822),
    ]
    assert result == answer

    result = flanking_primer_pairs(pIL68, pl, target=(550, 1200), automaton=atm)

    answer = [
        amplicon_tuple(fp=1215, rp=594, fposition=266, rposition=1689, size=1474),
        amplicon_tuple(fp=1215, rp=607, fposition=266, rposition=2157, size=1946),
    ]

    # CTCACTTGAAGTAATG
    # ||||||||||||||||
    # CTCACTTGAAGTAATGTATCGTGCACCTACCAAACCTCT
    # GAGTGAACTTCATTACATAGCACGTGGATGGTTTGGAGA
    #                        ||||||||||||||||
    #                        GTGGATGGTTTGGAGA

    f = Primer("CTCACTTGAAGTAATG", name="1")
    r = Primer("AGAGGTTTGGTAGGTG", name="2")
    t = Dseqrecord("CTCACTTGAAGTAATGTATCGTGCACCTACCAAACCTCT")
    automaton2 = make_automaton([f, r])

    assert flanking_primer_pairs(t, [f, r], target=(16, 23), automaton=automaton2) == [
        amplicon_tuple(fp=0, rp=1, fposition=16, rposition=23, size=39)
    ]
    assert flanking_primer_pairs(t, [f, r], target=(17, 23), automaton=automaton2) == [
        amplicon_tuple(fp=0, rp=1, fposition=16, rposition=23, size=39)
    ]
    assert flanking_primer_pairs(t, [f, r], target=(16, 22), automaton=automaton2) == [
        amplicon_tuple(fp=0, rp=1, fposition=16, rposition=23, size=39)
    ]
    assert flanking_primer_pairs(t, [f, r], target=(15, 23), automaton=automaton2) == []
    assert flanking_primer_pairs(t, [f, r], target=(16, 24), automaton=automaton2) == []

    #                        CTCACTTGAAGTAATG>
    #                        ||||||||||||||||
    # TATCGTGCACCTACCAAACCTCTCTCACTTGAAGTAATG   Linear template
    # ATAGCACGTGGATGGTTTGGAGAGAGTGAACTTCATTAC   No product
    #        ||||||||||||||||
    #       <GTGGATGGTTTGGAGA

    t = Dseqrecord("TATCGTGCACCTACCAAACCTCTCTCACTTGAAGTAATG")

    assert flanking_primer_pairs(t, [f, r], target=(0, 7), automaton=automaton2) == []
    assert flanking_primer_pairs(t, [f, r], target=(1, 6), automaton=automaton2) == []

    #                        CTCACTTGAAGTAATG>
    #                        ||||||||||||||||
    # TATCGTGCACCTACCAAACCTCTCTCACTTGAAGTAATG   Circular template
    # ATAGCACGTGGATGGTTTGGAGAGAGTGAACTTCATTAC   On product across the ori
    #        ||||||||||||||||
    #       <GTGGATGGTTTGGAGA

    t = Dseqrecord("TATCGTGCACCTACCAAACCTCTCTCACTTGAAGTAATG", circular=True)

    assert flanking_primer_pairs(t, [f, r], target=(1, 6), automaton=automaton2) == [
        amplicon_tuple(fp=0, rp=1, fposition=0, rposition=7, size=39)
    ]

    #                       CTCACTTGAAGTAATG>
    #                       ||||||||||||||||
    # ATCGTGCACCTACCAAACCTCTCTCACTTGAAGTAATGT   Circular template
    # TAGCACGTGGATGGTTTGGAGAGAGTGAACTTCATTACA  On product across the ori
    #       ||||||||||||||||
    #      <GTGGATGGTTTGGAGA

    t = Dseqrecord("ATCGTGCACCTACCAAACCTCTCTCACTTGAAGTAATGT", circular=True)

    answer = amplicon_tuple(fp=0, rp=1, fposition=38, rposition=6, size=39)

    assert flanking_primer_pairs(t, [f, r], target=(1, 5), automaton=automaton2) == [
        answer
    ]
    assert flanking_primer_pairs(t, [f, r], target=(0, 5), automaton=automaton2) == [
        answer
    ]
    assert flanking_primer_pairs(t, [f, r], target=(39, 5), automaton=automaton2) == [
        answer
    ]
    assert flanking_primer_pairs(t, [f, r], target=(38, 6), automaton=automaton2) == [
        answer
    ]
    assert flanking_primer_pairs(t, [f, r], target=(37, 6), automaton=automaton2) == []
    assert flanking_primer_pairs(t, [f, r], target=(38, 7), automaton=automaton2) == []

    #                           CTCACTTGAAGTAATG>
    #                           ||||||||||||||||           Circular template
    # ATGTATCGTGCACCTACCAAACCTCTCTCACTTGAAGTA              On product across the ori
    # TACATAGCACGTGGATGGTTTGGAGAGAGTGAACTTCAT
    #                                        ATGTATCGTGCACCTACCAAACCTCTCTCACTTGAAGTA
    #                                        TACATAGCACGTGGATGGTTTGGAGAGAGTGAACTTCAT
    #           ||||||||||||||||
    #          <GTGGATGGTTTGGAGA

    t = Dseqrecord("ATGTATCGTGCACCTACCAAACCTCTCTCACTTGAAGTA", circular=True)

    assert flanking_primer_pairs(t, [f, r], target=(0, 6), automaton=automaton2) == []
    assert flanking_primer_pairs(t, [f, r], target=(3, 9), automaton=automaton2) == [
        amplicon_tuple(fp=0, rp=1, fposition=3, rposition=10, size=39)
    ]


def test_diff_primer_pairs():

    results = diff_primer_pairs((nat, kan), pl, automaton=atm)

    assert results == [
        (
            primer_tuple(seq=nat, fp=82, rp=149, size=944),
            primer_tuple(seq=kan, fp=82, rp=149, size=1181),
        ),
    ]

    f = Primer("CTCACTTGAAGTAATG", name="1")
    r = Primer("AGAGGTTTGGTAGGTG", name="2")

    t1 = Dseqrecord("CTCACTTGAAGTAATGtaTCGTGCACCTACCAAACCTCT")
    t2 = Dseqrecord("CTCACTTGAAGTAATGccTCGTGCACCTACCAAACCTCT")
    automaton2 = make_automaton([f, r])

    # Returns emty list, since the two products are identical, so Callback is false
    assert diff_primer_pairs((t1, t2), [f, r], automaton=automaton2, short=0) == []

    s = primers[0] + "aaa" + primers[1].rc()

    s = "CAGATGCGAAGTTAAGTGCGaaaCGTCAAGACTGTCAAGGA"

    assert primer_pairs(Dseqrecord(s), pl, short=0) == [
        amplicon_tuple(fp=51, rp=82, fposition=20, rposition=23, size=41)
    ]


def test_diff_primer_triplets():

    triplets1 = diff_primer_triplets((wt, kan), pl, automaton=atm)
    assert len(triplets1) == 1
    results = [
        (
            primer_tuple(seq=wt, fp=701, rp=700, size=724),
            primer_tuple(seq=kan, fp=701, rp=1564, size=1450),
        ),
    ]
    assert triplets1 == results

    triplets2 = diff_primer_triplets([pIL68, pIL75], pl)
    assert len(triplets2) == 4
    answer = [
        (
            primer_tuple(seq=pIL68, fp=1215, rp=594, size=1474),
            primer_tuple(seq=pIL75, fp=51, rp=594, size=548),
        ),
        (
            primer_tuple(seq=pIL68, fp=595, rp=607, size=546),
            primer_tuple(seq=pIL75, fp=255, rp=607, size=1467),
        ),
        (
            primer_tuple(seq=pIL68, fp=1215, rp=594, size=1474),
            primer_tuple(seq=pIL75, fp=255, rp=594, size=1005),
        ),
        (
            primer_tuple(seq=pIL68, fp=595, rp=607, size=546),
            primer_tuple(seq=pIL75, fp=51, rp=607, size=1010),
        ),
    ]

    assert triplets2 == answer

    assert len(pcr(pl[1215], pl[594], pIL68)) == 1474
    assert len(pcr(pl[51], pl[594], pIL75)) == 548

    with pytest.raises(ValueError, match="No PCR product!"):
        pcr(pl[51], pl[594], pIL68)

    with pytest.raises(ValueError, match="No PCR product!"):
        pcr(pl[1215], pl[594], pIL75)

    f = Primer("CTCACTTGAAGTAATG", name="1")
    r = Primer("AGAGGTTTGGTAGGTG", name="2")
    r2 = Primer("AAGCATAACACACGTA", name="2a")
    t1 = Dseqrecord("CTCACTTGAAGTAATGtaTCGTGCACCTACCAAACCTCT")
    t3 = Dseqrecord("CTCACTTGAAGTAATGccAGCTATACGTGTGTTATGCTT")
    automaton3 = make_automaton([f, r, r2])

    assert (
        diff_primer_triplets((t1, t3), [f, r, r2], automaton=automaton3, short=0) == []
    )
