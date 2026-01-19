#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa: B950

from Bio.Restriction import (
    AatII,
    AjiI,
    AgeI,
    EcoRV,
    ZraI,
    SalI,
    EcoRI,
    RgaI,
    BsaI,
    BsrI,
    DrdI,
    HindIII,
    FauI,
    NotI,
)
from pydna.amplify import pcr
from pydna.dseq import Dseq
from pydna.readers import read
import pydna.assembly2 as assembly
from Bio.SeqFeature import ExactPosition, FeatureLocation, SeqFeature, SimpleLocation
from pydna.dseqrecord import Dseqrecord
from pydna.parsers import parse
from pydna.utils import eq
import pytest
from Bio.Seq import reverse_complement
from pydna.primer import Primer
import os

test_files = os.path.join(os.path.dirname(__file__))


def dseqrecord_list_to_dseq_list(dseqrecord_list):
    """Utility function for tests, to convert records to sequences and check equality."""
    return [i.seq for i in dseqrecord_list]


example_fragments = (
    Dseqrecord("AacgatCAtgctcc", name="a"),
    Dseqrecord("TtgctccTAAattctgc", name="b"),
    Dseqrecord("CattctgcGAGGacgatG", name="c"),
)


linear_results = (
    Dseqrecord("AacgatCAtgctccTAAattctgcGAGGacgatG", name="abc"),
    Dseqrecord("ggagcaTGatcgtCCTCgcagaatG", name="ac_rc"),
    Dseqrecord("AacgatG", name="ac"),
)

circular_results = (
    Dseqrecord("acgatCAtgctccTAAattctgcGAGG", name="abc", circular=True),
)


def test_built():

    asm = assembly.Assembly(example_fragments, limit=5)
    lin = sorted(asm.assemble_linear(), key=len)
    crc = asm.assemble_circular()

    assert [dseqr.seq for dseqr in lin] == [
        dseqr.seq for dseqr in sorted(linear_results, key=len)
    ]
    assert [c.seq.seguid() for c in crc] == [c.seq.seguid() for c in circular_results]


def test_new_assembly():

    #                   0000000000111111111222222222233333333334444444444555555555566
    #                   0123456780123456789012345678901234567890123456789012345678901
    #                   acgatgctatactggCCCCCtgtgctgtgctctaTTTTTtattctggctgtatctGGGGGT
    # ACTACGGCCTTCTCTCCCCCtgtgctgtgctctaTTTTTtattctggctgtatctGGGGGTacgatgctatactgg
    a = Dseqrecord("ACTACGGCCTTCTCTCCCCCtgtgctgtgctcta", name="one34")
    # ||||||||||||||
    # tgtgctgtgctcta 14
    # ||||||||||||||
    b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatct", name="two35")
    # ||||||||||||||||
    # tattctggctgtatct 16
    # ||||||||||||||||
    c = Dseqrecord("tattctggctgtatctGGGGGTacgatgctatactgg", name="three37")
    # ACTACGGCCTTCTCTCCCCCtgtgctgtgctctaTTTTTtattctggctgtatctGGGGGTacgatgctatactgg

    a.features = [
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(34), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(1), ExactPosition(33), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(1), ExactPosition(20), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(20), ExactPosition(33), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(21), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(19), ExactPosition(34), strand=1), type="misc"
        ),
    ]

    b.features = [
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(35), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(1), ExactPosition(34), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(19), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(19), ExactPosition(35), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(20), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(18), ExactPosition(35), strand=1), type="misc"
        ),
    ]

    c.features = [
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(37), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(1), ExactPosition(36), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(16), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(16), ExactPosition(37), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(17), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(15), ExactPosition(37), strand=1), type="misc"
        ),
    ]

    ln0 = assembly.Assembly((a, b, c), limit=14)
    dseqr = ln0.assemble_linear()[0]

    assert (
        str(dseqr.seq)
        == "ACTACGGCCTTCTCTCCCCCtgtgctgtgctctaTTTTTtattctggctgtatctGGGGGTacgatgctatactgg"
    )

    feature_seqs = (
        [f.extract(a).seq for f in a.features]
        + [f.extract(b).seq for f in b.features]
        + [f.extract(c).seq for f in c.features]
    )

    assembled_feature_seqs = [f.extract(dseqr).seq for f in dseqr.features]
    for f1, f2 in zip(feature_seqs, assembled_feature_seqs):
        assert f1 == f2

    a = Dseqrecord("ACTACGGCCTTCTCTCCCCCtgtgctgtgctcta", name="one34")
    a.add_feature(1, 33, label="first")  # ||||||||||||||
    # h0FkKjiojpATbfkOM0C3u4zg9hs
    # tgtgctgtgctcta 14
    # --------------
    brc = Dseqrecord("agatacagccagaataAAAAAtagagcacagcaca", name="twoArc35")
    # tgtgctgtgctctaTTTTTtattctggctgtatct
    brc.add_feature(1, 34, label="scnd")  # --------------
    # gmfjuQLVSPP4ayjJMPuig1jxxmE
    # tattctggctgtatct 16
    c = Dseqrecord("tattctggctgtatctGGGGGTacgatgctatactgg", name="three37")
    c.add_feature(1, 36, label="third")
    ln1 = assembly.Assembly((a, brc, c), limit=14)

    assert (
        str(ln1.assemble_linear()[0].seq)
        == "ACTACGGCCTTCTCTCCCCCtgtgctgtgctctaTTTTTtattctggctgtatctGGGGGTacgatgctatactgg"
    )

    a = Dseqrecord("ACTACGGCCTTCTCTCCCCCtgtgctgtgctcta", name="one")
    a.add_feature(1, 33, label="first")
    # h0FkKjiojpATbfkOM0C3u4zg9hs
    # tgtgctgtgctcta 14
    b2 = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatct", name="twoA")
    b2.add_feature(1, 34, label="scnd")
    b3 = Dseqrecord("tgtgctgtgctctaCCtattctggctgtatct", name="twoB")
    # gmfjuQLVSPP4ayjJMPuig1jxxmE
    # tattctggctgtatct 16
    c = Dseqrecord("tattctggctgtatctGGGGGTacgatgctatactgg", name="three")
    c.add_feature(1, 36, label="third")
    ln2 = assembly.Assembly((a, b2, b3, c), limit=14)
    linprods = ln2.assemble_linear()
    assert (
        str(linprods[0].seq)
        == "ACTACGGCCTTCTCTCCCCCtgtgctgtgctctaTTTTTtattctggctgtatctGGGGGTacgatgctatactgg"
    )
    assert (
        str(linprods[1].seq)
        == "ACTACGGCCTTCTCTCCCCCtgtgctgtgctctaCCtattctggctgtatctGGGGGTacgatgctatactgg"
    )

    # acgatgctatactgg 15
    a = Dseqrecord("acgatgctatactggCCCCCtgtgctgtgctctaGG", name="one36")
    # tgtgctgtgctcta 14
    b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatct", name="two35")
    # tattctggctgtatct 16
    c = Dseqrecord("tattctggctgtatctGGGGGTacgatgctatactgg", name="three37")
    # acgatgctatactgg 15
    a.features = [
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(34), strand=1), type="misc1"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(1), ExactPosition(33), strand=1), type="misc2"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(20), strand=1), type="misc3"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(20), ExactPosition(34), strand=1),
            type="misc4",
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(21), strand=1), type="misc5"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(19), ExactPosition(34), strand=1),
            type="misc6",
        ),
    ]

    b.features = [
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(35), strand=1), type="misc7"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(1), ExactPosition(34), strand=1), type="misc8"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(19), strand=1), type="misc9"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(19), ExactPosition(35), strand=1),
            type="misc10",
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(20), strand=1),
            type="misc11",
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(18), ExactPosition(35), strand=1),
            type="misc12",
        ),
    ]

    c.features = [
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(37), strand=1),
            type="misc13",
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(1), ExactPosition(36), strand=1),
            type="misc14",
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(16), strand=1),
            type="misc15",
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(16), ExactPosition(37), strand=1),
            type="misc16",
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(17), strand=1),
            type="misc17",
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(15), ExactPosition(37), strand=1),
            type="misc18",
        ),
    ]
    c1 = assembly.Assembly((a, b, c), limit=14)
    result = c1.assemble_circular()[0]
    assert result.seguid() == "cdseguid=CRIbOfddcwCZbvVOOOU4uJYP-So"
    assert (
        str(result.seq)
        == "acgatgctatactggCCCCCtgtgctgtgctctaTTTTTtattctggctgtatctGGGGGT"
    )
    # acgatgctatactggCCCCCtgtgctgtgctctaGG
    feature_seqs = (
        [f.extract(a).seq for f in a.features]
        + [f.extract(b).seq for f in b.features]
        + [f.extract(c).seq for f in c.features]
    )

    assembled_feature_seqs = [f.extract(result).seq for f in result.features]

    for f1, f2 in zip(feature_seqs, assembled_feature_seqs):
        assert f1 == f2

    a = Dseqrecord("acgatgctatactggCCCCCtgtgctgtgctctaGG", name="oneC")
    a.add_feature(1, 33, label="first")
    # h0FkKjiojpATbfkOM0C3u4zg9hs
    # tgtgctgtgctcta 14
    b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatct", name="twoA")
    b.add_feature(1, 34, label="scnd")
    b2 = Dseqrecord("tgtgctgtgctctaCCtattctggctgtatct", name="twoB")
    # gmfjuQLVSPP4ayjJMPuig1jxxmE
    # tattctggctgtatct 16
    c = Dseqrecord("tattctggctgtatctGGGGGTacgatgctatactgg", name="three")
    c.add_feature(1, 36, label="third")
    # VJtsIfDO2DkKXbW-sLF3nJ-AEe4
    # acgatgctatactgg 15

    c2 = assembly.Assembly((a, b, b2, c), limit=14)
    circprods = c2.assemble_circular()

    # Changed
    assert circprods[0].seguid() == "cdseguid=CRIbOfddcwCZbvVOOOU4uJYP-So"
    # assert circprods[1].seguid() == "cdseguid=CRIbOfddcwCZbvVOOOU4uJYP-So"
    assert circprods[1].seguid() == "cdseguid=zFIq5LWXL_YSxrSF2Q5hbzO0BPw"
    # assert circprods[1].seguid() == "cdseguid=zFIq5LWXL_YSxrSF2Q5hbzO0BPw"
    assert (
        str(circprods[0].seq)
        == "acgatgctatactggCCCCCtgtgctgtgctctaTTTTTtattctggctgtatctGGGGGT"
    )
    assert (
        str(circprods[1].seq)
        == "acgatgctatactggCCCCCtgtgctgtgctctaCCtattctggctgtatctGGGGGT"
    )

    # VJtsIfDO2DkKXbW-sLF3nJ-AEe4
    # acgatgctatactgg 15
    a = Dseqrecord("acgatgctatactggCCCCCtgtgctgtgctcta", name="oneC")

    # h0FkKjiojpATbfkOM0C3u4zg9hs
    # tgtgctgtgctcta 14
    b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatct", name="twoA")
    b.add_feature(1, 34, label="scnd")
    # gmfjuQLVSPP4ayjJMPuig1jxxmE
    # tattctggctgtatct 16
    c = Dseqrecord("tattctggctgtatctGGGGGTacgatgctatactgg", name="three")
    # VJtsIfDO2DkKXbW-sLF3nJ-AEe4
    # acgatgctatactgg 15

    c3 = assembly.Assembly((a, b, c), limit=14)
    assert (
        str(c3.assemble_circular()[0].seq)
        == "acgatgctatactggCCCCCtgtgctgtgctctaTTTTTtattctggctgtatctGGGGGT"
    )

    text1 = """
    >A_AgTEFp_b_631 NP+geg/4Ykv2pIwEqiLylYKPYOE
    TCCTGACGGGTAATTTTGATTTGCATGCCGTCCGGGTGAGTCATAGCGTCTGGTTGTTTTGCCAGATTCAGCAGAGTCTGTGCAATGCGGCCGCTGACcacatacgatttaggtgacactatagaacGCGGCCGCCAGCTGAAGCTTCGTACGCTGCAGGTCGACGGATCCCCGGGTTAATTAAGGCGCGCCAGATCTGTTTAGCTTGCCTTGTCCCCGCCGGGTCACCCGGCCAGCGACATGGAGGCCCAGAATACCCTCCTTGACAGTCTTGACGTGCGCAGCTCAGGGGCATGATGTGACTGTCGCCCGTACATTTAGCCCATACATCCCCATGTATAATCATTTGCATCCATACATTTTGATGGCCGCACGGCGCGAAGCAAAAATTACGGCTCCTCGCTGCAGACCTGCGAGCAGGGAAACGCTCCCCTCACAGACGCGTTGAATTGTCCCCACGCCGCGCCCCTGTAGAGAAATATAAAAGGTTAGGATTTGCCACTGAGGTTCTTCTTTCATATACTTCCTTTTAAAATCTTGCTAGGATACAGTTCTCACATcacatccgaacataaacaaccGTCGAGGAACGCCAGGTTGCCCACTTTctcactagtgacctgcagccgAC

    >B_hph_c_740 k2IHgW5Q2jg2iEJq/l7jCnq2mKM
    gtcgaggaacgccaggttgcccactTTCTCACTAGTGACCTGCAGCCGACCCAAtcacatcacatccgaacataaacaaccatgGGTACCACTCTTGACGACACGGCTTACCGGTACCGCACCAGTGTCCCGGGGGACGCCGAGGCCATCGAGGCACTGGATGGGTCCTTCACCACCGACACCGTCTTCCGCGTCACCGCCACCGGGGACGGCTTCACCCTGCGGGAGGTGCCGGTGGACCCGCCCCTGACCAAGGTGTTCCCCGACGACGAATCGGACGACGAATCGGACGACGGGGAGGACGGCGACCCGGACTCCCGGACGTTCGTCGCGTACGGGGACGACGGCGACCTGGCGGGCTTCGTGGTCGTCTCGTACTCCGGCTGGAACCGCCGGCTGACCGTCGAGGACATCGAGGTCGCCCCGGAGCACCGGGGGCACGGGGTCGGGCGCGCGTTGATGGGGCTCGCGACGGAGTTCGCCCGCGAGCGGGGCGCCGGGCACCTCTGGCTGGAGGTCACCAACGTCAACGCACCGGCGATCCACGCGTACCGGCGGATGGGGTTCACCCTCTGCGGCCTGGACACCGCCCTGTACGACGGCACCGCCTCGGACGGCGAGCAGGCGCTCTACATGAGCATGCCCTGCCCCtaatcagtactgacaataaaaagattcttgTAGGGATAACAGGGTAATCGGAGTGCCAtctgtgcagacaaacgcatcaggatTTAAAT

    >C_KlLEU2tt_d_650 8lAwzHM60BkV1hhzx3/bBsfAIYo
    CGCGCCATCTGTGCAGACAAACGCATCAGGATagagtcttttgtaacgacccCGTCTCCACCAACTTGGTATGCTTGAAATCTCAAGGCCATTACACATTCAGTTATGTGAACGAAAGGTCTTTATTTAACGTAGCATAAACTAAATAATACAGGTTCCGGTTAGCCTGCAATGTGTTAAATCTAAAGGAGCATACCCAAAATGAACTGAAGACAAGGAAATTTGCTTGTCCAGATGTGATTGAGCATTTGAACGTTAATAACATAACATTTTTATACTTAACTATAGAAAGACTTGTATAAAAACTGGCAAACGAGATATTCTGAATATTGGTGCATATTTCAGGTAGAAAAGCTTACAAAACAATCTAATCATAATATTGAGATGAAGAGAAAGATAAAAGAAAAAACGATAAGTCAGATGAGATTATGATTGTACTTTGAAATCGAGGAACAAAGTATATACGGTAGTAGTTCCCCGAGTTATAACGGGAGATCATGTAAATTGAGAAACCAGATAAAGATTTGGTATGCACTCTAGCAAGAAAATAAAATGATGAATCTATGAtatagatcacttttgttccagcATCCGGATTTACCTGAATCAATTGGCGAAATTTTTTGTACGAAATTTCAGCCACTTCACAG
    """

    list_of_formatted_seq_records = parse(text1)
    a = assembly.Assembly(list_of_formatted_seq_records, limit=25)

    candidate = a.assemble_linear()[0]
    correct = "TCCTGACGGGTAATTTTGATTTGCATGCCGTCCGGGTGAGTCATAGCGTCTGGTTGTTTTGCCAGATTCAGCAGAGTCTGTGCAATGCGGCCGCTGACcacatacgatttaggtgacactatagaacGCGGCCGCCAGCTGAAGCTTCGTACGCTGCAGGTCGACGGATCCCCGGGTTAATTAAGGCGCGCCAGATCTGTTTAGCTTGCCTTGTCCCCGCCGGGTCACCCGGCCAGCGACATGGAGGCCCAGAATACCCTCCTTGACAGTCTTGACGTGCGCAGCTCAGGGGCATGATGTGACTGTCGCCCGTACATTTAGCCCATACATCCCCATGTATAATCATTTGCATCCATACATTTTGATGGCCGCACGGCGCGAAGCAAAAATTACGGCTCCTCGCTGCAGACCTGCGAGCAGGGAAACGCTCCCCTCACAGACGCGTTGAATTGTCCCCACGCCGCGCCCCTGTAGAGAAATATAAAAGGTTAGGATTTGCCACTGAGGTTCTTCTTTCATATACTTCCTTTTAAAATCTTGCTAGGATACAGTTCTCACATcacatccgaacataaacaaccGTCGAGGAACGCCAGGTTGCCCACTTTctcactagtgacctgcagccgACCCAAtcacatcacatccgaacataaacaaccatgGGTACCACTCTTGACGACACGGCTTACCGGTACCGCACCAGTGTCCCGGGGGACGCCGAGGCCATCGAGGCACTGGATGGGTCCTTCACCACCGACACCGTCTTCCGCGTCACCGCCACCGGGGACGGCTTCACCCTGCGGGAGGTGCCGGTGGACCCGCCCCTGACCAAGGTGTTCCCCGACGACGAATCGGACGACGAATCGGACGACGGGGAGGACGGCGACCCGGACTCCCGGACGTTCGTCGCGTACGGGGACGACGGCGACCTGGCGGGCTTCGTGGTCGTCTCGTACTCCGGCTGGAACCGCCGGCTGACCGTCGAGGACATCGAGGTCGCCCCGGAGCACCGGGGGCACGGGGTCGGGCGCGCGTTGATGGGGCTCGCGACGGAGTTCGCCCGCGAGCGGGGCGCCGGGCACCTCTGGCTGGAGGTCACCAACGTCAACGCACCGGCGATCCACGCGTACCGGCGGATGGGGTTCACCCTCTGCGGCCTGGACACCGCCCTGTACGACGGCACCGCCTCGGACGGCGAGCAGGCGCTCTACATGAGCATGCCCTGCCCCtaatcagtactgacaataaaaagattcttgTAGGGATAACAGGGTAATCGGAGTGCCATCTGTGCAGACAAACGCATCAGGATagagtcttttgtaacgacccCGTCTCCACCAACTTGGTATGCTTGAAATCTCAAGGCCATTACACATTCAGTTATGTGAACGAAAGGTCTTTATTTAACGTAGCATAAACTAAATAATACAGGTTCCGGTTAGCCTGCAATGTGTTAAATCTAAAGGAGCATACCCAAAATGAACTGAAGACAAGGAAATTTGCTTGTCCAGATGTGATTGAGCATTTGAACGTTAATAACATAACATTTTTATACTTAACTATAGAAAGACTTGTATAAAAACTGGCAAACGAGATATTCTGAATATTGGTGCATATTTCAGGTAGAAAAGCTTACAAAACAATCTAATCATAATATTGAGATGAAGAGAAAGATAAAAGAAAAAACGATAAGTCAGATGAGATTATGATTGTACTTTGAAATCGAGGAACAAAGTATATACGGTAGTAGTTCCCCGAGTTATAACGGGAGATCATGTAAATTGAGAAACCAGATAAAGATTTGGTATGCACTCTAGCAAGAAAATAAAATGATGAATCTATGAtatagatcacttttgttccagcATCCGGATTTACCTGAATCAATTGGCGAAATTTTTTGTACGAAATTTCAGCCACTTCACAG"
    assert len(correct) == 1933
    assert eq(correct, candidate, circular=False)

    # Assembly:
    # Sequences........................: [632] [740] [650]
    # Sequences with shared homologies.: [632] [740] [650]
    # Homology limit (bp)..............: 25
    # Number of overlaps...............: 3
    # Nodes in graph(incl. 5' & 3')....: 5
    # Only terminal overlaps...........: No
    # Circular products................: [81]
    # Linear products..................: [1933] [1852] [1352] [1321] [1240]
    # [821] [713] [132] [51] [38]


def test_assembly():

    text1 = """
    >A_AgTEFp_b_631 NP+geg/4Ykv2pIwEqiLylYKPYOE
    TCCTGACGGGTAATTTTGATTTGCATGCCGTCCGGGTGAGTCATAGCGTCTGGTTGTTTTGCCAGATTCAGCAGAGTCTGTGCAATGCGGCCGCTGACcacatacgatttaggtgacactatagaacGCGGCCGCCAGCTGAAGCTTCGTACGCTGCAGGTCGACGGATCCCCGGGTTAATTAAGGCGCGCCAGATCTGTTTAGCTTGCCTTGTCCCCGCCGGGTCACCCGGCCAGCGACATGGAGGCCCAGAATACCCTCCTTGACAGTCTTGACGTGCGCAGCTCAGGGGCATGATGTGACTGTCGCCCGTACATTTAGCCCATACATCCCCATGTATAATCATTTGCATCCATACATTTTGATGGCCGCACGGCGCGAAGCAAAAATTACGGCTCCTCGCTGCAGACCTGCGAGCAGGGAAACGCTCCCCTCACAGACGCGTTGAATTGTCCCCACGCCGCGCCCCTGTAGAGAAATATAAAAGGTTAGGATTTGCCACTGAGGTTCTTCTTTCATATACTTCCTTTTAAAATCTTGCTAGGATACAGTTCTCACATcacatccgaacataaacaaccGTCGAGGAACGCCAGGTTGCCCACTTTctcactagtgacctgcagccgAC

    >B_hph_c_740 k2IHgW5Q2jg2iEJq/l7jCnq2mKM
    gtcgaggaacgccaggttgcccactTTCTCACTAGTGACCTGCAGCCGACCCAAtcacatcacatccgaacataaacaaccatgGGTACCACTCTTGACGACACGGCTTACCGGTACCGCACCAGTGTCCCGGGGGACGCCGAGGCCATCGAGGCACTGGATGGGTCCTTCACCACCGACACCGTCTTCCGCGTCACCGCCACCGGGGACGGCTTCACCCTGCGGGAGGTGCCGGTGGACCCGCCCCTGACCAAGGTGTTCCCCGACGACGAATCGGACGACGAATCGGACGACGGGGAGGACGGCGACCCGGACTCCCGGACGTTCGTCGCGTACGGGGACGACGGCGACCTGGCGGGCTTCGTGGTCGTCTCGTACTCCGGCTGGAACCGCCGGCTGACCGTCGAGGACATCGAGGTCGCCCCGGAGCACCGGGGGCACGGGGTCGGGCGCGCGTTGATGGGGCTCGCGACGGAGTTCGCCCGCGAGCGGGGCGCCGGGCACCTCTGGCTGGAGGTCACCAACGTCAACGCACCGGCGATCCACGCGTACCGGCGGATGGGGTTCACCCTCTGCGGCCTGGACACCGCCCTGTACGACGGCACCGCCTCGGACGGCGAGCAGGCGCTCTACATGAGCATGCCCTGCCCCtaatcagtactgacaataaaaagattcttgTAGGGATAACAGGGTAATCGGAGTGCCAtctgtgcagacaaacgcatcaggatTTAAAT

    >C_KlLEU2tt_d_650 8lAwzHM60BkV1hhzx3/bBsfAIYo
    CGCGCCATCTGTGCAGACAAACGCATCAGGATagagtcttttgtaacgacccCGTCTCCACCAACTTGGTATGCTTGAAATCTCAAGGCCATTACACATTCAGTTATGTGAACGAAAGGTCTTTATTTAACGTAGCATAAACTAAATAATACAGGTTCCGGTTAGCCTGCAATGTGTTAAATCTAAAGGAGCATACCCAAAATGAACTGAAGACAAGGAAATTTGCTTGTCCAGATGTGATTGAGCATTTGAACGTTAATAACATAACATTTTTATACTTAACTATAGAAAGACTTGTATAAAAACTGGCAAACGAGATATTCTGAATATTGGTGCATATTTCAGGTAGAAAAGCTTACAAAACAATCTAATCATAATATTGAGATGAAGAGAAAGATAAAAGAAAAAACGATAAGTCAGATGAGATTATGATTGTACTTTGAAATCGAGGAACAAAGTATATACGGTAGTAGTTCCCCGAGTTATAACGGGAGATCATGTAAATTGAGAAACCAGATAAAGATTTGGTATGCACTCTAGCAAGAAAATAAAATGATGAATCTATGAtatagatcacttttgttccagcATCCGGATTTACCTGAATCAATTGGCGAAATTTTTTGTACGAAATTTCAGCCACTTCACAG
    """

    text2 = """
    >693
    ttctagaactagtggatcccccgggctgcagatgagtgaaggccccgtcaaattcgaaaaaaataccgtcatatctgtctttggtgcgtcaggtgatctggcaaagaagaagacttttcccgccttatttgggcttttcagagaaggttaccttgatccatctaccaagatcttcggttatgcccggtccaaattgtccatggaggaggacctgaagtcccgtgtcctaccccacttgaaaaaacctcacggtgaagccgatgactctaaggtcgaacagttcttcaagatggtcagctacatttcgggaaattacgacacagatgaaggcttcgacgaattaagaacgcagatcgagaaattcgagaaaagtgccaacgtcgatgtcccacaccgtctcttctatctggccttgccgccaagcgtttttttgacggtggccaagcagatcaagagtcgtgtgtacgcagagaatggcatcacccgtgtaatcgtagagaaacctttcggccacgacctggcctctgccagggagctgcaaaaaaacctggggcccctctttaaagaagaagagttgtacagaattgaccattacttgggtaaagagttggtcaagaatcttttagtcttgaggttcggtaaccagtttttgaatgcctcgtggaatagagacaacattcaaagcgttcagat

    >934
    tatcgataagcttgatatcgaattcctgcagctaattatccttcgtatcttctggcttagtcacgggccaagcgtaagggtgcttttcgggcataacatacttgtgtttttgcatatattccttcaatccctttggacctcttgatccgtaggggtaaatttccggtgttggaccgtccggacgctctatgtgcttcagtaatggggtgaatatgccccaactgatatccaattcgtcatctctgacaaagttggaatggtcacccagtagggcgtctcttatcaacacctcgtaagcctctggaatccaaaagtcttggtacctgcttgcgtaagttagattcagatctgtgacttgggtagcatttgacagaccaggggtcttagcattaaactttaggtacacagcggcatcgggctgcactctgatgaccagttcgttatttggaatgtctttgaagacacccgatgcgaccgctttgtactgcagtctgatctccaccttggactcattcaaagccttaccggcacgcatcatgatggggacgccctcccaacgctcgttttcgatgttgaaagtcattgctgcaaaagtgacacatttagagtccttgtctacagtgtcatcatccacgtaggcgggcttagacccgtcctcagatttaccgtactggcccaagaggacgtcgtccgtgtcgatgggggccacggcctttagaaccttaaccttttcgtcacgaatagattccgggtcaaaagacaccggtctttccatagtcaagagagtcatgatttgtaacagatggttctgcatcacgtctctgattatgcctatagagtcgaaatagccgccacggccttcggtgccgaacctctctttaaacgaaatctgaacgctttgaatgttgtctctattccacgaggcattcaaaaa

    >7729
    gaattcgatatcaagcttatcgataccgtcgacctcgagtcatgtaattagttatgtcacgcttacattcacgccctccccccacatccgctctaaccgaaaaggaaggagttagacaacctgaagtctaggtccctatttatttttttatagttatgttagtattaagaacgttatttatatttcaaatttttcttttttttctgtacagacgcgtgtacgcatgtaacattatactgaaaaccttgcttgagaaggttttgggacgctcgaaggctttaatttgcggccggtacccaattcgccctatagtgagtcgtattacgcgcgctcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaatggcgaatggcgcgacgcgccctgtagcggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgctcctttcgctttcttcccttcctttctcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccctttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttgattagggtgatggttcacgtagtgggccatcgccctgatagacggtttttcgccctttgacgttggagtccacgttctttaatagtggactcttgttccaaactggaacaacactcaaccctatctcggtctattcttttgatttataagggattttgccgatttcggcctattggttaaaaaatgagctgatttaacaaaaatttaacgcgaattttaacaaaatattaacgtttacaatttcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatcgacggtcgaggagaacttctagtatatccacatacctaatattattgccttattaaaaatggaatcccaacaattacatcaaaatccacattctcttcaaaatcaattgtcctgtacttccttgttcatgtgtgttcaaaaacgttatatttataggataattatactctatttctcaacaagtaattggttgtttggccgagcggtctaaggcgcctgattcaagaaatatcttgaccgcagttaactgtgggaatactcaggtatcgtaagatgcaagagttcgaatctcttagcaaccattatttttttcctcaacataacgagaacacacaggggcgctatcgcacagaatcaaattcgatgactggaaattttttgttaatttcagaggtcgcctgacgcatatacctttttcaactgaaaaattgggagaaaaaggaaaggtgagaggccggaaccggcttttcatatagaatagagaagcgttcatgactaaatgcttgcatcacaatacttgaagttgacaatattatttaaggacctattgttttttccaataggtggttagcaatcgtcttactttctaacttttcttaccttttacatttcagcaatatatatatatatttcaaggatataccattctaatgtctgcccctatgtctgcccctaagaagatcgtcgttttgccaggtgaccacgttggtcaagaaatcacagccgaagccattaaggttcttaaagctatttctgatgttcgttccaatgtcaagttcgatttcgaaaatcatttaattggtggtgctgctatcgatgctacaggtgtcccacttccagatgaggcgctggaagcctccaagaaggttgatgccgttttgttaggtgctgtggctggtcctaaatggggtaccggtagtgttagacctgaacaaggtttactaaaaatccgtaaagaacttcaattgtacgccaacttaagaccatgtaactttgcatccgactctcttttagacttatctccaatcaagccacaatttgctaaaggtactgacttcgttgttgtcagagaattagtgggaggtatttactttggtaagagaaaggaagacgatggtgatggtgtcgcttgggatagtgaacaatacaccgttccagaagtgcaaagaatcacaagaatggccgctttcatggccctacaacatgagccaccattgcctatttggtccttggataaagctaatcttttggcctcttcaagattatggagaaaaactgtggaggaaaccatcaagaacgaattccctacattgaaggttcaacatcaattgattgattctgccgccatgatcctagttaagaacccaacccacctaaatggtattataatcaccagcaacatgtttggtgatatcatctccgatgaagcctccgttatcccaggttccttgggtttgttgccatctgcgtccttggcctctttgccagacaagaacaccgcatttggtttgtacgaaccatgccacggttctgctccagatttgccaaagaataaggttgaccctatcgccactatcttgtctgctgcaatgatgttgaaattgtcattgaacttgcctgaagaaggtaaggccattgaagatgcagttaaaaaggttttggatgcaggtatcagaactggtgatttaggtggttccaacagtaccaccgaagtcggtgatgctgtcgccgaagaagttaagaaaatccttgcttaaaaagattctctttttttatgatatttgtacataaactttataaatgaaattcataatagaaacgacacgaaattacaaaatggaatatgttcatagggtagacgaaactatatacgcaatctacatacatttatcaagaaggagaaaaaggaggatagtaaaggaatacaggtaagcaaattgatactaatggctcaacgtgataaggaaaaagaattgcactttaacattaatattgacaaggaggagggcaccacacaaaaagttaggtgtaacagaaaatcatgaaactacgattcctaatttgatattggaggattttctctaaaaaaaaaaaaatacaacaaataaaaaacactcaatgacctgaccatttgatggagtttaagtcaataccttcttgaagcatttcccataatggtgaaagttccctcaagaattttactctgtcagaaacggccttacgacgtagtcgatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagtatgatccaatatcaaaggaaatgatagcattgaaggatgagactaatccaattgaggagtggcagcatatagaacagctaaagggtagtgctgaaggaagcatacgataccccgcatggaatgggataatatcacaggaggtactagactacctttcatcctacataaatagacgcatataagtacgcatttaagcataaacacgcactatgccgttcttctcatgtatatatatatacaggcaacacgcagatataggtgcgacgtgaacagtgagctgtatgtgcgcagctcgcgttgcattttcggaagcgctcgttttcggaaacgctttgaagttcctattccgaagttcctattctctagaaagtataggaacttcagagcgcttttgaaaaccaaaagcgctctgaagacgcactttcaaaaaaccaaaaacgcaccggactgtaacgagctactaaaatattgcgaataccgcttccacaaacattgctcaaaagtatctctttgctatatatctctgtgctatatccctatataacctacccatccacctttcgctccttgaacttgcatctaaactcgacctctacattttttatgtttatctctagtattactctttagacaaaaaaattgtagtaagaactattcatagagtgaatcgaaaacaatacgaaaatgtaaacatttcctatacgtagtatatagagacaaaatagaagaaaccgttcataattttctgaccaatgaagaatcatcaacgctatcactttctgttcacaaagtatgcgcaatccacatcggtatagaatataatcggggatgcctttatcttgaaaaaatgcacccgcagcttcgctagtaatcagtaaacgcgggaagtggagtcaggctttttttatggaagagaaaatagacaccaaagtagccttcttctaaccttaacggacctacagtgcaaaaagttatcaagagactgcattatagagcgcacaaaggagaaaaaaagtaatctaagatgctttgttagaaaaatagcgctctcgggatgcatttttgtagaacaaaaaagaagtatagattctttgttggtaaaatagcgctctcgcgttgcatttctgttctgtaaaaatgcagctcagattctttgtttgaaaaattagcgctctcgcgttgcatttttgttttacaaaaatgaagcacagattcttcgttggtaaaatagcgctttcgcgttgcatttctgttctgtaaaaatgcagctcagattctttgtttgaaaaattagcgctctcgcgttgcatttttgttctacaaaatgaagcacagatgcttcgttcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgtagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaagcgggcagtgagcgcaacgcaattaatgtgagttacctcactcattaggcaccccaggctttacactttatgcttccggctcctatgttgtgtggaattgtgagcggataacaatttcacacaggaaacagctatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctcagtttatcattatcaatactcgccatttcaaagaatacgtaaataattaatagtagtgattttcctaactttatttagtcaaaaaattagccttttaattctgctgtaacccgtacatgcccaaaatagggggcgggttacacagaatatataacatcgtaggtgtctgggtgaacagtttattcctggcatccactaaatataatggagcccgctttttaagctggcatccagaaaaaaaaagaatcccagcaccaaaatattgttttcttcaccaaccatcagttcataggtccattctcttagcgcaactacagagaacaggggcacaaacaggcaaaaaacgggcacaacctcaatggagtgatgcaacctgcctggagtaaatgatgacacaaggcaattgacccacgcatgtatctatctcattttcttacaccttctattaccttctgctctctctgatttggaaaaagctgaaaaaaaaggttgaaaccagttccctgaaattattcccctacttgactaataagtatataaagacggtaggtattgattgtaattctgtaaatctatttcttaaacttcttaaattctacttttatagttagtcttttttttagttttaaaacaccagaacttagtttcgacggattctagaactagtggatcccccgggctgcag

    """

    text3 = """
    >1671bp ZdMoT7KL0/EgUFC01DddjNNfI/E (rc) WcvvNc76t2UYI2Py4Novcdu7bCY
    TTCTATAGACACACAAACACAAATACACACACTAAATTAATAATGGATGTCCACGAGGTCTCTATATCGGGATCAGCCTGCCTCGTACGCTGCAGGTCGACGGATCCCCGGGTTAATTAAGGCGCGCCAGATCTGTTTAGCTTGCCTCGTCCCCGCCGGGTCACCCGGCCAGCGACATGGAGGCCCAGAATACCCTCCTTGACAGTCTTGACGTGCGCAGCTCAGGGGCATGATGTGACTGTCGCCCGTACATTTAGCCCATACATCCCCATGTATAATCATTTGCATCCATACATTTTGATGGCCGCACGGCGCGAAGCAAAAATTACGGCTCCTCGCTGCAGACCTGCGAGCAGGGAAACGCTCCCCTCACAGACGCGTTGAATTGTCCCCACGCCGCGCCCCTGTAGAGAAATATAAAAGGTTAGGATTTGCCACTGAGGTTCTTCTTTCATATACTTCCTTTTAAAATCTTGCTAGGATACAGTTCTCACATCACATCCGAACATAAACAACCATGGGTAAGGAAAAGACTCACGTTTCGAGGCCGCGATTAAATTCCAACATGGATGCTGATTTATATGGGTATAAATGGGCTCGCGATAATGTCGGGCAATCAGGTGCGACAATCTATCGATTGTATGGGAAGCCCGATGCGCCAGAGTTGTTTCTGAAACATGGCAAAGGTAGCGTTGCCAATGATGTTACAGATGAGATGGTCAGACTAAACTGGCTGACGGAATTTATGCCTCTTCCGACCATCAAGCATTTTATCCGTACTCCTGATGATGCATGGTTACTCACCACTGCGATCCCCGGCAAAACAGCATTCCAGGTATTAGAAGAATATCCTGATTCAGGTGAAAATATTGTTGATGCGCTGGCAGTGTTCCTGCGCCGGTTGCATTCGATTCCTGTTTGTAATTGTCCTTTTAACAGCGATCGCGTATTTCGTCTCGCTCAGGCGCAATCACGAATGAATAACGGTTTGGTTGATGCGAGTGATTTTGATGACGAGCGTAATGGCTGGCCTGTTGAACAAGTCTGGAAAGAAATGCATAAGCTTTTGCCATTCTCACCGGATTCAGTCGTCACTCATGGTGATTTCTCACTTGATAACCTTATTTTTGACGAGGGGAAATTAATAGGTTGTATTGATGTTGGACGAGTCGGAATCGCAGACCGATACCAGGATCTTGCCATCCTATGGAACTGCCTCGGTGAGTTTTCTCCTTCATTACAGAAACGGCTTTTTCAAAAATATGGTATTGATAATCCTGATATGAATAAATTGCAGTTTCATTTGATGCTCGATGAGTTTTTCTAATCAGTACTGACAATAAAAAGATTCTTGTTTTCAAGAACTTGTCATTTGTATAGTTTTTTTATATTGTAGTTGTTCTATTTTAATCAAATGTTAGCGTGATTTATATTTTTTTTCGCCTCGACATCATCTGCCCAGATGCGAAGTTAAGTGCGCAGAAAGTAATATCATGCGTCAATCGTATGTGAATGCTGGTCGCTATACTGCTGTCGATTCGATACTAACGCCGCCATCCAGTGTCGAAAACGAGCTCGAATTCATCGATGAGCATCCTTATAGCCTCTTCTACGAGACCGACACCGTAAACAGGCCCCTTTTCCTTTGTCGATATCATGTAATTAGTTATG

    >YJR048W
    GAGGCACCAGCGTCAGCATTTTCAAAGGTGTGTTCTTCGTCAGACATGTTTTAGTGTGTGAATGAAATAGGTGTATGTTTTCTTTTTGCTAGACAATAATTAGGAACAAGGTAAGGGAACTAAAGTGTAGAATAAGATTAAAAAAGAAGAACAAGTTGAAAAGGCAAGTTGAAATTTCAAGAAAAAAGTCAATTGAAGTACAGTAAATTGACCTGAATATATCTGAGTTCCGACAACAATGAGTTTACCAAAGAGAACAATGGAATAGGAAACTTTGAACGAAGAAAGGAAAGCAGGAAAGGAAAAAATTTTTAGGCTCGAGAACAATAGGGCGAAAAAACAGGCAACGAACGAACAATGGAAAAACGAAAAAAAAAAAAAAAAACACAGAAAAGAATGCAGAAAGATGTCAACTGAAAAAAAAAAAGGTGAACACAGGAAAAAAAATAAAAAAAAAAAAAAAAAAAGGAGGACGAAACAAAAAAGTGAAAAAAAATGAAAATTTTTTTGGAAAACCAAGAAATGAATTATATTTCCGTGTGAGACGACATCGTCGAATATGATTCAGGGTAACAGTATTGATGTAATCAATTTCCTACCTGAATCTAAAATTCCCGGGAGCAAGATCAAGATGTTTTCACCGATCTTTCCGGTCTCTTTGGCCGGGGTTTACGGACGATGGCAGAAGACCAAAGCGCCAGTTCATTTGGCGAGCGTTGGTTGGTGGATCAAGCCCACGCGTAGGCAATCCTCGAGCAGATCCGCCAGGCGTGTATATATAGCGTGGATGGCCAGGCAACTTTAGTGCTGACACATACAGGCATATATATATGTGTGCGACGACACATGATCATATGGCATGCATGTGCTCTGTATGTATATAAAACTCTTGTTTTCTTCTTTTCTCTAAATATTCTTTCCTTATACATTAGGACCTTTGCAGCATAAATTACTATACTTCTATAGACACACAAACACAAATACACACACTAAATTAATAATGACTGAATTCAAGGCCGGTTCTGCTAAGAAAGGTGCTACACTTTTCAAGACTAGATGTCTACAATGCCACACCGTGGAAAAGGGTGGCCCACATAAGGTTGGTCCAAACTTGCATGGTATCTTTGGCAGACACTCTGGTCAAGCTGAAGGGTATTCGTACACAGATGCCAATATCAAGAAAAACGTGTTGTGGGACGAAAATAACATGTCAGAGTACTTGACTAACCCAAAGAAATATATTCCTGGTACCAAGATGGCCTTTGGTGGGTTGAAGAAGGAAAAAGACAGAAACGACTTAATTACCTACTTGAAAAAAGCCTGTGAGTAAACAGGCCCCTTTTCCTTTGTCGATATCATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCCTCCCACATCCGCTCTAACCGAAAAGGAAGGAGTTAGACAACCTGAAGTCTAGGTCCCTATTTATTTTTTTTAATAGTTATGTTAGTATTAAGAACGTTATTTATATTTCAAATTTTTCTTTTTTTTCTGTACAAACGCGTGTACGCATGTAACATTATACTGAAAACCTTGCTTGAGAAGGTTTTGGGACGCTCGAAGGCTTTAATTTGCAAGCTTCGCAGTTTACACTCTCATCGTCGCTCTCATCATCGCTTCCGTTGTTGTTTTCCTTAGTAGCGTCTGCTTCCAGAGAGTATTTATCTCTTATTACCTCTAAAGGTTCTGCTTGATTTCTGACTTTGTTCGCCTCATGTGCATATTTTTCTTGGTTCTTTTGGGACAAAATATGCGTAAAGGACTTTTGTTGTTCCCTCACATTCCAGTTTAGTTGTCGACTGATACTGTTAATAAACTCATCGGGCGAGGCTTCCACGGTTGGAAAAGCATATGGGCTGGCGCATATGGTTATAAAATCACCTTTTTGCAATTCAATTCTATCTTTCCCATCAAAAGCCGCCCATGCTGGAGCCCTTGACTTCATCGAGACTTTCACTTTTAAATTTATACTTTCTGGTAAGATGATGGGTCTGAAACTCAATGCATGTGGACAAATGGGTGTTAAAGCGATTGCATTGACGGTTGGGCATACCAATGACCCACCTGCACTCAAAGAATAGGCCGTGGACCCAGTCGGAGTAGCAGCAATCAGTCCGTCCGCCTGCGCAACGGTCATTAATGAGCCGTCACCATACAATTCTAACATGGATAGAAAAGGACTTGGACCACGATCGATGGTCACTTCGTTCAAAATGTGGTGTGTGCTTAGTTTTTCCACCACACATATTTTCTTCCCCGTGTTTGGGTCTACTTCAGGGCGGTGTCTACGATAAATTGTG

    """

    text4 = """
    >2389
    ctgtgaagtggctgaaatttcgtacaaaaaatttcgccaattgattcaggtaaatccggatgctggaacaaaagtgatctatatcatagattcatcattttattttcttgctagagtgcataccaaatctttatctggtttctcaatttacatgatctcccgttataactcggggaactactaccgtatatactttgttcctcgatttcaaagtacaatcataatctcatctgacttatcgttttttcttttatctttctcttcatctcaatattatgattagattgttttgtaagcttttctacctgaaatatgcaccaatattcagaatatctcgtttgccagtttttatacaagtctttctatagttaagtataaaaatgttatgttattaacgttcaaatgctcaatcacatctggacaagcaaatttccttgtcttcagttcattttgggtatgctcctttagatttaacacattgcaggctaaccggaacctgtattatttagtttatgctacgttaaataaagacctttcgttcacataactgaatgtgtaatggccttgagatttcaagcataccaagttggtggagacggggtcgttacaaaagactctatcctgatgcgtttgtctgcacagatggcactccgattaccctgttatccctacaagaatctttttattgtcagtactgattattcctttgccctcggacgagtgctggggcgtcggtttccactatcggcgagtacttctacacagccatcggtccagacggccgcgcttctgcgggcgatttgtgtacgcccgacagtcccggctccggatcggacgattgcgtcgcatcgaccctgcgcccaagctgcatcatcgaaattgccgtcaaccaagctctgatagagttggtcaagaccaatgcggagcatatacgcccggagccgcggcgatcctgcaagctccggatgcctccgctcgaagtagcgcgtctgctgctccatacaagccaaccacggcctccagaagaagatgttggcgacctcgtattgggaatccccgaacatcgcctcgctccagtcaatgaccgctgttatgcggccattgtccgtcaggacattgttggagccgaaatccgcgtgcacgaggtgccggacttcggggcagtcctcggcccaaagcatcagctcatcgagagcctgcgcgacggacgcactgacggtgtcgtccatcacagtttgccagtgatacacatggggatcagcaatcgcgcatatgaaatcacgccatgtagtgtattgaccgattccttgcggtccgaatgggccgaacccgctcgtctggctaagatcggccgcagcgatcgcatccatggcctccgcgaccggctgcagaacagcgggcagttcggtttcaggcaggtcttgcaacgtgacaccctgtgcacggcgggagatgcaataggtcaggctctcgctgaattccccaatgtcaagcacttccggaatcgggagcgcggccgatgcaaagtgccgataaacataacgatctttgtagaaaccatcggcgcagctatttacccgcaggacatatccacgccctcctacatcgaagctgaaagcacgagattcttcgccctccgagagctgcatcaggtcggagacgctgtcgaacttttcgatcagaaacttctcgacagacgtcgcggtgagttcaggctttttacccatggttgtttatgttcggatgtgatgtgattgggtcggctgcaggtcactagtgagaaagtgggcaacctggcgttcctcgacggttgtttatgttcggatgtgatgtgagaactgtatcctagcaagattttaaaaggaagtatatgaaagaagaacctcagtggcaaatcctaaccttttatatttctctacaggggcgcggcgtggggacaattcaacgcgtctgtgaggggagcgtttccctgctcgcaggtctgcagcgaggagccgtaatttttgcttcgcgccgtgcggccatcaaaatgtatggatgcaaatgattatacatggggatgtatgggctaaatgtacgggcgacagtcacatcatgcccctgagctgcgcacgtcaagactgtcaaggagggtattctgggcctccatgtcgctggccgggtgacccggcggggacaaggcaagctaaacagatctggcgcgccttaattaacccggggatccgtcgacctgcagcgtacgaagcttcagctggcggccgcgttctatagtgtcacctaaatcgtatgtggtcagcggccgcattgcacagactctgctgaatctggcaaaacaaccagacgctatgactcacccggacggcatgcaaatcaaaattacccgtcagga

    >2577
    tcctgacgggtaattttgatttgcatgccgtccgggtgagtcatagcgtctggttgttttgccagattcagcagagtctgtgcaatgcggccgctgacagagtcttttgtaacgaccccgtctccaccaacttggtatgcttgaaatctcaaggccattacacattcagttatgtgaacgaaaggtctttatttaacgtagcataaactaaataatacaggttccggttagcctgcaatgtgttaaatctaaaggagcatacccaaaatgaactgaagacaaggaaatttgcttgtccagatgtgattgagcatttgaacgttaataacataacatttttatacttaactatagaaagacttgtataaaaactggcaaacgagatattctgaatattggtgcatatttcaggtagaaaagcttacaaaacaatctaatcataatattgagatgaagagaaagataaaagaaaaaacgataagtcagatgagattatgattgtactttgaaatcgaggaacaaagtatatacggtagtagttccccgagttataacgggagatcatgtaaattgagaaaccagataaagatttggtatgcactctagcaagaaaataaaatgatgaatctatgatatagatcacttttgttccagcgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgacgatcagatctttcaggaaagtttcggaggagatagtgttcggcagtttgtacatcatctgcgggatcaggtacggtttgatcaggttgtagaagatcaggtaagacatagaatcgatgtagatgatcggtttgtttttgttgatttttacgtaacagttcagttggaatttgttacgcagacccttaaccaggtattctacttcttcgaaagtgaaagactgggtgttcagtacgatcgatttgttggtagagtttttgttgtaatcccatttaccaccatcatccatgaaccagtatgccagagacatcggggtcaggtagttttcaaccaggttgttcgggatggtttttttgttgttaacgatgaacaggttagccagtttgttgaaagcttggtgtttgaaagtctgggcgccccaggtgattaccaggttacccaggtggttaacacgttcttttttgtgcggcggggacagtacccactgatcgtacagcagacatacgtggtccatgtatgctttgtttttccactcgaactgcatacagtaggttttaccttcatcacgagaacggatgtaagcatcacccaggatcagaccgatacctgcttcgaactgttcgatgttcagttcgatcagctgggatttgtattctttcagcagtttagagttcggacccaggttcattacctggttttttttgatgtttttcatatgcatggatccggggttttttctccttgacgttaaagtatagaggtatattaacaattttttgttgatacttttattacatttgaataagaagtaatacaaaccgaaaatgttgaaagtattagttaaagtggttatgcagtttttgcatttatatatctgttaatagatcaaaaatcatcgcttcgctgattaattaccccagaaataaggctaaaaaactaatcgcattatcatcctatggttgttaatttgattcgttcatttgaaggtttgtggggccaggttactgccaatttttcctcttcataaccataaaagctagtattgtagaatctttattgttcggagcagtgcggcgcgaggcacatctgcgtttcaggaacgcgaccggtgaagacgaggacgcacggaggagagtcttccttcggagggctgtcacccgctcggcggcttctaatccgtacttcaatatagcaatgagcagttaagcgtattactgaaagttccaaagagaaggtttttttaggctaagataatggggctctttacatttccacaacatataagtaagattagatatggatatgtatatggatatgtatatggtggtaatgccatgtaatatgattattaaacttctttgcgtccatccaacgagatctggcgcgccttaattaacccaacctgcattaatgaatcggccaacgcgcggattaccctgttatccctacatattgttgtgccatctgtgcagacaaacgcatcaggattcagtactgacaataaaaagattcttgttttcaagaacttgtcatttgtatagtttttttatattgtagttgttctattttaatcaaatgttagcgtgatttatattttttttcgcctcgacatcatctgcccagatgcgaagttaagtgcgcagaaagtaatatcatgcgtcaatcgtatgtgaatgctggtcgctatactgctgtcgattcgatactaacgccgccatccagtgtcgaaaacgagctcgaattcatcgatgatatcagatccactagtggcctatgcggccgcggatctgccggtctccctatagtgagtcgatccggatttacctgaatcaattggcgaaattttttgtacgaaatttcagccacttcacag
    """

    text5 = """

    >2577
    tcctgacgggtaattttgatttgcatgccgtccgggtgagtcatagcgtctggttgttttgccagattcagcagagtctgtgcaatgcggccgctgacagagtcttttgtaacgaccccgtctccaccaacttggtatgcttgaaatctcaaggccattacacattcagttatgtgaacgaaaggtctttatttaacgtagcataaactaaataatacaggttccggttagcctgcaatgtgttaaatctaaaggagcatacccaaaatgaactgaagacaaggaaatttgcttgtccagatgtgattgagcatttgaacgttaataacataacatttttatacttaactatagaaagacttgtataaaaactggcaaacgagatattctgaatattggtgcatatttcaggtagaaaagcttacaaaacaatctaatcataatattgagatgaagagaaagataaaagaaaaaacgataagtcagatgagattatgattgtactttgaaatcgaggaacaaagtatatacggtagtagttccccgagttataacgggagatcatgtaaattgagaaaccagataaagatttggtatgcactctagcaagaaaataaaatgatgaatctatgatatagatcacttttgttccagcgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgacgatcagatctttcaggaaagtttcggaggagatagtgttcggcagtttgtacatcatctgcgggatcaggtacggtttgatcaggttgtagaagatcaggtaagacatagaatcgatgtagatgatcggtttgtttttgttgatttttacgtaacagttcagttggaatttgttacgcagacccttaaccaggtattctacttcttcgaaagtgaaagactgggtgttcagtacgatcgatttgttggtagagtttttgttgtaatcccatttaccaccatcatccatgaaccagtatgccagagacatcggggtcaggtagttttcaaccaggttgttcgggatggtttttttgttgttaacgatgaacaggttagccagtttgttgaaagcttggtgtttgaaagtctgggcgccccaggtgattaccaggttacccaggtggttaacacgttcttttttgtgcggcggggacagtacccactgatcgtacagcagacatacgtggtccatgtatgctttgtttttccactcgaactgcatacagtaggttttaccttcatcacgagaacggatgtaagcatcacccaggatcagaccgatacctgcttcgaactgttcgatgttcagttcgatcagctgggatttgtattctttcagcagtttagagttcggacccaggttcattacctggttttttttgatgtttttcatatgcatggatccggggttttttctccttgacgttaaagtatagaggtatattaacaattttttgttgatacttttattacatttgaataagaagtaatacaaaccgaaaatgttgaaagtattagttaaagtggttatgcagtttttgcatttatatatctgttaatagatcaaaaatcatcgcttcgctgattaattaccccagaaataaggctaaaaaactaatcgcattatcatcctatggttgttaatttgattcgttcatttgaaggtttgtggggccaggttactgccaatttttcctcttcataaccataaaagctagtattgtagaatctttattgttcggagcagtgcggcgcgaggcacatctgcgtttcaggaacgcgaccggtgaagacgaggacgcacggaggagagtcttccttcggagggctgtcacccgctcggcggcttctaatccgtacttcaatatagcaatgagcagttaagcgtattactgaaagttccaaagagaaggtttttttaggctaagataatggggctctttacatttccacaacatataagtaagattagatatggatatgtatatggatatgtatatggtggtaatgccatgtaatatgattattaaacttctttgcgtccatccaacgagatctggcgcgccttaattaacccaacctgcattaatgaatcggccaacgcgcggattaccctgttatccctacatattgttgtgccatctgtgcagacaaacgcatcaggattcagtactgacaataaaaagattcttgttttcaagaacttgtcatttgtatagtttttttatattgtagttgttctattttaatcaaatgttagcgtgatttatattttttttcgcctcgacatcatctgcccagatgcgaagttaagtgcgcagaaagtaatatcatgcgtcaatcgtatgtgaatgctggtcgctatactgctgtcgattcgatactaacgccgccatccagtgtcgaaaacgagctcgaattcatcgatgatatcagatccactagtggcctatgcggccgcggatctgccggtctccctatagtgagtcgatccggatttacctgaatcaattggcgaaattttttgtacgaaatttcagccacttcacag

    >2389
    ctgtgaagtggctgaaatttcgtacaaaaaatttcgccaattgattcaggtaaatccggatgctggaacaaaagtgatctatatcatagattcatcattttattttcttgctagagtgcataccaaatctttatctggtttctcaatttacatgatctcccgttataactcggggaactactaccgtatatactttgttcctcgatttcaaagtacaatcataatctcatctgacttatcgttttttcttttatctttctcttcatctcaatattatgattagattgttttgtaagcttttctacctgaaatatgcaccaatattcagaatatctcgtttgccagtttttatacaagtctttctatagttaagtataaaaatgttatgttattaacgttcaaatgctcaatcacatctggacaagcaaatttccttgtcttcagttcattttgggtatgctcctttagatttaacacattgcaggctaaccggaacctgtattatttagtttatgctacgttaaataaagacctttcgttcacataactgaatgtgtaatggccttgagatttcaagcataccaagttggtggagacggggtcgttacaaaagactctatcctgatgcgtttgtctgcacagatggcactccgattaccctgttatccctacaagaatctttttattgtcagtactgattattcctttgccctcggacgagtgctggggcgtcggtttccactatcggcgagtacttctacacagccatcggtccagacggccgcgcttctgcgggcgatttgtgtacgcccgacagtcccggctccggatcggacgattgcgtcgcatcgaccctgcgcccaagctgcatcatcgaaattgccgtcaaccaagctctgatagagttggtcaagaccaatgcggagcatatacgcccggagccgcggcgatcctgcaagctccggatgcctccgctcgaagtagcgcgtctgctgctccatacaagccaaccacggcctccagaagaagatgttggcgacctcgtattgggaatccccgaacatcgcctcgctccagtcaatgaccgctgttatgcggccattgtccgtcaggacattgttggagccgaaatccgcgtgcacgaggtgccggacttcggggcagtcctcggcccaaagcatcagctcatcgagagcctgcgcgacggacgcactgacggtgtcgtccatcacagtttgccagtgatacacatggggatcagcaatcgcgcatatgaaatcacgccatgtagtgtattgaccgattccttgcggtccgaatgggccgaacccgctcgtctggctaagatcggccgcagcgatcgcatccatggcctccgcgaccggctgcagaacagcgggcagttcggtttcaggcaggtcttgcaacgtgacaccctgtgcacggcgggagatgcaataggtcaggctctcgctgaattccccaatgtcaagcacttccggaatcgggagcgcggccgatgcaaagtgccgataaacataacgatctttgtagaaaccatcggcgcagctatttacccgcaggacatatccacgccctcctacatcgaagctgaaagcacgagattcttcgccctccgagagctgcatcaggtcggagacgctgtcgaacttttcgatcagaaacttctcgacagacgtcgcggtgagttcaggctttttacccatggttgtttatgttcggatgtgatgtgattgggtcggctgcaggtcactagtgagaaagtgggcaacctggcgttcctcgacggttgtttatgttcggatgtgatgtgagaactgtatcctagcaagattttaaaaggaagtatatgaaagaagaacctcagtggcaaatcctaaccttttatatttctctacaggggcgcggcgtggggacaattcaacgcgtctgtgaggggagcgtttccctgctcgcaggtctgcagcgaggagccgtaatttttgcttcgcgccgtgcggccatcaaaatgtatggatgcaaatgattatacatggggatgtatgggctaaatgtacgggcgacagtcacatcatgcccctgagctgcgcacgtcaagactgtcaaggagggtattctgggcctccatgtcgctggccgggtgacccggcggggacaaggcaagctaaacagatctggcgcgccttaattaacccggggatccgtcgacctgcagcgtacgaagcttcagctggcggccgcgttctatagtgtcacctaaatcgtatgtggtcagcggccgcattgcacagactctgctgaatctggcaaaacaaccagacgctatgactcacccggacggcatgcaaatcaaaattacccgtcagga

    >5681
    atccggatttacctgaatcaattggcgaaattttttgtacgaaatttcagccacttcacaggcggttttcgcacgtacccatgcgctacgttcctggccctcttcaaacaggcccagttcgccaataaaatcaccctgattcagataggagaggatcatttctttaccctcttcgtctttgatcagcactgccacagagcctttaacgatgtagtacagcgtttccgctttttcaccctggtgaataagcgtgctcttggatgggtacttatgaatgtggcaatgagacaagaaccattcgagagtaggatccgtttgaggtttaccaagtaccataagatccttaaatttttattatctagctagatgataatattatatcaagaattgtacctgaaagcaaataaattttttatctggcttaactatgcggcatcagagcagattgtactgagagtgcaccatatgcggtgtgaaataccgcacagatgcgtaaggagaaaataccgcatcaggcgctcttccgcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatacggttatccacagaatcaggggataacgcaggaaagaacatgtgagcaaaaggccagcaaaagcccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaaggacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgaggggtaataactgatataattaaattgaagctctaatttgtgagtttagtatacatgcatttacttataatacagttttttagttttgctggccgcatcttctcaaatatgcttcccagcctgcttttctgtaacgttcaccctctaccttagcatcccttccctttgcaaatagtcctcttccaacaataataatgtcagatcctgtagagaccacatcatccacggttctatactgttgacccaatgcgtctcccttgtcatctaaacccacaccgggtgtcataatcaaccaatcgtaaccttcatctcttccacccatgtctctttgagcaataaagccgataacaaaatctttgtcgctcttcgcaatgtcaacagtacccttagtatattctccagtagatagggagcccttgcatgacaattctgctaacatcaaaaggcctctaggttcctttgttacttcttctgccgcctgcttcaaaccgctaacaatacctgggcccaccacaccgtgtgcattcgtaatgtctgcccattctgctattctgtatacacccgcagagtactgcaatttgactgtattaccaatgtcagcaaattttctgtcttcgaagagtaaaaaattgtacttggcggataatgcctttagcggcttaactgtgccctccatggaaaaatcagtcaaaatatccacatgtgtttttagtaaacaaattttgggacctaatgcttcaactaactccagtaattccttggtggtacgaacatccaatgaagcacacaagtttgtttgcttttcgtgcatgatattaaatagcttggcagcaacaggactaggatgagtagcagcacgttccttatatgtagctttcgacatgatttatcttcgtttcctgcaggtttttgttctgtgcagttgggttaagaatactgggcaatttcatgtttcttcaacactacatatgcgtatatataccaatctaagtctgtgctccttccttcgttcttccttctgttcggagattaccgaatcaaaaaaatttcaaagaaaccgaaatcaaaaaaaagaataaaaaaaaaatgatgaattgaattgaaaagctagcttatcgatgataagctgtcaaagatgagaattaattccacggactatagactatactagatactccgtctactgtacgatacacttccgctcaggtccttgtcctttaacgaggccttaccactcttttgttactctattgatccagctcagcaaaggcagtgtgatctaagattctatcttcgcgatgtagtaaaactagctagaccgagaaagagactagaaatgcaaaaggcacttctacaatggctgccatcattattatccgatgtgacgctgcagcttctcaatgatattcgaatacgctttgaggagatacagcctaatatccgacaaactgttttacagatttacgatcgtacttgttacccatcattgaattttgaacatccgaacctgggagttttccctgaaacagatagtatatttgaacctgtataataatatatagtctagcgctttacggaagacaatgtatgtatttcggttcctggagaaactattgcatctattgcataggtaatcttgcacgtcgcatccccggttcattttctgcgtttccatcttgcacttcaatagcatatctttgttaacgaagcatctgtgcttcattttgtagaacaaaaatgcaacgcgagagcgctaatttttcaaacaaagaatctgagctgcatttttacagaacagaaatgcaacgcgaaagcgctattttaccaacgaagaatctgtgcttcatttttgtaaaacaaaaatgcaacgcgacgagagcgctaatttttcaaacaaagaatctgagctgcatttttacagaacagaaatgcaacgcgagagcgctattttaccaacaaagaatctatacttcttttttgttctacaaaaatgcatcccgagagcgctatttttctaacaaagcatcttagattactttttttctcctttgtgcgctctataatgcagtctcttgataactttttgcactgtaggtccgttaaggttagaagaaggctactttggtgtctattttctcttccataaaaaaagcctgactccacttcccgcgtttactgattactagcgaagctgcgggtgcattttttcaagataaaggcatccccgattatattctataccgatgtggattgcgcatactttgtgaacagaaagtgatagcgttgatgattcttcattggtcagaaaattatgaacggtttcttctattttgtctctatatactacgtataggaaatgtttacattttcgtattgttttcgattcactctatgaatagttcttactacaatttttttgtctaaagagtaatactagagataaacataaaaaatgtagaggtcgagtttagatgcaagttcaaggagcgaaaggtggatgggtaggttatatagggatatagcacagagatatatagcaaagagatacttttgagcaatgtttgtggaagcggtattcgcaatgggaagctccaccccggttgataatcagaaaagccccaaaaacaggaagattattatcaaaaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacagttaccaatgcttaatcagtgaggcacctatctcagcgatctgtctatttcgttcatccatagttgcctgactccccgtcgtgtagataactacgatacgggagcgcttaccatctggccccagtgctgcaatgataccgcgagacccacgctcaccggctccagatttatcagcaataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttggcattgctacaggcatcgtggtgtcactctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggataatagtgtatcacatagcagaactttaaaagtgctcatcattggaaaacgttcttcggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatgtaacccactcgtgcacccaactgatcttcagcatcttttactttcaccagcgtttctgggtgagcaaaaacaggaaggcaaaatgccgcaaaaaagggaataagggcgacacggaaatgttgaatactcatactcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaaaagtgccacctgctaagaaaccattattatcatgacattaacctataaaaataggcgtatcacgaggccctttcgtctcgcgcgtttcggtgatgacggtgaaaacctctgacacatgcagctcccggagacggtcacagcttgtctgtaagcggatgccgggagcagacaagcccgtcagggcgcgtcagcgggtgttggcgggtgtcggggctggcttaactatgcggcatcagagcagattgtactgagagtgcaccatagatcctgaggatcggggtgataaatcagtctgcgccacatcgggggaaacaaaatggcgcgagatctaaaaaaaaaggctccaaaaggagcctttcgcgctaccaggtaacgcgccactccgacgggattaacgagtgccgtaaacgacgatggttttaccgtgtgcggagatcaggttctgatcctcgagcatcttaagaattcgtcccacggtttgtctagagcagccgacaatctggccaatttcctgacgggtaattttgatttgcatgccgtccgggtgagtcatagcgtctggttgttttgccagattcagcagagtctgtgcaatgcggccgctgac
    """

    text6 = """
    >1671
    gaccaaattaactctatccttccattgcacaatttgcccagtatggatgtccacgaggtctctgctatggacacatttacgcccgtacgctgcaggtcgacggatccccgggttaattaaggcgcgccagatctgtttagcttgcctcgtccccgccgggtcacccggccagcgacatggaggcccagaataccctccttgacagtcttgacgtgcgcagctcaggggcatgatgtgactgtcgcccgtacatttagcccatacatccccatgtataatcatttgcatccatacattttgatggccgcacggcgcgaagcaaaaattacggctcctcgctgcagacctgcgagcagggaaacgctcccctcacagacgcgttgaattgtccccacgccgcgcccctgtagagaaatataaaaggttaggatttgccactgaggttcttctttcatatacttccttttaaaatcttgctaggatacagttctcacatcacatccgaacataaacaaccatgggtaaggaaaagactcacgtttcgaggccgcgattaaattccaacatggatgctgatttatatgggtataaatgggctcgcgataatgtcgggcaatcaggtgcgacaatctatcgattgtatgggaagcccgatgcgccagagttgtttctgaaacatggcaaaggtagcgttgccaatgatgttacagatgagatggtcagactaaactggctgacggaatttatgcctcttccgaccatcaagcattttatccgtactcctgatgatgcatggttactcaccactgcgatccccggcaaaacagcattccaggtattagaagaatatcctgattcaggtgaaaatattgttgatgcgctggcagtgttcctgcgccggttgcattcgattcctgtttgtaattgtccttttaacagcgatcgcgtatttcgtctcgctcaggcgcaatcacgaatgaataacggtttggttgatgcgagtgattttgatgacgagcgtaatggctggcctgttgaacaagtctggaaagaaatgcataagcttttgccattctcaccggattcagtcgtcactcatggtgatttctcacttgataaccttatttttgacgaggggaaattaataggttgtattgatgttggacgagtcggaatcgcagaccgataccaggatcttgccatcctatggaactgcctcggtgagttttctccttcattacagaaacggctttttcaaaaatatggtattgataatcctgatatgaataaattgcagtttcatttgatgctcgatgagtttttctaatcagtactgacaataaaaagattcttgttttcaagaacttgtcatttgtatagtttttttatattgtagttgttctattttaatcaaatgttagcgtgatttatattttttttcgcctcgacatcatctgcccagatgcgaagttaagtgcgcagaaagtaatatcatgcgtcaatcgtatgtgaatgctggtcgctatactgctgtcgattcgatactaacgccgccatccagtgtcgaaaacgagctcgaattcatcgatgaatgacgttatagacgagcctacgagaccgacaccgtaatcaataagcaattttttgcgaacatgattaagtttattttac

    >3527
    aaattgacctcgtttgattgctctcaggtttccaacaacttgtttgttaccggtttcagtacgggtattatcaagttatgggatgccagagcggctgaagctgccaccactgatttgacttatagacaaaacggtgaagatccaattcagaatgaaatcgctaacttctatcatgccggcggtgattctgttgtcgacgtccaattttctgccacttcaagctctgaattctttactgtgggcggcactggtaatatttaccactggaacactgactattccttgtcaaagtacaatcctgacgataccattgctcctcctcaagatgccactgaagaatcacaaacaaaatctctgagattcttgcacaagggtggaagcagaagatctccaaaacagattggaagaagaaacaccgccgcttggcatccagttattgaaaacttggttggtactgtcgatgacgatagtttagtcagcatctacaagccctataccgaggaaagcgaatagtcgccatgctaaacgcgcggaacaaggccatatttatatatttaatgcttttaactattattagtttttctcacccagggctggtttctttaccctgtgtgatgaaagtgcgcgcctaaagttatgtatgagctactgttagtatatttattattcaatatttaattaaacaatacttactactgatgaaagataccgtacacctgctggcgaataagtcactatacacgagcaattatatttctcaggtcattgcgacaattgaaaaatcggatcgactgattatcattaaggcgtgtgacggctcagaataaggaagttctattttaaagggcagcaaaaatttggcacaatttaaaaagaaagttggtctctggtagggttggatccggtttttcataacttgcgtttcatttatggttgcgtatttctcgcagtgtacatagaccaaattaactctatccttccattgcacaatttgcccagtatggcttcactgttcagacccccagaatctgcgaaatgcaacccaaactctcctagacttaaactgcccctcttacgaaataatcaggtagatgaaaataatatatacttgacgtcgaatgggagttccaccacagcttacagtagccacaccccagaaccactgacctcttccacatcgacacttttctcccaaactcgacttcatcctagcgactcttcaatgactttaaatacaatgaagaagaggcctgcaccgccatctttaccttcgctcagcataaattcacagtctaagtgcaaaacactacccgaactcgtacccatcgccgatgtgtctgatggtaaacatgatttaggattgaaacagcgtgtgatagccgaaaatgagttgtctggtaatagtgacttaaccccttcatcgatggcaagccccttttcacatacaaacacctcttctccctacctcagaaatgatctgagcaattcagtgggatctgacttttcaaatttgatatcggcatatgaacaaagttcaagtcccatcaagtcatcgtcccagcctaaatcatcttctgaatcgtacatagacttaaacagtgtacgagatgttgatcaattggatgaaaatggttggaaatatgcaaatttaaaagataggatcgagacattaggcattctaggagaaggagccggtggctctgtttccaagtgtaaattgaaaaatggatcaaaaatattcgctttaaaagtgataaacacattaaatacagatcccgagtatcagaagcaaatattcagagaattacagtttaataggagtttccaatccgaatatatcgtacgatattatggaatgtttacggatgacgaaaactcttcaatttatattgctatggagtacatgggtgggcgatcgttggatgctatttataaaaatttgttagagcgtggtggtaggatcagtgaaaaagtcctggggaagattgcagaagcggtactaagaggactatcatatttgcatgaaaaaaaagttattcatagagatattaagccccagaatattttactgaatgaaaatggtcaggtgaaactttgtgattttggggtcagtggagaagccgttaactcgctagccacaacattcacggggacgtcattctatatggctccagaaaggattcaaggtcaaccatacagtgtcacatctgatgtatggtcacttggattaacgattttggaagtggcgaacgggaaatttccatgctcttccgagaagatggcagccaatatagctccctttgaattgttaatgtggattttaacatttactcctgaattaaaggatgaacccgaatctaatatcatatggagtccatcattcaaatcctttatcgactactgtctgaaaaaagatagtcgtgaacggccatctccaagacaaatgatcaatcatccttggataaagggtcaaatgaaaaaaaatgtcaatatggaaaaattcgtgaggaagtgctggaaagattaatcaataagcaattttttgcgaacatgattaagtttattttactttatttcgcattttccataaaaaaattttcctccatacaagattcatacccaggatagcctaactaaaatgtatggctttatacagcttgacgacaaacttaatttgaatatatagatatattaatttaatcagcaggattagcatttaatgtcttacatgtctttgtatttggcatacgtattcagtgatattttaatacctcgcaacgcaatttccagcgagtttacgtttgaagcaaacgtaactcccatgagttatacataacgcctgctgcagctgaaccagaacaatataccgtggtattaccgtaagcatttgaaaacaataataagtttccaagtgccgctttagatccaaagtacaaggacatattgtggaagagtaagttgaatatttgaaaatgagagctttttcagcagccaccgttagggccacaactaggaagtcgttcatcccaatggcaccaagaactccttttgtgactccatcatttacaaagaatgtaggctcaatgagaagaatgagattttattctgatgaagccaaaagtgaagaatccaaagaaaacaatgaagatttgactgaagagcaatcagaaatcaagaaattagagagccagttaagcgcgaagactaaagaagcttctgaactcaaggacagattattaagatctgtggcagatttcagaaatttacaacaagtcacaaagaaggatattcagaaagctaaggactttgctttacagaagtttgcaaaggatttattggaatctgtagataactttggtcatgctttgaatgcttttaaagaggaagacttacaaaagtccaaggaaattagtgatttgtatacaggggttagaatgacaagagatgtttttgaaaacaccctaagaaagcacggtattgaaaaattagacccattgggagaaccatttgatccaaataaacacgaa
    """

    text7 = """
    >630
    tcctgacgggtaattttgatttgcatgccgtccgggtgagtcatagcgtctggttgttttgccagattcagcagagtctgtgcaatgcggccgctgacgtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaaggacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgagattatcaaaaaggatcttcacctagatgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgg

    >1196
    gtcgaggaacgccaggttgcccactttctcactagtgacctgcagccggccaatcacatcacatccgaacataaacaaccatgggtaaaaagcctgaactcaccgcgacgtctgtcgagaagtttctgatcgaaaagttcgacagcgtctccgacctgatgcagctctcggagggcgaagaatctcgtgctttcagcttcgatgtaggagggcgtggatatgtcctgcgggtaaatagctgcgccgatggtttctacaaagatcgttatgtttatcggcactttgcatcggccgcgctcccgattccggaagtgcttgacattggggaattcagcgagagcctgacctattgcatctcccgccgtgcacagggtgtcacgttgcaagacctgcctgaaaccgaactgcccgctgttctgcagccggtcgcggaggccatggatgcgatcgctgcggccgatcttagccagacgagcgggttcggcccattcggaccgcaaggaatcggtcaatacactacatggcgtgatttcatatgcgcgattgctgatccccatgtgtatcactggcaaactgtgatggacgacaccgtcagtgcgtccgtcgcgcaggctctcgatgagctgatgctttgggccgaggactgccccgaagtccggcacctcgtgcacgcggatttcggctccaacaatgtcctgacggacaatggccgcataacagcggtcattgactggagcgaggcgatgttcggggattcccaatacgaggtcgccaacatcttcttctggaggccgtggttggcttgtatggagcagcagacgcgctacttcgagcggaggcatccggagcttgcaggatcgccgcggctccgggcgtatatgctccgcattggtcttgaccaactctatcagagcttggttgacggcaatttcgatgatgcagcttgggcgcagggtcgatgcgacgcaatcgtccgatccggagccgggactgtcgggcgtacacaaatcgcccgcagaagcgcggccgtctggaccgatggctgtgtagaagtactcgccgatagtggaaaccgacgccccagcactcgtccgagggcaaaggaataatcagtactgacaataaaaagattcttgtagggataacagggtaatcggacgcgccatctgtgcagacaaacgcatcaggatttaaat

    >650
    cgcgccatctgtgcagacaaacgcatcaggatagagtcttttgtaacgaccccgtctccaccaacttggtatgcttgaaatctcaaggccattacacattcagttatgtgaacgaaaggtctttatttaacgtagcataaactaaataatacaggttccggttagcctgcaatgtgttaaatctaaaggagcatacccaaaatgaactgaagacaaggaaatttgcttgtccagatgtgattgagcatttgaacgttaataacataacatttttatacttaactatagaaagacttgtataaaaactggcaaacgagatattctgaatattggtgcatatttcaggtagaaaagcttacaaaacaatctaatcataatattgagatgaagagaaagataaaagaaaaaacgataagtcagatgagattatgattgtactttgaaatcgaggaacaaagtatatacggtagtagttccccgagttataacgggagatcatgtaaattgagaaaccagataaagatttggtatgcactctagcaagaaaataaaatgatgaatctatgatatagatcacttttgttccagcatccggatttacctgaatcaattggcgaaattttttgtacgaaatttcagccacttcacag

    >5681
    gtcagcggccgcattgcacagactctgctgaatctggcaaaacaaccagacgctatgactcacccggacggcatgcaaatcaaaattacccgtcaggaaattggccagattgtcggctgctctagacaaaccgtgggacgaattcttaagatgctcgaggatcagaacctgatctccgcacacggtaaaaccatcgtcgtttacggcactcgttaatcccgtcggagtggcgcgttacctggtagcgcgaaaggctccttttggagcctttttttttagatctcgcgccattttgtttcccccgatgtggcgcagactgatttatcaccccgatcctcaggatctatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtgatacactattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagagtgacaccacgatgcctgtagcaatgccaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagcgctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataataatcttcctgtttttggggcttttctgattatcaaccggggtggagcttcccattgcgaataccgcttccacaaacattgctcaaaagtatctctttgctatatatctctgtgctatatccctatataacctacccatccacctttcgctccttgaacttgcatctaaactcgacctctacattttttatgtttatctctagtattactctttagacaaaaaaattgtagtaagaactattcatagagtgaatcgaaaacaatacgaaaatgtaaacatttcctatacgtagtatatagagacaaaatagaagaaaccgttcataattttctgaccaatgaagaatcatcaacgctatcactttctgttcacaaagtatgcgcaatccacatcggtatagaatataatcggggatgcctttatcttgaaaaaatgcacccgcagcttcgctagtaatcagtaaacgcgggaagtggagtcaggctttttttatggaagagaaaatagacaccaaagtagccttcttctaaccttaacggacctacagtgcaaaaagttatcaagagactgcattatagagcgcacaaaggagaaaaaaagtaatctaagatgctttgttagaaaaatagcgctctcgggatgcatttttgtagaacaaaaaagaagtatagattctttgttggtaaaatagcgctctcgcgttgcatttctgttctgtaaaaatgcagctcagattctttgtttgaaaaattagcgctctcgtcgcgttgcatttttgttttacaaaaatgaagcacagattcttcgttggtaaaatagcgctttcgcgttgcatttctgttctgtaaaaatgcagctcagattctttgtttgaaaaattagcgctctcgcgttgcatttttgttctacaaaatgaagcacagatgcttcgttaacaaagatatgctattgaagtgcaagatggaaacgcagaaaatgaaccggggatgcgacgtgcaagattacctatgcaatagatgcaatagtttctccaggaaccgaaatacatacattgtcttccgtaaagcgctagactatatattattatacaggttcaaatatactatctgtttcagggaaaactcccaggttcggatgttcaaaattcaatgatgggtaacaagtacgatcgtaaatctgtaaaacagtttgtcggatattaggctgtatctcctcaaagcgtattcgaatatcattgagaagctgcagcgtcacatcggataataatgatggcagccattgtagaagtgccttttgcatttctagtctctttctcggtctagctagttttactacatcgcgaagatagaatcttagatcacactgcctttgctgagctggatcaatagagtaacaaaagagtggtaaggcctcgttaaaggacaaggacctgagcggaagtgtatcgtacagtagacggagtatctagtatagtctatagtccgtggaattaattctcatctttgacagcttatcatcgataagctagcttttcaattcaattcatcattttttttttattcttttttttgatttcggtttctttgaaatttttttgattcggtaatctccgaacagaaggaagaacgaaggaaggagcacagacttagattggtatatatacgcatatgtagtgttgaagaaacatgaaattgcccagtattcttaacccaactgcacagaacaaaaacctgcaggaaacgaagataaatcatgtcgaaagctacatataaggaacgtgctgctactcatcctagtcctgttgctgccaagctatttaatatcatgcacgaaaagcaaacaaacttgtgtgcttcattggatgttcgtaccaccaaggaattactggagttagttgaagcattaggtcccaaaatttgtttactaaaaacacatgtggatattttgactgatttttccatggagggcacagttaagccgctaaaggcattatccgccaagtacaattttttactcttcgaagacagaaaatttgctgacattggtaatacagtcaaattgcagtactctgcgggtgtatacagaatagcagaatgggcagacattacgaatgcacacggtgtggtgggcccaggtattgttagcggtttgaagcaggcggcagaagaagtaacaaaggaacctagaggccttttgatgttagcagaattgtcatgcaagggctccctatctactggagaatatactaagggtactgttgacattgcgaagagcgacaaagattttgttatcggctttattgctcaaagagacatgggtggaagagatgaaggttacgattggttgattatgacacccggtgtgggtttagatgacaagggagacgcattgggtcaacagtatagaaccgtggatgatgtggtctctacaggatctgacattattattgttggaagaggactatttgcaaagggaagggatgctaaggtagagggtgaacgttacagaaaagcaggctgggaagcatatttgagaagatgcggccagcaaaactaaaaaactgtattataagtaaatgcatgtatactaaactcacaaattagagcttcaatttaattatatcagttattacccctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagataaaaaatttatttgctttcaggtacaattcttgatataatattatcatctagctagataataaaaatttaaggatcttatggtacttggtaaacctcaaacggatcctactctcgaatggttcttgtctcattgccacattcataagtacccatccaagagcacgcttattcaccagggtgaaaaagcggaaacgctgtactacatcgttaaaggctctgtggcagtgctgatcaaagacgaagagggtaaagaaatgatcctctcctatctgaatcagggtgattttattggcgaactgggcctgtttgaagagggccaggaacgtagcgcatgggtacgtgcgaaaaccgcctgtgaagtggctgaaatttcgtacaaaaaatttcgccaattgattcaggtaaatccggat
    """

    text8 = """
    >5681 SEGUID JPjjayC2rXeAi4Jlju1WUlfqAvE
    gtcagcggccgcattgcacagactctgctgaatctggcaaaacaaccagacgctatgactcacccggacggcatgcaaatcaaaattacccgtcaggaaattggccagattgtcggctgctctagacaaaccgtgggacgaattcttaagatgctcgaggatcagaacctgatctccgcacacggtaaaaccatcgtcgtttacggcactcgttaatcccgtcggagtggcgcgttacctggtagcgcgaaaggctccttttggagcctttttttttagatctcgcgccattttgtttcccccgatgtggcgcagactgatttatcaccccgatcctcaggatctatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtgatacactattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagagtgacaccacgatgcctgtagcaatgccaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagcgctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataataatcttcctgtttttggggcttttctgattatcaaccggggtggagcttcccattgcgaataccgcttccacaaacattgctcaaaagtatctctttgctatatatctctgtgctatatccctatataacctacccatccacctttcgctccttgaacttgcatctaaactcgacctctacattttttatgtttatctctagtattactctttagacaaaaaaattgtagtaagaactattcatagagtgaatcgaaaacaatacgaaaatgtaaacatttcctatacgtagtatatagagacaaaatagaagaaaccgttcataattttctgaccaatgaagaatcatcaacgctatcactttctgttcacaaagtatgcgcaatccacatcggtatagaatataatcggggatgcctttatcttgaaaaaatgcacccgcagcttcgctagtaatcagtaaacgcgggaagtggagtcaggctttttttatggaagagaaaatagacaccaaagtagccttcttctaaccttaacggacctacagtgcaaaaagttatcaagagactgcattatagagcgcacaaaggagaaaaaaagtaatctaagatgctttgttagaaaaatagcgctctcgggatgcatttttgtagaacaaaaaagaagtatagattctttgttggtaaaatagcgctctcgcgttgcatttctgttctgtaaaaatgcagctcagattctttgtttgaaaaattagcgctctcgtcgcgttgcatttttgttttacaaaaatgaagcacagattcttcgttggtaaaatagcgctttcgcgttgcatttctgttctgtaaaaatgcagctcagattctttgtttgaaaaattagcgctctcgcgttgcatttttgttctacaaaatgaagcacagatgcttcgttaacaaagatatgctattgaagtgcaagatggaaacgcagaaaatgaaccggggatgcgacgtgcaagattacctatgcaatagatgcaatagtttctccaggaaccgaaatacatacattgtcttccgtaaagcgctagactatatattattatacaggttcaaatatactatctgtttcagggaaaactcccaggttcggatgttcaaaattcaatgatgggtaacaagtacgatcgtaaatctgtaaaacagtttgtcggatattaggctgtatctcctcaaagcgtattcgaatatcattgagaagctgcagcgtcacatcggataataatgatggcagccattgtagaagtgccttttgcatttctagtctctttctcggtctagctagttttactacatcgcgaagatagaatcttagatcacactgcctttgctgagctggatcaatagagtaacaaaagagtggtaaggcctcgttaaaggacaaggacctgagcggaagtgtatcgtacagtagacggagtatctagtatagtctatagtccgtggaattaattctcatctttgacagcttatcatcgataagctagcttttcaattcaattcatcattttttttttattcttttttttgatttcggtttctttgaaatttttttgattcggtaatctccgaacagaaggaagaacgaaggaaggagcacagacttagattggtatatatacgcatatgtagtgttgaagaaacatgaaattgcccagtattcttaacccaactgcacagaacaaaaacctgcaggaaacgaagataaatcatgtcgaaagctacatataaggaacgtgctgctactcatcctagtcctgttgctgccaagctatttaatatcatgcacgaaaagcaaacaaacttgtgtgcttcattggatgttcgtaccaccaaggaattactggagttagttgaagcattaggtcccaaaatttgtttactaaaaacacatgtggatattttgactgatttttccatggagggcacagttaagccgctaaaggcattatccgccaagtacaattttttactcttcgaagacagaaaatttgctgacattggtaatacagtcaaattgcagtactctgcgggtgtatacagaatagcagaatgggcagacattacgaatgcacacggtgtggtgggcccaggtattgttagcggtttgaagcaggcggcagaagaagtaacaaaggaacctagaggccttttgatgttagcagaattgtcatgcaagggctccctatctactggagaatatactaagggtactgttgacattgcgaagagcgacaaagattttgttatcggctttattgctcaaagagacatgggtggaagagatgaaggttacgattggttgattatgacacccggtgtgggtttagatgacaagggagacgcattgggtcaacagtatagaaccgtggatgatgtggtctctacaggatctgacattattattgttggaagaggactatttgcaaagggaagggatgctaaggtagagggtgaacgttacagaaaagcaggctgggaagcatatttgagaagatgcggccagcaaaactaaaaaactgtattataagtaaatgcatgtatactaaactcacaaattagagcttcaatttaattatatcagttattacccctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagataaaaaatttatttgctttcaggtacaattcttgatataatattatcatctagctagataataaaaatttaaggatcttatggtacttggtaaacctcaaacggatcctactctcgaatggttcttgtctcattgccacattcataagtacccatccaagagcacgcttattcaccagggtgaaaaagcggaaacgctgtactacatcgttaaaggctctgtggcagtgctgatcaaagacgaagagggtaaagaaatgatcctctcctatctgaatcagggtgattttattggcgaactgggcctgtttgaagagggccaggaacgtagcgcatgggtacgtgcgaaaaccgcctgtgaagtggctgaaatttcgtacaaaaaatttcgccaattgattcaggtaaatccggat

    >630 SEGUID 3GE7aiThG66A-DHF-12fiDtQJ0k
    TCCTGACGGGTAATTTTGATTTGCATGCCGTCCGGGTGAGTCATAGCGTCTGGTTGTTTTGCCAGATTCAGCAGAGTCTGTGCAATGCGGCCGCTGACCACATACGATTTAGGTGACACTATAGAACGCGGCCGCCAGCTGAAGCTTCGTACGCTGCAGGTCGACGGATCCCCGGGTTAATTAAGGCGCGCCAGATCTGTTTAGCTTGCCTTGTCCCCGCCGGGTCACCCGGCCAGCGACATGGAGGCCCAGAATACCCTCCTTGACAGTCTTGACGTGCGCAGCTCAGGGGCATGATGTGACTGTCGCCCGTACATTTAGCCCATACATCCCCATGTATAATCATTTGCATCCATACATTTTGATGGCCGCACGGCGCGAAGCAAAAATTACGGCTCCTCGCTGCAGACCTGCGAGCAGGGAAACGCTCCCCTCACAGACGCGTTGAATTGTCCCCACGCCGCGCCCCTGTAGAGAAATATAAAAGGTTAGGATTTGCCACTGAGGTTCTTCTTTCATATACTTCCTTTTAAAATCTTGCTAGGATACAGTTCTCACATCACATCCGAACATAAACAACCGTCGAGGAACGCCAGGTTGCCCACTTTCTCACTAGTGACCTGCAGCCGG

    >1196 SEGUID Bl-131D4eiLsj4bR3f6fue0gXik
    gtcgaggaacgccaggttgcccactttctcactagtgacctgcagccggccaatcacatcacatccgaacataaacaaccatgggtaaaaagcctgaactcaccgcgacgtctgtcgagaagtttctgatcgaaaagttcgacagcgtctccgacctgatgcagctctcggagggcgaagaatctcgtgctttcagcttcgatgtaggagggcgtggatatgtcctgcgggtaaatagctgcgccgatggtttctacaaagatcgttatgtttatcggcactttgcatcggccgcgctcccgattccggaagtgcttgacattggggaattcagcgagagcctgacctattgcatctcccgccgtgcacagggtgtcacgttgcaagacctgcctgaaaccgaactgcccgctgttctgcagccggtcgcggaggccatggatgcgatcgctgcggccgatcttagccagacgagcgggttcggcccattcggaccgcaaggaatcggtcaatacactacatggcgtgatttcatatgcgcgattgctgatccccatgtgtatcactggcaaactgtgatggacgacaccgtcagtgcgtccgtcgcgcaggctctcgatgagctgatgctttgggccgaggactgccccgaagtccggcacctcgtgcacgcggatttcggctccaacaatgtcctgacggacaatggccgcataacagcggtcattgactggagcgaggcgatgttcggggattcccaatacgaggtcgccaacatcttcttctggaggccgtggttggcttgtatggagcagcagacgcgctacttcgagcggaggcatccggagcttgcaggatcgccgcggctccgggcgtatatgctccgcattggtcttgaccaactctatcagagcttggttgacggcaatttcgatgatgcagcttgggcgcagggtcgatgcgacgcaatcgtccgatccggagccgggactgtcgggcgtacacaaatcgcccgcagaagcgcggccgtctggaccgatggctgtgtagaagtactcgccgatagtggaaaccgacgccccagcactcgtccgagggcaaaggaataatcagtactgacaataaaaagattcttgtagggataacagggtaatcggacgcgccatctgtgcagacaaacgcatcaggatttaaat

    >650 SEGUID 8lAwzHM60BkV1hhzx3_bBsfAIYo
    cgcgccatctgtgcagacaaacgcatcaggatagagtcttttgtaacgaccccgtctccaccaacttggtatgcttgaaatctcaaggccattacacattcagttatgtgaacgaaaggtctttatttaacgtagcataaactaaataatacaggttccggttagcctgcaatgtgttaaatctaaaggagcatacccaaaatgaactgaagacaaggaaatttgcttgtccagatgtgattgagcatttgaacgttaataacataacatttttatacttaactatagaaagacttgtataaaaactggcaaacgagatattctgaatattggtgcatatttcaggtagaaaagcttacaaaacaatctaatcataatattgagatgaagagaaagataaaagaaaaaacgataagtcagatgagattatgattgtactttgaaatcgaggaacaaagtatatacggtagtagttccccgagttataacgggagatcatgtaaattgagaaaccagataaagatttggtatgcactctagcaagaaaataaaatgatgaatctatgatatagatcacttttgttccagcatccggatttacctgaatcaattggcgaaattttttgtacgaaatttcagccacttcacag
    """

    # text1

    list_of_formatted_seq_records = parse(text1)
    a = assembly.Assembly(list_of_formatted_seq_records, limit=25)

    assert (
        repr(a)
        == """Assembly
fragments..: 631bp 740bp 650bp
limit(bp)..: 25
G.nodes....: 6
algorithm..: common_sub_strings"""
    )

    candidate = a.assemble_linear()[0]
    correct = "TCCTGACGGGTAATTTTGATTTGCATGCCGTCCGGGTGAGTCATAGCGTCTGGTTGTTTTGCCAGATTCAGCAGAGTCTGTGCAATGCGGCCGCTGACcacatacgatttaggtgacactatagaacGCGGCCGCCAGCTGAAGCTTCGTACGCTGCAGGTCGACGGATCCCCGGGTTAATTAAGGCGCGCCAGATCTGTTTAGCTTGCCTTGTCCCCGCCGGGTCACCCGGCCAGCGACATGGAGGCCCAGAATACCCTCCTTGACAGTCTTGACGTGCGCAGCTCAGGGGCATGATGTGACTGTCGCCCGTACATTTAGCCCATACATCCCCATGTATAATCATTTGCATCCATACATTTTGATGGCCGCACGGCGCGAAGCAAAAATTACGGCTCCTCGCTGCAGACCTGCGAGCAGGGAAACGCTCCCCTCACAGACGCGTTGAATTGTCCCCACGCCGCGCCCCTGTAGAGAAATATAAAAGGTTAGGATTTGCCACTGAGGTTCTTCTTTCATATACTTCCTTTTAAAATCTTGCTAGGATACAGTTCTCACATcacatccgaacataaacaaccGTCGAGGAACGCCAGGTTGCCCACTTTctcactagtgacctgcagccgACCCAAtcacatcacatccgaacataaacaaccatgGGTACCACTCTTGACGACACGGCTTACCGGTACCGCACCAGTGTCCCGGGGGACGCCGAGGCCATCGAGGCACTGGATGGGTCCTTCACCACCGACACCGTCTTCCGCGTCACCGCCACCGGGGACGGCTTCACCCTGCGGGAGGTGCCGGTGGACCCGCCCCTGACCAAGGTGTTCCCCGACGACGAATCGGACGACGAATCGGACGACGGGGAGGACGGCGACCCGGACTCCCGGACGTTCGTCGCGTACGGGGACGACGGCGACCTGGCGGGCTTCGTGGTCGTCTCGTACTCCGGCTGGAACCGCCGGCTGACCGTCGAGGACATCGAGGTCGCCCCGGAGCACCGGGGGCACGGGGTCGGGCGCGCGTTGATGGGGCTCGCGACGGAGTTCGCCCGCGAGCGGGGCGCCGGGCACCTCTGGCTGGAGGTCACCAACGTCAACGCACCGGCGATCCACGCGTACCGGCGGATGGGGTTCACCCTCTGCGGCCTGGACACCGCCCTGTACGACGGCACCGCCTCGGACGGCGAGCAGGCGCTCTACATGAGCATGCCCTGCCCCtaatcagtactgacaataaaaagattcttgTAGGGATAACAGGGTAATCGGAGTGCCATCTGTGCAGACAAACGCATCAGGATagagtcttttgtaacgacccCGTCTCCACCAACTTGGTATGCTTGAAATCTCAAGGCCATTACACATTCAGTTATGTGAACGAAAGGTCTTTATTTAACGTAGCATAAACTAAATAATACAGGTTCCGGTTAGCCTGCAATGTGTTAAATCTAAAGGAGCATACCCAAAATGAACTGAAGACAAGGAAATTTGCTTGTCCAGATGTGATTGAGCATTTGAACGTTAATAACATAACATTTTTATACTTAACTATAGAAAGACTTGTATAAAAACTGGCAAACGAGATATTCTGAATATTGGTGCATATTTCAGGTAGAAAAGCTTACAAAACAATCTAATCATAATATTGAGATGAAGAGAAAGATAAAAGAAAAAACGATAAGTCAGATGAGATTATGATTGTACTTTGAAATCGAGGAACAAAGTATATACGGTAGTAGTTCCCCGAGTTATAACGGGAGATCATGTAAATTGAGAAACCAGATAAAGATTTGGTATGCACTCTAGCAAGAAAATAAAATGATGAATCTATGAtatagatcacttttgttccagcATCCGGATTTACCTGAATCAATTGGCGAAATTTTTTGTACGAAATTTCAGCCACTTCACAG"
    assert len(correct) == 1933
    assert eq(correct, candidate, circular=True)

    # text2

    list_of_formatted_seq_records = parse(text2)
    a = assembly.Assembly(list_of_formatted_seq_records, limit=25)

    candidate = a.assemble_circular()[0]
    correct = "TTCTAGAACTAGTGGATCCCCCGGGCTGCAGATGAGTGAAGGCCCCGTCAAATTCGAAAAAAATACCGTCATATCTGTCTTTGGTGCGTCAGGTGATCTGGCAAAGAAGAAGACTTTTCCCGCCTTATTTGGGCTTTTCAGAGAAGGTTACCTTGATCCATCTACCAAGATCTTCGGTTATGCCCGGTCCAAATTGTCCATGGAGGAGGACCTGAAGTCCCGTGTCCTACCCCACTTGAAAAAACCTCACGGTGAAGCCGATGACTCTAAGGTCGAACAGTTCTTCAAGATGGTCAGCTACATTTCGGGAAATTACGACACAGATGAAGGCTTCGACGAATTAAGAACGCAGATCGAGAAATTCGAGAAAAGTGCCAACGTCGATGTCCCACACCGTCTCTTCTATCTGGCCTTGCCGCCAAGCGTTTTTTTGACGGTGGCCAAGCAGATCAAGAGTCGTGTGTACGCAGAGAATGGCATCACCCGTGTAATCGTAGAGAAACCTTTCGGCCACGACCTGGCCTCTGCCAGGGAGCTGCAAAAAAACCTGGGGCCCCTCTTTAAAGAAGAAGAGTTGTACAGAATTGACCATTACTTGGGTAAAGAGTTGGTCAAGAATCTTTTAGTCTTGAGGTTCGGTAACCAGTTTTTGAATGCCTCGTGGAATAGAGACAACATTCAAAGCGTTCAGATTTCGTTTAAAGAGAGGTTCGGCACCGAAGGCCGTGGCGGCTATTTCGACTCTATAGGCATAATCAGAGACGTGATGCAGAACCATCTGTTACAAATCATGACTCTCTTGACTATGGAAAGACCGGTGTCTTTTGACCCGGAATCTATTCGTGACGAAAAGGTTAAGGTTCTAAAGGCCGTGGCCCCCATCGACACGGACGACGTCCTCTTGGGCCAGTACGGTAAATCTGAGGACGGGTCTAAGCCCGCCTACGTGGATGATGACACTGTAGACAAGGACTCTAAATGTGTCACTTTTGCAGCAATGACTTTCAACATCGAAAACGAGCGTTGGGAGGGCGTCCCCATCATGATGCGTGCCGGTAAGGCTTTGAATGAGTCCAAGGTGGAGATCAGACTGCAGTACAAAGCGGTCGCATCGGGTGTCTTCAAAGACATTCCAAATAACGAACTGGTCATCAGAGTGCAGCCCGATGCCGCTGTGTACCTAAAGTTTAATGCTAAGACCCCTGGTCTGTCAAATGCTACCCAAGTCACAGATCTGAATCTAACTTACGCAAGCAGGTACCAAGACTTTTGGATTCCAGAGGCTTACGAGGTGTTGATAAGAGACGCCCTACTGGGTGACCATTCCAACTTTGTCAGAGATGACGAATTGGATATCAGTTGGGGCATATTCACCCCATTACTGAAGCACATAGAGCGTCCGGACGGTCCAACACCGGAAATTTACCCCTACGGATCAAGAGGTCCAAAGGGATTGAAGGAATATATGCAAAAACACAAGTATGTTATGCCCGAAAAGCACCCTTACGCTTGGCCCGTGACTAAGCCAGAAGATACGAAGGATAATTAGCTGCAGGAATTCGATATCAAGCTTATCGATACCGTCGACCTCGAGTCATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCCCCCCACATCCGCTCTAACCGAAAAGGAAGGAGTTAGACAACCTGAAGTCTAGGTCCCTATTTATTTTTTTATAGTTATGTTAGTATTAAGAACGTTATTTATATTTCAAATTTTTCTTTTTTTTCTGTACAGACGCGTGTACGCATGTAACATTATACTGAAAACCTTGCTTGAGAAGGTTTTGGGACGCTCGAAGGCTTTAATTTGCGGCCGGTACCCAATTCGCCCTATAGTGAGTCGTATTACGCGCGCTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCGACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTCTCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTAGGGTTCCGATTTAGTGCTTTACGGCACCTCGACCCCAAAAAACTTGATTAGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATAGACGGTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGTCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAACGCGAATTTTAACAAAATATTAACGTTTACAATTTCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATCGACGGTCGAGGAGAACTTCTAGTATATCCACATACCTAATATTATTGCCTTATTAAAAATGGAATCCCAACAATTACATCAAAATCCACATTCTCTTCAAAATCAATTGTCCTGTACTTCCTTGTTCATGTGTGTTCAAAAACGTTATATTTATAGGATAATTATACTCTATTTCTCAACAAGTAATTGGTTGTTTGGCCGAGCGGTCTAAGGCGCCTGATTCAAGAAATATCTTGACCGCAGTTAACTGTGGGAATACTCAGGTATCGTAAGATGCAAGAGTTCGAATCTCTTAGCAACCATTATTTTTTTCCTCAACATAACGAGAACACACAGGGGCGCTATCGCACAGAATCAAATTCGATGACTGGAAATTTTTTGTTAATTTCAGAGGTCGCCTGACGCATATACCTTTTTCAACTGAAAAATTGGGAGAAAAAGGAAAGGTGAGAGGCCGGAACCGGCTTTTCATATAGAATAGAGAAGCGTTCATGACTAAATGCTTGCATCACAATACTTGAAGTTGACAATATTATTTAAGGACCTATTGTTTTTTCCAATAGGTGGTTAGCAATCGTCTTACTTTCTAACTTTTCTTACCTTTTACATTTCAGCAATATATATATATATTTCAAGGATATACCATTCTAATGTCTGCCCCTATGTCTGCCCCTAAGAAGATCGTCGTTTTGCCAGGTGACCACGTTGGTCAAGAAATCACAGCCGAAGCCATTAAGGTTCTTAAAGCTATTTCTGATGTTCGTTCCAATGTCAAGTTCGATTTCGAAAATCATTTAATTGGTGGTGCTGCTATCGATGCTACAGGTGTCCCACTTCCAGATGAGGCGCTGGAAGCCTCCAAGAAGGTTGATGCCGTTTTGTTAGGTGCTGTGGCTGGTCCTAAATGGGGTACCGGTAGTGTTAGACCTGAACAAGGTTTACTAAAAATCCGTAAAGAACTTCAATTGTACGCCAACTTAAGACCATGTAACTTTGCATCCGACTCTCTTTTAGACTTATCTCCAATCAAGCCACAATTTGCTAAAGGTACTGACTTCGTTGTTGTCAGAGAATTAGTGGGAGGTATTTACTTTGGTAAGAGAAAGGAAGACGATGGTGATGGTGTCGCTTGGGATAGTGAACAATACACCGTTCCAGAAGTGCAAAGAATCACAAGAATGGCCGCTTTCATGGCCCTACAACATGAGCCACCATTGCCTATTTGGTCCTTGGATAAAGCTAATCTTTTGGCCTCTTCAAGATTATGGAGAAAAACTGTGGAGGAAACCATCAAGAACGAATTCCCTACATTGAAGGTTCAACATCAATTGATTGATTCTGCCGCCATGATCCTAGTTAAGAACCCAACCCACCTAAATGGTATTATAATCACCAGCAACATGTTTGGTGATATCATCTCCGATGAAGCCTCCGTTATCCCAGGTTCCTTGGGTTTGTTGCCATCTGCGTCCTTGGCCTCTTTGCCAGACAAGAACACCGCATTTGGTTTGTACGAACCATGCCACGGTTCTGCTCCAGATTTGCCAAAGAATAAGGTTGACCCTATCGCCACTATCTTGTCTGCTGCAATGATGTTGAAATTGTCATTGAACTTGCCTGAAGAAGGTAAGGCCATTGAAGATGCAGTTAAAAAGGTTTTGGATGCAGGTATCAGAACTGGTGATTTAGGTGGTTCCAACAGTACCACCGAAGTCGGTGATGCTGTCGCCGAAGAAGTTAAGAAAATCCTTGCTTAAAAAGATTCTCTTTTTTTATGATATTTGTACATAAACTTTATAAATGAAATTCATAATAGAAACGACACGAAATTACAAAATGGAATATGTTCATAGGGTAGACGAAACTATATACGCAATCTACATACATTTATCAAGAAGGAGAAAAAGGAGGATAGTAAAGGAATACAGGTAAGCAAATTGATACTAATGGCTCAACGTGATAAGGAAAAAGAATTGCACTTTAACATTAATATTGACAAGGAGGAGGGCACCACACAAAAAGTTAGGTGTAACAGAAAATCATGAAACTACGATTCCTAATTTGATATTGGAGGATTTTCTCTAAAAAAAAAAAAATACAACAAATAAAAAACACTCAATGACCTGACCATTTGATGGAGTTTAAGTCAATACCTTCTTGAAGCATTTCCCATAATGGTGAAAGTTCCCTCAAGAATTTTACTCTGTCAGAAACGGCCTTACGACGTAGTCGATATGGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAGTTAAGCCAGCCCCGACACCCGCCAACACCCGCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGGCATCCGCTTACAGACAAGCTGTGACCGTCTCCGGGAGCTGCATGTGTCAGAGGTTTTCACCGTCATCACCGAAACGCGCGAGACGAAAGGGCCTCGTGATACGCCTATTTTTATAGGTTAATGTCATGATAATAATGGTTTCTTAGTATGATCCAATATCAAAGGAAATGATAGCATTGAAGGATGAGACTAATCCAATTGAGGAGTGGCAGCATATAGAACAGCTAAAGGGTAGTGCTGAAGGAAGCATACGATACCCCGCATGGAATGGGATAATATCACAGGAGGTACTAGACTACCTTTCATCCTACATAAATAGACGCATATAAGTACGCATTTAAGCATAAACACGCACTATGCCGTTCTTCTCATGTATATATATATACAGGCAACACGCAGATATAGGTGCGACGTGAACAGTGAGCTGTATGTGCGCAGCTCGCGTTGCATTTTCGGAAGCGCTCGTTTTCGGAAACGCTTTGAAGTTCCTATTCCGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAGAGCGCTTTTGAAAACCAAAAGCGCTCTGAAGACGCACTTTCAAAAAACCAAAAACGCACCGGACTGTAACGAGCTACTAAAATATTGCGAATACCGCTTCCACAAACATTGCTCAAAAGTATCTCTTTGCTATATATCTCTGTGCTATATCCCTATATAACCTACCCATCCACCTTTCGCTCCTTGAACTTGCATCTAAACTCGACCTCTACATTTTTTATGTTTATCTCTAGTATTACTCTTTAGACAAAAAAATTGTAGTAAGAACTATTCATAGAGTGAATCGAAAACAATACGAAAATGTAAACATTTCCTATACGTAGTATATAGAGACAAAATAGAAGAAACCGTTCATAATTTTCTGACCAATGAAGAATCATCAACGCTATCACTTTCTGTTCACAAAGTATGCGCAATCCACATCGGTATAGAATATAATCGGGGATGCCTTTATCTTGAAAAAATGCACCCGCAGCTTCGCTAGTAATCAGTAAACGCGGGAAGTGGAGTCAGGCTTTTTTTATGGAAGAGAAAATAGACACCAAAGTAGCCTTCTTCTAACCTTAACGGACCTACAGTGCAAAAAGTTATCAAGAGACTGCATTATAGAGCGCACAAAGGAGAAAAAAAGTAATCTAAGATGCTTTGTTAGAAAAATAGCGCTCTCGGGATGCATTTTTGTAGAACAAAAAAGAAGTATAGATTCTTTGTTGGTAAAATAGCGCTCTCGCGTTGCATTTCTGTTCTGTAAAAATGCAGCTCAGATTCTTTGTTTGAAAAATTAGCGCTCTCGCGTTGCATTTTTGTTTTACAAAAATGAAGCACAGATTCTTCGTTGGTAAAATAGCGCTTTCGCGTTGCATTTCTGTTCTGTAAAAATGCAGCTCAGATTCTTTGTTTGAAAAATTAGCGCTCTCGCGTTGCATTTTTGTTCTACAAAATGAAGCACAGATGCTTCGTTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTCCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAACGCAATTAATGTGAGTTACCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCCTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCACACAGGAAACAGCTATGACCATGATTACGCCAAGCGCGCAATTAACCCTCACTAAAGGGAACAAAAGCTGGAGCTCAGTTTATCATTATCAATACTCGCCATTTCAAAGAATACGTAAATAATTAATAGTAGTGATTTTCCTAACTTTATTTAGTCAAAAAATTAGCCTTTTAATTCTGCTGTAACCCGTACATGCCCAAAATAGGGGGCGGGTTACACAGAATATATAACATCGTAGGTGTCTGGGTGAACAGTTTATTCCTGGCATCCACTAAATATAATGGAGCCCGCTTTTTAAGCTGGCATCCAGAAAAAAAAAGAATCCCAGCACCAAAATATTGTTTTCTTCACCAACCATCAGTTCATAGGTCCATTCTCTTAGCGCAACTACAGAGAACAGGGGCACAAACAGGCAAAAAACGGGCACAACCTCAATGGAGTGATGCAACCTGCCTGGAGTAAATGATGACACAAGGCAATTGACCCACGCATGTATCTATCTCATTTTCTTACACCTTCTATTACCTTCTGCTCTCTCTGATTTGGAAAAAGCTGAAAAAAAAGGTTGAAACCAGTTCCCTGAAATTATTCCCCTACTTGACTAATAAGTATATAAAGACGGTAGGTATTGATTGTAATTCTGTAAATCTATTTCTTAAACTTCTTAAATTCTACTTTTATAGTTAGTCTTTTTTTTAGTTTTAAAACACCAGAACTTAGTTTCGACGGA"
    assert len(correct) == 9253
    assert eq(correct, candidate, circular=True)

    # text3
    y, x = parse(text3)
    a = assembly.Assembly((x, y, x), limit=25)

    candidate = a.assemble_linear()[0]
    correct = "GAGGCACCAGCGTCAGCATTTTCAAAGGTGTGTTCTTCGTCAGACATGTTTTAGTGTGTGAATGAAATAGGTGTATGTTTTCTTTTTGCTAGACAATAATTAGGAACAAGGTAAGGGAACTAAAGTGTAGAATAAGATTAAAAAAGAAGAACAAGTTGAAAAGGCAAGTTGAAATTTCAAGAAAAAAGTCAATTGAAGTACAGTAAATTGACCTGAATATATCTGAGTTCCGACAACAATGAGTTTACCAAAGAGAACAATGGAATAGGAAACTTTGAACGAAGAAAGGAAAGCAGGAAAGGAAAAAATTTTTAGGCTCGAGAACAATAGGGCGAAAAAACAGGCAACGAACGAACAATGGAAAAACGAAAAAAAAAAAAAAAAACACAGAAAAGAATGCAGAAAGATGTCAACTGAAAAAAAAAAAGGTGAACACAGGAAAAAAAATAAAAAAAAAAAAAAAAAAAGGAGGACGAAACAAAAAAGTGAAAAAAAATGAAAATTTTTTTGGAAAACCAAGAAATGAATTATATTTCCGTGTGAGACGACATCGTCGAATATGATTCAGGGTAACAGTATTGATGTAATCAATTTCCTACCTGAATCTAAAATTCCCGGGAGCAAGATCAAGATGTTTTCACCGATCTTTCCGGTCTCTTTGGCCGGGGTTTACGGACGATGGCAGAAGACCAAAGCGCCAGTTCATTTGGCGAGCGTTGGTTGGTGGATCAAGCCCACGCGTAGGCAATCCTCGAGCAGATCCGCCAGGCGTGTATATATAGCGTGGATGGCCAGGCAACTTTAGTGCTGACACATACAGGCATATATATATGTGTGCGACGACACATGATCATATGGCATGCATGTGCTCTGTATGTATATAAAACTCTTGTTTTCTTCTTTTCTCTAAATATTCTTTCCTTATACATTAGGACCTTTGCAGCATAAATTACTATACTTCTATAGACACACAAACACAAATACACACACTAAATTAATAATGGATGTCCACGAGGTCTCTATATCGGGATCAGCCTGCCTCGTACGCTGCAGGTCGACGGATCCCCGGGTTAATTAAGGCGCGCCAGATCTGTTTAGCTTGCCTCGTCCCCGCCGGGTCACCCGGCCAGCGACATGGAGGCCCAGAATACCCTCCTTGACAGTCTTGACGTGCGCAGCTCAGGGGCATGATGTGACTGTCGCCCGTACATTTAGCCCATACATCCCCATGTATAATCATTTGCATCCATACATTTTGATGGCCGCACGGCGCGAAGCAAAAATTACGGCTCCTCGCTGCAGACCTGCGAGCAGGGAAACGCTCCCCTCACAGACGCGTTGAATTGTCCCCACGCCGCGCCCCTGTAGAGAAATATAAAAGGTTAGGATTTGCCACTGAGGTTCTTCTTTCATATACTTCCTTTTAAAATCTTGCTAGGATACAGTTCTCACATCACATCCGAACATAAACAACCATGGGTAAGGAAAAGACTCACGTTTCGAGGCCGCGATTAAATTCCAACATGGATGCTGATTTATATGGGTATAAATGGGCTCGCGATAATGTCGGGCAATCAGGTGCGACAATCTATCGATTGTATGGGAAGCCCGATGCGCCAGAGTTGTTTCTGAAACATGGCAAAGGTAGCGTTGCCAATGATGTTACAGATGAGATGGTCAGACTAAACTGGCTGACGGAATTTATGCCTCTTCCGACCATCAAGCATTTTATCCGTACTCCTGATGATGCATGGTTACTCACCACTGCGATCCCCGGCAAAACAGCATTCCAGGTATTAGAAGAATATCCTGATTCAGGTGAAAATATTGTTGATGCGCTGGCAGTGTTCCTGCGCCGGTTGCATTCGATTCCTGTTTGTAATTGTCCTTTTAACAGCGATCGCGTATTTCGTCTCGCTCAGGCGCAATCACGAATGAATAACGGTTTGGTTGATGCGAGTGATTTTGATGACGAGCGTAATGGCTGGCCTGTTGAACAAGTCTGGAAAGAAATGCATAAGCTTTTGCCATTCTCACCGGATTCAGTCGTCACTCATGGTGATTTCTCACTTGATAACCTTATTTTTGACGAGGGGAAATTAATAGGTTGTATTGATGTTGGACGAGTCGGAATCGCAGACCGATACCAGGATCTTGCCATCCTATGGAACTGCCTCGGTGAGTTTTCTCCTTCATTACAGAAACGGCTTTTTCAAAAATATGGTATTGATAATCCTGATATGAATAAATTGCAGTTTCATTTGATGCTCGATGAGTTTTTCTAATCAGTACTGACAATAAAAAGATTCTTGTTTTCAAGAACTTGTCATTTGTATAGTTTTTTTATATTGTAGTTGTTCTATTTTAATCAAATGTTAGCGTGATTTATATTTTTTTTCGCCTCGACATCATCTGCCCAGATGCGAAGTTAAGTGCGCAGAAAGTAATATCATGCGTCAATCGTATGTGAATGCTGGTCGCTATACTGCTGTCGATTCGATACTAACGCCGCCATCCAGTGTCGAAAACGAGCTCGAATTCATCGATGAGCATCCTTATAGCCTCTTCTACGAGACCGACACCGTAAACAGGCCCCTTTTCCTTTGTCGATATCATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCCTCCCACATCCGCTCTAACCGAAAAGGAAGGAGTTAGACAACCTGAAGTCTAGGTCCCTATTTATTTTTTTTAATAGTTATGTTAGTATTAAGAACGTTATTTATATTTCAAATTTTTCTTTTTTTTCTGTACAAACGCGTGTACGCATGTAACATTATACTGAAAACCTTGCTTGAGAAGGTTTTGGGACGCTCGAAGGCTTTAATTTGCAAGCTTCGCAGTTTACACTCTCATCGTCGCTCTCATCATCGCTTCCGTTGTTGTTTTCCTTAGTAGCGTCTGCTTCCAGAGAGTATTTATCTCTTATTACCTCTAAAGGTTCTGCTTGATTTCTGACTTTGTTCGCCTCATGTGCATATTTTTCTTGGTTCTTTTGGGACAAAATATGCGTAAAGGACTTTTGTTGTTCCCTCACATTCCAGTTTAGTTGTCGACTGATACTGTTAATAAACTCATCGGGCGAGGCTTCCACGGTTGGAAAAGCATATGGGCTGGCGCATATGGTTATAAAATCACCTTTTTGCAATTCAATTCTATCTTTCCCATCAAAAGCCGCCCATGCTGGAGCCCTTGACTTCATCGAGACTTTCACTTTTAAATTTATACTTTCTGGTAAGATGATGGGTCTGAAACTCAATGCATGTGGACAAATGGGTGTTAAAGCGATTGCATTGACGGTTGGGCATACCAATGACCCACCTGCACTCAAAGAATAGGCCGTGGACCCAGTCGGAGTAGCAGCAATCAGTCCGTCCGCCTGCGCAACGGTCATTAATGAGCCGTCACCATACAATTCTAACATGGATAGAAAAGGACTTGGACCACGATCGATGGTCACTTCGTTCAAAATGTGGTGTGTGCTTAGTTTTTCCACCACACATATTTTCTTCCCCGTGTTTGGGTCTACTTCAGGGCGGTGTCTACGATAAATTGTG"
    assert eq(correct, candidate, circular=False)

    # text4
    x, y = parse(text4)
    a = assembly.Assembly((x, y), limit=557)
    assert a.assemble_linear()[1].seguid() == "ldseguid=3ay9lXqyTHRr348YDQ4JrtFVx40"
    a = assembly.Assembly((y, x), limit=557)
    assert a.assemble_linear()[1].seguid() == "ldseguid=3ay9lXqyTHRr348YDQ4JrtFVx40"
    a = assembly.Assembly((x.rc(), y), limit=557)
    assert a.assemble_linear()[0].seguid() == "ldseguid=3ay9lXqyTHRr348YDQ4JrtFVx40"
    a = assembly.Assembly((x, y.rc()), limit=557)
    assert a.assemble_linear()[1].seguid() == "ldseguid=3ay9lXqyTHRr348YDQ4JrtFVx40"
    a = assembly.Assembly((x.rc(), y.rc()), limit=557)
    assert a.assemble_linear()[0].seguid() == "ldseguid=3ay9lXqyTHRr348YDQ4JrtFVx40"

    candidate = a.assemble_linear()[0]
    correct = "tcctgacgggtaattttgatttgcatgccgtccgggtgagtcatagcgtctggttgttttgccagattcagcagagtctgtgcaatgcggccgctgaccacatacgatttaggtgacactatagaacgcggccgccagctgaagcttcgtacgctgcaggtcgacggatccccgggttaattaaggcgcgccagatctgtttagcttgccttgtccccgccgggtcacccggccagcgacatggaggcccagaataccctccttgacagtcttgacgtgcgcagctcaggggcatgatgtgactgtcgcccgtacatttagcccatacatccccatgtataatcatttgcatccatacattttgatggccgcacggcgcgaagcaaaaattacggctcctcgctgcagacctgcgagcagggaaacgctcccctcacagacgcgttgaattgtccccacgccgcgcccctgtagagaaatataaaaggttaggatttgccactgaggttcttctttcatatacttccttttaaaatcttgctaggatacagttctcacatcacatccgaacataaacaaccgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgacccaatcacatcacatccgaacataaacaaccatgggtaaaaagcctgaactcaccgcgacgtctgtcgagaagtttctgatcgaaaagttcgacagcgtctccgacctgatgcagctctcggagggcgaagaatctcgtgctttcagcttcgatgtaggagggcgtggatatgtcctgcgggtaaatagctgcgccgatggtttctacaaagatcgttatgtttatcggcactttgcatcggccgcgctcccgattccggaagtgcttgacattggggaattcagcgagagcctgacctattgcatctcccgccgtgcacagggtgtcacgttgcaagacctgcctgaaaccgaactgcccgctgttctgcagccggtcgcggaggccatggatgcgatcgctgcggccgatcttagccagacgagcgggttcggcccattcggaccgcaaggaatcggtcaatacactacatggcgtgatttcatatgcgcgattgctgatccccatgtgtatcactggcaaactgtgatggacgacaccgtcagtgcgtccgtcgcgcaggctctcgatgagctgatgctttgggccgaggactgccccgaagtccggcacctcgtgcacgcggatttcggctccaacaatgtcctgacggacaatggccgcataacagcggtcattgactggagcgaggcgatgttcggggattcccaatacgaggtcgccaacatcttcttctggaggccgtggttggcttgtatggagcagcagacgcgctacttcgagcggaggcatccggagcttgcaggatcgccgcggctccgggcgtatatgctccgcattggtcttgaccaactctatcagagcttggttgacggcaatttcgatgatgcagcttgggcgcagggtcgatgcgacgcaatcgtccgatccggagccgggactgtcgggcgtacacaaatcgcccgcagaagcgcggccgtctggaccgatggctgtgtagaagtactcgccgatagtggaaaccgacgccccagcactcgtccgagggcaaaggaataatcagtactgacaataaaaagattcttgtagggataacagggtaatcggagtgccatctgtgcagacaaacgcatcaggatagagtcttttgtaacgaccccgtctccaccaacttggtatgcttgaaatctcaaggccattacacattcagttatgtgaacgaaaggtctttatttaacgtagcataaactaaataatacaggttccggttagcctgcaatgtgttaaatctaaaggagcatacccaaaatgaactgaagacaaggaaatttgcttgtccagatgtgattgagcatttgaacgttaataacataacatttttatacttaactatagaaagacttgtataaaaactggcaaacgagatattctgaatattggtgcatatttcaggtagaaaagcttacaaaacaatctaatcataatattgagatgaagagaaagataaaagaaaaaacgataagtcagatgagattatgattgtactttgaaatcgaggaacaaagtatatacggtagtagttccccgagttataacgggagatcatgtaaattgagaaaccagataaagatttggtatgcactctagcaagaaaataaaatgatgaatctatgatatagatcacttttgttccagcgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgacgatcagatctttcaggaaagtttcggaggagatagtgttcggcagtttgtacatcatctgcgggatcaggtacggtttgatcaggttgtagaagatcaggtaagacatagaatcgatgtagatgatcggtttgtttttgttgatttttacgtaacagttcagttggaatttgttacgcagacccttaaccaggtattctacttcttcgaaagtgaaagactgggtgttcagtacgatcgatttgttggtagagtttttgttgtaatcccatttaccaccatcatccatgaaccagtatgccagagacatcggggtcaggtagttttcaaccaggttgttcgggatggtttttttgttgttaacgatgaacaggttagccagtttgttgaaagcttggtgtttgaaagtctgggcgccccaggtgattaccaggttacccaggtggttaacacgttcttttttgtgcggcggggacagtacccactgatcgtacagcagacatacgtggtccatgtatgctttgtttttccactcgaactgcatacagtaggttttaccttcatcacgagaacggatgtaagcatcacccaggatcagaccgatacctgcttcgaactgttcgatgttcagttcgatcagctgggatttgtattctttcagcagtttagagttcggacccaggttcattacctggttttttttgatgtttttcatatgcatggatccggggttttttctccttgacgttaaagtatagaggtatattaacaattttttgttgatacttttattacatttgaataagaagtaatacaaaccgaaaatgttgaaagtattagttaaagtggttatgcagtttttgcatttatatatctgttaatagatcaaaaatcatcgcttcgctgattaattaccccagaaataaggctaaaaaactaatcgcattatcatcctatggttgttaatttgattcgttcatttgaaggtttgtggggccaggttactgccaatttttcctcttcataaccataaaagctagtattgtagaatctttattgttcggagcagtgcggcgcgaggcacatctgcgtttcaggaacgcgaccggtgaagacgaggacgcacggaggagagtcttccttcggagggctgtcacccgctcggcggcttctaatccgtacttcaatatagcaatgagcagttaagcgtattactgaaagttccaaagagaaggtttttttaggctaagataatggggctctttacatttccacaacatataagtaagattagatatggatatgtatatggatatgtatatggtggtaatgccatgtaatatgattattaaacttctttgcgtccatccaacgagatctggcgcgccttaattaacccaacctgcattaatgaatcggccaacgcgcggattaccctgttatccctacatattgttgtgccatctgtgcagacaaacgcatcaggattcagtactgacaataaaaagattcttgttttcaagaacttgtcatttgtatagtttttttatattgtagttgttctattttaatcaaatgttagcgtgatttatattttttttcgcctcgacatcatctgcccagatgcgaagttaagtgcgcagaaagtaatatcatgcgtcaatcgtatgtgaatgctggtcgctatactgctgtcgattcgatactaacgccgccatccagtgtcgaaaacgagctcgaattcatcgatgatatcagatccactagtggcctatgcggccgcggatctgccggtctccctatagtgagtcgatccggatttacctgaatcaattggcgaaattttttgtacgaaatttcagccacttcacag"
    assert eq(correct, candidate, circular=False)

    # text5
    list_of_formatted_seq_records = parse(text5)
    a = assembly.Assembly(list_of_formatted_seq_records, limit=60)
    candidate = a.assemble_circular()[-1]

    correct = "tcgcgcgtttcggtgatgacggtgaaaacctctgacacatgcagctcccggagacggtcacagcttgtctgtaagcggatgccgggagcagacaagcccgtcagggcgcgtcagcgggtgttggcgggtgtcggggctggcttaactatgcggcatcagagcagattgtactgagagtgcaccatagatcctgaggatcggggtgataaatcagtctgcgccacatcgggggaaacaaaatggcgcgagatctaaaaaaaaaggctccaaaaggagcctttcgcgctaccaggtaacgcgccactccgacgggattaacgagtgccgtaaacgacgatggttttaccgtgtgcggagatcaggttctgatcctcgagcatcttaagaattcgtcccacggtttgtctagagcagccgacaatctggccaatttcctgacgggtaattttgatttgcatgccgtccgggtgagtcatagcgtctggttgttttgccagattcagcagagtctgtgcaatgcggccgctgaccacatacgatttaggtgacactatagaacgcggccgccagctgaagcttcgtacgctgcaggtcgacggatccccgggttaattaaggcgcgccagatctgtttagcttgccttgtccccgccgggtcacccggccagcgacatggaggcccagaataccctccttgacagtcttgacgtgcgcagctcaggggcatgatgtgactgtcgcccgtacatttagcccatacatccccatgtataatcatttgcatccatacattttgatggccgcacggcgcgaagcaaaaattacggctcctcgctgcagacctgcgagcagggaaacgctcccctcacagacgcgttgaattgtccccacgccgcgcccctgtagagaaatataaaaggttaggatttgccactgaggttcttctttcatatacttccttttaaaatcttgctaggatacagttctcacatcacatccgaacataaacaaccgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgacccaatcacatcacatccgaacataaacaaccatgggtaaaaagcctgaactcaccgcgacgtctgtcgagaagtttctgatcgaaaagttcgacagcgtctccgacctgatgcagctctcggagggcgaagaatctcgtgctttcagcttcgatgtaggagggcgtggatatgtcctgcgggtaaatagctgcgccgatggtttctacaaagatcgttatgtttatcggcactttgcatcggccgcgctcccgattccggaagtgcttgacattggggaattcagcgagagcctgacctattgcatctcccgccgtgcacagggtgtcacgttgcaagacctgcctgaaaccgaactgcccgctgttctgcagccggtcgcggaggccatggatgcgatcgctgcggccgatcttagccagacgagcgggttcggcccattcggaccgcaaggaatcggtcaatacactacatggcgtgatttcatatgcgcgattgctgatccccatgtgtatcactggcaaactgtgatggacgacaccgtcagtgcgtccgtcgcgcaggctctcgatgagctgatgctttgggccgaggactgccccgaagtccggcacctcgtgcacgcggatttcggctccaacaatgtcctgacggacaatggccgcataacagcggtcattgactggagcgaggcgatgttcggggattcccaatacgaggtcgccaacatcttcttctggaggccgtggttggcttgtatggagcagcagacgcgctacttcgagcggaggcatccggagcttgcaggatcgccgcggctccgggcgtatatgctccgcattggtcttgaccaactctatcagagcttggttgacggcaatttcgatgatgcagcttgggcgcagggtcgatgcgacgcaatcgtccgatccggagccgggactgtcgggcgtacacaaatcgcccgcagaagcgcggccgtctggaccgatggctgtgtagaagtactcgccgatagtggaaaccgacgccccagcactcgtccgagggcaaaggaataatcagtactgacaataaaaagattcttgtagggataacagggtaatcggagtgccatctgtgcagacaaacgcatcaggatagagtcttttgtaacgaccccgtctccaccaacttggtatgcttgaaatctcaaggccattacacattcagttatgtgaacgaaaggtctttatttaacgtagcataaactaaataatacaggttccggttagcctgcaatgtgttaaatctaaaggagcatacccaaaatgaactgaagacaaggaaatttgcttgtccagatgtgattgagcatttgaacgttaataacataacatttttatacttaactatagaaagacttgtataaaaactggcaaacgagatattctgaatattggtgcatatttcaggtagaaaagcttacaaaacaatctaatcataatattgagatgaagagaaagataaaagaaaaaacgataagtcagatgagattatgattgtactttgaaatcgaggaacaaagtatatacggtagtagttccccgagttataacgggagatcatgtaaattgagaaaccagataaagatttggtatgcactctagcaagaaaataaaatgatgaatctatgatatagatcacttttgttccagcgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgacgatcagatctttcaggaaagtttcggaggagatagtgttcggcagtttgtacatcatctgcgggatcaggtacggtttgatcaggttgtagaagatcaggtaagacatagaatcgatgtagatgatcggtttgtttttgttgatttttacgtaacagttcagttggaatttgttacgcagacccttaaccaggtattctacttcttcgaaagtgaaagactgggtgttcagtacgatcgatttgttggtagagtttttgttgtaatcccatttaccaccatcatccatgaaccagtatgccagagacatcggggtcaggtagttttcaaccaggttgttcgggatggtttttttgttgttaacgatgaacaggttagccagtttgttgaaagcttggtgtttgaaagtctgggcgccccaggtgattaccaggttacccaggtggttaacacgttcttttttgtgcggcggggacagtacccactgatcgtacagcagacatacgtggtccatgtatgctttgtttttccactcgaactgcatacagtaggttttaccttcatcacgagaacggatgtaagcatcacccaggatcagaccgatacctgcttcgaactgttcgatgttcagttcgatcagctgggatttgtattctttcagcagtttagagttcggacccaggttcattacctggttttttttgatgtttttcatatgcatggatccggggttttttctccttgacgttaaagtatagaggtatattaacaattttttgttgatacttttattacatttgaataagaagtaatacaaaccgaaaatgttgaaagtattagttaaagtggttatgcagtttttgcatttatatatctgttaatagatcaaaaatcatcgcttcgctgattaattaccccagaaataaggctaaaaaactaatcgcattatcatcctatggttgttaatttgattcgttcatttgaaggtttgtggggccaggttactgccaatttttcctcttcataaccataaaagctagtattgtagaatctttattgttcggagcagtgcggcgcgaggcacatctgcgtttcaggaacgcgaccggtgaagacgaggacgcacggaggagagtcttccttcggagggctgtcacccgctcggcggcttctaatccgtacttcaatatagcaatgagcagttaagcgtattactgaaagttccaaagagaaggtttttttaggctaagataatggggctctttacatttccacaacatataagtaagattagatatggatatgtatatggatatgtatatggtggtaatgccatgtaatatgattattaaacttctttgcgtccatccaacgagatctggcgcgccttaattaacccaacctgcattaatgaatcggccaacgcgcggattaccctgttatccctacatattgttgtgccatctgtgcagacaaacgcatcaggattcagtactgacaataaaaagattcttgttttcaagaacttgtcatttgtatagtttttttatattgtagttgttctattttaatcaaatgttagcgtgatttatattttttttcgcctcgacatcatctgcccagatgcgaagttaagtgcgcagaaagtaatatcatgcgtcaatcgtatgtgaatgctggtcgctatactgctgtcgattcgatactaacgccgccatccagtgtcgaaaacgagctcgaattcatcgatgatatcagatccactagtggcctatgcggccgcggatctgccggtctccctatagtgagtcgatccggatttacctgaatcaattggcgaaattttttgtacgaaatttcagccacttcacaggcggttttcgcacgtacccatgcgctacgttcctggccctcttcaaacaggcccagttcgccaataaaatcaccctgattcagataggagaggatcatttctttaccctcttcgtctttgatcagcactgccacagagcctttaacgatgtagtacagcgtttccgctttttcaccctggtgaataagcgtgctcttggatgggtacttatgaatgtggcaatgagacaagaaccattcgagagtaggatccgtttgaggtttaccaagtaccataagatccttaaatttttattatctagctagatgataatattatatcaagaattgtacctgaaagcaaataaattttttatctggcttaactatgcggcatcagagcagattgtactgagagtgcaccatatgcggtgtgaaataccgcacagatgcgtaaggagaaaataccgcatcaggcgctcttccgcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatacggttatccacagaatcaggggataacgcaggaaagaacatgtgagcaaaaggccagcaaaagcccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaaggacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgaggggtaataactgatataattaaattgaagctctaatttgtgagtttagtatacatgcatttacttataatacagttttttagttttgctggccgcatcttctcaaatatgcttcccagcctgcttttctgtaacgttcaccctctaccttagcatcccttccctttgcaaatagtcctcttccaacaataataatgtcagatcctgtagagaccacatcatccacggttctatactgttgacccaatgcgtctcccttgtcatctaaacccacaccgggtgtcataatcaaccaatcgtaaccttcatctcttccacccatgtctctttgagcaataaagccgataacaaaatctttgtcgctcttcgcaatgtcaacagtacccttagtatattctccagtagatagggagcccttgcatgacaattctgctaacatcaaaaggcctctaggttcctttgttacttcttctgccgcctgcttcaaaccgctaacaatacctgggcccaccacaccgtgtgcattcgtaatgtctgcccattctgctattctgtatacacccgcagagtactgcaatttgactgtattaccaatgtcagcaaattttctgtcttcgaagagtaaaaaattgtacttggcggataatgcctttagcggcttaactgtgccctccatggaaaaatcagtcaaaatatccacatgtgtttttagtaaacaaattttgggacctaatgcttcaactaactccagtaattccttggtggtacgaacatccaatgaagcacacaagtttgtttgcttttcgtgcatgatattaaatagcttggcagcaacaggactaggatgagtagcagcacgttccttatatgtagctttcgacatgatttatcttcgtttcctgcaggtttttgttctgtgcagttgggttaagaatactgggcaatttcatgtttcttcaacactacatatgcgtatatataccaatctaagtctgtgctccttccttcgttcttccttctgttcggagattaccgaatcaaaaaaatttcaaagaaaccgaaatcaaaaaaaagaataaaaaaaaaatgatgaattgaattgaaaagctagcttatcgatgataagctgtcaaagatgagaattaattccacggactatagactatactagatactccgtctactgtacgatacacttccgctcaggtccttgtcctttaacgaggccttaccactcttttgttactctattgatccagctcagcaaaggcagtgtgatctaagattctatcttcgcgatgtagtaaaactagctagaccgagaaagagactagaaatgcaaaaggcacttctacaatggctgccatcattattatccgatgtgacgctgcagcttctcaatgatattcgaatacgctttgaggagatacagcctaatatccgacaaactgttttacagatttacgatcgtacttgttacccatcattgaattttgaacatccgaacctgggagttttccctgaaacagatagtatatttgaacctgtataataatatatagtctagcgctttacggaagacaatgtatgtatttcggttcctggagaaactattgcatctattgcataggtaatcttgcacgtcgcatccccggttcattttctgcgtttccatcttgcacttcaatagcatatctttgttaacgaagcatctgtgcttcattttgtagaacaaaaatgcaacgcgagagcgctaatttttcaaacaaagaatctgagctgcatttttacagaacagaaatgcaacgcgaaagcgctattttaccaacgaagaatctgtgcttcatttttgtaaaacaaaaatgcaacgcgacgagagcgctaatttttcaaacaaagaatctgagctgcatttttacagaacagaaatgcaacgcgagagcgctattttaccaacaaagaatctatacttcttttttgttctacaaaaatgcatcccgagagcgctatttttctaacaaagcatcttagattactttttttctcctttgtgcgctctataatgcagtctcttgataactttttgcactgtaggtccgttaaggttagaagaaggctactttggtgtctattttctcttccataaaaaaagcctgactccacttcccgcgtttactgattactagcgaagctgcgggtgcattttttcaagataaaggcatccccgattatattctataccgatgtggattgcgcatactttgtgaacagaaagtgatagcgttgatgattcttcattggtcagaaaattatgaacggtttcttctattttgtctctatatactacgtataggaaatgtttacattttcgtattgttttcgattcactctatgaatagttcttactacaatttttttgtctaaagagtaatactagagataaacataaaaaatgtagaggtcgagtttagatgcaagttcaaggagcgaaaggtggatgggtaggttatatagggatatagcacagagatatatagcaaagagatacttttgagcaatgtttgtggaagcggtattcgcaatgggaagctccaccccggttgataatcagaaaagccccaaaaacaggaagattattatcaaaaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacagttaccaatgcttaatcagtgaggcacctatctcagcgatctgtctatttcgttcatccatagttgcctgactccccgtcgtgtagataactacgatacgggagcgcttaccatctggccccagtgctgcaatgataccgcgagacccacgctcaccggctccagatttatcagcaataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttggcattgctacaggcatcgtggtgtcactctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggataatagtgtatcacatagcagaactttaaaagtgctcatcattggaaaacgttcttcggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatgtaacccactcgtgcacccaactgatcttcagcatcttttactttcaccagcgtttctgggtgagcaaaaacaggaaggcaaaatgccgcaaaaaagggaataagggcgacacggaaatgttgaatactcatactcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaaaagtgccacctgctaagaaaccattattatcatgacattaacctataaaaataggcgtatcacgaggccctttcgtc"

    assert len(correct) == 9772
    assert eq(
        Dseqrecord(correct, circular=True).seguid(),
        candidate.seq.seguid(),
        circular=True,
    )

    # text6

    y, x = parse(text6)
    a = assembly.Assembly((x, y, x), limit=40)
    candidate = a.assemble_linear()[0]
    correct = "aaattgacctcgtttgattgctctcaggtttccaacaacttgtttgttaccggtttcagtacgggtattatcaagttatgggatgccagagcggctgaagctgccaccactgatttgacttatagacaaaacggtgaagatccaattcagaatgaaatcgctaacttctatcatgccggcggtgattctgttgtcgacgtccaattttctgccacttcaagctctgaattctttactgtgggcggcactggtaatatttaccactggaacactgactattccttgtcaaagtacaatcctgacgataccattgctcctcctcaagatgccactgaagaatcacaaacaaaatctctgagattcttgcacaagggtggaagcagaagatctccaaaacagattggaagaagaaacaccgccgcttggcatccagttattgaaaacttggttggtactgtcgatgacgatagtttagtcagcatctacaagccctataccgaggaaagcgaatagtcgccatgctaaacgcgcggaacaaggccatatttatatatttaatgcttttaactattattagtttttctcacccagggctggtttctttaccctgtgtgatgaaagtgcgcgcctaaagttatgtatgagctactgttagtatatttattattcaatatttaattaaacaatacttactactgatgaaagataccgtacacctgctggcgaataagtcactatacacgagcaattatatttctcaggtcattgcgacaattgaaaaatcggatcgactgattatcattaaggcgtgtgacggctcagaataaggaagttctattttaaagggcagcaaaaatttggcacaatttaaaaagaaagttggtctctggtagggttggatccggtttttcataacttgcgtttcatttatggttgcgtatttctcgcagtgtacatagaccaaattaactctatccttccattgcacaatttgcccagtatggatgtccacgaggtctctgctatggacacatttacgcccgtacgctgcaggtcgacggatccccgggttaattaaggcgcgccagatctgtttagcttgcctcgtccccgccgggtcacccggccagcgacatggaggcccagaataccctccttgacagtcttgacgtgcgcagctcaggggcatgatgtgactgtcgcccgtacatttagcccatacatccccatgtataatcatttgcatccatacattttgatggccgcacggcgcgaagcaaaaattacggctcctcgctgcagacctgcgagcagggaaacgctcccctcacagacgcgttgaattgtccccacgccgcgcccctgtagagaaatataaaaggttaggatttgccactgaggttcttctttcatatacttccttttaaaatcttgctaggatacagttctcacatcacatccgaacataaacaaccatgggtaaggaaaagactcacgtttcgaggccgcgattaaattccaacatggatgctgatttatatgggtataaatgggctcgcgataatgtcgggcaatcaggtgcgacaatctatcgattgtatgggaagcccgatgcgccagagttgtttctgaaacatggcaaaggtagcgttgccaatgatgttacagatgagatggtcagactaaactggctgacggaatttatgcctcttccgaccatcaagcattttatccgtactcctgatgatgcatggttactcaccactgcgatccccggcaaaacagcattccaggtattagaagaatatcctgattcaggtgaaaatattgttgatgcgctggcagtgttcctgcgccggttgcattcgattcctgtttgtaattgtccttttaacagcgatcgcgtatttcgtctcgctcaggcgcaatcacgaatgaataacggtttggttgatgcgagtgattttgatgacgagcgtaatggctggcctgttgaacaagtctggaaagaaatgcataagcttttgccattctcaccggattcagtcgtcactcatggtgatttctcacttgataaccttatttttgacgaggggaaattaataggttgtattgatgttggacgagtcggaatcgcagaccgataccaggatcttgccatcctatggaactgcctcggtgagttttctccttcattacagaaacggctttttcaaaaatatggtattgataatcctgatatgaataaattgcagtttcatttgatgctcgatgagtttttctaatcagtactgacaataaaaagattcttgttttcaagaacttgtcatttgtatagtttttttatattgtagttgttctattttaatcaaatgttagcgtgatttatattttttttcgcctcgacatcatctgcccagatgcgaagttaagtgcgcagaaagtaatatcatgcgtcaatcgtatgtgaatgctggtcgctatactgctgtcgattcgatactaacgccgccatccagtgtcgaaaacgagctcgaattcatcgatgaatgacgttatagacgagcctacgagaccgacaccgtaatcaataagcaattttttgcgaacatgattaagtttattttactttatttcgcattttccataaaaaaattttcctccatacaagattcatacccaggatagcctaactaaaatgtatggctttatacagcttgacgacaaacttaatttgaatatatagatatattaatttaatcagcaggattagcatttaatgtcttacatgtctttgtatttggcatacgtattcagtgatattttaatacctcgcaacgcaatttccagcgagtttacgtttgaagcaaacgtaactcccatgagttatacataacgcctgctgcagctgaaccagaacaatataccgtggtattaccgtaagcatttgaaaacaataataagtttccaagtgccgctttagatccaaagtacaaggacatattgtggaagagtaagttgaatatttgaaaatgagagctttttcagcagccaccgttagggccacaactaggaagtcgttcatcccaatggcaccaagaactccttttgtgactccatcatttacaaagaatgtaggctcaatgagaagaatgagattttattctgatgaagccaaaagtgaagaatccaaagaaaacaatgaagatttgactgaagagcaatcagaaatcaagaaattagagagccagttaagcgcgaagactaaagaagcttctgaactcaaggacagattattaagatctgtggcagatttcagaaatttacaacaagtcacaaagaaggatattcagaaagctaaggactttgctttacagaagtttgcaaaggatttattggaatctgtagataactttggtcatgctttgaatgcttttaaagaggaagacttacaaaagtccaaggaaattagtgatttgtatacaggggttagaatgacaagagatgtttttgaaaacaccctaagaaagcacggtattgaaaaattagacccattgggagaaccatttgatccaaataaacacgaa"
    assert len(correct) == 3587
    assert eq(correct, candidate, circular=False)

    # text7 seguid G2drVQAIaRfFXBfqEc5Kddac36A
    list_of_formatted_seq_records = parse(text7)
    a = assembly.Assembly(list_of_formatted_seq_records, limit=28)

    candidate = a.assemble_circular()[1]

    correct = "aattggccagattgtcggctgctctagacaaaccgtgggacgaattcttaagatgctcgaggatcagaacctgatctccgcacacggtaaaaccatcgtcgtttacggcactcgttaatcccgtcggagtggcgcgttacctggtagcgcgaaaggctccttttggagcctttttttttagatctcgcgccattttgtttcccccgatgtggcgcagactgatttatcaccccgatcctcaggatctatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtgatacactattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagagtgacaccacgatgcctgtagcaatgccaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagcgctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataataatcttcctgtttttggggcttttctgattatcaaccggggtggagcttcccattgcgaataccgcttccacaaacattgctcaaaagtatctctttgctatatatctctgtgctatatccctatataacctacccatccacctttcgctccttgaacttgcatctaaactcgacctctacattttttatgtttatctctagtattactctttagacaaaaaaattgtagtaagaactattcatagagtgaatcgaaaacaatacgaaaatgtaaacatttcctatacgtagtatatagagacaaaatagaagaaaccgttcataattttctgaccaatgaagaatcatcaacgctatcactttctgttcacaaagtatgcgcaatccacatcggtatagaatataatcggggatgcctttatcttgaaaaaatgcacccgcagcttcgctagtaatcagtaaacgcgggaagtggagtcaggctttttttatggaagagaaaatagacaccaaagtagccttcttctaaccttaacggacctacagtgcaaaaagttatcaagagactgcattatagagcgcacaaaggagaaaaaaagtaatctaagatgctttgttagaaaaatagcgctctcgggatgcatttttgtagaacaaaaaagaagtatagattctttgttggtaaaatagcgctctcgcgttgcatttctgttctgtaaaaatgcagctcagattctttgtttgaaaaattagcgctctcgtcgcgttgcatttttgttttacaaaaatgaagcacagattcttcgttggtaaaatagcgctttcgcgttgcatttctgttctgtaaaaatgcagctcagattctttgtttgaaaaattagcgctctcgcgttgcatttttgttctacaaaatgaagcacagatgcttcgttaacaaagatatgctattgaagtgcaagatggaaacgcagaaaatgaaccggggatgcgacgtgcaagattacctatgcaatagatgcaatagtttctccaggaaccgaaatacatacattgtcttccgtaaagcgctagactatatattattatacaggttcaaatatactatctgtttcagggaaaactcccaggttcggatgttcaaaattcaatgatgggtaacaagtacgatcgtaaatctgtaaaacagtttgtcggatattaggctgtatctcctcaaagcgtattcgaatatcattgagaagctgcagcgtcacatcggataataatgatggcagccattgtagaagtgccttttgcatttctagtctctttctcggtctagctagttttactacatcgcgaagatagaatcttagatcacactgcctttgctgagctggatcaatagagtaacaaaagagtggtaaggcctcgttaaaggacaaggacctgagcggaagtgtatcgtacagtagacggagtatctagtatagtctatagtccgtggaattaattctcatctttgacagcttatcatcgataagctagcttttcaattcaattcatcattttttttttattcttttttttgatttcggtttctttgaaatttttttgattcggtaatctccgaacagaaggaagaacgaaggaaggagcacagacttagattggtatatatacgcatatgtagtgttgaagaaacatgaaattgcccagtattcttaacccaactgcacagaacaaaaacctgcaggaaacgaagataaatcatgtcgaaagctacatataaggaacgtgctgctactcatcctagtcctgttgctgccaagctatttaatatcatgcacgaaaagcaaacaaacttgtgtgcttcattggatgttcgtaccaccaaggaattactggagttagttgaagcattaggtcccaaaatttgtttactaaaaacacatgtggatattttgactgatttttccatggagggcacagttaagccgctaaaggcattatccgccaagtacaattttttactcttcgaagacagaaaatttgctgacattggtaatacagtcaaattgcagtactctgcgggtgtatacagaatagcagaatgggcagacattacgaatgcacacggtgtggtgggcccaggtattgttagcggtttgaagcaggcggcagaagaagtaacaaaggaacctagaggccttttgatgttagcagaattgtcatgcaagggctccctatctactggagaatatactaagggtactgttgacattgcgaagagcgacaaagattttgttatcggctttattgctcaaagagacatgggtggaagagatgaaggttacgattggttgattatgacacccggtgtgggtttagatgacaagggagacgcattgggtcaacagtatagaaccgtggatgatgtggtctctacaggatctgacattattattgttggaagaggactatttgcaaagggaagggatgctaaggtagagggtgaacgttacagaaaagcaggctgggaagcatatttgagaagatgcggccagcaaaactaaaaaactgtattataagtaaatgcatgtatactaaactcacaaattagagcttcaatttaattatatcagttattacccctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagataaaaaatttatttgctttcaggtacaattcttgatataatattatcatctagctagataataaaaatttaaggatcttatggtacttggtaaacctcaaacggatcctactctcgaatggttcttgtctcattgccacattcataagtacccatccaagagcacgcttattcaccagggtgaaaaagcggaaacgctgtactacatcgttaaaggctctgtggcagtgctgatcaaagacgaagagggtaaagaaatgatcctctcctatctgaatcagggtgattttattggcgaactgggcctgtttgaagagggccaggaacgtagcgcatgggtacgtgcgaaaaccgcctgtgaagtggctgaaatttcgtacaaaaaatttcgccaattgattcaggtaaatccggatgctggaacaaaagtgatctatatcatagattcatcattttattttcttgctagagtgcataccaaatctttatctggtttctcaatttacatgatctcccgttataactcggggaactactaccgtatatactttgttcctcgatttcaaagtacaatcataatctcatctgacttatcgttttttcttttatctttctcttcatctcaatattatgattagattgttttgtaagcttttctacctgaaatatgcaccaatattcagaatatctcgtttgccagtttttatacaagtctttctatagttaagtataaaaatgttatgttattaacgttcaaatgctcaatcacatctggacaagcaaatttccttgtcttcagttcattttgggtatgctcctttagatttaacacattgcaggctaaccggaacctgtattatttagtttatgctacgttaaataaagacctttcgttcacataactgaatgtgtaatggccttgagatttcaagcataccaagttggtggagacggggtcgttacaaaagactctatcctgatgcgtttgtctgcacagatggcgcgtccgattaccctgttatccctacaagaatctttttattgtcagtactgattattcctttgccctcggacgagtgctggggcgtcggtttccactatcggcgagtacttctacacagccatcggtccagacggccgcgcttctgcgggcgatttgtgtacgcccgacagtcccggctccggatcggacgattgcgtcgcatcgaccctgcgcccaagctgcatcatcgaaattgccgtcaaccaagctctgatagagttggtcaagaccaatgcggagcatatacgcccggagccgcggcgatcctgcaagctccggatgcctccgctcgaagtagcgcgtctgctgctccatacaagccaaccacggcctccagaagaagatgttggcgacctcgtattgggaatccccgaacatcgcctcgctccagtcaatgaccgctgttatgcggccattgtccgtcaggacattgttggagccgaaatccgcgtgcacgaggtgccggacttcggggcagtcctcggcccaaagcatcagctcatcgagagcctgcgcgacggacgcactgacggtgtcgtccatcacagtttgccagtgatacacatggggatcagcaatcgcgcatatgaaatcacgccatgtagtgtattgaccgattccttgcggtccgaatgggccgaacccgctcgtctggctaagatcggccgcagcgatcgcatccatggcctccgcgaccggctgcagaacagcgggcagttcggtttcaggcaggtcttgcaacgtgacaccctgtgcacggcgggagatgcaataggtcaggctctcgctgaattccccaatgtcaagcacttccggaatcgggagcgcggccgatgcaaagtgccgataaacataacgatctttgtagaaaccatcggcgcagctatttacccgcaggacatatccacgccctcctacatcgaagctgaaagcacgagattcttcgccctccgagagctgcatcaggtcggagacgctgtcgaacttttcgatcagaaacttctcgacagacgtcgcggtgagttcaggctttttacccatggttgtttatgttcggatgtgatgtgattggccggctgcaggtcactagtgagaaagtgggcaacctggcgttcctcgacatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacgtcagcggccgcattgcacagactctgctgaatctggcaaaacaaccagacgctatgactcacccggacggcatgcaaatcaaaattacccgtcagga"
    assert len(correct) == 7911
    assert eq(
        Dseqrecord(correct, circular=True).seguid(),
        candidate.seq.seguid(),
        circular=True,
    )

    # Contig not implemented
    # assert repr(candidate) == "Contig(o7911)"

    h = """\
 -|630|49
|      \\/
|      /\\
|      49|1196|32
|              \\/
|              /\\
|              32|650|61
|                     \\/
|                     /\\
|                     61|5681_rc|98
|                                \\/
|                                /\\
|                                98-
|                                   |
 -----------------------------------"""  # noqa: F841

    # assert h == candidate.small_figure()

    # assert candidate.detailed_figure() == ("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n"
    #                                        "tcctgacgggtaattttgatttgcatgccgtccgggtgagtcatagcgtctggttgttttgccagattcagcagagtctgtgcaatgcggccgctgacgtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaaggacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgagattatcaaaaaggatcttcacctagatgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgg\n"
    #                                        "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     GTCGAGGAACGCCAGGTTGCCCACTTTCTCACTAGTGACCTGCAGCCGG\n"
    #                                        "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     gtcgaggaacgccaggttgcccactttctcactagtgacctgcagccggccaatcacatcacatccgaacataaacaaccatgggtaaaaagcctgaactcaccgcgacgtctgtcgagaagtttctgatcgaaaagttcgacagcgtctccgacctgatgcagctctcggagggcgaagaatctcgtgctttcagcttcgatgtaggagggcgtggatatgtcctgcgggtaaatagctgcgccgatggtttctacaaagatcgttatgtttatcggcactttgcatcggccgcgctcccgattccggaagtgcttgacattggggaattcagcgagagcctgacctattgcatctcccgccgtgcacagggtgtcacgttgcaagacctgcctgaaaccgaactgcccgctgttctgcagccggtcgcggaggccatggatgcgatcgctgcggccgatcttagccagacgagcgggttcggcccattcggaccgcaaggaatcggtcaatacactacatggcgtgatttcatatgcgcgattgctgatccccatgtgtatcactggcaaactgtgatggacgacaccgtcagtgcgtccgtcgcgcaggctctcgatgagctgatgctttgggccgaggactgccccgaagtccggcacctcgtgcacgcggatttcggctccaacaatgtcctgacggacaatggccgcataacagcggtcattgactggagcgaggcgatgttcggggattcccaatacgaggtcgccaacatcttcttctggaggccgtggttggcttgtatggagcagcagacgcgctacttcgagcggaggcatccggagcttgcaggatcgccgcggctccgggcgtatatgctccgcattggtcttgaccaactctatcagagcttggttgacggcaatttcgatgatgcagcttgggcgcagggtcgatgcgacgcaatcgtccgatccggagccgggactgtcgggcgtacacaaatcgcccgcagaagcgcggccgtctggaccgatggctgtgtagaagtactcgccgatagtggaaaccgacgccccagcactcgtccgagggcaaaggaataatcagtactgacaataaaaagattcttgtagggataacagggtaatcggacgcgccatctgtgcagacaaacgcatcaggatttaaat\n"
    #                                        "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           CGCGCCATCTGTGCAGACAAACGCATCAGGAT\n"
    #                                        "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           cgcgccatctgtgcagacaaacgcatcaggatagagtcttttgtaacgaccccgtctccaccaacttggtatgcttgaaatctcaaggccattacacattcagttatgtgaacgaaaggtctttatttaacgtagcataaactaaataatacaggttccggttagcctgcaatgtgttaaatctaaaggagcatacccaaaatgaactgaagacaaggaaatttgcttgtccagatgtgattgagcatttgaacgttaataacataacatttttatacttaactatagaaagacttgtataaaaactggcaaacgagatattctgaatattggtgcatatttcaggtagaaaagcttacaaaacaatctaatcataatattgagatgaagagaaagataaaagaaaaaacgataagtcagatgagattatgattgtactttgaaatcgaggaacaaagtatatacggtagtagttccccgagttataacgggagatcatgtaaattgagaaaccagataaagatttggtatgcactctagcaagaaaataaaatgatgaatctatgatatagatcacttttgttccagcatccggatttacctgaatcaattggcgaaattttttgtacgaaatttcagccacttcacag\n"
    #                                        "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ATCCGGATTTACCTGAATCAATTGGCGAAATTTTTTGTACGAAATTTCAGCCACTTCACAG\n"
    #                                        "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        atccggatttacctgaatcaattggcgaaattttttgtacgaaatttcagccacttcacaggcggttttcgcacgtacccatgcgctacgttcctggccctcttcaaacaggcccagttcgccaataaaatcaccctgattcagataggagaggatcatttctttaccctcttcgtctttgatcagcactgccacagagcctttaacgatgtagtacagcgtttccgctttttcaccctggtgaataagcgtgctcttggatgggtacttatgaatgtggcaatgagacaagaaccattcgagagtaggatccgtttgaggtttaccaagtaccataagatccttaaatttttattatctagctagatgataatattatatcaagaattgtacctgaaagcaaataaattttttatctggcttaactatgcggcatcagagcagattgtactgagagtgcaccatatgcggtgtgaaataccgcacagatgcgtaaggagaaaataccgcatcaggcgctcttccgcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatacggttatccacagaatcaggggataacgcaggaaagaacatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaaggacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgaggggtaataactgatataattaaattgaagctctaatttgtgagtttagtatacatgcatttacttataatacagttttttagttttgctggccgcatcttctcaaatatgcttcccagcctgcttttctgtaacgttcaccctctaccttagcatcccttccctttgcaaatagtcctcttccaacaataataatgtcagatcctgtagagaccacatcatccacggttctatactgttgacccaatgcgtctcccttgtcatctaaacccacaccgggtgtcataatcaaccaatcgtaaccttcatctcttccacccatgtctctttgagcaataaagccgataacaaaatctttgtcgctcttcgcaatgtcaacagtacccttagtatattctccagtagatagggagcccttgcatgacaattctgctaacatcaaaaggcctctaggttcctttgttacttcttctgccgcctgcttcaaaccgctaacaatacctgggcccaccacaccgtgtgcattcgtaatgtctgcccattctgctattctgtatacacccgcagagtactgcaatttgactgtattaccaatgtcagcaaattttctgtcttcgaagagtaaaaaattgtacttggcggataatgcctttagcggcttaactgtgccctccatggaaaaatcagtcaaaatatccacatgtgtttttagtaaacaaattttgggacctaatgcttcaactaactccagtaattccttggtggtacgaacatccaatgaagcacacaagtttgtttgcttttcgtgcatgatattaaatagcttggcagcaacaggactaggatgagtagcagcacgttccttatatgtagctttcgacatgatttatcttcgtttcctgcaggtttttgttctgtgcagttgggttaagaatactgggcaatttcatgtttcttcaacactacatatgcgtatatataccaatctaagtctgtgctccttccttcgttcttccttctgttcggagattaccgaatcaaaaaaatttcaaagaaaccgaaatcaaaaaaaagaataaaaaaaaaatgatgaattgaattgaaaagctagcttatcgatgataagctgtcaaagatgagaattaattccacggactatagactatactagatactccgtctactgtacgatacacttccgctcaggtccttgtcctttaacgaggccttaccactcttttgttactctattgatccagctcagcaaaggcagtgtgatctaagattctatcttcgcgatgtagtaaaactagctagaccgagaaagagactagaaatgcaaaaggcacttctacaatggctgccatcattattatccgatgtgacgctgcagcttctcaatgatattcgaatacgctttgaggagatacagcctaatatccgacaaactgttttacagatttacgatcgtacttgttacccatcattgaattttgaacatccgaacctgggagttttccctgaaacagatagtatatttgaacctgtataataatatatagtctagcgctttacggaagacaatgtatgtatttcggttcctggagaaactattgcatctattgcataggtaatcttgcacgtcgcatccccggttcattttctgcgtttccatcttgcacttcaatagcatatctttgttaacgaagcatctgtgcttcattttgtagaacaaaaatgcaacgcgagagcgctaatttttcaaacaaagaatctgagctgcatttttacagaacagaaatgcaacgcgaaagcgctattttaccaacgaagaatctgtgcttcatttttgtaaaacaaaaatgcaacgcgacgagagcgctaatttttcaaacaaagaatctgagctgcatttttacagaacagaaatgcaacgcgagagcgctattttaccaacaaagaatctatacttcttttttgttctacaaaaatgcatcccgagagcgctatttttctaacaaagcatcttagattactttttttctcctttgtgcgctctataatgcagtctcttgataactttttgcactgtaggtccgttaaggttagaagaaggctactttggtgtctattttctcttccataaaaaaagcctgactccacttcccgcgtttactgattactagcgaagctgcgggtgcattttttcaagataaaggcatccccgattatattctataccgatgtggattgcgcatactttgtgaacagaaagtgatagcgttgatgattcttcattggtcagaaaattatgaacggtttcttctattttgtctctatatactacgtataggaaatgtttacattttcgtattgttttcgattcactctatgaatagttcttactacaatttttttgtctaaagagtaatactagagataaacataaaaaatgtagaggtcgagtttagatgcaagttcaaggagcgaaaggtggatgggtaggttatatagggatatagcacagagatatatagcaaagagatacttttgagcaatgtttgtggaagcggtattcgcaatgggaagctccaccccggttgataatcagaaaagccccaaaaacaggaagattattatcaaaaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacagttaccaatgcttaatcagtgaggcacctatctcagcgatctgtctatttcgttcatccatagttgcctgactccccgtcgtgtagataactacgatacgggagcgcttaccatctggccccagtgctgcaatgataccgcgagacccacgctcaccggctccagatttatcagcaataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttggcattgctacaggcatcgtggtgtcactctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggataatagtgtatcacatagcagaactttaaaagtgctcatcattggaaaacgttcttcggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatgtaacccactcgtgcacccaactgatcttcagcatcttttactttcaccagcgtttctgggtgagcaaaaacaggaaggcaaaatgccgcaaaaaagggaataagggcgacacggaaatgttgaatactcatactcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaaaagtgccacctgctaagaaaccattattatcatgacattaacctataaaaataggcgtatcacgaggccctttcgtctcgcgcgtttcggtgatgacggtgaaaacctctgacacatgcagctcccggagacggtcacagcttgtctgtaagcggatgccgggagcagacaagcccgtcagggcgcgtcagcgggtgttggcgggtgtcggggctggcttaactatgcggcatcagagcagattgtactgagagtgcaccatagatcctgaggatcggggtgataaatcagtctgcgccacatcgggggaaacaaaatggcgcgagatctaaaaaaaaaggctccaaaaggagcctttcgcgctaccaggtaacgcgccactccgacgggattaacgagtgccgtaaacgacgatggttttaccgtgtgcggagatcaggttctgatcctcgagcatcttaagaattcgtcccacggtttgtctagagcagccgacaatctggccaatttcctgacgggtaattttgatttgcatgccgtccgggtgagtcatagcgtctggttgttttgccagattcagcagagtctgtgcaatgcggccgctgac\n"
    #                                        "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       TCCTGACGGGTAATTTTGATTTGCATGCCGTCCGGGTGAGTCATAGCGTCTGGTTGTTTTGCCAGATTCAGCAGAGTCTGTGCAATGCGGCCGCTGAC\n")

    # text8 seguid wM7nM6oJer3bB6RV81IH78e02j4 7911bp
    list_of_formatted_seq_records = parse(text8)
    a = assembly.Assembly(list_of_formatted_seq_records, limit=28)
    candidate = a.assemble_circular()[0]
    assert candidate.seguid() == "cdseguid=-mVwekticpAYIT9C4JcXmOGFkRo"


def test_MXblaster1():
    """test MXblaster1"""

    primer = parse(f"{test_files}/primers.fas", ds=False)
    primer = primer[::-1]
    primer = primer[37:]

    for i, p in enumerate(primer):
        assert int(p.id.split("_")[0]) == i

    """ These are PCRs to get the genes and the terminator-promoters """
    AgTEFp = pcr(primer[524], primer[523], read(f"{test_files}/pAG25.gb"))
    hph = pcr(primer[502], primer[501], read(f"{test_files}/pAG32.gb"))
    KlLEU2tt = pcr(primer[520], primer[519], read(f"{test_files}/KlLEU2tt.gb"))

    """ The Gal1 promoter-ISceI fragment is made in two steps """
    gal1_ISceI_1 = pcr(primer[234], primer[316], read(f"{test_files}/pGSHU 7180bp.gb"))
    gal1_ISceI_2 = pcr(primer[562], primer[234], gal1_ISceI_1)
    AgTEFt = pcr(primer[522], primer[521], read(f"{test_files}/pAG25.gb"))

    """ load pCAPs and pCAPs-pSU0 sequences as Dseqrecord objects """
    pCAPs = read(f"{test_files}/pCAPs-AjiI.gb")
    pCAPs_pSU0 = read(f"{test_files}/pCAPs-pSU0.gb")

    # cut the pCAPs vectors for cloning

    pCAPs_ZraI = pCAPs.linearize(ZraI)
    pCAPs_PCR_prod = pcr(primer[492], primer[493], pCAPs)
    pCAPs_EcoRV = pCAPs.linearize(EcoRV)
    stuffer, pCAPs_pSU0_E_Z = pCAPs_pSU0.cut(EcoRV, ZraI)

    # make the pCAPs clones, six altoghether
    pCAPs_ZraI_AgTEFp = (pCAPs_ZraI + AgTEFp).looped()

    pCAPs_PCR_prod_hph = (pCAPs_PCR_prod + hph).looped()
    pCAPs_EcoRV_KlLEU2tt = (pCAPs_EcoRV + KlLEU2tt).looped()

    pCAPs_ZraI_KlLEU2tt = (pCAPs_ZraI + KlLEU2tt).looped()
    pCAPs_PCR_prod_gal1_ISceI_2 = (pCAPs_PCR_prod + gal1_ISceI_2).looped()
    pCAPs_EcoRV_AgTEFt = (pCAPs_EcoRV + AgTEFt).looped()

    # PCR each clone for the assembly in pCAPs

    A_AgTEFp_b = pcr([primer[167], primer[493]], pCAPs_ZraI_AgTEFp)
    B_hph_c = pcr([primer[467], primer[468]], pCAPs_PCR_prod_hph)
    C_KlLEU2tt_d = pcr([primer[492], primer[166]], pCAPs_EcoRV_KlLEU2tt)

    # Homologous recombination of the two tp-gene-tp building blocks

    a = assembly.Assembly((pCAPs_pSU0_E_Z, A_AgTEFp_b, B_hph_c, C_KlLEU2tt_d), limit=28)
    candidate = a.assemble_circular()[0]
    assert candidate.seguid() == "cdseguid=-mVwekticpAYIT9C4JcXmOGFkRo"
    assert len(candidate) == 7911
    YPK0_AgTEFp_hph_KlLEU2tt = candidate

    x = YPK0_AgTEFp_hph_KlLEU2tt

    AgTEFp_hph_KlLEU2tt_2 = pcr(primer[166], primer[167], YPK0_AgTEFp_hph_KlLEU2tt)

    A_KlLEU2tt_b = pcr([primer[167], primer[567]], pCAPs_ZraI_KlLEU2tt)
    B_gal1_ISceI_c = pcr([primer[467], primer[468]], pCAPs_PCR_prod_gal1_ISceI_2)
    C_AgTEFt_d = pcr([primer[568], primer[166]], pCAPs_EcoRV_AgTEFt)

    a = assembly.Assembly(
        (pCAPs_pSU0_E_Z, A_KlLEU2tt_b, B_gal1_ISceI_c, C_AgTEFt_d), limit=25
    )
    candidate = a.assemble_circular()[0]
    assert candidate.seguid() == "cdseguid=GGV2uIPTRD0lcV7RJPbVTItXNII"
    assert len(candidate) == 8099
    YPK0_KlLEU2tt_gal1_ISceI_AgTEFt = candidate

    feats = {}

    for f in YPK0_AgTEFp_hph_KlLEU2tt.features:
        feats[f.qualifiers["label"][0]] = f.extract(YPK0_AgTEFp_hph_KlLEU2tt).seq

    oldfeats = {}

    for x in (pCAPs_pSU0_E_Z, A_AgTEFp_b, B_hph_c, C_KlLEU2tt_d):
        for f in x.features:
            oldfeats[f.qualifiers["label"][0]] = f.extract(x).seq

    KlLEU2tt_gal1_ISceI_AgTEFt_2 = pcr(
        primer[166], primer[167], YPK0_KlLEU2tt_gal1_ISceI_AgTEFt
    )

    a = assembly.Assembly(
        (AgTEFp_hph_KlLEU2tt_2, KlLEU2tt_gal1_ISceI_AgTEFt_2, pCAPs_pSU0_E_Z), limit=61
    )
    candidate = a.assemble_circular()[0]
    assert len(candidate) == 9772
    assert candidate.seguid() == "cdseguid=bUl04KTp5LpAulZX3UHdejwnuIQ"
    pCAPs_MX4blaster1 = candidate

    pCAPs_MX4blaster1 = pCAPs_MX4blaster1.synced("tcgcgcgtttcggtgatgacggtgaaaacc")

    assert pCAPs_MX4blaster1.seguid() == "cdseguid=bUl04KTp5LpAulZX3UHdejwnuIQ"

    AX023560 = read(f"{test_files}/AX023560.gb")

    GAL10prom_slice = slice(
        AX023560.features[1].location.start, AX023560.features[1].location.end
    )

    GAL10prom = AX023560[GAL10prom_slice]

    assert GAL10prom.seq == AX023560.features[1].extract(AX023560).seq

    GIN11M86 = read(f"{test_files}/GIN11M86.gb")

    GAL_GIN = pcr(primer[592], primer[593], GAL10prom + GIN11M86)

    assert GAL_GIN.seguid() == "ldseguid=E4I5zDE-86ifdL8J9Wd9K827oLo"

    assert pCAPs.seguid() == "cdseguid=IFLMpKpCGZio0R0YkGSfqPaKKiw"

    pCAPs_GAL_GIN = (pCAPs.cut(AjiI)[0] + GAL_GIN).looped()

    assert pCAPs_GAL_GIN.seguid() == "cdseguid=WLRtLIqKRV3iM_le_IL7kDUqN2I"

    GAL_GIN2 = pcr(primer[592], primer[467], pCAPs_GAL_GIN)

    assert GAL_GIN2.seguid() == "ldseguid=VMbnoWQyowa92XOf2wbKlsM26f8"

    assert (
        pCAPs_MX4blaster1.seguid() == "cdseguid=bUl04KTp5LpAulZX3UHdejwnuIQ"
    )  # 9772bp__a

    pCAPs_MX4blaster1_AgeI = pCAPs_MX4blaster1.cut(AgeI)[0]

    pCAPs_MX4blaster1_AgeI.seq = pCAPs_MX4blaster1_AgeI.seq.fill_in()

    a = assembly.Assembly([GAL_GIN2, pCAPs_MX4blaster1_AgeI], limit=30)

    # Changed
    candidates = a.assemble_circular()
    candidate = candidates[2]

    assert len(candidate) == 10566

    assert candidate.seguid() == "cdseguid=c48cBUb3wF-Sdhzh0Tlprp-0CEg"

    pCAPs_MX4blaster2 = candidate

    pCAPs_MX4blaster2 = pCAPs_MX4blaster2.synced("tcgcgcgtttcggtgatgacggtgaaaacc")

    assert len(pCAPs_MX4blaster2) == 10566
    pCAPs_MX4blaster2_old = read(f"{test_files}/pMX4blaster2_old.gb")

    assert len(pCAPs_MX4blaster2_old) == 10566
    assert pCAPs_MX4blaster2_old.seguid() == "cdseguid=c48cBUb3wF-Sdhzh0Tlprp-0CEg"
    assert eq(pCAPs_MX4blaster2, pCAPs_MX4blaster2_old)
    assert pCAPs_MX4blaster2.seguid() == "cdseguid=c48cBUb3wF-Sdhzh0Tlprp-0CEg"


def test_assemble_pGUP1():

    GUP1rec1sens = read(f"{test_files}/GUP1rec1sens.txt")
    GUP1rec2AS = read(f"{test_files}/GUP1rec2AS.txt")
    GUP1_locus = read(f"{test_files}/GUP1_locus.gb")
    pGREG505 = read(f"{test_files}/pGREG505.gb")

    insert = pcr(GUP1rec1sens, GUP1rec2AS, GUP1_locus)

    his3, lin_vect = pGREG505.cut(SalI)

    ab = assembly.Assembly([insert, lin_vect], limit=28)

    pGUP1 = ab.assemble_circular()[0]

    pGUP1 = pGUP1.synced(pGREG505.seq[:50])

    pGUP1_correct = read(f"{test_files}/pGUP1_correct.gb")

    assert len(pGUP1_correct) == 9981
    assert len(pGUP1) == 9981
    assert eq(pGUP1, pGUP1_correct)
    assert pGUP1_correct.seguid() == "cdseguid=QiK2pH9yioTPfSobUTLz4CPiNzY"
    assert pGUP1.seguid() == "cdseguid=QiK2pH9yioTPfSobUTLz4CPiNzY"


def test_pYPK7_TDH3_GAL2_PGI1():

    pMEC1142 = read(f"{test_files}/pYPK0_TDH3_GAL2_PGI1.gb")

    pYPKp7 = read(f"{test_files}/pYPKp7.gb")

    pYPKp7_AatII = pYPKp7.linearize(AatII)

    asm = assembly.Assembly((pYPKp7_AatII, pMEC1142), limit=300)

    (result,) = asm.assemble_circular()

    pYPK7_TDH3_GAL2_PGI1 = read(f"{test_files}/pYPK7_TDH3_GAL2_PGI1.gb")

    assert (
        result.seguid()
        == pYPK7_TDH3_GAL2_PGI1.seguid()
        == "cdseguid=DeflrptvvS6m532WogvxQSgVKpk"
    )

    assert len(result) == len(pYPK7_TDH3_GAL2_PGI1) == 9780


def test_marker_replacement_on_plasmid():

    f, r, _, _ = parse(
        """

    >807_pYPK0_hygfwd2
    tctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgagCACATACGATTTAGGTGACAC

    >806_pYPK0_hygrev2
    atagtccgtggaattaattctcatctttgacagcttatcatcgataagctCGACTCACTATAGGGAGAC

    >678_pYPK0_hygfwd: (77-mer)
    ctcacgttaagggattttggtcatgagCACATACGATTTAGGTGACACTATAGAAC

    >666_pYPK0_hygrev (70-mer)
    catctttgacagcttatcatcgataagctCGACTCACTATAGGGAGACC
    """
    )

    pAG32 = read(f"{test_files}/pAG32.gb")
    pMEC1135 = read(f"{test_files}/pMEC1135.gb")

    hygromycin_product = pcr(f, r, pAG32)
    # This is an homologous recombination, so constrains should be applied
    asm_hyg = assembly.Assembly((hygromycin_product, pMEC1135), limit=50)
    candidate = asm_hyg.assemble_circular()[0]

    plasmid_feat = next(
        f
        for f in pMEC1135.features
        if "label" in f.qualifiers and f.qualifiers["label"][0] == "AmpR"
    )
    candidate_feat = next(
        f
        for f in candidate.features
        if "label" in f.qualifiers and f.qualifiers["label"][0] == "AmpR"
    )
    assert plasmid_feat.extract(pMEC1135).seq == candidate_feat.extract(candidate).seq


@pytest.mark.xfail(reason="contig not implemented")
def test_linear_with_annotations2():
    # Thanks to James Bagley for finding this bug
    # https://github.com/JamesBagley

    a = Dseqrecord("acgatgctatactgtgCCNCCtgtgctgtgctcta")
    a.add_feature(0, 10, label="a_feat")
    a_feat_seq = a.features[0].extract(a)
    # 12345678901234
    b = Dseqrecord("tgtgctgtgctctaTTTTTTTtattctggctgtatcCCCCCC")
    b.add_feature(0, 10, label="b_feat")
    b_feat_seq = b.features[0].extract(b)

    # 123456789012345
    c = Dseqrecord("GtattctggctgtatcGGGGGtacgatgctatactgtg")
    c.add_feature(0, 10, label="c_feat")
    c_feat_seq = c.features[0].extract(c)

    feature_sequences = {
        "a_feat": a_feat_seq,
        "b_feat": b_feat_seq,
        "c_feat": c_feat_seq,
    }

    a.name = "aaa"  # 1234567890123456
    b.name = "bbb"
    c.name = "ccc"
    asm = assembly.Assembly((a, b, c), limit=14)
    x = asm.assemble_linear()[0]
    answer = "aaa|14\n    \\/\n    /\\\n    14|bbb|15\n           \\/\n           /\\\n           15|ccc"

    assert x.figure() == answer.strip()
    answer = "acgatgctatactgtgCCNCCtgtgctgtgctcta\n                     TGTGCTGTGCTCTA\n                     tgtgctgtgctctaTTTTTTTtattctggctgtatc\n                                          TATTCTGGCTGTATC\n                                          tattctggctgtatcGGGGGtacgatgctatactgtg\n"
    assert x.detailed_figure()
    for feat in x.features:
        try:
            assert (
                feat.extract(x).seq == feature_sequences[feat.qualifiers["label"]].seq
            )
        except AssertionError:
            assert (
                feat.extract(x).seq == feature_sequences[feat.qualifiers["label"]].seq
            )


def test_sticky_ligation_algorithm_and_assembly():

    # Test full overlap
    seqrA = Dseqrecord(Dseq.from_full_sequence_and_overhangs("AAAGAT", 0, 3))
    seqrB = Dseqrecord(Dseq.from_full_sequence_and_overhangs("GATAAA", 3, 0))

    # AAAGAT   AAA
    # TTT   CTATTT
    assert assembly.sticky_end_sub_strings(seqrA, seqrB, 0) == [(3, 0, 3)]
    asm = assembly.Assembly(
        [seqrA, seqrB], limit=False, algorithm=assembly.sticky_end_sub_strings
    )
    assert str(asm.assemble_linear()[0].seq) == "AAAGATAAA"

    # Test partial overlap
    # AAAGAT AAA
    # TTT  TATTT
    seqrB = Dseqrecord(Dseq.from_full_sequence_and_overhangs("ATAAA", 2, 0))
    assert assembly.sticky_end_sub_strings(seqrA, seqrB, 1) == [(4, 0, 2)]

    # We test partial overlap settings here
    asm = assembly.Assembly(
        [seqrA, seqrB], limit=False, algorithm=assembly.sticky_end_sub_strings
    )
    assert asm.assemble_linear() == []
    asm = assembly.Assembly(
        [seqrA, seqrB], limit=True, algorithm=assembly.sticky_end_sub_strings
    )
    assert str(asm.assemble_linear()[0].seq) == "AAAGATAAA"

    # Test no overlap
    # AAAGAT CCC
    # TTT  GGGGG
    seqrB = Dseqrecord(Dseq.from_full_sequence_and_overhangs("CCCCC", 2, 0))
    assert assembly.sticky_end_sub_strings(seqrA, seqrB, 1) == []

    # Same examples, with opposite overhangs
    seqrA = Dseqrecord(Dseq.from_full_sequence_and_overhangs("AAAGAT", 0, -3))
    seqrB = Dseqrecord(Dseq.from_full_sequence_and_overhangs("GATAAA", -3, 0))
    assert assembly.sticky_end_sub_strings(seqrA, seqrB, 0) == [(3, 0, 3)]
    asm = assembly.Assembly(
        [seqrA, seqrB], limit=False, algorithm=assembly.sticky_end_sub_strings
    )
    assert str(asm.assemble_linear()[0].seq) == "AAAGATAAA"

    seqrB = Dseqrecord(Dseq.from_full_sequence_and_overhangs("ATAAA", -2, 0))
    assert assembly.sticky_end_sub_strings(seqrA, seqrB, 1) == [(4, 0, 2)]
    asm = assembly.Assembly(
        [seqrA, seqrB], limit=False, algorithm=assembly.sticky_end_sub_strings
    )
    assert asm.assemble_linear() == []
    asm = assembly.Assembly(
        [seqrA, seqrB], limit=True, algorithm=assembly.sticky_end_sub_strings
    )
    assert str(asm.assemble_linear()[0].seq) == "AAAGATAAA"

    seqrB = Dseqrecord(Dseq.from_full_sequence_and_overhangs("CCCCC", -2, 0))
    assert assembly.sticky_end_sub_strings(seqrA, seqrB, 1) == []


def test_restriction_ligation_algorithm():

    # Full overlap, negative ovhg
    seqrA = Dseqrecord("AAAGAATTCAAA")
    seqrB = Dseqrecord("CCCCGAATTCCCCGAATTC")
    assert assembly.restriction_ligation_overlap(seqrA, seqrB, [EcoRI], False) == [
        (4, 5, 4),
        (4, 14, 4),
    ]
    assert assembly.restriction_ligation_overlap(seqrB, seqrA, [EcoRI], False) == [
        (5, 4, 4),
        (14, 4, 4),
    ]

    # Full overlap, positive ovhg
    seqrA = Dseqrecord("TTGCGATCGCTT")
    seqrB = Dseqrecord("AAGCGATCGCAAGCGATCGCAA")
    assert assembly.restriction_ligation_overlap(seqrA, seqrB, [RgaI], False) == [
        (5, 5, 2),
        (5, 15, 2),
    ]
    assert assembly.restriction_ligation_overlap(seqrB, seqrA, [RgaI], False) == [
        (5, 5, 2),
        (15, 5, 2),
    ]

    # Full overlap, using two different enzymes
    seqrA = Dseqrecord("GACTAATGGGTC")
    seqrB = Dseqrecord("AAGCGATCGCAAGCGATCGCAA")
    assert assembly.restriction_ligation_overlap(seqrA, seqrB, [RgaI, DrdI], False) == [
        (5, 5, 2),
        (5, 15, 2),
    ]
    assert assembly.restriction_ligation_overlap(seqrB, seqrA, [RgaI, DrdI], False) == [
        (5, 5, 2),
        (15, 5, 2),
    ]

    # Partial overlap, positive ovhg
    seqrA = Dseqrecord("GACTAAAGGGTC")
    seqrB = Dseqrecord("AAGCGATCGCAAGCGATCGCAA")
    assert assembly.restriction_ligation_overlap(seqrA, seqrB, [RgaI, DrdI], True) == [
        (6, 5, 1),
        (6, 15, 1),
    ]
    assert assembly.restriction_ligation_overlap(seqrB, seqrA, [RgaI, DrdI], True) == []

    # Partial overlap, negative ovhg
    seqrA = Dseqrecord("CCCGCAAAAAAAAA")
    seqrB = Dseqrecord("AAGCTT")
    assert assembly.restriction_ligation_overlap(
        seqrA, seqrB, [HindIII, FauI], True
    ) == [(10, 1, 1)]
    assert (
        assembly.restriction_ligation_overlap(seqrB, seqrA, [HindIII, FauI], True) == []
    )

    # Circular molecule
    seqrA = Dseqrecord("GAATTCaaa")
    seqrB = Dseqrecord("CCCCGAATTCCCC", circular=True)
    for shift in range(len(seqrB)):
        index_in_circular = 5 - shift
        if index_in_circular < 0:
            index_in_circular = len(seqrB) + index_in_circular
        assert assembly.restriction_ligation_overlap(
            seqrA, seqrB.shifted(shift), [EcoRI], False
        ) == [(1, index_in_circular, 4)]


def test_fill_dseq():

    solution = Dseq("ACGT")
    for query in [
        Dseq("ACGT", "T", 0),
        Dseq("ACGT", "G", -1),
        Dseq("ACGT", "C", -2),
        Dseq("ACGT", "A", -3),
    ]:
        assert assembly.fill_dseq(query) == solution


def test_pcr_assembly_normal():

    primer1 = Primer("ACGTACGT")
    primer2 = Primer(reverse_complement("GCGCGCGC"))

    seq = Dseqrecord("ccccACGTACGTAAAAAAGCGCGCGCcccc")

    asm = assembly.PCRAssembly([primer1, seq, primer2], limit=8)

    prods = asm.assemble_linear()

    assert len(prods) == 1
    assert str(prods[0].seq) == "ACGTACGTAAAAAAGCGCGCGC"

    # Try to pass the reverse complement
    asm = assembly.PCRAssembly([primer1, seq.reverse_complement(), primer2], limit=8)
    prods = asm.assemble_linear()

    assert len(prods) == 1
    assert str(prods[0].seq) == "ACGTACGTAAAAAAGCGCGCGC"

    # When the product exactly matches the template
    asm = assembly.PCRAssembly(
        [primer1, Dseqrecord(Dseq("ACGTACGTAAAAAAGCGCGCGC")), primer2], limit=8
    )
    prods = asm.assemble_linear()

    assert len(prods) == 1
    assert str(prods[0].seq) == "ACGTACGTAAAAAAGCGCGCGC"

    # Primers with overhangs work
    primer1 = Primer("TTTACGTACGT")
    primer2 = Primer(reverse_complement("GCGCGCGCTTT"))

    asm = assembly.PCRAssembly([primer1, seq, primer2], limit=8)
    prods = asm.assemble_linear()

    assert len(prods) == 1
    assert str(prods[0].seq) == "TTTACGTACGTAAAAAAGCGCGCGCTTT"


def test_pcr_assembly_uracil():

    primer1 = Primer("AUUA")
    primer2 = Primer("UUAA")

    seq = Dseqrecord(Dseq("aaATTAggccggTTAAaa"))
    asm = assembly.PCRAssembly([primer1, seq, primer2], limit=4)

    assert (
        str(asm.assemble_linear()[0].seq) == "AUUAggccggTTOO"
    )  # FIXME: This is the expected result, see def in utils
    assert asm.assemble_linear()[0].seq.crick.startswith("UUAA")

    primer1 = Primer("ATAUUA")
    primer2 = Primer("ATUUAA")
    asm = assembly.PCRAssembly([primer1, seq, primer2], limit=6, mismatches=1)
    assert str(asm.assemble_linear()[0].seq) == "ATAUUAggccggTTOOAT"


def test_pcr_with_mistmaches():
    primer1 = Primer("atcttGagacgtgtattt")
    primer2 = Primer(reverse_complement("ccaagtgcctccGtttta"))

    seq = Dseqrecord("atcatcagacgtgtatttcacaagccagaagtgcatttggatccaagtgcctccatttta")

    # Works for 1 and 2 mismatches
    for i in [1, 2]:
        asm = assembly.PCRAssembly([primer1, seq, primer2], limit=14, mismatches=i)
        prods = asm.assemble_linear()

        assert len(prods) == 1
        assert (
            str(prods[0].seq)
            == "atcttGagacgtgtatttcacaagccagaagtgcatttggatccaagtgcctccGtttta"
        )

    asm = assembly.PCRAssembly([primer1, seq, primer2], limit=14, mismatches=0)
    prods = asm.assemble_linear()
    assert len(prods) == 0


# @pytest.mark.xfail(reason="This should not work.")
def test_pcrs_with_overlapping_primers_circular_templates():

    # ACGTTCGTGCGTTTTGC
    # |||||||||||||||||
    #       CACGCAAAACG
    #
    # ACGTTCGTGC
    # |||||||||||||||||
    # TGCAAGCACGCAAAACG

    seq = Dseqrecord(Dseq("ACGTTCGTGCGTTTTGC", circular=True))

    # Overlapping 5', edge case for circular extract_subfragment
    primer1 = Primer("ACGTTCGTGC")
    primer2 = Primer(reverse_complement("GTGCGTTTTGC"))  # CACGCAAAACG

    for shift in range(len(seq)):
        seq_shifted = seq.shifted(shift)
        asm = assembly.PCRAssembly([primer1, seq_shifted, primer2], limit=8)
        asm.assemble_linear()
        assert str(asm.assemble_linear()[0].seq) == "ACGTTCGTGCGTTTTGC"

    # Overlapping 5' and 3'
    primer1 = Primer("GCACGTTCGTG")
    for shift in range(len(seq)):
        seq_shifted = seq.shifted(shift)
        asm = assembly.PCRAssembly([primer1, seq_shifted, primer2], limit=8)
        assert str(asm.assemble_linear()[0].seq) == "GCACGTTCGTGCGTTTTGC"


# @pytest.mark.xfail(reason="This should not work.")
def test_pcrs_with_overlapping_primers_linear_templates():

    # ACGTTCGTGCGTTTTGC
    # |||||||||||||||||
    #       CACGCAAAACG
    #
    # ACGTTCGTGC
    # |||||||||||||||||
    # TGCAAGCACGCAAAACG

    seq = Dseqrecord(Dseq("ACGTTCGTGCGTTTTGC", circular=False))

    # Overlapping 5', do as normal
    primer1 = Primer("ACGTTCGTGC")
    primer2 = Primer(reverse_complement("GTGCGTTTTGC"))  # CACGCAAAACG

    asm = assembly.PCRAssembly([primer1, seq, primer2], limit=8)
    assert str(asm.assemble_linear()[0].seq) == "ACGTTCGTGCGTTTTGC"

    # Overlapping 3', should not work
    primer1 = Primer("GTGCGTTTTGC")
    primer2 = Primer(reverse_complement("ACGTTCGTGCG"))

    asm = assembly.PCRAssembly([primer1, seq, primer2], limit=8)
    assert len(asm.assemble_linear()) == 0


def test_pcr_assembly_invalid():

    primer1 = Primer("ACGTACGT")
    primer2 = Primer(reverse_complement("GCGCGCGC"))

    seq = Dseqrecord(Dseq("TTTTACGTACGTAAAAAAGCGCGCGCTTTTT"))

    # Limit too high
    asm = assembly.PCRAssembly([primer1, seq, primer2], limit=15)
    prods = asm.assemble_linear()

    assert len(prods) == 0

    # Clashing primers

    # ACGTTCGTGCGTTTTGC
    # |||||||||||||||||
    #       CACGCAAAACG
    #
    # ACGTTCGTGC
    # |||||||||||||||||
    # TGCAAGCACGCAAAACG

    seq = Dseqrecord(Dseq("ACGTTCGTGCGTTTTGC"))
    primer1 = Primer("ACGTTCGTGC")
    primer2 = Primer(reverse_complement("GTGCGTTTTGC"))

    asm = assembly.PCRAssembly([primer1, seq, primer2], limit=8)

    try:
        asm.assemble_linear()
        AssertionError("Clashing primers should give ValueError")
    except ValueError:
        pass

    # PCR on circular assemblies should not give an error if primers clash
    seq = seq.looped()
    for shift in range(len(seq)):
        seq_shifted = seq.shifted(shift)
        asm = assembly.PCRAssembly([primer1, seq_shifted, primer2], limit=8)
        asm.assemble_linear()

    # PCR assembly, raises error in case of circular assembly and insertion
    with pytest.raises(NotImplementedError):
        asm.get_circular_assemblies()

    # Error if length is not multiple of 3
    with pytest.raises(ValueError):
        assembly.PCRAssembly([primer1, seq], limit=8)

    # Error if not primer, seq, primer
    with pytest.raises(ValueError):
        assembly.PCRAssembly([seq, seq, seq], limit=8)
    with pytest.raises(ValueError):
        assembly.PCRAssembly([primer1, primer1, primer1], limit=8)

    # PCR assembly does not support only_adjacent_edges
    with pytest.raises(NotImplementedError):
        asm.get_linear_assemblies(only_adjacent_edges=True)

    # PCR assembly does not support insertion
    with pytest.raises(NotImplementedError):
        asm.get_insertion_assemblies()


def test_annotate_primer_binding_sites():
    primer1 = Primer("ACGTACGT", name="primer1")
    primer2 = Primer(reverse_complement("GCGCGCGC"), name="primer2")

    seq = Dseqrecord("ccccACGTACGTAAAAAAGCGCGCGCcccc")

    asm = assembly.PCRAssembly([primer1, seq, primer2], limit=8)

    prods = asm.assemble_linear()

    assert len(prods) == 1
    annotated = assembly.annotate_primer_binding_sites(
        prods[0], [primer1, seq, primer2]
    )
    assert annotated.features[0].type == "primer_bind"
    assert annotated.features[0].qualifiers["label"] == [primer1.name]
    assert annotated.features[0].qualifiers["note"] == ["sequence: " + str(primer1.seq)]
    assert annotated.features[0].location.start == 0
    assert annotated.features[0].location.end == len(primer1)
    assert annotated.features[1].type == "primer_bind"
    assert annotated.features[1].qualifiers["label"] == [primer2.name]
    assert annotated.features[1].qualifiers["note"] == ["sequence: " + str(primer2.seq)]
    assert annotated.features[1].location.start == len(prods[0]) - len(primer2)
    assert annotated.features[1].location.end == len(prods[0])

    # With tails
    primer1 = Primer("aaACGTACGT", name="primer1")
    primer2 = Primer(reverse_complement("aaGCGCGCGC"), name="primer2")

    seq = Dseqrecord("ccccACGTACGTAAAAAAGCGCGCGCcccc")

    asm = assembly.PCRAssembly([primer1, seq, primer2], limit=8)

    prods = asm.assemble_linear()

    assert len(prods) == 1
    annotated = assembly.annotate_primer_binding_sites(
        prods[0], [primer1, seq, primer2]
    )
    assert annotated.features[0].type == "primer_bind"
    assert annotated.features[0].qualifiers["label"] == [primer1.name]
    assert annotated.features[0].qualifiers["note"] == ["sequence: " + str(primer1.seq)]
    assert annotated.features[0].location.start == 0
    assert annotated.features[0].location.end == len(primer1)
    assert annotated.features[1].type == "primer_bind"
    assert annotated.features[1].qualifiers["label"] == [primer2.name]
    assert annotated.features[1].qualifiers["note"] == ["sequence: " + str(primer2.seq)]
    assert annotated.features[1].location.start == len(prods[0]) - len(primer2)
    assert annotated.features[1].location.end == len(prods[0])


def test_fragments_only_once():

    fragments = [
        Dseqrecord("TTTTacgatAAtgctccCCCC", circular=False),
        Dseqrecord("CCCCtcatGGGG", circular=False),
        Dseqrecord("GGGGatataTTTT", circular=False),
    ]

    asm = assembly.Assembly(
        fragments,
        limit=4,
        algorithm=assembly.gibson_overlap,
        use_all_fragments=True,
        use_fragment_order=False,
    )
    for a in asm.get_linear_assemblies():
        nodes_used = [
            f[0]
            for f in assembly.edge_representation2subfragment_representation(a, False)
        ]
        assert len(nodes_used) == len(set(nodes_used))

    for a in asm.get_circular_assemblies():
        nodes_used = [
            f[0]
            for f in assembly.edge_representation2subfragment_representation(a, True)
        ]
        assert len(nodes_used) == len(set(nodes_used))


def test_both_representations():
    # Linear examples
    assembly_linear = (
        (1, 2, "loc_1_r", "loc_2_l"),
        (2, 3, "loc_2_r", "loc_3_l"),
    )

    subf = assembly.edge_representation2subfragment_representation(
        assembly_linear, False
    )
    assert subf == (
        (1, None, "loc_1_r"),
        (2, "loc_2_l", "loc_2_r"),
        (3, "loc_3_l", None),
    )

    back_to_edge = assembly.subfragment_representation2edge_representation(subf, False)
    assert back_to_edge == assembly_linear

    assembly_linear = (
        (1, -2, "loc_1_r", "loc_2_l"),
        (-2, -3, "loc_2_r", "loc_3_l"),
    )

    subf = assembly.edge_representation2subfragment_representation(
        assembly_linear, False
    )
    assert subf == (
        (1, None, "loc_1_r"),
        (-2, "loc_2_l", "loc_2_r"),
        (-3, "loc_3_l", None),
    )

    back_to_edge = assembly.subfragment_representation2edge_representation(subf, False)
    assert back_to_edge == assembly_linear

    # Circular example

    assembly_circular = (
        (1, 2, "loc_1_r", "loc_2_l"),
        (2, 3, "loc_2_r", "loc_3_l"),
        (3, 1, "loc_3_r", "loc_1_l"),
    )

    subf = assembly.edge_representation2subfragment_representation(
        assembly_circular, True
    )
    assert subf == (
        (1, "loc_1_l", "loc_1_r"),
        (2, "loc_2_l", "loc_2_r"),
        (3, "loc_3_l", "loc_3_r"),
    )

    back_to_edge = assembly.subfragment_representation2edge_representation(subf, True)
    assert back_to_edge == assembly_circular


def test_ends_from_cutsite():

    seqs = (
        Dseq("AAAGAATTCAAA", circular=True),
        Dseq("CCGCGATCGCCC", circular=True),
        Dseq("AAGATATCAA", circular=True),
    )
    enzymes = (EcoRI, RgaI, EcoRV)
    expected_results = (("5'", "aatt"), ("3'", "at"), ("blunt", ""))

    for seq, enz, res in zip(seqs, enzymes, expected_results):
        for shift in range(len(seqs)):
            seq_shifted = seq.shifted(shift)
            cut = seq_shifted.get_cutsites(enz)[0]
            assert assembly.ends_from_cutsite(cut, seq_shifted) == (res, res)

    a = Dseq("AAGGTCTCAacccAA", circular=False)

    cut = a.get_cutsites([BsaI])[0]
    a1, a2 = a.cut([BsaI])
    assert assembly.ends_from_cutsite(cut, a) == (
        a1.three_prime_end(),
        a2.five_prime_end(),
    )

    a = Dseq("ACcTGatacACTGGATactA", circular=False)

    cut = a.get_cutsites([BsrI])[0]
    a1, a2 = a.cut([BsrI])
    assert assembly.ends_from_cutsite(cut, a) == (
        a1.three_prime_end(),
        a2.five_prime_end(),
    )

    with pytest.raises(ValueError):
        assembly.ends_from_cutsite(None, Dseq("ATCG"))


def test_restriction_ligation_assembly():

    seq_pairs = (
        (Dseqrecord("AAAGAATTCAAA"), Dseqrecord("CCCCGAATTCCCC")),
        (Dseqrecord("AAAGCGATCGCAAA"), Dseqrecord("CCCCGCGATCGCCCC")),
        (Dseqrecord("acaGATATCtta"), Dseqrecord("aaaGATATCata")),
    )
    enzymes = [EcoRI, RgaI, EcoRV]
    for (a, b), enz in zip(seq_pairs, enzymes):
        a1, a2 = a.cut([enz])
        b1, b2 = b.cut([enz])

        products = [
            a1 + b2,
            a1 + b1.reverse_complement(),
            b1 + a2,
            a2.reverse_complement() + b2,
        ]

        products2 = assembly.restriction_ligation_assembly([a, b], [enz])

        assert sorted(dseqrecord_list_to_dseq_list(products)) == sorted(
            dseqrecord_list_to_dseq_list(products2)
        )

    # Insertion in a vector
    f1 = Dseqrecord("GAATTCaaaGAATTC")
    f2 = Dseqrecord("CCCCGAATTCCCC", circular=True)

    _, a2, a3 = f1.cut([EcoRI])
    (b1,) = f2.cut([EcoRI])
    result_seguids = sorted(
        [
            (b1.seq + a2.seq).looped().seguid(),
            (b1.seq + a2.seq.reverse_complement()).looped().seguid(),
        ]
    )

    # We shift
    for shift in range(len(f2)):
        f2_shifted = f2.shifted(shift)
        products2 = assembly.restriction_ligation_assembly(
            [f1, f2_shifted], [EcoRI], circular_only=True
        )
        observed_seguids = sorted(x.seguid() for x in products2)
        assert len(result_seguids) == len(observed_seguids)
        assert result_seguids == observed_seguids

    # Cutting from plasmid
    f1 = Dseqrecord("aaGAATTCaaaGTCGACaa", circular=True)
    f2 = Dseqrecord("ccGAATTCccGTCGACc")

    f1_0, f1_1 = f1.cut([EcoRI, SalI])
    f2_0, f2_1, f2_3 = f2.cut([EcoRI, SalI])

    result_seguids = sorted(
        [
            (f1_1 + f2_1).looped().seguid(),
            (f1_0 + f2_1.reverse_complement()).looped().seguid(),
        ]
    )

    # We shift
    for shift in range(len(f1)):
        f1_shifted = f1.shifted(shift)
        products2 = assembly.restriction_ligation_assembly(
            [f1_shifted, f2], [EcoRI, SalI], circular_only=True
        )
        observed_seguids = sorted(x.seguid() for x in products2)
        assert len(result_seguids) == len(observed_seguids)
        assert result_seguids == observed_seguids

    # TODO: Check if features are transferred properly

    # Partial overlaps -> enzyme with negative overhang
    fragments = [Dseqrecord("GGTCTCCCCAATT"), Dseqrecord("GGTCTCCAACCAA")]

    # TODO: Needs fixing, related to https://github.com/pydna-group/pydna/issues/426
    # Not allowing partial overlaps
    def algo(x, y, _l):
        return assembly.restriction_ligation_overlap(x, y, [BsaI], False)

    f = assembly.Assembly(fragments, use_fragment_order=False, algorithm=algo)
    assert len(f.get_linear_assemblies(only_adjacent_edges=True)) == 0

    # Allowing partial overlaps
    def algo(x, y, _l):
        return assembly.restriction_ligation_overlap(x, y, [BsaI], partial=True)

    f = assembly.Assembly(fragments, algorithm=algo, use_fragment_order=False)
    assert len(f.get_linear_assemblies(only_adjacent_edges=True)) == 2
    p1, p2 = f.assemble_linear(only_adjacent_edges=True)
    assert str(p1.seq) == "GGTCTCCCCAACCAA"
    assert str(p2.seq) == "GGTCTCCAACCAATT"

    # Combining both partial and normal overlaps, to ensure that only_adjacent_edges keeps both.
    # In this case use_all_fragments=2, it should return 2 assemblies for 1 + 2, and
    fragments = [
        Dseqrecord("GGTCTCCCCAATT"),
        Dseqrecord("GGTCTCCAACCAA"),
        Dseqrecord("GGTCTCCCCAATT"),
    ]
    f = assembly.Assembly(
        fragments, algorithm=algo, use_fragment_order=False, use_all_fragments=False
    )
    products = f.assemble_linear(only_adjacent_edges=True)
    assert len(products) == 6
    products_seguid = set(p.seq.seguid() for p in products)
    assert products_seguid == set(
        [p1.seq.seguid(), p2.seq.seguid(), Dseqrecord("GGTCTCCCCAATT").seguid()]
    )

    # Partial overlaps -> enzyme with positive overhang
    fragments = [Dseqrecord("GACACCAGAGTC"), Dseqrecord("GACTAACGGGTC")]

    assert len(assembly.restriction_ligation_assembly(fragments, [DrdI])) == 0

    # Allowing partial overlaps
    def algo(x, y, _l):
        return assembly.restriction_ligation_overlap(x, y, [DrdI], True)

    f = assembly.Assembly(fragments, algorithm=algo, use_fragment_order=False)
    products = f.assemble_linear()
    assert str(products[0].seq) == "GACACCACGGGTC"
    assert str(products[1].seq) == "GACTAACAGAGTC"

    # Single fragment assemblies

    f1 = Dseqrecord("aaGAATTCtttGAATTCaa", circular=True)
    products = assembly.restriction_ligation_assembly([f1], [EcoRI], circular_only=True)
    assert len(products) == 2
    assert str(products[0].seq) == "AATTCaaaaG"
    assert str(products[1].seq) == "AATTCtttG"

    f1 = Dseqrecord("aaGAATTCtttGAATTCaa", circular=False)
    products = assembly.restriction_ligation_assembly(
        [f1], [EcoRI], circular_only=False
    )
    assert len(products) == 2
    assert str(products[1].seq) == "aaGAATTCaa"
    assert str(products[0].seq) == "AATTCtttG"

    # Mixing blunt and normal overhangs
    fragments = [Dseqrecord("aaaGATATCccGAATTCaa"), Dseqrecord("cgcGATATCataGAATTCtta")]
    products = assembly.restriction_ligation_assembly(
        fragments, [EcoRI, EcoRV], circular_only=True
    )
    assert len(products) == 1
    assert str(products[0].seq) == "ATCccGAATTCtatGAT"


@pytest.mark.xfail(reason="See https://github.com/pydna-group/pydna/issues/426")
def test_restriction_ligation_partial_overlaps():

    # Partial overlaps
    fragments = [Dseqrecord("GGTCTCCCCAATT"), Dseqrecord("GGTCTCCAACCAA")]
    products = assembly.restriction_ligation_assembly(
        fragments, [BsaI], allow_partial_overlap=True
    )

    assert len(products) == 2
    assert str(products[0].seq) == "GGTCTCCCCAACCAA"
    assert str(products[1].seq) == "GGTCTCCAACCAATT"


def test_only_adjacent_edges():
    """Tests that partially digested fragments are not used in the assembly"""
    # Partial overlaps -> enzyme with negative overhang
    fragments = [Dseqrecord("acaGAATTCcccGAATTCtta"), Dseqrecord("aaaGAATTCata")]

    def algo(x, y, _l):
        return assembly.restriction_ligation_overlap(x, y, [EcoRI])

    f = assembly.Assembly(fragments, algorithm=algo, use_fragment_order=False)

    assert len(f.get_linear_assemblies(only_adjacent_edges=True)) == 4


def test_golden_gate():
    """This also tests that no partial deletionsar"""

    # Circular assembly
    insert1 = Dseqrecord("GGTCTCAattaAAAAAttaaAGAGACC")
    insert2 = Dseqrecord("GGTCTCAttaaCCCCCatatAGAGACC")
    insert3 = Dseqrecord("GGTCTCAatatGGGGGccggAGAGACC")

    vector = Dseqrecord("TTTTattaAGAGACCTTTTTGGTCTCAccggTTTT", circular=True)

    i1_pre, i1, i1_post = insert1.cut(BsaI)
    _, i2, _ = insert2.cut(BsaI)
    i3_pre, i3, i3_post = insert3.cut(BsaI)
    _, v = vector.cut(BsaI)

    assembly_output = assembly.golden_gate_assembly(
        [insert1, insert2, insert3, vector], [BsaI], circular_only=True
    )

    assert len(assembly_output) == 1
    assert assembly_output[0].seguid() == (i1 + i2 + i3 + v).looped().seguid()


def test_gibson_assembly():
    # For the test we pass input fragments + expected output
    test_cases = [
        (["AGAGACCaaaAGAGACC"], ["AGAGACCaaa"]),
        (["GTCGACTaaaAGAGACC", "AGAGACCcgcGTCGACT"], ["GTCGACTaaaAGAGACCcgc"]),
    ]

    # Should return the same thing for gibson and equivalent functions:
    for gibson_like_function in [
        assembly.gibson_assembly,
        assembly.in_fusion_assembly,
        assembly.fusion_pcr_assembly,
    ]:
        for fragments_str, expected_outputs in test_cases:
            for mode in range(3):
                if mode == 0:
                    # No overhangs
                    fragments = [Dseqrecord(f) for f in fragments_str]
                elif mode == 1:
                    # 3' overhangs (should give the same results as no overhangs)
                    fragments = [
                        Dseqrecord(Dseq.from_full_sequence_and_overhangs(f, 3, 3))
                        for f in fragments_str
                    ]
                else:
                    # Add 5' overhangs that will be removed in Gibson, so should give same results as no overhangs
                    fragments = [
                        Dseqrecord(
                            Dseq.from_full_sequence_and_overhangs(
                                "aaa" + f + "aaa", -3, -3
                            )
                        )
                        for f in fragments_str
                    ]
                products = gibson_like_function(fragments, 7, circular_only=True)
                products_str = [str(p.seq) for p in products]
                assert products_str == expected_outputs


def test_insertion_assembly():

    # Insertion of linear sequence into linear sequence (like
    # homologous recombination of PCR product with homology arms in genome)
    a = Dseqrecord("CGTACGCACAbbbbCGTACGCACAC")
    b = Dseqrecord("CGTACGCACArrrrCGTACGCACAT")

    f = assembly.Assembly([a, b], use_fragment_order=False, limit=10)

    # All possibilities, including the single insertions
    results = [
        "CGTACGCACArrrrCGTACGCACAbbbbCGTACGCACAC",
        "CGTACGCACArrrrCGTACGCACAC",
        "CGTACGCACAbbbbCGTACGCACArrrrCGTACGCACAT",
        "CGTACGCACAbbbbCGTACGCACArrrrCGTACGCACAC",
        "CGTACGCACAbbbbCGTACGCACAT",
        "CGTACGCACArrrrCGTACGCACAbbbbCGTACGCACAT",
    ]

    assembly_products = [
        str(assembly.assemble([a, b], assem).seq)
        for assem in f.get_insertion_assemblies()
    ]
    assert sorted(assembly_products) == sorted(results)

    # TODO: debatable whether this kind of homologous recombination should happen, or how
    # the overlap restrictions should be applied.

    a = Dseqrecord("CGTACGCACAbbbbC")
    b = Dseqrecord("CGTACGCACArrrrCGTACGCACAT")
    f = assembly.Assembly([a, b], use_fragment_order=False, limit=10)
    results = ["CGTACGCACArrrrCGTACGCACAbbbbC"]
    for assem, result in zip(f.get_insertion_assemblies(), results):
        assert result == str(assembly.assemble([a, b], assem).seq)

    a = Dseqrecord("CGTACGCACAbbbbC")
    b = Dseqrecord("CGTACGCACArrrrT")
    f = assembly.Assembly([a, b], use_fragment_order=False, limit=10)
    assert len(f.get_insertion_assemblies()) == 0

    # Does not work for circular molecules
    a = Dseqrecord("CGTACGCACAbbbbCGTACGCACAC", circular=True)
    b = Dseqrecord("CGTACGCACArrrrCGTACGCACAT", circular=True)
    assert (
        assembly.Assembly(
            [a, b], use_fragment_order=False, limit=10
        ).get_insertion_assemblies()
        == []
    )

    a = Dseqrecord("CGTACGCACAbbbbC", circular=True)
    b = Dseqrecord("CGTACGCACArrrrCGTACGCACAT", circular=True)
    assert (
        assembly.Assembly(
            [a, b], use_fragment_order=False, limit=10
        ).get_insertion_assemblies()
        == []
    )

    a = Dseqrecord("CGTACGCACAbbbbC", circular=True)
    b = Dseqrecord("CGTACGCACArrrrT", circular=True)
    assert (
        assembly.Assembly(
            [a, b], use_fragment_order=False, limit=10
        ).get_insertion_assemblies()
        == []
    )

    # Only the right order is returned
    a = Dseqrecord("ttACGTTCGTccccTTAATTAAcc", circular=False)
    b = Dseqrecord("ttACGTTCGTttttCGGGCGCGaa", circular=True)
    c = Dseqrecord("aaCGGGCGCGggggTTAATTAAaa", circular=True)
    fragments = [a, b, c]
    for i in range(3):
        asm = assembly.Assembly(
            fragments[i:] + fragments[:i],
            use_fragment_order=False,
            limit=8,
            use_all_fragments=True,
        )
        prods = asm.assemble_insertion()

        assert len(prods) == 1
        assert str(prods[0].seq) == "ttACGTTCGTttttCGGGCGCGggggTTAATTAAcc"

    # Insertion / Excision of linear and circular molecules
    f = Dseqrecord("aaCGGGCGCGggggCGGGCGCGac", circular=False)
    asm = assembly.SingleFragmentAssembly([f], limit=8)
    prods = asm.assemble_insertion()
    assert len(prods) == 1
    assert str(prods[0].seq) == "aaCGGGCGCGac"
    prods = asm.assemble_circular()
    assert len(prods) == 1
    assert str(prods[0].seq) == "CGGGCGCGgggg"


def circles_assembly():
    a = Dseqrecord("xxxACGTAyyy", circular=True)
    b = Dseqrecord("bbACGTAbb", circular=True)

    f = assembly.Assembly([a, b], use_fragment_order=False, limit=5)
    circular_assemblies = f.get_circular_assemblies()
    assert len(circular_assemblies) == 1
    assert (
        str(assembly.assemble([a, b], circular_assemblies[0])) == "ACGTAyyyxxxACGTAbbbb"
    )

    # When more than two are provided, sequential homologous recombinations are returned
    a = Dseqrecord("aaACGTAACGTAaa", circular=True)
    b = Dseqrecord("ccACGTAACGTAcc", circular=True)
    c = Dseqrecord("ggACGTAACGTAgg", circular=True)

    f = assembly.Assembly(
        [a, b, c], use_fragment_order=False, limit=10, use_all_fragments=True
    )
    # All possibilities, including the single insertions
    results = [
        "ACGTAACGTAaaaaACGTAACGTAccccACGTAACGTAgggg",
        "ACGTAACGTAaaaaACGTAACGTAggggACGTAACGTAcccc",
    ]
    for assem, result in zip(f.get_insertion_assemblies(), results):
        assert result == str(assembly.assemble([a, b, c], assem).seq)

    # TODO: debatable whether this kind of homologous recombination should happen, same as above
    a = Dseqrecord("xxxACGTAyyy", circular=True)
    b = Dseqrecord("bbACGTAbbACGTAbb", circular=False)

    f = assembly.Assembly([a, b], use_fragment_order=False, limit=5)
    assert len(circular_assemblies) == 1
    assert (
        str(assembly.assemble([a, b], circular_assemblies[0])) == "ACGTAyyyxxxACGTAbb"
    )


def test_assemble_function():
    """A more granular test of the assemble function, independent of the experimental
    setup to make sure it works for all topologies"""

    f1 = Dseqrecord("aaaTTTctaGGGccc", circular=True)
    f2 = Dseqrecord("ccccTTTatgGGGaaa")

    # TTT features
    f1_feat1 = SeqFeature(SimpleLocation(3, 6))
    f2_feat1 = SeqFeature(SimpleLocation(4, 7))

    # GGG features
    f1_feat2 = SeqFeature(SimpleLocation(9, 12))
    f2_feat2 = SeqFeature(SimpleLocation(10, 13))

    f1.features = [f1_feat1, f1_feat2]
    f2.features = [f2_feat1, f2_feat2]

    for shift in range(len(f1)):

        f1_shifted = f1.shifted(shift)

        # Re-order the features so that TTT is first
        if str(f1_shifted.features[0].location.extract(f1_shifted.seq)) != "TTT":
            f1_shifted.features = f1_shifted.features[::-1]

        # Linear assembly 2 - 1 - 2 (ccccTTTctaGGGaaa)
        assembly_plan = [
            (2, 1, f2.features[0].location, f1_shifted.features[0].location),
            (1, 2, f1_shifted.features[1].location, f2.features[1].location),
        ]

        result = assembly.assemble([f1_shifted, f2], assembly_plan)
        assert str(result.seq) == "ccccTTTctaGGGaaa"
        assert len(result.features) == 4
        assert set(str(f.location.extract(result.seq)) for f in result.features) == {
            "TTT",
            "GGG",
        }

        # Circular assembly 1 - 2 (ccccTTTctaGGGaaa)
        assembly_plan = [
            (1, 2, f1_shifted.features[0].location, f2.features[0].location),
            (2, 1, f2.features[1].location, f1_shifted.features[1].location),
        ]

        result = assembly.assemble([f1_shifted, f2], assembly_plan)
        assert str(result.seq) == "GGGcccaaaTTTatg"
        assert len(result.features) == 4
        assert set(str(f.location.extract(result.seq)) for f in result.features) == {
            "TTT",
            "GGG",
        }

        # TODO: This type of assembly should maybe raise an error, no
        # linear assembly should start or finish with a circular sequence
        # assembly_plan = [
        #     (2, 1, f2.features[1].location, f1_shifted.features[1].location),
        # ]
        # result = assembly.assemble([f1_shifted, f2], assembly_plan)

    # Now both are circular, using a single insertion site
    f1 = Dseqrecord("aaaTTTcta", circular=True)
    f2 = Dseqrecord("ccccTTTatg", circular=True)
    f1.features = [SeqFeature(SimpleLocation(3, 6))]
    f2.features = [SeqFeature(SimpleLocation(4, 7))]

    for shift_1 in range(len(f1)):
        f1_shifted = f1.shifted(shift_1)

        for shift_2 in range(len(f2)):
            f2_shifted = f2.shifted(shift_2)
            # Linear assembly 2 - 1 - 2 (ccccTTTctaGGGaaa)
            assembly_plan = [
                (
                    1,
                    2,
                    f1_shifted.features[0].location,
                    f2_shifted.features[0].location,
                ),
                (
                    2,
                    1,
                    f2_shifted.features[0].location,
                    f1_shifted.features[0].location,
                ),
            ]

            result = assembly.assemble([f1_shifted, f2_shifted], assembly_plan)
            assert (
                result.seq.seguid()
                == Dseq("aaaTTTatgccccTTTcta", circular=True).seguid()
            )
            assert len(result.features) == 4
            assert set(
                str(f.location.extract(result.seq)) for f in result.features
            ) == {"TTT"}

    # Blunt assemblies
    fragments = [Dseqrecord("aaaTTTctaGGGccc"), Dseqrecord("ccccTTTatgGGGaa")]
    loc_end = SimpleLocation(15, 15)
    loc_start = SimpleLocation(0, 0)

    # A linear assembly
    assembly_plan = [
        (1, 2, loc_end, loc_start),
    ]

    assert (fragments[0] + fragments[1]).seq == assembly.assemble(
        fragments, assembly_plan
    ).seq

    # A circular assembly
    assembly_plan = [
        (1, 2, loc_end, loc_start),
        (2, 1, loc_end, loc_start),
    ]

    assert (fragments[0] + fragments[1]).looped().seq == assembly.assemble(
        fragments, assembly_plan
    ).seq


def test_assembly_is_valid():

    # fragments are only used for topology, their sequence is not used
    fragments = [Dseqrecord(""), Dseqrecord(""), Dseqrecord("")]

    # Normal assembly plan
    #
    # 1 ------
    #      |||
    # 2    ------
    #         |||
    # 3       ------
    #
    assembly_plan = [
        (1, 2, SimpleLocation(3, 6), SimpleLocation(0, 3)),
        (2, 3, SimpleLocation(3, 6), SimpleLocation(0, 3)),
    ]

    assert assembly.Assembly.assembly_is_valid(fragments, assembly_plan, False, True)

    # Partial overlap should be allowed, a fragment could act like a "bridge"
    # 1 ------
    #     ||||
    # 2   ------
    #        ||||
    # 3      ------
    # TODO: check meaning of this for non-homology assemblies like restriction ligation
    assembly_plan = [
        (1, 2, SimpleLocation(2, 6), SimpleLocation(0, 4)),
        (2, 3, SimpleLocation(2, 6), SimpleLocation(0, 4)),
    ]

    assert assembly.Assembly.assembly_is_valid(fragments, assembly_plan, False, True)

    # Complete overlap in linear assemblies should be discarded as is redundant
    # 1 ------
    #     ||||
    # 2   ------
    #     ||||
    # 3   ------
    # TODO: check meaning of this when 2 is circular molecule

    assembly_plan = [
        (1, 2, SimpleLocation(2, 6), SimpleLocation(0, 4)),
        (2, 3, SimpleLocation(0, 4), SimpleLocation(0, 4)),
    ]

    assert not assembly.Assembly.assembly_is_valid(
        fragments, assembly_plan, False, True
    )

    # Invalid assembly
    # 1   ------
    #         ||
    # 2   ------
    #     ||
    # 3   ------

    assembly_plan = [
        (1, 2, SimpleLocation(4, 6), SimpleLocation(4, 6)),
        (2, 3, SimpleLocation(0, 2), SimpleLocation(0, 2)),
    ]

    assert not assembly.Assembly.assembly_is_valid(
        fragments, assembly_plan, False, True
    )

    # Assembly plan including two fragments extracted from circular molecules:
    f1 = Dseqrecord("ccTTTc")
    f2 = Dseqrecord("TTTAAA", circular=True)
    f3 = Dseqrecord("AAACCC", circular=True)
    f4 = Dseqrecord("ggCCCg")
    f1.features = [SeqFeature(SimpleLocation(2, 5), id="f1_f2")]
    f2.features = [
        SeqFeature(SimpleLocation(0, 3), id="f1_f2"),
        SeqFeature(SimpleLocation(3, 6), id="f2_f3"),
    ]
    f3.features = [
        SeqFeature(SimpleLocation(0, 3), id="f2_f3"),
        SeqFeature(SimpleLocation(3, 6), id="f3_f4"),
    ]
    f4.features = [SeqFeature(SimpleLocation(2, 5), id="f3_f4")]

    def find_feature_by_id(f: Dseqrecord, id: str) -> SeqFeature:
        return next(f for f in f.features if f.id == id)

    for shift_2 in range(len(f2)):
        f2_shifted = f2.shifted(shift_2)

        for shift_3 in range(len(f3)):
            f3_shifted = f3.shifted(shift_3)
            fragments = [f1, f2_shifted, f3_shifted, f4]
            assembly_plan = [
                (
                    1,
                    2,
                    f1.features[0].location,
                    find_feature_by_id(f2_shifted, "f1_f2").location,
                ),
                (
                    2,
                    3,
                    find_feature_by_id(f2_shifted, "f2_f3").location,
                    find_feature_by_id(f3_shifted, "f2_f3").location,
                ),
                (
                    3,
                    4,
                    find_feature_by_id(f3_shifted, "f3_f4").location,
                    f4.features[0].location,
                ),
            ]
            assert assembly.Assembly.assembly_is_valid(
                fragments, assembly_plan, False, True
            )
            # Does not really belong here, but
            assert (
                str(assembly.assemble(fragments, assembly_plan).seq) == "ccTTTAAACCCg"
            )

    # is_circular must be set
    assert not assembly.Assembly.assembly_is_valid(fragments, assembly_plan, None, True)

    # In a circular assembly, first and last fragment must be the same
    assembly_plan[0] = (1, 2, SimpleLocation(0, 3), SimpleLocation(0, 3))
    assert not assembly.Assembly.assembly_is_valid(fragments, assembly_plan, True, True)

    # In a linear assembly, the first and last fragments cannot be circular
    fragments = [Dseqrecord("", circular=True), Dseqrecord("")]
    assembly_plan = [(1, 2, SimpleLocation(0, 3), SimpleLocation(0, 3))]
    assert not assembly.Assembly.assembly_is_valid(
        fragments, assembly_plan, False, True
    )


def test_extract_subfragment():
    def find_feature_by_id(f: Dseqrecord, id: str) -> SeqFeature:
        return next(f for f in f.features if f.id == id)

    f1 = Dseqrecord("aaTTTcccTTTaa", circular=True)
    f1.features = [
        SeqFeature(SimpleLocation(2, 5), id="left"),
        SeqFeature(SimpleLocation(8, 11), id="right"),
    ]

    for shift in range(len(f1)):
        f1_shifted = f1.shifted(shift)
        left = find_feature_by_id(f1_shifted, "left").location
        right = find_feature_by_id(f1_shifted, "right").location
        subfragment = assembly.extract_subfragment(f1_shifted, left, right)
        assert str(subfragment.seq) == "TTTcccTTT"
        assert len(subfragment.features) == 2

    # Old edge case, in circular molecules, if you extract the entire sequence
    # it leads to a [0:0] getitem call, before it gave an empty sequence, see https://github.com/BjornFJohansson/pydna/issues/161
    f1 = Dseqrecord("AAATTT", circular=True)
    f1.features = [
        SeqFeature(SimpleLocation(0, 3), id="left"),
        SeqFeature(SimpleLocation(3, 6), id="right"),
    ]

    for shift in range(len(f1)):
        f1_shifted = f1.shifted(shift)
        left = find_feature_by_id(f1_shifted, "left").location
        right = find_feature_by_id(f1_shifted, "right").location
        subfragment = assembly.extract_subfragment(f1_shifted, left, right)
        assert str(subfragment.seq) == "AAATTT"
        assert len(subfragment.features) == 2

    # In circular molecules, the same feature twice should extract the "opened up" sequence, as
    # when you cut open a plasmid with a restriction enzyme. This is useful to represent an integration
    # in which the integration seq is duplicated, or for digestion/ligation representing the opening up
    # of the plasmid.
    f1 = Dseqrecord("ATTTA", circular=True)
    f1.features = [SeqFeature(SimpleLocation(1, 4))]

    for shift in range(len(f1)):
        f1_shifted = f1.shifted(shift)
        loc = f1_shifted.features[0].location
        subfragment = assembly.extract_subfragment(f1_shifted, loc, loc)
        assert str(subfragment.seq) == "TTTAATTT"
        # The feature marking the overlap should be copied
        assert len(subfragment.features) == 2

    # Blunt ends assembly
    f1 = Dseqrecord("ATTTA")
    loc_start = SimpleLocation(0, 0)
    loc_end = SimpleLocation(5, 5)
    # Here we compare the seq because getitem changes the id of the sequence
    assert f1.seq == assembly.extract_subfragment(f1, loc_start, loc_end).seq

    f1 = Dseqrecord("ATTTA", circular=True)
    # Returns error if start or end location is None in circular sequences
    pytest.raises(ValueError, assembly.extract_subfragment, f1, None, None)
    pytest.raises(
        ValueError, assembly.extract_subfragment, f1, None, SimpleLocation(0, 5)
    )
    pytest.raises(
        ValueError, assembly.extract_subfragment, f1, SimpleLocation(0, 5), None
    )


def test_sticky_end_sub_strings():

    # Full overlap, negative ovhg
    a, b = Dseqrecord("AAAGAATTCAAA").cut(EcoRI)
    assert assembly.sticky_end_sub_strings(a, b) == [(4, 0, 4)]
    assert (
        assembly.sticky_end_sub_strings(a.reverse_complement(), b.reverse_complement())
        == []
    )
    assert assembly.sticky_end_sub_strings(b, a) == []

    # Full overlap, positive ovhg
    a, b = Dseqrecord("TTGCGATCGCTT").cut(RgaI)
    assert assembly.sticky_end_sub_strings(a, b) == [(5, 0, 2)]
    assert assembly.sticky_end_sub_strings(b, a) == []
    assert (
        assembly.sticky_end_sub_strings(a.reverse_complement(), b.reverse_complement())
        == []
    )

    # Blunt ends do not work either
    assert assembly.sticky_end_sub_strings(Dseqrecord("TT"), Dseqrecord("TT")) == []

    # Partial overlaps
    a = Dseqrecord(Dseq.from_full_sequence_and_overhangs("AAAGAA", 0, 3))
    b = Dseqrecord(Dseq.from_full_sequence_and_overhangs("AAAGAA", 3, 0))

    assert assembly.sticky_end_sub_strings(a, b) == []
    # Only when limit == True -> TODO: change this to not be an assembly parameter, but
    # functional instead.
    assert assembly.sticky_end_sub_strings(a, b, True) == [(4, 0, 2)]


def test_blunt_overlap():
    a = Dseqrecord(
        Dseq.from_full_sequence_and_overhangs("AAAAA", crick_ovhg=0, watson_ovhg=1)
    )
    b = Dseqrecord(
        Dseq.from_full_sequence_and_overhangs("CCCC", crick_ovhg=1, watson_ovhg=0)
    )

    assert assembly.blunt_overlap(a, b) == []
    assert assembly.blunt_overlap(a, b.reverse_complement()) == []
    assert assembly.blunt_overlap(a.reverse_complement(), b) == []
    assert assembly.blunt_overlap(b, a) == [(4, 0, 0)]
    assert assembly.blunt_overlap(a.reverse_complement(), b.reverse_complement()) == [
        (5, 0, 0)
    ]
    assert assembly.blunt_overlap(b.reverse_complement(), a.reverse_complement()) == []


def test_ligation_assembly():

    fragments = Dseqrecord("AAAGAATTCAAA").cut(EcoRI)
    assert assembly.ligation_assembly(fragments)[0].seq == Dseq("AAAGAATTCAAA")

    fragments = Dseqrecord("TTGCGATCGCTT").cut(RgaI)
    assert assembly.ligation_assembly(fragments)[0].seq == Dseq("TTGCGATCGCTT")

    # Circular ligation
    fragments = Dseqrecord("AAGAATTCTTGAATTCCC", circular=True).cut(EcoRI)
    expected_result = [
        (fragments[0] + fragments[1]).looped(),
        (fragments[0] + fragments[1].reverse_complement()).looped(),
    ]
    products = assembly.ligation_assembly(fragments, circular_only=True)
    assert sorted(dseqrecord_list_to_dseq_list(products)) == sorted(
        dseqrecord_list_to_dseq_list(expected_result)
    )

    # Partial ligation
    a = Dseqrecord(Dseq.from_full_sequence_and_overhangs("AAAGAA", 0, 3))
    b = Dseqrecord(Dseq.from_full_sequence_and_overhangs("AAAGAA", 3, 0))
    assert assembly.ligation_assembly([a, b]) == []
    assert (
        str(assembly.ligation_assembly([a, b], allow_partial_overlap=True)[0].seq)
        == "AAAGAAAGAA"
    )

    # Single fragment assemblies
    fragments = Dseqrecord("AAGAATTCTTGAATTCCC").cut(EcoRI)
    result = assembly.ligation_assembly([fragments[1]])
    assert len(result) == 1
    assert result[0].seq == fragments[1].looped().seq

    # Blunt ligation combined with sticky end
    fragments = Dseqrecord("AAAGAATTCAAA").cut(EcoRI)
    result = assembly.ligation_assembly(fragments, allow_blunt=True)
    result_str = [str(x.seq) for x in result]
    assert sorted(result_str) == sorted(["AAAGAATTCAAA"])
    assert result[0].circular


def test_blunt_assembly():
    # Linear assembly
    a = Dseqrecord(
        Dseq.from_full_sequence_and_overhangs("AAAAA", crick_ovhg=0, watson_ovhg=1)
    )
    b = Dseqrecord(
        Dseq.from_full_sequence_and_overhangs("CCCC", crick_ovhg=1, watson_ovhg=0)
    )

    asm = assembly.Assembly(
        [a, b],
        algorithm=assembly.blunt_overlap,
        use_all_fragments=True,
        use_fragment_order=False,
    )

    assert dseqrecord_list_to_dseq_list(asm.assemble_linear()) == [(b + a).seq]
    assert asm.assemble_circular() == []

    # Circular assembly
    a = Dseqrecord("ATT")
    b = Dseqrecord("CCCC")

    asm = assembly.Assembly(
        [a, b],
        algorithm=assembly.blunt_overlap,
        use_all_fragments=True,
        use_fragment_order=False,
    )

    assert dseqrecord_list_to_dseq_list(
        asm.assemble_linear()
    ) == dseqrecord_list_to_dseq_list(
        [a + b, a + b.reverse_complement(), b + a, a.reverse_complement() + b]
    )
    assert dseqrecord_list_to_dseq_list(
        asm.assemble_circular()
    ) == dseqrecord_list_to_dseq_list(
        [(a + b).looped(), (a + b.reverse_complement()).looped()]
    )

    # Circularisation
    asm = assembly.SingleFragmentAssembly(
        [Dseqrecord("AATT")], algorithm=assembly.blunt_overlap
    )
    result = asm.assemble_circular()
    assert len(result) == 1
    assert result[0].seq == Dseq("AATT", circular=True)


def test_format_insertion_assembly():

    loc1_a = SimpleLocation(2, 6)
    loc1_b = SimpleLocation(8, 12)
    loc2_a = SimpleLocation(0, 4)
    loc2_b = SimpleLocation(6, 10)

    seq1 = Dseqrecord("aaTTTTccTTTTaa")
    seq2 = Dseqrecord("TTTTccTTTT")

    fragments = [seq1, seq2]

    # This is just a dummy assembly planner, we only care for it containing 'fragments', which is used
    # by the class method
    assembly_planner = assembly.Assembly(fragments)
    asm_correct = [(1, 2, loc1_a, loc2_a), (2, 1, loc2_b, loc1_b)]
    assert asm_correct == assembly_planner.format_insertion_assembly(asm_correct)
    assert asm_correct == assembly_planner.format_insertion_assembly(asm_correct[::-1])

    # More fragments
    fragments = [seq1, seq2, seq2]
    assembly_planner = assembly.Assembly(fragments)
    asm_correct = [
        (1, 2, loc1_a, loc2_a),
        (2, 3, loc2_b, loc2_a),
        (3, 1, loc2_b, loc1_b),
    ]
    assert asm_correct == assembly_planner.format_insertion_assembly(asm_correct)
    assert asm_correct == assembly_planner.format_insertion_assembly(
        asm_correct[1:] + asm_correct[:1]
    )
    assert asm_correct == assembly_planner.format_insertion_assembly(
        asm_correct[2:] + asm_correct[:2]
    )

    asm_wrong = [(1, 2, loc1_b, loc2_a), (2, 3, loc2_b, loc2_a), (3, 1, loc2_b, loc1_a)]
    assert assembly_planner.format_insertion_assembly(asm_wrong) is None


def test_zip_leftwards():
    seq = Dseqrecord("AAAAACGTCCCGT")
    primer = Dseqrecord("ACGTCCCGT")
    match = (13, 9, 0)  # an empty match at the end of each
    assert (4, 0, 9) == assembly.zip_match_leftwards(seq, primer, match)

    # Works in circular molecules if the match spans the origin:
    seq = Dseqrecord("TCCCGTAAAAACG", circular=True)
    primer = Dseqrecord("ACGTCCCGT")
    match = (6, 9, 0)
    assert (10, 0, 9) == assembly.zip_match_leftwards(seq, primer, match)


def test_zip_rightwards():
    seq = Dseqrecord("AAAAACGTCCCGT")
    primer = Dseqrecord("ACGTCCCGT")
    match = (4, 0, 0)  # an empty match at the end of each
    assert (4, 0, 9) == assembly.zip_match_rightwards(seq, primer, match)

    # Works in circular molecules if the match spans the origin:
    seq = Dseqrecord("TCCCGTAAAAACG", circular=True)
    primer = Dseqrecord("ACGTCCCGT")
    match = (10, 0, 0)
    assert (10, 0, 9) == assembly.zip_match_rightwards(seq, primer, match)


def test_primer_template_overlap():
    template = Dseqrecord("AATTAGCAGCGATCGAGT", circular=True)
    for shift in range(len(template)):

        template_shifted = template.shifted(-shift)

        primer = Primer("TTAGCAGC")
        assert [
            (0, (2 + shift) % len(template), 8)
        ] == assembly.primer_template_overlap(primer, template_shifted, 8, 0)

        # The alignment is zipped if more bases align than the arg limit
        assert [
            (0, (2 + shift) % len(template), 8)
        ] == assembly.primer_template_overlap(primer, template_shifted, 6, 0)

        # Extra primer on the left
        primer_extra_left = Primer("GGGCTTTAGCAGC")
        assert [
            (5, (2 + shift) % len(template), 8)
        ] == assembly.primer_template_overlap(primer_extra_left, template_shifted, 8, 0)

        # Extra primer on the right does not work
        primer_extra_right = Primer("TTAGCAGCA")
        assert [] == assembly.primer_template_overlap(
            primer_extra_right, template_shifted, 8, 0
        )

        # Too short primer does not work either
        primer2short = Primer("TAGCAGC")
        assert [] == assembly.primer_template_overlap(
            primer2short, template_shifted, 8, 0
        )

        # Mismatch
        primer = Primer("TaAGCAGC")
        assert [
            (2, (4 + shift) % len(template), 6)
        ] == assembly.primer_template_overlap(primer, template_shifted, 8, 1)
        primer = Primer("aaAGCAGC")
        assert [
            (2, (4 + shift) % len(template), 6)
        ] == assembly.primer_template_overlap(primer, template_shifted, 8, 2)

        # Too many mismatches for argument
        primer = Primer("AAAGCAGC")
        assert [] == assembly.primer_template_overlap(primer, template_shifted, 8, 1)

        # Reverse primer
        primer = Primer("GCGATCGA")
        assert [
            ((8 + shift) % len(template), 0, 8)
        ] == assembly.primer_template_overlap(template_shifted, primer, 8, 0)

        # The alignment is zipped if more bases align than the arg limit
        assert [
            ((8 + shift) % len(template), 0, 8)
        ] == assembly.primer_template_overlap(template_shifted, primer, 6, 0)

        # Reverse primer with 5' extension
        primer = Primer("GCGATCGAAAAA")
        assert [
            ((8 + shift) % len(template), 0, 8)
        ] == assembly.primer_template_overlap(template_shifted, primer, 8, 0)

        # Extra bases on the 3' should not work
        primer = Primer("AAGCGATCGA")
        assert [] == assembly.primer_template_overlap(template_shifted, primer, 8, 0)

        # Mismatches
        primer = Primer("GCGtTCGA")
        assert [
            ((8 + shift) % len(template), 0, 3)
        ] == assembly.primer_template_overlap(template_shifted, primer, 8, 1)

    # Multiple matches
    primer = Primer("AATTAGCA")
    template = Dseqrecord("AATTAGCAGCGATCAATTAGCA")
    assert [(0, 0, 8), (0, 14, 8)] == assembly.primer_template_overlap(
        primer, template, 8, 1
    )

    # Gives the right error when passed sequences are not primer and sequence
    with pytest.raises(ValueError):
        assembly.primer_template_overlap(primer, primer, 8, 1)

    with pytest.raises(ValueError):
        assembly.primer_template_overlap(template, template, 8, 1)


def test_too_many_assemblies():
    fragments = [
        Dseqrecord("aaGCGGCCGCaaGCGGCCGC", circular=True),
        Dseqrecord("aaGCGGCCGCaaGCGGCCGC", circular=True),
        Dseqrecord("aaGCGGCCGCaaGCGGCCGCaaGCGGCCGCaaGCGGCCGC", circular=True),
    ]

    def algo(x, y, _l):
        return assembly.restriction_ligation_overlap(x, y, [NotI])  # noqa: B023

    asm = assembly.Assembly(fragments, algorithm=algo)

    with pytest.raises(ValueError):
        asm.assemble_circular()

    with pytest.raises(ValueError):
        asm.assemble_linear()

    with pytest.raises(ValueError):
        asm.get_insertion_assemblies()


def test_assembly_str():
    loc1_a = SimpleLocation(2, 6)
    loc1_b = SimpleLocation(8, 12)
    loc2_a = SimpleLocation(0, 4)
    loc2_b = SimpleLocation(6, 10)

    st = assembly.assembly2str([(1, 2, loc1_a, loc2_a), (2, 1, loc2_b, loc1_b)])
    assert st, "('1[2:6]:2[0:4]', '2[6:10]:1[8:12]')"

    st = assembly.assembly2str_tuple([(1, 2, loc1_a, loc2_a), (2, 1, loc2_b, loc1_b)])
    assert st, "((1, 2, '[2:6]', '[0:4]'), (2, 1, '[6:10]', '[8:12]'))"


def test_assembly_has_mismatches():
    fragments = [
        Dseqrecord("AAGAATTCTTGAATTCCC", circular=True),
        Dseqrecord("TCCCTTGAATTCCC", circular=True),
    ]

    for i in range(2):
        fragments[0].features = []
        fragments[1].features = []
        fragments[0].add_feature(14 - i, 18)
        fragments[1].add_feature(0, 4 + i)

        for shift in range(len(fragments[0])):
            seq1 = fragments[0].shifted(shift)
            for shift2 in range(len(fragments[1])):
                seq2 = fragments[1].shifted(shift2)
                asm = [(1, 2, seq1.features[0].location, seq2.features[0].location)]
                if i == 0:
                    assert not assembly.assembly_has_mismatches([seq1, seq2], asm)
                else:
                    assert assembly.assembly_has_mismatches([seq1, seq2], asm)


def test_single_fragment_assembly_error():
    fragments = [
        Dseqrecord("AAGAATTCTTGA"),
        Dseqrecord("TCCCTTGAATTCCC", circular=True),
    ]
    with pytest.raises(ValueError):
        assembly.SingleFragmentAssembly(
            fragments, algorithm=assembly.primer_template_overlap, limit=False
        )

    def algo(x, y, _l):
        return assembly.restriction_ligation_overlap(x, y, [EcoRI], False)

    f1 = Dseqrecord("aaGAATTCtttGAATTCaa", circular=False)
    f = assembly.SingleFragmentAssembly([f1], algorithm=algo)
    with pytest.raises(NotImplementedError):
        f.get_insertion_assemblies(only_adjacent_edges=True)
    with pytest.raises(NotImplementedError):
        f.get_linear_assemblies()

    f1 = Dseqrecord("GAATTCcatGAATTC", circular=False)
    f = assembly.SingleFragmentAssembly([f1], limit=5)

    assert len(f.assemble_insertion()) == 0
    assert len(f.assemble_circular()) == 1


def test_insertion_edge_case():
    a = Dseqrecord("cccgaggggaatcgaa")
    b = Dseqrecord("Acccgagggggaatc")

    asm = assembly.Assembly(
        [a, b], limit=5, use_all_fragments=True, use_fragment_order=False
    )

    product_seqs = set(str(prod.seq) for prod in asm.assemble_insertion())
    expected_seqs = {"cccgagggggaatcgaa", "Acccgaggggaatc"}
    assert product_seqs == expected_seqs

    b_circ = b.looped()

    for shift in range(len(b_circ)):
        b_shifted = b_circ.shifted(shift)
        asm = assembly.Assembly(
            [a, b_shifted], limit=5, use_all_fragments=True, use_fragment_order=False
        )

        product_seqs = [str(prod.seq) for prod in asm.assemble_insertion()]
        assert len(product_seqs) == 4
        assert "cccgagggggaatcgaa" in product_seqs


def test_common_sub_strings():

    a = Dseqrecord("RYBDKM", circular=True)
    for shift_1 in range(len(a)):
        a_shifted = a.shifted(shift_1)
        for shift_2 in range(len(a)):
            a_shifted_2 = a.shifted(shift_2)
            result = assembly.common_sub_strings(a_shifted, a_shifted_2, 3)
            assert len(result) == 1
            assert result[0][2] == 6


def test_in_vivo_assembly():
    # For the test we pass input fragments + expected output
    test_cases = [
        (["cccAGAGACCaaaAGAGACCttt"], ["AGAGACCaaa"]),
        (
            ["cccGTCGACTaaaAGAGACCttt", "cggAGAGACCcgcGTCGACTtta"],
            ["cGTCGACTaaaAGAGACCcg"],
        ),
    ]

    for fragments_str, expected_outputs in test_cases:
        for mode in range(3):
            if mode == 0:
                # No overhangs
                fragments = [Dseqrecord(f) for f in fragments_str]
            elif mode == 1:
                # 3' overhangs (should give the same results as no overhangs)
                fragments = [
                    Dseqrecord(Dseq.from_full_sequence_and_overhangs(f, 3, 3))
                    for f in fragments_str
                ]
            else:
                # Add 5' overhangs that will be removed in Gibson, so should give same results as no overhangs
                fragments = [
                    Dseqrecord(
                        Dseq.from_full_sequence_and_overhangs("aaa" + f + "aaa", -3, -3)
                    )
                    for f in fragments_str
                ]
            products = assembly.in_vivo_assembly(fragments, 7, circular_only=True)
            products_str = [str(p.seq) for p in products]
            assert products_str == expected_outputs


def test_gateway_assembly():

    attB1 = "ACAACTTTGTACAAAAAAGCAGAAG"
    attB2 = "ACAACTTTGTACAAGAAAGCTGGGC"
    attP1 = "AAAATAATGATTTTATTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCTGAACGAGAAGCGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA"
    attP2 = "AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTCTTATAATGCCAACTTTGTACAAGAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGACTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA"
    attR1 = "ACAACTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATGCAGTCACTATG"
    attL1 = "CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTCTTATAATGCCAACTTTGTACAAAAAAGCAGGCT"

    seq1 = Dseqrecord("aaa" + attB1 + "ccc")
    seq2 = Dseqrecord("aaa" + attP1 + "ccc")
    product_BP = "aaaACAACTTTGTACAAAAAAGCTGAACGAGAAGCGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTAccc"
    product_PB = "aaaAAAATAATGATTTTATTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCAGAAGccc"

    products_BP = assembly.gateway_assembly([seq1, seq2], "BP")
    assert [product_BP, product_PB] == [str(s.seq) for s in products_BP]

    seq3 = Dseqrecord("aaa" + attR1 + "ccc")
    seq4 = Dseqrecord("aaa" + attL1 + "ccc")
    product_RL = "aaaACAACTTTGTACAAAAAAGCAGGCTccc"
    product_LR = "aaaCAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTCTTATAATGCCAACTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATGCAGTCACTATGccc"

    products_LR = assembly.gateway_assembly([seq3, seq4], "LR")
    assert [product_RL, product_LR] == [str(s.seq) for s in products_LR]

    # Products are not valid inputs for the parent reaction
    with pytest.raises(ValueError):
        assembly.gateway_assembly(products_BP, "BP")
    with pytest.raises(ValueError):
        assembly.gateway_assembly(products_LR, "LR")

    # Test that greedy is being used (finds)
    with pytest.raises(ValueError) as e:
        assembly.gateway_assembly(products_LR, "LR", greedy=True)
    assert "fragment 2: attB1, attL1, attR1, attP1" in str(e.value)

    # Test multi site only
    seq1 = Dseqrecord("aaa" + attB1 + "ggg" + attB2 + "ccc", circular=True)
    seq2 = Dseqrecord("aaa" + attP1 + "ggg" + attP2 + "ccc", circular=True)

    products = assembly.gateway_assembly([seq1, seq2], "BP")
    assert len(products) == 4
    products = assembly.gateway_assembly([seq1, seq2], "BP", multi_site_only=True)
    assert len(products) == 2

    # Test single input
    seq1 = Dseqrecord("aaa" + attB1 + "ccc" + attP1)
    products = assembly.gateway_assembly([seq1], "BP")
    assert len(products) == 2
    assert (
        str(products[0].seq)
        == "GTACAAAAAAGCAGAAGcccAAAATAATGATTTTATTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTT"
    )
    assert products[0].circular
    assert (
        str(products[1].seq)
        == "aaaACAACTTTGTACAAAAAAGCTGAACGAGAAGCGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA"
    )
    assert not products[1].circular


def test_common_handle_insertion_fragments():

    genome = Dseqrecord("aaa")
    inserts = [Dseqrecord("bbb")]
    with pytest.raises(ValueError) as e:
        assembly.common_handle_insertion_fragments([genome], [])
    assert "Genome must be a Dseqrecord object" in str(e)

    with pytest.raises(ValueError) as e:
        assembly.common_handle_insertion_fragments(genome, [])
    assert "Inserts must be a non-empty list of Dseqrecord objects" in str(e)

    with pytest.raises(ValueError) as e:
        assembly.common_handle_insertion_fragments(genome, inserts[0])
    assert "Inserts must be a list of Dseqrecord objects" in str(e)

    assert [genome] + inserts == assembly.common_handle_insertion_fragments(
        genome, inserts
    )


def test_crispr_integration():
    genome = Dseqrecord(
        "AAaaccggttcaatgcaaacagtaatgatggatgacattcaaagcacTT", name="genome"
    )
    insert = Dseqrecord("aaccggttAAAAAAAAAttcaaagcac", name="insert")
    guide = Primer("ttcaatgcaaacagtaatga", name="guide")
    product, *_ = assembly.crispr_integration(genome, [insert], [guide], 8)
    assert product.seq == Dseq("AAaaccggttAAAAAAAAAttcaaagcacTT")

    # Same result if primer in caps
    guide_caps = Primer(str(guide.seq).upper())
    product, *_ = assembly.crispr_integration(genome, [insert], [guide_caps], 8)
    assert product.seq == Dseq("AAaaccggttAAAAAAAAAttcaaagcacTT")

    # Value error if no guides passed
    with pytest.raises(ValueError):
        assembly.crispr_integration(genome, [insert], [], 8)

    guide2 = Primer("ttcaatgcTTTcagtaatga")
    # Value error if it does not cut
    with pytest.raises(ValueError):
        assembly.crispr_integration(genome, [insert], [guide, guide2], 8)

    # If all guides cut outside, no products
    genome = Dseqrecord(
        "GTGGTCTGAACTCGGTGTTGAGCCCTGGCTGATGTAAaaccggttcaatgcaaacagtaatgatggatgacattcaaagcacTT",
        name="genome",
    )
    guide_outside = Primer("CTGAACTCGGTGTTGAGCCC")
    with pytest.warns(UserWarning):
        products = assembly.crispr_integration(genome, [insert], [guide_outside], 8)
    # Empty products
    assert len(products) == 0

    # If some guides cut outside and others don't, empty list and warning
    with pytest.warns(UserWarning):
        products = assembly.crispr_integration(
            genome, [insert], [guide_outside, guide], 8
        )
    assert len(products) == 0

    # Case where there is offtargets
    homology1 = "AAGTCCGTTCGTTTTACCTG"
    homology2 = "ctgaaatttcgtacaaaaaa"
    guide = Primer("ttcaatgcaaacagtaatga")
    genome = Dseqrecord(
        f"aaattcaatgcaaacagtaatgaTGGaaaaaaaa{homology1}aaattcaatgcaaacagtaatgaTGGaa{homology2}aaaaaa"
    )
    insert = Dseqrecord(f"{homology1}gggg{homology2}", circular=True)

    with pytest.raises(ValueError):
        assembly.crispr_integration(genome, [insert], [guide], 20)

    # Off-targets for the single integrations, but not for the one that uses both homology sites
    genome = Dseqrecord(
        f"aaaaaaaa{homology1}aaattcaatgcaaacagtaatgaTGGaa{homology2}aaaaaa"
    )
    assembly.crispr_integration(genome, [insert], [guide], 20)

    with pytest.warns(UserWarning):
        products = assembly.crispr_integration(genome, [insert], [guide], 20)
    assert len(products) == 1


def test_crispr_integration_edge_case():
    homology1 = "ATGCAAACAGTAATGATGGATGACATTCAAAGCACTGATT"
    guide = Primer("CATTCAAAGCACTGATTaat")
    genome = Dseqrecord(f"aaaaaa{homology1}aattggaa{homology1}tttttttt")

    insert = Dseqrecord(f"{homology1}acaa{homology1}")

    products = assembly.crispr_integration(genome, [insert], [guide], 40)
    assert len(products) == 2


def test_pcr_assembly():

    template = Dseqrecord("TTTTACGTACGTAAAAAAGCGCGCGCTTTTT")
    primer_fwd = Primer("ACGTACGT", name="forward")
    primer_rvs = Primer("GCGCGCGC", name="reverse")

    products = assembly.pcr_assembly(template, primer_fwd, primer_rvs, limit=8)
    assert len(products) == 1
    assert str(products[0].seq) == "ACGTACGTAAAAAAGCGCGCGC"

    # Works with the reverse complement as well, and directionality
    # is determined by forward primer
    products = assembly.pcr_assembly(
        template.reverse_complement(), primer_fwd, primer_rvs, limit=8
    )
    assert len(products) == 1
    assert str(products[0].seq) == "ACGTACGTAAAAAAGCGCGCGC"
    # No annotation
    assert len(products[0].features) == 0

    # By default, this should be false
    assert products[0].source.add_primer_features is False
    # Now with add_primer_features
    product, *_ = assembly.pcr_assembly(
        template, primer_fwd, primer_rvs, limit=8, add_primer_features=True
    )
    assert product.source.add_primer_features is True
    assert len(product.features) == 2
    assert product.features[0].type == "primer_bind"
    assert product.features[0].qualifiers["label"] == [primer_fwd.name]
    assert product.features[1].type == "primer_bind"
    assert product.features[1].qualifiers["label"] == [primer_rvs.name]
    assert str(product.seq) == "ACGTACGTAAAAAAGCGCGCGC"

    # Test with two identical primers
    template = Dseqrecord("TTTTACGTACGTAAAAAAACGTACGTTTTTT")
    primer_fwd = Primer("ACGTACGT")
    products = assembly.pcr_assembly(template, primer_fwd, primer_fwd, limit=8)
    assert len(products) == 1
