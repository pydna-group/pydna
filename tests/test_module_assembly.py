#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa: B950

from pydna import assembly
from pydna.dseqrecord import Dseqrecord
from pydna.parsers import parse
from pydna.utils import eq
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
from pydna.amplify import pcr
from pydna.readers import read
from Bio.Restriction import EcoRV, ZraI
from Bio.Restriction import AjiI, AgeI
from Bio.Restriction import SalI
from pydna.assembly import Assembly
from Bio.Restriction import AatII


def test_built():

    asm = assembly.Assembly(assembly.example_fragments, limit=5)
    lin = asm.assemble_linear()
    crc = asm.assemble_circular()

    assert [l.seq for l in lin] == [l.seq for l in assembly.linear_results]
    assert [c.seq for c in crc] == [c.seq for c in assembly.circular_results]


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
    l = ln0.assemble_linear()[0]

    assert (
        str(l.seq)
        == "ACTACGGCCTTCTCTCCCCCtgtgctgtgctctaTTTTTtattctggctgtatctGGGGGTacgatgctatactgg"
    )

    feature_seqs = (
        [f.extract(a).seq for f in a.features]
        + [f.extract(b).seq for f in b.features]
        + [f.extract(c).seq for f in c.features]
    )

    assembled_feature_seqs = [f.extract(l).seq for f in l.features]
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
            FeatureLocation(ExactPosition(0), ExactPosition(34), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(1), ExactPosition(33), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(0), ExactPosition(20), strand=1), type="misc"
        ),
        SeqFeature(
            FeatureLocation(ExactPosition(20), ExactPosition(34), strand=1), type="misc"
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
    assert circprods[0].seguid() == "cdseguid=CRIbOfddcwCZbvVOOOU4uJYP-So"
    assert circprods[1].seguid() == "cdseguid=CRIbOfddcwCZbvVOOOU4uJYP-So"
    assert circprods[2].seguid() == "cdseguid=zFIq5LWXL_YSxrSF2Q5hbzO0BPw"
    assert circprods[3].seguid() == "cdseguid=zFIq5LWXL_YSxrSF2Q5hbzO0BPw"
    assert (
        str(circprods[0].seq)
        == "acgatgctatactggCCCCCtgtgctgtgctctaTTTTTtattctggctgtatctGGGGGT"
    )
    assert (
        str(circprods[2].seq)
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
    assert eq(correct, candidate, linear=True)

    # text4
    x, y = parse(text4)
    a = assembly.Assembly((x, y), limit=557)

    assert a.assemble_linear()[0].seguid() == "ldseguid=3ay9lXqyTHRr348YDQ4JrtFVx40"
    a = assembly.Assembly((y, x), limit=557)
    assert a.assemble_linear()[0].seguid() == "ldseguid=3ay9lXqyTHRr348YDQ4JrtFVx40"
    a = assembly.Assembly((x.rc(), y), limit=557)
    assert a.assemble_linear()[0].seguid() == "ldseguid=3ay9lXqyTHRr348YDQ4JrtFVx40"
    a = assembly.Assembly((x, y.rc()), limit=557)
    assert a.assemble_linear()[0].seguid() == "ldseguid=3ay9lXqyTHRr348YDQ4JrtFVx40"
    a = assembly.Assembly((x.rc(), y.rc()), limit=557)
    assert a.assemble_linear()[0].seguid() == "ldseguid=3ay9lXqyTHRr348YDQ4JrtFVx40"

    candidate = a.assemble_linear()[0]
    correct = "tcctgacgggtaattttgatttgcatgccgtccgggtgagtcatagcgtctggttgttttgccagattcagcagagtctgtgcaatgcggccgctgaccacatacgatttaggtgacactatagaacgcggccgccagctgaagcttcgtacgctgcaggtcgacggatccccgggttaattaaggcgcgccagatctgtttagcttgccttgtccccgccgggtcacccggccagcgacatggaggcccagaataccctccttgacagtcttgacgtgcgcagctcaggggcatgatgtgactgtcgcccgtacatttagcccatacatccccatgtataatcatttgcatccatacattttgatggccgcacggcgcgaagcaaaaattacggctcctcgctgcagacctgcgagcagggaaacgctcccctcacagacgcgttgaattgtccccacgccgcgcccctgtagagaaatataaaaggttaggatttgccactgaggttcttctttcatatacttccttttaaaatcttgctaggatacagttctcacatcacatccgaacataaacaaccgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgacccaatcacatcacatccgaacataaacaaccatgggtaaaaagcctgaactcaccgcgacgtctgtcgagaagtttctgatcgaaaagttcgacagcgtctccgacctgatgcagctctcggagggcgaagaatctcgtgctttcagcttcgatgtaggagggcgtggatatgtcctgcgggtaaatagctgcgccgatggtttctacaaagatcgttatgtttatcggcactttgcatcggccgcgctcccgattccggaagtgcttgacattggggaattcagcgagagcctgacctattgcatctcccgccgtgcacagggtgtcacgttgcaagacctgcctgaaaccgaactgcccgctgttctgcagccggtcgcggaggccatggatgcgatcgctgcggccgatcttagccagacgagcgggttcggcccattcggaccgcaaggaatcggtcaatacactacatggcgtgatttcatatgcgcgattgctgatccccatgtgtatcactggcaaactgtgatggacgacaccgtcagtgcgtccgtcgcgcaggctctcgatgagctgatgctttgggccgaggactgccccgaagtccggcacctcgtgcacgcggatttcggctccaacaatgtcctgacggacaatggccgcataacagcggtcattgactggagcgaggcgatgttcggggattcccaatacgaggtcgccaacatcttcttctggaggccgtggttggcttgtatggagcagcagacgcgctacttcgagcggaggcatccggagcttgcaggatcgccgcggctccgggcgtatatgctccgcattggtcttgaccaactctatcagagcttggttgacggcaatttcgatgatgcagcttgggcgcagggtcgatgcgacgcaatcgtccgatccggagccgggactgtcgggcgtacacaaatcgcccgcagaagcgcggccgtctggaccgatggctgtgtagaagtactcgccgatagtggaaaccgacgccccagcactcgtccgagggcaaaggaataatcagtactgacaataaaaagattcttgtagggataacagggtaatcggagtgccatctgtgcagacaaacgcatcaggatagagtcttttgtaacgaccccgtctccaccaacttggtatgcttgaaatctcaaggccattacacattcagttatgtgaacgaaaggtctttatttaacgtagcataaactaaataatacaggttccggttagcctgcaatgtgttaaatctaaaggagcatacccaaaatgaactgaagacaaggaaatttgcttgtccagatgtgattgagcatttgaacgttaataacataacatttttatacttaactatagaaagacttgtataaaaactggcaaacgagatattctgaatattggtgcatatttcaggtagaaaagcttacaaaacaatctaatcataatattgagatgaagagaaagataaaagaaaaaacgataagtcagatgagattatgattgtactttgaaatcgaggaacaaagtatatacggtagtagttccccgagttataacgggagatcatgtaaattgagaaaccagataaagatttggtatgcactctagcaagaaaataaaatgatgaatctatgatatagatcacttttgttccagcgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgacgatcagatctttcaggaaagtttcggaggagatagtgttcggcagtttgtacatcatctgcgggatcaggtacggtttgatcaggttgtagaagatcaggtaagacatagaatcgatgtagatgatcggtttgtttttgttgatttttacgtaacagttcagttggaatttgttacgcagacccttaaccaggtattctacttcttcgaaagtgaaagactgggtgttcagtacgatcgatttgttggtagagtttttgttgtaatcccatttaccaccatcatccatgaaccagtatgccagagacatcggggtcaggtagttttcaaccaggttgttcgggatggtttttttgttgttaacgatgaacaggttagccagtttgttgaaagcttggtgtttgaaagtctgggcgccccaggtgattaccaggttacccaggtggttaacacgttcttttttgtgcggcggggacagtacccactgatcgtacagcagacatacgtggtccatgtatgctttgtttttccactcgaactgcatacagtaggttttaccttcatcacgagaacggatgtaagcatcacccaggatcagaccgatacctgcttcgaactgttcgatgttcagttcgatcagctgggatttgtattctttcagcagtttagagttcggacccaggttcattacctggttttttttgatgtttttcatatgcatggatccggggttttttctccttgacgttaaagtatagaggtatattaacaattttttgttgatacttttattacatttgaataagaagtaatacaaaccgaaaatgttgaaagtattagttaaagtggttatgcagtttttgcatttatatatctgttaatagatcaaaaatcatcgcttcgctgattaattaccccagaaataaggctaaaaaactaatcgcattatcatcctatggttgttaatttgattcgttcatttgaaggtttgtggggccaggttactgccaatttttcctcttcataaccataaaagctagtattgtagaatctttattgttcggagcagtgcggcgcgaggcacatctgcgtttcaggaacgcgaccggtgaagacgaggacgcacggaggagagtcttccttcggagggctgtcacccgctcggcggcttctaatccgtacttcaatatagcaatgagcagttaagcgtattactgaaagttccaaagagaaggtttttttaggctaagataatggggctctttacatttccacaacatataagtaagattagatatggatatgtatatggatatgtatatggtggtaatgccatgtaatatgattattaaacttctttgcgtccatccaacgagatctggcgcgccttaattaacccaacctgcattaatgaatcggccaacgcgcggattaccctgttatccctacatattgttgtgccatctgtgcagacaaacgcatcaggattcagtactgacaataaaaagattcttgttttcaagaacttgtcatttgtatagtttttttatattgtagttgttctattttaatcaaatgttagcgtgatttatattttttttcgcctcgacatcatctgcccagatgcgaagttaagtgcgcagaaagtaatatcatgcgtcaatcgtatgtgaatgctggtcgctatactgctgtcgattcgatactaacgccgccatccagtgtcgaaaacgagctcgaattcatcgatgatatcagatccactagtggcctatgcggccgcggatctgccggtctccctatagtgagtcgatccggatttacctgaatcaattggcgaaattttttgtacgaaatttcagccacttcacag"
    assert eq(correct, candidate, linear=True)

    # text5
    list_of_formatted_seq_records = parse(text5)
    a = assembly.Assembly(list_of_formatted_seq_records, limit=60)
    candidate = a.assemble_circular()[0]
    correct = "tcgcgcgtttcggtgatgacggtgaaaacctctgacacatgcagctcccggagacggtcacagcttgtctgtaagcggatgccgggagcagacaagcccgtcagggcgcgtcagcgggtgttggcgggtgtcggggctggcttaactatgcggcatcagagcagattgtactgagagtgcaccatagatcctgaggatcggggtgataaatcagtctgcgccacatcgggggaaacaaaatggcgcgagatctaaaaaaaaaggctccaaaaggagcctttcgcgctaccaggtaacgcgccactccgacgggattaacgagtgccgtaaacgacgatggttttaccgtgtgcggagatcaggttctgatcctcgagcatcttaagaattcgtcccacggtttgtctagagcagccgacaatctggccaatttcctgacgggtaattttgatttgcatgccgtccgggtgagtcatagcgtctggttgttttgccagattcagcagagtctgtgcaatgcggccgctgaccacatacgatttaggtgacactatagaacgcggccgccagctgaagcttcgtacgctgcaggtcgacggatccccgggttaattaaggcgcgccagatctgtttagcttgccttgtccccgccgggtcacccggccagcgacatggaggcccagaataccctccttgacagtcttgacgtgcgcagctcaggggcatgatgtgactgtcgcccgtacatttagcccatacatccccatgtataatcatttgcatccatacattttgatggccgcacggcgcgaagcaaaaattacggctcctcgctgcagacctgcgagcagggaaacgctcccctcacagacgcgttgaattgtccccacgccgcgcccctgtagagaaatataaaaggttaggatttgccactgaggttcttctttcatatacttccttttaaaatcttgctaggatacagttctcacatcacatccgaacataaacaaccgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgacccaatcacatcacatccgaacataaacaaccatgggtaaaaagcctgaactcaccgcgacgtctgtcgagaagtttctgatcgaaaagttcgacagcgtctccgacctgatgcagctctcggagggcgaagaatctcgtgctttcagcttcgatgtaggagggcgtggatatgtcctgcgggtaaatagctgcgccgatggtttctacaaagatcgttatgtttatcggcactttgcatcggccgcgctcccgattccggaagtgcttgacattggggaattcagcgagagcctgacctattgcatctcccgccgtgcacagggtgtcacgttgcaagacctgcctgaaaccgaactgcccgctgttctgcagccggtcgcggaggccatggatgcgatcgctgcggccgatcttagccagacgagcgggttcggcccattcggaccgcaaggaatcggtcaatacactacatggcgtgatttcatatgcgcgattgctgatccccatgtgtatcactggcaaactgtgatggacgacaccgtcagtgcgtccgtcgcgcaggctctcgatgagctgatgctttgggccgaggactgccccgaagtccggcacctcgtgcacgcggatttcggctccaacaatgtcctgacggacaatggccgcataacagcggtcattgactggagcgaggcgatgttcggggattcccaatacgaggtcgccaacatcttcttctggaggccgtggttggcttgtatggagcagcagacgcgctacttcgagcggaggcatccggagcttgcaggatcgccgcggctccgggcgtatatgctccgcattggtcttgaccaactctatcagagcttggttgacggcaatttcgatgatgcagcttgggcgcagggtcgatgcgacgcaatcgtccgatccggagccgggactgtcgggcgtacacaaatcgcccgcagaagcgcggccgtctggaccgatggctgtgtagaagtactcgccgatagtggaaaccgacgccccagcactcgtccgagggcaaaggaataatcagtactgacaataaaaagattcttgtagggataacagggtaatcggagtgccatctgtgcagacaaacgcatcaggatagagtcttttgtaacgaccccgtctccaccaacttggtatgcttgaaatctcaaggccattacacattcagttatgtgaacgaaaggtctttatttaacgtagcataaactaaataatacaggttccggttagcctgcaatgtgttaaatctaaaggagcatacccaaaatgaactgaagacaaggaaatttgcttgtccagatgtgattgagcatttgaacgttaataacataacatttttatacttaactatagaaagacttgtataaaaactggcaaacgagatattctgaatattggtgcatatttcaggtagaaaagcttacaaaacaatctaatcataatattgagatgaagagaaagataaaagaaaaaacgataagtcagatgagattatgattgtactttgaaatcgaggaacaaagtatatacggtagtagttccccgagttataacgggagatcatgtaaattgagaaaccagataaagatttggtatgcactctagcaagaaaataaaatgatgaatctatgatatagatcacttttgttccagcgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgacgatcagatctttcaggaaagtttcggaggagatagtgttcggcagtttgtacatcatctgcgggatcaggtacggtttgatcaggttgtagaagatcaggtaagacatagaatcgatgtagatgatcggtttgtttttgttgatttttacgtaacagttcagttggaatttgttacgcagacccttaaccaggtattctacttcttcgaaagtgaaagactgggtgttcagtacgatcgatttgttggtagagtttttgttgtaatcccatttaccaccatcatccatgaaccagtatgccagagacatcggggtcaggtagttttcaaccaggttgttcgggatggtttttttgttgttaacgatgaacaggttagccagtttgttgaaagcttggtgtttgaaagtctgggcgccccaggtgattaccaggttacccaggtggttaacacgttcttttttgtgcggcggggacagtacccactgatcgtacagcagacatacgtggtccatgtatgctttgtttttccactcgaactgcatacagtaggttttaccttcatcacgagaacggatgtaagcatcacccaggatcagaccgatacctgcttcgaactgttcgatgttcagttcgatcagctgggatttgtattctttcagcagtttagagttcggacccaggttcattacctggttttttttgatgtttttcatatgcatggatccggggttttttctccttgacgttaaagtatagaggtatattaacaattttttgttgatacttttattacatttgaataagaagtaatacaaaccgaaaatgttgaaagtattagttaaagtggttatgcagtttttgcatttatatatctgttaatagatcaaaaatcatcgcttcgctgattaattaccccagaaataaggctaaaaaactaatcgcattatcatcctatggttgttaatttgattcgttcatttgaaggtttgtggggccaggttactgccaatttttcctcttcataaccataaaagctagtattgtagaatctttattgttcggagcagtgcggcgcgaggcacatctgcgtttcaggaacgcgaccggtgaagacgaggacgcacggaggagagtcttccttcggagggctgtcacccgctcggcggcttctaatccgtacttcaatatagcaatgagcagttaagcgtattactgaaagttccaaagagaaggtttttttaggctaagataatggggctctttacatttccacaacatataagtaagattagatatggatatgtatatggatatgtatatggtggtaatgccatgtaatatgattattaaacttctttgcgtccatccaacgagatctggcgcgccttaattaacccaacctgcattaatgaatcggccaacgcgcggattaccctgttatccctacatattgttgtgccatctgtgcagacaaacgcatcaggattcagtactgacaataaaaagattcttgttttcaagaacttgtcatttgtatagtttttttatattgtagttgttctattttaatcaaatgttagcgtgatttatattttttttcgcctcgacatcatctgcccagatgcgaagttaagtgcgcagaaagtaatatcatgcgtcaatcgtatgtgaatgctggtcgctatactgctgtcgattcgatactaacgccgccatccagtgtcgaaaacgagctcgaattcatcgatgatatcagatccactagtggcctatgcggccgcggatctgccggtctccctatagtgagtcgatccggatttacctgaatcaattggcgaaattttttgtacgaaatttcagccacttcacaggcggttttcgcacgtacccatgcgctacgttcctggccctcttcaaacaggcccagttcgccaataaaatcaccctgattcagataggagaggatcatttctttaccctcttcgtctttgatcagcactgccacagagcctttaacgatgtagtacagcgtttccgctttttcaccctggtgaataagcgtgctcttggatgggtacttatgaatgtggcaatgagacaagaaccattcgagagtaggatccgtttgaggtttaccaagtaccataagatccttaaatttttattatctagctagatgataatattatatcaagaattgtacctgaaagcaaataaattttttatctggcttaactatgcggcatcagagcagattgtactgagagtgcaccatatgcggtgtgaaataccgcacagatgcgtaaggagaaaataccgcatcaggcgctcttccgcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatacggttatccacagaatcaggggataacgcaggaaagaacatgtgagcaaaaggccagcaaaagcccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaaggacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgaggggtaataactgatataattaaattgaagctctaatttgtgagtttagtatacatgcatttacttataatacagttttttagttttgctggccgcatcttctcaaatatgcttcccagcctgcttttctgtaacgttcaccctctaccttagcatcccttccctttgcaaatagtcctcttccaacaataataatgtcagatcctgtagagaccacatcatccacggttctatactgttgacccaatgcgtctcccttgtcatctaaacccacaccgggtgtcataatcaaccaatcgtaaccttcatctcttccacccatgtctctttgagcaataaagccgataacaaaatctttgtcgctcttcgcaatgtcaacagtacccttagtatattctccagtagatagggagcccttgcatgacaattctgctaacatcaaaaggcctctaggttcctttgttacttcttctgccgcctgcttcaaaccgctaacaatacctgggcccaccacaccgtgtgcattcgtaatgtctgcccattctgctattctgtatacacccgcagagtactgcaatttgactgtattaccaatgtcagcaaattttctgtcttcgaagagtaaaaaattgtacttggcggataatgcctttagcggcttaactgtgccctccatggaaaaatcagtcaaaatatccacatgtgtttttagtaaacaaattttgggacctaatgcttcaactaactccagtaattccttggtggtacgaacatccaatgaagcacacaagtttgtttgcttttcgtgcatgatattaaatagcttggcagcaacaggactaggatgagtagcagcacgttccttatatgtagctttcgacatgatttatcttcgtttcctgcaggtttttgttctgtgcagttgggttaagaatactgggcaatttcatgtttcttcaacactacatatgcgtatatataccaatctaagtctgtgctccttccttcgttcttccttctgttcggagattaccgaatcaaaaaaatttcaaagaaaccgaaatcaaaaaaaagaataaaaaaaaaatgatgaattgaattgaaaagctagcttatcgatgataagctgtcaaagatgagaattaattccacggactatagactatactagatactccgtctactgtacgatacacttccgctcaggtccttgtcctttaacgaggccttaccactcttttgttactctattgatccagctcagcaaaggcagtgtgatctaagattctatcttcgcgatgtagtaaaactagctagaccgagaaagagactagaaatgcaaaaggcacttctacaatggctgccatcattattatccgatgtgacgctgcagcttctcaatgatattcgaatacgctttgaggagatacagcctaatatccgacaaactgttttacagatttacgatcgtacttgttacccatcattgaattttgaacatccgaacctgggagttttccctgaaacagatagtatatttgaacctgtataataatatatagtctagcgctttacggaagacaatgtatgtatttcggttcctggagaaactattgcatctattgcataggtaatcttgcacgtcgcatccccggttcattttctgcgtttccatcttgcacttcaatagcatatctttgttaacgaagcatctgtgcttcattttgtagaacaaaaatgcaacgcgagagcgctaatttttcaaacaaagaatctgagctgcatttttacagaacagaaatgcaacgcgaaagcgctattttaccaacgaagaatctgtgcttcatttttgtaaaacaaaaatgcaacgcgacgagagcgctaatttttcaaacaaagaatctgagctgcatttttacagaacagaaatgcaacgcgagagcgctattttaccaacaaagaatctatacttcttttttgttctacaaaaatgcatcccgagagcgctatttttctaacaaagcatcttagattactttttttctcctttgtgcgctctataatgcagtctcttgataactttttgcactgtaggtccgttaaggttagaagaaggctactttggtgtctattttctcttccataaaaaaagcctgactccacttcccgcgtttactgattactagcgaagctgcgggtgcattttttcaagataaaggcatccccgattatattctataccgatgtggattgcgcatactttgtgaacagaaagtgatagcgttgatgattcttcattggtcagaaaattatgaacggtttcttctattttgtctctatatactacgtataggaaatgtttacattttcgtattgttttcgattcactctatgaatagttcttactacaatttttttgtctaaagagtaatactagagataaacataaaaaatgtagaggtcgagtttagatgcaagttcaaggagcgaaaggtggatgggtaggttatatagggatatagcacagagatatatagcaaagagatacttttgagcaatgtttgtggaagcggtattcgcaatgggaagctccaccccggttgataatcagaaaagccccaaaaacaggaagattattatcaaaaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacagttaccaatgcttaatcagtgaggcacctatctcagcgatctgtctatttcgttcatccatagttgcctgactccccgtcgtgtagataactacgatacgggagcgcttaccatctggccccagtgctgcaatgataccgcgagacccacgctcaccggctccagatttatcagcaataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttggcattgctacaggcatcgtggtgtcactctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggataatagtgtatcacatagcagaactttaaaagtgctcatcattggaaaacgttcttcggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatgtaacccactcgtgcacccaactgatcttcagcatcttttactttcaccagcgtttctgggtgagcaaaaacaggaaggcaaaatgccgcaaaaaagggaataagggcgacacggaaatgttgaatactcatactcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaaaagtgccacctgctaagaaaccattattatcatgacattaacctataaaaataggcgtatcacgaggccctttcgtc"
    assert len(correct) == 9772
    assert eq(correct, candidate, circular=True)

    # text6

    y, x = parse(text6)
    a = assembly.Assembly((x, y, x), limit=40)
    candidate = a.assemble_linear()[0]
    correct = "aaattgacctcgtttgattgctctcaggtttccaacaacttgtttgttaccggtttcagtacgggtattatcaagttatgggatgccagagcggctgaagctgccaccactgatttgacttatagacaaaacggtgaagatccaattcagaatgaaatcgctaacttctatcatgccggcggtgattctgttgtcgacgtccaattttctgccacttcaagctctgaattctttactgtgggcggcactggtaatatttaccactggaacactgactattccttgtcaaagtacaatcctgacgataccattgctcctcctcaagatgccactgaagaatcacaaacaaaatctctgagattcttgcacaagggtggaagcagaagatctccaaaacagattggaagaagaaacaccgccgcttggcatccagttattgaaaacttggttggtactgtcgatgacgatagtttagtcagcatctacaagccctataccgaggaaagcgaatagtcgccatgctaaacgcgcggaacaaggccatatttatatatttaatgcttttaactattattagtttttctcacccagggctggtttctttaccctgtgtgatgaaagtgcgcgcctaaagttatgtatgagctactgttagtatatttattattcaatatttaattaaacaatacttactactgatgaaagataccgtacacctgctggcgaataagtcactatacacgagcaattatatttctcaggtcattgcgacaattgaaaaatcggatcgactgattatcattaaggcgtgtgacggctcagaataaggaagttctattttaaagggcagcaaaaatttggcacaatttaaaaagaaagttggtctctggtagggttggatccggtttttcataacttgcgtttcatttatggttgcgtatttctcgcagtgtacatagaccaaattaactctatccttccattgcacaatttgcccagtatggatgtccacgaggtctctgctatggacacatttacgcccgtacgctgcaggtcgacggatccccgggttaattaaggcgcgccagatctgtttagcttgcctcgtccccgccgggtcacccggccagcgacatggaggcccagaataccctccttgacagtcttgacgtgcgcagctcaggggcatgatgtgactgtcgcccgtacatttagcccatacatccccatgtataatcatttgcatccatacattttgatggccgcacggcgcgaagcaaaaattacggctcctcgctgcagacctgcgagcagggaaacgctcccctcacagacgcgttgaattgtccccacgccgcgcccctgtagagaaatataaaaggttaggatttgccactgaggttcttctttcatatacttccttttaaaatcttgctaggatacagttctcacatcacatccgaacataaacaaccatgggtaaggaaaagactcacgtttcgaggccgcgattaaattccaacatggatgctgatttatatgggtataaatgggctcgcgataatgtcgggcaatcaggtgcgacaatctatcgattgtatgggaagcccgatgcgccagagttgtttctgaaacatggcaaaggtagcgttgccaatgatgttacagatgagatggtcagactaaactggctgacggaatttatgcctcttccgaccatcaagcattttatccgtactcctgatgatgcatggttactcaccactgcgatccccggcaaaacagcattccaggtattagaagaatatcctgattcaggtgaaaatattgttgatgcgctggcagtgttcctgcgccggttgcattcgattcctgtttgtaattgtccttttaacagcgatcgcgtatttcgtctcgctcaggcgcaatcacgaatgaataacggtttggttgatgcgagtgattttgatgacgagcgtaatggctggcctgttgaacaagtctggaaagaaatgcataagcttttgccattctcaccggattcagtcgtcactcatggtgatttctcacttgataaccttatttttgacgaggggaaattaataggttgtattgatgttggacgagtcggaatcgcagaccgataccaggatcttgccatcctatggaactgcctcggtgagttttctccttcattacagaaacggctttttcaaaaatatggtattgataatcctgatatgaataaattgcagtttcatttgatgctcgatgagtttttctaatcagtactgacaataaaaagattcttgttttcaagaacttgtcatttgtatagtttttttatattgtagttgttctattttaatcaaatgttagcgtgatttatattttttttcgcctcgacatcatctgcccagatgcgaagttaagtgcgcagaaagtaatatcatgcgtcaatcgtatgtgaatgctggtcgctatactgctgtcgattcgatactaacgccgccatccagtgtcgaaaacgagctcgaattcatcgatgaatgacgttatagacgagcctacgagaccgacaccgtaatcaataagcaattttttgcgaacatgattaagtttattttactttatttcgcattttccataaaaaaattttcctccatacaagattcatacccaggatagcctaactaaaatgtatggctttatacagcttgacgacaaacttaatttgaatatatagatatattaatttaatcagcaggattagcatttaatgtcttacatgtctttgtatttggcatacgtattcagtgatattttaatacctcgcaacgcaatttccagcgagtttacgtttgaagcaaacgtaactcccatgagttatacataacgcctgctgcagctgaaccagaacaatataccgtggtattaccgtaagcatttgaaaacaataataagtttccaagtgccgctttagatccaaagtacaaggacatattgtggaagagtaagttgaatatttgaaaatgagagctttttcagcagccaccgttagggccacaactaggaagtcgttcatcccaatggcaccaagaactccttttgtgactccatcatttacaaagaatgtaggctcaatgagaagaatgagattttattctgatgaagccaaaagtgaagaatccaaagaaaacaatgaagatttgactgaagagcaatcagaaatcaagaaattagagagccagttaagcgcgaagactaaagaagcttctgaactcaaggacagattattaagatctgtggcagatttcagaaatttacaacaagtcacaaagaaggatattcagaaagctaaggactttgctttacagaagtttgcaaaggatttattggaatctgtagataactttggtcatgctttgaatgcttttaaagaggaagacttacaaaagtccaaggaaattagtgatttgtatacaggggttagaatgacaagagatgtttttgaaaacaccctaagaaagcacggtattgaaaaattagacccattgggagaaccatttgatccaaataaacacgaa"
    assert len(correct) == 3587
    assert eq(correct, candidate, linear=True)

    # text7 cSEGUID G2drVQAIaRfFXBfqEc5Kddac36A
    list_of_formatted_seq_records = parse(text7)
    a = assembly.Assembly(list_of_formatted_seq_records, limit=28)

    candidate = a.assemble_circular()[0]
    correct = "aattggccagattgtcggctgctctagacaaaccgtgggacgaattcttaagatgctcgaggatcagaacctgatctccgcacacggtaaaaccatcgtcgtttacggcactcgttaatcccgtcggagtggcgcgttacctggtagcgcgaaaggctccttttggagcctttttttttagatctcgcgccattttgtttcccccgatgtggcgcagactgatttatcaccccgatcctcaggatctatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtgatacactattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagagtgacaccacgatgcctgtagcaatgccaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagcgctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataataatcttcctgtttttggggcttttctgattatcaaccggggtggagcttcccattgcgaataccgcttccacaaacattgctcaaaagtatctctttgctatatatctctgtgctatatccctatataacctacccatccacctttcgctccttgaacttgcatctaaactcgacctctacattttttatgtttatctctagtattactctttagacaaaaaaattgtagtaagaactattcatagagtgaatcgaaaacaatacgaaaatgtaaacatttcctatacgtagtatatagagacaaaatagaagaaaccgttcataattttctgaccaatgaagaatcatcaacgctatcactttctgttcacaaagtatgcgcaatccacatcggtatagaatataatcggggatgcctttatcttgaaaaaatgcacccgcagcttcgctagtaatcagtaaacgcgggaagtggagtcaggctttttttatggaagagaaaatagacaccaaagtagccttcttctaaccttaacggacctacagtgcaaaaagttatcaagagactgcattatagagcgcacaaaggagaaaaaaagtaatctaagatgctttgttagaaaaatagcgctctcgggatgcatttttgtagaacaaaaaagaagtatagattctttgttggtaaaatagcgctctcgcgttgcatttctgttctgtaaaaatgcagctcagattctttgtttgaaaaattagcgctctcgtcgcgttgcatttttgttttacaaaaatgaagcacagattcttcgttggtaaaatagcgctttcgcgttgcatttctgttctgtaaaaatgcagctcagattctttgtttgaaaaattagcgctctcgcgttgcatttttgttctacaaaatgaagcacagatgcttcgttaacaaagatatgctattgaagtgcaagatggaaacgcagaaaatgaaccggggatgcgacgtgcaagattacctatgcaatagatgcaatagtttctccaggaaccgaaatacatacattgtcttccgtaaagcgctagactatatattattatacaggttcaaatatactatctgtttcagggaaaactcccaggttcggatgttcaaaattcaatgatgggtaacaagtacgatcgtaaatctgtaaaacagtttgtcggatattaggctgtatctcctcaaagcgtattcgaatatcattgagaagctgcagcgtcacatcggataataatgatggcagccattgtagaagtgccttttgcatttctagtctctttctcggtctagctagttttactacatcgcgaagatagaatcttagatcacactgcctttgctgagctggatcaatagagtaacaaaagagtggtaaggcctcgttaaaggacaaggacctgagcggaagtgtatcgtacagtagacggagtatctagtatagtctatagtccgtggaattaattctcatctttgacagcttatcatcgataagctagcttttcaattcaattcatcattttttttttattcttttttttgatttcggtttctttgaaatttttttgattcggtaatctccgaacagaaggaagaacgaaggaaggagcacagacttagattggtatatatacgcatatgtagtgttgaagaaacatgaaattgcccagtattcttaacccaactgcacagaacaaaaacctgcaggaaacgaagataaatcatgtcgaaagctacatataaggaacgtgctgctactcatcctagtcctgttgctgccaagctatttaatatcatgcacgaaaagcaaacaaacttgtgtgcttcattggatgttcgtaccaccaaggaattactggagttagttgaagcattaggtcccaaaatttgtttactaaaaacacatgtggatattttgactgatttttccatggagggcacagttaagccgctaaaggcattatccgccaagtacaattttttactcttcgaagacagaaaatttgctgacattggtaatacagtcaaattgcagtactctgcgggtgtatacagaatagcagaatgggcagacattacgaatgcacacggtgtggtgggcccaggtattgttagcggtttgaagcaggcggcagaagaagtaacaaaggaacctagaggccttttgatgttagcagaattgtcatgcaagggctccctatctactggagaatatactaagggtactgttgacattgcgaagagcgacaaagattttgttatcggctttattgctcaaagagacatgggtggaagagatgaaggttacgattggttgattatgacacccggtgtgggtttagatgacaagggagacgcattgggtcaacagtatagaaccgtggatgatgtggtctctacaggatctgacattattattgttggaagaggactatttgcaaagggaagggatgctaaggtagagggtgaacgttacagaaaagcaggctgggaagcatatttgagaagatgcggccagcaaaactaaaaaactgtattataagtaaatgcatgtatactaaactcacaaattagagcttcaatttaattatatcagttattacccctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagataaaaaatttatttgctttcaggtacaattcttgatataatattatcatctagctagataataaaaatttaaggatcttatggtacttggtaaacctcaaacggatcctactctcgaatggttcttgtctcattgccacattcataagtacccatccaagagcacgcttattcaccagggtgaaaaagcggaaacgctgtactacatcgttaaaggctctgtggcagtgctgatcaaagacgaagagggtaaagaaatgatcctctcctatctgaatcagggtgattttattggcgaactgggcctgtttgaagagggccaggaacgtagcgcatgggtacgtgcgaaaaccgcctgtgaagtggctgaaatttcgtacaaaaaatttcgccaattgattcaggtaaatccggatgctggaacaaaagtgatctatatcatagattcatcattttattttcttgctagagtgcataccaaatctttatctggtttctcaatttacatgatctcccgttataactcggggaactactaccgtatatactttgttcctcgatttcaaagtacaatcataatctcatctgacttatcgttttttcttttatctttctcttcatctcaatattatgattagattgttttgtaagcttttctacctgaaatatgcaccaatattcagaatatctcgtttgccagtttttatacaagtctttctatagttaagtataaaaatgttatgttattaacgttcaaatgctcaatcacatctggacaagcaaatttccttgtcttcagttcattttgggtatgctcctttagatttaacacattgcaggctaaccggaacctgtattatttagtttatgctacgttaaataaagacctttcgttcacataactgaatgtgtaatggccttgagatttcaagcataccaagttggtggagacggggtcgttacaaaagactctatcctgatgcgtttgtctgcacagatggcgcgtccgattaccctgttatccctacaagaatctttttattgtcagtactgattattcctttgccctcggacgagtgctggggcgtcggtttccactatcggcgagtacttctacacagccatcggtccagacggccgcgcttctgcgggcgatttgtgtacgcccgacagtcccggctccggatcggacgattgcgtcgcatcgaccctgcgcccaagctgcatcatcgaaattgccgtcaaccaagctctgatagagttggtcaagaccaatgcggagcatatacgcccggagccgcggcgatcctgcaagctccggatgcctccgctcgaagtagcgcgtctgctgctccatacaagccaaccacggcctccagaagaagatgttggcgacctcgtattgggaatccccgaacatcgcctcgctccagtcaatgaccgctgttatgcggccattgtccgtcaggacattgttggagccgaaatccgcgtgcacgaggtgccggacttcggggcagtcctcggcccaaagcatcagctcatcgagagcctgcgcgacggacgcactgacggtgtcgtccatcacagtttgccagtgatacacatggggatcagcaatcgcgcatatgaaatcacgccatgtagtgtattgaccgattccttgcggtccgaatgggccgaacccgctcgtctggctaagatcggccgcagcgatcgcatccatggcctccgcgaccggctgcagaacagcgggcagttcggtttcaggcaggtcttgcaacgtgacaccctgtgcacggcgggagatgcaataggtcaggctctcgctgaattccccaatgtcaagcacttccggaatcgggagcgcggccgatgcaaagtgccgataaacataacgatctttgtagaaaccatcggcgcagctatttacccgcaggacatatccacgccctcctacatcgaagctgaaagcacgagattcttcgccctccgagagctgcatcaggtcggagacgctgtcgaacttttcgatcagaaacttctcgacagacgtcgcggtgagttcaggctttttacccatggttgtttatgttcggatgtgatgtgattggccggctgcaggtcactagtgagaaagtgggcaacctggcgttcctcgacatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacgtcagcggccgcattgcacagactctgctgaatctggcaaaacaaccagacgctatgactcacccggacggcatgcaaatcaaaattacccgtcagga"
    assert len(correct) == 7911
    assert eq(correct, candidate, circular=True)

    assert repr(candidate) == "Contig(o7911)"

    h = """\
 -|5681|61
|       \\/
|       /\\
|       61|650_rc|32
|                 \\/
|                 /\\
|                 32|1196_rc|49
|                            \\/
|                            /\\
|                            49|630_rc|98
|                                      \\/
|                                      /\\
|                                      98-
|                                         |
 -----------------------------------------"""

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
 -----------------------------------"""

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

    # text8 cSEGUID wM7nM6oJer3bB6RV81IH78e02j4 7911bp
    list_of_formatted_seq_records = parse(text8)
    a = assembly.Assembly(list_of_formatted_seq_records, limit=28)
    candidate = a.assemble_circular()[0]

    assert candidate.seguid() == "cdseguid=-mVwekticpAYIT9C4JcXmOGFkRo"


def test_MXblaster1():
    """test MXblaster1"""

    primer = parse("primers.fas", ds=False)
    primer = primer[::-1]
    primer = primer[37:]

    for i, p in enumerate(primer):
        assert int(p.id.split("_")[0]) == i

    """ These are PCRs to get the genes and the terminator-promoters """
    AgTEFp = pcr(primer[524], primer[523], read("pAG25.gb"))
    hph = pcr(primer[502], primer[501], read("pAG32.gb"))
    KlLEU2tt = pcr(primer[520], primer[519], read("KlLEU2tt.gb"))

    """ The Gal1 promoter-ISceI fragment is made in two steps """
    gal1_ISceI_1 = pcr(primer[234], primer[316], read("pGSHU 7180bp.gb"))
    gal1_ISceI_2 = pcr(primer[562], primer[234], gal1_ISceI_1)
    AgTEFt = pcr(primer[522], primer[521], read("pAG25.gb"))

    """ load pCAPs and pCAPs-pSU0 sequences as Dseqrecord objects """
    pCAPs = read("pCAPs-AjiI.gb")
    pCAPs_pSU0 = read("pCAPs-pSU0.gb")

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

    AX023560 = read("AX023560.gb")

    GAL10prom_slice = slice(
        AX023560.features[1].location.start, AX023560.features[1].location.end
    )

    GAL10prom = AX023560[GAL10prom_slice]

    assert GAL10prom.seq == AX023560.features[1].extract(AX023560).seq

    GIN11M86 = read("GIN11M86.gb")

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

    candidates = a.assemble_circular()

    candidate = candidates[1]
    assert len(candidate) == 10566

    assert candidate.seguid() == "cdseguid=c48cBUb3wF-Sdhzh0Tlprp-0CEg"

    pCAPs_MX4blaster2 = candidate

    pCAPs_MX4blaster2 = pCAPs_MX4blaster2.synced("tcgcgcgtttcggtgatgacggtgaaaacc")

    assert len(pCAPs_MX4blaster2) == 10566
    pCAPs_MX4blaster2_old = read("pMX4blaster2_old.gb")

    assert len(pCAPs_MX4blaster2_old) == 10566
    assert pCAPs_MX4blaster2_old.seguid() == "cdseguid=c48cBUb3wF-Sdhzh0Tlprp-0CEg"
    assert eq(pCAPs_MX4blaster2, pCAPs_MX4blaster2_old)
    assert pCAPs_MX4blaster2.seguid() == "cdseguid=c48cBUb3wF-Sdhzh0Tlprp-0CEg"


def test_assemble_pGUP1():

    GUP1rec1sens = read("GUP1rec1sens.txt")
    GUP1rec2AS = read("GUP1rec2AS.txt")
    GUP1_locus = read("GUP1_locus.gb")
    pGREG505 = read("pGREG505.gb")

    insert = pcr(GUP1rec1sens, GUP1rec2AS, GUP1_locus)

    his3, lin_vect = pGREG505.cut(SalI)

    ab = assembly.Assembly([insert, lin_vect], limit=28)

    pGUP1 = ab.assemble_circular()[0]

    pGUP1 = pGUP1.synced(pGREG505.seq[:50])

    pGUP1_correct = read("pGUP1_correct.gb")

    assert len(pGUP1_correct) == 9981
    assert len(pGUP1) == 9981
    assert eq(pGUP1, pGUP1_correct)
    assert pGUP1_correct.seguid() == "cdseguid=QiK2pH9yioTPfSobUTLz4CPiNzY"
    assert pGUP1.seguid() == "cdseguid=QiK2pH9yioTPfSobUTLz4CPiNzY"


# def test_35_36():
#    import sys
#    from pydna.assembly import _od
#    if sys.version_info < (3, 6):
#        from collections import OrderedDict
#        assert _od==OrderedDict
#    else:
#        assert _od==dict


def test_pYPK7_TDH3_GAL2_PGI1():

    pMEC1142 = read("pYPK0_TDH3_GAL2_PGI1.gb")

    pYPKp7 = read("pYPKp7.gb")

    pYPKp7_AatII = pYPKp7.linearize(AatII)

    z = Assembly((pYPKp7_AatII, pMEC1142), limit=300)

    assert z.assemble_circular()[1].seguid() == "cdseguid=DeflrptvvS6m532WogvxQSgVKpk"


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

    pAG32 = read("pAG32.gb")
    pMEC1135 = read("pMEC1135.gb")

    hygromycin_product = pcr(f, r, pAG32)

    asm_hyg = Assembly((pMEC1135, hygromycin_product, pMEC1135))

    candidate, other = asm_hyg.assemble_linear()

    # AmpR feature
    assert (
        pMEC1135.features[-1].extract(pMEC1135).seq
        == candidate.features[-1].extract(candidate).seq
    )


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
    asm = Assembly((a, b, c), limit=14)
    x = asm.assemble_linear()[0]
    # print(x.features)
    # print(x)
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
            print(feat.qualifiers["label"])
            print(feat.extract(x).seq, "extracted feat")
            print(feature_sequences[feat.qualifiers["label"]].seq, "original feat")
            assert (
                feat.extract(x).seq == feature_sequences[feat.qualifiers["label"]].seq
            )


# acgatgctatactgtgCCNCCtgtgctgtgctcta
#                      TGTGCTGTGCTCTA
#                      tgtgctgtgctctaTTTTTTTtattctggctgtatcCCCCCC
#                                           TATTCTGGCTGTATC
#                                          GtattctggctgtatcGGGGGtacgatgctatactgtg
