#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import textwrap

# from pydna.seq import Seq
from pydna.dseq import Dseq
from pydna.utils import eq

from Bio.Restriction import (
    Acc65I, ApaI, BamHI, BglII, BsaI, Bsp120I, EcoRI, EcoRV,
    KpnI, MaeII, NmeDI, NotI, PacI, PstI, RestrictionBatch, TaiI,
    BspLI, SmaI
)

from seguid import ldseguid
from seguid import cdseguid

def test_dseq():

    x = Dseq( "gGGATCC",
             "  CCTAGG"[::-1], 1)

    y = Dseq(" gGGATCC",
             "  CCTAGG"[::-1], 0)

    z = Dseq("  gGGATCC",
              "  CCTAGG"[::-1], -1)

    assert x == y == z

    x = Dseq( " GGATCC",
             " gCCTAGG"[::-1], 1)

    y = Dseq("  GGATCC",
             " gCCTAGG"[::-1], 0)

    z = Dseq("   GGATCC",
              " gCCTAGG"[::-1], -1)

    assert x == y == z

def test_cut1():


    """
    Acc65I.search(_Seq("GGTACC"))
    Out  [11]: [2]


        012345
        GGTACC
        CCATGG


    KpnI.search(_Seq("GGTACC"))
    Out  [12]: [6]
    """



    """


    |--||-------||-------||----|
    aaaGGTACCcccGGGCCCgggACGTaaa
    tttCCATGGgggCCCGGGcccTGCAttt
    |------||-------||-----||--|

    aaaG         GGCCCgggA
    tttCCATG         GcccTGC         Acc65I Bsp120I MaeII

        GTACCcccG         CGTaaa
            GgggCCCGG       Attt

    |   |        |        |           0  4  13  22
    0123456789111111111122222222
              012345678901234567


    aaaGGTAC         CgggACGT
    tttC         CCGGGccc             KpnI ApaI TaiI
                                      0 4 13 21
            CcccGGGCC        aaa
        CATGGgggC        TGCAttt
    |------||------||------||--|
    aaaGGTACCcccGGGCCCgggACGTaaa
    tttCCATGGgggCCCGGGcccTGCAttt
    |--||-------||-------||----|

    """

    lds = Dseq("aaaGGTACCcccGGGCCCgggACGTaaa")

    frags = lds.cut((Acc65I, Bsp120I, MaeII))

    first, second, third, fourth = frags

    assert first + second + third + fourth == lds

    # TODO: remove
    # assert (first.pos, second.pos, third.pos, fourth.pos) == (0, 4, 13, 22)

    frags2 = lds.cut((KpnI, ApaI, TaiI))

    first2, second2, third2, fourth2 = frags2

    assert first2 + second2 + third2 + fourth2 == lds

    # TODO: remove
    # assert (first2.pos, second2.pos, third2.pos, fourth2.pos) == (0, 4, 13, 21)

    # G      GTACC
    # CCATG      G
    first = Dseq.from_representation(
    """
    Dseq(-5)
    G
    CCATG
    """)
    second = Dseq.from_representation(
    """
    Dseq(-5)
    GTACC
        G
    """)

    assert (first, second) == Dseq("GGTACC").cut(Acc65I)


    # GGTAC      C
    # C      CATGG
    first = Dseq.from_representation(
    """
    Dseq(-5)
    GGTAC
    C
    """)
    second = Dseq.from_representation(
    """
    Dseq(-5)
        C
    CATGG
    """)
    assert (first, second) == Dseq("GGTACC").cut(KpnI)


    # GGT      ACC
    # CCA      TGG
    first = Dseq.from_representation(
    """
    Dseq(-3)
    GGT
    CCA
    """)
    second = Dseq.from_representation(
    """
    Dseq(-3)
    ACC
    TGG
    """)
    assert (first, second) == Dseq("GGTACC").cut(BspLI)

    #  GTACCG
    #      GCCATG
    lin = Dseq.from_representation(
    """
    Dseq(-10)
    GTACCG
        GCCATG
    """)
    assert (lin,) == Dseq("GGTACC", circular=True).cut(Acc65I)

    #      CGGTAC
    #  CATGGC
    lin = Dseq.from_representation(
    """
    Dseq(-10)
        CGGTAC
    CATGGC
    """)
    assert (lin,) == Dseq("GGTACC", circular=True).cut(KpnI)

    # ACCGGT
    # TGGCCA
    lin = Dseq.from_representation(
    """
    Dseq(-6)
    ACCGGT
    TGGCCA
    """)

    assert (lin,) == Dseq("GGTACC", circular=True).cut(BspLI)



def test_cas9():

    s = Dseq("gattcatgcatgtagcttacgtagtct")

    RNA = "catgcatgtagcttacgtag"

    assert slice(0, 21, 1), slice(21, 27, 1) == s.cas9(RNA)


def test_initialization():

    obj = Dseq(b"aaa")
    assert obj.crick == "ttt"
    assert obj.watson == "aaa"

    obj = Dseq("a", "t", 0)
    assert obj * 3 == Dseq("aaa", "ttt", 0)
    assert not obj == 123
    assert obj * 0 == Dseq("")

    with pytest.raises(TypeError):
        obj * 2.3

    assert obj.seguid() == "ldseguid=ydezQsYTZgUCcb3-adxMaq_Xf8g"

    assert obj == Dseq("a", "t", circular=False)

    with pytest.raises(ValueError):
        Dseq("a", ovhg=0)

    with pytest.raises(ValueError):
        Dseq("ttt", "tt")

    with pytest.raises(ValueError):
        Dseq("ttt", "aa")

    obj2 = Dseq("gata")

    assert not obj2.circular

    lin_seq = Dseq("gt")
    cir_seq = lin_seq.looped()

    assert not lin_seq.circular
    assert cir_seq.circular

    assert Dseq("gt", circular=False) == lin_seq
    assert Dseq("gt", circular=True) == cir_seq

    assert Dseq.quick(b"A") == Dseq("A")
    assert Dseq.quick(b"A", circular=True) == Dseq("A", circular=True)

    obj3 = Dseq.from_representation(
        """
                                   GGATCC
                                   CCTAGG
                                   """
    )
    assert obj3.ovhg == 0

    obj3 = Dseq.from_representation(
        """
                                   aGGATCC
                                    CCTAGGg
                                   """
    )
    assert obj3.ovhg == -1

    obj3 = Dseq.from_representation(
        """
                                    GGATCCg
                                   aCCTAGG
                                   """
    )
    assert obj3.ovhg == 1

    assert Dseq("g", "", 0) == Dseq("p")
    assert Dseq("a", "", 0) == Dseq("e")
    assert Dseq("t", "", 0) == Dseq("x")
    assert Dseq("c", "", 0) == Dseq("i")

    assert Dseq("", "g", 0) == Dseq("j")
    assert Dseq("", "a", 0) == Dseq("z")
    assert Dseq("", "t", 0) == Dseq("f")
    assert Dseq("", "c", 0) == Dseq("q")

    s = Dseq.from_representation(
    """
    Dseq(-6)
    G A C
    CcTaGg
    """)
    assert Dseq("G A C ", "CcTaGg"[::-1], 0) == s

    s = Dseq.from_representation(
    """
    Dseq(-6)
    G A C
    C T G
    """)
    assert Dseq("G A C ", "C T G"[::-1], 0) == s # TODO: Discuss if this should raise an exeption.

    s = Dseq.from_representation(
    """
    Dseq(-6)
    GA T CC
    CTAACGG
    """)
    assert Dseq("GA T CC", "CTAACGG"[::-1], 0) == s



def test_cut_around_and_religate():

    def cut_and_religate_Dseq(seq_string, enz, top):
        ds = Dseq(seq_string, circular=not top)
        frags = list(ds.cut(enz))
        if not frags:
            return
        a = frags.pop(0)

        for f in frags:
            a += f
        if not top:
            a = a.looped()
        assert eq(a, ds)

    seqs = [
        ("aaaGGTACCcccGGTACCCgggGGTACCttt", BamHI, False),
        ("aaaGGTACCcccGGTACCCgggGGTACCttt", Acc65I, False),
        ("aaaGGTACCcccGGTACCCgggGGTACCttt", KpnI, False),
        ("aaaGGTACCcccGGTACCCgggGGTACCttt", Acc65I, True),
        ("aaaGGTACCcccGGTACCCgggGGTACCttt", KpnI, True),
        ("aaaGGTACCcccGGTACCCgggGGTACCttt", BamHI, True),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [Acc65I, BamHI], False),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [KpnI, BamHI], False),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [BamHI, Acc65I], False),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [BamHI, KpnI], False),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [Acc65I, BamHI], True),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [KpnI, BamHI], True),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [BamHI, Acc65I], True),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [BamHI, KpnI], True),
    ]

    for sek, enz, lin in seqs:
        for i in range(len(sek)):
            zek = sek[i:] + sek[:i]
            cut_and_religate_Dseq(zek, enz, lin)


def test_Dseq_cutting_adding():

    a = Dseq(
        "GGATCCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtGGATCC",
        "CCTAGGagtagatgatagtagcatcgcatgactagataagacgacgagtagtagccatgagagatattaatatatatatacgcgcaCCTAGG"[::-1],
        ovhg=0,
    )

    b = a.cut(BamHI)[1]

    assert b.watson == "GATCCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtG"
    assert b.crick == "GATCCacgcgcatatatatataattatagagagtaccgatgatgagcagcagaatagatcagtacgctacgatgatagtagatgaG"
    c = Dseq(
        "nCTGCAGtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtGAATTCn",
        "nGACGTCagtagatgatagtagcatcgcatgactagataagacgacgagtagtagccatgagagatattaatatatatatacgcgcaCTTAAGn"[::-1],
        ovhg=0,
    )

    f, d, _ = c.cut((EcoRI, PstI))

    assert d.watson == "GtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtG"
    assert d.crick == "AATTCacgcgcatatatatataattatagagagtaccgatgatgagcagcagaatagatcagtacgctacgatgatagtagatgaCTGCA"

    e = Dseq(
        "nGAATTCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtCTGCAGn",
        "nCTTAAGagtagatgatagtagcatcgcatgactagataagacgacgagtagtagccatgagagatattaatatatatatacgcgcaGACGTCn"[::-1],
        ovhg=0,
    )

    f = e.cut((EcoRI, PstI))[1]

    assert f.watson == "AATTCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtCTGCA"
    assert f.crick == "GacgcgcatatatatataattatagagagtaccgatgatgagcagcagaatagatcagtacgctacgatgatagtagatgaG"


def test_dseq():

    obj1 = Dseq("a", "t", circular=True)
    obj2 = Dseq("a", "t")

    with pytest.raises(TypeError):
        obj1 + obj2

    with pytest.raises(TypeError):
        obj2 + obj1

    with pytest.raises(TypeError):
        obj1 + ""


    # with pytest.raises(AttributeError):
    #     obj2 + ""

    assert obj2 + "" == obj2
    assert "" + obj2 == obj2

    obj1 = Dseq("at", "t")
    obj2 = Dseq("a", "t")

    with pytest.raises(TypeError):
        obj1 + obj2

    obj = Dseq("aaa", "ttt", circular=True)
    assert obj[1:2] == Dseq("a", "t", 0)

    assert obj[:] == Dseq("aaa", "ttt", circular=False)

    obj = Dseq("atg", "cat", 0, circular=False)

    assert obj[1:2]._data == b"atg"[1:2]

    assert obj[2:1]._data == b"atg"[2:1]

    assert obj.reverse_complement() == obj.rc() == Dseq("cat", "atg", 0)

    obj = Dseq("atg", "cat", circular=True)

    assert obj.looped() == obj

    assert obj[:] == Dseq("atg", "cat", 0, circular=False)

    assert obj[1:2]._data == b"atg"[1:2]

    assert obj[2:1]._data == b"ga"

    obj = Dseq("G", "", 0)
    assert obj.five_prime_end() == ("single", "g")
    assert obj.three_prime_end() == ("single", "g")
    obj = Dseq("", "C", 0)
    assert obj.five_prime_end() == ("single", "c")
    assert obj.three_prime_end() == ("single", "c")




    obj = Dseq("ccGGATCC", "aaggatcc", -2)
    # assert obj._data == b"ccGGATCCtt"
    assert obj._data == b"iiGGATCCzz"
    assert str(obj.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-10)
    ccGGATCC
      CCTAGGaa
    """
    ).strip()


    assert repr(obj) == rpr

    assert obj[3] == Dseq("G", "c", 0)

    assert obj.fill_in() == Dseq("ccGGATCCtt", "aaggatccgg", 0)

    assert obj + Dseq("") == obj
    assert Dseq("") + obj == obj

    obj = Dseq("gatcAAAAAA", "gatcTTTTTT")
    assert obj.fill_in("gatc") == Dseq("gatcAAAAAAgatc", "gatcTTTTTTgatc")
    assert obj.fill_in("atc") == obj
    assert obj.fill_in("ac") == obj
    assert obj.fill_in("at") == obj

    obj = Dseq("AAAAAAgatc", "TTTTTTgatc")
    assert obj.fill_in("gatc") == obj
    assert obj.fill_in("atc") == obj
    assert obj.fill_in("ac") == obj
    assert obj.fill_in("at") == obj

    obj = Dseq("gatcAAAAAA", "gatcTTTTTT")
    assert obj.t4() == Dseq("gatcAAAAAAgatc", "gatcTTTTTTgatc")

    assert obj.t4("at") == obj
    assert obj.t4("atg") == Dseq("gatcAAAAAAgat", "gatcTTTTTTgat")
    assert obj.t4("atgc") == Dseq("gatcAAAAAAgatc", "gatcTTTTTTgatc")
    assert obj.mung() == Dseq("AAAAAA", "TTTTTT")

    obj = Dseq("AAAAAAgatc", "TTTTTTgatc")
    assert obj.t4() == obj.t4("at") == Dseq("AAAAAA")
    assert obj.t4("atc") == obj.t4("atg") == obj.t4("atcg") == Dseq("AAAAAA")

    assert Dseq("GGATCC", "GGATCC").t4() == Dseq("GGATCC", "GGATCC")
    assert Dseq("GGATCCa", "GGATCC").t4() == Dseq("GGATCC", "GGATCC")
    assert Dseq("aGGATCC", "GGATCC").t4() == Dseq("aGGATCC", "GGATCCt")
    assert Dseq("aGGATCCa", "GGATCC").t4() == Dseq("aGGATCC", "GGATCCt")
    assert Dseq("GGATCC", "aGGATCC").t4() == Dseq("GGATCCt", "aGGATCC")
    assert Dseq("GGATCC", "GGATCCa").t4() == Dseq("GGATCC", "GGATCC")
    assert Dseq("GGATCC", "aGGATCCa").t4() == Dseq("GGATCCt", "aGGATCC")

    assert Dseq("GGATCC", "ATCC").t4("g") == Dseq("gg", "", ovhg=0)
    assert Dseq("GGATCC", "GGATCC").t4("gat") == Dseq("ggat", "ggat", ovhg=-2)

    a2 = Dseq("ccGGATCCaa", "ggatcc", -2)
    # assert a2._data == b"ccGGATCCaa"
    assert a2._data == b"iiGGATCCee"
    assert str(a2.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-10)
    ccGGATCCaa
      CCTAGG
    """
    ).strip()

    assert repr(a2) == rpr

    a3 = Dseq("ccGGATCC", "ggatcc", -2)
    # assert a3._data == b"ccGGATCC"
    assert a3._data == b"iiGGATCC"
    assert str(a3.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-8)
    ccGGATCC
      CCTAGG
    """
    ).strip()
    assert repr(a3) == rpr

    b = Dseq("GGATCC", "aaggatcccc", 2)
    # assert b._data == b"ggGGATCCtt"
    assert b._data == b"qqGGATCCzz"
    assert str(b.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-10)
      GGATCC
    ccCCTAGGaa
    """
    ).strip()
    assert repr(b) == rpr

    b2 = Dseq("GGATCCaa", "ggatcccc", 2)
    # assert b2._data == b"ggGGATCCaa"
    assert b2._data == b"qqGGATCCee"
    assert str(b2.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-10)
      GGATCCaa
    ccCCTAGG
    """
    ).strip()
    assert repr(b2) == rpr
    assert b2.seguid() == "ldseguid=F0z-LxHZqAK3HvqQiqjM7A28daE"
    assert b2.rc().seguid() == "ldseguid=F0z-LxHZqAK3HvqQiqjM7A28daE"

    b3 = Dseq("GGATCC", "ggatcccc", 2)
    # assert b3._data == b"ggGGATCC"
    assert b3._data == b"qqGGATCC"
    assert str(b3.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-8)
      GGATCC
    ccCCTAGG
    """
    ).strip()
    assert repr(b3) == rpr

    c = Dseq("GGATCCaaa", "ggatcc", 0)
    # assert c._data == b"GGATCCaaa"
    assert c._data == b"GGATCCeee"
    assert str(c.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-9)
    GGATCCaaa
    CCTAGG
    """
    ).strip()
    assert repr(c) == rpr

    d = Dseq("GGATCC", "aaaggatcc", 0)
    # assert d._data == b"GGATCCttt"
    assert d._data == b"GGATCCzzz"
    assert str(d.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-9)
    GGATCC
    CCTAGGaaa
    """
    ).strip()
    assert repr(d) == rpr

    obj = Dseq("GGATCCaaa", "ggatcc", 0)


    frag1 = Dseq("G", "gatcc", 0)
    frag2 = Dseq("GATCCaaa", "g", -4)

    assert obj.cut(BamHI) == (frag1, frag2)

    assert frag1 + frag2 == obj

    assert obj.seguid() == "ldseguid=qvssQpZe_4SlasGZYdKJSkuvQtc"

    assert frag1.seguid() == "ldseguid=jcVhCJ9Aa8aIQdBlkSU_XHTWmDc"
    assert frag2.seguid() == "ldseguid=SO1HxaZPDpcj-QffzS-mfF6_eag"

    assert frag1.rc().seguid() == "ldseguid=jcVhCJ9Aa8aIQdBlkSU_XHTWmDc"
    assert frag2.rc().seguid() == "ldseguid=SO1HxaZPDpcj-QffzS-mfF6_eag"

    obj = Dseq("tagcgtagctgtagtatgtgatctggtcta", "tagaccagatcacatactacagctacgcta")
    assert repr(obj) == "Dseq(-30)\ntagcgtagctgtagtatgtgatctggtcta\natcgcatcgacatcatacactagaccagat"

    obj2 = Dseq("tagcgtagctgtagtatgtgatctggtcta")

    obj3 = obj = Dseq("tagcgtagctgtagtatgtgatctggtcta", "tagaccagatcacatactacagctacgcta", 0)

    assert obj == obj2 == obj3

    assert obj.find("ggatcc") == -1

    assert obj.find("tgtagta") == 9

    obj = Dseq("tagcgtagctgtagtatgtgatctggtctaa", "ttagaccagatcacatactacagctacgcta")

    obj = Dseq("tagcgtagctgtagtatgtgatctggtctaa", "CCCttagaccagatcacatactacagctacgcta")

    assert repr(obj) == "Dseq(-34)\ntagc..ctaa\natcg..gattCCC"

    obj = Dseq("tagcgtagctgtagtatgtgatctggtctaaCCC", "ttagaccagatcacatactacagctacgcta")

    assert repr(obj) == "Dseq(-34)\ntagc..ctaaCCC\natcg..gatt"

    obj = Dseq("agcgtagctgtagtatgtgatctggtctaa", "ttagaccagatcacatactacagctacgcta")
    assert repr(obj) == "Dseq(-31)\n agcg..ctaa\natcgc..gatt"

    obj = Dseq("Atagcgtagctgtagtatgtgatctggtctaa", "ttagaccagatcacatactacagctacgcta")
    assert repr(obj) == "Dseq(-32)\nAtagc..ctaa\n atcg..gatt"

    obj = Dseq("tagcgtagctgtagtatgtgatctggtctaa", "tatcgcatcgacatcatacactagaccagatt"[::-1])

    assert repr(obj) == "Dseq(-32)\n tagc..ctaa\ntatcg..gatt"

    obj1 = Dseq(
        "tagcgtagctgtagtatgtgatctggtcta",
        "tagaccagatcacatactacagctacgcta",
        circular=True,
    )
    obj2 = Dseq(
        "tagcgtagctgtagtatgtgatctggtcta",
        "tagaccagatcacatactacagctacgcta",
        circular=True,
    )
    obj3 = Dseq(
        "tagcgtagctgtagtatgtgatctggtcta",
        "tagaccagatcacatactacagctacgcta",
        circular=True,
    )

    assert obj1 == obj2 == obj3

    assert obj1.find("ggatcc") == -1

    assert obj1.find("tgtagta") == 9

    assert obj1.find("gtcta" "tag") == 25 # find substring over origin

    assert Dseq("tagcgtagctgtagtatgtgatctggtcta", "tagaccagatcacatactacagctacgcta").looped() == obj1



    obj = Dseq("ggatcc")

    assert BglII in obj.no_cutters()
    assert BamHI not in obj.no_cutters()

    assert BamHI in obj.unique_cutters()

    assert BamHI in obj.once_cutters()

    assert BamHI in (obj + obj).twice_cutters()
    assert BamHI not in obj.twice_cutters()

    assert BamHI in obj.n_cutters(1)
    assert BamHI in obj.cutters()



    rb = RestrictionBatch((BamHI, BglII))

    assert obj.cut(rb) == obj.cut(BamHI, BglII) == obj.cut(BglII, BamHI)

    obj = Dseq("ggatccAGATCT")

    assert obj.cut(rb) == obj.cut(BamHI, BglII) == obj.cut(BglII, BamHI)

    obj = Dseq("AGATCTggatcc")

    assert obj.cut(rb) == obj.cut(BamHI, BglII) == obj.cut(BglII, BamHI)

    obj = Dseq("ggatccAGATCT", circular=True)
    # TODO: address this test change Related to https://github.com/pydna-group/pydna/issues/78
    assert obj.cut(rb) == obj.cut(BamHI, BglII) == obj.cut(BglII, BamHI)

    obj = Dseq("AGATCTggatcc", circular=True)

    assert obj.cut(rb) == obj.cut(BglII, BamHI) == obj.cut(BamHI, BglII)


def test_Dseq_slicing():

    a = Dseq("ggatcc", "ggatcc", 0)

    assert a[:].watson == a.watson
    assert a[:].crick == a.crick
    assert a.ovhg == a[:].ovhg
    b, c = a.cut(BamHI)
    d = b[1:5]
    e = d.rc()
    # assert  d+e == Dseq("gatc","gatc",0)
    assert e + d == Dseq("gatc", "gatc", 0)


def test_Dseq_slicing2():



    a = Dseq("aaGGATCCnnnnnnnnnGAATTCccc", circular=True)
    # TODO: address this test change Related to https://github.com/pydna-group/pydna/issues/78

    assert a.cut(
        EcoRI,
        BamHI,
        KpnI,
    ) == a.cut(
        BamHI,
        EcoRI,
        KpnI,
    )


def test_Dseq___getitem__():
    # test the slicing

    assert str(Dseq("gatc")[1:3]) == "gatc"[1:3]
    assert str(Dseq("gatc")[:]) == "gatc"[:]
    assert str(Dseq("gatc")[0:0]) == "gatc"[0:0]
    assert str(Dseq("gatc")[:0]) == "gatc"[:0]
    assert str(Dseq("gatc")[0:]) == "gatc"[0:]
    assert str(Dseq("gatc")[None:None]) == "gatc"[None:None]
    assert str(Dseq("gatc")[None:0]) == "gatc"[None:0]
    assert str(Dseq("gatc")[0:None]) == "gatc"[0:None]

    s = Dseq("GGATCC", circular=False)
    assert s[1:-1] == Dseq("GATC", circular=False)
    t = Dseq("GGATCC", circular=True)
    assert t[1:5] == Dseq("GATC")
    assert t[1:5].__dict__ == Dseq("GATC").__dict__
    assert s[1:5] == Dseq("GATC")
    assert s[1:5] == Dseq("GATC", circular=False)
    assert s[5:1:-1] == Dseq("CCTA")

    assert t[5:1] == Dseq("CG")
    assert s[9:1] == Dseq("")
    assert t[9:1] == Dseq("TCCG")  # XXX: This is important!
    assert t[1:9] == Dseq("GATCCGGA")  # XXX: This is important!


    # Indexing of full circular molecule (https://github.com/pydna-group/pydna/issues/161)
    s = Dseq("GGATCC", circular=True)
    str_seq = str(s)
    for shift in range(len(s)):
        assert str(s[shift:shift]) == str_seq[shift:] + str_seq[:shift]


def test_cut_circular():

    test = "aaaaaaGGTACCggtctcaaaa"

    for i in range(len(test)):
        nt = test[i:] + test[:i]

        a = Dseq(nt, circular=True).cut(Acc65I)[0]  # G^GTACC

        assert a.watson.upper() == "GTACCGGTCTCAAAAAAAAAAG"
        assert a.crick.upper() == "GTACCTTTTTTTTTTGAGACCG"
        assert a.ovhg == -4  # CggtctcaaaaaaaaaaGGTAC
        b = Dseq(nt, circular=True).cut(KpnI)[0]  # GGTAC^C
        assert b.watson.upper() == "CGGTCTCAAAAAAAAAAGGTAC"
        assert b.crick.upper() == "CTTTTTTTTTTGAGACCGGTAC"
        assert b.ovhg == 4
        c = Dseq(nt, circular=True).cut(BsaI)[0]  # ggtctcnnn
        assert c.watson.upper() == "AAAAAAAAAGGTACCGGTCTCA"
        assert c.crick.upper() == "TTTTTGAGACCGGTACCTTTTT"
        assert c.ovhg == -4
        d = Dseq(nt, circular=True).cut(NotI)
        assert d == ()


def test_repr():


    a = Dseq("gattcgtatgctgatcgtacgtactgaaaac")

    assert repr(a) == "Dseq(-31)\ngatt..aaac\nctaa..tttg"

    b = Dseq("gattcgtatgctgatcgtacgtactgaaaac", "gactagcatgcatgacttttg"[::-1])

    assert repr(b) == "Dseq(-31)\ngattcgtatgctga..aaac\n          gact..tttg"

    c = Dseq("gattcgtatgctgatcgtacgtactgaaaac", "actagcatgcatgacttttg"[::-1])

    assert repr(c) == "Dseq(-31)\ngatt..atgctgat..aaac\n          acta..tttg"

    d = Dseq("gattcgtatgctgatcgtacg", "gactagcatgc"[::-1])

    assert repr(d) == "Dseq(-21)\ngattcgtatgctgatcgtacg\n          gactagcatgc"

    e = Dseq("gactagcatgcatgacttttc", "gattcgtatgctgatcgtacgtactgaaaag"[::-1])

    assert repr(e) == "Dseq(-31)\n          gact..tttc\ngattcgtatgctga..aaag"

    f = Dseq("Cgactagcatgcatgacttttc", "gattcgtatgctgatcgtacgtactgaaaag"[::-1])

    assert repr(f) == "Dseq(-31)\n         Cgac..tttc\ngattcgtatGctg..aaag"

    g = Dseq("gattcgtatgctgatcgtacgtactgaaaac", "ctaagcatacgactagc"[::-1])

    assert repr(g) == "Dseq(-31)\ngatt..atcgtacg..aaac\nctaa..tagc"

    h = Dseq("cgtatgctgatcgtacgtactgaaaac", "gcatacgactagc"[::-1])

    assert repr(h) == "Dseq(-27)\ncgtatgctgatcgtacgtactgaaaac\ngcatacgactagc"

    i = Dseq("cgtatgctgatcgtacgtactgaaaacagact", "gcatacgactagc"[::-1])

    assert repr(i) == "Dseq(-32)\ncgta..atcgtacg..gact\ngcat..tagc"

    j = Dseq("gtttcgtatgctgatcgtacgtactgaaaac", "acAAGGAGAGAtg", ovhg=11)

    assert repr(j) == "Dseq(-42)\n          gtttcg..aaac\ngtAG..GGAAca"

    k = Dseq("g", "gattcgtatgctgatcgtacgtactgaaaac", ovhg=0)

    assert repr(k) == "Dseq(-31)\ng\ncaaaa..ttag"

    x = Dseq("gattcgtatgctgatcgtacgtactgaaaa")

    assert repr(x) == "Dseq(-30)\ngattcgtatgctgatcgtacgtactgaaaa\nctaagcatacgactagcatgcatgactttt"

    y = Dseq("gattcgtatgctgatcgtacgtactgaaaa", "gactagcatgcatgactttt"[::-1])

    assert repr(y) == "Dseq(-30)\ngattcgtatgctgatcgtacgtactgaaaa\n          gactagcatgcatgactttt"

    z = Dseq("gattcgtatgctgatcgtacgtactgaaaa", "actagcatgcatgactttt"[::-1])

    assert repr(z) == "Dseq(-30)\ngattcgtatgctgatcgtacgtactgaaaa\n           actagcatgcatgactttt"


def test_shifted():

    a = Dseq("gatc", circular=True)

    assert a.shifted(1) == Dseq("atcg", circular=True)

    assert a.shifted(4) == a

    b = Dseq("gatc", circular=False)
    with pytest.raises(TypeError):
        b.shifted(1)

    # Shifted with zero gives a copy of the sequence, not the same sequence
    assert a.shifted(0) == a
    assert a.shifted(0) is not a


def test_looped():

    # Looping a circular sequence should return a copy of the sequence
    # not the same sequence



    a = Dseq("gatc", circular=True)

    assert a.looped() == a
    assert a.looped() is not a


def test_misc():


    x = Dseq("ctcgGCGGCCGCcagcggccg", circular=True)



    a, b = x.cut(NotI)

    z = (a + b).looped()
    # TODO: address this test change Related to https://github.com/pydna-group/pydna/issues/78
    assert z.shifted(-6) == x


def test_cut_missing_enzyme():


    x = Dseq("ctcgGCGGCCGCcagcggccg")

    assert x.cut(PstI) == ()

    x = Dseq("ctcgGCGGCCGCcagcggccg", circular=True)

    assert x.cut(PstI) == ()


def test_cut_with_no_enzymes():


    x = Dseq("ctcgGCGGCCGCcagcggccg")

    assert x.cut([]) == ()

    x = Dseq("ctcgGCGGCCGCcagcggccg", circular=True)

    assert x.cut([]) == ()


def test_transcribe():


    x = Dseq("ATGAAATAA")

    assert str(x.transcribe()) == "AUGAAAUAA"

    assert str(x.reverse_complement().transcribe()) == "UUAUUUCAU"


def test_translate():


    x = Dseq("ATGAAATAA")

    assert str(x.translate()) == "MK*"

    assert str(x.reverse_complement().translate()) == "LFH"


def test_from_full_sequence_and_overhangs():


    test_cases = [
        (2, 2, "AAAA", "TTTT"),
        (-2, 2, "AAAAAA", "TT"),
        (2, -2, "AA", "TTTTTT"),
        (-2, -2, "AAAA", "TTTT"),
        (0, 0, "AAAAAA", "TTTTTT"),
    ]
    for crick_ovhg, watson_ovhg, watson, crick in test_cases:
        dseq_1 = Dseq.from_full_sequence_and_overhangs("AAAAAA", crick_ovhg=crick_ovhg, watson_ovhg=watson_ovhg)
        dseq_2 = Dseq(watson, crick, ovhg=crick_ovhg, circular=False)

        assert dseq_1 == dseq_2
        assert dseq_2.watson_ovhg == watson_ovhg


def test_right_end_position():

    test_cases = [
        ("AAA", "TT", "AAE", (3, 2)),
        ("AA", "TTT", "AAF", (2, 3)),
        ("AAA", "TTT", "AAA", (3, 3)),
    ]
    for watson, crick, dscode, expected in test_cases:
        dseq = Dseq(dscode, circular=False)
        assert dseq.right_end_position() == expected


def test_left_end_position():

    test_cases = [
        ("AAA", "TT", "EAA", (0, 1), -1),
        ("AA", "TTT", "FAA", (1, 0), 1),
        ("AAA", "TTT", "AAA", (0, 0), 0),
    ]
    for watson, crick, dscode, expected, ovhg in test_cases:
        dseq = Dseq(dscode, circular=False)
        assert dseq.left_end_position() == expected


def test_apply_cut():

    seq = Dseq("aaGAATTCaa", circular=False)

    # A cut where both sides are None returns the same sequence
    assert seq.apply_cut(None, None) == seq

    # A cut where one side is None leaves that side intact
    EcoRI_cut = ((3, -4), None)

    assert seq.apply_cut(None, EcoRI_cut) == Dseq.from_full_sequence_and_overhangs(
        "aaGAATT", watson_ovhg=-4, crick_ovhg=0
    )
    assert seq.apply_cut(EcoRI_cut, None) == Dseq.from_full_sequence_and_overhangs(
        "AATTCaa", watson_ovhg=0, crick_ovhg=-4
    )

    # It respects the original overhang
    seq = Dseq.from_full_sequence_and_overhangs("aaGAATTCaa", watson_ovhg=1, crick_ovhg=1)
    assert seq.apply_cut(None, EcoRI_cut) == Dseq.from_full_sequence_and_overhangs(
        "aaGAATT", watson_ovhg=-4, crick_ovhg=1
    )
    assert seq.apply_cut(EcoRI_cut, None) == Dseq.from_full_sequence_and_overhangs(
        "AATTCaa", watson_ovhg=1, crick_ovhg=-4
    )

    seq = Dseq.from_full_sequence_and_overhangs("aaGAATTCaa", watson_ovhg=-1, crick_ovhg=-1)
    assert seq.apply_cut(None, EcoRI_cut) == Dseq.from_full_sequence_and_overhangs(
        "aaGAATT", watson_ovhg=-4, crick_ovhg=-1
    )
    assert seq.apply_cut(EcoRI_cut, None) == Dseq.from_full_sequence_and_overhangs(
        "AATTCaa", watson_ovhg=-1, crick_ovhg=-4
    )

    # A repeated cut in a circular molecule opens it up
    seq = Dseq("aaGAATTCaa", circular=True)
    assert seq.apply_cut(EcoRI_cut, EcoRI_cut) == Dseq.from_full_sequence_and_overhangs(
        "AATTCaaaaGAATT", watson_ovhg=-4, crick_ovhg=-4
    )

    # Two cuts extract a subsequence
    seq = Dseq("aaGAATTCaaGAATTCaa", circular=True)
    EcoRI_cut_2 = ((11, -4), None)
    assert seq.apply_cut(EcoRI_cut, EcoRI_cut_2) == Dseq.from_full_sequence_and_overhangs(
        "AATTCaaGAATT", watson_ovhg=-4, crick_ovhg=-4
    )

    # Overlapping cuts should return an error
    seq = Dseq("aaGAATTCaa", circular=True)
    first_cuts = [
        ((3, -4), BamHI),
        ((7, 4), BamHI),
        # Spanning the origin
        ((9, -8), BamHI),
        ((8, 8), BamHI),
    ]

    overlapping_cuts = [
        ((4, -4), EcoRI),
        ((2, -4), EcoRI),
        ((2, -6), EcoRI),
        ((8, 4), EcoRI),
        ((6, 4), EcoRI),
        ((8, 6), EcoRI),
        # Spanning the origin
        ((7, -8), EcoRI),
        ((6, 8), EcoRI),
    ]

    for first_cut in first_cuts:
        for second_cut in overlapping_cuts:
            try:
                seq.apply_cut(first_cut, second_cut)
            except ValueError as e:
                assert e.args[0] == "Cuts by BamHI EcoRI overlap."
            else:
                print(first_cut, second_cut)
                assert False, "Expected ValueError"

    # Rotating the sequence, apply the same cut
    seq = Dseq("acgtATGaatt", circular=True)
    for shift in range(len(seq)):

        seq_shifted = seq.shifted(shift)

        start = 4 - shift
        if start < 0:
            start += len(seq)
        # Cut with negative ovhg
        new_cut = ((start, -3), None)
        out = seq_shifted.apply_cut(new_cut, new_cut)
        assert out.to_blunt_string() == "ATGaattacgtATG"

        # Cut with positive ovhg
        start = (start + 3) % len(seq)
        new_cut = ((start, 3), None)
        out = seq_shifted.apply_cut(new_cut, new_cut)
        assert out.to_blunt_string() == "ATGaattacgtATG"

        # A blunt cut
        start = 4 - shift
        new_cut = ((start, 0), None)
        out = seq_shifted.apply_cut(new_cut, new_cut)
        assert out.to_blunt_string() == "ATGaattacgt"


def test_cutsite_is_valid():

    # Works for circular case
    seqs = ["GAATTC", "TTAATTAAC", "GATATC"]
    enzs = [EcoRI, PacI, EcoRV]
    for seq, enz in zip(seqs, enzs):
        dseq = Dseq(seq, circular=True)
        for shift in range(len(seq)):
            dseq_shifted = dseq.shifted(shift)
            (cutsite,) = dseq_shifted.get_cutsites([enz])

            assert dseq_shifted.cutsite_is_valid(cutsite)

    # Works for overhangs
    seqs = ["GAATTC", "TTAATTAA", "GATATC"]
    for seq, enz in zip(seqs, enzs):
        for ovhg in [-1, 0, 1]:
            dseq = Dseq.from_full_sequence_and_overhangs(seq, ovhg, 0)
            if ovhg != 0:
                assert len(dseq.get_cutsites([enz])) == 0
            else:
                assert len(dseq.get_cutsites([enz])) == 1

            dseq = Dseq.from_full_sequence_and_overhangs(seq, 0, ovhg)
            if ovhg != 0:
                assert len(dseq.get_cutsites([enz])) == 0
            else:
                assert len(dseq.get_cutsites([enz])) == 1

    # Special cases:
    dseq = Dseq.from_full_sequence_and_overhangs("AAAAAATTTTTTTGCCGGCAAAAAAAATTTTT", 0, 0)
    assert len(dseq.get_cutsites([NmeDI])) == 2
    # Remove left cutting place
    assert len(dseq[2:].get_cutsites([NmeDI])) == 1
    # Remove right cutting place
    assert len(dseq[:-2].get_cutsites([NmeDI])) == 1
    # Remove both cutting places
    assert len(dseq[2:-2].get_cutsites([NmeDI])) == 0

    # overhang left side
    dseq = Dseq.from_full_sequence_and_overhangs("AAAAAATTTTTTTGCCGGCAAAAAAAATTTTT", -2, 0)
    assert len(dseq.get_cutsites([NmeDI])) == 1
    dseq = Dseq.from_full_sequence_and_overhangs("AAAAAATTTTTTTGCCGGCAAAAAAAATTTTT", 2, 0)
    assert len(dseq.get_cutsites([NmeDI])) == 1

    # overhang right side
    dseq = Dseq.from_full_sequence_and_overhangs("AAAAAATTTTTTTGCCGGCAAAAAAAATTTTT", 0, 2)
    assert len(dseq.get_cutsites([NmeDI])) == 1
    dseq = Dseq.from_full_sequence_and_overhangs("AAAAAATTTTTTTGCCGGCAAAAAAAATTTTT", 0, -2)
    assert len(dseq.get_cutsites([NmeDI])) == 1

    # overhang both sides
    dseq = Dseq.from_full_sequence_and_overhangs("AAAAAATTTTTTTGCCGGCAAAAAAAATTTTT", 2, 2)
    assert len(dseq.get_cutsites([NmeDI])) == 0
    dseq = Dseq.from_full_sequence_and_overhangs("AAAAAATTTTTTTGCCGGCAAAAAAAATTTTT", -2, -2)
    assert len(dseq.get_cutsites([NmeDI])) == 0

    # overhang on recognition site removes both cutting places
    dseq = Dseq.from_full_sequence_and_overhangs("AAAAAATTTTTTTGCCGGCAAAAAAAATTTTT", 16, 0)
    assert len(dseq.get_cutsites([NmeDI])) == 0
    dseq = Dseq.from_full_sequence_and_overhangs("AAAAAATTTTTTTGCCGGCAAAAAAAATTTTT", 0, 16)
    assert len(dseq.get_cutsites([NmeDI])) == 0


def test_get_cutsite_pairs():


    # in the test, we replace cuts by integers for clarity.

    dseq = Dseq("A")

    # Empty returns empty list
    assert dseq.get_cutsite_pairs([]) == []

    # Single cut on linear seq returns two fragments
    assert dseq.get_cutsite_pairs([1]) == [(None, 1), (1, None)]

    # Two cuts on linear seq return three fragments
    assert dseq.get_cutsite_pairs([1, 2]) == [(None, 1), (1, 2), (2, None)]

    dseq = Dseq("A", circular=True)

    # Empty returns empty list
    assert dseq.get_cutsite_pairs([]) == []

    # Single cut on circular seq returns opened molecule
    assert dseq.get_cutsite_pairs([1]) == [(1, 1)]

    # Two cuts on circular seq return 2 fragments
    assert dseq.get_cutsite_pairs([1, 2]) == [(1, 2), (2, 1)]


def test_get_cut_parameters():



    dseq = Dseq.from_full_sequence_and_overhangs("aaaACGTaaa", 3, 3)
    assert dseq.get_cut_parameters(None, True) == (*dseq.left_end_position(), dseq.ovhg)
    assert dseq.get_cut_parameters(None, False) == (*dseq.right_end_position(), dseq.watson_ovhg)

    assert dseq.get_cut_parameters(((4, -2), None), True) == (4, 6, -2)
    assert dseq.get_cut_parameters(((4, -2), None), False) == (4, 6, -2)
    assert dseq.get_cut_parameters(((6, 2), None), True) == (6, 4, 2)
    assert dseq.get_cut_parameters(((6, 2), None), False) == (6, 4, 2)

    dseq = Dseq("aaaACGTaaa", circular=True)

    # None cannot be used on circular molecules
    try:
        assert dseq.get_cut_parameters(None, True) == (*dseq.left_end_position(), dseq.ovhg)
    except AssertionError as e:
        assert e.args[0] == "Circular sequences should not have None cuts"
    else:
        assert False, "Expected AssertionError"

    try:
        assert dseq.get_cut_parameters(None, False) == (*dseq.right_end_position(), dseq.watson_ovhg)
    except AssertionError as e:
        assert e.args[0] == "Circular sequences should not have None cuts"
    else:
        assert False, "Expected AssertionError"

    # "Normal" cuts
    assert dseq.get_cut_parameters(((4, -2), None), True) == (4, 6, -2)
    assert dseq.get_cut_parameters(((4, -2), None), False) == (4, 6, -2)
    assert dseq.get_cut_parameters(((6, 2), None), True) == (6, 4, 2)
    assert dseq.get_cut_parameters(((6, 2), None), False) == (6, 4, 2)

    # Origin-spannign cuts
    assert dseq.get_cut_parameters(((9, -2), None), True) == (9, 1, -2)
    assert dseq.get_cut_parameters(((9, -2), None), False) == (9, 1, -2)
    assert dseq.get_cut_parameters(((1, 2), None), True) == (1, 9, 2)
    assert dseq.get_cut_parameters(((1, 2), None), False) == (1, 9, 2)


def test_checksums():




    # AT
    # TA

    dlDNA_ldseguid = "odgytmQKSOnFEUorGIWK3NDjqUA"
    truth = f"ldseguid={dlDNA_ldseguid}"
    assert ldseguid("AT", "AT") == truth == Dseq("AT", "AT").seguid()

    #  -AT
    #  AT-

    dlDNA2_ldseguid = "-9xkp3UfucL4bSPxYODh8i9KFEE"
    truth = f"ldseguid={dlDNA2_ldseguid}"

    assert ldseguid("-AT", "-TA") == truth == Dseq("AT", "TA", 1).seguid()

    # TA-
    # -TA

    dlDNA3_ldseguid = "kI9qYVNRPF8epm2xem0ZUP8J-CI"
    truth = f"ldseguid={dlDNA3_ldseguid}"
    assert ldseguid("TA-", "AT-") == truth == Dseq("TA", "AT", -1).seguid()

    # CTATAG
    # --TA--

    dlDNA4_ldseguid = "ToSxUXWMCIKz-FYdXJ3Qq-bS_8o"
    truth = f"ldseguid={dlDNA4_ldseguid}"
    assert ldseguid("CTATAG", "--AT--") == truth == Dseq("CTATAG", "AT", -2).seguid()

    # --AT--
    # GATATC

    assert ldseguid("--AT--", "CTATAG") == truth == Dseq("AT", "CTATAG", 2).seguid()

    truth = "cdseguid=5fHMG19IbYxn7Yr7_sOCkvaaw7U"
    assert cdseguid("ACGTT", "AACGT") == truth == Dseq("ACGTT", "AACGT", circular=True).seguid()
    assert cdseguid("AACGT", "ACGTT") == truth == Dseq("AACGT", "ACGTT", circular=True).seguid()


def test_ovhg():
    # No overhang
    assert Dseq("AAAA").ovhg == 0
    assert Dseq("AAAA", circular=True).ovhg == 0
    # Sticky ends
    assert Dseq("FFAA").ovhg == 2
    assert Dseq("EEAA").ovhg == -2

    # Sticky end on the other hang does not matter
    assert Dseq("AAFF").ovhg == 0
    assert Dseq("AAEE").ovhg == 0

    #
    assert Dseq("FFAAFF").ovhg == 2
    assert Dseq("FFAAEE").ovhg == 2
    assert Dseq("EEAAEE").ovhg == -2
    assert Dseq("EEAAFF").ovhg == -2

    # Single strand
    assert Dseq("EEEE").ovhg is None
    assert Dseq("FFFF").ovhg is None


def test_watson_ovhg():
    # No overhang
    for seq in [
        "AAAA",
        "AAAA",
        "FFAA",
        "EEAA",
        "AAFF",
        "AAEE",
        "FFAAFF",
        "FFAAEE",
        "EEAAEE",
        "EEAAFF",
    ]:
        assert (
            Dseq(seq).watson_ovhg == Dseq(seq).reverse_complement().ovhg
        ), f"error for {seq}"

    # Single strand
    assert Dseq("EEEE").watson_ovhg is None
    assert Dseq("FFFF").watson_ovhg is None





def test_melt():

    assert Dseq("AGJGaGEg").melt(2) == (Dseq("EP"), Dseq("FQJGaGEp"), Dseq("q"))

    assert Dseq("AGIGaGFg").melt(2) == (Dseq("FQ"), Dseq("EPIGaGFq"), Dseq("p"))
    assert Dseq("AGJGaGFg").melt(2) == (Dseq("EP"), Dseq("FQJGaGFq"), Dseq("p"))

    assert Dseq("AGIGaGEg").melt(2) == (Dseq("FQ"), Dseq("EPIGaGEp"), Dseq("q"))
    assert Dseq("GATPGCPGCA").melt(2) == (Dseq("QJ"), Dseq("GATPPIPGCA"))
    assert Dseq("GATQGCQGCA").melt(2) == (Dseq("PI"), Dseq("GATQQJQGCA"))
    assert Dseq("PEXIGAQFZJ").melt(2) == (Dseq("PEXIPE"), Dseq("QFQFZJ"))
    assert Dseq("QFZJGAPEXI").melt(2) == (Dseq("QFZJQF"), Dseq("PEPEXI"))
    assert Dseq("AGJGaGEgGATC").melt(2) == (Dseq("EP"), Dseq("FQJGaGEgGATC"))
    assert Dseq("AGIGaGFgGATC").melt(2) == (Dseq("FQ"), Dseq("EPIGaGFgGATC"))
    assert Dseq("AGJGaGFgGATC").melt(2) == (Dseq("EP"), Dseq("FQJGaGFgGATC"))
    assert Dseq("AGIGaGEgGATC").melt(2) == (Dseq("FQ"), Dseq("EPIGaGEgGATC"))
    assert Dseq("GATPGGPGCAGATC").melt(2) == (Dseq("QQ"), Dseq("GATPPPPGCAGATC"))
    assert Dseq("GATQGGQGCAGATC").melt(2) == (Dseq("PP"), Dseq("GATQQQQGCAGATC"))
    assert Dseq("PEXIGAQFZJGATC").melt(2) == (Dseq("PEXIPE"), Dseq("QFQFZJGATC"))
    assert Dseq("QFZJGAPEXIGATC").melt(2) == (Dseq("QFZJQF"), Dseq("PEPEXIGATC"))
    assert Dseq("GATCAGJGaGEg").melt(2) == (Dseq("GATCAGJGaGEp"), Dseq("q"))
    assert Dseq("GATCAGIGaGFg").melt(2) == (Dseq("GATCAGIGaGFq"), Dseq("p"))
    assert Dseq("GATCAGJGaGFg").melt(2) == (Dseq("GATCAGJGaGFq"), Dseq("p"))
    assert Dseq("GATCAGIGaGEg").melt(2) == (Dseq("GATCAGIGaGEp"), Dseq("q"))
    assert Dseq("GATCGATPGGPGCA").melt(2) == (Dseq("QQ"), Dseq("GATCGATPPPPGCA"))
    assert Dseq("GATCGATQGGQGCA").melt(2) == (Dseq("PP"), Dseq("GATCGATQQQQGCA"))
    assert Dseq("GATCPEXIGAQFZJ").melt(2) == (Dseq("GATCPEXIPE"), Dseq("QFQFZJ"))
    assert Dseq("GATCQFZJGAPEXI").melt(2) == (Dseq("GATCQFZJQF"), Dseq("PEPEXI"))
    assert Dseq("GATCAGJGaGEgGATC").melt(2) == ()
    assert Dseq("GATCAGIGaGFgGATC").melt(2) == ()
    assert Dseq("GATCAGJGaGFgGATC").melt(2) == ()
    assert Dseq("GATCAGIGaGEgGATC").melt(2) == ()
    assert Dseq("GATCGATPGGPGCAGATC").melt(2) == (Dseq("QQ"), Dseq("GATCGATPPPPGCAGATC"))
    assert Dseq("GATCGATQGGQGCAGATC").melt(2) == (Dseq("PP"), Dseq("GATCGATQQQQGCAGATC"))
    assert Dseq("GATCPEXIGAQFZJGATC").melt(2) == (Dseq("GATCPEXIPE"), Dseq("QFQFZJGATC"))
    assert Dseq("GATCQFZJGAPEXIGATC").melt(2) == (Dseq("GATCQFZJQF"), Dseq("PEPEXIGATC"))
    assert Dseq("GATCPEXIGAQFZJGATC").melt(2) == (Dseq("GATCPEXIPE"), Dseq("QFQFZJGATC"))
    assert Dseq("GATCQFZJGAPEXIGATC").melt(2) == (Dseq("GATCQFZJQF"), Dseq("PEPEXIGATC"))


def test__get_ds_meltsites():

    assert Dseq("AGJGaGEg").get_ds_meltsites(2) == [((2, 2), None), ((8, 1), None)]
    assert Dseq("AGIGaGFg").get_ds_meltsites(2) == [((0, -2), None), ((7, -1), None)]
    assert Dseq("AGJGaGFg").get_ds_meltsites(2) == [((2, 2), None), ((7, -1), None)]
    assert Dseq("AGIGaGEg").get_ds_meltsites(2) == [((0, -2), None), ((8, 1), None)]

    assert Dseq("PEXIGAQFZJ").get_ds_meltsites(2) == [((6, 2), None)]
    assert Dseq("QFZJGAPEXI").get_ds_meltsites(2) == [((4, -2), None)]

    assert Dseq("GATCPEXIGAQFZJGATC").get_ds_meltsites(2) == [((10, 2), None)]
    assert Dseq("GATCQFZJGAPEXIGATC").get_ds_meltsites(2) == [((8, -2), None)]

    assert Dseq("AGCPAGQGAT", circular=True).get_ds_meltsites(2) == [((6, 2), None)]
    assert Dseq("AGCQAGPGAT", circular=True).get_ds_meltsites(2) == [((4, -2), None)]

def test__get_ds_meltsites():

    assert Dseq("AGJGaGEg").get_ds_meltsites(2) == [((2, 2), None), ((8, 1), None)]
    assert Dseq("AGIGaGFg").get_ds_meltsites(2) == [((0, -2), None), ((7, -1), None)]
    assert Dseq("AGJGaGFg").get_ds_meltsites(2) == [((2, 2), None), ((7, -1), None)]
    assert Dseq("AGIGaGEg").get_ds_meltsites(2) == [((0, -2), None), ((8, 1), None)]

    assert Dseq("PEXIGAQFZJ").get_ds_meltsites(2) == [((6, 2), None)]
    assert Dseq("QFZJGAPEXI").get_ds_meltsites(2) == [((4, -2), None)]

    assert Dseq("GATCPEXIGAQFZJGATC").get_ds_meltsites(2) == [((10, 2), None)]
    assert Dseq("GATCQFZJGAPEXIGATC").get_ds_meltsites(2) == [((8, -2), None)]

    assert Dseq("AGCPAGQGAT", circular=True).get_ds_meltsites(2) == [((6, 2), None)]
    assert Dseq("AGCQAGPGAT", circular=True).get_ds_meltsites(2) == [((4, -2), None)]


def test_nibble():

    assert Dseq("pgatc").nibble_five_prime_left() == Dseq("gatc")
    assert Dseq("pgatc").nibble_five_prime_left(2) == Dseq("qatc")
    assert Dseq("pgatc").nibble_three_prime_left() == Dseq("ppatc")
    assert Dseq("pgatc").nibble_five_prime_right() == Dseq("pgati")
    assert Dseq("pgatc").nibble_three_prime_right() == Dseq("pgatj")

    assert Dseq("qgatc").nibble_five_prime_left() == Dseq("qqatc")
    assert Dseq("qgatc").nibble_five_prime_left(2) == Dseq("qqftc")
    assert Dseq("qgatc").nibble_three_prime_left() == Dseq("gatc")
    assert Dseq("qgatc").nibble_five_prime_right() == Dseq("qgati")
    assert Dseq("qgatc").nibble_three_prime_right() == Dseq("qgatj")

    # Dseq(-2)
    # 5-gg-3
    #   ||
    # 3-  -5
    assert Dseq("pp").nibble_five_prime_left() == Dseq("p")
    assert Dseq("pp").nibble_three_prime_left() == Dseq("pp")
    assert Dseq("pp").nibble_five_prime_right() == Dseq("pp")
    assert Dseq("pp").nibble_three_prime_right() == Dseq("p")

    # Dseq(-2)
    # 5-  -3
    #   ||
    # 3-cc-5
    assert Dseq("qq").nibble_five_prime_left() == Dseq("qq")
    assert Dseq("qq").nibble_three_prime_left() == Dseq("q")
    assert Dseq("qq").nibble_five_prime_right() == Dseq("q")
    assert Dseq("qq").nibble_three_prime_right() == Dseq("qq")

@pytest.mark.xfail(reason="Not implemented.")
def test_anneal():

    # Dseq(-9)
    # GGATC   G
    # C   GTAGC

    assert Dseq("GPEXI") / Dseq("JFZJG") == Dseq("GPEXCFZJG")

    with pytest.raises(TypeError):
        Dseq("GPEXI") + Dseq("JFZJG")

    # Dseq(-9)
    # GGATC   G
    # C   GTAGC
    assert Dseq("GPEXI") / Dseq("JFZJG") == Dseq("GPEXCFZJG")

    # GGACT       G
    # C       GTAGC
    assert Dseq("GPEIX") / Dseq("JFZJG") == None

    # Dseq(-8)
    # GGACT  G
    # C  GATGC
    assert Dseq("GPEIX") / Dseq("JZFJG") == Dseq("GPECTFJG")

    # Dseq(-8)
    # GGACTA G
    # C  GATGC
    Dseq("GPEIXE") / Dseq("JZFJG") == Dseq("GPECTAJG")

    # Dseq(-8)
    # GGACTA G
    # C TGATGC
    assert Dseq("GPEIXE") / Dseq("FJZFJG") == Dseq("GPACTAJG")

    # Dseq(-8)
    # GGACTA G
    # C  GATGC
    assert Dseq("GPEIXE") / Dseq("JZFJG") == Dseq("GPECTAJG")

    # Dseq(-8)
    # GGACTA G
    # C TGATGC
    assert Dseq("GPEIXE") / Dseq("FJZFJG") == Dseq("GPACTAJG")

    # Dseq(-8)
    # GGACTACG
    # C TGATGC
    assert Dseq("GPEIXEI") / Dseq("FJZFJG") == Dseq("GPACTACG")

    # Dseq(-8)
    # GGACTACG
    # CCTGATGC
    assert Dseq("GPEIXEI") / Dseq("QFJZFJG") == Dseq("GGACTACG")


def test_mw():

    from Bio.Data.IUPACData import unambiguous_dna_weights
    from Bio.Data.IUPACData import unambiguous_rna_weights
    from Bio.Data.IUPACData import atom_weights

    # The molecular weight values for a short DNA molecule agrees very well
    # with a the online tool https://molbiotools.com/dnacalculator.php (*)
    # accessed December 20, 2025

    double_strand_linear = Dseq("GATTACA")

    assert round(double_strand_linear.mw(), 1) == 4359.8 # 4359.81 Da (*)

    double_strand_circular = Dseq("GATTACA", circular = True)

    assert round(double_strand_circular.mw(), 1) == 4323.8 # 4323.78 Da (*)

    single_strand_linear = Dseq("PEXXEIE")

    assert round(single_strand_linear.mw(), 1) == 2184.4 # 2184.41 Da (*)

    single_strand_circular = Dseq("PEXXEIE", circular = True)

    assert round(single_strand_circular.mw(), 1) == 2166.4 # 2166.39 Da (*)

    ds_lin_obj2 = Dseq("GATZFCA")
    assert repr(ds_lin_obj2) == textwrap.dedent("""\
                                Dseq(-7)
                                GAT  CA
                                CTAATGT""")
    # ds_lin_obj2 is missing a T and an A compared to double_strand_linear
    mw = round(ds_lin_obj2.mw(), 1)

    h2o = atom_weights["H"]*2 + atom_weights["O"]

    T = unambiguous_dna_weights["T"]
    A = unambiguous_dna_weights["A"]

    assert round(double_strand_linear.mw() - A - T + h2o, 1) == round(mw, 1)


def test_misc2():

    target = Dseq("GGTCTCtggatccccc", circular=True)
    targetrc = target.rc()

    assert target.cut(BsaI) == tuple([x.rc() for x in targetrc.cut(BsaI)])

    target = Dseq("GGTCTCtGGATCCaCACGTCtGACGTCcGGTACC", circular=True)
    targetrc = target.rc()

    assert target.cut(BsaI) == tuple([x.rc() for x in targetrc.cut(BsaI)])

    target = Dseq("GGTCTCtttCCCCttttGAGACC")
    targetrc = target.rc()

    assert target.cut(BsaI) == tuple([x.rc() for x in targetrc.cut(BsaI)])[::-1]

    puc19mcs = Dseq("gaattcgagctcggtacccggggatcctctagagtcgacctgcaggcatgcaagctt")


    #           10        20        30        40        50
    #       *    *    *    *    *    *    *    *    *    *    *
    #                   17 XmaI
    #               13 KpnI  22 BamHI          40 PstI
    #         7 SacI    17 SmaI          34 SalI     46 SphI
    #   1 EcoRI     13 Acc65I      28 XbaI    39 SbfI      52 HindIII
    #   |     |     |   |    |     |     |    ||     |     |
    # 1 gaattcgagctcggtacccggggatcctctagagtcgacctgcaggcatgcaagctt 57
    #   E  F  E  L  G  T  R  G  S  S  R  V  D  L  Q  A  C  K  L
    #    N  S  S  S  V  P  G  D  P  L  E  S  T  C  R  H  A  S
    #     I  R  A  R  Y  P  G  I  L  *  S  R  P  A  G  M  Q  A
    #   cttaagctcgagccatgggcccctaggagatctcagctggacgtccgtacgttcgaa

    assert puc19mcs.cut(SmaI) == (
    Dseq.from_representation("""
    Dseq(-19)
    gaattcgagctcggtaccc
    cttaagctcgagccatggg
    """),
    Dseq.from_representation("""
    Dseq(-38)
    ggggatcctctagagtcgacctgcaggcatgcaagctt
    cccctaggagatctcagctggacgtccgtacgttcgaa"""))

    puc19mcs.cut(BamHI) == (
    Dseq.from_representation("""
    Dseq(-19)
    gaattcgagctcggtacccggg
    cttaagctcgagccatgggcccctag
    """),
    Dseq.from_representation("""
    Dseq(-38)
    gatcctctagagtcgacctgcaggcatgcaagctt
        gagatctcagctggacgtccgtacgttcgaa
    """))

    puc19mcs.cut(SmaI, BamHI) == (
    Dseq.from_representation("""
    Dseq(-19)
    gaattcgagctcggtaccc
    cttaagctcgagccatggg
    """),
    Dseq.from_representation("""
    Dseq(-7)
    ggg
    cccctag
    """),
    Dseq.from_representation("""
    gatcctctagagtcgacctgcaggcatgcaagctt
        gagatctcagctggacgtccgtacgttcgaa
    """))

    with pytest.raises(ValueError) as err:
        puc19mcs.cut(Acc65I, KpnI)
    assert err.match("overlap")

@pytest.mark.xfail(reason="Can not yet deal with uracil.")
def test_seguid_with_uracil():

    x = Dseq("PEXIaUaOaQFZJ")

    assert x.seguid() == "ldseguid=AA-pG6wsGMu2J-AuWBWlOi3iDc4"

    y = Dseq("QFZJaUaOaPEXI")

    assert y.seguid() == "ldseguid=10e3S-WXX1z0BHqpbYTPeb2MwUY"

    # Dseq(-5)
    # CUGOT
    # GACOA

# def test_pair():

#     from pydna.dseq import Dseq, pair

#     x = Dseq("PEXIaUaOaQFZJ")

#     x_wco = Dseq(pair("GATCaUaAa", "tAtUtCTAG"[::-1], -4))

#     x_wco2 = Dseq(pair(x.watson, x.crick, x.ovhg))

#     assert x == x_wco == x_wco2

#     y = Dseq("QFZJaUaOaPEXI")

#     y_wco = Dseq(pair("aUaAaGATC", "CTAGtAtUt"[::-1], 4))

#     y_wco2 = Dseq(pair(y.watson, y.crick, y.ovhg))

#     assert y == y_wco == y_wco2

#     assert pair("GATT", "CTAA"[::-1]) == b"GATT"

#     assert pair("GATTa", "aCTAA"[::-1]) == b"zGATTe"

#     assert pair("GATCaUaAa", "tAtUtCTAG"[::-1], -4) == b"PEXIaUaOaQFZJ"

def test_dscode():

    assert str(Dseq("PEXIaaaQFZJ")) == "GATCaaaGATC"
    assert str(Dseq("QFZJaaaPEXI")) == "GATCaaaGATC"

def cut2():
    from pydna.dseq import Dseq
    from Bio.Restriction import BamHI, AjiI, ZraI, KpnI, Acc65I, BsaI

    target = Dseq("GGTCTCtGGATCCaCACGTCtGACGTCcGGTACC")
    targetrc = target.rc()

    for enz in (BamHI, AjiI, ZraI, KpnI, Acc65I, BsaI):
        assert target.cut(enz) == tuple([x.rc() for x in targetrc.cut(enz)])[::-1]
















def test_cut_subtypes():
    from pydna.dseq import Dseq

    # http://rebase.neb.com/rebase/sublist.html

    # from Bio.Restriction import HaeIV  #  ImportError: cannot import name 'HaeIV' from 'Bio.Restriction'
    from Bio.Restriction import AciI
    from Bio.Restriction import AhdI
    from Bio.Restriction import BcgI
    from Bio.Restriction import Bpu10I
    from Bio.Restriction import BsgI
    from Bio.Restriction import BslI
    from Bio.Restriction import DpnI
    from Bio.Restriction import Eco57I
    from Bio.Restriction import EcoRI
    from Bio.Restriction import EcoRII
    from Bio.Restriction import FokI
    from Bio.Restriction import GsuI
    from Bio.Restriction import MmeI
    from Bio.Restriction import NaeI
    from Bio.Restriction import PpuMI
    from Bio.Restriction import SfiI
    from Bio.Restriction import SgrAI

    enz = (
        AciI,
        AhdI,
        BcgI,
        Bpu10I,
        BsgI,
        BslI,
        DpnI,
        Eco57I,
        EcoRI,
        EcoRII,
        FokI,
        GsuI,
        MmeI,
        NaeI,
        PpuMI,
        SfiI,
        SgrAI,
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtACCGCGtcaagtgtctatgcgtagatcgtA").cut(AciI) == (
        Dseq("GtcaagtgtctatgcgtagatcgtACJQ"),
        Dseq("IPCGtcaagtgtctatgcgtagatcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtAGACAAAAAGTCGtcaagtgtctatgcgtagatcgtA").cut(AhdI) == (
        Dseq("GtcaagtgtctatgcgtagatcgtAGACAAE"),
        Dseq("FAAGTCGtcaagtgtctatgcgtagatcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtACGAAAAAAATGCGtcaagtgtctatgcgtagatcgtA").cut(BcgI) == (
        Dseq("Gtcaagtgtctatpi"),
        Dseq("qjgtagatcgtACGAAAAAAATGCGtcaagtgtcxe"),
        Dseq("zftgcgtagatcgtA"),
    )
    assert Dseq("GtcaagtgtctatgcgtagatcgtACGAAAAAAATGCGtcaagtgtctatgcgtagatcgtA").rc().cut(BcgI) == (
        Dseq("Tacgatctacgcaxe"),
        Dseq("zfgacacttgaCGCATTTTTTTCGTacgatctacpi"),
        Dseq("qjatagacacttgaC"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtACCTAAGCGtcaagtgtctatgcgtagatcgtA").cut(Bpu10I) == (
        Dseq("GtcaagtgtctatgcgtagatcgtACCZFF"),
        Dseq("XEEGCGtcaagtgtctatgcgtagatcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtAGTGCAGGtcaagtgtctatgcgtagatcgtA").cut(BsgI) == (
        Dseq("GtcaagtgtctatgcgtagatcgtAGTGCAGGtcaagtgtctatgip"),
        Dseq("jqtagatcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtACCAAAAAAAGGGtcaagtgtctatgcgtagatcgtA").cut(BslI) == (
        Dseq("GtcaagtgtctatgcgtagatcgtACCAAEEE"),
        Dseq("FFFAAGGGtcaagtgtctatgcgtagatcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtAGATCGtcaagtgtctatgcgtagatcgtA").cut(DpnI) == (
        Dseq("Gtcaagtgtctatgcgtaga"),
        Dseq("tcgtAGA"),
        Dseq("TCGtcaagtgtctatgcgtaga"),
        Dseq("tcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtACTGAAGGtcaagtgtctatgcgtagatcgtA").cut(Eco57I) == (
        Dseq("GtcaagtgtctatgcgtagatcgtACTGAAGGtcaagtgtctatgip"),
        Dseq("jqtagatcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtAGAATTCGtcaagtgtctatgcgtagatcgtA").cut(EcoRI) == (
        Dseq("GtcaagtgtctatgcgtagatcgtAGFFZZ"),
        Dseq("EEXXCGtcaagtgtctatgcgtagatcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtACCAGGGtcaagtgtctatgcgtagatcgtA").cut(EcoRII) == (
        Dseq("GtcaagtgtctatgcgtagatcgtAJJFQQ"),
        Dseq("IIEPPGtcaagtgtctatgcgtagatcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtAGGATGGtcaagtgtctatgcgtagatcgtA").cut(FokI) == (
        Dseq("GtcaagtgtctatgcgtagatcgtAGGATGGtcaagtgtjzfz"),
        Dseq("ixexgcgtagatcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtACTGGAGGtcaagtgtctatgcgtagatcgtA").cut(GsuI) == (
        Dseq("GtcaagtgtctatgcgtagatcgtACTGGAGGtcaagtgtctatgip"),
        Dseq("jqtagatcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtATCCAACGtcaagtgtctatgcgtagatcgtA").cut(MmeI) == (
        Dseq("GtcaagtgtctatgcgtagatcgtATCCAACGtcaagtgtctatgcgtape"),
        Dseq("qftcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtAGCCGGCGtcaagtgtctatgcgtagatcgtA").cut(NaeI) == (
        Dseq("GtcaagtgtctatgcgtagatcgtAGCC"),
        Dseq("GGCGtcaagtgtctatgcgtagatcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtAAGGACCTGtcaagtgtctatgcgtagatcgtA").cut(PpuMI) == (
        Dseq("GtcaagtgtctatgcgtagatcgtAAGQFJ"),
        Dseq("PEICTGtcaagtgtctatgcgtagatcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtAGGCCAAAAAGGCCGtcaagtgtctatgcgtagatcgtA").cut(SfiI) == (
        Dseq("GtcaagtgtctatgcgtagatcgtAGGCCAEEE"),
        Dseq("FFFAGGCCGtcaagtgtctatgcgtagatcgtA"),
    )

    assert Dseq("GtcaagtgtctatgcgtagatcgtACACCGGTGGtcaagtgtctatgcgtagatcgtA").cut(SgrAI) == (
        Dseq("GtcaagtgtctatgcgtagatcgtACAJJQQ"),
        Dseq("IIPPTGGtcaagtgtctatgcgtagatcgtA"),
    )


def test_bcgi_circular():
    from pydna.dseq import Dseq
    from Bio.Restriction import BcgI

    t = "GAGttttttttttCGAaaaaaaTGCttttttttttTCT"
    BcgI.search(Dseq(t))

    # GAGttttttttttCGAaaaaaaTGCttttttttttTCT
    # CTCaaaaaaaaaaGCTttttttACGaaaaaaaaaaAGA

    # GAG   tttt..ttttTC   T
    # C   TCaaaa..aaaa   AGA
    for topology in [False]:
        obj = Dseq(t, circular=topology)
        for e in [BcgI]:
            for i in range(len(obj) + 1):
                cs = obj.get_cutsites(e)
                assert cs == [((3, 2), BcgI), ((37, 2), BcgI)]




def test_cut_circular2():
    from pydna.dseq import Dseq
    from Bio.Restriction import BsaI, KpnI, Acc65I, NotI

    test = "aaaaaaGGTACCggtctcaaaa"

    for i in range(len(test)):
        nt = test[i:] + test[:i]

        a = Dseq(nt, circular=True).cut(Acc65I)[0]  # G^GTACC
        assert a.watson.upper() == "GTACCGGTCTCAAAAAAAAAAG"
        assert a.crick.upper() == "GTACCTTTTTTTTTTGAGACCG"
        assert a.ovhg == -4  # CggtctcaaaaaaaaaaGGTAC
        b = Dseq(nt, circular=True).cut(KpnI)[0]  # GGTAC^C
        assert b.watson.upper() == "CGGTCTCAAAAAAAAAAGGTAC"
        assert b.crick.upper() == "CTTTTTTTTTTGAGACCGGTAC"
        assert b.ovhg == 4
        c = Dseq(nt, circular=True).cut(BsaI)[0]  # ggtctcnnn
        assert c.watson.upper() == "AAAAAAAAAGGTACCGGTCTCA"
        assert c.crick.upper() == "TTTTTGAGACCGGTACCTTTTT"
        assert c.ovhg == -4
        d = Dseq(nt, circular=True).cut(NotI)
        assert d == ()

    from pydna.dseq import Dseq

    # from pydna.dseq_old import oDseq

    a = Dseq("gattcgtatgctgatcgtacgtactgaaaac")

    assert repr(a) == "Dseq(-31)\ngatt..aaac\nctaa..tttg"

    b = Dseq("pexxipxexpctgatcgtacgtactgaaaac")

    assert repr(b) == "Dseq(-31)\ngattcgtatgctga..aaac\n          gact..tttg"

    c = Dseq("pexxipxexpitgatcgtacgtactgaaaac")

    assert repr(c) == "Dseq(-31)\ngatt..atgctgat..aaac\n          acta..tttg"

    d = Dseq("pexxipxexpctgatcgtacg")

    assert repr(d) == "Dseq(-21)\ngattcgtatgctgatcgtacg\n          gactagcatgc"

    e = Dseq("jzffqjfzfjgactagcatgcatgacttttc")

    assert repr(e) == "Dseq(-31)\n          gact..tttc\ngattcgtatgctga..aaag"

    f = Dseq("JzffqjfzfCgactagcatgcatgacttttc")

    assert repr(f) == "Dseq(-31)\n         Cgac..tttc\nGattcgtatGctg..aaag"

    g = Dseq("gattcgtatgctgatcgxeipxeixpeeeei")

    assert repr(g) == "Dseq(-31)\ngatt..atcgtacg..aaac\nctaa..tagc"

    h = Dseq("cgtatgctgatcgxeipxeixpeeeei")

    assert repr(h) == "Dseq(-27)\ncgtatgctgatcgtacgtactgaaaac\ngcatacgactagc"

    i = Dseq("cgtatgctgatcgxeipxeixpeeeeiepeix")

    assert repr(i) == "Dseq(-32)\ncgta..atcgtacg..gact\ngcat..tagc"

    j = Dseq("jfZJZJZJJZZgaxxipxexpixpexipxeipxeixpeeeei")

    assert repr(j) == "Dseq(-42)\n          gattcg..aaac\ngtAG..GGAAct"

    k = Dseq("gzzzzjfqzfjqzfjqfzjfqjfzfjqffzj")

    assert repr(k) == "Dseq(-31)\ng\ncaaaa..ttag"

    x = Dseq("gattcgtatgctgatcgtacgtactgaaaa")

    assert repr(x) == "Dseq(-30)\ngattcgtatgctgatcgtacgtactgaaaa\nctaagcatacgactagcatgcatgactttt"

    y = Dseq("pexxipxexpctgatcgtacgtactgaaaa")

    assert repr(y) == "Dseq(-30)\ngattcgtatgctgatcgtacgtactgaaaa\n          gactagcatgcatgactttt"

    z = Dseq("pexxipxexpitgatcgtacgtactgaaaa")

    assert repr(z) == "Dseq(-30)\ngattcgtatgctgatcgtacgtactgaaaa\n           actagcatgcatgactttt"


@pytest.mark.xfail(reason="Shifting circular sequences may fail with missing nucleotides.")
def test_shifted2():
    from pydna.dseq import Dseq

    a = Dseq("gatc", circular=True)

    assert a.shifted(1) == Dseq("atcg", circular=True)

    assert a.shifted(4) == a

    b = Dseq("gatc", circular=False)
    with pytest.raises(TypeError):
        b.shifted(1)

    # Shifted with zero gives a copy of the sequence, not the same sequence
    assert a.shifted(0) == a
    assert a.shifted(0) is not a

    x = Dseq("TPGJG", circular=True)

    """
    TGG G
    A CGC

    GG GT
     CGCA

    G GTG
    CGCA

     GTGG
    GCA C

    GTGG
    CA CG

    TGG G
    A CGC
    """

    for i in range(len(x)):
        assert x.shifted(i)._data == x._data[i:] + x._data[:i]















@pytest.mark.xfail(reason="Fail with repeated enzymes")
def test_overlapping_cuts():

    import pytest
    from pydna.dseq import Dseq
    from Bio.Restriction import PstI, SbfI, BamHI, DpnI, BcgI

    s = Dseq("CCTGCAGG")

    assert s.cut(PstI) == (Dseq(b"CCXPIE"), Dseq(b"ZQJFGG"))
    assert s.cut(SbfI) == (Dseq(b"CCXPIE"), Dseq(b"ZQJFGG"))

    with pytest.raises(ValueError):
        s.cut(PstI, SbfI)

    t = Dseq("GGATCC")

    assert t.cut(BamHI) == (Dseq("GQFZJ"), Dseq("PEXIC"))
    assert t.cut(DpnI) == (Dseq("GGA"), Dseq("TCC"))

    with pytest.raises(ValueError):
        t.cut(BamHI, DpnI)

    q = Dseq("TgcgtagatcgtACGAggatccTGCGtcaagtgtctat")

    with pytest.raises(ValueError):
        q.cut(BamHI, BamHI, BamHI)

    assert q.cut(BamHI, BcgI) == (Dseq("Tpi"), Dseq("qjgtagatcgtACGAgqfzj"), Dseq("pexicTGCGtcaagtgtcxe"), Dseq("zft"))














def test_new():
    from pydna.dseq import Dseq
    from Bio.Restriction import KpnI, Acc65I, BsaI, XmaI, SmaI, BamHI

    fiveoh = Dseq("PEXIaaaQFZJ")
    assert str(fiveoh + fiveoh) == "GATCaaaGATCaaaGATC"
    threeoh = Dseq("QFZJQtttPEXIP")
    assert str(threeoh + threeoh) == "GATCGtttGATCGtttGATCG"

    assert repr(Dseq("AIXEP") + Dseq("JZFQA")) == "Dseq(-6)\nACTAGA\nTGATCT"

    with pytest.raises(TypeError):
        # assert repr(Dseq("AIXEP") + Dseq("ZFQA")) == "Dseq(-6)\nACTAGA\nT ATCT"
        Dseq("AIXEP") + Dseq("ZFQA")

    with pytest.raises(TypeError):
        # assert repr(Dseq("AIXE") + Dseq("JZFQA")) == "Dseq(-6)\nACTA A\nTGATCT"
        Dseq("AIXE") + Dseq("JZFQA")

    assert repr(Dseq("AP") + Dseq("QA")) == "Dseq(-3)\nAGA\nTCT"
    assert repr(Dseq("AE") + Dseq("FA")) == "Dseq(-3)\nAAA\nTTT"
    assert repr(Dseq("AX") + Dseq("ZA")) == "Dseq(-3)\nATA\nTAT"
    assert repr(Dseq("AI") + Dseq("JA")) == "Dseq(-3)\nACA\nTGT"

    with pytest.raises(TypeError):
        # assert repr(Dseq("APP") + Dseq("QA")) == "Dseq(-4)\nAGGA\nT CT"
        Dseq("APP") + Dseq("QA")

    with pytest.raises(TypeError):
        # assert repr(Dseq("AEE") + Dseq("FA")) == "Dseq(-4)\nAAAA\nT TT"
        Dseq("AEE") + Dseq("FA")

    with pytest.raises(TypeError):
        # assert repr(Dseq("AXX") + Dseq("ZA")) == "Dseq(-4)\nATTA\nT AT"
        Dseq("AXX") + Dseq("ZA")

    with pytest.raises(TypeError):
        # assert repr(Dseq("AII") + Dseq("JA")) == "Dseq(-4)\nACCA\nT GT"
        Dseq("AXX") + Dseq("ZA")

    with pytest.raises(TypeError):
        # assert repr(Dseq("AP") + Dseq("QQA")) == "Dseq(-4)\nAG A\nTCCT"
        Dseq("AP") + Dseq("QQA")

    with pytest.raises(TypeError):
        # assert repr(Dseq("AE") + Dseq("FFA")) == "Dseq(-4)\nAA A\nTTTT"
        Dseq("AE") + Dseq("FFA")

    with pytest.raises(TypeError):
        # assert repr(Dseq("AX") + Dseq("ZZA")) == "Dseq(-4)\nAT A\nTAAT"
        Dseq("AX") + Dseq("ZZA")

    with pytest.raises(TypeError):
        # assert repr(Dseq("AI") + Dseq("JJA")) == "Dseq(-4)\nAC A\nTGGT"
        Dseq("AI") + Dseq("JJA")

    assert Dseq("QAP").looped()[:] == Dseq("GA")
    assert Dseq("PAQ").looped()[:] == Dseq("GA")
    assert Dseq("EAF").looped()[:] == Dseq("AA")
    assert Dseq("FAE").looped()[:] == Dseq("AA")
    assert Dseq("XAZ").looped()[:] == Dseq("TA")
    assert Dseq("ZAX").looped()[:] == Dseq("TA")
    assert Dseq("IAJ").looped()[:] == Dseq("CA")
    assert Dseq("JAI").looped()[:] == Dseq("CA")

    s = Dseq("PEXIAAAQFZJ")

    assert s.looped()._data == b"GATCAAA"

    s = Dseq("PEXIAAAQFZJ")

    # assert s._fill_in_right("gatc") == "GATCAAAQFZJ"
    # assert s._fill_in_right("gat") == "PATCAAAQFZJ"
    # assert s._fill_in_right("ga") == "PETCAAAQFZJ"
    # assert s._fill_in_right("g") == "PEXCAAAQFZJ"
    # assert s._fill_in_right("") == s._data

    # assert s._fill_in_right("atc")._data == s._data
    # assert s._fill_in_right("at")._data == s._data
    # assert s._fill_in_right("a")._data == s._data
    # assert s._fill_in_right("")._data == s._data

    # assert s._fill_in_left("gatc")._data == b"PEXIAAAGATC"
    # assert s._fill_in_left("gat")._data == b"PEXIAAAGATJ"
    # assert s._fill_in_left("ga")._data == b"PEXIAAAGAZJ"
    # assert s._fill_in_left("g")._data == b"PEXIAAAGFZJ"
    # assert s._fill_in_left("")._data == s._data

    # assert s._fill_in_left("atc")._data == s._data
    # assert s._fill_in_left("at")._data == s._data
    # assert s._fill_in_left("a")._data == s._data
    # assert s._fill_in_left("")._data == s._data

    # assert s.fill_in("gatc")._data == b"GATCAAAGATC"
    # assert s.fill_in("gat")._data == b"PATCAAAGATJ"
    # assert s.fill_in("ga")._data == b"PETCAAAGAZJ"
    # assert s.fill_in("g")._data == b"PEXCAAAGFZJ"
    # assert s.fill_in("")._data == s._data

    assert s == Dseq.from_representation(repr(s))

    assert Dseq.from_full_sequence_and_overhangs("AAAAAA", crick_ovhg=2, watson_ovhg=2)._data == Dseq("FFAAEE")._data
    assert Dseq.from_full_sequence_and_overhangs("AAAAAA", crick_ovhg=-2, watson_ovhg=2)._data == Dseq("EEAAEE")._data
    assert Dseq.from_full_sequence_and_overhangs("AAAAAA", crick_ovhg=2, watson_ovhg=-2)._data == Dseq("FFAAFF")._data
    assert Dseq.from_full_sequence_and_overhangs("AAAAAA", crick_ovhg=-2, watson_ovhg=-2)._data == Dseq("EEAAFF")._data

    assert repr(Dseq("XIpexiPEXIpexiAqfzjQFZJqfzjQF")) == (
        "Dseq(-29)\n" "TCgatcGATCgatcA\n" "              TctagCTAGctagCT"
    )

    assert repr(Dseq("XIpexiPEXIpexiAAqfzjQFZJqfzjQF")) == (
        "Dseq(-30)\n" "TCgatcGATCgatcAA\n" "              TTctagCTAGctagCT"
    )

    assert repr(Dseq("EXIpexiPEXIpexiAqfzjQFZJqfzjQFZ")) == ("Dseq(-31)\n" "ATCg..gatcA\n" "          Tctag..gCTA")

    assert repr(Dseq("XIpexiPEXIpexiAaAqfzjQFZJqfzjQF")) == ("Dseq(-31)\n" "TCga..gatcAaA\n" "          TtTctag..agCT")

    dsdna = """
               Dseq(
                GGATCC
               aCCTAGGg"""

    assert Dseq.from_representation(dsdna)._data == b"zGGATCCj"

    dsdna = """
               Dseq(
               aGGATCCg
                CCTAGG"""

    assert Dseq.from_representation(dsdna)._data == b"eGGATCCp"

    dsdna = """
               Dseq(
                 GGATCCaa
               cccctagg"""

    assert Dseq.from_representation(dsdna)._data == b"qqGGATCCee"

    dsdna = """
               Dseq(
                 ggtaccaa
               ccccatgg"""

    assert Dseq.from_representation(dsdna)._data == b"qqggtaccee"

    assert Dseq("AP") + Dseq("") == Dseq("AP")
    assert Dseq("AE") + Dseq("") == Dseq("AE")
    assert Dseq("AX") + Dseq("") == Dseq("AX")
    assert Dseq("AI") + Dseq("") == Dseq("AI")
    assert Dseq("Ap") + Dseq("") == Dseq("Ap")
    assert Dseq("Ax") + Dseq("") == Dseq("Ax")
    assert Dseq("Ai") + Dseq("") == Dseq("Ai")

    assert Dseq("AQ") + Dseq("") == Dseq("AQ")
    assert Dseq("AF") + Dseq("") == Dseq("AF")
    assert Dseq("AZ") + Dseq("") == Dseq("AZ")
    assert Dseq("AJ") + Dseq("") == Dseq("AJ")
    assert Dseq("Aq") + Dseq("") == Dseq("Aq")
    assert Dseq("Af") + Dseq("") == Dseq("Af")
    assert Dseq("Az") + Dseq("") == Dseq("Az")
    assert Dseq("Ai") + Dseq("") == Dseq("Ai")

    # assert Dseq("PEXIGGATCCQFZJ").T4(b"CG") == Dseq(b"PEXCGGATCCGFZJ")
    # assert Dseq("PAGAJ").T4(b"ATG") == Dseq(b"PAGAJ")

    from pydna import _PydnaWarning

    # with pytest.warns(_PydnaWarning):
    #     assert Dseq("QFZJGGATCCPEXI").T4("") == Dseq("")

    # assert Dseq("pexiAqfzj").T4(b"gatc") == Dseq("gatcAgatc")

    assert Dseq("qqGGATCCee").seguid() == "ldseguid=F0z-LxHZqAK3HvqQiqjM7A28daE"

    ##################################################################################

    from Bio.SeqFeature import SimpleLocation, ExactPosition, CompoundLocation

    assert Dseq("ggtctcAAgcTT", circular=False).get_cutsites(BsaI) == [((7, -4), BsaI)]
    # assert Dseq("TggtctcAAgcT", circular=False).get_cutsites(BsaI) == [((8, -4), BsaI)]

    assert Dseq("TggtctcAAgcT", circular=False).get_cutsites(BsaI) == []
    assert Dseq("TTggtctcAAgc", circular=False).get_cutsites(BsaI) == []
    assert Dseq("cTTggtctcAAg", circular=False).get_cutsites(BsaI) == []
    assert Dseq("gcTTggtctcAA", circular=False).get_cutsites(BsaI) == []

    assert Dseq("gtctcAAgcTTg", circular=True).get_cutsites(BsaI) == [((6, -4), BsaI)]
    assert Dseq("ggtctcAAgcTT", circular=True).get_cutsites(BsaI) == [((7, -4), BsaI)]
    assert Dseq("TggtctcAAgcT", circular=True).get_cutsites(BsaI) == [((8, -4), BsaI)]
    assert Dseq("TTggtctcAAgc", circular=True).get_cutsites(BsaI) == [((9, -4), BsaI)]
    assert Dseq("cTTggtctcAAg", circular=True).get_cutsites(BsaI) == [((10, -4), BsaI)]
    assert Dseq("gcTTggtctcAA", circular=True).get_cutsites(BsaI) == [((11, -4), BsaI)]
    assert Dseq("AgcTTggtctcA", circular=True).get_cutsites(BsaI) == [((0, -4), BsaI)]
    assert Dseq("AAgcTTggtctc", circular=True).get_cutsites(BsaI) == [((1, -4), BsaI)]
    assert Dseq("cAAgcTTggtct", circular=True).get_cutsites(BsaI) == [((2, -4), BsaI)]
    assert Dseq("tcAAgcTTggtc", circular=True).get_cutsites(BsaI) == [((3, -4), BsaI)]
    assert Dseq("ctcAAgcTTggt", circular=True).get_cutsites(BsaI) == [((4, -4), BsaI)]
    assert Dseq("tctcAAgcTTgg", circular=True).get_cutsites(BsaI) == [((5, -4), BsaI)]
    assert Dseq("gtctcAAgcTTg", circular=True).get_cutsites(BsaI) == [((6, -4), BsaI)]

    assert Dseq("gcTTggtctcAA", circular=True).get_cutsites(BsaI)

    assert Dseq("gtctcAAgcTTggtctcAA", circular=True).get_cutsites(BsaI)
    assert Dseq("gtctcAAgcTTggtctcAA", circular=True).get_cutsites(BsaI)

    assert Dseq("ggatcc").cut(BamHI) == (Dseq("gqfzj"), Dseq("pexic"))
    assert Dseq("GGATCC").cut(BamHI) == (Dseq("GQFZJ"), Dseq("PEXIC"))
    assert Dseq("GGATCc").cut(BamHI) == (Dseq("GQFZJ"), Dseq("PEXIc"))
    assert Dseq("gGATCC").cut(BamHI) == (Dseq("gQFZJ"), Dseq("PEXIC"))
    assert Dseq("gGATCc").cut(BamHI) == (Dseq("gQFZJ"), Dseq("PEXIc"))

    assert Dseq("ggtacc").cut(Acc65I) == (Dseq("gqzfj"), Dseq("pxeic"))
    assert Dseq("GGTACC").cut(Acc65I) == (Dseq("GQZFJ"), Dseq("PXEIC"))
    assert Dseq("GGTACc").cut(Acc65I) == (Dseq("GQZFJ"), Dseq("PXEIc"))
    assert Dseq("gGTACC").cut(Acc65I) == (Dseq("gQZFJ"), Dseq("PXEIC"))
    assert Dseq("gGTACc").cut(Acc65I) == (Dseq("gQZFJ"), Dseq("PXEIc"))

    assert Dseq("ggtacc").cut(KpnI) == (Dseq("gpxei"), Dseq("qzfjc"))
    assert Dseq("GGTACC").cut(KpnI) == (Dseq("GPXEI"), Dseq("QZFJC"))
    assert Dseq("GGTACc").cut(KpnI) == (Dseq("GPXEI"), Dseq("QZFJc"))
    assert Dseq("gGTACC").cut(KpnI) == (Dseq("gPXEI"), Dseq("QZFJC"))
    assert Dseq("gGTACc").cut(KpnI) == (Dseq("gPXEI"), Dseq("QZFJc"))

    s = Dseq("CCCGGGGCATCGTAGTGATCGGTACC")

    a, b, c = s.cut([XmaI, Acc65I])

    assert a + b + c == s
    assert repr(a + b + c) == repr(s)

    a, b, c = s.cut([SmaI, Acc65I])

    assert a + b + c == s
    assert repr(a + b + c) == repr(s)

    a, b, c = s.cut([SmaI, KpnI])

    assert a + b + c == s
    assert repr(a + b + c) == repr(s)

    s = Dseq("PPCCCGGGGCATCGTAGTGATCGGTACC")

    a, b, c = s.cut([XmaI, Acc65I])

    assert a + b + c == s
    assert repr(a + b + c) == repr(s)

    s = Dseq("QQCCCGGGGCATCGTAGTGATCGGTACC")

    a, b, c = s.cut([XmaI, Acc65I])

    assert a + b + c == s
    assert repr(a + b + c) == repr(s)

    ##################################################################################

    assert Dseq("aa").get_cutsites(BamHI) == []

    assert Dseq("GGTACCa", circular=True).get_cutsites(Acc65I) == [((1, -4), Acc65I)]
    assert Dseq("GTACCaG", circular=True).get_cutsites(Acc65I) == [((0, -4), Acc65I)]
    assert Dseq("TACCaGG", circular=True).get_cutsites(Acc65I) == [((6, -4), Acc65I)]
    assert Dseq("ACCaGGT", circular=True).get_cutsites(Acc65I) == [((5, -4), Acc65I)]
    assert Dseq("CCaGGTA", circular=True).get_cutsites(Acc65I) == [((4, -4), Acc65I)]
    assert Dseq("CaGGTAC", circular=True).get_cutsites(Acc65I) == [((3, -4), Acc65I)]
    assert Dseq("aGGTACC", circular=True).get_cutsites(Acc65I) == [((2, -4), Acc65I)]
    assert Dseq("GGTACCa", circular=True).get_cutsites(Acc65I) == [((1, -4), Acc65I)]

    assert Dseq("GGTACCa", circular=True).get_cutsites(KpnI) == [((5, 4), KpnI)]
    assert Dseq("GTACCaG", circular=True).get_cutsites(KpnI) == [((4, 4), KpnI)]
    assert Dseq("TACCaGG", circular=True).get_cutsites(KpnI) == [((3, 4), KpnI)]
    assert Dseq("ACCaGGT", circular=True).get_cutsites(KpnI) == [((2, 4), KpnI)]
    assert Dseq("CCaGGTA", circular=True).get_cutsites(KpnI) == [((1, 4), KpnI)]
    assert Dseq("CaGGTAC", circular=True).get_cutsites(KpnI) == [((0, 4), KpnI)]
    assert Dseq("aGGTACC", circular=True).get_cutsites(KpnI) == [((6, 4), KpnI)]
    assert Dseq("GGTACCa", circular=True).get_cutsites(KpnI) == [((5, 4), KpnI)]

    assert Dseq("GGTACCa", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"),)
    assert Dseq("GTACCaG", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"),)
    assert Dseq("TACCaGG", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"),)
    assert Dseq("ACCaGGT", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"),)
    assert Dseq("CCaGGTA", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"),)
    assert Dseq("CaGGTAC", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"),)
    assert Dseq("aGGTACC", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"),)

    assert Dseq("GGTACCa", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"),)
    assert Dseq("GTACCaG", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"),)
    assert Dseq("TACCaGG", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"),)
    assert Dseq("ACCaGGT", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"),)
    assert Dseq("CCaGGTA", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"),)
    assert Dseq("CaGGTAC", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"),)
    assert Dseq("aGGTACC", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"),)

    assert Dseq("GGTACCaGGTACCa", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"), Dseq("PXEICaGQZFJ"))
    assert Dseq("GTACCaGGTACCaG", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"), Dseq("PXEICaGQZFJ"))
    assert Dseq("TACCaGGTACCaGG", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"), Dseq("PXEICaGQZFJ"))
    assert Dseq("ACCaGGTACCaGGT", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"), Dseq("PXEICaGQZFJ"))
    assert Dseq("CCaGGTACCaGGTA", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"), Dseq("PXEICaGQZFJ"))
    assert Dseq("CaGGTACCaGGTAC", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"), Dseq("PXEICaGQZFJ"))
    assert Dseq("aGGTACCaGGTACC", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"), Dseq("PXEICaGQZFJ"))
    assert Dseq("GGTACCaGGTACCa", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"), Dseq("PXEICaGQZFJ"))
    assert Dseq("GTACCaGGTACCaG", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"), Dseq("PXEICaGQZFJ"))
    assert Dseq("TACCaGGTACCaGG", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"), Dseq("PXEICaGQZFJ"))
    assert Dseq("ACCaGGTACCaGGT", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"), Dseq("PXEICaGQZFJ"))
    assert Dseq("CCaGGTACCaGGTA", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"), Dseq("PXEICaGQZFJ"))
    assert Dseq("CaGGTACCaGGTAC", circular=True).cut(Acc65I) == (Dseq("PXEICaGQZFJ"), Dseq("PXEICaGQZFJ"))

    assert Dseq("GGTACCaGGTACCa", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"), Dseq("QZFJCaGPXEI"))
    assert Dseq("GTACCaGGTACCaG", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"), Dseq("QZFJCaGPXEI"))
    assert Dseq("TACCaGGTACCaGG", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"), Dseq("QZFJCaGPXEI"))
    assert Dseq("ACCaGGTACCaGGT", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"), Dseq("QZFJCaGPXEI"))
    assert Dseq("CCaGGTACCaGGTA", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"), Dseq("QZFJCaGPXEI"))
    assert Dseq("CaGGTACCaGGTAC", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"), Dseq("QZFJCaGPXEI"))
    assert Dseq("aGGTACCaGGTACC", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"), Dseq("QZFJCaGPXEI"))
    assert Dseq("GGTACCaGGTACCa", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"), Dseq("QZFJCaGPXEI"))
    assert Dseq("GTACCaGGTACCaG", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"), Dseq("QZFJCaGPXEI"))
    assert Dseq("TACCaGGTACCaGG", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"), Dseq("QZFJCaGPXEI"))
    assert Dseq("ACCaGGTACCaGGT", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"), Dseq("QZFJCaGPXEI"))
    assert Dseq("CCaGGTACCaGGTA", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"), Dseq("QZFJCaGPXEI"))
    assert Dseq("CaGGTACCaGGTAC", circular=True).cut(KpnI) == (Dseq("QZFJCaGPXEI"), Dseq("QZFJCaGPXEI"))

    dig01 = Dseq("AGGTACCggtctcaAAA", circular=True).cut(BsaI, Acc65I)
    dig02 = Dseq("GGTACCggtctcaAAAA", circular=True).cut(BsaI, Acc65I)
    dig03 = Dseq("GTACCggtctcaAAAAG", circular=True).cut(BsaI, Acc65I)
    dig04 = Dseq("TACCggtctcaAAAAGG", circular=True).cut(BsaI, Acc65I)
    dig05 = Dseq("ACCggtctcaAAAAGGT", circular=True).cut(BsaI, Acc65I)
    dig06 = Dseq("CCggtctcaAAAAGGTA", circular=True).cut(BsaI, Acc65I)
    dig07 = Dseq("CggtctcaAAAAGGTAC", circular=True).cut(BsaI, Acc65I)
    dig08 = Dseq("ggtctcaAAAAGGTACC", circular=True).cut(BsaI, Acc65I)
    dig09 = Dseq("gtctcaAAAAGGTACCg", circular=True).cut(BsaI, Acc65I)
    dig10 = Dseq("tctcaAAAAGGTACCgg", circular=True).cut(BsaI, Acc65I)
    dig11 = Dseq("ctcaAAAAGGTACCggt", circular=True).cut(BsaI, Acc65I)
    dig12 = Dseq("tcaAAAAGGTACCggtc", circular=True).cut(BsaI, Acc65I)
    dig13 = Dseq("caAAAAGGTACCggtct", circular=True).cut(BsaI, Acc65I)
    dig14 = Dseq("AAAAGGTACCggtctca", circular=True).cut(BsaI, Acc65I)
    dig15 = Dseq("AAAGGTACCggtctcaA", circular=True).cut(BsaI, Acc65I)
    dig16 = Dseq("AAGGTACCggtctcaAA", circular=True).cut(BsaI, Acc65I)

    a = (
        dig01,
        dig02,
        dig03,
        dig04,
        dig05,
        dig06,
        dig07,
        dig08,
        dig09,
        dig10,
        dig11,
        dig12,
        dig13,
        dig14,
        dig15,
        dig16,
    )

    # x, y = set(a)

    # assert sorted(x) == sorted(y)


def test_BsaI():
    from pydna.dseq import Dseq
    from Bio.Restriction import BsaI

    assert Dseq("ggtctcAAgcTT").get_cutsites(BsaI) == [((7, -4), BsaI)]
    assert Dseq("ggtctcAAgcTT").cut(BsaI) == (Dseq("ggtctcAFqjZ"), Dseq("EpiXT"))
    assert Dseq("TggtctcAAgcT").get_cutsites(BsaI) == []
    assert Dseq("TggtctcAAgcT").cut(BsaI) == ()  # (Dseq("TggtctcAFqjZ"), Dseq("EpiX"))
    assert Dseq("TTggtctcAAgc").get_cutsites(BsaI) == []
    assert Dseq("TTggtctcAAgc").cut(BsaI) == ()

    assert Dseq("gtctcAAgcTTg", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)
    assert Dseq("ggtctcAAgcTT", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)
    assert Dseq("TggtctcAAgcT", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)
    assert Dseq("TTggtctcAAgc", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)
    assert Dseq("cTTggtctcAAg", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)



def test__add_fail():

    from pydna.dseq import Dseq

    assert Dseq("G") + Dseq("A") == Dseq("GA")

    with pytest.raises(TypeError):
        Dseq("GP") + Dseq("A")

    with pytest.raises(TypeError):
        Dseq("GE") + Dseq("A")

    with pytest.raises(TypeError):
        Dseq("GX") + Dseq("A")

    with pytest.raises(TypeError):
        Dseq("GI") + Dseq("A")

    with pytest.raises(TypeError):
        Dseq("GQ") + Dseq("A")

    with pytest.raises(TypeError):
        Dseq("GF") + Dseq("A")

    with pytest.raises(TypeError):
        Dseq("GZ") + Dseq("A")

    with pytest.raises(TypeError):
        Dseq("GJ") + Dseq("A")

    with pytest.raises(TypeError):
        Dseq("A") + Dseq("PG")

    with pytest.raises(TypeError):
        Dseq("A") + Dseq("EG")

    with pytest.raises(TypeError):
        Dseq("A") + Dseq("XG")

    with pytest.raises(TypeError):
        Dseq("A") + Dseq("IG")

    with pytest.raises(TypeError):
        Dseq("A") + Dseq("QG")

    with pytest.raises(TypeError):
        Dseq("A") + Dseq("FG")

    with pytest.raises(TypeError):
        Dseq("A") + Dseq("ZG")

    with pytest.raises(TypeError):
        Dseq("A") + Dseq("JG")


@pytest.mark.xfail(reason="Bug in Dseq.join")
def tests_misc():
    from pydna.dseq import Dseq
    from Bio.Restriction import Acc65I, NlaIV, KpnI, NotI

    assert Dseq("aa").cut(NotI) == ()

    obj = Dseq("GGTACC")

    obj_Acc65I = obj.cut(Acc65I)
    obj_KpnI = obj.cut(KpnI)
    obj_NlaIV = obj.cut(NlaIV)

    assert obj_Acc65I == (Dseq("GQZFJ"), Dseq("PXEIC"))
    assert obj_KpnI == (Dseq("GPXEI"), Dseq("QZFJC"))
    assert obj_NlaIV == (Dseq("GGT"), Dseq("ACC"))

    Dseq("").join(obj_Acc65I[::-1])
    cobj = obj.looped()

    for i in range(len(cobj) + 1):
        shifted_obj = cobj.shifted(i)

        assert shifted_obj.cut(Acc65I) == (Dseq("").join(obj_Acc65I[::-1]),)
        assert shifted_obj.cut(KpnI) == (Dseq("").join(obj_KpnI[::-1]),)
        assert shifted_obj.cut(NlaIV) == (Dseq("").join(obj_NlaIV[::-1]),)

    from Bio.Restriction import RestrictionBatch
    from Bio.Restriction import BamHI, BglII, BsaI, KpnI, Acc65I, XhoII, DpnI, BcgI, NlaIV
    from Bio.Restriction.Restriction_Dictionary import rest_dict, typedict

    enzymes = RestrictionBatch((BamHI, BglII, BsaI, KpnI, Acc65I, XhoII))

    obj = Dseq("GGTACCnnnGGtaCC")

    obj_Acc65I = obj.cut(Acc65I)
    obj_KpnI = obj.cut(KpnI)
    obj_NlaIV = obj.cut(NlaIV)

    assert obj_Acc65I == (Dseq("GQZFJ"), Dseq("PXEICnnnGQzfJ"), Dseq("PxeIC"))
    assert obj_KpnI == (Dseq("GPXEI"), Dseq("QZFJCnnnGPxeI"), Dseq("QzfJC"))
    assert obj_NlaIV == (Dseq("GGT"), Dseq("ACCnnnGGt"), Dseq("aCC"))


def test_doublecut_on_circular():
    pass
    # from Bio.Restriction import MaeII, MboI

    # enzymes = [MaeII, MboI]

    # aseq = Dseq("gatcacgt", circular=True)
    # bseq = Dseq("atcacgtg", circular=True)
    # cseq = Dseq("tcacgtga", circular=True)
    # dseq = Dseq("cacgtgat", circular=True)
    # eseq = Dseq("acgtgatc", circular=True)
    # fseq = Dseq("cgtgatca", circular=True)
    # gseq = Dseq("gtgatcac", circular=True)
    # hseq = Dseq("gatcacgt", circular=True)

    # a = aseq.get_cutsites(*enzymes)
    # b = bseq.get_cutsites(*enzymes)
    # c = cseq.get_cutsites(*enzymes)
    # d = dseq.get_cutsites(*enzymes)
    # e = eseq.get_cutsites(*enzymes)
    # f = fseq.get_cutsites(*enzymes)
    # g = gseq.get_cutsites(*enzymes)
    # h = hseq.get_cutsites(*enzymes)

    # print(a)
    # print(b)
    # print(c)
    # print(d)
    # print(e)
    # print(f)
    # print(g)
    # print(h)

    # print()

    # aa, *_ = _get_cutsite_pairs(a, True, 8)
    # bb, *_ = _get_cutsite_pairs(b, True, 8)
    # cc, *_ = _get_cutsite_pairs(c, True, 8)
    # dd, *_ = _get_cutsite_pairs(d, True, 8)
    # ee, *_ = _get_cutsite_pairs(e, True, 8)
    # ff, *_ = _get_cutsite_pairs(f, True, 8)
    # gg, *_ = _get_cutsite_pairs(g, True, 8)
    # hh, *_ = _get_cutsite_pairs(h, True, 8)

    # print(aa)
    # print(bb)
    # print(cc)
    # print(dd)
    # print(ee)
    # print(ff)
    # print(gg)
    # print(hh)

    # tuple(aseq.apply_cut(*cutsite_pair, 0, 4) for cutsite_pair in aa)
    # tuple(bseq.apply_cut(*cutsite_pair, 4, 2) for cutsite_pair in bb)
    # tuple(cseq.apply_cut(*cutsite_pair, 3, 2) for cutsite_pair in cc)
    # tuple(dseq.apply_cut(*cutsite_pair, 2, 2) for cutsite_pair in dd)
    # tuple(eseq.apply_cut(*cutsite_pair, 1, 2) for cutsite_pair in ee)
    # tuple(fseq.apply_cut(*cutsite_pair, 0, 2) for cutsite_pair in ff)
    # tuple(gseq.apply_cut(*cutsite_pair, 2, 4) for cutsite_pair in gg)
    # tuple(hseq.apply_cut(*cutsite_pair, 0, 4) for cutsite_pair in hh)
