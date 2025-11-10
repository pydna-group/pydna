#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import textwrap

from pydna.dseq import Dseq
from pydna.utils import eq

from Bio.Restriction import (
    Acc65I, ApaI, BamHI, BglII, BsaI, Bsp120I, EcoRI, EcoRV,
    KpnI, MaeII, NmeDI, NotI, PacI, PstI, RestrictionBatch, TaiI
)

from seguid import ldseguid
from seguid import cdseguid

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

    assert round(obj.mw(), 1) == 19535.6

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
        ("AAA", "TT", (3, 2)),
        ("AA", "TTT", (2, 3)),
        ("AAA", "TTT", (3, 3)),
    ]
    for watson, crick, expected in test_cases:
        dseq = Dseq(watson, crick, ovhg=0, circular=False)
        assert dseq.right_end_position() == expected


def test_left_end_position():



    test_cases = [
        ("AAA", "TT", (0, 1), -1),
        ("AA", "TTT", (1, 0), 1),
        ("AAA", "TTT", (0, 0), 0),
    ]
    for watson, crick, expected, ovhg in test_cases:
        dseq = Dseq(watson, crick, ovhg=ovhg, circular=False)
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
