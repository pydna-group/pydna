#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_dsiupac():
    from pydna.dseq import Dseq

    assert str(Dseq("PEXIaaaQFZJ")) == "PEXIaaaQFZJ"
    assert str(Dseq("QFZJaaaPEXI")) == "QFZJaaaPEXI"


def test_cut1():
    from pydna.dseq import Dseq

    """
    Acc65I.search(_Seq("GGTACC"))
    Out  [11]: [2]


        012345
        GGTACC
        CCATGG


    KpnI.search(_Seq("GGTACC"))
    Out  [12]: [6]
    """

    from Bio.Restriction import Acc65I, Bsp120I, KpnI, ApaI, TaiI, MaeII

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
    from pydna.dseq import Dseq

    s = Dseq("gattcatgcatgtagcttacgtagtct")

    RNA = "catgcatgtagcttacgtag"

    assert slice(0, 21, 1), slice(21, 27, 1) == s.cas9(RNA)


def test_initialization():
    import pytest
    from pydna.dseq import Dseq

    obj = Dseq("a")
    assert obj.crick == "t"
    assert obj.watson == "a"
    assert obj * 2 == Dseq("aa")
    assert not obj == 123
    assert obj * 0 == Dseq("")

    with pytest.raises(TypeError):
        obj * 2.3

    assert obj.seguid() == "ldseguid=ydezQsYTZgUCcb3-adxMaq_Xf8g"

    assert obj == Dseq("a", circular=False)

    # with pytest.raises(ValueError):
    #     DseXXXq("a", ovhg=0)

    # with pytest.raises(ValueError):
    #     DseXXXq("ttt", "tt")

    # with pytest.raises(ValueError):
    #     DseXXXq("ttt", "aa")

    ln = Dseq("gt")
    ci = ln.looped()

    assert not ln.circular
    assert ci.circular

    assert Dseq("gt", circular=False) == ln
    assert Dseq("gt", circular=True) == ci

    # assert Dseq.from_string("A") == Dseq("A")
    # assert Dseq.from_string("A", circular=True) == Dseq("A", circular=True)

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

    # assert Dseq.from_dsiupac("PEXIaaaQFZJ") == Dseq("GATCaaa", "CTAGttt", -4)
    # assert Dseq.from_dsiupac("QFZJaaaPEXI") == Dseq("aaaGATC", "tttCTAG", 4)


def test_cut_around_and_religate():
    from pydna.dseq import Dseq
    from pydna.utils import eq
    from Bio.Restriction import KpnI, BamHI, Acc65I

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

    for s in seqs:
        sek, enz, lin = s
        for i in range(len(sek)):
            zek = sek[i:] + sek[:i]
            cut_and_religate_Dseq(zek, enz, lin)


def test_Dseq_cutting_adding():
    from pydna.dseq import Dseq
    from Bio.Restriction import BamHI, PstI, EcoRI

    a = Dseq("GGATCCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtGGATCC")

    b = a.cut(BamHI)[1]

    assert b.watson == "GATCCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtG"
    assert b.crick == "GATCCacgcgcatatatatataattatagagagtaccgatgatgagcagcagaatagatcagtacgctacgatgatagtagatgaG"
    c = Dseq("nCTGCAGtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtGAATTCn")

    f, d, l_ = c.cut((EcoRI, PstI))

    assert d.watson == "GtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtG"
    assert d.crick == "AATTCacgcgcatatatatataattatagagagtaccgatgatgagcagcagaatagatcagtacgctacgatgatagtagatgaCTGCA"

    e = Dseq("nGAATTCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtCTGCAGn")

    f = e.cut((EcoRI, PstI))[1]

    assert f.watson == "AATTCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtCTGCA"
    assert f.crick == "GacgcgcatatatatataattatagagagtaccgatgatgagcagcagaatagatcagtacgctacgatgatagtagatgaG"


def test_dseq():
    import pytest
    import textwrap
    from pydna.dseq import Dseq

    obj1 = Dseq("a", circular=True)
    obj2 = Dseq("a")

    with pytest.raises(TypeError):
        obj1 + obj2

    with pytest.raises(TypeError):
        obj2 + obj1

    with pytest.raises(TypeError):
        obj1 + ""

    # with pytest.raises(AttributeError):
    #     obj2 + "" # TODO: discuss this

    obj1 = Dseq("ax")
    obj2 = Dseq("a")

    with pytest.raises(TypeError):
        obj1 + obj2

    obj = Dseq("gat", circular=True)

    assert obj[1:2] == Dseq("a")

    assert obj[:] == Dseq("gat", circular=False)

    obj = Dseq("atg", circular=False)

    assert obj[1:2]._data == b"atg"[1:2]

    assert obj[2:1]._data == b"atg"[2:1]

    assert obj.reverse_complement() == obj.rc() == Dseq("cat")

    obj = Dseq("atg", circular=True)

    assert obj.looped() == obj

    assert obj[:] == Dseq("atg", circular=False)

    assert obj[1:2]._data == b"atg"[1:2]

    assert obj[2:1]._data == b""  # TODO: change this?

    # obj = Dseq("G")
    # assert obj.five_prime_end() == ("5'", "g")
    # obj = Dseq("", "C", 0)
    # assert obj.five_prime_end() == ("3'", "c")

    # obj = Dseq("ccGGATCC", "aaggatcc", -2)
    obj = Dseq("iiGGATCCzz")
    assert obj._data == b"iiGGATCCzz"
    assert str(obj.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """\
    Dseq(-10)
    ccGGATCC
      CCTAGGaa
    """
    ).strip()
    assert repr(obj) == rpr

    assert obj[3] == Dseq("G")

    assert obj.fill_in("gt") == Dseq("ccGGATCCtt")

    assert obj + Dseq("") == obj  # TODO: fix this
    assert Dseq("") + obj == obj  # TODO: fix this

    obj = Dseq("pexiAAAAAAqfzj")
    assert obj.fill_in("gatc") == Dseq("gatcAAAAAAgatc")
    assert obj.fill_in("atc") == obj
    assert obj.fill_in("ac") == obj
    assert obj.fill_in("at") == obj

    obj = Dseq("qfzjAAAAAApexi")
    assert obj.fill_in("gatc") == obj
    assert obj.fill_in("atc") == obj
    assert obj.fill_in("ac") == obj
    assert obj.fill_in("at") == obj

    obj = Dseq("pexiAqfzj")
    assert obj.T4(b"GATC") == Dseq("gatcAgatc")
    assert obj.T4("at") == Dseq("pexiAqfzj")
    assert obj.T4("atg") == Dseq("patcAgatj")
    assert obj.T4("atgc") == Dseq("gatcAgatc")
    assert obj.mung() == Dseq("A")

    obj = Dseq("qfzjApexi")
    assert obj.T4("a") == Dseq("")
    assert obj.T4("t") == Dseq("")
    assert obj.T4("at") == Dseq("")
    assert obj.T4("atc") == Dseq("")
    assert obj.T4("atg") ==  Dseq("")
    assert obj.T4("atcg") == Dseq("A")

    assert Dseq("GGATCC").t4("gatc") == Dseq("GGATCC")
    assert Dseq("GGATCCe").t4("gatc") == Dseq("GGATCC")
    assert Dseq("eGGATCC").t4("gatc") == Dseq("aGGATCC")
    assert Dseq("eGGATCCe").t4("gatc") == Dseq("aGGATCC")
    assert Dseq("GGATCCz").t4("gatc") == Dseq("GGATCCt")
    assert Dseq("zGGATCC").t4("gatc") == Dseq("GGATCC")
    assert Dseq("zGGATCCz").t4("gatc") == Dseq("GGATCCt")

    assert Dseq("GGATII").t4("g") == Dseq("") # Dseq("gg", "", ovhg=0)
    assert Dseq("GGATCC").t4("gat") == Dseq('PPATJJ')

    a2 = Dseq("iiGGATCCee")
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

    a3 = Dseq("iiGGATCC")
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

    b = Dseq("qqGGATCCzz")
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

    b2 = Dseq("qqGGATCCee")
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

    b3 = Dseq("qqGGATCC")
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

    c = Dseq("GGATCCeee")
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

    d = Dseq("GGATCCzzz")
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

    obj = Dseq("GGATCCeee")

    from Bio.Restriction import BamHI

    frag1 = Dseq("GQFZJ")
    frag2 = Dseq("PEXICeee")

    assert obj.cut(BamHI) == (frag1, frag2)

    assert frag1 + frag2 == obj

    assert obj.seguid() == "ldseguid=qvssQpZe_4SlasGZYdKJSkuvQtc"

    assert frag1.seguid() == "ldseguid=jcVhCJ9Aa8aIQdBlkSU_XHTWmDc"
    assert frag2.seguid() == "ldseguid=SO1HxaZPDpcj-QffzS-mfF6_eag"

    assert frag1.rc().seguid() == "ldseguid=jcVhCJ9Aa8aIQdBlkSU_XHTWmDc"
    assert frag2.rc().seguid() == "ldseguid=SO1HxaZPDpcj-QffzS-mfF6_eag"

    obj = Dseq("tagcgtagctgtagtatgtgatctggtcta")

    assert repr(obj) == "Dseq(-30)\ntagcgtagctgtagtatgtgatctggtcta\natcgcatcgacatcatacactagaccagat"

    obj2 = Dseq("tagcgtagctgtagtatgtgatctggtcta")

    obj3 = obj = Dseq("tagcgtagctgtagtatgtgatctggtcta")

    assert obj == obj2 == obj3

    assert obj.find("ggatcc") == -1

    assert obj.find("tgtagta") == 9

    obj = Dseq("tagcgtagctgtagtatgtgatctggtctaa")

    obj = Dseq("tagcgtagctgtagtatgtgatctggtctaaQQQ")

    assert repr(obj) == "Dseq(-34)\ntagc..ctaa\natcg..gattCCC"

    obj = Dseq("tagcgtagctgtagtatgtgatctggtctaaIII")

    assert repr(obj) == "Dseq(-34)\ntagc..ctaaCCC\natcg..gatt"

    obj = Dseq("zagcgtagctgtagtatgtgatctggtctaa")
    assert repr(obj) == "Dseq(-31)\n agcg..ctaa\natcgc..gatt"

    obj = Dseq("Etagcgtagctgtagtatgtgatctggtctaa")
    assert repr(obj) == "Dseq(-32)\nAtagc..ctaa\n atcg..gatt"

    obj = Dseq("ftagcgtagctgtagtatgtgatctggtctaa")

    assert repr(obj) == "Dseq(-32)\n tagc..ctaa\ntatcg..gatt"

    assert round(obj.mw(), 1) == 19535.6

    obj1 = Dseq(
        "tagcgtagctgtagtatgtgatctggtcta",
        circular=True,
    )
    obj2 = Dseq(
        "tagcgtagctgtagtatgtgatctggtcta",
        circular=True,
    )
    obj3 = Dseq(
        "tagcgtagctgtagtatgtgatctggtcta",
        circular=True,
    )

    assert obj1 == obj2 == obj3

    assert obj1.find("ggatcc") == -1

    assert obj1.find("tgtagta") == 9

    assert Dseq("tagcgtagctgtagtatgtgatctggtcta").looped() == obj1

    from Bio.Restriction import BglII, BamHI

    obj = Dseq("ggatcc")

    assert BglII in obj.no_cutters()
    assert BamHI not in obj.no_cutters()

    assert BamHI in obj.unique_cutters()

    assert BamHI in obj.once_cutters()

    assert BamHI in (obj + obj).twice_cutters()
    assert BamHI not in obj.twice_cutters()

    assert BamHI in obj.n_cutters(1)
    assert BamHI in obj.cutters()

    from Bio.Restriction import RestrictionBatch

    rb = RestrictionBatch((BamHI, BglII))

    assert obj.cut(rb) == obj.cut(BamHI, BglII) == obj.cut(BglII, BamHI)

    obj = Dseq("ggatccAGATCT")

    assert obj.cut(rb) == obj.cut(BamHI, BglII) == obj.cut(BglII, BamHI)

    obj = Dseq("AGATCTggatcc")

    assert obj.cut(rb) == obj.cut(BamHI, BglII) == obj.cut(BglII, BamHI)

    obj = Dseq("ggatccAGATCT", circular=True)
    # TODO: address this test change Related to https://github.com/BjornFJohansson/pydna/issues/78
    assert obj.cut(rb) == obj.cut(BamHI, BglII) == obj.cut(BglII, BamHI)

    obj = Dseq("AGATCTggatcc", circular=True)

    assert obj.cut(rb) == obj.cut(BglII, BamHI) == obj.cut(BamHI, BglII)


def test_Dseq_slicing():
    from pydna.dseq import Dseq

    # from pydna.readers import read
    # from pydna.utils import eq

    # from Bio.Seq import Seq
    # from Bio.SeqRecord import SeqRecord as Srec
    from Bio.Restriction import BamHI

    a = Dseq("ggatcc")

    # assert a[:].watson == a.watson
    # assert a[:].crick == a.crick
    # assert a.ovhg == a[:].ovhg
    # b, c = a.cut(BamHI)
    # d = b[1:5]   # TODO: fix adding two ssdna together
    # e = d.rc()
    # assert e + d == Dseq("gatc")


def test_Dseq_slicing2():
    from pydna.dseq import Dseq
    from Bio.Restriction import BamHI, EcoRI, KpnI

    a = Dseq("aaGGATCCnnnnnnnnnGAATTCccc", circular=True)
    # TODO: address this test change Related to https://github.com/BjornFJohansson/pydna/issues/78

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
    from pydna.dseq import Dseq

    s = Dseq("GGATCC", circular=False)
    assert s[1:-1] == Dseq("GATC", circular=False)
    t = Dseq("GGATCC", circular=True)
    assert t[1:5] == Dseq("GATC")
    assert t[1:5].__dict__ == Dseq("GATC").__dict__
    assert s[1:5] == Dseq("GATC")
    assert s[1:5] == Dseq("GATC", circular=False)
    assert s[5:1:-1] == Dseq("CCTA")

    assert t[5:1] == Dseq("") # TODO: discuss this
    assert s[9:1] == Dseq("")
    assert t[9:1] == Dseq("")

    # Indexing of full circular molecule (https://github.com/BjornFJohansson/pydna/issues/161)
    s = Dseq("GGATCC", circular=True)
    str_seq = str(s)
    for shift in range(len(s)):
        assert str(s[shift:shift]) == str_seq[shift:] + str_seq[:shift]


def test_cut_circular():
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


def test_shifted():
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


def test_looped():

    # Looping a circular sequence should return a copy of the sequence
    # not the same sequence

    from pydna.dseq import Dseq

    a = Dseq("gatc", circular=True)

    assert a.looped() == a
    assert a.looped() is not a


def test_misc():
    from pydna.dseq import Dseq

    x = Dseq("ctcgGCGGCCGCcagcggccg", circular=True)

    from Bio.Restriction import NotI

    a, b = x.cut(NotI)

    z = (a + b).looped()
    # TODO: address this test change Related to https://github.com/BjornFJohansson/pydna/issues/78
    assert z.shifted(-6) == x


def test_cut_missing_enzyme():
    from pydna.dseq import Dseq

    x = Dseq("ctcgGCGGCCGCcagcggccg")

    from Bio.Restriction import PstI

    assert x.cut(PstI) == ()

    x = Dseq("ctcgGCGGCCGCcagcggccg", circular=True)

    assert x.cut(PstI) == ()


def test_cut_with_no_enzymes():
    from pydna.dseq import Dseq

    x = Dseq("ctcgGCGGCCGCcagcggccg")

    assert x.cut([]) == ()

    x = Dseq("ctcgGCGGCCGCcagcggccg", circular=True)

    assert x.cut([]) == ()


def test_transcribe():
    from pydna.dseq import Dseq

    x = Dseq("ATGAAATAA")

    assert str(x.transcribe()) == "AUGAAAUAA"

    assert str(x.reverse_complement().transcribe()) == "UUAUUUCAU"


def test_translate():
    from pydna.dseq import Dseq

    x = Dseq("ATGAAATAA")

    assert str(x.translate()) == "MK*"

    assert str(x.reverse_complement().translate()) == "LFH"


def test_from_full_sequence_and_overhangs():
    from pydna.dseq import Dseq

    test_cases = [
        (2, 2, "AAAA", "TTTT", "FFAAEE"),
        (-2, 2, "AAAAAA", "TT", "EEAAEE"),
        (2, -2, "AA", "TTTTTT", "FFAAFF"),
        (-2, -2, "AAAA", "TTTT", "EEAAFF"),
        (0, 0, "AAAAAA", "TTTTTT", "AAAAAA"),
    ]
    for crick_ovhg, watson_ovhg, watson, crick, dsiupac in test_cases:
        dseq_1 = Dseq.from_full_sequence_and_overhangs("AAAAAA", crick_ovhg=crick_ovhg, watson_ovhg=watson_ovhg)
        dseq_2 = Dseq(dsiupac, circular=False)

        assert dseq_1 == dseq_2
        assert dseq_2.watson_ovhg() == watson_ovhg


def test_right_end_position():

    from pydna.dseq import Dseq

    test_cases = [
        ("AAA", "TT", "AAE", (3, 2)),
        ("AA", "TTT", "AAF", (2, 3)),
        ("AAA", "TTT", "AAA", (3, 3)),
    ]
    for watson, crick, dsiupac, expected in test_cases:
        dseq = Dseq(dsiupac, circular=False)
        assert dseq.right_end_position() == expected


def test_left_end_position():

    from pydna.dseq import Dseq

    test_cases = [
        ("AAA", "TT",  "EAA", (0, 1), -1),
        ("AA", "TTT", "FAA", (1, 0), 1),
        ("AAA", "TTT", "AAA", (0, 0), 0),
    ]
    for watson, crick, dsiupac, expected, ovhg in test_cases:
        dseq = Dseq(dsiupac, circular=False)
        assert dseq.left_end_position() == expected


# def test_apply_cut():
#     from pydna.dseq import Dseq
#     from Bio.Restriction import EcoRI, BamHI

#     seq = Dseq("aaGAATTCaa", circular=False)

#     # A cut where both sides are None returns the same sequence
#     assert seq.apply_cut(None, None) == seq

#     # A cut where one side is None leaves that side intact
#     EcoRI_cut = ((3, -4), None)

#     assert seq.apply_cut(None, EcoRI_cut) == Dseq.from_full_sequence_and_overhangs(
#         "aaGAATT", watson_ovhg=-4, crick_ovhg=0
#     )
#     assert seq.apply_cut(EcoRI_cut, None) == Dseq.from_full_sequence_and_overhangs(
#         "AATTCaa", watson_ovhg=0, crick_ovhg=-4
#     )

#     # It respects the original overhang
#     seq = Dseq.from_full_sequence_and_overhangs("aaGAATTCaa", watson_ovhg=1, crick_ovhg=1)
#     assert seq.apply_cut(None, EcoRI_cut) == Dseq.from_full_sequence_and_overhangs(
#         "aaGAATT", watson_ovhg=-4, crick_ovhg=1
#     )
#     assert seq.apply_cut(EcoRI_cut, None) == Dseq.from_full_sequence_and_overhangs(
#         "AATTCaa", watson_ovhg=1, crick_ovhg=-4
#     )

#     seq = Dseq.from_full_sequence_and_overhangs("aaGAATTCaa", watson_ovhg=-1, crick_ovhg=-1)
#     assert seq.apply_cut(None, EcoRI_cut) == Dseq.from_full_sequence_and_overhangs(
#         "aaGAATT", watson_ovhg=-4, crick_ovhg=-1
#     )
#     assert seq.apply_cut(EcoRI_cut, None) == Dseq.from_full_sequence_and_overhangs(
#         "AATTCaa", watson_ovhg=-1, crick_ovhg=-4
#     )

#     # A repeated cut in a circular molecule opens it up
#     seq = Dseq("aaGAATTCaa", circular=True)
#     assert seq.apply_cut(EcoRI_cut, EcoRI_cut) == Dseq.from_full_sequence_and_overhangs(
#         "AATTCaaaaGAATT", watson_ovhg=-4, crick_ovhg=-4
#     )

#     # Two cuts extract a subsequence
#     seq = Dseq("aaGAATTCaaGAATTCaa", circular=True)
#     EcoRI_cut_2 = ((11, -4), None)
#     assert seq.apply_cut(EcoRI_cut, EcoRI_cut_2) == Dseq.from_full_sequence_and_overhangs(
#         "AATTCaaGAATT", watson_ovhg=-4, crick_ovhg=-4
#     )

#     # Overlapping cuts should return an error
#     seq = Dseq("aaGAATTCaa", circular=True)
#     first_cuts = [
#         ((3, -4), BamHI),
#         ((7, 4), BamHI),
#         # Spanning the origin
#         ((9, -8), BamHI),
#         ((8, 8), BamHI),
#     ]

#     overlapping_cuts = [
#         ((4, -4), EcoRI),
#         ((2, -4), EcoRI),
#         ((2, -6), EcoRI),
#         ((8, 4), EcoRI),
#         ((6, 4), EcoRI),
#         ((8, 6), EcoRI),
#         # Spanning the origin
#         ((7, -8), EcoRI),
#         ((6, 8), EcoRI),
#     ]

#     for first_cut in first_cuts:
#         for second_cut in overlapping_cuts:
#             try:
#                 seq.apply_cut(first_cut, second_cut)
#             except ValueError as e:
#                 assert e.args[0] == "Cuts by BamHI EcoRI overlap."
#             else:
#                 print(first_cut, second_cut)
#                 assert False, "Expected ValueError"

#     # Rotating the sequence, apply the same cut
#     seq = Dseq("acgtATGaatt", circular=True)
#     for shift in range(len(seq)):
#         seq_shifted = seq.shifted(shift)
#         start = 4 - shift
#         if start < 0:
#             start += len(seq)
#         # Cut with negative ovhg
#         new_cut = ((start, -3), None)
#         out = seq_shifted.apply_cut(new_cut, new_cut)
#         assert str(out) == "ATGaattacgtATG"

#         # Cut with positive ovhg
#         start = (start + 3) % len(seq)
#         new_cut = ((start, 3), None)
#         out = seq_shifted.apply_cut(new_cut, new_cut)
#         assert str(out) == "ATGaattacgtATG"

#         # A blunt cut
#         start = 4 - shift
#         new_cut = ((start, 0), None)
#         out = seq_shifted.apply_cut(new_cut, new_cut)
#         assert str(out) == "ATGaattacgt"


# def test_cutsite_is_valid():

#     from pydna.dseq import Dseq
#     from Bio.Restriction import EcoRI, PacI, NmeDI, EcoRV

#     # Works for circular case
#     seqs = ["GAATTC", "TTAATTAAC", "GATATC"]
#     enzs = [EcoRI, PacI, EcoRV]
#     for seq, enz in zip(seqs, enzs):
#         dseq = Dseq(seq, circular=True)
#         for shift in range(len(seq)):
#             dseq_shifted = dseq.shifted(shift)
#             (cutsite,) = dseq_shifted.get_cutsites([enz])

#             assert dseq_shifted.cutsite_is_valid(cutsite)

#     # Works for overhangs
#     seqs = ["GAATTC", "TTAATTAA", "GATATC"]
#     for seq, enz in zip(seqs, enzs):
#         for ovhg in [-1, 0, 1]:
#             dseq = Dseq.from_full_sequence_and_overhangs(seq, ovhg, 0)
#             if ovhg != 0:
#                 assert len(dseq.get_cutsites([enz])) == 0
#             else:
#                 assert len(dseq.get_cutsites([enz])) == 1

#             dseq = Dseq.from_full_sequence_and_overhangs(seq, 0, ovhg)
#             if ovhg != 0:
#                 assert len(dseq.get_cutsites([enz])) == 0
#             else:
#                 assert len(dseq.get_cutsites([enz])) == 1

#     # Special cases:
#     dseq = Dseq.from_full_sequence_and_overhangs("AAAAAAAAAAAAAGCCGGCAAAAAAAAAAAA", 0, 0)
#     assert len(dseq.get_cutsites([NmeDI])) == 2
#     # Remove left cutting place
#     assert len(dseq[2:].get_cutsites([NmeDI])) == 1
#     # Remove right cutting place
#     assert len(dseq[:-2].get_cutsites([NmeDI])) == 1
#     # Remove both cutting places
#     assert len(dseq[2:-2].get_cutsites([NmeDI])) == 0

#     # overhang left side
#     dseq = Dseq.from_full_sequence_and_overhangs("AAAAAAAAAAAAAGCCGGCAAAAAAAAAAAA", -2, 0)
#     assert len(dseq.get_cutsites([NmeDI])) == 1
#     dseq = Dseq.from_full_sequence_and_overhangs("AAAAAAAAAAAAAGCCGGCAAAAAAAAAAAA", 2, 0)
#     assert len(dseq.get_cutsites([NmeDI])) == 1

#     # overhang right side
#     dseq = Dseq.from_full_sequence_and_overhangs("AAAAAAAAAAAAAGCCGGCAAAAAAAAAAAA", 0, 2)
#     assert len(dseq.get_cutsites([NmeDI])) == 1
#     dseq = Dseq.from_full_sequence_and_overhangs("AAAAAAAAAAAAAGCCGGCAAAAAAAAAAAA", 0, -2)
#     assert len(dseq.get_cutsites([NmeDI])) == 1

#     # overhang both sides
#     dseq = Dseq.from_full_sequence_and_overhangs("AAAAAAAAAAAAAGCCGGCAAAAAAAAAAAA", 2, 2)
#     assert len(dseq.get_cutsites([NmeDI])) == 0
#     dseq = Dseq.from_full_sequence_and_overhangs("AAAAAAAAAAAAAGCCGGCAAAAAAAAAAAA", -2, -2)
#     assert len(dseq.get_cutsites([NmeDI])) == 0

#     # overhang on recognition site removes both cutting places
#     dseq = Dseq.from_full_sequence_and_overhangs("AAAAAAAAAAAAAGCCGGCAAAAAAAAAAAA", 16, 0)
#     assert len(dseq.get_cutsites([NmeDI])) == 0
#     dseq = Dseq.from_full_sequence_and_overhangs("AAAAAAAAAAAAAGCCGGCAAAAAAAAAAAA", 0, 16)
#     assert len(dseq.get_cutsites([NmeDI])) == 0


# def test_get_cutsite_pairs():
#     from pydna.dseq import Dseq

#     # in the test, we replace cuts by integers for clarity.

#     dseq = Dseq("A")

#     # Empty returns empty list
#     assert dseq.get_cutsite_pairs([]) == []

#     # Single cut on linear seq returns two fragments
#     assert dseq.get_cutsite_pairs([1]) == [(None, 1), (1, None)]

#     # Two cuts on linear seq return three fragments
#     assert dseq.get_cutsite_pairs([1, 2]) == [(None, 1), (1, 2), (2, None)]

#     dseq = Dseq("A", circular=True)

#     # Empty returns empty list
#     assert dseq.get_cutsite_pairs([]) == []

#     # Single cut on circular seq returns opened molecule
#     assert dseq.get_cutsite_pairs([1]) == [(1, 1)]

#     # Two cuts on circular seq return 2 fragments
#     assert dseq.get_cutsite_pairs([1, 2]) == [(1, 2), (2, 1)]


# def test_get_cut_parameters():

#     from pydna.dseq import Dseq

#     dseq = Dseq.from_full_sequence_and_overhangs("aaaACGTaaa", 3, 3)
#     assert dseq.get_cut_parameters(None, True) == (*dseq.left_end_position(), dseq.ovhg)
#     assert dseq.get_cut_parameters(None, False) == (*dseq.right_end_position(), dseq.watson_ovhg())

#     assert dseq.get_cut_parameters(((4, -2), None), True) == (4, 6, -2)
#     assert dseq.get_cut_parameters(((4, -2), None), False) == (4, 6, -2)
#     assert dseq.get_cut_parameters(((6, 2), None), True) == (6, 4, 2)
#     assert dseq.get_cut_parameters(((6, 2), None), False) == (6, 4, 2)

#     dseq = Dseq("aaaACGTaaa", circular=True)

#     # None cannot be used on circular molecules
#     try:
#         assert dseq.get_cut_parameters(None, True) == (*dseq.left_end_position(), dseq.ovhg)
#     except AssertionError as e:
#         assert e.args[0] == "Circular sequences should not have None cuts"
#     else:
#         assert False, "Expected AssertionError"

#     try:
#         assert dseq.get_cut_parameters(None, False) == (*dseq.right_end_position(), dseq.watson_ovhg())
#     except AssertionError as e:
#         assert e.args[0] == "Circular sequences should not have None cuts"
#     else:
#         assert False, "Expected AssertionError"

#     # "Normal" cuts
#     assert dseq.get_cut_parameters(((4, -2), None), True) == (4, 6, -2)
#     assert dseq.get_cut_parameters(((4, -2), None), False) == (4, 6, -2)
#     assert dseq.get_cut_parameters(((6, 2), None), True) == (6, 4, 2)
#     assert dseq.get_cut_parameters(((6, 2), None), False) == (6, 4, 2)

#     # Origin-spannign cuts
#     assert dseq.get_cut_parameters(((9, -2), None), True) == (9, 1, -2)
#     assert dseq.get_cut_parameters(((9, -2), None), False) == (9, 1, -2)
#     assert dseq.get_cut_parameters(((1, 2), None), True) == (1, 9, 2)
#     assert dseq.get_cut_parameters(((1, 2), None), False) == (1, 9, 2)


def test_checksums():

    from seguid import ldseguid, cdseguid
    from pydna.dseq import Dseq

    # AT
    # TA

    dlDNA_ldseguid = "odgytmQKSOnFEUorGIWK3NDjqUA"
    truth = f"ldseguid={dlDNA_ldseguid}"
    assert ldseguid("AT", "AT") == truth == Dseq("AT").seguid()

    #  -AT
    #  AT-

    dlDNA2_ldseguid = "-9xkp3UfucL4bSPxYODh8i9KFEE"
    truth = f"ldseguid={dlDNA2_ldseguid}"
    assert ldseguid("-AT", "-TA") == truth == Dseq("ZAX").seguid()

    # TA-
    # -TA

    dlDNA3_ldseguid = "kI9qYVNRPF8epm2xem0ZUP8J-CI"
    truth = f"ldseguid={dlDNA3_ldseguid}"
    assert ldseguid("TA-", "AT-") == truth == Dseq("XAZ").seguid()

    # CTATAG
    # --TA--

    dlDNA4_ldseguid = "ToSxUXWMCIKz-FYdXJ3Qq-bS_8o"
    truth = f"ldseguid={dlDNA4_ldseguid}"
    assert ldseguid("CTATAG", "--AT--") == truth == Dseq("IXATEP").seguid()

    # --AT--
    # GATATC

    assert ldseguid("--AT--", "CTATAG") == truth == Dseq("JZATFQ").seguid()

    truth = "cdseguid=5fHMG19IbYxn7Yr7_sOCkvaaw7U"
    assert cdseguid("ACGTT", "AACGT") == truth == Dseq("ACGTT", circular=True).seguid()
    assert cdseguid("AACGT", "ACGTT") == truth == Dseq("AACGT", circular=True).seguid()


def test_new():
    from pydna.dseq import Dseq
    from Bio.Restriction import KpnI, Acc65I, BsaI, XmaI, SmaI, BamHI

    fiveoh = Dseq('PEXIaaaQFZJ')
    assert str(fiveoh + fiveoh) == 'PEXIaaaGATCaaaQFZJ'
    threeoh = Dseq('QFZJQtttPEXIP')
    assert str(threeoh + threeoh) == "QFZJQtttGATCGtttPEXIP"

    assert repr(Dseq("AIXEP") + Dseq("JZFQA")) == "Dseq(-6)\nACTAGA\nTGATCT"
    assert repr(Dseq("AIXEP") + Dseq("ZFQA")) == "Dseq(-6)\nACTAGA\nT ATCT"
    assert repr(Dseq("AIXE") + Dseq("JZFQA")) == "Dseq(-6)\nACTA A\nTGATCT"
    assert repr(Dseq("AP") + Dseq("QA")) == "Dseq(-3)\nAGA\nTCT"
    assert repr(Dseq("AE") + Dseq("FA")) == "Dseq(-3)\nAAA\nTTT"
    assert repr(Dseq("AX") + Dseq("ZA")) == "Dseq(-3)\nATA\nTAT"
    assert repr(Dseq("AI") + Dseq("JA")) == "Dseq(-3)\nACA\nTGT"
    assert repr(Dseq("APP") + Dseq("QA")) == "Dseq(-4)\nAGGA\nT CT"
    assert repr(Dseq("AEE") + Dseq("FA")) == "Dseq(-4)\nAAAA\nT TT"
    assert repr(Dseq("AXX") + Dseq("ZA")) == "Dseq(-4)\nATTA\nT AT"
    assert repr(Dseq("AII") + Dseq("JA")) == "Dseq(-4)\nACCA\nT GT"

    assert repr(Dseq("AP") + Dseq("QQA")) == "Dseq(-4)\nAG A\nTCCT"
    assert repr(Dseq("AE") + Dseq("FFA")) == "Dseq(-4)\nAA A\nTTTT"
    assert repr(Dseq("AX") + Dseq("ZZA")) == "Dseq(-4)\nAT A\nTAAT"
    assert repr(Dseq("AI") + Dseq("JJA")) == "Dseq(-4)\nAC A\nTGGT"

    assert Dseq("QAP").looped()[:] == Dseq("GA")
    assert Dseq("PAQ").looped()[:] == Dseq("GA")
    assert Dseq("EAF").looped()[:] == Dseq("AA")
    assert Dseq("FAE").looped()[:] == Dseq("AA")
    assert Dseq("XAZ").looped()[:] == Dseq("TA")
    assert Dseq("ZAX").looped()[:] == Dseq("TA")
    assert Dseq("IAJ").looped()[:] == Dseq("CA")
    assert Dseq("JAI").looped()[:] == Dseq("CA")

    s = Dseq("PEXIAAAQFZJ")

    assert s.looped()._data == b'GATCAAA'

    s = Dseq("PEXIAAAQFZJ")

    assert s._fill_in_three_prime("gatc")._data == b'GATCAAAQFZJ'
    assert s._fill_in_three_prime("gat")._data == b'PATCAAAQFZJ'
    assert s._fill_in_three_prime("ga")._data == b'PETCAAAQFZJ'
    assert s._fill_in_three_prime("g")._data == b'PEXCAAAQFZJ'
    assert s._fill_in_three_prime("")._data == s._data

    assert s._fill_in_three_prime("atc")._data == s._data
    assert s._fill_in_three_prime("at")._data == s._data
    assert s._fill_in_three_prime("a")._data == s._data
    assert s._fill_in_three_prime("")._data == s._data

    assert s._fill_in_five_prime("gatc")._data == b'PEXIAAAGATC'
    assert s._fill_in_five_prime("gat")._data == b'PEXIAAAGATJ'
    assert s._fill_in_five_prime("ga")._data == b'PEXIAAAGAZJ'
    assert s._fill_in_five_prime("g")._data == b'PEXIAAAGFZJ'
    assert s._fill_in_five_prime("")._data == s._data

    assert s._fill_in_five_prime("atc")._data == s._data
    assert s._fill_in_five_prime("at")._data == s._data
    assert s._fill_in_five_prime("a")._data == s._data
    assert s._fill_in_five_prime("")._data == s._data

    assert s.fill_in("gatc")._data == b'GATCAAAGATC'
    assert s.fill_in("gat")._data == b'PATCAAAGATJ'
    assert s.fill_in("ga")._data == b'PETCAAAGAZJ'
    assert s.fill_in("g")._data == b'PEXCAAAGFZJ'
    assert s.fill_in("")._data == s._data

    assert s == Dseq.from_representation(repr(s))

    assert Dseq.from_full_sequence_and_overhangs('AAAAAA', crick_ovhg=2, watson_ovhg=2)._data == Dseq("FFAAEE")._data
    assert Dseq.from_full_sequence_and_overhangs('AAAAAA', crick_ovhg=-2, watson_ovhg=2)._data == Dseq("EEAAEE")._data
    assert Dseq.from_full_sequence_and_overhangs('AAAAAA', crick_ovhg=2, watson_ovhg=-2)._data == Dseq("FFAAFF")._data
    assert Dseq.from_full_sequence_and_overhangs('AAAAAA', crick_ovhg=-2, watson_ovhg=-2)._data == Dseq("EEAAFF")._data

    assert repr(Dseq("XIpexiPEXIpexiAqfzjQFZJqfzjQF")) == (
        "Dseq(-29)\n" "TCgatcGATCgatcA\n" "              TctagCTAGctagCT"
    )

    assert repr(Dseq("XIpexiPEXIpexiAAqfzjQFZJqfzjQF")) == ( "Dseq(-30)\n"
                                                             "TCgatcGATCgatcAA\n"
                                                             "              TTctagCTAGctagCT"
    )

    assert repr(Dseq("EXIpexiPEXIpexiAqfzjQFZJqfzjQFZ")) == ("Dseq(-31)\n" "ATCg..gatcA\n" "          Tctag..gCTA")

    assert repr(Dseq("XIpexiPEXIpexiAaAqfzjQFZJqfzjQF")) == ("Dseq(-31)\n" "TCga..gatcAaA\n" "          TtTctag..agCT")

    dsdna = """
               Dseq(
                GGATCC
               aCCTAGGg"""

    assert Dseq.from_representation(dsdna)._data == b"zGGATCCj"

    dsdna = """Dseq(
               aGGATCCg
                CCTAGG"""

    assert Dseq.from_representation(dsdna)._data == b"eGGATCCp"

    dsdna = """Dseq(
               GGATCCaa
             cccctagg"""

    assert Dseq.from_representation(dsdna)._data == b"qqGGATCCee"

    dsdna = """Dseq(
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

    assert Dseq('PEXIGGATCCQFZJ').T4(b"CG") == Dseq(b'PEXCGGATCCGFZJ')
    assert Dseq('PAGAJ').T4(b"ATG") == Dseq(b'PAGAJ')
    assert Dseq('QFZJGGATCCPEXI').T4("") == Dseq("")

    assert Dseq("pexiAqfzj").T4(b"gatc") == Dseq("gatcAgatc")

    assert Dseq("qqGGATCCee").seguid() == 'ldseguid=F0z-LxHZqAK3HvqQiqjM7A28daE'

    ##################################################################################

    assert Dseq("ggtctcAAgcTT", circular=False).get_cutsites(BsaI) == [(7, BsaI)]
    assert Dseq("TggtctcAAgcT", circular=False).get_cutsites(BsaI) == [(8, BsaI)]
    assert Dseq("TTggtctcAAgc", circular=False).get_cutsites(BsaI) == []
    assert Dseq("cTTggtctcAAg", circular=False).get_cutsites(BsaI) == []
    assert Dseq("gcTTggtctcAA", circular=False).get_cutsites(BsaI) == []

    assert Dseq("gtctcAAgcTTg", circular=True).get_cutsites(BsaI) == [(6, BsaI)]
    assert Dseq("ggtctcAAgcTT", circular=True).get_cutsites(BsaI) == [(7, BsaI)]
    assert Dseq("TggtctcAAgcT", circular=True).get_cutsites(BsaI) == [(8, BsaI)]
    assert Dseq("TTggtctcAAgc", circular=True).get_cutsites(BsaI) == [(9, BsaI)]
    assert Dseq("cTTggtctcAAg", circular=True).get_cutsites(BsaI) == [(10, BsaI)]
    assert Dseq("gtctcAAgcTTg", circular=True).get_cutsites(BsaI) == [(6, BsaI)]

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

    assert Dseq("GGTACCa", circular = True).get_cutsites(Acc65I) == [(1, Acc65I)]
    assert Dseq("GTACCaG", circular = True).get_cutsites(Acc65I) == [(0, Acc65I)]
    assert Dseq("TACCaGG", circular = True).get_cutsites(Acc65I) == [(6, Acc65I)]
    assert Dseq("ACCaGGT", circular = True).get_cutsites(Acc65I) == [(5, Acc65I)]
    assert Dseq("CCaGGTA", circular = True).get_cutsites(Acc65I) == [(4, Acc65I)]
    assert Dseq("CaGGTAC", circular = True).get_cutsites(Acc65I) == [(3, Acc65I)]
    assert Dseq("aGGTACC", circular = True).get_cutsites(Acc65I) == [(2, Acc65I)]

    assert Dseq("GGTACCa", circular = True).get_cutsites(KpnI) == [(5, KpnI)]
    assert Dseq("GTACCaG", circular = True).get_cutsites(KpnI) == [(4, KpnI)]
    assert Dseq("TACCaGG", circular = True).get_cutsites(KpnI) == [(3, KpnI)]
    assert Dseq("ACCaGGT", circular = True).get_cutsites(KpnI) == [(2, KpnI)]
    assert Dseq("CCaGGTA", circular = True).get_cutsites(KpnI) == [(1, KpnI)]
    assert Dseq("CaGGTAC", circular = True).get_cutsites(KpnI) == [(0, KpnI)]
    assert Dseq("aGGTACC", circular = True).get_cutsites(KpnI) == [(6, KpnI)]

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

    a = (dig01, dig02, dig03, dig04, dig05, dig06, dig07, dig08, dig09, dig10,
         dig11, dig12, dig13, dig14, dig15, dig16,)

    x, y = set(a)

    assert sorted(x) == sorted(y)

    assert Dseq("ggtctcAAgcTT").get_cutsites(BsaI) == [(7, BsaI)]
    assert Dseq("ggtctcAAgcTT").cut(BsaI) == (Dseq("ggtctcAFqjZ"), Dseq("EpiXT"))
    assert Dseq("TggtctcAAgcT").get_cutsites(BsaI) == [(8, BsaI)]
    assert Dseq("TggtctcAAgcT").cut(BsaI) == (Dseq("TggtctcAFqjZ"), Dseq("EpiX"))
    assert Dseq("TTggtctcAAgc").get_cutsites(BsaI) == []
    assert Dseq("TTggtctcAAgc").cut(BsaI) == ()

    assert Dseq("gtctcAAgcTTg", circular=True).get_cutsites(BsaI) == [(6, BsaI)]
    assert Dseq("ggtctcAAgcTT", circular=True).get_cutsites(BsaI) == [(7, BsaI)]
    assert Dseq("TggtctcAAgcT", circular=True).get_cutsites(BsaI) == [(8, BsaI)]
    assert Dseq("TTggtctcAAgc", circular=True).get_cutsites(BsaI) == [(9, BsaI)]
    assert Dseq("cTTggtctcAAg", circular=True).get_cutsites(BsaI) == [(10, BsaI)]
    assert Dseq("gtctcAAgcTTg", circular=True).get_cutsites(BsaI) == [(6, BsaI)]

    assert Dseq("gtctcAAgcTTg", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)
    assert Dseq("ggtctcAAgcTT", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)
    assert Dseq("TggtctcAAgcT", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)
    assert Dseq("TTggtctcAAgc", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)
    assert Dseq("cTTggtctcAAg", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
