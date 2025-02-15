#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_misc():

    from Bio.Restriction import BamHI, AjiI, ZraI, KpnI, Acc65I, BsaI

    from pydna.dseq import Dseq

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

    from Bio.Restriction import Acc65I, BamHI, EcoRI, HindIII, KpnI, PstI, SacI, SalI, SbfI, SmaI, SphI, XbaI, XmaI

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

    puc19mcs.cut(SmaI)
    puc19mcs.cut(BamHI)
    puc19mcs.cut(SmaI, BamHI)

    puc19mcs.cut(Acc65I, KpnI)

    x = Dseq("PEXIaUaOaQFZJ")

    assert x.seguid() == "ldseguid=AA-pG6wsGMu2J-AuWBWlOi3iDc4"

    y = Dseq("QFZJaUaOaPEXI")

    assert y.seguid() == "ldseguid=10e3S-WXX1z0BHqpbYTPeb2MwUY"

    # Dseq(-5)
    # CUGOT
    # GACOA


def test_pair():

    from pydna.dseq import Dseq, pair

    x = Dseq("PEXIaUaOaQFZJ")

    x_wco = Dseq(pair("GATCaUaAa", "tAtUtCTAG"[::-1], -4))

    x_wco2 = Dseq(pair(x.watson, x.crick, x.ovhg))

    assert x == x_wco == x_wco2

    y = Dseq("QFZJaUaOaPEXI")

    y_wco = Dseq(pair("aUaAaGATC", "CTAGtAtUt"[::-1], 4))

    y_wco2 = Dseq(pair(y.watson, y.crick, y.ovhg))

    assert y == y_wco == y_wco2

    assert pair("GATT", "CTAA"[::-1]) == b"GATT"

    assert pair("GATTa", "aCTAA"[::-1]) == b"zGATTe"

    assert pair("GATCaUaAa", "tAtUtCTAG"[::-1], -4) == b"PEXIaUaOaQFZJ"


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

    from Bio.Restriction import BamHI

    # x = Dseq("TCCtGAATTCcGGTCTCGGTACCaGGA", circular=True)

    # assert x.get_cutsites(BamHI)[0].recsite.extract(x) == Dseq("GGATCC")


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

    Dseq.from_representation("GTACCTTTG\n" "    GAAACCTAG")
    Dseq.from_representation("\nGTACCTTTG\n" "    GAAACCTAG")

    Dseq.from_representation(
        """
                             GTACCTTTG
                                 GAAACCTAG"""
    )

    Dseq.from_representation("\nGTACCTTTG\n    GAAACCTAG")

    with pytest.raises(KeyError):
        Dseq.from_representation(
            """GTACCTTTG
                                        GAAACCTAG"""
        )


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

    from pydna import _PydnaWarning

    with pytest.warns(_PydnaWarning):
        obj = Dseq("qfzjApexi")
        assert obj.T4("a") == Dseq("")
        assert obj.T4("t") == Dseq("")
        assert obj.T4("at") == Dseq("")
        assert obj.T4("atc") == Dseq("")
        assert obj.T4("atg") == Dseq("")
        assert obj.T4("atcg") == Dseq("A")

    assert Dseq("GGATCC").t4("gatc") == Dseq("GGATCC")
    assert Dseq("GGATCCe").t4("gatc") == Dseq("GGATCC")
    assert Dseq("eGGATCC").t4("gatc") == Dseq("aGGATCC")
    assert Dseq("eGGATCCe").t4("gatc") == Dseq("aGGATCC")
    assert Dseq("GGATCCz").t4("gatc") == Dseq("GGATCCt")
    assert Dseq("zGGATCC").t4("gatc") == Dseq("GGATCC")
    assert Dseq("zGGATCCz").t4("gatc") == Dseq("GGATCCt")

    with pytest.warns(_PydnaWarning):
        assert Dseq("GGATII").t4("g") == Dseq("")  # Dseq("gg", "", ovhg=0)

    assert Dseq("GGATCC").t4("gat") == Dseq("PPATJJ")

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

    assert t[5:1] == Dseq("")  # TODO: discuss this
    assert s[9:1] == Dseq("")
    assert t[9:1] == Dseq("")

    # Indexing of full circular molecule (https://github.com/BjornFJohansson/pydna/issues/161)
    # s = Dseq("GGATCC", circular=True) # TODO: discuss this
    # str_seq = str(s)
    # for shift in range(len(s)):
    #     assert str(s[shift:shift]) == str_seq[shift:] + str_seq[:shift]


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

    x = Dseq("TPGJG", circular=True)

    for i in range(len(x)):
        assert x.shifted(i)._data == x._data[i:] + x._data[:i]


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
        ("AAA", "TT", "EAA", (0, 1), -1),
        ("AA", "TTT", "FAA", (1, 0), 1),
        ("AAA", "TTT", "AAA", (0, 0), 0),
    ]
    for watson, crick, dsiupac, expected, ovhg in test_cases:
        dseq = Dseq(dsiupac, circular=False)
        assert dseq.left_end_position() == expected


def test_overlapping_cuts():

    import pytest
    from pydna.dseq import Dseq
    from Bio.Restriction import PstI, SbfI, BamHI, DpnI, BcgI

    s = Dseq("CCTGCAGG")

    assert s.cut(PstI) == (Dseq(b'CCXPIE'), Dseq(b'ZQJFGG'))
    assert s.cut(SbfI) == (Dseq(b'CCXPIE'), Dseq(b'ZQJFGG'))

    with pytest.raises(ValueError):
        s.cut(PstI, SbfI)

    t = Dseq("GGATCC")

    assert t.cut(BamHI) == (Dseq("GQFZJ"), Dseq("PEXIC"))
    assert t.cut(DpnI) == (Dseq("GGA"), Dseq("TCC"))

    with pytest.raises(ValueError):
        t.cut(BamHI, DpnI)

    q = Dseq("TgcgtagatcgtACGAggatccTGCGtcaagtgtctat")

    with pytest.raises(ValueError):
        q.cut(BamHI, BamHI)

    assert q.cut(BamHI, BcgI) == (Dseq("Tpi"), Dseq("qjgtagatcgtACGAgqfzj"), Dseq("pexicTGCGtcaagtgtcxe"), Dseq("zft"))


def test_apply_cut():

    from pydna.dseq import Dseq

    from Bio.Restriction import EcoRI, BamHI

    seq = Dseq("aaGAATTCaa", circular=False)

    # An empty cut returns an empty Dseq
    assert seq.apply_cut(((0, 0), None), ((0, 0), None), shift=0, stuffer=0) == Dseq("")

    # A cut where one side is None leaves that side intact
    EcoRI_cut = (3, -4), None

    assert seq.apply_cut(((0, 0), None), EcoRI_cut, shift=0, stuffer=0) == Dseq.from_full_sequence_and_overhangs(
        "aaGAATT", watson_ovhg=-4, crick_ovhg=0
    )
    assert seq.apply_cut(EcoRI_cut, ((10, 0), None), shift=0, stuffer=0) == Dseq.from_full_sequence_and_overhangs(
        "AATTCaa", watson_ovhg=0, crick_ovhg=-4
    )

    # It respects the original overhang
    seq = Dseq.from_full_sequence_and_overhangs("aaGAATTCaa", watson_ovhg=1, crick_ovhg=1)
    assert seq.apply_cut(((0, 0), None), EcoRI_cut, shift=0, stuffer=0) == Dseq.from_full_sequence_and_overhangs(
        "aaGAATT", watson_ovhg=-4, crick_ovhg=1
    )
    assert seq.apply_cut(EcoRI_cut, ((10, 0), None), shift=0, stuffer=0) == Dseq.from_full_sequence_and_overhangs(
        "AATTCaa", watson_ovhg=1, crick_ovhg=-4
    )

    seq = Dseq.from_full_sequence_and_overhangs("aaGAATTCaa", watson_ovhg=-1, crick_ovhg=-1)
    assert seq.apply_cut(((0, 0), None), EcoRI_cut, shift=0, stuffer=0) == Dseq.from_full_sequence_and_overhangs(
        "aaGAATT", watson_ovhg=-4, crick_ovhg=-1
    )
    assert seq.apply_cut(EcoRI_cut, ((10, 0), None), shift=0, stuffer=0) == Dseq.from_full_sequence_and_overhangs(
        "AATTCaa", watson_ovhg=-1, crick_ovhg=-4
    )

    # A repeated cut in a circular molecule opens it up
    seq = Dseq("aaGAATTCaa", circular=True)
    from pydna.utils import get_cutsite_pairs

    pair, shift, stuffer = get_cutsite_pairs([EcoRI_cut], circular=True, length=len(seq))
    assert pair == [(((0, -4), None), ((10, -4), None))]
    assert shift == 3
    assert stuffer == 4
    assert seq.apply_cut(*pair[0], shift=shift, stuffer=stuffer) == Dseq.from_full_sequence_and_overhangs(
        "AATTCaaaaGAATT", watson_ovhg=-4, crick_ovhg=-4
    )

    # Two cuts extract a subsequence
    seq = Dseq("aaGAATTCaaGAATTCaa", circular=True)
    EcoRI_cut_2 = ((11, -4), None)
    from pydna.utils import get_cutsite_pairs

    pair, shift, stuffer = get_cutsite_pairs([EcoRI_cut, EcoRI_cut_2], circular=True, length=len(seq))
    assert shift == 3
    assert stuffer == 4
    assert pair == [(((0, -4), None), ((8, -4), None)), (((8, -4), None), ((18, -4), None))]
    assert seq.apply_cut(*pair[0], shift=shift, stuffer=stuffer) == Dseq.from_full_sequence_and_overhangs(
        "AATTCaaGAATT", watson_ovhg=-4, crick_ovhg=-4
    )

    # Blunt cuts never overlap
    seq = Dseq("gat")
    first_cuts = [
        ((0, 0), None),
        ((1, 0), None),
        ((2, 0), None),
        ((3, 0), None),
    ]

    second_cuts = [
        ((0, 0), None),
        ((1, 0), None),
        ((2, 0), None),
        ((3, 0), None),
    ]

    for first_cut in first_cuts:
        for second_cut in second_cuts:
            try:
                seq.apply_cut(first_cut, second_cut, shift=0, stuffer=0)
            except Exception as e:
                print(e)
                pytest.fail("Unexpected error ..")

    assert seq.apply_cut(((0, 0), None), ((1, -1), None), shift=0, stuffer=0) == Dseq("gf")
    assert seq.apply_cut(((1, -1), None), ((3, 0), None), shift=0, stuffer=0) == Dseq("et")

    assert seq.apply_cut(((0, 0), None), ((1, -1), None), shift=0, stuffer=0) == Dseq("gf")
    assert seq.apply_cut(((1, -1), None), ((2, 0), None), shift=0, stuffer=0) == Dseq("e")
    assert seq.apply_cut(((2, 0), None), ((3, 0), None), shift=0, stuffer=0) == Dseq("t")

    assert seq.apply_cut(((0, 0), None), ((1, -1), None), shift=0, stuffer=0) == Dseq("gf")
    with pytest.raises(ValueError):
        seq.apply_cut(((1, -1), None), ((2, -1), None), shift=0, stuffer=0)
    assert seq.apply_cut(((2, -1), None), ((3, 0), None), shift=0, stuffer=0) == Dseq("x")

    with pytest.raises(ValueError):
        seq.apply_cut(((1, -1), None), ((2, 1), None), shift=0, stuffer=0)

    # # Rotating the sequence, apply the same cut # FIXME: reintroduce these tests
    # seq = Dseq("acgtATGaatt", circular=True)

    # for shift in range(len(seq)):
    #     seq_shifted = seq.shifted(shift)
    #     start = 4 - shift
    #     if start < 0:
    #         start += len(seq)
    #     # Cut with negative ovhg
    #     new_cut = ((start, -3), None)

    #     out = seq_shifted.apply_cut(new_cut, new_cut)
    #     assert str(out) == "ATGaattacgtATG"

    #     # Cut with positive ovhg
    #     start = (start + 3) % len(seq)
    #     new_cut = ((start, 3), None)
    #     out = seq_shifted.apply_cut(new_cut, new_cut)
    #     assert str(out) == "ATGaattacgtATG"

    #     # A blunt cut
    #     start = 4 - shift
    #     new_cut = ((start, 0), None)
    #     out = seq_shifted.apply_cut(new_cut, new_cut)
    #     assert str(out) == "ATGaattacgt"


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


def test_get_cutsite_pairs():

    from pydna.dseq import Dseq

    # in the test, we replace cuts by integers for clarity.

    from pydna.utils import get_cutsite_pairs

    # Empty returns empty list
    assert get_cutsite_pairs([], circular=None, length=None) == ([], None, None)

    # Single cut on linear seq returns two sites
    assert get_cutsite_pairs([((0, 0), None)], circular=True, length=0) == ([(((0, 0), None), ((0, 0), None))], 0, 0)

    # Two cuts on linear seq return three fragments
    assert get_cutsite_pairs([((0, 0), None), ((1, 0), None)], circular=True, length=0) == (
        [(((0, 0), None), ((1, 0), None)), (((1, 0), None), ((0, 0), None))],
        0,
        0,
    )


# def test_get_cut_parameters():  # FIXME: obsolete?

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

    fiveoh = Dseq("PEXIaaaQFZJ")
    assert str(fiveoh + fiveoh) == "PEXIaaaGATCaaaQFZJ"
    threeoh = Dseq("QFZJQtttPEXIP")
    assert str(threeoh + threeoh) == "QFZJQtttGATCGtttPEXIP"

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

    assert s._fill_in_three_prime("gatc")._data == b"GATCAAAQFZJ"
    assert s._fill_in_three_prime("gat")._data == b"PATCAAAQFZJ"
    assert s._fill_in_three_prime("ga")._data == b"PETCAAAQFZJ"
    assert s._fill_in_three_prime("g")._data == b"PEXCAAAQFZJ"
    assert s._fill_in_three_prime("")._data == s._data

    assert s._fill_in_three_prime("atc")._data == s._data
    assert s._fill_in_three_prime("at")._data == s._data
    assert s._fill_in_three_prime("a")._data == s._data
    assert s._fill_in_three_prime("")._data == s._data

    assert s._fill_in_five_prime("gatc")._data == b"PEXIAAAGATC"
    assert s._fill_in_five_prime("gat")._data == b"PEXIAAAGATJ"
    assert s._fill_in_five_prime("ga")._data == b"PEXIAAAGAZJ"
    assert s._fill_in_five_prime("g")._data == b"PEXIAAAGFZJ"
    assert s._fill_in_five_prime("")._data == s._data

    assert s._fill_in_five_prime("atc")._data == s._data
    assert s._fill_in_five_prime("at")._data == s._data
    assert s._fill_in_five_prime("a")._data == s._data
    assert s._fill_in_five_prime("")._data == s._data

    assert s.fill_in("gatc")._data == b"GATCAAAGATC"
    assert s.fill_in("gat")._data == b"PATCAAAGATJ"
    assert s.fill_in("ga")._data == b"PETCAAAGAZJ"
    assert s.fill_in("g")._data == b"PEXCAAAGFZJ"
    assert s.fill_in("")._data == s._data

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

    assert Dseq("PEXIGGATCCQFZJ").T4(b"CG") == Dseq(b"PEXCGGATCCGFZJ")
    assert Dseq("PAGAJ").T4(b"ATG") == Dseq(b"PAGAJ")

    from pydna import _PydnaWarning

    with pytest.warns(_PydnaWarning):
        assert Dseq("QFZJGGATCCPEXI").T4("") == Dseq("")

    assert Dseq("pexiAqfzj").T4(b"gatc") == Dseq("gatcAgatc")

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

    x, y = set(a)

    assert sorted(x) == sorted(y)

    assert Dseq("ggtctcAAgcTT").get_cutsites(BsaI) == [((7, -4), BsaI)]
    assert Dseq("ggtctcAAgcTT").cut(BsaI) == (Dseq("ggtctcAFqjZ"), Dseq("EpiXT"))
    assert Dseq("TggtctcAAgcT").get_cutsites(BsaI) == [((8, -4), BsaI)]
    assert Dseq("TggtctcAAgcT").cut(BsaI) == (Dseq("TggtctcAFqjZ"), Dseq("EpiX"))
    assert Dseq("TTggtctcAAgc").get_cutsites(BsaI) == []
    assert Dseq("TTggtctcAAgc").cut(BsaI) == ()

    assert Dseq("gtctcAAgcTTg", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)
    assert Dseq("ggtctcAAgcTT", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)
    assert Dseq("TggtctcAAgcT", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)
    assert Dseq("TTggtctcAAgc", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)
    assert Dseq("cTTggtctcAAg", circular=True).cut(BsaI) == (Dseq("EpiXTggtctcAFqjZ"),)


def test_melt():
    from pydna.seq import Seq
    from pydna.dseq import Dseq

    assert Dseq(b"AGJGaGEg").melt(2) == (Dseq(b"EP"), Dseq(b"FQJGaGEp"), Dseq(b"q"))
    assert Dseq(b"AGIGaGFg").melt(2) == (Dseq(b"FQ"), Dseq(b"EPIGaGFq"), Dseq(b"p"))
    assert Dseq(b"AGJGaGFg").melt(2) == (Dseq(b"EP"), Dseq(b"FQJGaGFq"), Dseq(b"p"))
    assert Dseq(b"AGIGaGEg").melt(2) == (Dseq(b"FQ"), Dseq(b"EPIGaGEp"), Dseq(b"q"))
    assert Dseq(b"GATPGCPGCA").melt(2) == (Seq('CG'), Dseq(b"GATPPIPGCA"))
    assert Dseq(b"GATQGCQGCA").melt(2) == (Seq('CG'), Dseq(b"GATQQJQGCA"))
    assert Dseq(b"PEXIGAQFZJ").melt(2) == (Dseq(b"PEXIPE"), Dseq(b"QFQFZJ"))
    assert Dseq(b"QFZJGAPEXI").melt(2) == (Dseq(b"QFZJQF"), Dseq(b"PEPEXI"))
    assert Dseq(b"AGJGaGEgGATC").melt(2) == (Dseq(b"EP"), Dseq(b"FQJGaGEgGATC"))
    assert Dseq(b"AGIGaGFgGATC").melt(2) == (Dseq(b"FQ"), Dseq(b"EPIGaGFgGATC"))
    assert Dseq(b"AGJGaGFgGATC").melt(2) == (Dseq(b"EP"), Dseq(b"FQJGaGFgGATC"))
    assert Dseq(b"AGIGaGEgGATC").melt(2) == (Dseq(b"FQ"), Dseq(b"EPIGaGEgGATC"))
    assert Dseq(b"GATPGGPGCAGATC").melt(2) == (Seq("CC"), Dseq(b"GATPPPPGCAGATC"))
    assert Dseq(b"GATQGGQGCAGATC").melt(2) == (Seq("GG"), Dseq(b"GATQQQQGCAGATC"))
    assert Dseq(b"PEXIGAQFZJGATC").melt(2) == (Dseq(b"PEXIPE"), Dseq(b"QFQFZJGATC"))
    assert Dseq(b"QFZJGAPEXIGATC").melt(2) == (Dseq(b"QFZJQF"), Dseq(b"PEPEXIGATC"))
    assert Dseq(b"GATCAGJGaGEg").melt(2) == (Dseq(b"GATCAGJGaGEp"), Dseq(b"q"))
    assert Dseq(b"GATCAGIGaGFg").melt(2) == (Dseq(b"GATCAGIGaGFq"), Dseq(b"p"))
    assert Dseq(b"GATCAGJGaGFg").melt(2) == (Dseq(b"GATCAGJGaGFq"), Dseq(b"p"))
    assert Dseq(b"GATCAGIGaGEg").melt(2) == (Dseq(b"GATCAGIGaGEp"), Dseq(b"q"))
    assert Dseq(b"GATCGATPGGPGCA").melt(2) == (Seq("CC"), Dseq(b"GATCGATPPPPGCA"))
    assert Dseq(b"GATCGATQGGQGCA").melt(2) == (Seq("GG"), Dseq(b"GATCGATQQQQGCA"))
    assert Dseq(b"GATCPEXIGAQFZJ").melt(2) == (Dseq(b"GATCPEXIPE"), Dseq(b"QFQFZJ"))
    assert Dseq(b"GATCQFZJGAPEXI").melt(2) == (Dseq(b"GATCQFZJQF"), Dseq(b"PEPEXI"))
    assert Dseq(b"GATCAGJGaGEgGATC").melt(2) == ()
    assert Dseq(b"GATCAGIGaGFgGATC").melt(2) == ()
    assert Dseq(b"GATCAGJGaGFgGATC").melt(2) == ()
    assert Dseq(b"GATCAGIGaGEgGATC").melt(2) == ()
    assert Dseq(b"GATCGATPGGPGCAGATC").melt(2) == (Seq("CC"), Dseq(b"GATCGATPPPPGCAGATC"))
    assert Dseq(b"GATCGATQGGQGCAGATC").melt(2) == (Seq("GG"), Dseq(b"GATCGATQQQQGCAGATC"))
    assert Dseq(b"GATCPEXIGAQFZJGATC").melt(2) == (Dseq(b"GATCPEXIPE"), Dseq(b"QFQFZJGATC"))
    assert Dseq(b"GATCQFZJGAPEXIGATC").melt(2) == (Dseq(b"GATCQFZJQF"), Dseq(b"PEPEXIGATC"))
    assert Dseq(b"GATCPEXIGAQFZJGATC").melt(2) == (Dseq(b"GATCPEXIPE"), Dseq(b"QFQFZJGATC"))
    assert Dseq(b"GATCQFZJGAPEXIGATC").melt(2) == (Dseq(b"GATCQFZJQF"), Dseq(b"PEPEXIGATC"))


def test__get_ds_meltsites():

    from pydna.dseq import Dseq

    assert Dseq(b"AGJGaGEg").get_ds_meltsites(2) == [((2, 2), None), ((8, 1), None)]
    assert Dseq(b"AGIGaGFg").get_ds_meltsites(2) == [((0, -2), None), ((7, -1), None)]
    assert Dseq(b"AGJGaGFg").get_ds_meltsites(2) == [((2, 2), None), ((7, -1), None)]
    assert Dseq(b"AGIGaGEg").get_ds_meltsites(2) == [((0, -2), None), ((8, 1), None)]

    assert Dseq(b"PEXIGAQFZJ").get_ds_meltsites(2) == [((6, 2), None)]
    assert Dseq(b"QFZJGAPEXI").get_ds_meltsites(2) == [((4, -2), None)]

    assert Dseq(b"GATCPEXIGAQFZJGATC").get_ds_meltsites(2) == [((10, 2), None)]
    assert Dseq(b"GATCQFZJGAPEXIGATC").get_ds_meltsites(2) == [((8, -2), None)]

    assert Dseq("AGCPAGQGAT", circular=True).get_ds_meltsites(2) == [((6, 2), None)]
    assert Dseq("AGCQAGPGAT", circular=True).get_ds_meltsites(2) == [((4, -2), None)]


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


def test_anneal():
    from pydna.dseq import Dseq

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


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
