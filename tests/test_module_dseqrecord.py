#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from pydna import _PydnaWarning
from Bio import BiopythonWarning

# from Bio.SeqIO.InsdcIO import BiopythonWarning

import warnings

warnings.simplefilter("always", BiopythonWarning)


def test_new():
    from Bio.Restriction import Acc65I, KpnI
    from pydna.dseqrecord import Dseqrecord
    from Bio.SeqFeature import ExactPosition
    from Bio.SeqFeature import SimpleLocation

    GGTACC = Dseqrecord("GGTACC")
    GGTACC.add_feature(0, 5)
    GGTACC.add_feature(1, 5)
    GGTACC.add_feature(1, 2)
    GGTACC.add_feature(4, 5)
    GGTACC.add_feature(1, 6)

    GQZFJ, PXEIC = GGTACC.cut(Acc65I)
    GPXEI, QZFJC = GGTACC.cut(KpnI)

    assert GQZFJ.features == GGTACC.features[:4]
    assert GPXEI.features == GGTACC.features[:4]
    assert [f.location for f in PXEIC.features] == [f.location - 1 for f in GGTACC.features[1:]]
    assert [f.location for f in QZFJC.features] == [f.location - 1 for f in GGTACC.features[1:]]

    oGGTACC = Dseqrecord("GGTACC", circular=True)
    oGGTACC.add_feature(1, 5)
    (PXEICGQZFJ,) = oGGTACC.cut(Acc65I)
    assert PXEICGQZFJ.features[0].location == SimpleLocation(ExactPosition(0), ExactPosition(4), strand=1)

    oGGTACC = Dseqrecord("GGTACC", circular=True)
    oGGTACC.add_feature(1, 2)
    (PXEICGQZFJ,) = oGGTACC.cut(Acc65I)
    assert PXEICGQZFJ.features[0].location == SimpleLocation(ExactPosition(0), ExactPosition(1), strand=1)

    oGGTACC = Dseqrecord("GGTACC", circular=True)
    oGGTACC.add_feature(4, 5)
    (PXEICGQZFJ,) = oGGTACC.cut(Acc65I)
    assert PXEICGQZFJ.features[0].location == SimpleLocation(ExactPosition(3), ExactPosition(4), strand=1)

    oGGTACC = Dseqrecord("GGTACC", circular=True)
    oGGTACC.add_feature(0, 6)
    (PXEICGQZFJ,) = oGGTACC.cut(Acc65I)
    assert PXEICGQZFJ.features == []

    oGGTACC = Dseqrecord("GGTACC", circular=True)
    oGGTACC.add_feature(1, 5)
    (PXEICGQZFJ,) = oGGTACC.cut(KpnI)
    assert PXEICGQZFJ.features[0].location == SimpleLocation(ExactPosition(0), ExactPosition(4), strand=1)

    oGGTACC = Dseqrecord("GGTACC", circular=True)
    oGGTACC.add_feature(1, 2)
    (PXEICGQZFJ,) = oGGTACC.cut(KpnI)
    assert PXEICGQZFJ.features[0].location == SimpleLocation(ExactPosition(0), ExactPosition(1), strand=1)

    oGGTACC = Dseqrecord("GGTACC", circular=True)
    oGGTACC.add_feature(4, 5)
    (PXEICGQZFJ,) = oGGTACC.cut(KpnI)
    assert PXEICGQZFJ.features[0].location == SimpleLocation(ExactPosition(3), ExactPosition(4), strand=1)

    oGGTACC = Dseqrecord("GGTACC", circular=True)
    oGGTACC.add_feature(0, 6)
    (PXEICGQZFJ,) = oGGTACC.cut(KpnI)
    assert PXEICGQZFJ.features == []

    oGGTACC = Dseqrecord("GGTACC", circular=True)
    oGGTACC.add_feature(0, 5)
    (PXEICGQZFJ,) = oGGTACC.cut(KpnI)
    assert PXEICGQZFJ.features[0].location == SimpleLocation(ExactPosition(5), ExactPosition(10), strand=1)

    s = """

    LOCUS       New_DNA                   13 bp    DNA     circular     22-DEC-2024
    DEFINITION  .
    ACCESSION
    VERSION
    SOURCE      .
      ORGANISM  .
    FEATURES             Location/Qualifiers
         misc_feature    2..8
                         /label="in"

         misc_feature    join(9..13,1..1)
                         /label="bb"
    ORIGIN
            1 CCGGaaaCCG Gct
    //

    """

    from pydna.readers import read

    dsr = read(s)

    dsr.seq

    from Bio.Restriction import HpaII

    a, b = dsr.cut(HpaII)


def test_orfs():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("AatgaaataaT")

    assert s.orfs(9)[0].seq == s.seq[1:-1]
    assert s.orfs(10) == tuple()

    s = Dseqrecord(
        "tctgcaataATGGGTAATGAAATCGATGAGAAAAATCAGGCCCCCGTGCAACAAGAATGCCTGAAAGAGATGATTCAGAATGGGCATGCTCGGCGTATGGGATCTGTTGAAGATCTGTATGTTGCTCTCAACAGACAAAACTTATATCGAAACTTCTGCACATATGGAGAATTGAGTGATTACTGTACTAGGGATCAGCTCACATTAGCTTTGAGGGAAATCTGCCTGAAAAATCCAACTCTTTTACATATTGTTCTACCAACAAGATGGCCAAATCATGAAAATTATTATCGCAGTTCCGAATACTATTCACGGCCACATCCAGTGCATGATTATATTTCAGTATTACAAGAATTGAAACTGAGTGGTGTGGTTCTCAATGAACAACCTGAGTACAGTGCAGTAATGAAGCAAATATTAGAAGAATTCAAAAATAGTAAGGGTTCCTATACTGCAAAAATTTTTAAACTTACTACCACTTTGACTATTCCTTACTTTGGACCAACAGGACCGAGTTGGCGGCTAATTTGTCTTCCAGAAGAGCACACAGAAAAGTGGAAAAAATTTATCTTTGTATCTAATCATTGCATGTCTGATGGTCGGTCTTCGATCCACTTTTTTCATGATTTAAGAGACGAATTAAATAATATTAAAACTCCACCAAAAAAATTAGATTACATTTTCAAGTACGAGGAGGATTACCAATTATTGAGGAAACTTCCAGAACCGATCGAAAAGGTGATAGACTTTAGACCACCGTACTTGTTTATTCCGAAGTCACTTCTTTCGGGTTTCATCTACAATCATTTGAGATTTTCTTCAAAAGGTGTCTGTATGAGAATGGATGATGTGGAAAAAACCGATGATGTTGTCACCGAGATCATCAATATTTCACCAACAGAATTTCAAGCGATTAAAGCAAATATTAAATCAAATATCCAAGGTAAGTGTACTATCACTCCGTTTTTACATGTTTGTTGGTTTGTATCTCTTCATAAATGGGGTAAATTTTTCAAACCATTGAACTTCGAATGGCTTACGGATATTTTTATCCCCGCAGATTGCCGCTCACAACTACCAGATGATGATGAAATGAGACAGATGTACAGATATGGCGCTAACGTTGGATTTATTGACTTCACCCCCTGGATAAGCGAATTTGACATGAATGATAACAAAGAAAATTTTTGGCCACTTATTGAGCACTACCATGAAGTAATTTCGGAAGCTTTAAGAAATAAAAAGCATCTCCATGGCTTAGGGTTCAATATACAAGGCTTCGTTCAAAAATATGTGAACATTGACAAGGTAATGTGCGATCGTGCCATCGGGAAAAGACGCGGAGGTACATTGTTAAGCAATGTAGGTCTGTTTAATCAGTTAGAGGAGCCCGATGCCAAATATTCTATATGCGATTTGGCATTTGGCCAATTTCAAGGATCCTGGCACCAAGCATTTTCCTTGGGTGTTTGTTCGACTAATGTAAAGGGGATGAATATTGTTGTTGCTTCAACAAAGAATGTTGTTGGTAGTCAAGAATCTCTCGAAGAGCTTTGCTCCATTTACAAAGCTCTCCTTTTAGGCCCTTAA"
    )

    lens = (1581, 1521, 1494, 1185, 1002)
    for orf, ln in zip(s.orfs(1002), lens):
        assert len(orf) == ln


def test_cas9():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("gattcatgcatgtagcttacgtagtct")

    RNA = "catgcatgtagcttacgtag"

    (f1, f2), (f3,) = s.cas9(RNA)

    assert f1.seq == Dseqrecord("gattcatgcatgtagcttacg").seq
    assert f2.seq == Dseqrecord("tagtct").seq
    assert f3.seq == Dseqrecord("gattcatgcatgtagcttacgtagtct").seq


def test_FadiBakoura():
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord

    dseq_record = Dseqrecord(Dseq("ACTCTTTCTCTCTCT", circular=True))
    dseq_record.features = [SeqFeature(FeatureLocation(start=2, end=4))]
    # below assert fails (not anymore)
    assert len(dseq_record[6:1].features) == 0


def test_IPython_missing(monkeypatch):
    import IPython
    from pydna import dseqrecord

    # assert dseqrecord._display is IPython.display.display
    assert dseqrecord._display_html == IPython.display.display_html
    import sys
    import copy

    fakesysmodules = copy.copy(sys.modules)
    fakesysmodules["IPython.display"] = None
    monkeypatch.delitem(sys.modules, "IPython.display")
    monkeypatch.setattr("sys.modules", fakesysmodules)
    from importlib import reload

    reload(dseqrecord)
    from pydna import dseqrecord

    assert dseqrecord._display_html("item") == "item"
    # assert dseqrecord._HTML("item") == "item"


def test_initialization():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec

    a = []

    a.append(Dseqrecord("attt", id="1"))
    a.append(Dseqrecord(Dseq("attt"), id="2"))
    a.append(Dseqrecord(Seq("attt"), id="3"))
    a.append(Dseqrecord(Srec(Seq("attt"), id="4")))
    a.append(Dseqrecord(Dseqrecord("attt"), id="5"))
    a.append(Dseqrecord(Dseqrecord(Dseq("attt", circular=True)), id="6", circular=False))

    for b in a:
        assert isinstance(b.seq, Dseq)
        assert str(b.seq.watson) == "attt"
        assert str(b.seq.crick) == "aaat"
        assert str(b.seq) == "attt"
        assert str(b.seq) == "attt"
        assert b.seq.circular == False
        assert b.circular == False

    a = []
    a.append(Dseqrecord("attt", circular=True))
    a.append(Dseqrecord(Dseq("attt"), circular=True))
    a.append(Dseqrecord(Seq("attt"), circular=True))
    a.append(Dseqrecord(Srec(Seq("attt")), circular=True))
    a.append(Dseqrecord(Dseqrecord("attt"), circular=True))

    for b in a:
        assert isinstance(b.seq, Dseq)
        assert str(b.seq.watson) == "attt"
        assert str(b.seq.crick) == "aaat"
        assert str(b.seq) == "attt"
        assert str(b.seq) == "attt"
        assert b.circular == True
        assert b.seq.circular == True

    a = []
    a.append(Dseqrecord(Dseq("attt", circular=True), circular=True))
    a.append(Dseqrecord(Dseq("attt", circular=False), circular=True))
    a.append(Dseqrecord(Dseq("attt", circular=True), circular=False))
    a.append(Dseqrecord(Dseq("attt", circular=False), circular=False))

    circular = [True, True, False, False]

    for b, ci in zip(a, circular):
        assert isinstance(b.seq, Dseq)
        assert str(b.seq.watson) == "attt"
        assert str(b.seq.crick) == "aaat"
        assert str(b.seq) == "attt"
        assert str(b.seq) == "attt"
        assert b.circular == ci
        assert b.seq.circular == ci

    a = []
    ds = Dseq("etttf")
    assert ds.ovhg == -1
    assert str(ds.watson) == "attt"
    assert str(ds.crick) == "taaa"

    #   attt
    #    aaat

    dsr = Dseqrecord(ds, circular=False)

    assert isinstance(dsr.seq, Dseq)
    assert dsr.seq.watson == "attt"
    assert dsr.seq.crick == "taaa"
    assert dsr.circular == False
    assert dsr.seq.circular == False
    assert str(dsr.seq) == "etttf"

    dsr = Dseqrecord(ds, circular=True)

    assert isinstance(dsr.seq, Dseq)
    assert dsr.seq.watson == "attt"
    assert dsr.seq.crick == "aaat"
    assert dsr.circular == True
    assert dsr.seq.circular == True

    a = []
    ds = Dseq("etttq")
    assert ds.circular == False
    assert ds.ovhg == -1

    a.append(Dseqrecord(ds, circular=False))

    with pytest.raises(TypeError):
        Dseqrecord(ds, circular=True)

    with pytest.raises(TypeError):
        Dseqrecord(ds, circular=True)

    with pytest.raises(ValueError):
        b = Dseqrecord([])

    with pytest.raises(ValueError):
        b = Dseqrecord(("a",))

    with pytest.raises(ValueError):
        b = Dseqrecord(0)

    input = """
            LOCUS       New_DNA                    4 bp ds-DNA     linear       30-MAR-2013
            DEFINITION  .
            ACCESSION
            VERSION
            SOURCE      .
              ORGANISM  .
            COMMENT
            COMMENT     ApEinfo:methylated:1
            FEATURES             Location/Qualifiers
                 misc_feature    2..3
                                 /label=NewFeature
                                 /ApEinfo_fwdcolor=cyan
                                 /ApEinfo_revcolor=green
                                 /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                 width 5 offset 0
            ORIGIN
                    1 acgt
            //
            """
    a = read(input)

    assert a.features[0].extract(a).seq.watson == "CG"

    b = a + a

    for f in b.features:
        assert b.features[0].extract(a).seq.watson == "CG"
    feature = a.features[0]

    s = Dseq("zagctx")
    b = Dseqrecord(s)
    b.features.append(feature)
    cb = Dseqrecord(b, circular=True)
    assert b.features[0].extract(b).seq.watson.lower() == cb.features[0].extract(b).seq.watson.lower()
    assert b.features[0].extract(b).seq.crick.lower() == cb.features[0].extract(b).seq.crick.lower()
    s = Dseq("eagctz")
    # print s.fig()
    # Dseq(-6)
    # aagct
    # tcgaa
    b = Dseqrecord(s)
    with pytest.raises(TypeError):
        cb = Dseqrecord(b, circular=True)

    s = Dseq("zagctx")
    # print s.fig()
    # Dseq(-6)
    # agcta
    # ttcga

    b = Dseqrecord(s)
    b.features.append(feature)
    cb = Dseqrecord(b, circular=True)
    assert b.features[0].extract(b).seq.watson.lower() == cb.features[0].extract(b).seq.watson.lower()
    assert b.features[0].extract(b).seq.crick.lower() == cb.features[0].extract(b).seq.crick.lower()


def test_linear_circular():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec

    """ test Dseqrecord linear & circular property"""
    a = Dseqrecord("attt")
    a.stamp()
    assert a.stamp()
    a = Dseqrecord("attt", circular=False)
    assert a.circular == False
    assert a.rc().circular == False

    a = Dseqrecord("attt", circular=True)

    assert a.circular == True
    assert a.rc().circular == True
    assert a.seq.circular == True

    a = Dseqrecord("attt", circular=True)

    assert a.circular == True
    assert a.rc().circular == True
    assert a.seq.circular == True

    a = Dseqrecord("attt", circular=False)

    assert a.circular == False
    assert a.rc().circular == False
    assert a.seq.circular == False


def test_stamp():
    from pydna.dseqrecord import Dseqrecord

    lin = Dseqrecord("attt")
    lin.stamp()
    first = lin.annotations["comment"]
    assert "ldseguid=BPR1_ojRL1ZJa6EgF02NLHAxr68" in first
    lin.stamp()
    assert first[:42] == lin.annotations["comment"][:42]

    crc = Dseqrecord("attt", circular=True)
    crc.stamp()
    first = crc.annotations["comment"]
    assert "cdseguid=BPR1_ojRL1ZJa6EgF02NLHAxr68" in first
    crc.stamp()
    assert first[:42] == crc.annotations["comment"][:42]
    assert len(first) == len(crc.annotations["comment"])

    from pydna.dseq import Dseq

    blunt = Dseqrecord(Dseq("aa"))

    assert blunt.stamp()[:42] == "ldseguid=TEwydy0ugvGXh3VJnVwgtxoyDQA"

    staggered = Dseqrecord(Dseq("zaa"))
    assert staggered.stamp()[:42] == "ldseguid=WPLhxEZErSzQmVMmVhZrQ5aSc78"

    staggered = Dseqrecord(Dseq("aaz"))
    assert staggered.stamp()[:42] == "ldseguid=Vma2bZhvSl9otSfAvTQP5eUsXYY"

    staggered = Dseqrecord(Dseq("zaaz"))
    assert staggered.stamp()[:42] == "ldseguid=8Fy5Jaz0IKJ_I4cvAFUj0XX718g"


def test_revcomp():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec

    #      ----
    #     attcccgggg
    #     taagggcccc
    #
    #      ttcc
    #
    #
    #     ccccgggaat
    #     ggggccctta
    #          ----

    a = Dseqrecord("attcccgggg")

    a.add_feature(1, 5)

    assert str(a.features[0].extract(a).seq) == "ttcc"

    assert a.features[0].location.strand == 1

    rc = a.rc()

    assert str(rc.features[0].extract(rc).seq) == "ttcc"

    assert rc.features[0].location.strand == -1

    a = Dseqrecord("attcccgggg")

    a.add_feature(1, 5)

    a.features[0].location.strand = None

    rc = a.rc()

    assert rc.features[0].location.strand is None


def test_m():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("A" * 5000)
    assert f"{s.m():.3e}" == "1.544e-07"


def test_extract_feature():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("tttGGATCCaaa")
    s.add_feature(3, 9)
    s.extract_feature(0).seq == Dseqrecord("GGATCC").seq


def test_find():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("tttGGATCCaaa")
    assert s.find(Dseqrecord("ggatcc")) == 3
    assert s.find(Dseqrecord("ggatTc")) == -1
    s = Dseqrecord("tttGGATCCaaa", circular=True)
    assert s.find(Dseqrecord("aaattt")) == 9


def test_find_aa():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("tttGGATCCaaa")
    assert s.find_aa("FGSK") == slice(0, 12, None)
    assert s.find_aa("MMMM") is None


def test_str():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("tttGGATCCaaa")
    s.annotations = {"date": "03-JAN-2018"}
    assert (
        str(s)
        == "Dseqrecord\ncircular: False\nsize: 12\nID: id\nName: name\nDescription: description\nNumber of features: 0\n/date=03-JAN-2018\nDseq(-12)\ntttGGATCCaaa\naaaCCTAGGttt"
    )
    s = s.looped()
    assert (
        str(s)
        == "Dseqrecord\ncircular: True\nsize: 12\nID: id\nName: name\nDescription: description\nNumber of features: 0\n/date=03-JAN-2018\nDseq(o12)\ntttGGATCCaaa\naaaCCTAGGttt"
    )


def test___contains__():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("tttGGATCCaaa")
    assert "ggatcc" in s


def test_seguid():
    from pydna.dseqrecord import Dseqrecord

    l = Dseqrecord("tttGGATCCaaa")
    assert l.seguid() == "ldseguid=jbGRr-Jhpl0tVyt0Bx5nmY9_G6E"
    c = Dseqrecord("tttGGATCCaaa", circular=True)
    assert c.seguid() == "cdseguid=r5dYgWx-W5r1KxnmIqMA19t9Hh8"


def test_format():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=True)
    s.format("gb")
    s.format("genbank")
    s = Dseqrecord("GGATCC", circular=False)
    s.format("gb")
    s.format("genbank")

    s = Dseqrecord("GGATCC", circular=True)
    s.format("fasta")
    s = Dseqrecord("GGATCC", circular=False)
    s.format("fasta")

    # A long sequence name is properly handled
    s.name = "A" * 45
    genbank_str = s.format("genbank")
    locus_line = genbank_str.split("\n")[0]

    assert (
        locus_line
        == "LOCUS       AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA 6 bp    DNA     linear   UNK 01-JAN-1980"
    )

    # with pytest.warns((BiopythonWarning,)):
    #     assert (locus_line
    #         == "LOCUS       AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA 6 bp    DNA     linear   UNK 01-JAN-1980")

    # with warnings.catch_warnings() as w:
    #     #warnings.simplefilter("always")
    #     assert (locus_line
    #         == "LOCUS       AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA 6 bp    DNA     linear   UNK 01-JAN-1980")
    #     # assert any(isinstance(warning.message, BiopythonWarning) for warning in w), "BiopythonWarning not emitted"


def test_write():
    from unittest.mock import patch
    from unittest.mock import mock_open
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=True)
    m = mock_open()
    with patch("pydna.dseqrecord.open", m):
        s.write(filename="AAA.gb")
    m.mock_calls
    m.assert_called_once_with("AAA.gb", "w", encoding="utf8")
    handle = m()
    handle.write.assert_called_once_with(Dseqrecord("GGATCC", circular=True).format())

    m = mock_open()
    with patch("pydna.dseqrecord.open", m):
        s.write()
    m.mock_calls
    m.assert_called_once_with("name.gb", "w", encoding="utf8")
    handle = m()
    handle.write.assert_called_once_with(Dseqrecord("GGATCC", circular=True).format())

    with pytest.raises(TypeError):
        s.write(filename=123)


def test_write_content_to_existing_file_with_same_content():

    from pytest import MonkeyPatch
    from pydna.dseqrecord import Dseqrecord
    from unittest.mock import mock_open

    m_open = mock_open()
    s = Dseqrecord("Ggatcc", circular=True)

    def mock_display_html(arg, raw=None):
        return arg, raw

    with MonkeyPatch.context() as mp:
        mp.setattr("pydna.dseqrecord._os.path.isfile", lambda x: True)
        mp.setattr("pydna.dseqrecord._os.path.getmtime", lambda x: 0.0)
        mp.setattr("pydna.readers.read", lambda x: s)
        mp.setattr("pydna.dseqrecord._display_html", mock_display_html)
        mp.setattr("builtins.open", m_open)
        msg, raw = s.write(filename="Ggatcc.gb")
        assert raw is True
        assert msg == "<font face=monospace><a href='Ggatcc.gb' target='_blank'>Ggatcc.gb</a></font><br>"
        m_open.assert_not_called()


def test_write_content_to_existing_file_with_different_content():

    from pytest import MonkeyPatch
    from pydna.dseqrecord import Dseqrecord
    from unittest.mock import mock_open

    m_open = mock_open()
    s = Dseqrecord("Ggatcc", circular=True)
    d = Dseqrecord("GgatcA", circular=True)

    def mock_display_html(arg, raw=None):
        return arg, raw

    data = "LOCUS       name                       6 bp    DNA     circular UNK 01-JAN-1980\nDEFINITION  description.\nACCESSION   id\nVERSION     id\nKEYWORDS    .\nSOURCE      .\n  ORGANISM  .\n            .\nFEATURES             Location/Qualifiers\nORIGIN\n        1 ggatcc\n//"

    with MonkeyPatch.context() as mp:
        mp.setattr("pydna.dseqrecord._os.path.isfile", lambda x: True)
        mp.setattr("pydna.dseqrecord._os.path.getmtime", lambda x: 0.0)
        mp.setattr("pydna.dseqrecord._os.rename", lambda x, y: True)
        mp.setattr("pydna.readers.read", lambda x: d)
        mp.setattr("pydna.dseqrecord._display_html", mock_display_html)
        mp.setattr("builtins.open", m_open)
        s.write(filename="Ggatcc.gb")
        msg, raw = s.write(filename="Ggatcc.gb")
        assert raw is True
        assert "Ggatcc_OLD_" in msg
        m_open.assert_called_with("Ggatcc.gb", "w", encoding="utf8")
        assert m_open().write.call_args[0][0] == data


def test_write_content_to_stamped_existing_file_with_different_content():

    from pytest import MonkeyPatch
    from pydna.dseqrecord import Dseqrecord
    from unittest.mock import mock_open
    from textwrap import dedent

    data = dedent(
        """\
    LOCUS       name                       6 bp    DNA     circular UNK 01-JAN-1980
    DEFINITION  description.
    ACCESSION   id
    VERSION     id
    KEYWORDS    .
    SOURCE      .
      ORGANISM  .
                .
    COMMENT     pydna cdseguid=NNNNNNNNNNNNNNNNNNNNNNNNNNN
                pydna cdseguid=-REF1S8RXPtqi-FtmqexGtGvvSo 2025-01-09T08:27:12
    FEATURES             Location/Qualifiers
    ORIGIN
            1 ggatcc
    //"""
    )

    new = Dseqrecord("Ggatcc", circular=True)
    old = Dseqrecord("Ggatcc", circular=True)

    old.annotations["comment"] = "pydna cdseguid=-REF1S8RXPtqi-FtmqexGtGvvSo 2025-01-09T08:27:12"
    new.annotations["comment"] = "pydna cdseguid=NNNNNNNNNNNNNNNNNNNNNNNNNNN"

    m_open = mock_open(read_data=old.format())

    def mock_display_html(arg, raw=None):
        return arg, raw

    with MonkeyPatch.context() as mp:
        mp.setattr("pydna.dseqrecord._os.path.isfile", lambda x: True)
        mp.setattr("pydna.dseqrecord._os.path.getmtime", lambda x: 0.0)
        mp.setattr("pydna.dseqrecord._os.rename", lambda x, y: True)
        mp.setattr("pydna.readers.read", lambda x: old)
        mp.setattr("pydna.dseqrecord._display_html", mock_display_html)
        mp.setattr("builtins.open", m_open)
        msg, raw = new.write(filename="Ggatcc.gb")
        assert raw is True
        assert msg == "<font face=monospace><a href='Ggatcc.gb' target='_blank'>Ggatcc.gb</a></font><br>"
        m_open.assert_called_with("Ggatcc.gb", "w", encoding="utf8")
        assert m_open().write.call_args[0][0] == data


def test_cut_args():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC")
    from Bio.Restriction import BamHI, BglII, RestrictionBatch

    rb = RestrictionBatch((BamHI, BglII))
    assert s.cut(BamHI)[0].seq == s.cut(BamHI + BglII)[0].seq == s.cut(rb)[0].seq
    assert s.cut(BamHI)[1].seq == s.cut(BamHI + BglII)[1].seq == s.cut(rb)[1].seq


def test_cut_circular():
    from pydna.dseqrecord import Dseqrecord
    from Bio.Restriction import BsaI, KpnI, Acc65I, NotI

    test = "aaGGTACCggtctcaaaa"

    for i in range(len(test)):
        nt = test[i:] + test[:i]

        a = Dseqrecord(nt, circular=True).cut(Acc65I)[0]
        assert a.seq.watson.upper() == "GTACCGGTCTCAAAAAAG"
        assert a.seq.crick.upper() == "GTACCTTTTTTGAGACCG"
        assert a.seq.ovhg == -4
        b = Dseqrecord(nt, circular=True).cut(KpnI)[0]
        assert b.seq.watson.upper() == "CGGTCTCAAAAAAGGTAC"
        assert b.seq.crick.upper() == "CTTTTTTGAGACCGGTAC"
        assert b.seq.ovhg == 4
        c = Dseqrecord(nt, circular=True).cut(BsaI)[0]
        assert c.seq.watson.upper() == "AAAAAGGTACCGGTCTCA"
        assert c.seq.crick.upper() == "TTTTTGAGACCGGTACCT"
        assert c.seq.ovhg == -4
        d = Dseqrecord(nt, circular=True).cut(NotI)
        assert d == ()


def test_cut_add():
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec

    from Bio.Seq import Seq
    from Bio.Restriction import BamHI, EcoRI, PstI, EcoRV, SmaI

    from Bio.SeqUtils.CheckSum import seguid

    a = Dseqrecord("GGATCCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtGGATCC").seq
    b = a.cut(BamHI)[1]
    c = Dseqrecord("nCTGCAGtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtGAATTCn").seq
    f, d, l = c.cut((EcoRI, PstI))

    pUC19 = read("pUC19.gb")

    assert pUC19.circular == True

    assert len(pUC19) == 2686
    assert len(pUC19.seq.watson) == 2686
    assert len(pUC19.seq.crick) == 2686

    assert pUC19.seq.circular == True

    pUC19_SmaI = pUC19.cut(SmaI)
    assert len(pUC19_SmaI) == 1
    pUC19_SmaI = pUC19_SmaI[0]

    assert not pUC19_SmaI.circular
    assert len(pUC19_SmaI) == 2686

    pUC19_SmaI_a = pUC19_SmaI.seq + a

    assert not pUC19_SmaI_a.circular
    assert pUC19_SmaI_a.circular == False

    pUC19_SmaI_a = pUC19_SmaI_a.looped()
    assert len(pUC19_SmaI_a) == 2778

    assert pUC19_SmaI_a.circular
    assert eq(pUC19_SmaI_a, read("pUC19-SmaI-a.gb"))

    """ sticky end cloning """

    pUC19_BamHI = pUC19.cut(BamHI)

    assert len(pUC19_BamHI) == 1

    pUC19_BamHI = pUC19_BamHI[0].seq

    assert len(pUC19_BamHI.watson) == len(pUC19_BamHI.crick) == 2686

    pUC19_BamHI_a = pUC19_BamHI + b

    assert len(pUC19_BamHI_a.watson) == len(pUC19_BamHI_a.crick) == 2772

    assert pUC19_BamHI_a.circular == False

    pUC19_BamHI_a = pUC19_BamHI_a.looped()

    assert pUC19_BamHI_a.circular == True

    assert eq(pUC19_BamHI_a, read("pUC19-BamHI-a.gb"))

    pUC19_BamHI_a_rc = pUC19_BamHI + b.rc()

    pUC19_BamHI_a_rc = pUC19_BamHI_a_rc.looped()

    assert pUC19_BamHI_a.circular == True

    assert eq(pUC19_BamHI_a_rc, read("pUC19-BamHI-a-rc.gb"))

    """ adding (ligating) dsDNA objects """
    with pytest.raises(TypeError) as excinfo:
        pUC19 + a
    assert "circular" in str(excinfo.value)
    with pytest.raises(TypeError) as excinfo:
        a + pUC19
    assert "circular" in str(excinfo.value)
    with pytest.raises(TypeError) as excinfo:
        a + b
    assert "compatible" in str(excinfo.value)
    with pytest.raises(TypeError) as excinfo:
        b + a
    assert "compatible" in str(excinfo.value)
    with pytest.raises(TypeError) as excinfo:
        d + d
    assert "compatible" in str(excinfo.value)

    """ directional cloning """

    pUC19_EcoRI_PstI = pUC19.cut(EcoRI, PstI)[1]

    with pytest.raises(TypeError) as excinfo:
        pUC19_EcoRI_PstI + d
    assert "compatible" in str(excinfo.value)

    pUC19_EcoRI_PstI_d = pUC19_EcoRI_PstI + d.rc()

    pUC19_EcoRI_PstI_d = pUC19_EcoRI_PstI_d.looped()

    assert eq(pUC19_EcoRI_PstI_d, read("pUC19-EcoRI_PstI-d-rc.gb"))
    assert eq(pUC19_EcoRI_PstI_d.rc(), read("pUC19-EcoRI_PstI-d-rc.gb"))

    from Bio.Restriction import Bsu36I, BstAPI

    pCAPs = read("pCAPs.gb")
    a, b = pCAPs.cut(Bsu36I, BstAPI)
    c = (a + b).looped()
    assert eq(c, pCAPs)


def test_Dseqrecord_cutting_adding_2():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec

    a = (
        Dseqrecord(Dseq("EEXXCACANGGTACCNGGTACCNGCGGATATC")),
        Dseqrecord(Dseq("CACANGGTACCNGGTACCNGCGGATATC")),
        Dseqrecord(Dseq("ZZFFCACANGGTACCNGGTACCNGCGGATATC")),
    )

    from Bio.Restriction import KpnI, Acc65I, NlaIV

    enzymes = [Acc65I, NlaIV, KpnI]

    for enz in enzymes:
        for f in a:
            b, c, d = f.cut(enz)
            # print(b.seq.__repr__())
            # print(c.seq.__repr__())
            # print(d.seq.__repr__())
            e = b + c + d
            assert str(e.seq).lower() == str(f.seq).lower()


def test_Dseqrecord_cutting_adding_3():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec
    from Bio.Restriction import KpnI, BamHI, Acc65I, NlaIV, EcoRI, EcoRV

    a = read(
        """

LOCUS       New_DNA                   10 bp ds-DNA     linear       02-APR-2013
DEFINITION
ACCESSION   New_DNA
VERSION     New_DNA
KEYWORDS    .
SOURCE
  ORGANISM  . .
COMMENT
COMMENT     ApEinfo:methylated:1
FEATURES             Location/Qualifiers
     misc_feature    1..1
                     /label=1
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    2..2
                     /label=2
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    3..3
                     /label=3
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    4..4
                     /label=4
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    5..5
                     /label=5
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    6..6
                     /label=6
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    7..7
                     /label=7
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    8..8
                     /label=8
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    9..9
                     /label=9
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    10..10
                     /label=10
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
ORIGIN
        1 ttGGTACCgg
//"""
    )

    b, c = a.cut(Acc65I)

    assert [f.qualifiers["label"] for f in b.features] == [
        ["1"],
        ["2"],
        ["3"],
        ["4"],
        ["5"],
        ["6"],
        ["7"],
    ]
    assert [f.qualifiers["label"] for f in c.features] == [
        ["4"],
        ["5"],
        ["6"],
        ["7"],
        ["8"],
        ["9"],
        ["10"],
    ]


def test_Dseqrecord_cutting_adding_4():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec
    from Bio.Restriction import KpnI, Acc65I, NlaIV, EcoRI, EcoRV

    a = read(
        """
    LOCUS       New_DNA                   33 bp ds-DNA     linear       08-NOV-2012
    DEFINITION  .
    ACCESSION   .
    VERSION     .
    SOURCE      .
      ORGANISM  .
    COMMENT
    COMMENT     ApEinfo:methylated:1
    FEATURES             Location/Qualifiers
         misc_feature    1..11
                         /label=Acc65I-1
                         /ApEinfo_fwdcolor=cyan
                         /ApEinfo_revcolor=green
                         /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                         width 5 offset 0
         misc_feature    12..18
                         /label=Acc65I-2
                         /ApEinfo_fwdcolor=cyan
                         /ApEinfo_revcolor=green
                         /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                         width 5 offset 0
         misc_feature    19..33
                         /label=Acc65I-3
                         /ApEinfo_fwdcolor=cyan
                         /ApEinfo_revcolor=green
                         /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                         width 5 offset 0
         misc_feature    1..15
                         /label=KpnI-1
                         /ApEinfo_fwdcolor=cyan
                         /ApEinfo_revcolor=green
                         /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                         width 5 offset 0
         misc_feature    16..22
                         /label=KpnI-2
                         /ApEinfo_fwdcolor=cyan
                         /ApEinfo_revcolor=green
                         /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                         width 5 offset 0
         misc_feature    23..33
                         /label=KpnI-3
                         /ApEinfo_fwdcolor=cyan
                         /ApEinfo_revcolor=green
                         /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                         width 5 offset 0
         misc_feature    1..13
                         /label=NlaIV-1
                         /ApEinfo_fwdcolor=cyan
                         /ApEinfo_revcolor=green
                         /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                         width 5 offset 0
         misc_feature    14..20
                         /label=NlaIV-2
                         /ApEinfo_fwdcolor=cyan
                         /ApEinfo_revcolor=green
                         /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                         width 5 offset 0
         misc_feature    21..33
                         /label=NlaIV-3
                         /ApEinfo_fwdcolor=cyan
                         /ApEinfo_revcolor=green
                         /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                         width 5 offset 0
    ORIGIN
            1 GAATTCacag ggtaccaGGT ACCagcgGAT ATC
    //

    """
    )

    assert a.seguid() == "ldseguid=O1VPUWEeJ093UfmqTE5ObrI_2Yw"

    assert [x.qualifiers["label"][0] for x in a.features] == [
        "Acc65I-1",
        "Acc65I-2",
        "Acc65I-3",
        "KpnI-1",
        "KpnI-2",
        "KpnI-3",
        "NlaIV-1",
        "NlaIV-2",
        "NlaIV-3",
    ]

    b, c, d = a.cut(Acc65I)

    assert [x.qualifiers["label"][0] for x in b.features] == [
        "Acc65I-1",
        "KpnI-1",
        "NlaIV-1",
    ]
    assert [x.qualifiers["label"][0] for x in c.features] == [
        "Acc65I-2",
        "KpnI-2",
        "NlaIV-2",
    ]
    assert [x.qualifiers["label"][0] for x in d.features] == [
        "Acc65I-3",
        "KpnI-3",
        "NlaIV-3",
    ]
    e = b + c + d
    assert sorted([x.qualifiers["label"][0] for x in e.features]) == [x.qualifiers["label"][0] for x in a.features]
    assert str(a.seq) == str(e.seq)

    b, c, d = a.cut(KpnI)
    assert [x.qualifiers["label"][0] for x in b.features] == [
        "Acc65I-1",
        "KpnI-1",
        "NlaIV-1",
    ]
    assert [x.qualifiers["label"][0] for x in c.features] == [
        "Acc65I-2",
        "KpnI-2",
        "NlaIV-2",
    ]  # ???
    assert [x.qualifiers["label"][0] for x in d.features] == [
        "Acc65I-3",
        "KpnI-3",
        "NlaIV-3",
    ]
    e = b + c + d
    assert sorted([x.qualifiers["label"][0] for x in e.features]) == [x.qualifiers["label"][0] for x in a.features]

    b, c, d = a.cut(NlaIV)
    assert [x.qualifiers["label"][0] for x in b.features] == ["Acc65I-1", "NlaIV-1"]
    assert [x.qualifiers["label"][0] for x in c.features] == ["NlaIV-2"]
    assert [x.qualifiers["label"][0] for x in d.features] == ["KpnI-3", "NlaIV-3"]
    e = b + c + d
    assert str(a.seq) == str(e.seq)

    b, c = a.cut(EcoRI)
    e = b + c
    assert str(a.seq) == str(e.seq)

    b, c = a.cut(EcoRV)
    e = b + c
    assert str(a.seq) == str(e.seq)

    b, c, d = a.cut(EcoRI, EcoRV)
    e = b + c + d

    assert str(a.seq) == str(e.seq)

    b, c, d, f = a.cut(Acc65I, EcoRI)
    e = b + c + d + f
    assert str(a.seq) == str(e.seq)

    b, c, d, f = a.cut(EcoRI, Acc65I)
    e = b + c + d + f
    assert str(a.seq) == str(e.seq)


def test_features_on_slice():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqFeature import SimpleLocation

    dseq_record = Dseqrecord(Dseq("ACTCTTTCTCTCTCT", circular=True))
    dseq_record.features = [SeqFeature(SimpleLocation(start=2, end=4))]
    assert len(dseq_record[6:1].features) == 0
    assert len(dseq_record[6:3].features) == 0
    assert len(dseq_record[6:4].features) == 0  # ??? This slice should be empty?
    assert len(dseq_record[6:5].features) == 0  # ??? This slice should be empty?
    assert len(dseq_record[:].features) == 1

    dseq_record2 = Dseqrecord(Dseq("ACTCTTTCTCTCTCT"))
    assert dseq_record2[6:1].seq == Dseq("")


def test_features_change_ori():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    # Shifted a sequence by zero returns a copy
    s = Dseqrecord("GGATCC", circular=True)
    assert s.shifted(0) == s
    assert s.shifted(0) is not s

    s1 = read(
        """
        LOCUS       New_DNA                   13 bp ds-DNA     circular     08-JAN-2018
        DEFINITION  .
        SOURCE      .
        COMMENT     feature copied in ApE version 3.1.2 yields: tcccgtttt
        COMMENT     ApEinfo:methylated:1
        FEATURES             Location/Qualifiers
             misc_feature    join(2..3,5..7,9..12)
                             /locus_tag="hej"
                             /label="hej"
                             /ApEinfo_label="hej"
                             /ApEinfo_fwdcolor="cyan"
                             /ApEinfo_revcolor="green"
                             /ApEinfo_graphicformat="arrow_data {{0 1 2 0 0 -1} {} 0}
                             width 5 offset 0"
        ORIGIN
                1 atcaccgatt tta
        //"""
    )

    for i in range(1, len(s1)):
        b = s1.shifted(i)
        assert str(b.features[0].extract(b).seq).lower() == "tcccgtttt"

    s2 = read(
        """
        LOCUS       New_DNA                   13 bp ds-DNA     circular     08-JAN-2018
        DEFINITION  .
        SOURCE      .
        COMMENT
        COMMENT     feature copied in ApE version 3.1.2 yields: aaaacggga
        COMMENT     should be tcccgtttt, since feature is on the other strand
        COMMENT     ApEinfo:methylated:1
        FEATURES             Location/Qualifiers
             misc_feature    complement(join(2..5,7..9,11..12))
                             /locus_tag="hej"
                             /label="hej"
                             /ApEinfo_label="hej"
                             /ApEinfo_fwdcolor="cyan"
                             /ApEinfo_revcolor="green"
                             /ApEinfo_graphicformat="arrow_data {{0 1 2 0 0 -1} {} 0}
                             width 5 offset 0"
        ORIGIN
                1 taaaatcggt gat
        //

        """
    )

    assert str(s2.features[0].extract(s2).seq).lower() == "tcccgtttt"

    for i in range(0, len(s2)):
        b = s2.shifted(i)
        assert str(b.features[0].extract(b).seq).lower() == "tcccgtttt"

    # for x in b.features[0].location.parts:
    #     print(x.extract(b.seq))

    # from pydna.utils import shift_location
    # from Bio.SeqFeature import SimpleLocation
    # from Bio.SeqFeature import CompoundLocation
    # from Bio.SeqFeature import ExactPosition

    # ml6 = s2.shifted(6).features[0].location
    # ml7 = s2.shifted(7).features[0].location

    # ccc = shift_location(s2.features[0].location, -6, 13)
    # ddd = shift_location(s2.features[0].location, -7, 13)

    # print(b.seq, b.features[0].extract(b).seq, b.features[0].location)

    s3 = read(
        """
            LOCUS       New_DNA                   21 bp ds-DNA     circular     03-APR-2013
            DEFINITION  a
            ACCESSION
            VERSION
            SOURCE      .
              ORGANISM  .
            COMMENT
            COMMENT     ApEinfo:methylated:1
            FEATURES             Location/Qualifiers
                 misc_feature    join(18..21,1..4)
                                 /label=bb
                                 /ApEinfo_fwdcolor=cyan
                                 /ApEinfo_revcolor=green
                                 /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                 width 5 offset 0
                 misc_feature    5..17
                                 /label=ins
                                 /ApEinfo_fwdcolor=#e03c2b
                                 /ApEinfo_revcolor=green
                                 /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                 width 5 offset 0
            ORIGIN
                    1 aaagGTACCT TTGGATCcggg
            //
        """
    )

    bbfeat = Dseq.from_representation("cgggaaag\n" "gccctttc")
    insfeat = Dseq.from_representation("GTACCTTTGGATC\n" "CATGGAAACCTAG")

    assert str(s3.features[0].extract(s3).seq) == str(bbfeat).upper()  # bb
    assert str(s3.features[1].extract(s3).seq) == str(insfeat).upper()  # ins
    assert s3.features[0].extract(s3).seq == bbfeat.upper()  # bb
    assert s3.features[1].extract(s3).seq == insfeat  # ins

    for i in range(1, len(s3)):
        b = s3.shifted(i)

        assert [str(f.extract(b).seq) for f in b.features if f.qualifiers["label"][0] == "ins"][0] == "GTACCTTTGGATC"
        assert [str(f.extract(b).seq) for f in b.features if f.qualifiers["label"][0] == "bb"][0] == "CGGGAAAG"

    from Bio.Restriction import Acc65I, BamHI

    inseq = Dseq.from_representation("GTACCTTTG\n" "    GAAACCTAG")

    bbseq = Dseq.from_representation("GATCCGGGAAAG\n" "    GCCCTTTCCATG")

    assert s3.seq.cut(Acc65I, BamHI) == (inseq, bbseq)

    bb1, ins1 = sorted(s3.cut(Acc65I, BamHI), key=len, reverse=True)

    assert str(bbfeat).upper() in bbseq
    assert str(insfeat).upper() in inseq

    for i in range(0, len(s3)):
        b = s3.shifted(i)
        feature1 = b.features[0].extract(b).seq
        feature2 = b.features[1].extract(b).seq
        assert sorted([feature1, feature2]) == [s3.features[0].extract(s3).seq, s3.features[1].extract(s3).seq]

        """        11
                   |
        ATCcgggaaagGTACCTTTGG
        taggccctttccatggaaacc

                   GTACCTTTGGATCCGGGAAAG
                       GAAACCTAGGCCCTTTCCATG

                   GTACCTTTG
                       GAAACCTAG

                            GATCCGGGAAAG
                                GCCCTTTCCATG


        """

        bb, ins = sorted(b.cut(Acc65I, BamHI), key=len, reverse=True)

        assert eq(bb1, bb)
        assert eq(ins1, ins)

        assert bb.features[0].extract(bb).seq == bbfeat.upper()

        assert ins.features[0].extract(ins).seq.full_sequence == insfeat


def test_amijalis():
    # Thanks to https://github.com/amijalis
    from pydna.dseqrecord import Dseqrecord

    test_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCG"

    # length of test_seq is 36.

    test_seq_dseqrecord = Dseqrecord(test_seq).looped()
    test_seq_dseqrecord.add_feature(0, 36, type_="test", label="test")

    # print(f'Features before shifting: {test_seq_dseqrecord.features}')

    f1 = test_seq_dseqrecord.features[0]

    test_seq_shifted = test_seq_dseqrecord.shifted(10)

    f2 = test_seq_shifted.features[0]

    assert f1.extract(test_seq_dseqrecord).seq == f2.extract(test_seq_shifted).seq

    # print(f'Features after shifting: {test_seq_shifted.features}')


def test_figure():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from Bio.Restriction import Acc65I, KpnI, ApaI, Bsp120I
    from Bio.SeqFeature import SimpleLocation
    from Bio.SeqFeature import SeqFeature

    # broken feature linear

    linearDseq = Dseq.from_representation(
        """
    gatcggtaccgatcATGAAATAAgatcGGGCCCgatc
    ctagccatggctagTACTTTATTctagCCCGGGctag
    """
    )

    linearDseqrecord = Dseqrecord(linearDseq)

    assert (
        linearDseqrecord.figure()
        == "Dseqrecord(-37)\n\x1b[48;5;11m\x1b[0mgatcggtaccgatcATGAAATAAgatcGGGCCCgatc\nctagccatggctagTACTTTATTctagCCCGGGctag"
    )

    linearDseqrecord.features.append(SeqFeature(SimpleLocation(14, 17, 1) + SimpleLocation(20, 23, 1), type="test"))

    expect = "Dseqrecord(-37)\ngatcggtaccgatc\x1b[48;5;11mATG\x1b[0mAAA\x1b[48;5;11mTAA\x1b[0mgatcGGGCCCgatc\nctagccatggctagTACTTTATTctagCCCGGGctag"
    assert linearDseqrecord.figure() == expect

    # short feature linear

    linearDseqrecord = Dseqrecord(linearDseq)
    linearDseqrecord.add_feature(14, 23)

    assert (
        linearDseqrecord.figure()
        == "Dseqrecord(-37)\ngatcggtaccgatc\x1b[48;5;11mATGAAATAA\x1b[0mgatcGGGCCCgatc\nctagccatggctagTACTTTATTctagCCCGGGctag"
    )

    feat = Dseq("ATGAAATAA")
    assert linearDseqrecord.features[0].extract(linearDseqrecord).seq == feat

    a1, b1, c1 = linearDseqrecord.cut(Acc65I, Bsp120I)
    assert (
        b1.figure() == "Dseqrecord(-27)\ngtaccgatc\x1b[48;5;11mATGAAATAA\x1b[0mgatcG    \n    gctagTACTTTATTctagCCCGG"
    )
    assert b1.extract_feature(0).seq == feat
    assert b1.rc().extract_feature(0).seq == feat
    assert (
        b1.rc().figure()
        == "Dseqrecord(-27)\nGGCCCgatcTTATTTCATgatcg    \n    Gctag\x1b[48;5;11mAATAAAGTA\x1b[0mctagccatg"
    )

    a2, b2, c2 = linearDseqrecord.cut(KpnI, Bsp120I)
    assert (
        b2.figure() == "Dseqrecord(-27)\n    cgatc\x1b[48;5;11mATGAAATAA\x1b[0mgatcG    \ncatggctagTACTTTATTctagCCCGG"
    )
    assert b2.extract_feature(0).seq == feat
    assert b2.rc().extract_feature(0).seq == feat
    assert (
        b2.rc().figure()
        == "Dseqrecord(-27)\nGGCCCgatcTTATTTCATgatcggtac\n    Gctag\x1b[48;5;11mAATAAAGTA\x1b[0mctagc    "
    )

    a3, b3, c3 = linearDseqrecord.cut(Acc65I, ApaI)
    assert (
        b3.figure() == "Dseqrecord(-27)\ngtaccgatc\x1b[48;5;11mATGAAATAA\x1b[0mgatcGGGCC\n    gctagTACTTTATTctagC    "
    )
    assert b3.extract_feature(0).seq == feat
    assert b3.rc().extract_feature(0).seq == feat
    assert (
        b3.rc().figure()
        == "Dseqrecord(-27)\n    CgatcTTATTTCATgatcg    \nCCGGGctag\x1b[48;5;11mAATAAAGTA\x1b[0mctagccatg"
    )

    a4, b4, c4 = linearDseqrecord.cut(KpnI, ApaI)
    assert (
        b4.figure() == "Dseqrecord(-27)\n    cgatc\x1b[48;5;11mATGAAATAA\x1b[0mgatcGGGCC\ncatggctagTACTTTATTctagC    "
    )
    assert b4.extract_feature(0).seq == feat
    assert b4.rc().extract_feature(0).seq == feat
    assert (
        b4.rc().figure()
        == "Dseqrecord(-27)\n    CgatcTTATTTCATgatcggtac\nCCGGGctag\x1b[48;5;11mAATAAAGTA\x1b[0mctagc    "
    )

    # short feature circular

    circularDseqrecord = linearDseqrecord.looped()
    assert (
        circularDseqrecord.figure()
        == "Dseqrecord(o37)\ngatcggtaccgatc\x1b[48;5;11mATGAAATAA\x1b[0mgatcGGGCCCgatc\nctagccatggctagTACTTTATTctagCCCGGGctag"
    )
    assert circularDseqrecord.features[0].extract(circularDseqrecord).seq == feat

    a5, b5 = circularDseqrecord.cut(Acc65I, Bsp120I)
    assert (
        a5.figure() == "Dseqrecord(-27)\ngtaccgatc\x1b[48;5;11mATGAAATAA\x1b[0mgatcG    \n    gctagTACTTTATTctagCCCGG"
    )
    assert a5.extract_feature(0).seq == feat

    a6, b6 = circularDseqrecord.cut(KpnI, Bsp120I)
    assert (
        a6.figure() == "Dseqrecord(-27)\n    cgatc\x1b[48;5;11mATGAAATAA\x1b[0mgatcG    \ncatggctagTACTTTATTctagCCCGG"
    )
    assert a6.extract_feature(0).seq == feat

    a7, b7 = circularDseqrecord.cut(Acc65I, ApaI)
    assert (
        a7.figure() == "Dseqrecord(-27)\ngtaccgatc\x1b[48;5;11mATGAAATAA\x1b[0mgatcGGGCC\n    gctagTACTTTATTctagC    "
    )
    assert a7.extract_feature(0).seq == feat

    a8, b8 = circularDseqrecord.cut(KpnI, ApaI)
    assert (
        a8.figure() == "Dseqrecord(-27)\n    cgatc\x1b[48;5;11mATGAAATAA\x1b[0mgatcGGGCC\ncatggctagTACTTTATTctagC    "
    )
    assert a8.extract_feature(0).seq == feat

    a9, b9 = circularDseqrecord.cut(Bsp120I, KpnI)
    assert (
        a9.figure() == "Dseqrecord(-27)\n    cgatc\x1b[48;5;11mATGAAATAA\x1b[0mgatcG    \ncatggctagTACTTTATTctagCCCGG"
    )
    assert a9.extract_feature(0).seq == feat

    # longer feature linear

    linearDseqrecord = Dseqrecord(linearDseq)
    linearDseqrecord.add_feature(5, 32)
    assert (
        linearDseqrecord.figure()
        == "Dseqrecord(-37)\ngatcg\x1b[48;5;11mgtaccgatcATGAAATAAgatcGGGCC\x1b[0mCgatc\nctagccatggctagTACTTTATTctagCCCGGGctag"
    )
    feat = Dseq.from_representation(
        """
    gtaccgatcATGAAATAAgatcGGGCC
    catggctagTACTTTATTctagCCCGG
    """
    )
    assert linearDseqrecord.features[0].extract(linearDseqrecord).seq == feat

    a10, b10, c10 = linearDseqrecord.cut(Acc65I, Bsp120I)
    assert (
        b10.figure() == "Dseqrecord(-27)\n\x1b[48;5;11mgtaccgatcATGAAATAAgatcG    \x1b[0m\n    gctagTACTTTATTctagCCCGG"
    )
    feat10 = Dseq.from_representation(
        """
    gtaccgatcATGAAATAAgatcG
        gctagTACTTTATTctagCCCGG
    """
    )
    assert b10.extract_feature(0).seq == feat10
    assert b10.rc().extract_feature(0).seq == feat10
    assert (
        b10.rc().figure()
        == "Dseqrecord(-27)\nGGCCCgatcTTATTTCATgatcg    \n\x1b[48;5;11m    GctagAATAAAGTActagccatg\x1b[0m"
    )

    a11, b11, c11 = linearDseqrecord.cut(KpnI, Bsp120I)
    assert (
        b11.figure() == "Dseqrecord(-27)\n\x1b[48;5;11m    cgatcATGAAATAAgatcG    \x1b[0m\ncatggctagTACTTTATTctagCCCGG"
    )
    feat11 = Dseq.from_representation(
        """
        cgatcATGAAATAAgatcG
    catggctagTACTTTATTctagCCCGG
    """
    )
    assert b11.extract_feature(0).seq == feat11
    assert b11.rc().extract_feature(0).seq == feat11
    assert (
        b11.rc().figure()
        == "Dseqrecord(-27)\nGGCCCgatcTTATTTCATgatcggtac\n\x1b[48;5;11m    GctagAATAAAGTActagc    \x1b[0m"
    )

    a12, b12, c12 = linearDseqrecord.cut(Acc65I, ApaI)
    assert (
        b12.figure() == "Dseqrecord(-27)\n\x1b[48;5;11mgtaccgatcATGAAATAAgatcGGGCC\x1b[0m\n    gctagTACTTTATTctagC    "
    )
    feat12 = Dseq.from_representation(
        """
    gtaccgatcATGAAATAAgatcGGGCC
        gctagTACTTTATTctagC
    """
    )
    assert b12.extract_feature(0).seq == feat12
    assert b12.rc().extract_feature(0).seq == feat12
    assert (
        b12.rc().figure()
        == "Dseqrecord(-27)\n    CgatcTTATTTCATgatcg    \n\x1b[48;5;11mCCGGGctagAATAAAGTActagccatg\x1b[0m"
    )

    a13, b13, c13 = linearDseqrecord.cut(KpnI, ApaI)
    assert (
        b13.figure() == "Dseqrecord(-27)\n\x1b[48;5;11m    cgatcATGAAATAAgatcGGGCC\x1b[0m\ncatggctagTACTTTATTctagC    "
    )
    feat13 = Dseq.from_representation(
        """
        cgatcATGAAATAAgatcGGGCC
    catggctagTACTTTATTctagC
    """
    )
    assert b13.extract_feature(0).seq == feat13
    assert b13.rc().extract_feature(0).seq == feat13
    assert (
        b13.rc().figure()
        == "Dseqrecord(-27)\n    CgatcTTATTTCATgatcggtac\n\x1b[48;5;11mCCGGGctagAATAAAGTActagc    \x1b[0m"
    )

    # longer feature circular

    circularDseqrecord = linearDseqrecord.looped()
    circularDseqrecord.figure()
    a14, b14 = circularDseqrecord.cut(KpnI, Bsp120I)
    assert (
        a14.figure() == "Dseqrecord(-27)\n\x1b[48;5;11m    cgatcATGAAATAAgatcG    \x1b[0m\ncatggctagTACTTTATTctagCCCGG"
    )
    feat14 = Dseq.from_representation(
        """
        cgatcATGAAATAAgatcG
    catggctagTACTTTATTctagCCCGG
    """
    )
    assert a14.extract_feature(0).seq == feat14

    a15, b15 = circularDseqrecord.cut(KpnI, ApaI)
    assert (
        a15.figure() == "Dseqrecord(-27)\n\x1b[48;5;11m    cgatcATGAAATAAgatcGGGCC\x1b[0m\ncatggctagTACTTTATTctagC    "
    )
    feat15 = Dseq.from_representation(
        """
        cgatcATGAAATAAgatcGGGCC
    catggctagTACTTTATTctagC
    """
    )
    assert a15.extract_feature(0).seq == feat15

    a16, b16 = circularDseqrecord.cut(Acc65I, Bsp120I)
    assert (
        a16.figure() == "Dseqrecord(-27)\n\x1b[48;5;11mgtaccgatcATGAAATAAgatcG    \x1b[0m\n    gctagTACTTTATTctagCCCGG"
    )
    feat16 = Dseq.from_representation(
        """
    gtaccgatcATGAAATAAgatcG
        gctagTACTTTATTctagCCCGG
    """
    )
    assert a16.extract_feature(0).seq == feat16

    a17, b17 = circularDseqrecord.cut(Acc65I, ApaI)
    assert (
        a17.figure() == "Dseqrecord(-27)\n\x1b[48;5;11mgtaccgatcATGAAATAAgatcGGGCC\x1b[0m\n    gctagTACTTTATTctagC    "
    )
    feat17 = Dseq.from_representation(
        """
    gtaccgatcATGAAATAAgatcGGGCC
        gctagTACTTTATTctagC
    """
    )
    assert a17.extract_feature(0).seq == feat17

    # Wrap around feature circular on watson

    circularDseqrecord = Dseqrecord(linearDseq, circular=True)
    circularDseqrecord.add_feature(32, 5)
    assert (
        circularDseqrecord.figure()
        == "Dseqrecord(o37)\n\x1b[48;5;11mgatcg\x1b[0mgtaccgatcATGAAATAAgatcGGGCC\x1b[48;5;11mCgatc\x1b[0m\nctagccatggctagTACTTTATTctagCCCGGGctag"
    )

    feat = Dseq.from_representation(
        """
    Cgatcgatcg
    Gctagctagc
    """
    )
    assert circularDseqrecord.extract_feature(0).seq == feat

    a18, b18 = circularDseqrecord.cut(KpnI, Bsp120I)
    assert b18.figure() == "Dseqrecord(-18)\nGGCC\x1b[48;5;11mCgatcgatcg\x1b[0mgtac\n    Gctagctagc    "
    assert b18.extract_feature(0).seq == feat

    a19, b19 = circularDseqrecord.cut(KpnI, ApaI)
    assert b19.figure() == "Dseqrecord(-18)\n    \x1b[48;5;11mCgatcgatcg\x1b[0mgtac\nCCGGGctagctagc    "
    assert b19.extract_feature(0).seq == feat

    a20, b20 = circularDseqrecord.cut(Acc65I, Bsp120I)
    assert b20.figure() == "Dseqrecord(-18)\nGGCC\x1b[48;5;11mCgatcgatcg\x1b[0m    \n    Gctagctagccatg"
    assert b20.extract_feature(0).seq == feat

    a21, b21 = circularDseqrecord.cut(Acc65I, ApaI)
    assert b21.figure() == "Dseqrecord(-18)\n    \x1b[48;5;11mCgatcgatcg\x1b[0m    \nCCGGGctagctagccatg"
    assert b21.extract_feature(0).seq == feat

    # Wrap around feature circular on crick

    circularDseqrecord = Dseqrecord(linearDseq, circular=True)
    circularDseqrecord.add_feature(32, 5)
    circularDseqrecord = circularDseqrecord.rc()
    assert (
        circularDseqrecord.figure()
        == "Dseqrecord(o37)\ngatcGGGCCCgatcTTATTTCATgatcggtaccgatc\n\x1b[48;5;11mctagC\x1b[0mCCGGGctagAATAAAGTActagccatg\x1b[48;5;11mgctag\x1b[0m"
    )

    feat = Dseq.from_representation(
        """
    Cgatcgatcg
    Gctagctagc
    """
    )
    assert circularDseqrecord.extract_feature(0).seq == feat

    a22, b22 = circularDseqrecord.cut(KpnI, Bsp120I)

    # Passes the tests if changed to "Dseqrecord(-18)\n    cgatcgatcG    \ncatg\x1b[48;5;11mgctagctagC\x1b[0mCCGG"
    assert b22.figure() == "Dseqrecord(-18)\n    cgatcgatcG    \ncatg\x1b[48;5;11mgctagctagC\x1b[0mCCGG"
    assert b22.extract_feature(0).seq == feat

    a23, b23 = circularDseqrecord.cut(KpnI, ApaI)
    assert b23.figure() == "Dseqrecord(-18)\n    cgatcgatcGGGCC\ncatg\x1b[48;5;11mgctagctagC\x1b[0m    "
    assert b23.extract_feature(0).seq == feat

    a24, b24 = circularDseqrecord.cut(Acc65I, Bsp120I)
    assert b24.figure() == "Dseqrecord(-18)\ngtaccgatcgatcG    \n    \x1b[48;5;11mgctagctagC\x1b[0mCCGG"
    assert b24.extract_feature(0).seq == feat

    a25, b25 = circularDseqrecord.cut(Acc65I, ApaI)
    assert b25.figure() == "Dseqrecord(-18)\ngtaccgatcgatcGGGCC\n    \x1b[48;5;11mgctagctagC\x1b[0m    "
    assert b25.extract_feature(0).seq == feat


# @pytest.mark.xfail(reason="issue #78")
def test_jan_glx():
    # Thanks to https://github.com/jan-glx
    from Bio.Restriction import NdeI, BamHI
    from pydna.readers import read

    # from pydna.genbank import Genbank
    # gb = Genbank("bjornjobb@gmail.com")
    # puc19 = gb.nucleotide("M77789.2")
    # assert puc19.seguid() == "n-NZfWfjHgA7wKoEBU6zfoXib_0"
    # puc19.write("pUC19_M77789.gb")
    puc19 = read("pUC19_M77789.gb")
    assert puc19.seguid() == "cdseguid=mCC0B3UMZfgLyh3Pl574MVjm30U"
    insert, bb = puc19.cut(NdeI, BamHI)  # Note the order !

    puc19_ = (bb + insert).looped().synced(puc19)
    assert puc19_.seguid() == "cdseguid=mCC0B3UMZfgLyh3Pl574MVjm30U"

    # Some features are lost because they spanned the cutting sites in puc19.
    assert puc19_.extract_feature(0).seq == puc19.extract_feature(2).seq
    assert puc19_.extract_feature(1).seq == puc19.extract_feature(4).seq
    assert puc19_.extract_feature(2).seq == puc19.extract_feature(6).seq
    assert puc19_.extract_feature(3).seq == puc19.extract_feature(7).seq


def test_synced():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read

    bb_ins = Dseqrecord("tcgcgcgtttcgATGTAgtgatgacggtgaA", circular=True)

    for i in range(1, len(bb_ins)):
        cand = bb_ins.shifted(i)
        assert str(cand.synced("tcgcgcgtttcg", limit=11).seq).upper() == str(bb_ins.seq.upper())
        assert str(cand.synced(Dseq("tcgcgcgtttcg"), limit=11).seq).upper() == str(bb_ins.seq.upper())
        assert str(cand.synced(Seq("tcgcgcgtttcg"), limit=11).seq).upper() == str(bb_ins.seq.upper())
        assert str(cand.synced(SeqRecord(Seq("tcgcgcgtttcg")), limit=11).seq).upper() == str(bb_ins.seq.upper())

    pGUP1 = read("pGUP1_correct.gb")
    pGREG505 = read("pGREG505.gb")
    pGUP1_not_synced = read("pGUP1_not_synced.gb")
    assert pGUP1_not_synced.synced(pGREG505).seguid() == "cdseguid=QiK2pH9yioTPfSobUTLz4CPiNzY" == pGUP1.seguid()

    bb_ins = Dseqrecord("tcgcgcgtttcgAgtgatgacggtgaA", circular=True)

    for i in range(1, len(bb_ins)):
        cand = bb_ins.shifted(i)
        assert str(cand.synced("tcgcgcgtttcg", limit=11).seq).upper() == str(bb_ins.seq.upper())

    a = Dseqrecord("gtatcgttctagctattggctagta", circular=True)
    b = Dseqrecord("ggttagtcagttatatcggcttatc", circular=True)

    with pytest.raises(TypeError):
        a.synced(b)

    a = Dseqrecord("gtatcgttctagctattggctagta", circular=False)
    with pytest.raises(TypeError):
        a.synced(b)

    a = Dseqrecord("gtatcgttctagctattggctagta", circular=True)
    b = Dseqrecord("acgatactactagccaatagctaga", circular=True)

    assert str(a.synced(b).seq) == "acgatactactagccaatagctaga"
    assert str(a.synced(a).seq) == "gtatcgttctagctattggctagta"
    assert str(a.rc().synced(a).seq) == "gtatcgttctagctattggctagta"


def test_map_pCR_MCT1_HA46():
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read

    pCR_MCT1_HA46 = read("pCR_MCT1_HA46.gb")

    slc = pCR_MCT1_HA46.find_aa("VFFKE YPYDVPDYA IEG".replace(" ", ""))

    pCR_MCT1_HA46.map_target = slc

    with pytest.raises(ValueError):
        pCR_MCT1_HA46.map_trace_files("/subfolder/*.ab1")

    map_ = pCR_MCT1_HA46.map_trace_files("*.ab1")

    assert set(map_) == set(["28-1rev_D04_026.ab1", "32-3rev_H04_018.ab1", "36-5rev_D05_041.ab1"])

    assert set([x.fname for x in pCR_MCT1_HA46.matching_reads]) == set(
        ["28-1rev_D04_026.ab1", "32-3rev_H04_018.ab1", "36-5rev_D05_041.ab1"]
    )

    assert set([x.fname for x in pCR_MCT1_HA46.not_matching_reads]) == set(["02-G1_B01_013.ab1"])

    assert pCR_MCT1_HA46.find_aa("YPYDVPDYA".replace(" ", "")) == slice(1088, 1115, None)

    assert pCR_MCT1_HA46.find_aa("VFFKE YPYDVPDYA IEG".replace(" ", "")) == slice(1073, 1124, None)

    pCR_MCT1_HA46.add_feature(1073, 1124)
    mt = pCR_MCT1_HA46.features[-1]
    pCR_MCT1_HA46.map_target = mt

    map_ = pCR_MCT1_HA46.map_trace_files("*.ab1")
    assert set(map_) == set(["28-1rev_D04_026.ab1", "32-3rev_H04_018.ab1", "36-5rev_D05_041.ab1"])

    pCR_MCT1_HA46.map_target = None

    map_ = pCR_MCT1_HA46.map_trace_files("*.ab1")

    assert set(map_) == set(["32-3rev_H04_018.ab1", "36-5rev_D05_041.ab1", "28-1rev_D04_026.ab1"])


def test_map_short():
    from pydna.dseqrecord import Dseqrecord

    t = Dseqrecord("AAGTTAAAATAAGGCTAGTCCGTTAT")
    t.map_target = slice(0, 26)
    t.map_trace_files("*.ab1")
    assert t.map_trace_files("*.ab1") == ["02-G1_B01_013.ab1"]


def test_map_too_short():
    from pydna.dseqrecord import Dseqrecord

    t = Dseqrecord("AAGTTAAAATAAGGCTAGTCCGTT")
    # t.map_target = slice(0,23)
    t.map_trace_files("*.ab1")
    assert t.map_trace_files("*.ab1") == []


def test_map_no_match():
    from pydna.dseqrecord import Dseqrecord

    t = Dseqrecord(
        "AGGAGGGGCCCACACCGAGGAAGTAGACTGTTATACGTCGGCGATGGTGGTAGCTAACTATGTTGCCTGCCACTACAACAGTATCTAAGCCGTGTAAAGG"
    )  # random DNA!
    t.map_trace_files("*.ab1")
    assert t.map_trace_files("*.ab1") == []


def test_slicing2():
    from pydna.dseqrecord import Dseqrecord
    from Bio.Restriction import BamHI, EcoRI, KpnI

    a = Dseqrecord("aaGGATCCnnnnnnnnnGAATTCccc", circular=True)
    assert (
        a.cut(
            EcoRI,
            BamHI,
            KpnI,
        )[0].seq
        == a.cut(
            BamHI,
            EcoRI,
            KpnI,
        )[0].seq
    )


def test_rogerstager():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.utils import eq
    from Bio.Seq import Seq
    from Bio.Restriction import BsaI

    answ = []
    answ.append(Dseq("eeeeaaaaaaaaggtctcaffff"))
    answ.append(Dseq("eeeeaaaaaggtctcaffff"))

    tests = [Seq("aaaaaaggtctcaaaaaaa"), Seq("aaaaaaggtctcaaaa")]

    for s in tests:
        d = Dseqrecord(s).looped()
        for f in d.cut(BsaI):
            a = answ.pop(0)
            assert f.seq.watson == a.watson
            assert f.seq.crick == a.crick
            assert f.seq.ovhg == a.ovhg
            assert eq(f.seq, a)


def test___mul__():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=False)
    assert s * 3 == Dseqrecord("GGATCCGGATCCGGATCC", circular=False)
    assert s * 0 == Dseqrecord("")
    s = Dseqrecord("GGATCC", circular=True)
    with pytest.raises(TypeError):
        s * 3
    with pytest.raises(TypeError):
        s * 3.1


def test___repr__():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=False)
    assert repr(s) == "Dseqrecord(-6)"
    s = Dseqrecord("GGATCC", circular=True)
    assert repr(s) == "Dseqrecord(o6)"


def test__repr_pretty_():
    from unittest.mock import MagicMock
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=False)
    pp = MagicMock()
    s._repr_pretty_(pp, None)
    pp.text.assert_called_with("Dseqrecord(-6)")
    s = Dseqrecord("GGATCC", circular=True)
    s._repr_pretty_(pp, None)
    pp.text.assert_called_with("Dseqrecord(o6)")


def test___getitem__():
    from pydna.dseqrecord import Dseqrecord
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    s = Dseqrecord("GGATCC", circular=False)
    assert s[1:-1].seq == Dseqrecord("GATC", circular=False).seq
    t = Dseqrecord("GGATCC", circular=True)
    assert t[1:5].seq == Dseqrecord("GATC").seq
    #    assert t[1:5].__dict__ == Dseqrecord("GATC").__dict__
    #    assert str(t[1:5].__dict__) == str(Dseqrecord("GATC").__dict__)
    assert s[1:5].seq == Dseqrecord("GATC").seq
    assert s[1:5].seq == Dseqrecord("GATC", circular=False).seq
    assert s[5:1:-1].seq == Dseqrecord("CCTA", circular=False).seq

    assert t[1:1].seq == Dseqrecord("GATCCG").seq
    assert t[5:1].seq == Dseqrecord("").seq
    assert t[9:1].seq == Dseqrecord("").seq  # ??? How should slicing work?
    assert t[1:9].seq == Dseqrecord("GATCC").seq
    assert t[9:10].seq == Dseqrecord("").seq
    assert t[10:9].seq == Dseqrecord("").seq

    # Test how slicing works with features (using sequence as in test_features_change_ori)
    seqRecord = Dseqrecord("aaagGTACCTTTGGATCcggg", circular=True)
    f1 = SeqFeature(SimpleLocation(4, 17, 1), type="misc_feature")
    f2 = SeqFeature(SimpleLocation(17, 21, 1) + SimpleLocation(0, 4, 1), type="misc_feature")
    seqRecord.features = [f1, f2]

    # Exact feature sliced for normal and origin-spanning features
    assert len(seqRecord[4:17].features) == 1
    assert len(seqRecord[17:4].features) == 0  # ??? how sould slicing work?

    # Partial feature sliced for normal and origin-spanning features
    assert len(seqRecord[2:20].features) == 1
    assert len(seqRecord[13:8].features) == 0  # ??? how sould slicing work?

    # Indexing of full circular molecule (https://github.com/BjornFJohansson/pydna/issues/161)
    # s = Dseqrecord("GGATCC", circular=True) # TODO: discuss this
    # str_seq = str(s.seq)
    # for shift in range(len(s)):
    #     assert str(s[shift:shift].seq) == str_seq[shift:] + str_seq[:shift]


def test___eq__():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=False)
    t = Dseqrecord("GGATCC", circular=False)
    u = Dseqrecord("GGATCC", circular=True)
    v = Dseqrecord("GGATCC", circular=True)

    assert s == t
    assert u == v
    assert s != u

    assert s != 123


def test___ne__():
    pass


def test___hash__():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=False)
    t = Dseqrecord("GGATCC", circular=False)
    u = Dseqrecord("GGATCC", circular=True)
    v = Dseqrecord("GGATCC", circular=True)

    assert hash(s) == hash(t)
    assert hash(u) == hash(v)
    assert hash(s) != hash(u)


def test_linearize():
    from Bio.Restriction import BamHI, BglII
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord

    u = Dseqrecord("GGATCC", circular=True)
    frag = u.linearize(BamHI)
    assert frag.seq == Dseq("PEXICGQFZJ")
    result = Dseqrecord(Dseq("PEXICGQFZJ"))
    assert frag.seq == result.seq

    t = Dseqrecord("GGATCC", circular=False)
    with pytest.raises(TypeError):
        t.linearize(BamHI)
    with pytest.raises(TypeError):
        u.linearize(BglII)
    u = Dseqrecord("GGATCCagatct", circular=True)
    u.linearize(BglII)
    u.linearize(BamHI)
    with pytest.raises(TypeError):
        u.linearize(BamHI, BglII)
    with pytest.raises(TypeError):
        u.linearize(BamHI + BglII)
    u = Dseqrecord("GGATCCGGATCC", circular=True)
    with pytest.raises(TypeError):
        u.linearize(BamHI)


def test_cutters():
    from pydna.dseqrecord import Dseqrecord

    obj = Dseqrecord("ggatcc")
    from Bio.Restriction import BglII, BamHI

    assert BglII in obj.no_cutters()
    assert BamHI not in obj.no_cutters()

    assert BamHI in obj.unique_cutters()

    assert BamHI in obj.once_cutters()

    assert BamHI in (obj + obj).twice_cutters()
    assert BamHI not in obj.twice_cutters()

    assert BamHI in obj.n_cutters(1)
    assert BamHI in obj.cutters()


def test_number_of_cuts():
    from Bio.Restriction import EcoRI, BamHI, BglII
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCCgaattc", circular=False)
    s.cut(BamHI)
    s.cut(EcoRI)
    s.cut(BglII)
    from Bio.Restriction import RestrictionBatch

    rb = RestrictionBatch((EcoRI, BamHI, BglII))
    assert s.number_of_cuts(BamHI) == 1
    assert s.number_of_cuts(EcoRI) == 1
    assert s.number_of_cuts(BglII) == 0
    assert s.number_of_cuts(rb) == 2
    assert s.number_of_cuts(EcoRI + BamHI) == 2


def test_reverse_complement():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=False)
    assert s.reverse_complement() == s.rc()
    assert s.rc().seq == s.seq
    assert s.reverse_complement().seq == s.seq

    from pydna.readers import read

    s = read(
        """
                LOCUS       New_DNA                   13 bp ds-DNA     circular     12-NOV-2013
                DEFINITION  .
                ACCESSION
                VERSION
                SOURCE      .
                ORGANISM  .
                COMMENT
                COMMENT     ApEinfo:methylated:1
                FEATURES             Location/Qualifiers
                     misc_feature    join(9..10,12..13,1..1,3..6)
                                     /label=hej
                                     /ApEinfo_fwdcolor=cyan
                                     /ApEinfo_revcolor=green
                                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                     width 5 offset 0
                ORIGIN
                        1 gattttaatc acc
                //"""
    )


def test_shifted():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCCgaattc", circular=False)
    with pytest.raises(TypeError):
        s.shifted(1)


def test_looped():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord

    from Bio import BiopythonDeprecationWarning
    import warnings

    warnings.simplefilter("always")

    a = Dseqrecord("aAAa")
    a.add_feature()
    b = a.looped()
    assert a.features[0].extract(a).seq == b.features[0].extract(b).seq
    assert str(a.features[0].extract(a).seq) == str(b.features[0].extract(b).seq)
    assert a.features == b.features

    a = Dseqrecord("aAAa")
    a.add_feature(0, 4)
    b = a.looped()
    assert a.features[0].extract(a).seq == b.features[0].extract(b).seq
    assert str(a.features[0].extract(a).seq) == str(b.features[0].extract(b).seq)
    assert a.features == b.features

    a = Dseqrecord("aAAa")
    a.add_feature(2, 4)
    b = a.looped()
    assert a.features[0].extract(a).seq == b.features[0].extract(b).seq
    assert str(a.features[0].extract(a).seq) == str(b.features[0].extract(b).seq)
    assert a.features == b.features

    a = Dseqrecord("aAAa")
    a.add_feature(2, 4)
    b = a.looped()
    assert a.features[0].extract(a).seq == b.features[0].extract(b).seq
    assert str(a.features[0].extract(a).seq) == str(b.features[0].extract(b).seq)
    assert a.features == b.features

    a = Dseqrecord("aAAa")
    a.add_feature(2, 4)
    b = a.looped()
    assert a.features[0].extract(a).seq == b.features[0].extract(b).seq
    assert str(a.features[0].extract(a).seq) == str(b.features[0].extract(b).seq)
    assert a.features == b.features

    a = Dseqrecord(Dseq("pAAaq"))
    a.add_feature(2, 4)
    b = a.looped()
    assert a.features[0].extract(a).seq == b.features[0].extract(b).seq
    assert str(a.features[0].extract(a).seq) == str(b.features[0].extract(b).seq)
    assert a.features == b.features

    a = Dseqrecord(Dseq("iaaaj"))
    a.add_feature(0, 5)
    b = a.looped()
    assert a.features == b.features

    a = Dseqrecord(Dseq("iaaaj"))
    a.add_feature(0, 5, strand=-1)
    b = a.looped()
    assert a.features == b.features

    a = Dseqrecord(Dseq("jaaai"))
    a.add_feature(2, 4)
    b = a.looped()
    assert a.features == b.features

    a = Dseqrecord(Dseq("faaae"))
    a.add_feature(0, 5)
    b = a.looped()
    assert a.features == b.features

    a = Dseqrecord(Dseq("faaae"))
    a.add_feature(0, 5, strand=-1)
    b = a.looped()
    assert a.features == b.features

    a = Dseqrecord(Dseq("eaaaf"))
    a.add_feature(0, 6)
    b = a.looped()
    assert a.features == b.features


def test_upper():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("Gc")
    s.annotations["sample"] = ["sample"]
    u = s.upper()
    assert s.upper() == u
    assert s.seq.upper() == u.seq.upper()
    del u.__dict__["_seq"]
    del s.__dict__["_seq"]
    assert u.__dict__ == s.__dict__


def test_lower():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("Gc")
    s.annotations["sample"] = ["sample"]
    l = s.lower()
    assert s.lower() == l
    assert s.seq.upper() == l.seq.upper()
    del l.__dict__["_seq"]
    del s.__dict__["_seq"]
    assert l.__dict__ == s.__dict__


def test_map():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec
    from Bio.SeqIO import read as abiread

    traces = []

    import glob

    for name in glob.glob("*.ab1"):
        traces.append(abiread(name, "abi"))

    for t in traces:
        d = Dseqrecord(t.seq)

        if "ITVFFKEYPYDVPDYAIEGIFHAT" in d:
            #                     Y   P   Y   D   V   P   D   Y   A
            s = Dseqrecord(Dseq("tat cca tat gac gtt cca gac tat gca".replace(" ", "")))
            sl = s.find_aa("YPYDVPDYA")
            assert str(s[sl].seq.translate()) == "YPYDVPDYA"
            assert "YPYDVPDYA" in s

            s = Dseqrecord(Dseq("EEE tat cca tat gac gtt cca gac tat gca".replace(" ", "")))
            sl = s.find_aa("YPYDVPDYA")
            assert str(s[sl].seq.translate()) == "YPYDVPDYA"
            assert "YPYDVPDYA" in s

            s = Dseqrecord(Dseq("FFF tat cca tat gac gtt cca gac tat gca".replace(" ", "")))
            sl = s.find_aa("YPYDVPDYA")
            assert str(s[sl].seq.translate()) == "YPYDVPDYA"
            assert "YPYDVPDYA" in s

            s = Dseqrecord(Dseq("FFF tat cca tat gac gtt cca gac tat gca".replace(" ", "")).rc())
            assert s.find_aa("YPYDVPDYA") == None

            assert str(s.rc()[sl].seq.translate(to_stop=False)) == "YPYDVPDYA"
            assert "YPYDVPDYA" in s.rc()

            s = Dseqrecord(Dseq("aaa tat cca tat gac gtt cca gac tat gca".replace(" ", ""), circular=True))
            sl = s.find_aa("YPYDVPDYA")
            assert str(s[sl].seq.translate()) == "YPYDVPDYA"
            assert "YPYDVPDYA" in s


def test_assemble_YEp24PGK_XK():
    """
    test YEp24PGK_XK
    """
    import pytest
    import sys

    from pydna.readers import read
    from pydna.utils import eq
    from pydna.amplify import pcr

    """ test YEp24PGK_XK"""
    p1 = read("primer1.txt", ds=False)
    p3 = read("primer3.txt", ds=False)
    XKS1 = read("XKS1_orf.txt")
    YEp24PGK = read("YEp24PGK.txt")

    PCR_prod = pcr(p1, p3, XKS1)

    from Bio.Restriction import BamHI

    stuffer1, insert, stuffer2 = PCR_prod.cut(BamHI)

    from Bio.Restriction import BglII

    YEp24PGK_BglII = YEp24PGK.cut(BglII)[0]

    YEp24PGK_XK = YEp24PGK_BglII + insert

    assert YEp24PGK_XK.seguid() == "ldseguid=hlyzwrknN5F_ATOvtTCUtQy6YJY"

    YEp24PGK_XK = YEp24PGK_XK.looped()

    assert YEp24PGK_XK.seguid() == "cdseguid=Rszaoc76OKSdw6Q78zj2RZzmR0I"

    YEp24PGK_XK = YEp24PGK_XK.synced("gaattctgaaccagtcctaaaacgagtaaataggaccggcaattc")  # YEp24PGK)

    assert YEp24PGK_XK.seguid() == "cdseguid=Rszaoc76OKSdw6Q78zj2RZzmR0I"

    YEp24PGK_XK_correct = read("YEp24PGK_XK_manually_assembled.txt")

    assert YEp24PGK_XK_correct.seguid() == "cdseguid=Rszaoc76OKSdw6Q78zj2RZzmR0I"
    assert eq(YEp24PGK_XK, YEp24PGK_XK_correct)


def test_melt():

    from pydna.dseqrecord import Dseqrecord
    from Bio.SeqFeature import ExactPosition, SimpleLocation

    dr = Dseqrecord("GATCGATPGAPGCAGATC")
    dr.add_feature(8, 10)

    # GATCGATGGAGGCAGATC
    # CTAGCTA CT CGTCTAG
    #         ==

    dr2, srs = dr.shed_ss_dna(2)
    (seqrec,) = srs
    seqrec.features[0].location == SimpleLocation(ExactPosition(0), ExactPosition(2), strand=1)

    dr = Dseqrecord("GATCAGJGaGEgatc")
    dr.add_feature(7, 10)

    # GATCAG GaGAgatc
    # CTAGTCGCtC ctag
    #        ===

    dr.figure()
    a, b = dr.melt(3)


# def test_apply_cut():

#     from pydna.dseqrecord import Dseqrecord
#     from Bio.SeqFeature import SeqFeature, SimpleLocation
#     from pydna.utils import location_boundaries as _location_boundaries

#     def find_feature_by_id(f: Dseqrecord, id: str) -> SeqFeature:
#         return next(f for f in f.features if f.id == id)

#     # Single cut case, check that features are transmitted correctly.
#     for strand in [1, -1, None]:
#         seq = Dseqrecord("acgtATGaatt", circular=True)
#         seq.features.append(SeqFeature(SimpleLocation(4, 7, strand), id="full_overlap"))
#         seq.features.append(SeqFeature(SimpleLocation(3, 7, strand), id="left_side"))
#         seq.features.append(SeqFeature(SimpleLocation(4, 8, strand), id="right_side"))
#         seq.features.append(SeqFeature(SimpleLocation(3, 10, strand), id="throughout"))
#         for shift in range(len(seq)):
#             seq_shifted = seq.shifted(shift)
#             cut_feature = find_feature_by_id(seq_shifted, "full_overlap")
#             start, end = _location_boundaries(cut_feature.location)
#             # Cut leaving + and - overhangs in the feature full_overlap
#             for dummy_cut in (((start, -3), None), ((end, 3), None)):
#                 open_seq = seq_shifted.apply_cut(dummy_cut, dummy_cut)
#                 assert len(open_seq.features) == 4
#                 new_locs = sorted(str(f.location) for f in open_seq.features)
#                 assert str(open_seq.seq) == "ATGaattacgtATG"
#                 if strand == 1:
#                     assert new_locs == sorted(["[0:3](+)", "[0:4](+)", "[11:14](+)", "[10:14](+)"])
#                 elif strand == -1:
#                     assert new_locs == sorted(["[0:3](-)", "[0:4](-)", "[11:14](-)", "[10:14](-)"])
#                 if strand == None:
#                     assert new_locs == sorted(["[0:3]", "[0:4]", "[11:14]", "[10:14]"])


# def test_apply_cut():

#     from pydna.dseqrecord import Dseqrecord
#     from Bio.SeqFeature import SeqFeature, SimpleLocation
#     from pydna.utils import location_boundaries as _location_boundaries

#     def find_feature_by_id(f: Dseqrecord, id: str) -> SeqFeature:
#         return next(f for f in f.features if f.id == id)

#     # Single cut case, check that features are transmitted correctly.
#     for strand in [1, -1, None]:
#         seq = Dseqrecord("acgtATGaatt", circular=True)
#         seq.features.append(SeqFeature(SimpleLocation(4, 7, strand), id="full_overlap"))
#         seq.features.append(SeqFeature(SimpleLocation(3, 7, strand), id="left_side"))
#         seq.features.append(SeqFeature(SimpleLocation(4, 8, strand), id="right_side"))
#         seq.features.append(SeqFeature(SimpleLocation(3, 10, strand), id="throughout"))
#         for shift in range(len(seq)):
#             seq_shifted = seq.shifted(shift)
#             cut_feature = find_feature_by_id(seq_shifted, "full_overlap")
#             start, end = _location_boundaries(cut_feature.location)
#             # Cut leaving + and - overhangs in the feature full_overlap
#             for dummy_cut in (((start, -3), None), ((end, 3), None)):
#                 open_seq = seq_shifted.apply_cut(dummy_cut, dummy_cut)
#                 assert len(open_seq.features) == 4
#                 new_locs = sorted(str(f.location) for f in open_seq.features)
#                 assert str(open_seq.seq) == "ATGaattacgtATG"
#                 if strand == 1:
#                     assert new_locs == sorted(["[0:3](+)", "[0:4](+)", "[11:14](+)", "[10:14](+)"])
#                 elif strand == -1:
#                     assert new_locs == sorted(["[0:3](-)", "[0:4](-)", "[11:14](-)", "[10:14](-)"])
#                 if strand == None:
#                     assert new_locs == sorted(["[0:3]", "[0:4]", "[11:14]", "[10:14]"])


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s", "-x"])
    # args = [
    #     __file__,
    #     "--cov=pydna",
    #     "--cov-append",
    #     "--cov-report=html:../htmlcov",
    #     "--cov-report=xml",
    #     "--capture=no",
    #     "--durations=10",
    #     "--import-mode=importlib",
    #     "--nbval",
    #     "--current-env",
    #     "--doctest-modules",
    #     "--vvv",
    #     # "-x",
    # ]
    # pytest.main(args)
