# -*- coding: utf-8 -*-
import pytest


def test_contig(monkeypatch):
    monkeypatch.setenv("pydna_cached_funcs", "")

    from pydna import contig
    from pydna.assembly import Assembly
    from pydna.dseqrecord import Dseqrecord

    a = Dseqrecord("acgatgctatactgCCCCCtgtgctgtgctcta", name="one")
    b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatc", name="two")
    c = Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactg", name="three")
    asm = Assembly((a, b, c), limit=14)

    cnt = asm.assemble_circular()[0]

    assert repr(cnt) == "Contig(o59)"

    assert cnt.detailed_figure() == str(
        "||||||||||||||\n"
        "acgatgctatactgCCCCCtgtgctgtgctcta\n"
        "                   TGTGCTGTGCTCTA\n"
        "                   tgtgctgtgctctaTTTTTtattctggctgtatc\n"
        "                                      TATTCTGGCTGTATC\n"
        "                                      tattctggctgtatcGGGGGtacgatgctatactg\n"
        "                                                           ACGATGCTATACTG\n"
    )

    from textwrap import indent

    fig = """ -|one|14
|      \\/
|      /\\
|      14|two|15
|             \\/
|             /\\
|             15|three|14
|                      \\/
|                      /\\
|                      14-
|                         |
 -------------------------"""

    cnt2 = asm.assemble_linear()[0]

    fig = (
        "one|14\n" "    \\/\n" "    /\\\n" "    14|two|15\n" "           \\/\n" "           /\\\n" "           15|three"
    )

    assert fig == cnt2.figure()

    assert repr(cnt2) == "Contig(-73)"

    # print(repr(cnt2._repr_html_()))

    assert (
        cnt2._repr_html_()
        == "<pre>one|14\n    \\/\n    /\\\n    14|two|15\n           \\/\n           /\\\n           15|three</pre>"
    )

    from unittest.mock import MagicMock

    pp = MagicMock()

    cnt2._repr_pretty_(pp, None)

    pp.text.assert_called_with("Contig(-73)")

    from Bio.Seq import Seq

    from pydna.seqrecord import SeqRecord

    arg = SeqRecord(Seq("aaa"))

    import networkx as nx

    x = contig.Contig.from_SeqRecord(arg, graph=nx.MultiDiGraph())


def test_reverse_complement(monkeypatch):
    from pydna._pretty import pretty_str
    from pydna.assembly import Assembly
    from pydna.dseqrecord import Dseqrecord

    a = Dseqrecord("acgatgctatactgtgCCNCCtgtgctgtgctcta")
    # 12345678901234
    b = Dseqrecord("tgtgctgtgctctaTTTTTTTtattctggctgtatc")
    # 123456789012345
    c = Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactgtg")
    a.name = "aaa"  # 1234567890123456
    b.name = "bbb"
    c.name = "ccc"
    asm = Assembly((a, b, c), limit=14)
    x = asm.assemble_circular()[0]
    y = x.rc()
    z = y.rc()
    assert x.figure() == z.figure()
    assert x.detailed_figure() == z.detailed_figure()

    xfig = """\
 -|aaa|14
|      \\/
|      /\\
|      14|bbb|15
|             \\/
|             /\\
|             15|ccc|16
|                    \\/
|                    /\\
|                    16-
|                       |
 -----------------------
     """.rstrip()

    xdfig = pretty_str(
        """\
||||||||||||||||
acgatgctatactgtgCCNCCtgtgctgtgctcta
                     TGTGCTGTGCTCTA
                     tgtgctgtgctctaTTTTTTTtattctggctgtatc
                                          TATTCTGGCTGTATC
                                          tattctggctgtatcGGGGGtacgatgctatactgtg
                                                               ACGATGCTATACTGTG
    """.rstrip()
        + "\n"
    )

    assert x.figure() == xfig
    assert x.detailed_figure() == xdfig

    yfig = """\
 -|ccc_rc|15
|         \\/
|         /\\
|         15|bbb_rc|14
|                   \\/
|                   /\\
|                   14|aaa_rc|16
|                             \\/
|                             /\\
|                             16-
|                                |
 --------------------------------
     """.rstrip()

    ydfig = (
        """\
||||||||||||||||
cacagtatagcatcgtaCCCCCgatacagccagaata
                      GATACAGCCAGAATA
                      gatacagccagaataAAAAAAAtagagcacagcaca
                                            TAGAGCACAGCACA
                                            tagagcacagcacaGGNGGcacagtatagcatcgt
                                                               CACAGTATAGCATCGT
    """.rstrip()
        + "\n"
    )

    assert y.figure() == yfig
    assert y.detailed_figure() == ydfig


def test_linear(monkeypatch):
    from pydna._pretty import pretty_str
    from pydna.assembly import Assembly
    from pydna.dseqrecord import Dseqrecord

    a = Dseqrecord("acgatgctatactgtgCCNCCtgtgctgtgctcta")
    # 12345678901234
    b = Dseqrecord("tgtgctgtgctctaTTTTTTTtattctggctgtatc")
    # 123456789012345
    c = Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactgtg")
    a.name = "aaa"  # 1234567890123456
    b.name = "bbb"
    c.name = "ccc"
    asm = Assembly((a, b, c), limit=14)
    x = asm.assemble_linear()[0]

    answer = "aaa|14\n    \\/\n    /\\\n    14|bbb|15\n           \\/\n           /\\\n           15|ccc"

    assert x.figure() == answer.strip()
    answer = "acgatgctatactgtgCCNCCtgtgctgtgctcta\n                     TGTGCTGTGCTCTA\n                     tgtgctgtgctctaTTTTTTTtattctggctgtatc\n                                          TATTCTGGCTGTATC\n                                          tattctggctgtatcGGGGGtacgatgctatactgtg\n"
    assert x.detailed_figure()


def test_rich(monkeypatch):

    from pydna.readers import read
    from pydna.amplify import pcr
    from pydna.parsers import parse_primers
    from pydna.myprimers import PrimerList

    p = {}

    p[1113], p[987], p[1196], p[1195], p[978], p[977], p[984], p[983], p[1804], p[1347] = parse_primers(
        '''

    >1113_Amp.fw.nw 55-mer
    GAAAAGCGTTTACCTCGGAACTCTATTGTAGAACCCCTATTTGTTTATTTTTCTA

    >987_Amp.REV 52-merAmp.REV
    AGAAAGTCTACACCTTACGCTGATTGGATTTGAAGTTTTAAATCAATCTAAA

    >1196_Pbr.FW 49-mer
    AATCCAATCAGCGTAAGGTGTAGACTTTCTCTGTCAGACCAAGTTTACT

    >1195_Pbr.REV 44-mer
    GTTGACTACTATTTACGCAGCAGATACGATCTCGTTTCATCGGT

    >978_Crp.FW 46-merCrp.FW
    AACTGTAAAATCAGGTATCTCGTAGTCCGTGTTCTGATCCTCGAGC

    >977_Crp.REV 47-merCrp.REV
    TACAATAGAGTTCCGAGGTAAACGCTTTTCGTTCTTGTCTCATTGCC

    >984_Micron.FW 49-merMicron.FW
    ATCGTATCTGCTGCGTAAATAGTAGTCAACGATCGTACTTGTTACCCAT

    >983_Micron.REV 48-merMicron.REV
    CAGAGCAGACAGTTCCTTTACGAGATTTTAGATCCAATATCAAAGGAA

    >1804_TRP1fp_pTA 51-merreplaces 1348 & 1680
    TAAAATCTCGTAAAGGAACTGTCTGCTCTGtataaaaataggcgtatcacg

    >1347_TRP1rp_pTA 51-mer
    ACGGACTACGAGATACCTGATTTTACAGTTGATCTTTTATGCTTGCTTTTC

    '''
    )

    lol = [
        ['pBR322.gb', 1113, 987, 'amp'],
        ['pBR322.gb', 1196, 1195, 'pBR'],
        ['pYPKpw.gb', 978, 977, 'Δcrp'],
        ['YEplac195_snapgene.gb', 984, 983, '2µ'],
        ['YIplac204_snapgene.gb', 1804, 1347, 'TRP1'],
    ]

    fragments = []

    for row in lol:
        template, fp, rp, target = row
        fp = p[int(fp)]
        rp = p[int(rp)]
        template = read(template.strip())
        fragment = pcr(fp, rp, template)
        fragment.name = target.strip()
        fragments.append(fragment)

    from pydna.assembly import Assembly

    asm = Assembly(fragments)

    cps = asm.assemble_circular()

    cp = cps[0]

    result = cp.graphic_figure()
    # TODO: test for graph content
    # cp.graphic_figure_plotly()


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
