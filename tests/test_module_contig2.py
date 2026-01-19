# -*- coding: utf-8 -*-
import pytest
from pydna.dseqrecord import Dseqrecord
from pydna.assembly import Assembly
from pydna import contig
from pydna._pretty import pretty_str
from pydna.amplify import pcr
from pydna.parsers import parse_primers
from pydna.readers import read


a = Dseqrecord("acgatgctatactgCCCCCtgtgctgtgctcta", name="one")
b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatc", name="two")
c = Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactg", name="three")


@pytest.mark.mpl_image_compare
def test_contig_linear():

    asm = Assembly((a, b, c), limit=14)

    cont = asm.assemble_linear()[0]

    fig = (
        "one|14\n"
        "    \\/\n"
        "    /\\\n"
        "    14|two|15\n"
        "           \\/\n"
        "           /\\\n"
        "           15|three"
    )

    assert fig == cont.figure()

    assert repr(cont) == "Contig(-73)"

    assert cont._repr_html_() == f"<pre>{fig}</pre>"

    assert cont.detailed_figure() == (
        "acgatgctatactgCCCCCtgtgctgtgctcta\n"
        "                   TGTGCTGTGCTCTA\n"
        "                   tgtgctgtgctctaTTTTTtattctggctgtatc\n"
        "                                      TATTCTGGCTGTATC\n"
        "                                      tattctggctgtatcGGGGGtacgatgctatactg\n"
    )

    return cont.figure_mpl()


@pytest.mark.mpl_image_compare
def test_contig_circular():

    asm = Assembly((a, b, c), limit=14)

    cont = asm.assemble_circular()[0]

    assert repr(cont) == "Contig(o59)"

    assert cont.detailed_figure() == (
        "||||||||||||||\n"
        "acgatgctatactgCCCCCtgtgctgtgctcta\n"
        "                   TGTGCTGTGCTCTA\n"
        "                   tgtgctgtgctctaTTTTTtattctggctgtatc\n"
        "                                      TATTCTGGCTGTATC\n"
        "                                      tattctggctgtatcGGGGGtacgatgctatactg\n"
        "                                                           ACGATGCTATACTG\n"
    )

    fig = (
        " -|one|14\n"
        "|      \\/\n"
        "|      /\\\n"
        "|      14|two|15\n"
        "|             \\/\n"
        "|             /\\\n"
        "|             15|three|14\n"
        "|                      \\/\n"
        "|                      /\\\n"
        "|                      14-\n"
        "|                         |\n"
        " -------------------------"
    )

    assert fig == cont.figure()

    return cont.figure_mpl()

    from unittest.mock import MagicMock

    pp = MagicMock()

    cont._repr_pretty_(pp, None)

    pp.text.assert_called_with("Contig(o59)")

    from Bio.Seq import Seq

    from pydna.seqrecord import SeqRecord

    arg = SeqRecord(Seq("aaa"))

    import networkx as nx

    x = contig.Contig.from_SeqRecord(arg, graph=nx.MultiDiGraph())


def test_reverse_complement():

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


def test_linear():

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


@pytest.mark.mpl_image_compare
def test_mpl1():

    a = Dseqrecord(
        "GCCGATTCTCATCCGGGCCT" "CACTATAGGACCATTCCGTT" "GCCCTGCGCTGCGCTGTATA", name="1"
    )

    b = Dseqrecord(
        "GCCCTGCGCTGCGCTGTATA" "TCCCCAGGAACAGACTTCCT" "GCCGATTCTCATCCGGGCCT", name="2"
    )

    asm = Assembly((a, b), limit=20)
    cps = asm.assemble_circular()
    cp = cps[0]
    return cp.figure_mpl()


@pytest.mark.mpl_image_compare
def test_mpl2():

    x = Dseqrecord(
        "CTAAGATATTCTTACGTGTA" "CACTATAGGACCATTCCGTT" "GCCCTGCGCTGCGCTGTATA", name="x"
    )

    y = Dseqrecord(
        "GCCCTGCGCTGCGCTGTATA" "TCCCCAGGAACAGACTTCCT" "GCCGATTCTCATCCGGGCCT", name="y"
    )

    z = Dseqrecord(
        "GCCGATTCTCATCCGGGCCT" "GGATGCAATGCGATCCTCCG" "CTAAGATATTCTTACGTGTA", name="z"
    )

    asm2 = Assembly((x, y, z), limit=20)
    cps2 = asm2.assemble_circular()
    cp2 = cps2[0]
    return cp2.figure_mpl()


@pytest.mark.mpl_image_compare
def test_mpl3():

    p = {}

    (
        p[1113],
        p[987],
        p[1196],
        p[1195],
        p[978],
        p[977],
        p[984],
        p[983],
        p[1804],
        p[1347],
    ) = parse_primers(
        """

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

    """
    )

    lol = [
        ["pBR322.gb", 1113, 987, "amp"],
        ["pBR322.gb", 1196, 1195, "pBR"],
        ["pYPKpw.gb", 978, 977, "Δcrp"],
        ["YEplac195_snapgene.gb", 984, 983, "2µ"],
        ["YIplac204_snapgene.gb", 1804, 1347, "TRP1"],
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

    asm = Assembly(fragments)

    cps = asm.assemble_circular()

    cp = cps[0]

    return cp.figure_mpl()


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
