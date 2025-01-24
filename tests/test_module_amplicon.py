# -*- coding: utf-8 -*-
import pytest


def test_amplicon():
    from pydna.amplify import Anneal
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer
    import textwrap

    template = Dseqrecord("AAAtacactcaccgtctatcattatctactatcgactgtatcatctgatagcacTTT")

    p1 = Primer("CCCtacactcaccgtctatcattatc")
    p2 = Primer("GGGgtgctatcagatgatacagtcg")

    ann = Anneal((p1, p2), template)

    prod = ann.products[0]

    assert str(prod.seq) == str(("CCC" + template[3:-3] + "CCC").seq)

    assert prod.forward_primer.seq == p1.seq
    assert prod.reverse_primer.seq == p2.seq

    prod.forward_primer.__dict__

    assert repr(prod) == "Amplicon(57)"

    assert prod._repr_html_() == "Amplicon(57)"

    from unittest.mock import MagicMock

    pp = MagicMock()

    prod._repr_pretty_(pp, None)

    pp.text.assert_called_with("Amplicon(57)")

    fig = """\
       5tacactcaccgtctatcattatc...cgactgtatcatctgatagcac3
                                  ||||||||||||||||||||||
                                 3gctgacatagtagactatcgtgGGG5
    5CCCtacactcaccgtctatcattatc3
        |||||||||||||||||||||||
       3atgtgagtggcagatagtaatag...gctgacatagtagactatcgtg5"""

    assert prod.figure() == textwrap.dedent(fig)

    ctemplate = template.looped()
    cann = Anneal((p1, p2), ctemplate)

    cprod = cann.products[0]
    assert repr(cprod) == "Amplicon(57)"

    assert cprod.figure() == textwrap.dedent(fig)

    # assert prod.program() == prod.taq_program()

    # assert prod.pfu_sso7d_program() == prod.dbd_program()

    from pydna.amplicon import Amplicon

    from Bio.Seq import Seq

    from pydna.seqrecord import SeqRecord

    arg = SeqRecord(Seq("aaa"))

    x = Amplicon.from_SeqRecord(arg)


def test_amplicon_dbd():
    from pydna.amplify import Anneal
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer
    from textwrap import dedent

    template = Dseqrecord("GCGTCCAGCGGCTGCCCGAGGCGCCAAGTG" + "GATC" * 360 + "CCCGGGCCGAGCCCGCATCTGAGGCCGCCGCGGGC")

    p1 = Primer("GCGTCCAGCGGCTGCCCGAGG")
    p2 = Primer("GCCCGCGGCGGCCTCAGATGCGG")

    ann = Anneal((p1, p2), template)

    prod = ann.products[0]

    assert repr(prod) == "Amplicon(1505)"

    fig = r"""
    |95°C|95°C               |    |tmf:77.1
    |____|_____          72°C|72°C|tmr:80.9
    |3min|30s  \ 64.6°C _____|____|45s/kb
    |    |      \______/ 1:07|5min|GC 51%
    |    |       30s         |    |1505bp
    """
    fig = dedent(fig).strip()
    assert str(prod.program()) == fig

    fig = r"""
    |98°C|98°C      |    |tmf:71.6
    |____|____      |    |tmr:75.3
    |30s |10s \ 72°C|72°C|15s/kb
    |    |     \____|____|GC 51%
    |    |      0:22|5min|1505bp
    """
    fig = dedent(fig).strip()
    assert str(prod.dbd_program()) == fig


def test_amplicon_dbd_low_gc():
    from pydna.amplify import Anneal
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer
    from textwrap import dedent

    template = Dseqrecord("AAAATATTTTTATACAT" + "GAAA" * 370 + "ATAAAAAATACGTCATC")

    p1 = Primer("AAAATATTTTTATACAT")
    p2 = Primer("GATGACGTATTTTTTAT")

    ann = Anneal((p1, p2), template)

    prod = ann.products[0]

    assert repr(prod) == "Amplicon(1514)"

    fig = r"""
    |98°C|98°C               |    |tmf:32.6
    |____|_____          72°C|72°C|tmr:39.6
    |30s |10s  \ 35.6°C _____|____|15s/kb
    |    |      \______/ 0:22|5min|GC 24%
    |    |       10s         |    |1514bp
    """

    fig = dedent(fig).strip()

    assert str(prod.dbd_program()) == fig


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
