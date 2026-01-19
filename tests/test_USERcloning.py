#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydna.dseq import Dseq
from pydna.parsers import parse, parse_primers
from pydna.amplify import pcr
from textwrap import dedent


def test_USER_cloning():

    # >a
    # GUGGATT

    primers = """
    >a
    GUGGATT

    >b
    cUCGCCG
    """
    template = """
    >temp
    gtGGATTaaaCGGCGag
    """
    fp, rp = parse_primers(primers)
    te, *rest = parse(template)
    te.add_feature()
    p = pcr((fp, rp, te), limit=5)
    assert p.seq._data == b"GUGGATTaaaCGGCGOg"

    figure = p.figure()

    correct_figure = dedent(
        """\
                               5gtGGATT...CGGCGag3
                                          |||||||
                                         3GCCGCUc5
                               5GUGGATT3
                                |||||||
                               3caCCTAA...GCCGCtc5"""
    )
    assert figure == correct_figure

    w = "GUGGATTaaaCGGCGAg"
    c = "CACCTAAtttGCCGCUc"

    assert p.seq.watson == w
    assert p.seq.crick == c[::-1]

    assert te.features[0] == p.features[0]

    USERprocessed = p.seq.user()

    correct_figure = dedent(
        """\
    Dseq(-17)
    G GGATTaaaCGGCGAg
    CACCTAAtttGCCGC c
    """
    ).strip()

    assert repr(USERprocessed) == correct_figure

    melted1 = USERprocessed.melt(1)

    melted12 = USERprocessed.melt(12)

    assert melted1 == melted12

    stuffer, insert, stuffer = melted1

    correct_figure = dedent(
        """\
    Dseq(-17)
      GGATTaaaCGGCGAg
    CACCTAAtttGCCGC
    """
    ).strip()

    assert repr(insert) == correct_figure

    plasmid = Dseq.from_representation(
        """
                                       Dseq(-7)
                                         aaaGT
                                       Tcttt
                                       """
    )

    plasmid_insert = (plasmid + insert).looped()

    assert plasmid_insert == Dseq("AgaaaGTGGATTaaaCGGCG", circular=True)
