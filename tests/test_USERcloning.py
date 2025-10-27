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
    assert p.seq.crick ==  c[::-1]

    assert te.features[0] == p.features[0]

    USERprocessed = p.seq.user()

    assert USERprocessed.melt(1) == USERprocessed.melt(12)

    stuffer, insert, stuffer = USERprocessed.melt(1)

    plasmid = Dseq.from_representation("""
                                       Dseq(-7)
                                         aaaGT
                                       Tcttt
                                       """)

    plasmid_insert = (plasmid + insert).looped()

    assert plasmid_insert == Dseq("AgaaaGTGGATTaaaCGGCG", circular=True)
