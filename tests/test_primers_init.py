# -*- coding: utf-8 -*-

from pydna.seq import Seq as pydnaSeq
from pydna.primer import Primer
from pydna.parsers import parse_primers
from Bio.SeqFeature import SimpleLocation
from pydna.dseqrecord import Dseqrecord


def test_primer_init():
    primer = Primer("AAAAA")

    parsed_primer, *_ = parse_primers(">primer\nAAAAA")

    loc = SimpleLocation(1, 3)

    assert type(loc.extract(primer).seq) == pydnaSeq

    assert type(loc.extract(parsed_primer).seq) == pydnaSeq

    assert type(Primer(Dseqrecord("AAAAA", id="primer")).seq) == pydnaSeq
