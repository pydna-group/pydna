# -*- coding: utf-8 -*-


def test_primer():
    from pydna import primer

    x = primer.Primer("gtatcatatctatctatcta", footprint=12)

    assert str(x.seq) == "gtatcatatctatctatcta"
    assert str(x.tail) == "gtatcata"
    assert str(x.footprint) == "tctatctatcta"
    assert repr(x) == "id 20-mer:5'-gtatcatatctatctatcta-3'"
    # assert x.tm() == 41.04993874467033

    assert x[0] == "g"

    w = x[2:15]

    assert repr(w) == "part_id 13-mer:5'-atcatatctatct-3'"
    assert str(w.seq) == "atcatatctatct"

    y = "AAA" + w

    assert isinstance(y, type(w))
    assert str(y.seq) == "AAAatcatatctatct"

    # print(str(x.seq).rjust(20))
    # print(str(x.tail))
    # print(str(x.footprint).rjust(20))
    #
    #
    #
    # print(str(w.seq).ljust(20))
    # print(str(w.tail))
    # print(str(w.footprint).rjust(15))
