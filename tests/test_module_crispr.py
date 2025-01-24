#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_crispr():

    from pydna.crispr import cas9, protospacer
    from pydna.dseqrecord import Dseqrecord

    a = Dseqrecord("GTTACTTTACCCGACGT")
    b = Dseqrecord("CCCaGG")

    target = a + b

    scaffold = "gttttagagctagaaatagcaagttaaaataagg"

    containing_sgRNA = target[:20] + scaffold

    (ps,) = protospacer(containing_sgRNA)

    assert ps == "GTTACTTTACCCGACGTCCC"

    c9 = cas9(ps)

    assert c9.__dict__ == cas9(ps.lower()).__dict__

    # assert target.cut(c9) == (a, b)

    assert c9.search(target) == [17]
    assert c9.search(target.seq) == [17]
    assert c9.search(str(target.seq)) == [17]
    assert c9.search(str(target.seq).lower()) == [17]

    containing_sgRNA_rc = containing_sgRNA.rc()

    (ps_rc,) = protospacer(containing_sgRNA_rc)

    assert ps_rc == "GTTACTTTACCCGACGTCCC"

    target_rc = target.rc()

    assert c9.search(target_rc) == [6]
    assert c9.search(target_rc.seq) == [6]
    assert c9.search(str(target_rc.seq)) == [6]
    assert c9.search(str(target_rc.seq).lower()) == [6]

    assert [s.seq for s in target_rc.cut(c9)] == [b.seq.rc(), a.seq.rc()]


# CCtGGGACGTCGGGTAAAGTAAC
# GGaCCCTGCAGCCCATTTCATTG

#    CCCTGCAGCCCATTTCATTG

# CCtGG
# GGaCC

#      GACGTCGGGTAAAGTAAC
#      CTGCAGCCCATTTCATTG

if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
