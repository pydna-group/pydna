# -*- coding: utf-8 -*-
import pytest
from pydna.ligate import ligate
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord

def test_ligate():


    a = Dseqrecord(
        Dseq.from_representation(
            """
                                            GATCaaa
                                                tttTTCC"""
        )
    )
    b = Dseqrecord(
        Dseq.from_representation(
            """
                                            AAGGatta
                                                taatAGGA"""
        )
    )
    c = Dseqrecord(
        Dseq.from_representation(
            """
                                            TCCTccact
                                                ggtgaCTAG"""
        )
    )

    d = Dseqrecord(
        Dseq.from_representation(
            """
                                               Tcgcgc
                                                gcgcgC"""
        )
    )

    e = Dseqrecord(
        Dseq.from_representation(
            """
                                                Gcaatt
                                                 gttaa"""
        )
    )

    fragments = [a, b, c, d, e]
    rcfragments = [a.rc(), b.rc(), c.rc(), d.rc(), e.rc()]

    def list_combinations(a, b):
        N = len(a)
        combinations = []
        for i in range(2**N):
            current = []
            for j in range(N):
                if i & (1 << j):
                    current.append(a[j])
                else:
                    current.append(b[j])
            combinations.append(current)

        return combinations

    # Example usage:
    combinations = list_combinations(fragments, rcfragments)

    for frgs in combinations:
        csequences, lsequences = ligate(frgs)
        for cs in csequences:
            assert cs.seguid() == "cdseguid=ogVOwoSZEk_4kUWnK1NCl5Dv1z4"
            assert len(cs) == 24
        for ss in lsequences:
            assert ss.seguid() == "ldseguid=vBQxmszgfR4b84O0wI7d_ya9uDA"
            assert len(ss) == 12
