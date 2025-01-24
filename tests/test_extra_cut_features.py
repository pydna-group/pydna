#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test cut
"""
import pytest


def test_cut_feat():
    from pydna.readers import read
    from pydna.amplify import pcr

    from pydna.design import primer_design
    from pydna.dseqrecord import Dseqrecord

    from Bio.Restriction import EcoRI

    puc19 = read("pUC19_MarkBudde.gb")

    assert len(puc19.features) == 19
    puc_lin = puc19[:]
    assert len(puc_lin.features) == 19

    ampl = primer_design(puc_lin)
    pf, pr = ampl.forward_primer, ampl.reverse_primer

    # pf
    # Out[3]: f2686 13-mer:5'-TCGCGCGTTTCGG-3'
    # pr
    # Out[4]: r2686 15-mer:5'-GACGAAAGGGCCTCG-3'

    pcrProd = pcr(pf, pr, puc19)
    assert len(pcrProd.features) == 21
    assert len(pcrProd.cut(EcoRI)[1].features) == 16

    # from pydna.editor import ape
    # ape(pcrProd)
    # ape(puc_lin)

    def amplicon_to_dseqrecord(a):
        d = Dseqrecord(a.seq)
        d.features = a.features
        return d

    pcrProdDseqrecord = amplicon_to_dseqrecord(pcrProd)
    assert len(pcrProdDseqrecord.cut(EcoRI)[1].features) == 16


if __name__ == "__main__":
    pytest.cmdline.main([__file__, "-v", "-s"])
