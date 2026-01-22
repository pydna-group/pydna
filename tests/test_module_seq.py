# -*- coding: utf-8 -*-
import pytest
from pydna.seq import Seq


def test_seguid():
    assert Seq("gattaca").seguid() == "lsseguid=tp2jzeCM2e3W4yxtrrx09CMKa_8"


def test_gc():
    assert Seq("atgtaa").gc() == 0.167


def test_cai():
    pytest.importorskip("cai2")
    assert Seq("atgtaa").cai() == 1.0


def test_rare_codons():
    pytest.importorskip("cai2")

    from pydna.codon import rare_codons

    lol = {}

    lol["sce"] = [
        [
            "cds",
            "len",
            "cai",
            "gc",
            "sta",
            "stp",
            "n-end",
            "CGA",
            "CGG",
            "CGC",
            "CCG",
            "CTC",
            "GCG",
            "rare",
        ],
        ["ATG...TAA", 8.0, 0.219, 0.708, 1.0, 0.47, "2 min", 1, 1, 1, 1, 1, 1, 0.75],
    ]

    lol["eco"] = [
        [
            "cds",
            "len",
            "cai",
            "gc",
            "sta",
            "stp",
            "n-end",
            "CGA",
            "CGG",
            "CGC",
            "CCG",
            "CTC",
            "GCG",
            "rare",
        ],
        ["ATG...TAA", 10.0, 0.387, 0.5, 1.0, 0.47, "2 min", 1, 1, 0, 0, 0, 0, 0.2],
    ]

    for organism, codons in rare_codons.items():
        s = Seq("atg" + "".join(codons) + "taa")
        slices = s.rarecodons(organism=organism)
        for slc in slices:
            assert s[slc].upper() in codons
        assert s.express().lol() == lol[organism]


def test_startcodon():
    assert Seq("atgtaa").startcodon() == 1.0


def test_stopcodon():
    assert Seq("atgtaa").stopcodon() == 0.47


def test_orf():
    s = Seq("atgaaattttaa")
    assert s.orfs(2) == [(0, 12)]


def test_no_orf():
    s2 = Seq("aaaaaaaaaaaaaaa")
    assert s2.orfs(2) == []


#    >>> import warnings
#    >>> from pydna import _PydnaWarning
#    >>> warnings.simplefilter('ignore', _PydnaWarning)

# obj.name = "1111111111111111"


#    def __init__(self, seq, id="<unknown id>", name="<unknown name>",
#                 description="<unknown description>", dbxrefs=None,
#                 features=None, annotations=None,
#                 letter_annotations=None):
#        """Create a SeqRecord.
