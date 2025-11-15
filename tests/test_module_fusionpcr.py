# -*- coding: utf-8 -*-
import pytest

from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.utils import eq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from pydna.fusionpcr import fuse_by_pcr


x = Dseqrecord(
    Dseq("tctggtcaagctgaagggtattc"), features=[SeqFeature(SimpleLocation(5, 10, strand=1), type="misc_feature")]
)
y = Dseqrecord(Dseq("tattcgtacacagatg"), features=[SeqFeature(SimpleLocation(5, 10, strand=1), type="misc_feature")])
z = Dseqrecord(Dseq("acagatgacgtgt"), features=[SeqFeature(SimpleLocation(5, 10, strand=1), type="misc_feature")])
n = Dseqrecord(Dseq("cgatcaaaaaaacc"))
r = Dseqrecord(Dseq("tctggtcaagctgaagggtattcgtacacagatgacgtgt"))


def test_simple_case():
    result, *rest = fuse_by_pcr((x, y, z), limit=5)
    assert eq(result, r)


def test_fusionpcr1():
    """
    tctggtcaagctgaagggtattcgtacacagatgacgtgt   (40bp)

    tctggtcaagctgaagggtattc
    agaccagttcgacttcccataag
                      tattcgtacacagatg
                      ataagcatgtgtctac
                               acagatgacgtgt
                               tgtctactgcaca
    """

    seqtuples = [
        (x, y, z),
        (x, z, y),
        (y, x, z),
        (y, z, x),
        (z, x, y),
        (z, y, x),
        (x.rc(), y.rc(), z.rc()),
        (x.rc(), z.rc(), y.rc()),
        (y.rc(), x.rc(), z.rc()),
        (y.rc(), z.rc(), x.rc()),
        (z.rc(), x.rc(), y.rc()),
        (z.rc(), y.rc(), x.rc()),
    ]

    for arg in seqtuples:
        result = fuse_by_pcr(arg, limit=5).pop()
        assert eq(result, r)


def test_fusionpcr2():
    seqtuples = [
        (x.rc(), y, z),
        (x, z.rc(), y),
        (y, x, z.rc()),
        (y.rc(), z, x),
        (z, x.rc(), y),
        (z, y, x.rc()),
    ]

    seqtuples = [(x, y, z.rc())]
    for arg in seqtuples:
        result = fuse_by_pcr(arg, limit=5).pop()
        assert eq(result, r)


def test_fusionpcr3():
    seqtuples = [
        (x, y, z, n),
    ]
    for arg in seqtuples:
        result = fuse_by_pcr(arg, limit=5).pop()
        assert eq(result, r)




x = Dseqrecord(
    Dseq("tctggtcaagctgaagggtattc"), features=[SeqFeature(SimpleLocation(5, 10, strand=1), type="misc_feature")]
)
y = Dseqrecord(Dseq("tattcgtacacagatg"), features=[SeqFeature(SimpleLocation(5, 10, strand=1), type="misc_feature")])
z = Dseqrecord(Dseq("acagatgacgtgt"), features=[SeqFeature(SimpleLocation(5, 10, strand=1), type="misc_feature")])


seqtuples = [(x, y, z)]
for arg in seqtuples:
    result = fuse_by_pcr(arg, limit=5).pop()
