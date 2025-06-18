#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Provides the Dseq class for handling double stranded DNA sequences.

Dseq is a subclass of :class:`Bio.Seq.Seq`. The Dseq class
is mostly useful as a part of the :class:`pydna.dseqrecord.Dseqrecord` class
which can hold more meta data.

The Dseq class support the notion of circular and linear DNA topology.
"""

from abc import ABC, abstractmethod
import re
from pydna.utils import rc


class _cas(ABC):
    scaffold = "ND"
    pam = "ND"
    size = 0
    fst5 = 0
    fst3 = 0

    def __init__(self, protospacer):
        self.protospacer = protospacer.upper()
        self.compsite = re.compile(
            f"(?=(?P<watson>{self.protospacer}{self.pam}))|(?=(?P<crick>{rc(self.pam)}{rc(self.protospacer)}))",
            re.UNICODE,
        )

    @abstractmethod
    def search(self, dna, linear=True):
        """To override in subclass."""
        pass

    def __repr__(self):
        return f"{type(self).__name__}({self.protospacer[:3]}..{self.protospacer[-3:]})"

    @abstractmethod
    def __str__(self):
        """To override in subclass."""
        pass


class cas9(_cas):
    """docstring.

    .. code-block::

            |----size----------|

            ---protospacer------
                            -fst3
            fst5             |-|
            |--------------|
                                PAM
        5-NNGGAAGAGTAATACACTA-AAANGGNN-3
        ||||||||||||||||||| ||||||||
        3-NNCCTTCTCATTATGTGAT-TTTNCCNN-5
            ||||||||||||||||| |||
        5-GGAAGAGTAATACACTA-AAAg-u-a-a-g-g  Scaffold
            ---gRNA spacer---    u-a
                                u-a
                                u-a
                                u-a
                                a-u
                                g-u-g
                                a    a
                                g-c-a
                                c-g
                                u-a
                                a-u
                                g   a  tetraloop
                                a-a
    """

    scaffold = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG"
    pam = ".GG"
    size = 20
    fst5 = 17
    fst3 = -3
    ovhg = fst5 - (size + fst3)

    def search(self, dna, linear=True):
        """docstring."""
        dna = str(dna).upper()
        if linear:
            dna = dna
        else:
            dna = dna + dna[1 : self.size]
        results = []
        for mobj in self.compsite.finditer(dna):
            w, c = mobj.groups()
            if w:
                results.append(mobj.start("watson") + 1 + self.fst5)
            if c:
                results.append(mobj.start("crick") + len(self.pam) + 1 - self.fst3)
        return results

    def __str__(self):
        """docstring."""
        return f">{type(self).__name__} protospacer scaffold\n{self.protospacer} {self.scaffold}"


def protospacer(guide_construct, cas=cas9):
    """docstring."""
    in_watson = [
        mobj.group("ps")
        for mobj in re.finditer(
            f"(?P<ps>.{{{cas.size}}})(?:{cas.scaffold})",
            str(guide_construct.seq).upper(),
        )
    ]
    in_crick = [
        rc(mobj.group("ps"))
        for mobj in re.finditer(
            f"(?:{rc(cas.scaffold)})(?P<ps>.{{{cas.size}}})",
            str(guide_construct.seq).upper(),
        )
    ]
    return in_watson + in_crick
