#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
# doctest: +NORMALIZE_WHITESPACE
# doctest: +SKIP
"""This module provides the :class:`Amplicon` class for PCR simulation.
This class is not meant to be use directly but is
used by the :mod:`amplify` module"""

from pydna.tm import dbd_program
from pydna.tm import program
from pydna.primer import Primer
from pydna._pretty import pretty_str
from pydna.dseqrecord import Dseqrecord
from pydna.seqrecord import SeqRecord
import textwrap
import copy


class Amplicon(Dseqrecord):
    """The Amplicon class holds information about a PCR reaction involving two
    primers and one template. This class is used by the Anneal class and is not
    meant to be instantiated directly.

    Parameters
    ----------
    forward_primer : SeqRecord(Biopython)
        SeqRecord object holding the forward (sense) primer

    reverse_primer : SeqRecord(Biopython)
        SeqRecord object holding the reverse (antisense) primer

    template : Dseqrecord
        Dseqrecord object holding the template (circular or linear)


    """

    def __init__(
        self,
        record,
        *args,
        template=None,
        forward_primer=None,
        reverse_primer=None,
        **kwargs,
    ):
        super().__init__(record, *args)
        self.template = template
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.__dict__.update(kwargs)

        # https://medium.com/@chipiga86/circular-references-without-memory-
        # leaks-and-destruction-of-objects-in-python-43da57915b8d

    @classmethod
    def from_SeqRecord(cls, record, *args, path=None, **kwargs):
        obj = super().from_SeqRecord(record, *args, **kwargs)
        obj.path = path
        return obj

    def __getitem__(self, sl):
        answer = copy.copy(self)
        answer.seq = answer.seq.__getitem__(sl)
        # answer.seq.alphabet = self.seq.alphabet
        sr = SeqRecord("n" * len(self))
        sr.features = self.features
        answer.features = SeqRecord.__getitem__(sr, sl).features
        return answer

    def __repr__(self):
        """returns a short string representation of the object"""
        return "Amplicon({})".format(self.__len__())

    def _repr_pretty_(self, p, cycle):
        p.text("Amplicon({})".format(self.__len__()))

    def _repr_html_(self):
        return "Amplicon({})".format(self.__len__())

    def reverse_complement(self):
        r = type(self)(super().reverse_complement())
        r.template = self.template.rc()
        r.forward_primer = copy.copy(self.reverse_primer)
        r.reverse_primer = copy.copy(self.forward_primer)
        r.forward_primer.position, r.reverse_primer.position = (
            r.reverse_primer.position,
            r.forward_primer.position,
        )
        return r

    rc = reverse_complement

    def figure(self):
        """
        This method returns a simple figure of the two primers binding
        to a part of the template.

        ::

         5tacactcaccgtctatcattatc...cgactgtatcatctgatagcac3
                                    ||||||||||||||||||||||
                                   3gctgacatagtagactatcgtg5
         5tacactcaccgtctatcattatc3
          |||||||||||||||||||||||
         3atgtgagtggcagatagtaatag...gctgacatagtagactatcgtg5



        Returns
        -------
        figure:string
             A string containing a text representation of the primers
             annealing on the template (see example above).
        """

        fp = self.forward_primer
        rp = self.reverse_primer
        tp = self.template
        ft = len(fp) - fp._fp  # forward tail length
        # rt = len(rp) - rp._fp  # reverse tail length
        faz = tp[fp.position - fp._fp : fp.position].seq
        raz = tp[rp.position : rp.position + rp._fp].seq
        sp3 = " " * (len(fp.seq) + 3)
        # breakpoint()
        fzc = tp.seq.crick[::-1][fp.position - fp._fp : fp.position]
        rzc = tp.seq.crick[::-1][rp.position : rp.position + rp._fp]
        f = f"""
            {" " * ft}5{faz}...{raz}3
             {sp3}{"|" * rp._fp}
            {sp3}3{rp.seq[::-1]}5
            5{fp.seq}3
             {"|" * fp._fp:>{len(fp)}}
            {" " * ft}3{fzc}...{rzc}5
            """
        # breakpoint()
        return pretty_str(textwrap.dedent(f).strip("\n"))

    def set_forward_primer_footprint(self, length):
        self.forward_primer = Primer(
            self.forward_primer.tail + self.seq[:length], footprint=length
        )

    def set_reverse_primer_footprint(self, length):
        self.reverse_primer = Primer(
            self.reverse_primer.tail + self.seq[:length], footprint=length
        )

    def program(self):
        return program(self)

    def dbd_program(self):
        return dbd_program(self)

    def primers(self):
        return self.forward_primer, self.reverse_primer
