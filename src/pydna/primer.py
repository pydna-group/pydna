#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""This module provide the Primer class that is a subclass of the biopython SeqRecord."""

from pydna.seq import Seq
from pydna.seqrecord import SeqRecord


class Primer(SeqRecord):
    """Primer and its position on a template, footprint and tail."""

    def __init__(
        self, record, *args, amplicon=None, position=None, footprint=0, **kwargs
    ):
        if hasattr(record, "features"):  # Seqrecord
            self.__dict__.update(record.__dict__)
            self.__dict__.update(kwargs)
        elif hasattr(record, "transcribe"):  # Seq
            super().__init__(record, *args, **kwargs)
        else:  # string?
            super().__init__(Seq(record), *args, **kwargs)
        self.amplicon = amplicon
        self.position = position
        self._fp = footprint or len(record)

    @property
    def footprint(self):
        return self.seq[-self._fp :] if self._fp else ""

    @property
    def tail(self):
        return self.seq[: -self._fp] if self._fp else ""

    def __repr__(self):
        s = min((self.seq, "{}..{}".format(self.seq[:15], self.seq[-3:])), key=len)
        return "{id} {len}-mer:5'-{seq}-3'".format(id=self.id, len=len(self), seq=s)

    def __radd__(self, other):
        new = super().__radd__(other)
        return Primer(new, position=self.position, footprint=self._fp)

    def __getitem__(self, index):
        result = super().__getitem__(index)
        if hasattr(index, "indices"):  # index is a slice
            i1, i2, i3 = index.indices(len(self))
            j1, j2, j3 = slice(-(self._fp or 0), None).indices(len(self))
            result._fp = self._fp - (i1 - j1 > 0) * abs(i1 - j1)
        return result

    def reverse_complement(self, *args, **kwargs):
        """Return the reverse complement of the sequence."""
        answer = super().reverse_complement(*args, **kwargs)
        answer.amplicon = None
        answer.position = None
        answer._fp = len(self)
        return answer
