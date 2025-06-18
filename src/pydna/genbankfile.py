#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
from pydna.dseqrecord import Dseqrecord as _Dseqrecord


class GenbankFile(_Dseqrecord):
    def __init__(self, record, *args, path=None, **kwargs):
        super().__init__(record, *args, **kwargs)
        self.path = path

    @classmethod
    def from_SeqRecord(cls, record, *args, path=None, **kwargs):
        obj = super().from_SeqRecord(record, *args, **kwargs)
        obj.path = path
        return obj

    def __repr__(self):
        """returns a short string representation of the object"""
        return "File({})({}{})".format(
            self.id, {True: "-", False: "o"}[not self.circular], len(self)
        )

    def _repr_pretty_(self, p, cycle):
        """returns a short string representation of the object"""
        p.text(
            "File({})({}{})".format(
                self.id, {True: "-", False: "o"}[not self.circular], len(self)
            )
        )

    def _repr_html_(self):
        return "<a href='{path}' target='_blank'>{path}</a><br>".format(path=self.path)

    def reverse_complement(self):
        answer = type(self)(super().reverse_complement(), path=self.path)
        return answer

    rc = reverse_complement
