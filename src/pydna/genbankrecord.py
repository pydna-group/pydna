#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from pydna._pretty import pretty_str as _ps
import os as _os


class GenbankRecord(_Dseqrecord):
    def __init__(
        self, record, *args, item="accession", start=None, stop=None, strand=1, **kwargs
    ):
        super().__init__(record, *args, **kwargs)
        self.item = item
        self.start = start
        self.stop = stop
        self.strand = strand
        self._repr = item
        if self.start is not None and self.stop is not None:
            self._repr += " {}-{}".format(self.start, self.stop)
        self._linktemplate = "<a href='https://www.ncbi.nlm.nih.gov/nuccore/{item}?from={start}&to={stop}&strand={strand}' target='_blank'>{text}</a>"
        self.hyperlink = _ps(
            self._linktemplate.format(
                item=self.item,
                start=self.start or "",
                stop=self.stop or "",
                strand=self.strand,
                text=self._repr,
            )
        )

    @classmethod
    def from_string(
        cls,
        record: str = "",
        *args,
        item="accession",
        start=None,
        stop=None,
        strand=1,
        **kwargs,
    ):
        """docstring."""
        obj = super().from_string(record, *args, **kwargs)
        obj.item = item
        obj.start = start
        obj.stop = stop
        obj.strand = strand
        obj._repr = item
        if obj.start is not None and obj.stop is not None:
            obj._repr += " {}-{}".format(obj.start, obj.stop)
        obj._linktemplate = "<a href='https://www.ncbi.nlm.nih.gov/nuccore/{item}?from={start}&to={stop}&strand={strand}' target='_blank'>{text}</a>"
        obj.hyperlink = _ps(
            obj._linktemplate.format(
                item=obj.item,
                start=obj.start or "",
                stop=obj.stop or "",
                strand=obj.strand,
                text=obj._repr,
            )
        )
        return obj

    @classmethod
    def from_SeqRecord(
        cls, record, *args, item="accession", start=None, stop=None, strand=1, **kwargs
    ):
        obj = super().from_SeqRecord(record, *args, **kwargs)
        obj.item = item
        obj.start = start
        obj.stop = stop
        obj.strand = strand
        obj._repr = item
        if obj.start is not None and obj.stop is not None:
            obj._repr += " {}-{}".format(obj.start, obj.stop)
        obj._linktemplate = "<a href='https://www.ncbi.nlm.nih.gov/nuccore/{item}?from={start}&to={stop}&strand={strand}' target='_blank'>{text}</a>"
        obj.hyperlink = _ps(
            obj._linktemplate.format(
                item=obj.item,
                start=obj.start or "",
                stop=obj.stop or "",
                strand=obj.strand,
                text=obj._repr,
            )
        )
        return obj

    def __getitem__(self, sl):
        answer = super().__getitem__(sl)
        answer.item = self.item
        answer.start = (self.start or 0) + (sl.start or 0)
        answer.stop = (self.start or 0) + (sl.stop or 0)
        answer.strand = self.strand
        return answer

    def __repr__(self):
        """returns a short string representation of the object"""
        return "Gbnk({}{} {})".format(
            {True: "-", False: "o"}[not self.circular], len(self), self._repr
        )

    def _repr_pretty_(self, p, cycle):
        """returns a short string representation of the object"""
        p.text(self.__repr__())

    def _repr_html_(self):
        return self.hyperlink

    def reverse_complement(self):
        answer = type(self)(
            super().reverse_complement(),
            item=self.item,
            start=self.start,
            stop=self.stop,
            strand={1: 2, 2: 1}[self.strand],
        )
        return answer

    rc = reverse_complement

    def pydna_code(self):
        """docstring."""  # FIXME

        code = (
            "from pydna.genbank import Genbank\n"
            f"gb = Genbank('{_os.getenv('pydna_email')}')\n"
            f"seq = gb.nucleotide('{self.item}'"
        )
        if self.start and self.start:
            code += (
                ",\n"
                f"                    seq_start={self.start},\n"
                f"                    seq_stop={self.stop},\n"
                f"                    strand={self.strand})"
            )
        else:
            code += ")"

        return _ps(code)

    def biopython_code(self):
        """docstring."""  # FIXME

        code = (
            "from Bio import Entrez, SeqIO\n"
            f"Entrez.email = '{_os.getenv('pydna_email')}'\n"
            "handle = Entrez.efetch(db='nuccore',\n"
            f"                       id='{self.item}',\n"
            "                       rettype='gbwithparts',\n"
            "                       retmode='text',"
        )
        if self.start and self.stop:
            code += (
                "\n"
                f"                       seq_start={self.start},\n"
                f"                       seq_stop={self.stop},\n"
                f"                       strand={self.strand})\n"
            )
        else:
            code += ")\n"

        code += "record = SeqIO.read(handle, 'genbank')"

        return _ps(code)
