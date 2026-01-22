#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""This module provide most pydna functionality in the local namespace.

Example
-------

>>> from pydna.all import *
>>> Dseq("aaa")
Dseq(-3)
aaa
ttt
>>> Dseqrecord("aaa")
Dseqrecord(-3)
>>> from pydna.all import __all__
>>> __all__
['Anneal', 'pcr', 'Assembly', 'genbank', 'Genbank', 'Dseqrecord',
'Dseq', 'read', 'read_primer', 'parse', 'parse_primers', 'primer_design', 'assembly_fragments', 'eq', 'gbtext_clean']
>>>
"""


__all__ = [
    "Anneal",
    "pcr",
    "Assembly",
    "genbank",
    "Genbank",
    "Dseqrecord",
    "Dseq",
    "read",
    "read_primer",
    "parse",
    "parse_primers",
    "primer_design",
    "assembly_fragments",
    "eq",
    "gbtext_clean",
]


from pydna.amplify import Anneal
from pydna.amplify import pcr
from pydna.assembly import Assembly
from pydna.genbank import genbank
from pydna.genbank import Genbank
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydna.readers import read
from pydna.readers import read_primer
from pydna.parsers import parse
from pydna.parsers import parse_primers
from pydna.design import primer_design
from pydna.design import assembly_fragments
from pydna.utils import eq
from pydna.genbankfixer import gbtext_clean
