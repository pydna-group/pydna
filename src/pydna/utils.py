#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Miscellaneous functions."""

from Bio.Data.IUPACData import ambiguous_dna_complement as _ambiguous_dna_complement
from Bio.Seq import _maketrans

# import shelve as _shelve
# import os as _os
import re as _re

# import logging as _logging
# import base64 as _base64
# import pickle as _pickle
# import hashlib as _hashlib
import keyword as _keyword
import collections as _collections
import itertools as _itertools
from copy import deepcopy as _deepcopy

import sys as _sys
import random
import subprocess as _subprocess
from bisect import bisect as _bisect
from math import ceil as _ceil

from pydna.codon import weights as _weights
from pydna.codon import rare_codons as _rare_codons

from Bio.SeqFeature import SimpleLocation as _sl
from Bio.SeqFeature import CompoundLocation as _cl
from Bio.SeqFeature import Location as _Location

from typing import Union as _Union, TypeVar as _TypeVar, List as _List

# For functions that take str or bytes as input and return str or bytes as output, matching the input type
StrOrBytes = _TypeVar("StrOrBytes", str, bytes)

# _module_logger = _logging.getLogger("pydna." + __name__)
# _ambiguous_dna_complement.update({"U": "A"})
_ambiguous_dna_complement.update(
    {  # dsIUPAC
        "P": "J",  # G in top strand, complementary strand empty.
        "E": "Z",  # A "
        "X": "F",  # T "
        "I": "Q",  # C "
        "U": "O",  # U in top strand, A in complementary strand.
        "O": "U",  # A in top strand, U in complementary strand.
        "J": "P",  # top strand empty, G in complementary strand.
        "Z": "E",  # "                 A "
        "F": "X",  # "                 T "
        "Q": "I",  # "                 C "
    }
)

_keys = "".join(_ambiguous_dna_complement.keys()).encode("ASCII")
_values = "".join(_ambiguous_dna_complement.values()).encode("ASCII")
_complement_table = bytes.maketrans(_keys + _keys.lower(), _values + _values.lower())
_complement_table = _maketrans(_ambiguous_dna_complement)

to_watson_table = bytes.maketrans(
    b"PEXIQFZJpexiqfzj" b"12" b"UuOoGATCgatc", b"GATC    gatc    " b" ." b"UuAaGATCgatc"
)

to_crick_table = bytes.maketrans(
    _keys + _keys.lower() + b"PEXIQFZJpexiqfzj" b"12" b"UuOoGATCgatc",
    _values + _values.lower() + b"    CTAG    ctag" b". " b"AaUuCTAGctag",
)

to_N = rd = bytes.maketrans(b"PEXIpexiQFZJqfzjUuOo", b"NNNNNNNNNNNNNNNNNNNN")
to_5tail_table = bytes.maketrans(b"GATCgatc", b"QFZJqfzj")
to_3tail_table = bytes.maketrans(b"GATCgatc", b"PEXIpexi")
to_full_sequence = bytes.maketrans(b"PEXIpexiQFZJqfzj", b"GATCgatcGATCgatc")


iupac_regex = {  # IUPAC Ambiguity Codes for Nucleotide Degeneracy and U for Uracile
    "A": "(?:A)",
    "C": "(?:C)",
    "G": "(?:G)",
    "T": "(?:T|U)",
    "U": "(?:U|T|U)",
    "R": "(?:A|G|R)",
    "Y": "(?:C|T|Y)",
    "S": "(?:G|C|S)",
    "W": "(?:A|T|W)",
    "K": "(?:G|T|K)",
    "M": "(?:A|C|M)",
    "B": "(?:C|G|T|B)",
    "D": "(?:A|G|T|D)",
    "H": "(?:A|C|T|H)",
    "V": "(?:A|C|G|V)",
    "N": "(?:A|G|C|T|N)",
}

iupac_compl_regex = {  # IUPAC Ambiguity Code complements
    "A": "(?:T|U)",
    "C": "(?:G)",
    "G": "(?:C)",
    "T": "(?:A)",
    "U": "(?:A)",
    "R": "(?:T|C|Y)",
    "Y": "(?:G|A|R)",
    "S": "(?:G|C|S)",
    "W": "(?:A|T|W)",
    "K": "(?:C|AM)",
    "M": "(?:T|G|K)",
    "B": "(?:C|G|A|V)",
    "D": "(?:A|C|T|H)",
    "H": "(?:A|G|T|D)",
    "V": "(?:T|C|G|B)",
    "N": "(?:A|G|C|T|N)",
}


bp_dict = {
    (b"P", b"Q"): b"G",  # P / Q  >>--->  G
    (b"E", b"F"): b"A",  # E / F  >>--->  A
    (b"X", b"Z"): b"T",  # X / Z  >>--->  T
    (b"I", b"J"): b"C",  # I / J  >>--->  C
    (b"p", b"q"): b"g",  # p / q  >>--->  g
    (b"e", b"f"): b"a",  # e / f  >>--->  a
    (b"x", b"z"): b"t",  # x / z  >>--->  t
    (b"i", b"j"): b"c",  # i / j  >>--->  c
    (b"p", b"Q"): b"g",  # p / Q  >>--->  g
    (b"e", b"F"): b"a",  # e / F  >>--->  a
    (b"x", b"Z"): b"t",  # x / Z  >>--->  t
    (b"i", b"J"): b"c",  # i / J  >>--->  c
    (b"P", b"q"): b"G",  # P / q  >>--->  G
    (b"E", b"f"): b"A",  # E / f  >>--->  A
    (b"X", b"z"): b"T",  # X / z  >>--->  T
    (b"I", b"j"): b"C",  # I / j  >>--->  C
    # (b"P", b" "): b"P",  # P / q  >>--->  G
    # (b"E", b" "): b"E",  # E / f  >>--->  A
    # (b"X", b" "): b"X",  # X / z  >>--->  T
    # (b"I", b" "): b"I",  # I / j  >>--->  C
    (b"Q", b"P"): b"G",  # Q / P  >>--->  G
    (b"F", b"E"): b"A",  # F / E  >>--->  A
    (b"Z", b"X"): b"T",  # Z / X  >>--->  T
    (b"J", b"I"): b"C",  # J / I  >>--->  C
    (b"q", b"p"): b"g",  # q / p  >>--->  g
    (b"f", b"e"): b"a",  # f / e  >>--->  a
    (b"z", b"x"): b"t",  # z / x  >>--->  t
    (b"j", b"i"): b"c",  # j / i  >>--->  c
    (b"q", b"P"): b"G",  # Q / P  >>--->  G
    (b"f", b"E"): b"A",  # F / E  >>--->  A
    (b"z", b"X"): b"T",  # Z / X  >>--->  T
    (b"j", b"I"): b"C",  # J / I  >>--->  C
    (b"Q", b"p"): b"g",  # q / p  >>--->  g
    (b"F", b"e"): b"a",  # f / e  >>--->  a
    (b"Z", b"x"): b"t",  # z / x  >>--->  t
    (b"J", b"i"): b"c",  # j / i  >>--->  c
    # (b" ", b"Q"): b"Q",  # Q / P  >>--->  G
    # (b" ", b"F"): b"F",  # F / E  >>--->  A
    # (b" ", b"Z"): b"Z",  # Z / X  >>--->  T
    # (b" ", b"J"): b"J",  # J / I  >>--->  C
    (b"G", b" "): b"P",
    (b"A", b" "): b"E",
    (b"T", b" "): b"X",
    (b"C", b" "): b"I",
    (b"g", b" "): b"p",
    (b"a", b" "): b"e",
    (b"t", b" "): b"x",
    (b"c", b" "): b"i",
    (b" ", b"G"): b"J",
    (b" ", b"A"): b"Z",
    (b" ", b"T"): b"F",
    (b" ", b"C"): b"Q",
    (b" ", b"g"): b"j",
    (b" ", b"a"): b"z",
    (b" ", b"t"): b"f",
    (b" ", b"c"): b"q",
    (b"G", b"C"): b"G",
    (b"A", b"T"): b"A",
    (b"T", b"A"): b"T",
    (b"C", b"G"): b"C",
    (b"g", b"c"): b"g",
    (b"a", b"t"): b"a",
    (b"t", b"a"): b"t",
    (b"c", b"g"): b"c",
    (b"G", b"c"): b"G",
    (b"A", b"t"): b"A",
    (b"T", b"a"): b"T",
    (b"C", b"g"): b"C",
    (b"g", b"C"): b"g",
    (b"a", b"T"): b"a",
    (b"t", b"A"): b"t",
    (b"c", b"G"): b"c",
    (b"U", b"O"): b"U",
    (b"O", b"U"): b"A",
    (b"u", b"o"): b"u",
    (b"o", b"u"): b"a",
    (b"u", b"O"): b"u",
    (b"O", b"u"): b"A",
    (b"U", b"o"): b"U",
    (b"o", b"U"): b"a",
    (b"U", b"A"): b"U",
    (b"A", b"U"): b"O",
    (b"u", b"a"): b"u",
    (b"a", b"u"): b"o",
    (b"u", b"A"): b"u",
    (b"A", b"u"): b"O",
    (b"U", b"a"): b"U",
    (b"a", b"U"): b"o",
    (b"M", b"T"): b"A",
    (b"m", b"t"): b"a",
    (b"M", b"t"): b"A",
    (b"m", b"T"): b"a",
    (b"T", b"M"): b"K",
    (b"t", b"m"): b"k",
    (b"T", b"m"): b"K",
    (b"t", b"M"): b"k",
    (b"M", b"G"): b"C",
    (b"m", b"g"): b"c",
    (b"M", b"g"): b"C",
    (b"m", b"G"): b"c",
    (b"G", b"M"): b"K",
    (b"g", b"m"): b"k",
    (b"G", b"m"): b"K",
    (b"g", b"M"): b"k",
    (b"R", b"T"): b"A",
    (b"r", b"t"): b"a",
    (b"R", b"t"): b"A",
    (b"r", b"T"): b"a",
    (b"T", b"R"): b"Y",
    (b"t", b"r"): b"y",
    (b"T", b"r"): b"Y",
    (b"t", b"R"): b"y",
    (b"R", b"C"): b"G",
    (b"r", b"c"): b"g",
    (b"R", b"c"): b"G",
    (b"r", b"C"): b"g",
    (b"C", b"R"): b"Y",
    (b"c", b"r"): b"y",
    (b"C", b"r"): b"Y",
    (b"c", b"R"): b"y",
    (b"W", b"T"): b"A",
    (b"w", b"t"): b"a",
    (b"W", b"t"): b"A",
    (b"w", b"T"): b"a",
    (b"T", b"W"): b"W",
    (b"t", b"w"): b"w",
    (b"T", b"w"): b"W",
    (b"t", b"W"): b"w",
    (b"W", b"A"): b"T",
    (b"w", b"a"): b"t",
    (b"W", b"a"): b"T",
    (b"w", b"A"): b"t",
    (b"A", b"W"): b"W",
    (b"a", b"w"): b"w",
    (b"A", b"w"): b"W",
    (b"a", b"W"): b"w",
    (b"S", b"G"): b"C",
    (b"s", b"g"): b"c",
    (b"S", b"g"): b"C",
    (b"s", b"G"): b"c",
    (b"G", b"S"): b"S",
    (b"g", b"s"): b"s",
    (b"G", b"s"): b"S",
    (b"g", b"S"): b"s",
    (b"S", b"C"): b"G",
    (b"s", b"c"): b"g",
    (b"S", b"c"): b"G",
    (b"s", b"C"): b"g",
    (b"C", b"S"): b"S",
    (b"c", b"s"): b"s",
    (b"C", b"s"): b"S",
    (b"c", b"S"): b"s",
    (b"Y", b"G"): b"C",
    (b"y", b"g"): b"c",
    (b"Y", b"g"): b"C",
    (b"y", b"G"): b"c",
    (b"G", b"Y"): b"R",
    (b"g", b"y"): b"r",
    (b"G", b"y"): b"R",
    (b"g", b"Y"): b"r",
    (b"Y", b"A"): b"T",
    (b"y", b"a"): b"t",
    (b"Y", b"a"): b"T",
    (b"y", b"A"): b"t",
    (b"A", b"Y"): b"R",
    (b"a", b"y"): b"r",
    (b"A", b"y"): b"R",
    (b"a", b"Y"): b"r",
    (b"K", b"C"): b"G",
    (b"k", b"c"): b"g",
    (b"K", b"c"): b"G",
    (b"k", b"C"): b"g",
    (b"C", b"K"): b"M",
    (b"c", b"k"): b"m",
    (b"C", b"k"): b"M",
    (b"c", b"K"): b"m",
    (b"K", b"A"): b"T",
    (b"k", b"a"): b"t",
    (b"K", b"a"): b"T",
    (b"k", b"A"): b"t",
    (b"A", b"K"): b"M",
    (b"a", b"k"): b"m",
    (b"A", b"k"): b"M",
    (b"a", b"K"): b"m",
    (b"V", b"T"): b"A",
    (b"v", b"t"): b"a",
    (b"V", b"t"): b"A",
    (b"v", b"T"): b"a",
    (b"T", b"V"): b"B",
    (b"t", b"v"): b"b",
    (b"T", b"v"): b"B",
    (b"t", b"V"): b"b",
    (b"V", b"G"): b"C",
    (b"v", b"g"): b"c",
    (b"V", b"g"): b"C",
    (b"v", b"G"): b"c",
    (b"G", b"V"): b"B",
    (b"g", b"v"): b"b",
    (b"G", b"v"): b"B",
    (b"g", b"V"): b"b",
    (b"V", b"C"): b"G",
    (b"v", b"c"): b"g",
    (b"V", b"c"): b"G",
    (b"v", b"C"): b"g",
    (b"C", b"V"): b"B",
    (b"c", b"v"): b"b",
    (b"C", b"v"): b"B",
    (b"c", b"V"): b"b",
    (b"H", b"T"): b"A",
    (b"h", b"t"): b"a",
    (b"H", b"t"): b"A",
    (b"h", b"T"): b"a",
    (b"T", b"H"): b"D",
    (b"t", b"h"): b"d",
    (b"T", b"h"): b"D",
    (b"t", b"H"): b"d",
    (b"H", b"G"): b"C",
    (b"h", b"g"): b"c",
    (b"H", b"g"): b"C",
    (b"h", b"G"): b"c",
    (b"G", b"H"): b"D",
    (b"g", b"h"): b"d",
    (b"G", b"h"): b"D",
    (b"g", b"H"): b"d",
    (b"H", b"A"): b"T",
    (b"h", b"a"): b"t",
    (b"H", b"a"): b"T",
    (b"h", b"A"): b"t",
    (b"A", b"H"): b"D",
    (b"a", b"h"): b"d",
    (b"A", b"h"): b"D",
    (b"a", b"H"): b"d",
    (b"D", b"T"): b"A",
    (b"d", b"t"): b"a",
    (b"D", b"t"): b"A",
    (b"d", b"T"): b"a",
    (b"T", b"D"): b"H",
    (b"t", b"d"): b"h",
    (b"T", b"d"): b"H",
    (b"t", b"D"): b"h",
    (b"D", b"C"): b"G",
    (b"d", b"c"): b"g",
    (b"D", b"c"): b"G",
    (b"d", b"C"): b"g",
    (b"C", b"D"): b"H",
    (b"c", b"d"): b"h",
    (b"C", b"d"): b"H",
    (b"c", b"D"): b"h",
    (b"D", b"A"): b"T",
    (b"d", b"a"): b"t",
    (b"D", b"a"): b"T",
    (b"d", b"A"): b"t",
    (b"A", b"D"): b"H",
    (b"a", b"d"): b"h",
    (b"A", b"d"): b"H",
    (b"a", b"D"): b"h",
    (b"B", b"G"): b"C",
    (b"b", b"g"): b"c",
    (b"B", b"g"): b"C",
    (b"b", b"G"): b"c",
    (b"G", b"B"): b"V",
    (b"g", b"b"): b"v",
    (b"G", b"b"): b"V",
    (b"g", b"B"): b"v",
    (b"B", b"C"): b"G",
    (b"b", b"c"): b"g",
    (b"B", b"c"): b"G",
    (b"b", b"C"): b"g",
    (b"C", b"B"): b"V",
    (b"c", b"b"): b"v",
    (b"C", b"b"): b"V",
    (b"c", b"B"): b"v",
    (b"B", b"A"): b"T",
    (b"b", b"a"): b"t",
    (b"B", b"a"): b"T",
    (b"b", b"A"): b"t",
    (b"A", b"B"): b"V",
    (b"a", b"b"): b"v",
    (b"A", b"b"): b"V",
    (b"a", b"B"): b"v",
    # (b"N", b"C"): b"G",
    # (b"n", b"c"): b"g",
    # (b"N", b"c"): b"G",
    # (b"n", b"C"): b"g",
    # (b"C", b"N"): b"N",
    # (b"c", b"n"): b"n",
    # (b"C", b"n"): b"N",
    # (b"c", b"N"): b"n",
    # (b"N", b"T"): b"A",
    # (b"n", b"t"): b"a",
    # (b"N", b"t"): b"A",
    # (b"n", b"T"): b"a",
    # (b"T", b"N"): b"N",
    # (b"t", b"n"): b"n",
    # (b"T", b"n"): b"N",
    # (b"t", b"N"): b"n",
    # (b"N", b"A"): b"T",
    # (b"n", b"a"): b"t",
    # (b"N", b"a"): b"T",
    # (b"n", b"A"): b"t",
    # (b"A", b"N"): b"N",
    # (b"a", b"n"): b"n",
    # (b"A", b"n"): b"N",
    # (b"a", b"N"): b"n",
    # (b"N", b"G"): b"C",
    # (b"n", b"g"): b"c",
    # (b"N", b"g"): b"C",
    # (b"n", b"G"): b"c",
    # (b"G", b"N"): b"N",
    # (b"g", b"n"): b"n",
    # (b"G", b"n"): b"N",
    # (b"g", b"N"): b"n",
    (b"N", b"N"): b"N",
    (b"n", b"N"): b"n",
    (b"N", b"n"): b"N",
    (b"n", b"n"): b"n",
    # (b" ", b"N"): b"N",
    # (b" ", b"n"): b"n",
}

bp_dict_str = {
    (x.decode("ascii"), y.decode("ascii")): z.decode("ascii")
    for (x, y), z in bp_dict.items()
}


def three_frame_orfs(
    dna: str,
    limit: int = 100,
    startcodons: tuple = ("ATG",),
    stopcodons: tuple = ("TAG", "TAA", "TGA"),
    # startcodons: tuple[str, ...] = ("ATG",),
    # stopcodons: tuple[str, ...] = ("TAG", "TAA", "TGA"),
):
    """Overlapping orfs in three frames."""
    # breakpoint()
    limit = _ceil(limit / 3) - 1
    dna = dna.upper()

    orfs = []

    for frame in (0, 1, 2):

        codons = [dna[i : i + 3] for i in range(frame, len(dna), 3)]

        startdindices = [i for i, cd in enumerate(codons) if cd in startcodons]
        stopdindices = [i for i, cd in enumerate(codons) if cd in stopcodons]

        for startindex in startdindices:
            try:
                stopindex = stopdindices[_bisect(stopdindices, startindex)]
            except IndexError:
                pass
            else:
                if stopindex - startindex >= limit:
                    orfs.append(
                        (frame, startindex * 3 + frame, (stopindex + 1) * 3 + frame)
                    )
                # print(stopindex, startindex, limit)
    return orfs


def shift_location(original_location, shift, lim):
    """docstring."""
    newparts = []
    strand = original_location.strand
    if lim is None:
        if min(original_location) + shift < 0:
            raise ValueError(
                "Shift moves location below zero, use a `lim` to loop around if sequence is circular."
            )
        lim = _sys.maxsize

    for part in original_location.parts:
        new_start = (part.start + shift) % lim
        new_end = (part.end + shift) % lim or lim
        old_start, old_end = (
            (newparts[-1].start, newparts[-1].end) if len(newparts) else (None, None)
        )

        # The "join with old" cases are for features with multiple parts
        # in which consecutive parts do not have any bases between them.
        # This type of feature is generated to represent a feature that
        # spans the origin of a circular sequence. See more details in
        # https://github.com/pydna-group/pydna/issues/195

        if len(part) == 0:
            newparts.append(_sl(new_start, new_start, strand))
            continue
        # Join with old, case 1
        elif strand != -1 and old_end == new_start:
            part = newparts.pop()
            part._end = new_end
            new_start = part.start
        # Join with old, case 2
        elif strand == -1 and old_start == new_end:
            part = newparts.pop()
            part._start = new_start
            new_end = part.end
        if new_start < new_end:
            newparts.append(_sl(new_start, new_end, strand))
        else:
            parttuple = (_sl(new_start, lim, strand), _sl(0, new_end, strand))
            newparts.extend(parttuple if strand != -1 else parttuple[::-1])
    try:
        newloc = _cl(newparts)
    except ValueError:
        newloc, *n = newparts
    assert len(newloc) == len(original_location)
    return newloc


# def shift_feature(feature, shift, lim):
#     """Return a new feature with shifted location."""
#     # TODO: Missing tests
#     new_location = shift_location(feature.location, shift, lim)
#     new_feature = _deepcopy(feature)
#     new_feature.location = new_location
#     return new_feature


def shift_feature(feature, shift, lim):
    """Return a new feature with shifted location."""
    # TODO: Missing tests
    new_location = shift_location(feature.location, shift, lim)
    new_feature = _deepcopy(feature)
    new_feature.location = new_location
    return new_feature


# def smallest_rotation(s):
#     """Smallest rotation of a string.

#     Algorithm described in Pierre Duval, Jean. 1983. Factorizing Words
#     over an Ordered Alphabet. Journal of Algorithms & Computational Technology
#     4 (4) (December 1): 363–381. and Algorithms on strings and sequences based
#     on Lyndon words, David Eppstein 2011.
#     https://gist.github.com/dvberkel/1950267

#     Examples
#     --------
#     >>> from pydna.utils import smallest_rotation
#     >>> smallest_rotation("taaa")
#     'aaat'

#     """
#     prev, rep = None, 0
#     ds = _array("u", 2 * s)
#     lens, lends = len(s), len(ds)
#     old = 0
#     k = 0
#     w = ""
#     while k < lends:
#         i, j = k, k + 1
#         while j < lends and ds[i] <= ds[j]:
#             i = (ds[i] == ds[j]) and i + 1 or k
#             j += 1
#         while k < i + 1:
#             k += j - i
#             prev = w
#             w = ds[old:k]
#             old = k
#             if w == prev:
#                 rep += 1
#             else:
#                 prev, rep = w, 1
#             if len(w) * rep == lens:
#                 return "".join(w * rep)


def smallest_rotation(s):
    """Smallest rotation of a string.

    Algorithm described in Pierre Duval, Jean. 1983. Factorizing Words
    over an Ordered Alphabet. Journal of Algorithms & Computational Technology
    4 (4) (December 1): 363–381. and Algorithms on strings and sequences based
    on Lyndon words, David Eppstein 2011.
    https://gist.github.com/dvberkel/1950267

    Examples
    --------
    >>> from pydna.utils import smallest_rotation
    >>> smallest_rotation("taaa")
    'aaat'
    """
    from pydivsufsort import min_rotation

    k = min_rotation(bytes(s, "ascii"))
    return s[k:] + s[:k]


def anneal_from_left(watson: str, crick: str) -> int:
    """
    The length of the common prefix shared by two strings.

    Args:
        str1 (str): The first string.
        str2 (str): The second string.

    Returns:
        int: The length of the common prefix.
    """

    result = len(
        list(
            _itertools.takewhile(
                lambda x: bp_dict_str.get((x[0], x[1])), zip(watson, crick[::-1])
            )
        )
    )

    return result


def cai(seq: str, organism: str = "sce", weights: dict = _weights):
    """docstring."""
    from cai2 import CAI as _CAI

    return round(_CAI(seq.upper(), weights=weights[organism]), 3)


def rarecodons(seq: str, organism="sce"):
    """docstring."""
    rare = _rare_codons[organism]
    s = seq.upper()
    slices = []
    for i in range(0, len(seq) // 3):
        x, y = i * 3, i * 3 + 3
        trip = s[x:y]
        if trip in rare:
            slices.append(slice(x, y, 1))
    return slices


def express(seq: str, organism="sce"):
    """docstring.

    **NOT IMPLEMENTED YET**
    """
    # x = _PrettyTable(["cds", "len", "cai", "gc", "sta", "stp", "n-end"] + _rare_codons[organism] + ["rare"])
    # val = []

    # val.append(f"{self._data.upper().decode('ASCII')[:3]}..." f"{self._data.upper().decode('ASCII')[-3:]}")
    # val.append(len(self) / 3)
    # val.append(cai(organism))
    # val.append(gc())
    # val.append(startcodon())
    # val.append(stopcodon())
    # val.append(_n_end[organism].get(_seq3(self[3:6].translate())))
    # s = self._data.upper().decode("ASCII")
    # trps = [s[i * 3 : i * 3 + 3] for i in range(0, len(s) // 3)]
    # tot = 0
    # for cdn in _rare_codons[organism]:
    #     cnt = trps.count(cdn)
    #     tot += cnt
    #     val.append(cnt)
    # val.append(round(tot / len(trps), 3))
    # x.add_row(val)
    # return x
    raise NotImplementedError


def open_folder(pth):
    """docstring."""
    if _sys.platform == "win32":
        _subprocess.run(["start", pth], shell=True)
    elif _sys.platform == "darwin":
        _subprocess.run(["open", pth])
    else:
        try:
            _subprocess.run(["xdg-open", pth])
        except OSError:
            return "no cache to open."


def rc(sequence: StrOrBytes) -> StrOrBytes:
    """Reverse complement.

    accepts mixed DNA/RNA
    """
    return sequence.translate(_complement_table)[::-1]


def complement(sequence: str):
    """Complement.

    accepts mixed DNA/RNA
    """
    return sequence.translate(_complement_table)


# def memorize(filename):
#     """Cache functions and classes.

#     see pydna.download
#     """

#     def decorator(f):
#         def wrappee(*args, **kwargs):
#             _module_logger.info("#### memorizer ####")
#             _module_logger.info("cache filename                   = %s", filename)
#             _module_logger.info(
#                 "os.environ['pydna_cached_funcs'] = %s",
#                 _os.getenv("pydna_cached_funcs", ""),
#             )
#             if filename not in _os.getenv("pydna_cached_funcs", ""):
#                 _module_logger.info("cache filename not among cached functions, made it new!")
#                 return f(*args, **kwargs)
#             key = _base64.urlsafe_b64encode(_hashlib.sha1(_pickle.dumps((args, kwargs))).digest()).decode("ascii")
#             _module_logger.info("key = %s", key)
#             cache = _shelve.open(
#                 _os.path.join(_os.environ["pydna_data_dir"], identifier_from_string(filename)),
#                 writeback=False,
#             )
#             try:
#                 result = cache[key]
#             except KeyError:
#                 _module_logger.info(
#                     "no result for key %s in shelve %s",
#                     key,
#                     identifier_from_string(filename),
#                 )
#                 result = f(*args, **kwargs)
#                 _module_logger.info("made it new!")
#                 cache[key] = result
#                 _module_logger.info("saved result under key %s", key)
#             else:
#                 _module_logger.info("found %s in cache", key)
#             cache.close()
#             return result

#         return wrappee

#     return decorator


def identifier_from_string(s: str) -> str:
    """Return a valid python identifier.

    based on the argument s or an empty string
    """
    s = s.strip()
    s = _re.sub(r"\s+", r"_", s)
    s.replace("-", "_")
    s = _re.sub("[^0-9a-zA-Z_]", "", s)
    if s and not s[0].isidentifier() or _keyword.iskeyword(s):
        s = "_{s}".format(s=s)
    assert s == "" or s.isidentifier()
    return s


def flatten(*args) -> _List:
    """Flattens an iterable of iterables.

    Down to str, bytes, bytearray or any of the pydna or Biopython seq objects
    """
    output = []
    args = list(args)
    while args:
        top = args.pop()
        if (
            isinstance(top, _collections.abc.Iterable)
            and not isinstance(top, (str, bytes, bytearray))
            and not hasattr(top, "reverse_complement")
        ):
            args.extend(top)
        else:
            output.append(top)
    return output[::-1]


def seq31(seq):
    """Turn a three letter code protein sequence into one with one letter code.

    The single input argument 'seq' should be a protein sequence using single
    letter codes, as a python string.

    This function returns the amino acid sequence as a string using the one
    letter amino acid codes. Output follows the IUPAC standard (including
    ambiguous characters B for "Asx", J for "Xle" and X for "Xaa", and also U
    for "Sel" and O for "Pyl") plus "Ter" for a terminator given as an
    asterisk.

    Any unknown
    character (including possible gap characters), is changed into 'Xaa'.

    Examples
    --------
    >>> from Bio.SeqUtils import seq3
    >>> seq3("MAIVMGRWKGAR*")
    'MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer'
    >>> from pydna.utils import seq31
    >>> seq31('MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer')
    'M  A  I  V  M  G  R  W  K  G  A  R  *'
    """
    threecode = {
        "Ala": "A",
        "Asx": "B",
        "Cys": "C",
        "Asp": "D",
        "Glu": "E",
        "Phe": "F",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Lys": "K",
        "Leu": "L",
        "Met": "M",
        "Asn": "N",
        "Pro": "P",
        "Gln": "Q",
        "Arg": "R",
        "Ser": "S",
        "Thr": "T",
        "Val": "V",
        "Trp": "W",
        "Tyr": "Y",
        "Glx": "Z",
        "Xaa": "X",
        "Ter": "*",
        "Sel": "U",
        "Pyl": "O",
        "Xle": "J",
    }

    nr_of_codons = int(len(seq) / 3)
    sequence = [seq[i * 3 : i * 3 + 3].title() for i in range(nr_of_codons)]
    padding = " " * 2
    return padding.join([threecode.get(aa, "X") for aa in sequence])


def randomRNA(length, maxlength=None):
    """docstring."""
    if maxlength and maxlength > length:
        length = int(round(random.triangular(length, maxlength)))
    return "".join([random.choice("GAUC") for x in range(length)])


def randomDNA(length, maxlength=None):
    """docstring."""
    if maxlength and maxlength > length:
        length = int(round(random.triangular(length, maxlength)))
    return "".join([random.choice("GATC") for x in range(length)])


def randomORF(length, maxlength=None):
    """docstring."""
    length -= 2
    if maxlength and maxlength > length:
        length = int(round(random.triangular(length, maxlength - 2)))

    cdns = (
        "TTT",
        "TTC",
        "TTA",
        "TTG",
        "TCT",
        "TCC",
        "TCA",
        "TCG",
        "TAT",
        "TAC",
        "TGT",
        "TGC",
        "TGG",
        "CTT",
        "CTC",
        "CTA",
        "CTG",
        "CCT",
        "CCC",
        "CCA",
        "CCG",
        "CAT",
        "CAC",
        "CAA",
        "CAG",
        "CGT",
        "CGC",
        "CGA",
        "CGG",
        "ATT",
        "ATC",
        "ATA",
        "ATG",
        "ACT",
        "ACC",
        "ACA",
        "ACG",
        "AAT",
        "AAC",
        "AAA",
        "AAG",
        "AGT",
        "AGC",
        "AGA",
        "AGG",
        "GTT",
        "GTC",
        "GTA",
        "GTG",
        "GCT",
        "GCC",
        "GCA",
        "GCG",
        "GAT",
        "GAC",
        "GAA",
        "GAG",
        "GGT",
        "GGC",
        "GGA",
        "GGG",
    )

    starts = ("ATG",)
    stops = ("TAA", "TAG", "TGA")

    return (
        random.choice(starts)
        + "".join([random.choice(cdns) for x in range(length)])
        + random.choice(stops)
    )


def randomprot(length, maxlength=None):
    """docstring."""
    if maxlength and maxlength > length:
        length = int(round(random.triangular(length, maxlength)))
    return "".join([random.choice("ACDEFGHIKLMNPQRSTVWY") for x in range(length)])


def eq(*args, **kwargs):
    """Compare two or more DNA sequences for equality.

    Compares two or more DNA sequences for equality i.e. if they
    represent the same double stranded DNA molecule.

    Parameters
    ----------
    args : iterable
        iterable containing sequences
        args can be strings, Biopython Seq or SeqRecord, Dseqrecord
        or dsDNA objects.
    circular : bool, optional
        Consider all molecules circular or linear
    linear : bool, optional
        Consider all molecules circular or linear

    Returns
    -------
    eq : bool
        Returns True or False

    Notes
    -----
    Compares two or more DNA sequences for equality i.e. if they
    represent the same DNA molecule.

    Two linear sequences are considiered equal if either:

    1. They have the same sequence (case insensitive)
    2. One sequence is the reverse complement of the other

    Two circular sequences are considered equal if they are circular
    permutations meaning that they have the same length and:

    1. One sequence can be found in the concatenation of the other sequence with itself.
    2. The reverse complement of one sequence can be found in the concatenation of the other sequence with itself.

    The topology for the comparison can be set using one of the keywords
    linear or circular to True or False.

    If circular or linear is not set, it will be deduced from the topology of
    each sequence for sequences that have a linear or circular attribute
    (like Dseq and Dseqrecord).

    Examples
    --------
    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.utils import eq
    >>> eq("aaa","AAA")
    True
    >>> eq("aaa","AAA","TTT")
    True
    >>> eq("aaa","AAA","TTT","tTt")
    True
    >>> eq("aaa","AAA","TTT","tTt", linear=True)
    True
    >>> eq("Taaa","aTaa", linear = True)
    False
    >>> eq("Taaa","aTaa", circular = True)
    True
    >>> a=Dseqrecord("Taaa")
    >>> b=Dseqrecord("aTaa")
    >>> eq(a,b)
    False
    >>> eq(a,b,circular=True)
    True
    >>> a=a.looped()
    >>> b=b.looped()
    >>> eq(a,b)
    True
    >>> eq(a,b,circular=False)
    False
    >>> eq(a,b,linear=True)
    False
    >>> eq(a,b,linear=False)
    True
    >>> eq("ggatcc","GGATCC")
    True
    >>> eq("ggatcca","GGATCCa")
    True
    >>> eq("ggatcca","tGGATCC")
    True
    """
    args = flatten(args)  # flatten

    topology = None

    if "linear" in kwargs:
        if kwargs["linear"] is True:
            topology = "linear"
        if kwargs["linear"] is False:
            topology = "circular"
    elif "circular" in kwargs:
        if kwargs["circular"] is True:
            topology = "circular"
        if kwargs["circular"] is False:
            topology = "linear"
    else:
        topology = set(
            [arg.circular if hasattr(arg, "circular") else None for arg in args]
        )

        if len(topology) != 1:
            raise ValueError("sequences have different topologies")
        topology = topology.pop()
        if topology in (False, None):
            topology = "linear"
        elif topology is True:
            topology = "circular"

    args = [arg.seq if hasattr(arg, "seq") else arg for arg in args]
    args_string_list = [
        arg.watson.lower() if hasattr(arg, "watson") else str(arg).lower()
        for arg in args
    ]

    length = set((len(s) for s in args_string_list))

    if len(length) != 1:
        return False
    same = True

    if topology == "circular":
        # force circular comparison of all given sequences
        for s1, s2 in _itertools.combinations(args_string_list, 2):
            if not (s1 in s2 + s2 or rc(s1) in s2 + s2):
                same = False
    elif topology == "linear":
        # force linear comparison of all given sequences
        for s1, s2 in _itertools.combinations(args_string_list, 2):
            if not (s1 == s2 or s1 == rc(s2)):
                same = False
    return same


# def cuts_overlap(left_cut, right_cut, seq_len):
#     # Special cases:
#     if left_cut is None or right_cut is None or left_cut == right_cut:
#         return False

#     # This block of code would not be necessary if the cuts were
#     # initially represented like this
#     (left_watson, left_ovhg), _ = left_cut
#     (right_watson, right_ovhg), _ = right_cut
#     # Position of the cut on the crick strands on the left and right
#     left_crick = left_watson - left_ovhg
#     right_crick = right_watson - right_ovhg
#     if left_crick >= seq_len:
#         left_crick -= seq_len
#         left_watson -= seq_len
#     if right_crick >= seq_len:
#         right_crick -= seq_len
#         right_watson -= seq_len

#     # Convert into ranges x and y and see if ranges overlap
#     x = sorted([left_watson, left_crick])
#     y = sorted([right_watson, right_crick])
#     return (x[1] > y[0]) != (y[1] < x[0])


# def location_boundaries(loc: _Union[_sl, _cl]):
#     if loc.strand == -1:
#         return loc.parts[-1].start, loc.parts[0].end
#     else:
#         return loc.parts[0].start, loc.parts[-1].end


def cuts_overlap(left_cut, right_cut, seq_len):
    # Special cases:
    if left_cut is None or right_cut is None or left_cut == right_cut:
        return False

    # This block of code would not be necessary if the cuts were
    # initially represented like this
    (left_watson, left_ovhg), _ = left_cut
    (right_watson, right_ovhg), _ = right_cut
    # Position of the cut on the crick strands on the left and right
    left_crick = left_watson - left_ovhg
    right_crick = right_watson - right_ovhg
    if left_crick >= seq_len:
        left_crick -= seq_len
        left_watson -= seq_len
    if right_crick >= seq_len:
        right_crick -= seq_len
        right_watson -= seq_len

    # Convert into ranges x and y and see if ranges overlap
    x = sorted([left_watson, left_crick])
    y = sorted([right_watson, right_crick])
    return (x[1] > y[0]) != (y[1] < x[0])


def location_boundaries(loc: _Union[_sl, _cl]):
    if loc.strand == -1:
        return loc.parts[-1].start, loc.parts[0].end
    else:
        return loc.parts[0].start, loc.parts[-1].end


def locations_overlap(loc1: _Union[_sl, _cl], loc2: _Union[_sl, _cl], seq_len):
    start1, end1 = location_boundaries(loc1)
    start2, end2 = location_boundaries(loc2)

    boundaries1 = [(start1, end1)]
    boundaries2 = [(start2, end2)]

    if start1 > end1:
        boundaries1 = [
            [start1, end1 + seq_len],
            [start1 - seq_len, end1],
        ]
    if start2 > end2:
        boundaries2 = [
            [start2, end2 + seq_len],
            [start2 - seq_len, end2],
        ]

    for b1, b2 in _itertools.product(boundaries1, boundaries2):
        if b1[0] < b2[1] and b1[1] > b2[0]:
            return True

    return False


def sum_is_sticky(
    three_prime_end: tuple[str, str],
    five_prime_end: tuple[str, str],
    partial: bool = False,
) -> int:
    """Return the overlap length if the 3' end of seq1 and 5' end of seq2 ends are sticky and compatible for ligation.
    Return 0 if they are not compatible."""
    type_seq1, sticky_seq1 = three_prime_end
    type_seq2, sticky_seq2 = five_prime_end

    if (
        "blunt" != type_seq2
        and type_seq2 == type_seq1
        and str(sticky_seq2) == str(rc(sticky_seq1))
    ):
        return len(sticky_seq1)

    if not partial:
        return 0

    if type_seq1 != type_seq2 or type_seq2 == "blunt":
        return 0
    elif type_seq2 == "5'":
        sticky_seq1 = str(rc(sticky_seq1))
    elif type_seq2 == "3'":
        sticky_seq2 = str(rc(sticky_seq2))

    ovhg_len = min(len(sticky_seq1), len(sticky_seq2))
    # [::-1] to try the longest overhangs first
    for i in range(1, ovhg_len + 1)[::-1]:
        if sticky_seq1[-i:] == sticky_seq2[:i]:
            return i
    else:
        return 0


def limit_iterator(iterator, limit):
    """
    Call the function with an iterator to raise an error if the number of items is greater than the limit.
    """
    for i, x in enumerate(iterator):
        if i >= limit:
            raise ValueError(f"Too many possible paths (more than {limit})")
        yield x


def create_location(
    start: int, end: int, lim: int, strand: int | None = None
) -> _Location:
    """
    Create a location object from a start and end position.
    If the end position is less than the start position, the location is circular. It handles negative positions.

    Parameters
    ----------
    start : int
        The start position of the location.
    end : int
        The end position of the location.
    lim : int
        The length of the sequence.
    strand : int, optional
        The strand of the location. None, 1 or -1.

    Returns
    -------
    location : Location
        The location object. Can be a SimpleLocation or a CompoundLocation if the feature spans the origin of
        a circular sequence.

    Examples
    --------
    >>> from pydna.utils import create_location
    >>> str(create_location(0, 5, 10,-1))
    '[0:5](-)'
    >>> str(create_location(0, 5, 10,+1))
    '[0:5](+)'
    >>> str(create_location(0, 5, 10))
    '[0:5]'
    >>> str(create_location(8, 2, 10))
    'join{[8:10], [0:2]}'
    >>> str(create_location(8, 2, 10,-1))
    'join{[0:2](-), [8:10](-)}'
    >>> str(create_location(-2, 2, 10))
    'join{[8:10], [0:2]}'

    Note this special case, 0 is the same as len(seq)
    >>> str(create_location(5, 0, 10))
    '[5:10]'

    Note the special case where if start and end are the same,
    the location spans the entire sequence (it's not empty).
    >>> str(create_location(5, 5, 10))
    'join{[5:10], [0:5]}'

    """
    while start < 0:
        start += lim
    while end < 0:
        end += lim
    if end > start:
        return _sl(start, end, strand)
    else:
        return shift_location(_sl(start, end + lim, strand), 0, lim)
