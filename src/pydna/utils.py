#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Miscellaneous functions."""

from Bio.Data.IUPACData import ambiguous_dna_complement as _ambiguous_dna_complement
from Bio.Seq import _maketrans
import shelve as _shelve
import os as _os
import re as _re
import logging as _logging
import base64 as _base64
import pickle as _pickle
import hashlib as _hashlib
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

from typing import Union as _Union, TypeVar as _TypeVar, List as _List

# For functions that take str or bytes as input and return str or bytes as output, matching the input type
StrOrBytes = _TypeVar("StrOrBytes", str, bytes)

_module_logger = _logging.getLogger("pydna." + __name__)
_ambiguous_dna_complement.update({"U": "A"})
_complement_table = _maketrans(_ambiguous_dna_complement)


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
                    orfs.append((frame, startindex * 3 + frame, (stopindex + 1) * 3 + frame))
                # print(stopindex, startindex, limit)
    return orfs


def shift_location(original_location, shift, lim):
    """docstring."""
    newparts = []
    strand = original_location.strand
    if lim is None:
        if min(original_location) + shift < 0:
            raise ValueError("Shift moves location below zero, use a `lim` to loop around if sequence is circular.")
        lim = _sys.maxsize

    for part in original_location.parts:
        new_start = (part.start + shift) % lim
        new_end = (part.end + shift) % lim or lim
        old_start, old_end = (newparts[-1].start, newparts[-1].end) if len(newparts) else (None, None)

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


def memorize(filename):
    """Cache functions and classes.

    see pydna.download
    """

    def decorator(f):
        def wrappee(*args, **kwargs):
            _module_logger.info("#### memorizer ####")
            _module_logger.info("cache filename                   = %s", filename)
            _module_logger.info(
                "os.environ['pydna_cached_funcs'] = %s",
                _os.getenv("pydna_cached_funcs", ""),
            )
            if filename not in _os.getenv("pydna_cached_funcs", ""):
                _module_logger.info("cache filename not among cached functions, made it new!")
                return f(*args, **kwargs)
            key = _base64.urlsafe_b64encode(_hashlib.sha1(_pickle.dumps((args, kwargs))).digest()).decode("ascii")
            _module_logger.info("key = %s", key)
            cache = _shelve.open(
                _os.path.join(_os.environ["pydna_data_dir"], identifier_from_string(filename)),
                writeback=False,
            )
            try:
                result = cache[key]
            except KeyError:
                _module_logger.info(
                    "no result for key %s in shelve %s",
                    key,
                    identifier_from_string(filename),
                )
                result = f(*args, **kwargs)
                _module_logger.info("made it new!")
                cache[key] = result
                _module_logger.info("saved result under key %s", key)
            else:
                _module_logger.info("found %s in cache", key)
            cache.close()
            return result

        return wrappee

    return decorator


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

    return random.choice(starts) + "".join([random.choice(cdns) for x in range(length)]) + random.choice(stops)


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
        topology = set([arg.circular if hasattr(arg, "circular") else None for arg in args])

        if len(topology) != 1:
            raise ValueError("sequences have different topologies")
        topology = topology.pop()
        if topology in (False, None):
            topology = "linear"
        elif topology is True:
            topology = "circular"

    args = [arg.seq if hasattr(arg, "seq") else arg for arg in args]
    args_string_list = [arg.watson.lower() if hasattr(arg, "watson") else str(arg).lower() for arg in args]

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


if __name__ == "__main__":
    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
