#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""Provides two functions, parse and parse_primers"""

# import os as _os
import re as _re
import io as _io
import textwrap as _textwrap

from Bio import SeqIO as _SeqIO
from pydna.genbankfile import GenbankFile as _GenbankFile
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from pydna.primer import Primer as _Primer

# from pydna.amplify import pcr as _pcr
# from copy import deepcopy as _deepcopy
# from Bio.SeqFeature import SeqFeature as _SeqFeature
# import xml.etree.ElementTree as _et

try:
    from itertools import pairwise as _pairwise
except ImportError:

    def _pairwise(iterable):
        # pairwise('ABCDEFG') → AB BC CD DE EF FG
        iterator = iter(iterable)
        a = next(iterator, None)
        for b in iterator:
            yield a, b
            a = b


# "^>.+?^(?=$|LOCUS|ID|>|\#)|^(?:LOCUS|ID).+?^//"
# "(?:^>.+\n^(?:^[^>]+?)(?=\n\n|>|^LOCUS|ID))|(?:(?:^LOCUS|ID)(?:(?:.|\n)+?)^//)"

# gb_fasta_embl_regex = r"(?:>.+\n^(?:^[^>]+?)(?=\n\n|>|LOCUS|ID))|(?:(?:LOCUS|ID)(?:(?:.|\n)+?)^//)"

gb_fasta_embl_regex = (
    r"(?:^>.+\n^(?:^[^>]+?)(?=\n\n|>|^LOCUS|^ID))|(?:(?:^LOCUS|^ID)(?:(?:.|\n)+?)^//)"
)

# The gb_fasta_embl_regex is meant to be able to extract sequences from
# text where sequences are mixed with other contents as well
# use https://regex101.com to get an idea how it works.


def extract_from_text(text):
    """docstring."""
    data = _textwrap.dedent(str(text))
    mos = list(_re.finditer(gb_fasta_embl_regex, data + "\n\n", flags=_re.MULTILINE))

    class Fakemo(object):
        def start(self):
            return len(data)

        def end(self):
            return 0

    mofirst = molast = Fakemo()

    gaps = []

    for mo1, mo2 in _pairwise([mofirst] + mos + [molast]):
        gaps.append(data[mo1.end() : mo2.start()])

    return tuple(mo.group(0) for mo in mos), tuple(gaps)


def embl_gb_fasta(text):
    """Parse embl, genbank or fasta format from text.

    Returns list of Bio.SeqRecord.SeqRecord

    annotations["molecule_type"]
    annotations["topology"]

    """
    chunks, gaps = extract_from_text(text)
    result_list = []
    # topology = "linear"

    for chunk in chunks:
        handle = _io.StringIO(chunk)
        # circular = False
        first_line = chunk.splitlines()[0].lower().split()
        try:
            parsed = _SeqIO.read(handle, "embl")
        except ValueError:
            handle.seek(0)
            try:
                parsed = _SeqIO.read(handle, "genbank")
            except ValueError:
                handle.seek(0)
                try:
                    parsed = _SeqIO.read(handle, "fasta-blast")
                except ValueError:
                    handle.close()
                    continue
                else:
                    # hack to pick up molecule_type from FASTA header line
                    if "protein" in first_line:
                        parsed.annotations["molecule_type"] = "protein"
                        parsed.annotations["topology"] = "linear"
                    else:
                        parsed.annotations["molecule_type"] = "DNA"
            # else:
            #     if _re.match(r"LOCUS\s+(\S+)\s+(\S+)\s+aa", " ".join(first_line)):
            #         parsed.annotations["molecule_type"] = "protein"
            #         parsed.annotations["topology"] = "linear"
        handle.close()
        # hack to pick up topology from FASTA and malformed gb files
        first_line = chunk.splitlines()[0].lower().split()
        parsed.annotations["topology"] = "linear"
        if "circular" in first_line:
            parsed.annotations["topology"] = "circular"
        # assert parsed.annotations.get("topology"), "topology must be set"
        assert parsed.annotations.get("molecule_type"), "molecule_type  must be set"
        if not parsed.annotations.get("molecule_type"):
            print(parsed)
        result_list.append(parsed)
    return tuple(result_list)


def parse(data, ds=True):
    """Return *all* DNA sequences found in data.

    If no sequences are found, an empty list is returned. This is a greedy
    function, use carefully.

    Parameters
    ----------
    data : string or iterable
        The data parameter is a string containing:

        1. an absolute path to a local file.
           The file will be read in text
           mode and parsed for EMBL, FASTA
           and Genbank sequences. Can be
           a string or a Path object.

        2. a string containing one or more
           sequences in EMBL, GENBANK,
           or FASTA format. Mixed formats
           are allowed.

        3. data can be a list or other iterable where the elements are 1 or 2

    ds : bool
        If True double stranded :class:`Dseqrecord` objects are returned.
        If False single stranded :class:`Bio.SeqRecord` [#]_ objects are
        returned.

    Returns
    -------
    list
        contains Dseqrecord or SeqRecord objects

    References
    ----------
    .. [#] http://biopython.org/wiki/SeqRecord

    See Also
    --------
    read

    """

    # a string is an iterable datatype but on Python2.x
    # it doesn't have an __iter__ method.
    if not hasattr(data, "__iter__") or isinstance(data, (str, bytes)):
        data = (data,)

    sequences = []

    for item in data:
        try:
            # item is a path to a utf-8 encoded text file?
            with open(item, "r", encoding="utf-8") as f:
                raw = f.read()
        except IOError:
            # item was not a path, add sequences parsed from item
            raw = item
            path = None
        else:
            # item was a readable text file, seqences are parsed from the file
            path = item
        finally:
            newsequences = embl_gb_fasta(raw)
            # nfs = [_SeqFeature() for f in parsed.features]
            # for f, nf in zip(parsed.features, nfs):
            #     nf.__dict__ = _deepcopy(f.__dict__)
            # parsed.features = nfs
            for s in newsequences:
                if ds and path:
                    sequences.append(_GenbankFile.from_SeqRecord(s, path=path))
                elif ds:
                    sequences.append(_Dseqrecord.from_SeqRecord(s))
                else:
                    sequences.append(s)
    return sequences


def parse_primers(data):
    """docstring."""
    return [_Primer(x) for x in parse(data, ds=False)]


def parse_snapgene(file_path: str) -> list[_Dseqrecord]:
    """Parse a SnapGene file and return a Dseqrecord object.

    Parameters
    ----------
    file_path : str
        The path to the SnapGene file to parse.

    Returns
    -------
    Dseqrecord
        The parsed SnapGene file as a Dseqrecord object.

    """
    with open(file_path, "rb") as f:
        parsed_seq = next(_SeqIO.parse(f, "snapgene"))
        circular = (
            "topology" in parsed_seq.annotations.keys()
            and parsed_seq.annotations["topology"] == "circular"
        )
        return [_Dseqrecord(parsed_seq, circular=circular)]
