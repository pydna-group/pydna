#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""Provides two functions, parse and parse_primers"""

import re
import io
import textwrap
import os
from Bio import SeqIO
from Bio.SeqIO.InsdcIO import GenBankScanner, GenBankIterator
import warnings

from pydna.dseqrecord import Dseqrecord
from Bio.SeqRecord import SeqRecord
from pydna.opencloning_models import UploadedFileSource
from pydna.primer import Primer


try:
    from itertools import pairwise
except ImportError:

    def pairwise(iterable):
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
    data = textwrap.dedent(str(text))
    mos = list(re.finditer(gb_fasta_embl_regex, data + "\n\n", flags=re.MULTILINE))

    class Fakemo(object):
        def start(self):
            return len(data)

        def end(self):
            return 0

    mofirst = molast = Fakemo()

    gaps = []

    for mo1, mo2 in pairwise([mofirst] + mos + [molast]):
        gaps.append(data[mo1.end() : mo2.start()])

    return tuple(mo.group(0) for mo in mos), tuple(gaps)


class CustomGenBankScanner(GenBankScanner):
    """
    A custom GenBank scanner that parses the LOCUS line and extracts the name, size, molecule_type, topology and date
    from malformed GenBank files.

    For example, the following line:

    ```
        LOCUS       pKM265       4536 bp    DNA   circular  SYN        21-JUN-2013
    ```
    """

    def _feed_first_line(self, consumer, line):
        # A regex for
        m = re.match(
            r"(?i)LOCUS\s+(?P<name>\S+)\s+(?P<size>\d+ bp)\s+(?P<molecule_type>\S+)(?:\s+(?P<topology>circular|linear))?(?:\s+.+\s+)?(?P<date>\d+-\w+-\d+)?",
            line,
        )
        if m is None:  # pragma: no cover - I don't think this can happen, but JIC
            raise ValueError("LOCUS line cannot be parsed")
        name, size, molecule_type, topology, date = m.groups()

        consumer.locus(name)
        consumer.size(size[:-3])
        consumer.molecule_type(molecule_type)
        consumer.topology(topology.lower() if topology is not None else None)
        consumer.date(date)


class CustomGenBankIterator(GenBankIterator):

    def __init__(self, source):
        super(GenBankIterator, self).__init__(source, fmt="GenBank")
        self.records = CustomGenBankScanner(debug=0).parse_records(self.stream)


def parse_genbank(handle: io.StringIO) -> SeqRecord:
    """
    Equivalent to SeqIO.read(handle, "genbank") but with more permissive parsing for the LOCUS line. It also
    removes BASE_COUNT lines that can be misplaced and anyway are not stored in the SeqRecord.annotations.
    """

    filtered_lines = list()
    for line in handle:
        if not line.lstrip().startswith("BASE COUNT"):
            filtered_lines.append(line)
    filtered_handle = io.StringIO("".join(filtered_lines))

    try:
        # This may raise a ValueError such as "LOCUS line does not contain space at position X"
        # then use the more permissive CustomGenBankIterator
        parsed = SeqIO.read(filtered_handle, "genbank")

        # Sometimes the line is enough to parse the record, but not enough to parse the topology,
        # then raise an error starting with "LOCUS line does not contain" to trigger the
        # more permissive CustomGenBankIterator
        if "topology" not in parsed.annotations.keys():
            raise ValueError("LOCUS line does not contain topology")
        return parsed
    except ValueError as e:
        if "LOCUS line does not contain" not in str(e):
            raise e
        filtered_handle.seek(0)
        warnings.warn(
            "LOCUS line is wrongly formatted, we used a more permissive parser.",
            stacklevel=2,
        )
        return next(CustomGenBankIterator(filtered_handle))


def embl_gb_fasta(text):
    """Parse embl, genbank or fasta format from text.

    Returns list of Bio.SeqRecord.SeqRecord

    annotations["molecule_type"]
    annotations["topology"]

    """
    chunks, gaps = extract_from_text(text)
    result_list = []

    for chunk in chunks:
        handle = io.StringIO(chunk)
        first_line = chunk.splitlines()[0].lower().split()
        try:
            parsed = SeqIO.read(handle, "embl")
            parsed.annotations["pydna_parse_sequence_file_format"] = "embl"
        except ValueError:
            handle.seek(0)
            try:
                parsed = parse_genbank(handle)
                parsed.annotations["pydna_parse_sequence_file_format"] = "genbank"
            except ValueError:
                handle.seek(0)
                try:
                    parsed = SeqIO.read(handle, "fasta-blast")
                    parsed.annotations["pydna_parse_sequence_file_format"] = "fasta"
                except ValueError:
                    handle.close()
                    continue
                else:
                    # molecule_type is not set by the Biopython FASTA parser
                    parsed.annotations["molecule_type"] = "DNA"
        handle.close()
        # hack to pick up topology from FASTA and malformed gb files
        first_line = chunk.splitlines()[0].lower().split()
        parsed.annotations["topology"] = "linear"
        if "circular" in first_line:
            parsed.annotations["topology"] = "circular"
        molecule_type = parsed.annotations.get("molecule_type")
        assert molecule_type, "molecule_type must be set"
        assert molecule_type != "protein", "molecule_type can not be 'protein'"
        result_list.append(parsed)
    return tuple(result_list)


def parse(data, ds=True, is_path=None) -> list[Dseqrecord | SeqRecord]:
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

    is_path : bool, optional
        If True, the data is treated as a path to a file. If False, the data is treated as a string (e.g. FASTA file content).
        If None, both are tried.

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
        if is_path is None:
            path = item if os.path.isfile(item) else None
        elif is_path:
            path = item
        else:
            path = None

        if path:
            with open(path, "r", encoding="utf-8") as f:
                raw = f.read()
        else:
            raw = item

        newsequences = embl_gb_fasta(raw)
        for s in newsequences:
            if ds and path:
                from pydna.opencloning_models import UploadedFileSource

                result = Dseqrecord.from_SeqRecord(s)
                result.source = UploadedFileSource(
                    file_name=str(path),  # we use str to handle PosixPath
                    sequence_file_format=s.annotations[
                        "pydna_parse_sequence_file_format"
                    ],
                    index_in_file=0,
                )
                sequences.append(result)
            elif ds:
                sequences.append(Dseqrecord.from_SeqRecord(s))
            else:
                sequences.append(s)
    return sequences


def parse_primers(data):
    """docstring."""
    return [Primer(x) for x in parse(data, ds=False)]


def parse_snapgene(file_path: str) -> list[Dseqrecord]:
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
        parsed_seq = next(SeqIO.parse(f, "snapgene"))
        circular = (
            "topology" in parsed_seq.annotations.keys()
            and parsed_seq.annotations["topology"] == "circular"
        )

        source = UploadedFileSource(
            file_name=str(file_path),
            sequence_file_format="snapgene",
            index_in_file=0,
        )
        return [Dseqrecord(parsed_seq, circular=circular, source=source)]
