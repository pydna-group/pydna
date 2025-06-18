#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""This module provides a class for downloading sequences from genbank
called Genbank and an function that does the same thing called genbank.

The function can be used if the environmental variable **pydna_email** has
been set to a valid email address. The easiest way to do this permanantly is to edit the
`pydna.ini` file. See the documentation of :func:`pydna.open_config_folder`"""

# from pydna.utils import memorize as _memorize
from pydna.genbankrecord import GenbankRecord as _GenbankRecord
from pydna.readers import read as _read

from Bio import Entrez as _Entrez
from typing import Literal as _Literal, Optional as _Optional
import re as _re
import os as _os

# import logging as _logging

# _module_logger = _logging.getLogger("pydna." + __name__)


# TODO http://httpbin.org/ use for testing?


class Genbank:
    """Class to facilitate download from genbank. It is easier and
    quicker to use the :func:`pydna.genbank.genbank` function directly.

    Parameters
    ----------
    users_email : string
        Has to be a valid email address. You should always tell
        Genbanks who you are, so that they can contact you.

    Examples
    --------

    >>> from pydna.genbank import Genbank
    >>> gb=Genbank("bjornjobb@gmail.com")
    >>> rec = gb.nucleotide("LP002422.1")   # <- entry from genbank
    >>> print(len(rec))
    1
    """

    def __init__(
        self,
        users_email: str,
        *,
        tool: str = "pydna",
    ) -> None:
        if not _re.match(
            r"[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}", users_email, _re.IGNORECASE
        ):
            raise ValueError("email address {} is not valid.".format(users_email))

        # _module_logger.info("#### Genbank ititiation ####")
        # _module_logger.info("Genbank initiated with email: %s", users_email)
        # _module_logger.info("Genbank initiated with tool : %s", tool)

        if users_email == "someone@example.com":
            raise ValueError(
                "you have to set your email address in order to download from Genbank"
            )
        self.email = users_email
        self.tool = tool

    def __repr__(self) -> str:
        """This method returns a short representation containing the email used to initiate."""
        return "GenbankConnection({})".format(self.email)

    # @_memorize("pydna.genbank.Genbank.nucleotide")
    def nucleotide(
        self,
        item: str,
        seq_start: _Optional[int] = None,
        seq_stop: _Optional[int] = None,
        strand: _Literal[1, 2] = 1,
    ) -> _GenbankRecord:
        """This method downloads a genbank nuclotide record from genbank. This method is
        cached by default. This can be controlled by editing the **pydna_cached_funcs** environment
        variable. The best way to do this permanently is to edit the edit the
        `pydna.ini` file. See the documentation of :func:`pydna.open_config_folder`

        Item is a string containing one genbank accession number
        for a nucleotide file. Genbank nucleotide accession numbers have this format:

        | A12345   = 1 letter  + 5 numerals
        | AB123456 = 2 letters + 6 numerals

        The accession number is sometimes followed by a point and version number

        | BK006936.2

        Item can also contain optional interval information in the following formats:

        | BK006936.2 REGION: complement(613900..615202)
        | NM_005546 REGION: 1..100
        | NM_005546 REGION: complement(1..100)
        | 21614549:1-100
        | 21614549:c100-1
        | 21614549 1-100
        | 21614549 c100-1

        It is useful to set an interval for large genbank records to limit the download time.
        The items above containing interval information and can be obtained directly by
        looking up an entry in Genbank and setting the `Change region shown` on the
        upper right side of the page. The `ACCESSION` line of the displayed Genbank
        file will have the formatting shown.

        Alternatively, seq_start and seq_stop can be set explicitly to the sequence intervals to be
        downloaded.

        If strand is 2. "c", "C", "crick", "Crick", "antisense","Antisense",
        "2", 2, "-" or "-1", the antisense (Crick) strand is returned, otherwise
        the sense (Watson) strand is returned.

        Result is returned as a :class:`pydna.genbankrecord.GenbankRecord` object.

        References
        ----------

        .. [#]   http://www.dsimb.inserm.fr/~fuchs/M2BI/AnalSeq/Annexes/Sequences/Accession_Numbers.htm
        .. [#]   http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
        """
        matches = (
            (1, _re.search(r"(REGION:\s(?P<start>\d+)\.\.(?P<stop>\d+))", item)),
            (
                2,
                _re.search(
                    r"(REGION: complement\((?P<start>\d+)\.\.(?P<stop>\d+)\))", item
                ),
            ),
            (1, _re.search(r"(:|\s)(?P<start>\d+)-(?P<stop>\d+)", item)),
            (2, _re.search(r"(:|\s)c(?P<start>\d+)-(?P<stop>\d+)", item)),
        )

        for strand_, match in matches:
            if match:
                seq_start = match.group("start")
                seq_stop = match.group("stop")
                item = item[: match.start()]
                strand = strand_
                break

        if strand not in [1, 2]:
            try:
                strand = {"c": 2, "crick": 2, "antisense": 2, "2": 2, "-": 2, "-1": 2}[
                    strand.lower()
                ]
            except (KeyError, AttributeError):
                strand = 1

        # _module_logger.info("#### Genbank download ####")
        # _module_logger.info("item  %s", item)
        # _module_logger.info("start %s", seq_start)
        # _module_logger.info("stop  %s", seq_stop)

        # _module_logger.info("strand  %s", str(strand))

        _Entrez.email = self.email
        _Entrez.tool = self.tool

        # _module_logger.info("Entrez.email  %s", self.email)
        text = _Entrez.efetch(
            db="nuccore",
            id=item,
            rettype="gbwithparts",
            seq_start=seq_start,
            seq_stop=seq_stop,
            strand=strand,
            retmode="text",
        ).read()

        # _module_logger.info("text[:160]  %s", text[:160])

        return _GenbankRecord(
            _read(text), item=item, start=seq_start, stop=seq_stop, strand=strand
        )


def genbank(accession: str = "CS570233.1", *args, **kwargs) -> _GenbankRecord:
    """
    Download a genbank nuclotide record.

    This function takes the same paramenters as the
    :func:pydna.genbank.Genbank.nucleotide method. The email address stored
    in the `pydna_email` environment variable is used. The easiest way set
    this permanantly is to edit the `pydna.ini` file.
    See the documentation of :func:`pydna.open_config_folder`

    if no accession is given, a `very short Genbank
    entry <https://www.ncbi.nlm.nih.gov/nuccore/CS570233.1>`_
    is used as an example (see below). This can be useful for testing the
    connection to Genbank.

    Please note that this result is also cached by default by settings in
    the pydna.ini file.
    See the documentation of :func:`pydna.open_config_folder`

    ::

        LOCUS       CS570233                  14 bp    DNA     linear   PAT 18-MAY-2007
        DEFINITION  Sequence 6 from Patent WO2007025016.
        ACCESSION   CS570233
        VERSION     CS570233.1
        KEYWORDS    .
        SOURCE      synthetic construct
          ORGANISM  synthetic construct
                    other sequences; artificial sequences.
        REFERENCE   1
          AUTHORS   Shaw,R.W. and Cottenoir,M.
          TITLE     Inhibition of metallo-beta-lactamase by double-stranded dna
          JOURNAL   Patent: WO 2007025016-A1 6 01-MAR-2007;
                    Texas Tech University System (US)
        FEATURES             Location/Qualifiers
             source          1..14
                             /organism="synthetic construct"
                             /mol_type="unassigned DNA"
                             /db_xref="taxon:32630"
                             /note="This is a 14bp aptamer inhibitor."
        ORIGIN
                1 atgttcctac atga
        //

    """
    email = _os.getenv("pydna_email")
    # _module_logger.info("#### genbank function called ####")
    # _module_logger.info("email      %s", email)
    # _module_logger.info("accession  %s", email)
    gb = Genbank(email)
    return gb.nucleotide(accession, *args, **kwargs)
