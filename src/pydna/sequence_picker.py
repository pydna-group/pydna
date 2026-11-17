#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2013-2026 Björn Johansson
# SPDX-FileCopyrightText: 2023-2026 The Project Contributors
# SPDX-License-Identifier: BSD-3-Clause

from pydna.dseqrecord import Dseqrecord
import os

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


email = os.getenv("pydna_email")
tool = "pydna"


def genbank_accession(s: str) -> Dseqrecord:
    """docstring."""
    s = Dseqrecord(s)

    NCBIWWW.email = email
    NCBIWWW.tool = tool

    result_handle = NCBIWWW.qblast(
        "blastn",
        "nt",
        str(s.seq),
        hitlist_size=1,
        alignments=1,
        descriptions=1,
        expect=1e-8,
        megablast=True,
        service="megablast",
        ungapped_alignment=True,
    )

    blast_records = NCBIXML.read(result_handle)
    best_alignment, *rest = blast_records.alignments
    best_hsp, *rest = best_alignment.hsps
    dbs = best_hsp.sbjct
    start, stop = sorted((best_hsp.sbjct_start, best_hsp.sbjct_end))
    result = Dseqrecord(
        dbs,
        circular=False,
        id=s.name,
        name=s.name,
        description=(f"{best_alignment.accession} REGION: {start}..{stop}"),
    )
    return result
