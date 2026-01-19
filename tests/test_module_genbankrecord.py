#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from textwrap import dedent


def test_genbankfile():
    from pydna import genbankrecord

    gbr = genbankrecord.GenbankRecord("aaa")

    assert (
        gbr.hyperlink
        == "<a href='https://www.ncbi.nlm.nih.gov/nuccore/accession?from=&to=&strand=1' target='_blank'>accession</a>"
    )

    assert repr(gbr) == "Gbnk(-3 accession)"

    assert (
        gbr._repr_html_()
        == "<a href='https://www.ncbi.nlm.nih.gov/nuccore/accession?from=&to=&strand=1' target='_blank'>accession</a>"
    )

    from unittest.mock import MagicMock

    pp = MagicMock()

    gbr._repr_pretty_(pp, None)

    pp.text.assert_called_with("Gbnk(-3 accession)")

    gbr = genbankrecord.GenbankRecord("aaa", start=1, stop=2)

    assert gbr.hyperlink == (
        "<a href='https://www.ncbi.nlm.nih.gov/nuccore/"
        "accession?from=1&to=2&strand=1'"
        " target='_blank'>accession 1-2</a>"
    )

    gbr_rc = gbr.rc()

    assert gbr_rc.strand == 2

    gbr2 = genbankrecord.GenbankRecord(
        "atgttcctacatg", item="CS570233.1", start=1, stop=13, strand=1
    )
    assert repr(gbr2) == "Gbnk(-13 CS570233.1 1-13)"

    assert gbr2.pydna_code() == dedent(
        "".join(
            [
                "from pydna.genbank import Genbank\n",
                "gb = Genbank('None')\n",
                "seq = gb.nucleotide('CS570233.1',\n",
                "                    seq_start=1,\n",
                "                    seq_stop=13,\n",
                "                    strand=1)",
            ]
        )
    )

    assert gbr2.biopython_code() == dedent(
        """\
    from Bio import Entrez, SeqIO
    Entrez.email = 'None'
    handle = Entrez.efetch(db='nuccore',
                           id='CS570233.1',
                           rettype='gbwithparts',
                           retmode='text',
                           seq_start=1,
                           seq_stop=13,
                           strand=1)
    record = SeqIO.read(handle, 'genbank')"""
    )

    from Bio.Seq import Seq
    from pydna.seqrecord import SeqRecord

    arg = SeqRecord(Seq("aaa"))

    genbankrecord.GenbankRecord.from_SeqRecord(arg)

    genbankrecord.GenbankRecord.from_SeqRecord(arg, start=0, stop=2)


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
