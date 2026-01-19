#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_pydna_gbtext_clean():
    pytest.importorskip("pyparsing")

    from pydna.readers import read
    from pydna.genbankfixer import gbtext_clean

    files = [
        ("sequence.gb", "ldseguid=_GHnVotdYSxAKX2bAUIa10ZZDhE"),
        ("NCBI_example.gb", "ldseguid=_GHnVotdYSxAKX2bAUIa10ZZDhE"),
        ("YEplac181.txt", "ldseguid=ktUgiF4Ah2i6zM8DbWewnTuwmkw"),
        ("pGADT7-Rec.gb", "cdseguid=L0wFG4Mx6rWMYej1cNpCqSYGOd0"),
        (
            "P30350%20(2013-10-11%2013_49_14).dna.txt",
            "cdseguid=cr9QvaW0QQWgiRNmG3W-ebeCQK0",
        ),
        ("ApE_example.gb", "ldseguid=2oLvuMsAc-fgAISC7b3XJBuejGI"),
        ("VectorNTI_example.gb", "cdseguid=fLRW6jTMqqwZCbIFOGyl1TdDbrs"),
        ("hej.txt", "ldseguid=ktUgiF4Ah2i6zM8DbWewnTuwmkw"),
        ("fakeGenBankFile.gb", "cdseguid=cJQ-aN4c5sLYq7FcNLKr5ZVSdSM"),
    ]

    for file_, seg in files:
        with open("broken_genbank_files/" + file_, "r") as f:
            infile = f.read()
        if file_ == "hej.txt":
            from Bio import BiopythonParserWarning

            with pytest.warns(BiopythonParserWarning):
                assert read(gbtext_clean(infile).gbtext).seguid() == seg
        else:
            assert read(gbtext_clean(infile).gbtext).seguid() == seg


def test_wrapstring():
    pytest.importorskip("pyparsing")

    from pydna.genbankfixer import wrapstring

    assert (
        wrapstring(
            "0123456789",
            0,
            5,
        )
        == "01234\n56789\n\n"
    )
    assert wrapstring("0123456789", 1, 5, padfirst=False) == "0123\n 4567\n 89\n"


if __name__ == "__main__":
    pytest.main([__file__, "-x", "-vvv", "-s"])
