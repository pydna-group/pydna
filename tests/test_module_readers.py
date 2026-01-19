#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from pydna.readers import read
from pydna.readers import read_primer
from pydna.parsers import parse
from pathlib import Path
from Bio.SeqIO import read as BPread
from Bio.SeqIO import parse as BPparse


def test_read_no_data():
    data = ""
    with pytest.raises(ValueError) as err:
        read(data)
    assert err.match("No sequence found in data")
    with pytest.raises(ValueError):
        read(data, ds=False)
    assert err.match("No sequence found in data")
    with pytest.raises(ValueError):
        read(data, ds=True)
    assert err.match("No sequence found in data")


def test_read_too_many_sequences():
    data = ">a\na\n>g\ng"
    with pytest.raises(ValueError) as err:
        read(data)
    assert err.match("More than one sequence found in data")
    with pytest.raises(ValueError) as err:
        read(data, ds=False)
    assert err.match("More than one sequence found in data")
    with pytest.raises(ValueError) as err:
        read(data, ds=True)
    assert err.match("More than one sequence found in data")


def test_read_primer_no_data():
    data = ""
    with pytest.raises(ValueError):
        read_primer(data)
    with pytest.raises(ValueError):
        read_primer(data)
    with pytest.raises(ValueError):
        read_primer(data)


def test_read_primer():
    data = ">pr\ngatc"
    pr = read_primer(data)
    assert str(pr.seq) == "gatc"


def test_pydna_read_test():
    assert read("pydna_read_test.txt").format("gb")[349:368] == '/label="2micron 2µ"'


def test_parse_and_read_with_biopython_and_pydna():

    q = BPread("read1.gb", "gb")
    w = BPread("read2.gb", "gb")
    e = BPread("read3.fasta", "fasta")
    r = BPread("read4.fasta", "fasta")

    # a, b = BPparse("pth1.txt", "gb")
    with open("pth1.txt", "r", encoding="utf-8") as f:
        a, b = BPparse(f, "gb")

    # print("|" + a.features[13].qualifiers["label"][0] + "|")
    # print("|" + a.format("gb")[3314:3324] + "|")

    assert a.features[13].qualifiers["label"][0] == "2micron 2µ"
    assert a.format("gb")[3268:3278] == "2micron 2µ"

    x, y = parse("pth1.txt")

    assert "".join(a.format("gb").splitlines()[1:]) == "".join(
        x.format("gb").splitlines()[1:]
    )
    assert "".join(b.format("gb").strip().splitlines()[4:]) == "".join(
        y.format("gb").splitlines()[4:]
    )


def test_read_from_string():

    input_ = """
            LOCUS       New_DNA                    4 bp ds-DNA     linear       30-MAR-2013
            DEFINITION  .
            ACCESSION
            VERSION
            SOURCE      .
              ORGANISM  .
            COMMENT
            COMMENT     ApEinfo:methylated:1
            FEATURES             Location/Qualifiers
                 misc_feature    2..3
                                 /label=NewFeature
                                 /ApEinfo_fwdcolor=cyan
                                 /ApEinfo_revcolor=green
                                 /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                 width 5 offset 0
            ORIGIN
                    1 acgt
            //
            """
    a = read(input_)

    assert str(a.seq) == "ACGT"

    input_ = """>hej
               acgt"""

    assert str(a.seq) == "ACGT"

    input_ = """
            LOCUS       New_DNA                    4 bp ds-DNA     linear       30-MAR-2013
            DEFINITION  .
            ACCESSION
            VERSION
            SOURCE      .
              ORGANISM  .
            COMMENT
            COMMENT     ApEinfo:methylated:1
            FEATURES             Location/Qualifiers
                 misc_feature    2..3
                                 /label=NewFeature
                                 /ApEinfo_fwdcolor=cyan
                                 /ApEinfo_revcolor=green
                                 /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                 width 5 offset 0
            ORIGIN
                    1 acgt
            //
            """
    a = read(input_)
    assert str(a.seq) == "ACGT"

    input_ = """>hej
                acgt"""
    assert str(a.seq) == "ACGT"

    input_ = """>hej öööh!
                acgt"""
    assert str(a.seq) == "ACGT"

    input_ = """
                LOCUS       New_DNA                    4 bp ds-DNA     linear       30-MAR-2013
                DEFINITION  öööh!
                ACCESSION
                VERSION
                SOURCE      .
                  ORGANISM  .
                COMMENT
                COMMENT     ApEinfo:methylated:1
                FEATURES             Location/Qualifiers
                     misc_feature    2..3
                                     /label=öööh!
                                     /ApEinfo_fwdcolor=cyan
                                     /ApEinfo_revcolor=green
                                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                     width 5 offset 0
                ORIGIN
                        1 acgt
                //
            """
    a = read(input_)
    assert str(a.seq) == "ACGT"


def test_read_from_unicode():

    with open("pth1.txt", "r", encoding="utf-8") as f:
        text = f.read()
    assert isinstance(text, str)
    x, y = parse(text)
    assert x.format()[3268:3278] == "2micron 2µ"


def test_read_from_file():

    a = read("read1.gb")
    b = read("read2.gb")
    c = read("read3.fasta")
    d = read("read4.fasta")
    e = read(Path("read1.gb"))

    a.format("gb")
    b.format("gb")
    c.format("gb")
    d.format("gb")
    e.format("gb")

    assert (
        str(a.seq).lower()
        == str(b.seq).lower()
        == str(c.seq).lower()
        == str(d.seq).lower()
    )


def test_read_with_feature_spanning_ori():

    test = """
    LOCUS       New_DNA                   10 bp ds-DNA     circular     23-AUG-2018
    DEFINITION  .
    ACCESSION
    VERSION
    SOURCE      .
      ORGANISM  .
    COMMENT
    COMMENT     ApEinfo:methylated:1
    FEATURES             Location/Qualifiers
         misc_feature    join(9..10,1..2)
                         /locus_tag="myfeature"
                         /label="myfeature"
                         /ApEinfo_label="myfeature"
                         /ApEinfo_fwdcolor="cyan"
                         /ApEinfo_revcolor="green"
                         /ApEinfo_graphicformat="arrow_data {{0 1 2 0 0 -1} {} 0}
                         width 5 offset 0"
    ORIGIN
            1 accgggtttt
    //
    """
    a = read(test)
    assert str(a.seq).lower() == "accgggtttt"
    assert str(a.features[0].extract(a).seq) == "TTAC"
    assert a.features[0].location.strand == 1

    b = a.rc()

    assert str(b.seq).lower() == "aaaacccggt"
    assert str(b.features[0].extract(a).seq) == "GTAA"
