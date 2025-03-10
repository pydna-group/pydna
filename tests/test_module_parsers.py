#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test parse
"""

import pytest


def test_extract_from_text():

    text = """\
            >a
            aaaa
            LOCUS
            //
            >b
            bbbbbb
            ID
            //
            """
    from pydna.parsers import extract_from_text

    seqs, gaps = extract_from_text(text)
    assert seqs == (">a\naaaa\n", "LOCUS\n//", ">b\nbbbbbb\n", "ID\n//")
    assert [g.strip() for g in gaps] == ["", "", "", "", ""]
    text = """\
    comment 0
    LOCUS a
    //
    comment 1
    LOCUS b
    //
    comment 2
    >c
    ccccc

    comment 3
    >ddd
    dddddd
    ID e
    //
    comment 4
    """
    seqs, gaps = extract_from_text(text)
    assert seqs == ("LOCUS a\n//", "LOCUS b\n//", ">c\nccccc", ">ddd\ndddddd\n", "ID e\n//")
    assert tuple(g.strip() for g in gaps) == ("comment 0", "comment 1", "comment 2", "comment 3", "", "comment 4")

    from pydna.parsers import embl_gb_fasta

    text = """\
            LOCUS       New_linear_DNA             2 bp    DNA     linear       29-MAR-2024
            DEFINITION  .
            ACCESSION
            VERSION
            SOURCE      .
              ORGANISM  .
            ORIGIN
                    1 aa
            //
            LOCUS       New_circular_DNA           2 bp    DNA     circular     29-MAR-2024
            DEFINITION  .
            ACCESSION
            VERSION
            SOURCE      .
              ORGANISM  .
            ORIGIN
                    1 aa
            //
            """

    lin, crc = embl_gb_fasta(text)

    assert lin.annotations.get("topology") == "linear"

    assert crc.annotations.get("topology") == "circular"

    text = """\
            >a
            aaa
            >c
            ccc
            >g
            ggg
            >t
            ttt
            """

    a, c, t, g = embl_gb_fasta(text)

    assert [x.annotations.get("topology") for x in (a, c, g, t)] == ["linear", "linear", "linear", "linear"]

    text = """\
            >a circular
            aaa
            >c circular
            ccc
            >g circular
            ggg
            >t circular
            ttt
            """

    a, c, t, g = embl_gb_fasta(text)

    assert [x.annotations.get("topology") for x in (a, c, g, t)] == ["circular", "circular", "circular", "circular"]


def test_parse1():
    from pydna.parsers import parse
    from pydna.readers import read

    """ test parsing fasta sequences from a text"""

    text = """
            points....: 1

            The sequence seq below represents a double stranded linear DNA molecule.

            >seq
            CTCCCCTATCACCAGGGTACCGATAGCCACGAATCT

            Give the sequence(s) of the fragment(s) formed after digesting seq
            with the restriction enzyme Acc65I in the order that they appear in seq.

            Use FASTA format and give the Watson strand(s) in 5'-3' direction below.
            Give the sequences the names frag1,frag2,... etc.
            >frag1
            CTCCCCTATCACCAGG

            >frag2
            GTACCGATAGCCACGAATCT

            *********** Question 4 ***********

            QuestionID:
            """
    result = parse(text)

    correct = [
        "CTCCCCTATCACCAGGGTACCGATAGCCACGAATCT",
        "CTCCCCTATCACCAGG",
        "GTACCGATAGCCACGAATCT",
    ]

    assert [str(s.seq) for s in result] == correct
    assert [s.circular for s in result] == [False, False, False]

    input = """
            LOCUS       ScCYC1                   330 bp    DNA              UNK 01-JAN-1980
            DEFINITION  ScCYC1
            ACCESSION   ScCYC1
            VERSION     ScCYC1
            KEYWORDS    .
            SOURCE      .
              ORGANISM  .
                        .
            FEATURES             Location/Qualifiers
            ORIGIN
                    1 ATGACTGAAT TCAAGGCCGG TTCTGCTAAG AAAGGTGCTA CACTTTTCAA GACTAGATGT
                   61 CTACAATGCC ACACCGTGGA AAAGGGTGGC CCACATAAGG TTGGTCCAAA CTTGCATGGT
                  121 ATCTTTGGCA GACACTCTGG TCAAGCTGAA GGGTATTCGT ACACAGATGC CAATATCAAG
                  181 AAAAACGTGT TGTGGGACGA AAATAACATG TCAGAGTACT TGACTAACCC AAAGAAATAT
                  241 ATTCCTGGTA CCAAGATGGC CTTTGGTGGG TTGAAGAAGG AAAAAGACAG AAACGACTTA
                  301 ATTACCTACT TGAAAAAAGC CTGTGAGTAA
            //
            """
    result = parse(input).pop()

    assert str(result.seq) == str(read(input).seq)
    correct = """ATGACTGAATTCAAGGCCGGTTCTGCTAAGAAAGGTGCTACACTTTTCAAGACTAGATGTCTACAATGCCACACCGTGGAAAAGGGTGGCCCACATAAGGTTGGTCCAAACTTGCATGGTATCTTTGGCAGACACTCTGGTCAAGCTGAAGGGTATTCGTACACAGATGCCAATATCAAGAAAAACGTGTTGTGGGACGAAAATAACATGTCAGAGTACTTGACTAACCCAAAGAAATATATTCCTGGTACCAAGATGGCCTTTGGTGGGTTGAAGAAGGAAAAAGACAGAAACGACTTAATTACCTACTTGAAAAAAGCCTGTGAGTAA"""

    assert str(result.seq) == correct
    assert result.circular == False

    seqs = parse("RefDataBjorn.fas")
    assert len(seqs) == 771
    assert list(set([len(a) for a in seqs])) == [901]

    pAG25 = read("pAG25.gb")
    assert pAG25.circular == True

    pCAPs = read("pCAPs.gb")
    assert pCAPs.circular == True

    pUC19 = read("pUC19.gb")
    assert pUC19.circular == True

    input = """
    ID   example    standard; DNA; UNC; 3 BP.
    SQ   Sequence 3 BP;
         aaa                                                                       3
    //
    """
    result = parse(input).pop()
    assert str(result.seq) == "AAA"

    input = """
    ID   name?      standard; circular DNA; UNK; 100 BP.
    XX
    DT   25-DEC-2017
    XX
    DE   description?.
    XX
    AC   id?;
    XX
    SV   id?
    XX
    KW   .
    XX
    OS   .
    OC   .
    OC   .
    XX
    FH   Key             Location/Qualifiers
    SQ   Sequence 100 BP;
         aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa        60
         aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa                             100
    //
    """
    result = parse(input).pop()
    assert str(result.seq) == "A" * 100


def test_parse2():
    from pydna.parsers import parse
    from pydna.readers import read

    seqs = parse("RefDataBjorn.fas")

    assert len(seqs) == 771
    assert list(set([len(a) for a in seqs])) == [901]

    for i, s in enumerate(seqs):
        a = s.description
        b = a.split()
        c = "|".join([b[0], b[1], b[3]])
        s.id = b[2].replace(" ", "_") + "_" + str(i)
        s.description = ""
        if b[3] == "Zenion hololepis":
            s.id = b[3].replace(" ", "_") + "_" + str(i)


def test_parse_primers():
    from pydna.parsers import parse_primers

    data = str(">1\n" "aaaa\n" ">2\n" "cccc\n")
    parse_primers(data)

    f0, r0 = parse_primers(
        """
                             >ForwardPrimer
                             gctactacacacgtactgactg

                             >ReversePrimer
                             tgtggttactgactctatcttg"""
    )
    assert str(f0.seq) == "gctactacacacgtactgactg"
    assert str(r0.seq) == "tgtggttactgactctatcttg"


def test_parse_error():
    from pydna.parsers import parse

    data = """
LOCUS
DATA_IS_NOT_A_SEQUENCE
//"""
    assert parse(data) == []


def test_parse_list():
    from pydna.parsers import parse_primers

    data = str(">1\n" "aaaa\n" ">2\n" "cccc\n")

    assert [str(x.seq) for x in parse_primers([data, data])] == ["aaaa", "cccc", "aaaa", "cccc"]


def test_misc_parse():

    from pydna.parsers import parse
    from Bio.SeqIO import read as BPread
    from Bio.SeqIO import parse as BPparse

    q = BPread("read1.gb", "gb")
    w = BPread("read2.gb", "gb")
    e = BPread("read3.fasta", "fasta")
    r = BPread("read4.fasta", "fasta")

    with open("pth1.txt", "r", encoding="utf-8") as f:
        a, b = BPparse(f, "gb")

    assert "|" + a.features[13].qualifiers["label"][0] + "|" == "|2micron 2µ|"
    assert "|" + a.format("gb")[3314:3324] + "|" == '|olor="gree|'

    assert a.features[13].qualifiers["label"][0] == "2micron 2µ"
    assert a.format("gb")[3268:3278] == "2micron 2µ"

    x, y = parse("pth1.txt")
    x.format("gb")
    y.format("gb")
    assert x.format()[3268:3278] == "2micron 2µ"
    assert x.features[13].qualifiers["label"][0] == "2micron 2µ"

    assert "".join(a.format("gb").splitlines()[1:]) == "".join(x.format("gb").splitlines()[1:])
    assert "".join(b.format("gb").strip().splitlines()[4:]) == "".join(y.format("gb").splitlines()[4:])


def test_dna2949():
    from pydna.parsers import parse

    with open("dna2943.gb") as f:
        f.read()
    seqlist = parse("dna2943.gb", ds=True)
    assert len(seqlist) == 1
    assert seqlist[0].seguid() == "ldseguid=ScLoSddUf2c0GIAGpvIi33nLvFY"


def proteins():
    from pydna.parsers import embl_gb_fasta

    proteins = """\
    >pdb|3VQM|V Chain V, C-terminal peptide from Small heat shock protein StHsp14.0
    VIKIE

    LOCUS       3VQM_W                     5 aa            linear   SYN 08-NOV-2023
    DEFINITION  Chain W, C-terminal peptide from Small heat shock protein
                StHsp14.0.
    ACCESSION   3VQM_W
    VERSION     3VQM_W
    DBSOURCE    pdb: molecule 3VQM, chain W, release Nov 8, 2023;
                deposition: Mar 26, 2012;
                class: CHAPERONE;
                source: Mmdb_id: 100300, Pdb_id 1: 3VQM;
                Exp. method: X-ray Diffraction.
    KEYWORDS    .
    SOURCE      synthetic construct
      ORGANISM  synthetic construct
                other sequences; artificial sequences.
    REFERENCE   1  (residues 1 to 5)
      AUTHORS   Hanazono,Y., Takeda,K., Yohda,M. and Miki,K.
      TITLE     Structural studies on the oligomeric transition of a small heat
                shock protein, StHsp14.0
      JOURNAL   J Mol Biol 422 (1), 100-108 (2012)
       PUBMED   22613762
    REFERENCE   2  (residues 1 to 5)
      AUTHORS   Hanazono,Y., Takeda,K. and Miki,K.
      TITLE     Direct Submission
      JOURNAL   Submitted (26-MAR-2012)
    COMMENT     Small heat shock protein hsp14.0 of C-terminal deletion variant
                with C-terminal peptide.
    FEATURES             Location/Qualifiers
         source          1..5
                         /organism="synthetic construct"
                         /db_xref="taxon:32630"
    ORIGIN
            1 vikie
    //
    """

    fa, gb = embl_gb_fasta(proteins)

    assert fa.annotations["molecule_type"] == "protein"
    assert gb.annotations["molecule_type"] == "protein"

    assert fa.annotations["topology"] == "linear"
    assert gb.annotations["topology"] == "linear"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
