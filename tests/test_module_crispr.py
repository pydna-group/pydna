#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from Bio.Restriction import BamHI
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from Bio.Restriction import SapI
from pydna.parsers import parse_primers
from pydna.crispr import cas9, protospacer
from pydna.assembly import Assembly
from textwrap import dedent
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import FormattedSeq


def test_crispr():

    a = Dseq.from_representation(
        """\
    GTTACTTTACCCGACGT
    CAATGAAATGGGCTGCA
    """
    )

    b = Dseq.from_representation(
        """\
    CCCaGG
    GGGtCC
    """
    )

    sgr_text = "GTTACTTTACCCGACGTCCCgttttagagctagaaatagcaagttaaaataagg"
    target_string = "GTTACTTTACCCGACGTCCCaGG"

    for _sg, _tgt in [
        (sgr_text, target_string),
        (sgr_text.upper(), target_string.lower()),
        (sgr_text.lower(), target_string.upper()),
    ]:
        containing_sgRNA = Dseqrecord(sgr_text)
        target = Dseqrecord(target_string)

        assert [
            f.seq
            for f in target.cut([cas9(ps) for ps in protospacer(containing_sgRNA)])
        ] == [a, b]
        assert [
            f.seq
            for f in target.cut([cas9(ps) for ps in protospacer(containing_sgRNA.rc())])
        ] == [a, b]
        assert [
            f.seq
            for f in target.rc().cut([cas9(ps) for ps in protospacer(containing_sgRNA)])
        ] == [b.rc(), a.rc()]
        assert [
            f.seq
            for f in target.rc().cut(
                [cas9(ps) for ps in protospacer(containing_sgRNA.rc())]
            )
        ] == [
            b.rc(),
            a.rc(),
        ]

    c9 = cas9("GTTACTTTACCCGACGTCCC")

    assert (
        str(c9)
        == ">cas9 protospacer scaffold\nGTTACTTTACCCGACGTCCC GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG"
    )

    assert target.cut(c9) == target.cut(cas9("GTTACTTTACCCGACGTCCC".lower()))

    with pytest.raises(TypeError):
        BamHI.search("GGATCC")

    assert BamHI.search(Seq("GGATCC")) == [2]
    assert BamHI.search(MutableSeq("GGATCC")) == [2]
    assert BamHI.search(Dseq("GGATCC")) == [2]
    with pytest.raises(TypeError):
        assert BamHI.search(SeqRecord("GGATCC")) == [2]
    with pytest.raises(TypeError):
        assert BamHI.search(FormattedSeq("GGATCC")) == [2]
    with pytest.raises(TypeError):
        assert BamHI.search(Dseqrecord("GGATCC")) == [2]
    assert BamHI.search(Seq("GATCCG"), linear=False) == [1]
    assert BamHI.search(MutableSeq("GATCCG"), linear=False) == [1]
    assert BamHI.search(Dseq("GATCCG"), linear=False) == [1]
    assert BamHI.search(Dseq("GATCCG", circular=True), linear=True) == []

    with pytest.raises(TypeError):
        c9.search(target_string)
    assert c9.search(Seq(target_string)) == [18]
    assert c9.search(MutableSeq(target_string)) == [18]
    assert c9.search(Dseq(target_string)) == [18]
    with pytest.raises(TypeError):
        assert c9.search(SeqRecord(target_string)) == [18]
    with pytest.raises(TypeError):
        assert c9.search(FormattedSeq(target_string)) == [18]
    with pytest.raises(TypeError):
        assert c9.search(Dseqrecord(target_string)) == [18]
    rotated_target_string = target_string[1:] + target_string[0]
    assert c9.search(Seq(rotated_target_string), linear=False) == [17]
    assert c9.search(MutableSeq(rotated_target_string), linear=False) == [17]
    assert c9.search(Dseq(rotated_target_string), linear=False) == [17]
    assert c9.search(Dseq("ggat", circular=True), linear=True) == []


def test_crispr_torulaspora():
    """
    1. Two single stranded oligos are annealed in-vitro
    2. They form a 26 bp dsDNA fragment with 5' overhangs
    3. This fragment is cloned in plasmid pTgCRISPRplanB cut with SapI (2 cuts)
       The resulting pTgCRISPR_guide has a complete guide sequrnce with scaffold
    4. The pTgCRISPR_guide is used to instantiate a Cas9 object (gr)
    5. torulaspora_locus is cut into two fragments (before & after)
    6. before, repair_template, after are assembled into a new deletion locus.

    Lombardi, L., Oliveira-Pacheco, J., & Butler, G. (2019).
    Plasmid-Based CRISPR-Cas9 Gene Editing in Multiple Candida
    Species. mSphere, 4(2), e00125-19.
    https://doi.org/10.1128/mSphere.00125-19

    """

    guides = """
    >NewGuide_Fw
    CCAAAGGCTGCTGCTTATGCCAG

    >NewGuide_Rv
    AACCTGGCATAAGCAGCAGCCTT
    """
    watson, crick = parse_primers(guides)

    guide_insert = Dseqrecord(Dseq(str(watson.seq), str(crick.seq)))

    assert repr(guide_insert.seq) == dedent(
        """\
                                            Dseq(-26)
                                            CCAAAGGCTGCTGCTTATGCCAG
                                               TTCCGACGACGAATACGGTCCAA"""
    )

    # This is only a partial sequence of the plasmid
    pTgCRISPRplanB = Dseqrecord(
        "gggcgtgtggcgtagttggtagcgcgttcccttagcatgg"
        "gaaaggtcatgagttcgattcttatctcgtccaggaagag"
        "ctcgcgagctcttccgttttagagctagaaatagcaagtt"
        "aaaataaggctagtccgttatcaacttgaaaaagtggcac"
        "cgagtcggtgc",
        name="plasmid",
        circular=True,
    )

    stuffer, linear_plasmid = pTgCRISPRplanB.cut(SapI)

    assert repr(linear_plasmid.seq) == dedent(
        """\
                                              Dseq(-149)
                                              gttttag..tcgt
                                                 aatc..agcaggt"""
    )

    pTgCRISPR_guide = (linear_plasmid + guide_insert).looped()

    gr = cas9(protospacer(pTgCRISPR_guide).pop())

    assert protospacer(pTgCRISPR_guide) == protospacer(pTgCRISPR_guide.rc())

    assert repr(gr) == "cas9(AAG..CAG)"

    repair_template = Dseqrecord(
        "aggaggcttacaacgaatcataatcccacgactactc"
        "tactttttacatagatgtaatatctccttatcgcact"
        "gtcgtcctcttgctgtttccgtaagcgatgagcctct"
        "cgcccaccttgtcgtcctgttggccggacccttgctt"
        "gtatttcaacagcaaactgcgccggactcttcctgat"
        "atgtacgtcgccgctctgttgtgacgtctaatagtct"
        "gggggtgtgattggttgatgactagacctccggaata"
        "acttcgaagcatttcttttaagagatggtctgctcat"
        "tcgcagaaggagtgtcattgagacgtttagctggtcg"
        "acgtcgagctatataaacgtgcgggtggtacaagatg"
        "aacattttctagattcagttcgagagcagcaacattc"
        "ttttaagctcacaaatgcatgctcacactcggtcata"
        "atgtaatgaaccaacctacctatttttttatcgttta"
        "tttaccaattacgttcggaaaaacttgttgttccgag"
        "ggacagtttgctgccagcaaccggaaaaactttgacc"
        "cgagggacaaggctgaagctcgaaatgccaaggctgt"
        "cgtcgagctgggcaaagaaaggtgcgatatttcaaat"
        "cacttattctacttttcttgttgttttagcacgttgg"
        "agcctaaatcttggaagcttgtaaactcaaaatcgca"
        "agatctagtaggaatattcggaaatttcgagttggag"
        "gttgtttctccccgctgaaattcgtagtctcttgtta"
        "ttcaggtctctcgtatgggttcgctgctggtatggtc"
        "ccgatcccaccttcttagtagactttggctggtagtt"
        "ctcgggaaaacttggttgcttaggaaaagaaagtttt"
        "tccgggccctcacaaggccctgcagaa",
        name="repair_template",
    )

    torulaspora_locus = Dseqrecord(
        "AGTCCGATGTCATTACTGAGGAGGCTTACAACGAAT"
        "CATAATCCCACGACTACTCTACTTTTTACATAGATG"
        "TAATATCTCCTTATCGCACTGTCGTCCTCTTGCTGT"
        "TTCCGTAAGCGATGAGCCTCTCGCCCACCTTGTCGT"
        "CCTGTTGGCCGGACCCTTGCTTGTATTTCAACAGCA"
        "AACTGCGCCGGACTCTTCCTGATATGTACGTCGCCG"
        "CTCTGTTGTGACGTCTAATAGTCTGGGGGTGTGATT"
        "GGTTGATGACTAGACCTCCGGAATAACTTCGAAGCA"
        "TTTCTTTTAAGAGATGGTCTGCTCATTCGCAGAAGG"
        "AGTGTCATTGAGACGTTTAGCTGGTCGACGTCGAGC"
        "TATATAAACGTGCGGGTGGTACAAGATGAACATTTT"
        "CTAGATTCAGTTCGAGAGCAGCAACATTCTTTTAAG"
        "CATGTCGGAAGGTCAAAATATAGCTGTTGCGACAAG"
        "CACtaattcaaattcaaattctaaCTCGACGTCCAT"
        "CGAGAAGACAGATAAGGCTGCTGCTTATGCCAGTGG"
        "TGGGACAGCTGCTGACGGGGCTGATTTGGAGCAGTT"
        "TAGAAGTCATGATGAAGGTGAGAAGGATGAGATGAC"
        "TCAGCTGGTCCGCGATGCGGACCAGCAGGCTGCTAC"
        "CATGAAGAAGTCTGACCTGATGTTTGTTTCGTTATG"
        "TTGTGTTATGGTTGCATTTGGGGGGTTTGTATTTGG"
        "TTGGGATACTGGTACTATATCTGGGTTTGTTAACCA"
        "GACCGACTGGATTAGACGGTTTGGGTCTTATAACCA"
        "AAATACTGGTCAGTACTACTTGTCTGATGTGCGTAC"
        "AGGTTTGCTGGTTTCGATCTTCAATATTGGGTGTGC"
        "CATTGGTGGTATTGTGCTCAGTAAAGCCGGTGACTT"
        "GTATGGTAGACGTTTGGGATTGGTTTTCGTTGTGGT"
        "CATTTATACCGTTGGTATATTAATTCAAATCTGTTC"
        "TTTCGATAAGTGGTACCAGTACTTTATTGGAAGGAT"
        "TATATCTGGTTTAGGTGTAGGTGGTATTACCGTTAT"
        "GTCGCCCATGTTGATATCTGAGGTGGCCCCCAAGCA"
        "AACTAGAGGTACCCTTATTTCGTGCTACCAGTTGAT"
        "GATCACTGGTGGTATCTTTTTAGGTTACTGTACCAA"
        "TTTTGGTACCAAGCAGTATGATGATTCTACGCAATG"
        "GAGAGTACCACTGGGTCTGTGCTTTGCCTGGGCTAT"
        "GTTTATGATCGGTGGTATGTTGTTCGTGCCAGAATC"
        "CCCACGTTACTTAGTAGAAGTTGGTAAGATGGAAGA"
        "GGCAAGACGTTCGTTGGCTAGAGCCAATAGGTGTCC"
        "CACAGAATCTCCATTGGTAACTTTGGAATTAGAAAA"
        "TTTCCAGGCTTCCGTAGAGGCTGAAAAGGTCGCGGG"
        "TAGCGCCAGTTGGTCCGAATTGATCACTGGTAAACC"
        "TGCTATCTTCAAGCGTACCCTGATGGGTGTTATGAT"
        "CCAGTCCTTGCAGCAATTGACTGGTGATAACTACTT"
        "CTTTTACTACGGTACTACGATTTTCAGAGCAGTCGG"
        "GTTGGAGGATTCGTTCGAGACTTCGATCGTTCTCGG"
        "TGTAGTCAACTTTGCCTCTACATTCCTTGCATTGTA"
        "CACTGTGGATCACTTTGGTCGTCGTAAATGTCTACT"
        "ATGGGGTTGTATTGGTATGGTCTGTTGTTATGTCGT"
        "TTATGCATCTGTCGGTGTCACTAGATTATGGCCAGA"
        "TGGACAAGATCAACCATCTTCCAAGGGTGCTGGTAA"
        "CTGTATGATTGTTTTCGCatgtttcttcatcttctg"
        "ttTCGCTACCACTTGGGCTCCAATTGCATACGTTAT"
        "CATCTCTGAATCTTTCCCATTGAGAATCAAGTCCAA"
        "GGCTATGTCCATTGCAAGTGCTTCTAACTGGCTATG"
        "GGGTTTCTTGATTGGTTTCTTCACTCCTTTCATCAC"
        "ATCCGCGATCAATTTCTATTACGGGTACGTGTTTAT"
        "GGGCTGTATGGTGTTTGCATTTTTCtatgtcttctt"
        "cttcgtcccAGAGACCAAAGGtttgactttggaaga"
        "agtgaaCATCATGTATGAAGAACGTGTCTTGCCTTG"
        "GAAGTCAGCATCTTGGTTACCTCCACACAGAAGGGG"
        "AGCCGACTTTGACGCAGATGCCTACACTAGGGACGA"
        "CCTCCCCTTCTACAAAAGAATGCTGCGTAGAAATTA"
        "ATCACAAATGCATGCTCACACTCGGTCATAATGTAA"
        "TGAACCAACCTACCTATTTTTTTATCGTTTATTTAC"
        "CAATTACGTTCGGAAAAACTTGTTGTTCCGAGGGAC"
        "AGTTTGCTGCCAGCAACCGGAAAAACTTTGACCCGA"
        "GGGACAAGGCTGAAGCTCGAAATGCCAAGGCTGTCG"
        "TCGAGCTGGGCAAAGAAAGGTGCGatatttcaaatc"
        "acttaTTCTacttttcttgttgtttTAGCACGTTGG"
        "AGCCTAAATCTTGGAAGCTTGTAAACTCAAAATCGC"
        "AAGATCTAGTAGGAATATTCGGAAATTTCGAGTTGG"
        "AGGTTGTTTCTCCCCGCTGAAATTCGTAGTCTCTTG"
        "TTATTCAGGTCTCTCGTATGGGTTCGCTGCTGGTAT"
        "GGTCCCGATCCCACCTTCTTAGTAGACTTTGGCTGG"
        "TAGTTCTCGGGAAAACTTGGTTGCTTAGGAAAAGAA"
        "AGTTTTTCCGGGCCCTCACAAGGCCCTGCAGAAGAG"
        "TTTT"
    )

    before, after = torulaspora_locus.cut(gr)

    before.name = "before"
    after.name = "after"

    assert len(before) == 534
    assert len(after) == 2170

    asm = Assembly((before, repair_template, after))

    deleted = asm.assemble_linear().pop()

    expected = (
        "before|415\n"
        "       \\/\n"
        "       /\\\n"
        "       415|repair_template|500\n"
        "                           \\/\n"
        "                           /\\\n"
        "                           500|after"
    )

    assert deleted.figure() == expected

    assert repr(deleted) == "Contig(-940)"

    assert deleted.seguid() == "ldseguid=HxGfGzivCuau8hBCOTAlzJq_FRo"

    final_deletion = Dseqrecord(
        "AGTCCGATGTCATTACTGaggaggcttacaacgaat"
        "cataatcccacgactactctactttttacatagatg"
        "taatatctccttatcgcactgtcgtcctcttgctgt"
        "ttccgtaagcgatgagcctctcgcccaccttgtcgt"
        "cctgttggccggacccttgcttgtatttcaacagca"
        "aactgcgccggactcttcctgatatgtacgtcgccg"
        "ctctgttgtgacgtctaatagtctgggggtgtgatt"
        "ggttgatgactagacctccggaataacttcgaagca"
        "tttcttttaagagatggtctgctcattcgcagaagg"
        "agtgtcattgagacgtttagctggtcgacgtcgagc"
        "tatataaacgtgcgggtggtacaagatgaacatttt"
        "ctagattcagttcgagagcagcaacattcttttaag"
        "cTCACAAATGCATGCTCACACTCGGTCATAATGTAA"
        "TGAACCAACCTACCTATTTTTTTATCGTTTATTTAC"
        "CAATTACGTTCGGAAAAACTTGTTGTTCCGAGGGAC"
        "AGTTTGCTGCCAGCAACCGGAAAAACTTTGACCCGA"
        "GGGACAAGGCTGAAGCTCGAAATGCCAAGGCTGTCG"
        "TCGAGCTGGGCAAAGAAAGGTGCGATATTTCAAATC"
        "ACTTATTCTACTTTTCTTGTTGTTTTAGCACGTTGG"
        "AGCCTAAATCTTGGAAGCTTGTAAACTCAAAATCGC"
        "AAGATCTAGTAGGAATATTCGGAAATTTCGAGTTGG"
        "AGGTTGTTTCTCCCCGCTGAAATTCGTAGTCTCTTG"
        "TTATTCAGGTCTCTCGTATGGGTTCGCTGCTGGTAT"
        "GGTCCCGATCCCACCTTCTTAGTAGACTTTGGCTGG"
        "TAGTTCTCGGGAAAACTTGGTTGCTTAGGAAAAGAA"
        "AGTTTTTCCGGGCCCTCACAAGGCCCTGCAGAAGAG"
        "TTTT"
    )

    assert deleted.seq == final_deletion.seq


def test_cut():
    seq = Dseqrecord("ttCATTACTGTTTGCATTGAACaGGCCC", circular=True)
    enz = cas9("CATTACTGTTTGCATTGAAC")
    expected = set([Dseq("AACAGGCCCTTCATTACTGTTTGCATTG").seguid()])

    for shift in range(len(seq)):
        seq_shifted = seq.shifted(shift)
        for rc in [False, True]:
            seq_shifted_rc = seq_shifted.reverse_complement() if rc else seq_shifted
            products = seq_shifted_rc.cut(enz)
            seguids = set(f.seq.seguid() for f in products)
            assert seguids == expected
