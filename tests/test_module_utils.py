#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_flatten():
    from Bio.Seq import Seq as bSeq
    from pydna.dseq import Dseq
    from Bio.SeqRecord import SeqRecord as bSeqRecord
    from pydna.seqrecord import SeqRecord
    from pydna.dseqrecord import Dseqrecord
    from pydna.utils import flatten

    testtuple = (
        [1, 2, 3],
        [bSeq("a"), bSeq("a"), bSeq("a")],
        [Dseq("a"), Dseq("a"), Dseq("a")],
        [bSeqRecord(bSeq("a")), bSeqRecord(bSeq("a")), bSeqRecord(bSeq("a"))],
        [SeqRecord(bSeq("a")), SeqRecord(bSeq("a")), SeqRecord(bSeq("a"))],
        [Dseqrecord(bSeq("a")), Dseqrecord(bSeq("a")), Dseqrecord(bSeq("a"))],
        [Dseqrecord("a"), Dseqrecord("a"), Dseqrecord("a")],
        ["a", "b", "c"],
    )

    for t in testtuple:
        assert flatten(t) == t
    testtuple2 = (([1, [2, 3]], [1, 2, 3]),)
    for argument, result in testtuple2:
        assert flatten(argument) == result


def test_eq():
    from pydna.dseqrecord import Dseqrecord
    from pydna.utils import eq
    from Bio.Seq import Seq
    from pydna.seq import Seq as Seq_
    from pydna.seqrecord import SeqRecord as SeqRecord_
    from Bio.SeqRecord import SeqRecord

    assert eq("AAA", "TTT", linear=True)
    assert eq("AAA", "TTT", linear=False)

    assert eq("aAA", "TtT", linear=True)
    assert eq("AAa", "TtT", linear=False)

    assert eq("ATA", "AAT", circular=True)
    assert not eq("ATA", "AAT", circular=False)
    assert eq("AAA", "AAA", linear=True)
    assert eq("AAA", "AAA", linear=False)

    assert eq("ATA", Seq("AAT"), circular=True)
    assert not eq("ATA", Seq("AAT"), circular=False)
    assert eq("AAA", Seq("AAA"), linear=True)
    assert eq("AAA", Seq("AAA"), linear=False)

    assert eq("ATA", SeqRecord_("AAT"), circular=True)
    assert not eq("ATA", SeqRecord_("AAT"), circular=False)
    assert eq("AAA", SeqRecord_("AAA"), linear=True)
    assert eq("AAA", SeqRecord_("AAA"), linear=False)

    assert eq(Seq_("ATA"), SeqRecord_(Seq_("AAT")), circular=True)
    assert not eq(Seq_("ATA"), SeqRecord_(Seq_("AAT")), circular=False)
    assert eq(Seq_("AAA"), SeqRecord_(Seq_("AAA")), linear=True)
    assert eq(Seq_("AAA"), SeqRecord_(Seq_("AAA")), linear=False)

    assert eq("ATA", Dseqrecord("AAT"), circular=True)
    assert not eq("ATA", Dseqrecord("AAT"), circular=False)
    assert eq("AAA", Dseqrecord("AAA"), linear=True)
    assert eq("AAA", Dseqrecord("AAA"), linear=False)

    assert eq(Seq("ATA"), SeqRecord(Seq("AAT")), circular=True)
    assert not eq(Seq("ATA"), SeqRecord(Seq("AAT")), circular=False)
    assert eq(Seq("AAA"), SeqRecord(Seq("AAA")), linear=True)
    assert eq(Seq("AAA"), SeqRecord(Seq("AAA")), linear=False)

    assert eq(Seq("ATA"), Dseqrecord("AAT"), circular=True)
    assert not eq(Seq("ATA"), Dseqrecord("AAT"), circular=False)
    assert eq(Seq("AAA"), Dseqrecord("AAA"), linear=True)
    assert eq(Seq("AAA"), Dseqrecord("AAA"), linear=False)

    assert eq(Dseqrecord("AAA", circular=False), Dseqrecord("AAA", circular=False))
    assert eq(Dseqrecord("AAA", circular=True), Dseqrecord("AAA", circular=True))
    assert not eq(Dseqrecord("ATA", circular=False), Dseqrecord("AAT", circular=False))
    assert eq(Dseqrecord("ATA", circular=True), Dseqrecord("AAT", circular=True))

    with pytest.raises(ValueError):
        eq(Dseqrecord("ATA", circular=True), Dseqrecord("ATA", circular=False))

    assert not eq(Dseqrecord("ATA", circular=True), Dseqrecord("ATAA", circular=True))

    assert eq(Dseqrecord("ATA"), Dseqrecord("ATA"), circular=True)
    assert not eq(Dseqrecord("ATA"), Dseqrecord("CCC"), circular=True)
    assert not eq(
        Dseqrecord("ATA"), Dseqrecord("ATA"), Dseqrecord("CCC"), circular=True
    )


# def test_shift_origin():
#    from pydna.readers import read
#    from pydna.dseqrecord import Dseqrecord
#
#
#    from pydna.utils import shift_origin, eq
#    from Bio.Seq import Seq
#    from Bio.SeqRecord import SeqRecord
#    pCAPs   = read("pCAPs.gb")
#    assert pCAPs.circular
#    pCAPs_b = shift_origin(pCAPs, 200)
#    assert len(pCAPs) == len(pCAPs_b)
#    assert pCAPs_b.circular
#    assert eq(pCAPs, pCAPs_b)
#    pCAPs_b_linear = pCAPs_b.tolinear()
#    assert eq(pCAPs, pCAPs_b_linear, circular=True)
#    pCAPs_c = pCAPs[200:]+pCAPs[:200]
#    assert eq(pCAPs, pCAPs_c, circular=True)
#    #with self.assertRaisesRegex(ValueError, "shift"):
#    #    pCAPs_b = shift_origin(pCAPs, 20000)


# def test_copy_features():
#    from pydna.readers import read
#    from pydna.utils import seguid, copy_features
#
#    a=read("pCAPs.gb")
#    b=read("pCAPs_fasta.txt")
#
#    for sh in [1,2,3,3127,3128,3129]:
#        newb = (b[sh:]+b[:sh]).looped()
#        copy_features(a, newb)
#        #print "a",[len(str(f.extract(a).seq.lower()) for f in a.features if len(f)>10]
#        #print "b",[len(str(f.extract(newb).seq).lower()) for f in newb.features]
#        x= sorted([str(f.extract(a).seq).lower() for f in a.features if len(f)>10],key=len)
#        y= sorted([str(f.extract(newb).seq).lower() for f in newb.features],key=len)
#        assert x==y
#
#    b=b.rc()
#
#    for sh in [1,2,3,3127,3128,3129]:
#        newb = b[sh:]+b[:sh]
#        copy_features(a, newb)
#
#
#        x = sorted([str(f.extract(a).seq).lower() for f in a.features if len(f)>10],key=len)
#        y = sorted([str(f.extract(newb).seq).lower() for f in newb.features],key=len)
#        assert x==y
#
#    seguid_bla = "riT98j2v4NxVS8sbw_Q8epCwQwo"
#    seguid_cre = "xLZ2xs2O8CUMmWh2OrhmNFp5ZLg"
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre, seguid_cre, seguid_bla, seguid_bla]
#
#    b=read("pCAPs_fasta.txt").looped()
#
#    b=b.synced("attaacgagtgccgtaaacgacgatggttttacc")
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,seguid_bla,seguid_bla]
#
#    b=read("pCAPs_fasta.txt").looped()
#    b=b.synced("ttaacgagtgccgtaaacgacgatggttttacc")
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,seguid_bla,seguid_bla]
#
#    b=read("pCAPs_fasta.txt").looped()
#    b=b.synced("taacgagtgccgtaaacgacgatggttttacc")
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_bla,seguid_bla]
#
#    b=read("pCAPs_fasta.txt").looped()
#    b=b.synced("gttaccaatgcttaatcagtgaggcacctatctcagc")
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,seguid_bla,seguid_bla]
#
#    b=read("pCAPs_fasta.txt").looped()
#    b=b.synced("ttaccaatgcttaatcagtgaggcacctatctcagc")
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,seguid_bla,seguid_bla]
#
#    b=read("pCAPs_fasta.txt").looped()
#    b=b.synced("taccaatgcttaatcagtgaggcacctatctcagc")
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,]


# def test_cseguid():
#     from pydna.utils import cseguid

#     x = "tcgcgcgtttcggtgatgacggtgAAAAcctctgacacatgcagctcccggattgtactgagagtgc"
#     assert cseguid(x) == cseguid(x.upper()) == cseguid(x.lower()) == "naaZmDzyMa58OsNXROe5SvjC7WU"


# def test_lseguid():
#     from pydna.utils import lseguid_blunt

#     x = "tcgcgcgtttcggtgatgacggtgAAAAcctctgacacatgcagctcccggattgtactgagagtgc"
#     assert lseguid_blunt(x) == lseguid_blunt(x.upper()) == lseguid_blunt(x.lower()) == "bHrqalTJ793oAigMQ5_qCttJRTk"


def test_rc():
    from pydna.utils import rc

    assert rc("gattc") == "gaatc"
    assert rc("GATTC") == "GAATC"


def test_complement():
    from pydna.utils import complement

    assert complement("gattc") == "ctaag"
    assert complement("GATTC") == "CTAAG"


def test_seq31():
    from pydna.utils import seq31

    assert (
        seq31("MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer")
        == "M  A  I  V  M  G  R  W  K  G  A  R  *"
    )


# def test_parse_text_table():
#     from pydna.utils import parse_text_table

#     table1 = "one     two   three\n" "  four   five    six\n" "seven   eight     nine\n"
#     (
#         formatted,
#         columnsplit,
#         rowsplit,
#         list_of_lists_rc,
#         list_of_lists_cr,
#     ) = parse_text_table(table1)

#     assert formatted == "one   two   three\nfour  five  six  \nseven eight nine "

#     cs = "one\n" "four\n" "seven\n" "|||\n" "two\n" "five\n" "eight\n" "|||\n" "three\n" "six\n" "nine"

#     assert cs == columnsplit

#     rs = "one\n" "two\n" "three\n" "---\n" "four\n" "five\n" "six\n" "---\n" "seven\n" "eight\n" "nine"

#     assert rs == rowsplit


# def test_join_list_to_table():
#     from pydna.utils import join_list_to_table

#     cs = "one\n" "four\n" "seven\n" "|||\n" "two\n" "five\n" "eight\n" "|||\n" "three\n" "six\n" "nine"

#     assert join_list_to_table(cs) == "one   two   three\nfour  five  six  \nseven eight nine "

#     cs = "one\n" "four\n" "seven\n" "|||\n" "two\n" "five\n" "\n" "|||\n" "three\n" "six\n" "nine"

#     answer = 'one   two  three\nfour  five six  \nseven "    nine '
#     assert join_list_to_table(cs) == answer

#     rs = "one\n" "two\n" "three\n" "---\n" "four\n" "five\n" "six\n" "---\n" "seven\n" "eight\n" "nine"

#     assert join_list_to_table(rs) == "one   two   three\nfour  five  six  \nseven eight nine "

#     assert join_list_to_table("somestring") is None

#     cs = "one\n" "four\n" " \n" "|||\n" "two\n" " \n" "eight"

#     assert join_list_to_table(cs) == 'one  two  \nfour "    \n"    eight'


# def test_expandtolist():
#     from pydna.utils import expandtolist

#     samples = "Sample [1..3] prepared according to [A..C]"

#     result = (
#         "Sample 1 prepared according to A\n" "Sample 2 prepared according to B\n" "Sample 3 prepared according to C\n"
#     )

#     assert result == expandtolist(samples)

#     samples = "Sample [1..4] prepared according to [A..B]"

#     result = (
#         "Sample 1 prepared according to A\n"
#         "Sample 2 prepared according to A\n"
#         "Sample 3 prepared according to B\n"
#         "Sample 4 prepared according to B\n"
#     )

#     assert result == expandtolist(samples)


def test_randomRNA():
    from pydna.utils import randomRNA

    rna = randomRNA(10)
    assert len(rna) == 10
    assert all([nt in "GAUC" for nt in rna])
    rna = randomRNA(10, 20)
    assert 10 <= len(rna) <= 20
    assert all([nt in "GAUC" for nt in rna])


def test_randomDNA():
    from pydna.utils import randomDNA

    dna = randomDNA(10)
    assert len(dna) == 10
    assert all([nt in "GATC" for nt in dna])
    dna = randomDNA(10, 20)
    assert 10 <= len(dna) <= 20
    assert all([nt in "GATC" for nt in dna])


def test_randomORF():
    from pydna.utils import randomORF

    orf = randomORF(3)
    assert orf.startswith("ATG")

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

    assert orf[3:6] in cdns

    assert orf.endswith(("TAA", "TAG", "TGA"))

    assert 9 <= len(randomORF(3, 6)) <= 18


def test_randomprot():
    from pydna.utils import randomprot

    protein = randomprot(10)
    assert len(protein) == 10
    assert all([aa in "ACDEFGHIKLMNPQRSTVWY" for aa in protein])
    protein = randomprot(10, 20)
    assert 10 <= len(protein) <= 20
    assert all([aa in "ACDEFGHIKLMNPQRSTVWY" for aa in protein])


def test_smallest_rotation():
    from pydna.utils import smallest_rotation as sr

    assert sr("tttaaa") == "aaattt"


def test_shift_location():
    from pydna.utils import shift_location
    from Bio.SeqFeature import SimpleLocation

    # TODO: more tests here

    # Shifting of locations should be reversible (https://github.com/pydna-group/pydna/issues/195)
    for strand in (1, -1, None):
        loc = SimpleLocation(0, 2, strand)
        assert shift_location(shift_location(loc, 1, 6), -1, 6) == loc

    # Shifting location on circular sequence
    for strand in (1, -1, None):
        loc = SimpleLocation(0, 4, strand)
        assert shift_location(loc, 1, 6) == SimpleLocation(1, 5, strand)
        if strand == -1:
            assert shift_location(loc, -1, 6) == SimpleLocation(
                0, 3, strand
            ) + SimpleLocation(5, 6, strand)
        else:
            assert shift_location(loc, -1, 6) == SimpleLocation(
                5, 6, strand
            ) + SimpleLocation(0, 3, strand)

    # Shifting ignoring the sequence length
    # See https://github.com/pydna-group/pydna/issues/281
    for strand in (1, -1, None):
        loc = SimpleLocation(4, 6, strand)
        assert shift_location(loc, 1000, None) == SimpleLocation(1004, 1006, strand)
        assert shift_location(loc, -4, None) == SimpleLocation(0, 2, strand)
        try:
            shift_location(loc, -1000, None)
            raise AssertionError("Shift below zero should raise ValueError")
        except ValueError:
            pass

        composed_loc = SimpleLocation(2, 4, strand) + SimpleLocation(5, 6, strand)
        assert shift_location(composed_loc, 1000, None) == SimpleLocation(
            1002, 1004, strand
        ) + SimpleLocation(1005, 1006, strand)
        assert shift_location(composed_loc, -2, None) == SimpleLocation(
            0, 2, strand
        ) + SimpleLocation(3, 4, strand)


def test_locations_overlap():
    from pydna.utils import locations_overlap, shift_location
    from Bio.SeqFeature import SimpleLocation

    # exact        =====          |
    # greater    =========        |
    # inner         ==            |
    # right         =====         |
    # left        =====           |
    # main         =====          |
    #         -----------------------------------------
    #         0123456789
    main_overlap = SimpleLocation(5, 10)
    inner_overlap = SimpleLocation(6, 8)
    right_overlap = SimpleLocation(6, 11)
    left_overlap = SimpleLocation(4, 9)
    exact_overlap = SimpleLocation(5, 10)
    greater_overlap = SimpleLocation(3, 12)
    no_overlap_left = SimpleLocation(0, 5)
    no_overlap_right = SimpleLocation(11, 15)

    overlapping_locations = [
        inner_overlap,
        right_overlap,
        left_overlap,
        exact_overlap,
        greater_overlap,
    ]
    non_overlapping_locations = [no_overlap_left, no_overlap_right]

    for shift in range(20):
        main_shifted = shift_location(main_overlap, shift, 20)
        for loc in overlapping_locations:
            loc_shifted = shift_location(loc, shift, 20)
            assert locations_overlap(main_shifted, loc_shifted, 20)
        for loc in non_overlapping_locations:
            loc_shifted = shift_location(loc, shift, 20)
            assert not locations_overlap(main_shifted, loc_shifted, 20)


def test_create_location():
    from pydna.utils import create_location

    # Basic location creation with different strands
    assert str(create_location(0, 5, 10, -1)) == "[0:5](-)"
    assert str(create_location(0, 5, 10, +1)) == "[0:5](+)"
    assert str(create_location(0, 5, 10)) == "[0:5]"

    # Circular location (end < start)
    assert str(create_location(8, 2, 10)) == "join{[8:10], [0:2]}"
    assert str(create_location(8, 2, 10, -1)) == "join{[0:2](-), [8:10](-)}"

    # Negative positions (should be normalized)
    assert str(create_location(-2, 2, 10)) == "join{[8:10], [0:2]}"

    # Special case: 0 is the same as len(seq)
    assert str(create_location(5, 0, 10)) == "[5:10]"

    # Special case: start and end are the same, spans entire sequence
    assert str(create_location(5, 5, 10)) == "join{[5:10], [0:5]}"


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
