# -*- coding: utf-8 -*-
import pytest
from Bio.SeqFeature import SimpleLocation
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydna.assembly2 import (
    Assembly,
    recombinase_excision,
    recombinase_integration,
    recombinase_excision_or_inversion,
    recombinase_assembly,
)
from pydna.recombinase import (
    _recombinase_homology_offset_and_length,
    Recombinase,
    RecombinaseCollection,
)
from Bio.Seq import reverse_complement

# ---------------------------------------------------------------------------
# _recombinase_homology_offset_and_length
# ---------------------------------------------------------------------------


def test_recombinase_homology_offset_and_length_basic():
    site1 = "ATGCCCTAAaaTT"  # lowercase 'aa' at positions 9-10
    site2 = "AAaaTTTTTTTCCCT"  # lowercase 'aa' at positions 2-3

    off1, len1 = _recombinase_homology_offset_and_length(site1)
    off2, len2 = _recombinase_homology_offset_and_length(site2)

    assert (off1, len1) == (9, 2)
    assert (off2, len2) == (2, 2)


def test_recombinase_homology_offset_and_length_wrong_pattern_raises():
    for pattern in ["ATGCCCTAAATT", "ATGaa", "aaATG", "ATGaaCCaa", "ATGaaCCaaTT"]:
        with pytest.raises(ValueError):
            _recombinase_homology_offset_and_length(pattern)


def test_recombinase_homology_offset_and_length_invalid_chars():
    with pytest.raises(ValueError):
        _recombinase_homology_offset_and_length("ATGaa123")
    with pytest.raises(ValueError):
        _recombinase_homology_offset_and_length("ATGaa-CC")


# ---------------------------------------------------------------------------
# Recombinase.overlap – basic (linear, exact)
# ---------------------------------------------------------------------------


def test_recombinase_init_errors():
    # Different cores
    with pytest.raises(ValueError) as e:
        Recombinase("ATGacCTAAATT", "AAaaTTTTTTTCCCT")
    assert "Recombinase recognition sites do not have matching" in str(e.value)


def test_recombinase_overlap_single_match():
    site1 = "ATGCCCTAAaaTT"
    site2 = "AAaaTTTTTTTCCCT"

    seqA = Dseqrecord("aaaATGCCCTAAaaTTtt")
    seqB = Dseqrecord("tataAAaaTTTTTTTCCCTaaa")

    rec = Recombinase(site1, site2)
    matches = rec.overlap(seqA, seqB, limit=None)
    assert matches == [(12, 6, 2)]


def test_recombinase_overlap_no_site_returns_empty():
    site1 = "ATGCCCTAAaaTT"
    site2 = "AAaaTTTTTTTCCCT"

    seqA = Dseqrecord("aaaATGCCCTAAaaTTtt")
    seqB = Dseqrecord("tataAAccTTTTTTTCCCTaaa")

    rec = Recombinase(site1, site2)
    matches = rec.overlap(seqA, seqB, limit=None)
    assert matches == []


def test_recombinase_algorithm_with_assembly():
    site1 = "ATGCCCTAAaaTT"
    site2 = "AAaaTTTTTTTCCCT"

    seqA = Dseqrecord("aaaATGCCCTAAaaTTtt")
    seqB = Dseqrecord("tataAAaaTTTTTTTCCCTaaa")

    rec = Recombinase(site1, site2)
    asm = Assembly([seqA, seqB], algorithm=rec.overlap, use_fragment_order=False)
    products = asm.assemble_linear()

    assert len(products) == 2
    assert str(products[0].seq) == "aaaATGCCCTAAaaTTTTTTTCCCTaaa"
    assert str(products[1].seq) == "tataAAaaTTtt"


# ---------------------------------------------------------------------------
# Circular sequences
# ---------------------------------------------------------------------------


def test_circular_site_spanning_origin():
    """Sites that span the origin of a circular sequence should be found."""
    site1 = "ATGCCCTAAaaTT"
    site2 = "AAaaTTTTTTTCCCT"

    circ_seq = Dseqrecord("aaTTgggggATGCCCTAA", circular=True)
    linear_seq = Dseqrecord("tataAAaaTTTTTTTCCCTaaa")

    rec = Recombinase(site1, site2)
    for shift in range(len(circ_seq)):
        circ_seq_shifted = circ_seq.shifted(shift)
        matches = rec.overlap(circ_seq_shifted, linear_seq, limit=None)
        assert len(matches) == 1
        start, _, length = matches[0]
        assert (
            str(circ_seq_shifted.seq[start : (start + length) % len(circ_seq_shifted)])
            == "aa"
        )


def test_circular_both_sequences():
    """Both sequences circular, sites should still be found."""
    site1 = "GGGaaaCCC"
    site2 = "TTaaaTT"

    circ_x = Dseqrecord("aaaCCCcaGGG", circular=True)
    circ_y = Dseqrecord("aaaTTgaTT", circular=True)

    rec = Recombinase(site1, site2)
    matches = rec.overlap(circ_x, circ_y, limit=None)
    assert len(matches) == 1
    assert matches[0] == (0, 0, 3)


# ---------------------------------------------------------------------------
# Degenerate sequences
# ---------------------------------------------------------------------------


def test_degenerate_site():
    """Sites with IUPAC ambiguity codes should match concrete sequences."""
    site1 = "ATGNNNaaTT"
    site2 = "CCNNaaTTGG"

    seqA = Dseqrecord("xxxATGCCCaaTTyyy")
    seqB = Dseqrecord("xxxCCGGaaTTGGyyy")

    rec = Recombinase(site1, site2)
    matches = rec.overlap(seqA, seqB, limit=None)
    assert len(matches) == 1
    assert matches[0] == (9, 7, 2)


# ---------------------------------------------------------------------------
# Reverse complement
# ---------------------------------------------------------------------------

# TODO: check tests below


def test_reverse_complement_sites():
    """Sites on the reverse strand should be detected."""

    site1 = "ATGCCCTAAaaTT"
    site2 = "AAaaTTTTTTTCCCT"

    rc_site1 = reverse_complement(site1)
    rc_site2 = reverse_complement(site2)
    seqA = Dseqrecord(f"ggg{rc_site1}ggg")
    seqB = Dseqrecord(f"ttt{rc_site2}ttt")

    rec = Recombinase(site1, site2)
    matches = rec.overlap(seqA, seqB, limit=None)
    assert len(matches) == 1
    assert matches[0] == (5, 14, 2)


# ---------------------------------------------------------------------------
# Site annotation
# ---------------------------------------------------------------------------


def test_find_recombinase_sites():
    site1 = "ATGCCCTAAaaTT"
    site2 = "AAaaTTTTTTTCCCT"
    rec = Recombinase(site1, site2, site1_name="s1", site2_name="s2")
    seq = Dseqrecord(f"ggg{site1.upper()}ttt")
    sites = rec.find(seq)
    assert "s1" in sites
    assert sites["s1"] == [SimpleLocation(3, 3 + len(site1), 1)]

    # Works origin-spanning sites
    seq_looped = seq.looped()
    for shift in range(len(seq_looped)):
        seq_looped_shifted = seq_looped.shifted(shift)
        sites = rec.find(seq_looped_shifted)
        assert str(sites["s1"][0].extract(seq_looped_shifted.seq)) == site1.upper()

    # Works with multiple sites
    seq = Dseqrecord(f"ggg{site1.upper()}ttt{site2.upper()}ttt{site1.upper()}ttt")
    sites = rec.find(seq)
    assert len(sites["s1"]) == 2
    assert len(sites["s2"]) == 1


def test_annotate_recombinase_sites():
    site1 = "ATGCCCTAAaaTT"
    site2 = "AAaaTTTTTTTCCCT"
    rec = Recombinase(site1, site2, site1_name="mysite1", site2_name="mysite2")
    seq = Dseqrecord(f"ggg{site1.upper()}ttt")

    annotated_seq = rec.annotate(seq)
    assert len(annotated_seq.features) == 1
    assert annotated_seq.features[0].qualifiers.get("label", []) == ["mysite1"]
    assert annotated_seq.features[0].location == SimpleLocation(3, 3 + len(site1), 1)

    # Does not change original sequence
    assert len(seq.features) == 0


# ---------------------------------------------------------------------------
# recombinase_integration and recombinase_excision
# ---------------------------------------------------------------------------


def test_recombinase_integration():
    """Integration of insert into genome via recombinase sites."""
    site1 = "ATGCCCTAAaaCT"
    site2 = "CAaaTTTTTTTCCCT"

    genome = Dseqrecord(f"cccccc{site1.upper()}aaaaa")
    insert = Dseqrecord(f"{site2.upper()}bbbbb", circular=True)
    rec = Recombinase(site1, site2)

    products = recombinase_integration(genome, [insert], rec)

    assert len(products) == 1
    assert len(products[0].features) == 0
    assert str(products[0].seq) == "ccccccATGCCCTAAAATTTTTTTCCCTbbbbbCAAACTaaaaa"


def test_recombinase_excision():
    """Excision of plasmid from genome with two recombinase sites."""
    site1 = "ATGCCCTAAaaCT"
    site2 = "AAaaTTTTTTTCCCT"

    rec = Recombinase(site1, site2)
    genome = Dseqrecord(f"cccccc{site1.upper()}tttt{site2.upper()}aaaaa")

    products = recombinase_excision_or_inversion(genome, rec)
    assert len(products) == 2
    assert products[0].seq.upper() == Dseq("AACTttttAA".upper())
    assert (
        products[1].seq.upper()
        == Dseq("ccccccATGCCCTAAAATTTTTTTCCCTaaaaa", circular=True).upper()
    )


def test_recombinase_excision_deprecated_alias_warns_and_matches_new_name():
    site1 = "ATGCCCTAAaaCT"
    site2 = "AAaaTTTTTTTCCCT"
    rec = Recombinase(site1, site2)
    genome = Dseqrecord(f"cccccc{site1.upper()}tttt{site2.upper()}aaaaa")

    with pytest.warns(DeprecationWarning, match="recombinase_excision_or_inversion"):
        products_old = recombinase_excision(genome, rec)
    products_new = recombinase_excision_or_inversion(genome, rec)

    assert [p.seq.seguid() for p in products_old] == [
        p.seq.seguid() for p in products_new
    ]


def test_recombinase_integration_excision_reversibility():
    """With same site on both sides (Cre-Lox style), integrate then excise returns originals."""
    site = "ATGaaGTA"

    genome = Dseqrecord(f"cccccc{site.upper()}aaaaa")
    insert = Dseqrecord(f"{site.upper()}bbbbb", circular=True)
    rec = Recombinase(site, site)

    products = recombinase_integration(genome, [insert], rec)
    assert len(products) == 1
    integrated = products[0]

    excised = recombinase_excision_or_inversion(integrated, rec)
    assert len(excised) == 2

    assert excised[1].seq.seguid() == genome.seq.seguid()
    assert excised[0].seq.seguid() == insert.seq.seguid()


def test_recombinase_reverse():
    site1 = "ATGCCCTAAaaTT"
    site2 = "AAaaTTTTTTTCCCT"
    rec = Recombinase(site1, site2)
    rev_rec = rec.get_reverse_recombinase()
    assert rev_rec.site1 == "ATGCCCTAAaaTTTTTTTCCCT"
    assert rev_rec.site2 == "AAaaTT"
    assert rev_rec.site1_name == "site12"
    assert rev_rec.site2_name == "site21"

    # Renaming convention
    for site1_name, site2_name in [
        ("attB", "attP"),
        ("attP", "attB"),
        ("attL", "attR"),
        ("attR", "attL"),
    ]:
        rec = Recombinase(site1, site2, site1_name=site1_name, site2_name=site2_name)
        rev_rec = rec.get_reverse_recombinase()
        if site1_name == "attB":
            expected = ("attR", "attL")
        elif site1_name == "attP":
            expected = ("attL", "attR")
        elif site1_name == "attL":
            expected = ("attP", "attB")
        elif site1_name == "attR":
            expected = ("attB", "attP")

        assert rev_rec.site1_name == expected[0]
        assert rev_rec.site2_name == expected[1]


# ---------------------------------------------------------------------------
# RecombinaseCollection
# ---------------------------------------------------------------------------


def test_recombinase_collection():
    site1 = "AAaaTTC"
    site2 = "CCaaTTC"
    site3 = "GAccACC"
    site4 = "TCccAAC"
    rec1 = Recombinase(site1, site2, site1_name="s1", site2_name="s2")
    rec2 = Recombinase(site3, site4, site1_name="s3", site2_name="s4")
    collection = RecombinaseCollection([rec1, rec2])
    seq1 = Dseqrecord(f"ggg{site1.upper()}ttt{site3.upper()}ttt")
    seq2 = Dseqrecord(f"gggc{site2.upper()}ttt{site4.upper()}ttt")
    assert collection.overlap(seq1, seq2) == [(5, 6, 2), (15, 16, 2)]


def test_recombinase_collection_find():
    # Here the important thing to test is that if
    # two recombinases have the same site name, the find method should return
    # both locations.
    site1 = "AAaaTTC"
    site2 = "CCaaTTC"
    site3 = "GAccACC"
    site4 = "TCccAAC"
    rec1 = Recombinase(site1, site2, site1_name="s1", site2_name="s2")
    rec2 = Recombinase(site3, site4, site1_name="s1", site2_name="s2")
    collection = RecombinaseCollection([rec1, rec2])
    seq = Dseqrecord(f"ggg{site1.upper()}ttt{site3.upper()}ttt")
    assert len(collection.find(seq)["s1"]) == 2


def test_recombinase_collection_annotate():
    # Same as above, checking that two sites with the same name are annotated.
    site1 = "AAaaTTC"
    site2 = "CCaaTTC"
    site3 = "GAccACC"
    site4 = "TCccAAC"
    rec1 = Recombinase(site1, site2, site1_name="s1", site2_name="s2")
    rec2 = Recombinase(site3, site4, site1_name="s1", site2_name="s2")
    collection = RecombinaseCollection([rec1, rec2])
    seq = Dseqrecord(f"ggg{site1.upper()}ttt{site3.upper()}ttt")
    annotated_seq = collection.annotate(seq)
    assert len(annotated_seq.features) == 2
    assert annotated_seq.features[0].qualifiers.get("label", []) == ["s1"]
    assert annotated_seq.features[1].qualifiers.get("label", []) == ["s1"]

    # It does not re-annotate if the sequence is already annotated.
    annotated_seq = collection.annotate(annotated_seq)
    assert len(annotated_seq.features) == 2


def test_recombinase_collection_init_errors():
    # Check that the init method raises errors for invalid inputs.
    with pytest.raises(ValueError):
        RecombinaseCollection(None)
    with pytest.raises(ValueError):
        RecombinaseCollection([])
    with pytest.raises(ValueError):
        RecombinaseCollection([Recombinase("AAaaTTC", "CCaaTTC"), "blah"])


def test_recombinase_collection_in_assembly_functions():
    site1 = "AAaaTTC"
    site2 = "CCaaTTC"
    site3 = "GAccACC"
    site4 = "TCccAAC"
    rec1 = Recombinase(site1, site2, site1_name="s1", site2_name="s2")
    rec2 = Recombinase(site3, site4, site1_name="s1", site2_name="s2")
    collection = RecombinaseCollection([rec1, rec2])
    seq = Dseqrecord(f"ggg{site1.upper()}ttt{site3.upper()}ttt")
    seq2 = Dseqrecord(f"ggg{site2.upper()}ttt{site4.upper()}ttt")
    products = recombinase_integration(seq, [seq2], collection)
    assert len(products) == 1


def test_inversion_recombinase():
    attL = "AAAttGCGC"
    attR = "TATttCCAA"
    attLR = "AAAttCCAA"
    attRL = "TATttGCGC"

    seq1 = Dseqrecord("ggg" + attL + "cca" + reverse_complement(attR) + "tttt")
    rec = Recombinase(attL, attR)

    prods = recombinase_excision_or_inversion(seq1, rec)
    expected = (
        "ggg" + attLR + reverse_complement("cca") + reverse_complement(attRL) + "tttt"
    )
    assert len(prods) == 1
    assert str(prods[0].seq) == expected
    prods[0].validate_history()

    # Same results with assembly function
    assert prods == recombinase_assembly([seq1], rec)

    seq1 = seq1.looped()
    prods = recombinase_excision_or_inversion(seq1, rec)
    assert len(prods) == 1
    assert prods[0].seq.seguid() == Dseq(expected, circular=True).seguid()
    # Same results with assembly function
    assert prods == recombinase_assembly([seq1], rec)
    prods[0].validate_history()


def test_recombinase_assembly():
    attL = "AAAttGCGC"
    attR = "TATttCCAA"
    attLR = "AAAttCCAA"
    attRL = "TATttGCGC"

    seq1 = Dseqrecord(f"ggg{attL}aaa")
    seq2 = Dseqrecord(f"ccc{attR}ttt")
    rec = Recombinase(attL, attR)

    prods = recombinase_assembly([seq1, seq2], rec)
    assert len(prods) == 2
    assert str(prods[0].seq) == f"ggg{attLR}ttt"
    assert str(prods[1].seq) == f"ccc{attRL}aaa"

    # Circular sequences should merge
    seq1 = seq1.looped()
    seq2 = seq2.looped()
    prods = recombinase_assembly([seq1, seq2], rec)
    assert len(prods) == 1
    assert (
        prods[0].seguid() == Dseq(f"ggg{attLR}tttccc{attRL}aaa", circular=True).seguid()
    )
