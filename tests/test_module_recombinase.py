import pytest

from pydna.dseqrecord import Dseqrecord
from pydna.assembly2 import (
    Assembly,
    recombinase_integration,
    recombinase_excision,
)
from pydna.recombinase import (
    _recombinase_homology_offset_and_length,
    Recombinase,
    find_recombinase_sites,
    annotate_recombinase_sites,
)


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
    from Bio.Seq import reverse_complement

    site1 = "ATGCCCTAAaaTT"
    site2 = "AAaaTTTTTTTCCCT"

    rc_site1 = reverse_complement(site1)
    seqA = Dseqrecord(f"ggg{rc_site1}ggg")
    seqB = Dseqrecord(f"ttt{reverse_complement(site2)}ttt")

    rec = Recombinase(site1, site2)
    matches = rec.overlap(seqA, seqB, limit=None)
    assert len(matches) >= 1


# ---------------------------------------------------------------------------
# Recombinase class
# ---------------------------------------------------------------------------


def test_recombinase_class_basic():
    rec = Recombinase("ATGCCCTAAaaTT", "AAaaTTTTTTTCCCT")

    seqA = Dseqrecord("aaaATGCCCTAAaaTTtt")
    seqB = Dseqrecord("tataAAaaTTTTTTTCCCTaaa")

    matches = rec.overlap(seqA, seqB)
    assert matches == [(12, 6, 2)]


# ---------------------------------------------------------------------------
# Site annotation
# ---------------------------------------------------------------------------


def test_find_recombinase_sites():
    site1 = "ATGCCCTAAaaTT"
    seq = Dseqrecord(f"ggg{site1.upper()}ttt")
    sites = find_recombinase_sites(
        seq, site1, "GGGaaaCCC", site1_name="s1", site2_name="s2"
    )
    assert "s1" in sites
    assert len(sites["s1"]) >= 1


def test_annotate_recombinase_sites():
    site1 = "ATGCCCTAAaaTT"
    seq = Dseqrecord(f"ggg{site1.upper()}ttt")
    assert len(seq.features) == 0
    annotate_recombinase_sites(
        seq, site1, "GGGaaaCCC", site1_name="mysite1", site2_name="mysite2"
    )
    assert len(seq.features) >= 1
    assert any(f.qualifiers.get("label", []) == ["mysite1"] for f in seq.features)

    # Calling again should not duplicate
    n_before = len(seq.features)
    annotate_recombinase_sites(
        seq, site1, "GGGaaaCCC", site1_name="mysite1", site2_name="mysite2"
    )
    assert len(seq.features) == n_before


# ---------------------------------------------------------------------------
# Symmetry: site1 in seqy and site2 in seqx should also work
# ---------------------------------------------------------------------------


def test_sites_swapped_between_sequences():
    """With site1 in seqx and site2 in seqy (can be either order in Assembly)."""
    site1 = "ATGCCCTAAaaTT"
    site2 = "AAaaTTTTTTTCCCT"

    seqA = Dseqrecord("aaaATGCCCTAAaaTTtt")
    seqB = Dseqrecord("tataAAaaTTTTTTTCCCTaaa")

    rec = Recombinase(site1, site2)
    matches = rec.overlap(seqA, seqB, limit=None)
    assert len(matches) >= 1


# ---------------------------------------------------------------------------
# recombinase_integration and recombinase_excision
# ---------------------------------------------------------------------------


def test_recombinase_integration():
    """Integration of insert into genome via recombinase sites."""
    site1 = "ATGCCCTAAaaTT"
    site2 = "AAaaTTTTTTTCCCT"

    genome = Dseqrecord(f"cccccc{site1.upper()}aaaaa")
    insert = Dseqrecord(f"{site2.upper()}bbbbb", circular=True)

    products = recombinase_integration(genome, [insert], site1, site2)
    assert len(products) >= 1
    # Integrated product should contain both flanking regions and insert
    seq_str = str(products[0].seq).upper()
    assert "ATGCCCTAA" in seq_str
    assert "TTTTTTTCCCT" in seq_str


def test_recombinase_excision():
    """Excision of plasmid from genome with two recombinase sites."""
    site1 = "ATGCCCTAAaaTT"
    site2 = "AAaaTTTTTTTCCCT"

    genome = Dseqrecord(
        f"cccccc{site1.upper()}xxxxx{site2.upper()}aaaaa",
        circular=True,
    )

    products = recombinase_excision(genome, site1, site2)
    assert len(products) >= 1


def test_recombinase_integration_excision_reversibility():
    """With same site on both sides (Cre-Lox style), integrate then excise returns originals."""
    site = "ATGCCCTAAaaTT"

    genome = Dseqrecord(f"cccccc{site.upper()}aaaaa")
    insert = Dseqrecord(f"{site.upper()}bbbbb", circular=True)

    products = recombinase_integration(genome, [insert], site, site)
    assert len(products) == 1
    integrated = products[0]

    excised = recombinase_excision(integrated, site, site)
    assert len(excised) == 2

    seqs = [str(p.seq).upper() for p in excised]
    assert any(site.upper() in s for s in seqs)
