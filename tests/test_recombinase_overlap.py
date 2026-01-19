import pytest

from pydna.dseqrecord import Dseqrecord
from pydna.assembly2 import (
    make_recombinase_algorithm,
    _recombinase_homology_offset_and_length,
    Assembly,
)


def test_recombinase_homology_offset_and_length_basic():
    """
    The helper should correctly identify the lowercase homology region.
    """
    attB = "ATGCCCTAAaaTT"         # lowercase 'aa' at positions 9–10
    attP = "AAaaTTTTTTTCCCT"       # lowercase 'aa' at positions 2–3

    offB, lenB = _recombinase_homology_offset_and_length(attB)
    offP, lenP = _recombinase_homology_offset_and_length(attP)

    assert (offB, lenB) == (9, 2)
    assert (offP, lenP) == (2, 2)


def test_recombinase_homology_offset_and_length_no_lowercase_raises():
    """
    If there is no lowercase region in the site, an error should be raised.
    """
    with pytest.raises(ValueError):
        _recombinase_homology_offset_and_length("ATGCCCTAAATT")  # all uppercase


def test_make_recombinase_algorithm_single_match_minimal_example():
    """
    For the example in the issue, the algorithm should return [(12, 6, 2)].
    """
    attB = "ATGCCCTAAaaTT"
    attP = "AAaaTTTTTTTCCCT"

    seqA = Dseqrecord("aaaATGCCCTAAaaTTtt")
    seqB = Dseqrecord("tataAAaaTTTTTTTCCCTaaa")

    algo = make_recombinase_algorithm(attB, attP)

    matches = algo(seqA, seqB, limit=1)
    assert matches == [(12, 6, 2)]


def test_make_recombinase_algorithm_no_site_returns_empty():
    """
    If either attB or attP is missing from the fragments, no overlaps should be found.
    """
    attB = "ATGCCCTAAaaTT"
    attP = "AAaaTTTTTTTCCCT"

    seqA = Dseqrecord("aaaATGCCCTAAaaTTtt")
    # Break attP homology in seqB
    seqB = Dseqrecord("tataAAxxTTTTTTTCCCTaaa")

    algo = make_recombinase_algorithm(attB, attP)

    matches = algo(seqA, seqB, limit=1)
    assert matches == []


def test_recombinase_algorithm_integration_with_assembly():
    """
    Integration test: the recombinase algorithm should work inside Assembly.
    We mainly check that Assembly can build at least one product using this
    algorithm in the minimal example.
    """
    attB = "ATGCCCTAAaaTT"
    attP = "AAaaTTTTTTTCCCT"

    seqA = Dseqrecord("aaaATGCCCTAAaaTTtt")
    seqB = Dseqrecord("tataAAaaTTTTTTTCCCTaaa")

    algo = make_recombinase_algorithm(attB, attP)

    asm = Assembly([seqA, seqB], algorithm=algo, use_fragment_order=False)
    products = asm.assemble_linear()

    # There should be at least one assembled product
    assert len(products) >= 1

    # Optionally: check something about the sequence length or content
    # (keeping it loose for now since Assembly semantics may evolve)
    seq_str = str(products[0].seq)
    assert "ATGCCCTAA" in seq_str
    assert "TTTTTTTCCCT" in seq_str
