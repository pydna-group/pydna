import regex
import pytest

from pydna.fuzzy import fuzzymatch


def test_fuzzy_match_alignment():
    needle = "GGTGCTGGTCGCCGCTT"
    haystack = "AGTGCTGGTCTCCGGCTTCGGTGCTGGTCGCCGCTTC"

    results = fuzzymatch(needle, haystack, 15, 1, 1, 1)
    first = next(result for result in results if result["start"] == 1)
    second = next(result for result in results if result["start"] == 19)

    assert first == {
        "fuzzy_counts": (1, 1, 1),
        "start": 1,
        "stop": 18,
        "cigar": "1=1I8=1X2=1D4=",
        "alignment": (
            "GGTGCTGGTCGCC-GCTT\n" "| ||||||||.|| ||||\n" "G-TGCTGGTCTCCGGCTT"
        ),
    }
    assert second["fuzzy_counts"] == (0, 0, 0)
    assert second["cigar"] == "17="


def test_limit_and_literal_needle():
    assert len(fuzzymatch("A", "AAA", limit=2)) == 2
    assert fuzzymatch("A.T", "xxA.Txx") == [
        {
            "fuzzy_counts": (0, 0, 0),
            "start": 2,
            "stop": 5,
            "cigar": "3=",
            "alignment": "A.T\n|||\nA.T",
        }
    ]


def test_overlapping_matches():
    results = fuzzymatch("TATA", "TATATA")

    assert [(result["start"], result["stop"]) for result in results] == [
        (0, 4),
        (2, 6),
    ]


@pytest.mark.parametrize(
    ("kwargs", "error"),
    [
        ({"limit": -1}, ValueError),
        ({"limit": 1.5}, TypeError),
        ({"substitutions": -1}, ValueError),
        ({"insertions": 1.5}, TypeError),
    ],
)
def test_invalid_limits(kwargs, error):
    with pytest.raises(error):
        fuzzymatch("A", "A", **kwargs)


def test_case_sensitive_flag():
    assert fuzzymatch("a", "A", flags=regex.ENHANCEMATCH) == []
