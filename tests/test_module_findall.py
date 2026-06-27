#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from pydna.findall import findall

<<<<<<< HEAD
=======

>>>>>>> 3052ee83 (added tests)
def test_arguments():
    with pytest.raises(TypeError, match="needle and haystack must be strings"):
        findall(1, 2)

    with pytest.raises(TypeError, match="needle and haystack must be strings"):
        findall("ACGT", 2)

    with pytest.raises(TypeError, match="needle and haystack must be strings"):
        findall(1, "ACGT")

    with pytest.raises(ValueError, match="needle must not be empty"):
        findall("", "ACGT")

    with pytest.raises(ValueError, match="haystack must not be empty"):
        findall("ACGT", "")

    with pytest.raises(TypeError, match="max_edits must be an integer"):
        findall("ACGT", "ACGT", max_edits=1.5)

    with pytest.raises(TypeError, match="max_edits must be an integer"):
        findall("ACGT", "ACGT", max_edits="1")

    with pytest.raises(ValueError, match="max_edits must be a positive integer"):
        findall("ACGT", "ACGT", max_edits=-1)

    assert findall("ACGT", "ACG") == []
    assert findall("ACGT", "ACG", circular=True) == []


def test_find_with_deletion():
    # ACGT-ACGT
    # ||||-||||
    # ACGTTACGT
    needle = "ACGTACGT"
    haystack = "TTTACGTTACGTAAA"

    results = findall(needle, haystack, 1)

    assert results == [
        {
            "distance": 1,
            "start": 3,
            "stop": 12,
            "cigar": "4=1D4=",
            "alignment": "ACGT-ACGT\n||||-||||\nACGTTACGT",
        }
    ]


def test_find_with_insert():
    # ACGTTACGT
    # ||||-||||
    # ACGT-ACGT
    needle = "ACGTTACGT"
    haystack = "TTTACGTACGTAAA"

    results = findall(needle, haystack, 1, circular=False)

    assert results == [
        {
            "distance": 1,
            "start": 3,
            "stop": 11,
            "cigar": "4=1I4=",
            "alignment": "ACGTTACGT\n||||-||||\nACGT-ACGT",
        }
    ]


def test_find_with_sub():
    # ACGTTACGT
    # ||||.||||
    # ACGTCACGT
    needle = "ACGTTACGT"
    haystack = "TTTACGTCACGTAAA"

    results = findall(needle, haystack, 1, circular=False)

    assert results == [
        {
            "distance": 1,
            "start": 3,
            "stop": 12,
            "cigar": "4=1X4=",
            "alignment": "ACGTTACGT\n||||.||||\nACGTCACGT",
        }
    ]


def test_find_with_insert_deletion():
    # ACGT-ACAGTTT
    # ||||-||-||||
    # ACGTTAC-GTTT
    needle = "ACGTACAGTTT"
    haystack = "TTTACGTTACGTTTAAA"

    results = findall(needle, haystack, 2, circular=False)

    assert results == [
        {
            "distance": 2,
            "start": 3,
            "stop": 14,
            "cigar": "4=1D2=1I4=",
            "alignment": "ACGT-ACAGTTT\n||||-||-||||\nACGTTAC-GTTT",
        }
    ]


def test_find_across_circular_origin():
    needle = "TAAC"
    haystack = "ACGTTA"

    results = findall(
        needle,
        haystack,
        max_edits=0,
        circular=True,
    )

    assert {
        "distance": 0,
        "start": 4,
        "stop": 2,
        "cigar": "4=",
        "alignment": "TAAC\n||||\nTAAC",
    } in results


def test_find_two():
    needle = "TAAC"
    haystack = "TAACACGTTAACTA"

    results = findall(needle, haystack, max_edits=0)

    assert results == [
        {
            "distance": 0,
            "start": 0,
            "stop": 4,
            "cigar": "4=",
            "alignment": "TAAC\n||||\nTAAC",
        },
        {
            "distance": 0,
            "start": 8,
            "stop": 12,
            "cigar": "4=",
            "alignment": "TAAC\n||||\nTAAC",
        },
    ]


def test_find_overlap():
    needle = "TATA"
    haystack = "TATATA"

    results = findall(needle, haystack, max_edits=0)

    assert results == [
        {
            "distance": 0,
            "start": 0,
            "stop": 4,
            "cigar": "4=",
            "alignment": "TATA\n||||\nTATA",
        },
        {
            "distance": 0,
            "start": 2,
            "stop": 6,
            "cigar": "4=",
            "alignment": "TATA\n||||\nTATA",
        },
    ]


def test_iupac_n_matches_all_bases():
    results = findall("ANT", "AAT AGT ACT ATT".replace(" ", ""), max_edits=0)

    assert [(r["start"], r["stop"], r["distance"], r["cigar"]) for r in results] == [
        (0, 3, 0, "3="),
        (3, 6, 0, "3="),
        (6, 9, 0, "3="),
        (9, 12, 0, "3="),
    ]


def test_iupac_purine_r_matches_a_or_g():
    results = findall("ART", "AATAGTACTATT", max_edits=0)

    assert [(r["start"], r["stop"]) for r in results] == [
        (0, 3),  # AAT
        (3, 6),  # AGT
    ]


def test_iupac_ambiguity_counts_as_zero_edits():
    (result,) = findall("ANT", "AGT", max_edits=0)

    assert result["distance"] == 0
    assert result["cigar"] == "3="
    assert result["alignment"] == "ANT\n|||\nAGT"


def test_iupac_still_allows_real_edits():
    (result,) = findall("ANT", "AG", max_edits=1)

    assert result["distance"] == 1
    assert result["start"] == 0
    assert result["stop"] == 2


def test_iupac_circular_match():
    results = findall("NTA", "AACGT", max_edits=0, circular=True)

    assert {
        "distance": 0,
        "start": 3,
        "stop": 1,
        "cigar": "3=",
        "alignment": "NTA\n|||\nGTA",
    } in results
