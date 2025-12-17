#!/usr/bin/env python
# -*- coding: utf-8 -*-

from textwrap import dedent

from pydna.dseq import Dseq
from pydna.alphabet import (
    # dictionaries
    basepair_dict,
    annealing_dict,
    annealing_dict_w_holes,
    complement_dict_for_dscode,

    # translation tables
    complement_table_for_dscode,
    dscode_to_watson_table,
    dscode_to_crick_table,
    dscode_to_watson_tail_table,
    dscode_to_crick_tail_table,
    dscode_to_full_sequence_table,

    # letter sets
    ds_letters,
    ss_letters_watson,
    ss_letters_crick,

    # regex-related
    iupac_compl_regex,
    regex_ss_melt_factory,
    regex_ds_melt_factory,

    # data structures
    DseqParts,

    # helpers
    get_parts,
    dsbreaks,
    representation_tuple,
    anneal_strands,
)


# ---------------------------------------------------------------------------
# Translation tables
# ---------------------------------------------------------------------------

def test_translation_tables():

    assert "GATC".translate(complement_table_for_dscode) == "CTAG"
    assert "GATC".translate(dscode_to_watson_table) == "GATC"
    assert "GATC".translate(dscode_to_crick_table) == "CTAG"

    assert "G".translate(dscode_to_watson_tail_table) == "P"
    assert "C".translate(dscode_to_crick_tail_table) == "J"

    assert "P".translate(dscode_to_full_sequence_table) == "G"
    assert "Q".translate(dscode_to_full_sequence_table) == "G"


# ---------------------------------------------------------------------------
# Alphabet letter sets
# ---------------------------------------------------------------------------

def test_letters():

    assert ds_letters == "GATCUORYMKSWHBVDNgatcuorymkswhbvdn"
    assert ds_letters == "GATCUORYMKSWHBVDNgatcuorymkswhbvdn"
    assert ss_letters_watson == "PEXI$pexi$"
    assert ss_letters_crick == "QFZJ%qfzj%"


# ---------------------------------------------------------------------------
# Core dictionaries
# ---------------------------------------------------------------------------

def test_basepair_dict():

    assert basepair_dict["G", "C"] == "G"
    assert basepair_dict["C", "G"] == "C"
    assert basepair_dict["A", "T"] == "A"
    assert basepair_dict["T", "A"] == "T"

    # mixed case must work
    assert basepair_dict["g", "c"] == "g"
    assert basepair_dict["G", "c"] == "G"


def test_complement_dict_for_dscode():

    assert complement_dict_for_dscode["G"] == "C"
    assert complement_dict_for_dscode["C"] == "G"
    assert complement_dict_for_dscode["A"] in {"T", "U"}


def test_annealing_dicts():

    # perfect annealing
    assert anneal_strands("TTA", "AAT"[::-1]) is True
    assert anneal_strands("UUA", "AAT"[::-1]) is True

    # imperfect annealing
    assert anneal_strands("TG", "AA") is False

    # hole-tolerant dict must be a superset
    assert len(annealing_dict_w_holes) >= len(annealing_dict)


# ---------------------------------------------------------------------------
# IUPAC complement regex
# ---------------------------------------------------------------------------

def test_iupac_compl_regex():

    assert iupac_compl_regex["A"] == "(?:T|U)"
    assert iupac_compl_regex["N"] == "(?:A|G|C|T|N)"
    assert set(iupac_compl_regex) == set("ACGTURYSWKMBDHVN")


# ---------------------------------------------------------------------------
# regex_ss_melt_factory
# ---------------------------------------------------------------------------

def test_regex_ss_melt_factory():

    regex = regex_ss_melt_factory(3)
    assert regex.groups == 8

    s = Dseq("GFTTAJA")
    assert repr(s) == dedent("""\
        Dseq(-7)
        G TTA A
        CTAATGT
    """).strip()

    m = regex.search(s._data)
    assert m is not None
    assert m.groupdict()["watson"] == b"TTA"


# ---------------------------------------------------------------------------
# regex_ds_melt_factory
# ---------------------------------------------------------------------------

def test_regex_ds_melt_factory():

    regex = regex_ds_melt_factory(3)
    assert regex.groups == 8

    s = Dseq("aaaGFTTAIAttt")
    assert repr(s) == dedent("""\
        Dseq(-13)
        aaaG TTACAttt
        tttCTAAT Taaa
    """).strip()

    m = regex.search(s._data)
    assert m is not None
    assert m.groupdict()["crick"] == b"TTA"


# ---------------------------------------------------------------------------
# DseqParts and helpers
# ---------------------------------------------------------------------------

def test_get_parts_and_DseqParts():

    parts = get_parts("eeATCGuggCCGgg")

    assert isinstance(parts, DseqParts)

    (
        sticky_left5,
        sticky_left3,
        middle,
        sticky_right3,
        sticky_right5,
        single_watson,
        single_crick,
    ) = parts

    assert sticky_left5 == "ee"
    assert middle == "ATCGuggCCGgg"
    assert single_watson == ""
    assert single_crick == ""


def test_representation_tuple():

    w, c = representation_tuple("GATC")
    assert w == "GATC"
    assert c == "CTAG"


def test_dsbreaks():

    breaks = dsbreaks("GATPFTAA")
    assert len(breaks) == 1
    assert "[0:8]" in breaks[0]
