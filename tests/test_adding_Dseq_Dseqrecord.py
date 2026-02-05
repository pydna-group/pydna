#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord

biopython_Seq = Seq("AT")
mutable_biopython_Seq = MutableSeq("AT")
biopython_SeqRecord = SeqRecord(biopython_Seq)


def test_adding_Dseq_to_strings_both_sides():
    """
    Strings can be added on both sides of a Dseq.
    The string will be cast as a Dseq with a single watson strand.
    After that, it will be treated in the same way.

    The joining of Dseqs can be understood as forming one or two phosphodiester
    bonds with compatible base pairing (or no base pairing in case of single
    strand sequences).

    ::

        g    tc    gtc
        |  + ||  = |||
        c    ag    cag

        g    tc    gtc
        |  +  |  = |||
        ca    g    cag

        gt    c    gtc
        |  +  |  = |||
        c    ag    cag

    A special case is the joining of two single stranded Dseqs.
    For single stranded sequences, either a phosphodiester bond is formed
    or a double stranded sequence depending on the orientation and
    complementarity of the sequences.

    ::

        gt + ca    gtca
        ||   ||  = ||||

        gt         gt
        || + ||  = ||
             ca    ca

        ||   ||  = ||||
        gt + ca    gtca


    """

    assert "gatc" + Dseq("GT") + "gatc" == Dseq.from_representation(
        """
    Dseq(-10)
    gatcGTgatc
        CA
    """
    )
    assert Dseq("pexi") + Dseq("GT") + Dseq("pexi") == Dseq.from_representation(
        """
    Dseq(-10)
    gatcGTgatc
        CA
    """
    )
    assert "gatc" + Dseq("PEXI") + "gatc" == Dseq.from_representation(
        """
    Dseq(-12)
    gatcGATCgatc
    ||||||||||||
    """
    )
    assert "gatc" + Dseq("QFZJ") == Dseq.from_representation(
        """
    Dseq(-4)
    gatc
    ctag
    """
    )
    assert Dseq("pexi") + Dseq("QFZJ") == Dseq.from_representation(
        """
    Dseq(-4)
    gatc
    ctag
    """
    )
    assert Dseq("QFZJ") + "gatc" == Dseq.from_representation(
        """
    Dseq(-4)
    GATC
    CTAG
    """
    )
    assert Dseq("QFZJ") + Dseq("pexi") == Dseq.from_representation(
        """
    Dseq(-4)
    GATC
    CTAG
    """
    )
    assert Dseq("QFZJ") + Dseq("qfzj") == Dseq.from_representation(
        """
    Dseq(-8)
    ||||||||
    CTAGctag
    """
    )

    assert Dseq("GP") + Dseq("QT") == Dseq("GGT")
    assert Dseq("GP") + Dseq("qT") == Dseq("GGT")
    assert Dseq("Gp") + Dseq("QT") == Dseq("GgT")
    assert Dseq("Gp") + Dseq("qT") == Dseq("GgT")

    assert Dseq("GP") + Dseq("Q") == Dseq("GG")
    assert Dseq("GP") + Dseq("q") == Dseq("GG")
    assert Dseq("Gp") + Dseq("Q") == Dseq("Gg")
    assert Dseq("Gp") + Dseq("q") == Dseq("Gg")

    with pytest.raises(TypeError) as e:
        Dseq("GP") + Dseq("G")
        assert e.match("sticky ends not compatible!")
    with pytest.raises(TypeError) as e:
        Dseq("GP") + Dseq("A")
        assert e.match("sticky ends not compatible!")
    with pytest.raises(TypeError) as e:
        Dseq("GP") + Dseq("A")
        assert e.match("sticky ends not compatible!")
    with pytest.raises(TypeError) as e:
        Dseq("GP") + Dseq("C")
        assert e.match("sticky ends not compatible!")

    assert Dseq("GP") + Dseq("QT") == Dseq("GGT")
    with pytest.raises(TypeError) as e:
        Dseq("GP") + Dseq("FT")
        assert e.match("sticky ends not compatible!")
    with pytest.raises(TypeError) as e:
        Dseq("GP") + Dseq("ZT")
        assert e.match("sticky ends not compatible!")
    with pytest.raises(TypeError) as e:
        Dseq("GP") + Dseq("JT")
        assert e.match("sticky ends not compatible!")


def test_adding_strings_to_Dseqrecord_on_both_sides():
    """ """
    assert ("aaa" + Dseqrecord("GT") + "aaa").seq == Dseq.from_representation(
        """
    Dseq(-8)
    aaaGTaaa
       CA
    """
    )
    assert ("gatc" + Dseqrecord("PI") + "gatc").seq == Dseq.from_representation(
        """
    Dseq(-10)
    gatcGCgatc
    ||||||||||
    """
    )
    assert ("gatc" + Dseqrecord("QFZJ")).seq == Dseq.from_representation(
        """
    Dseq(-4)
    gatc
    ctag
    """
    )
    assert (Dseqrecord("PEXI") + Dseqrecord("PEXI")).seq == Dseq.from_representation(
        """
    Dseq(-8)
    GATCGATC
    ||||||||
    """
    )
    assert (Dseqrecord("QFZJ") + Dseqrecord("QFZJ")).seq == Dseq.from_representation(
        """
    Dseq(-8)
    ||||||||
    CTAGCTAG
    """
    )

    with pytest.raises(TypeError) as e:
        assert (Dseqrecord("QFZJ") + "aaa").seq
        assert e.match("sticky ends not compatible!")
