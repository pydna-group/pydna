# -*- coding: utf-8 -*-
from pydna.cre_lox import cre_loxP_overlap, annotate_loxP_sites
from pydna.assembly2 import cre_lox_excision_or_inversion
from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord
from pydna.parsers import parse
import os
import unittest

# ATAACTTCGTATANNNTANNNTATACGAAGTTAT
loxP_sequence = "ATAACTTCGTATAGCATACATTATACGAAGTTAT"
loxP_sequence2 = "ATAACTTCGTATAGTTTACATTATACGAAGTTAT"

lox66_sequence = "ATAACTTCGTATAGTTTACATTATACGAACGGTA"
lox71_sequence = "TACCGTTCGTATAGTTTACATTATACGAAGTTAT"

test_files = os.path.join(os.path.dirname(__file__))


class TestCreLox(unittest.TestCase):
    def test_cre_loxP_overlap(self):
        # Works with consensus sequence and other sequence
        for loxP in [loxP_sequence, loxP_sequence2]:
            seqA = Dseqrecord("aa" + loxP + "acgt")
            seqB = Dseqrecord("ccc" + loxP + "acgt")
            self.assertEqual(cre_loxP_overlap(seqA, seqB), [(15, 16, 8)])
            self.assertEqual(
                cre_loxP_overlap(seqA.reverse_complement(), seqB.reverse_complement()),
                [(17, 17, 8)],
            )
            self.assertEqual(cre_loxP_overlap(seqA, seqB.reverse_complement()), [])
            self.assertEqual(cre_loxP_overlap(seqA.reverse_complement(), seqB), [])

        # It does not mix different loxP sites, even if they match the consensus
        seqA = Dseqrecord("aa" + loxP_sequence + "acgt")
        seqB = Dseqrecord("ccc" + loxP_sequence2 + "acgt")
        self.assertEqual(cre_loxP_overlap(seqA, seqB), [])
        self.assertEqual(
            cre_loxP_overlap(seqA.reverse_complement(), seqB.reverse_complement()), []
        )
        self.assertEqual(cre_loxP_overlap(seqA, seqB.reverse_complement()), [])
        self.assertEqual(cre_loxP_overlap(seqA.reverse_complement(), seqB), [])

        # Works for the lox66 and lox71
        seqA = Dseqrecord("aa" + lox66_sequence + "acgt")
        seqB = Dseqrecord("ccc" + lox71_sequence + "acgt")
        self.assertEqual(cre_loxP_overlap(seqA, seqB), [(15, 16, 8)])
        self.assertEqual(
            cre_loxP_overlap(seqA.reverse_complement(), seqB.reverse_complement()),
            [(17, 17, 8)],
        )
        self.assertEqual(cre_loxP_overlap(seqA, seqB.reverse_complement()), [])
        self.assertEqual(cre_loxP_overlap(seqA.reverse_complement(), seqB), [])

    def test_annotate_loxP_sites(self):
        seq = parse(os.path.join(test_files, "all_cre_lox.gb"))[0]

        # Make a copy without features
        seq_no_features = Dseqrecord(seq.seq)
        seq_annotated = annotate_loxP_sites(seq_no_features)
        self.assertEqual(len(seq_annotated.features), 6)

        # Check that re-annotating does not add more features
        seq_annotated2 = annotate_loxP_sites(seq)
        self.assertEqual(len(seq_annotated2.features), 6)

    def test_cre_lox_excision_or_inversion(self):
        seq1 = Dseqrecord(
            f"cccc{loxP_sequence}ttggca{reverse_complement(loxP_sequence)}aaaaa"
        )
        prods = cre_lox_excision_or_inversion(seq1)
        assert len(prods) == 1
        expected = f"cccc{loxP_sequence}tgccaa{reverse_complement(loxP_sequence)}aaaaa"
        assert str(prods[0].seq) == expected
        seq1 = seq1.looped()
        prods = cre_lox_excision_or_inversion(seq1)
        assert len(prods) == 1
        assert prods[0].seq.seguid() == Dseqrecord(expected, circular=True).seguid()
