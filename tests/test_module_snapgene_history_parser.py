#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Tests for pydna.snapgene_history_parser."""

import glob
from unittest import TestCase
import os
from pydna.opencloning_models import AddgeneIdSource
from pydna.snapgene_history_parser import (
    parse_snapgene_history,
    SnapgeneHistoryParserWarning,
)
import warnings

TEST_FILES = glob.glob(
    os.path.join(os.path.dirname(__file__), "snapgene_history_files", "*.dna")
)[::-1]

METHOD_NOT_SUPPORTED = [
    "topo_ta_cloning.dna",
    "topo_directional_cloning.dna",
    "gc_cloning_including_overhangs.dna",
    "gc_cloning_overhangs_after.dna",
    "ta_cloning_including_overhangs.dna",
    "ta_cloning_overhangs_after.dna",
    "topo_blunt_cloning.dna",
    "delete_restriction_fragment_without_compatible_overhangs.dna",
    "pcr_modifying_ends_insert_R.dna",
    "destroy_restriction_site.dna",
]

EXPECTED_VALUE_ERROR = [
    # The problem here is that it's a linear blunt ligation
    # of several sequences. pydna does not return this because
    # it would circularize.
    "blunt_linear_ligation.dna"
]


class TestSnapgeneHistoryParser(TestCase):

    def test_files_exist(self):
        self.assertEqual(len(TEST_FILES), 49)

    def test_correctly_parsed(self):
        seqr_dict = {}
        warning_dict = {}
        for file in TEST_FILES:
            filename = os.path.basename(file)
            try:

                with warnings.catch_warnings(record=True) as wlist:
                    warnings.simplefilter("always")
                    seqr_dict[filename] = parse_snapgene_history(file)
                    warning_dict[filename] = [
                        str(w.message)
                        for w in wlist
                        if issubclass(w.category, SnapgeneHistoryParserWarning)
                    ]
            except NotImplementedError as e:
                if filename not in METHOD_NOT_SUPPORTED:
                    raise AssertionError(f"File {filename} not supported") from e
            except ValueError as e:
                if (
                    filename not in EXPECTED_VALUE_ERROR
                    or "No product found for expected SEGUID" not in str(e)
                ):
                    raise AssertionError(f"File {filename} not supported") from e

        # Check special cases
        seqr = seqr_dict["import_addgene_then_clone.dna"]
        self.assertIsInstance(seqr.source.input[0].sequence.source, AddgeneIdSource)
        # Can't do this anymore, (see _source_from_metadata)
        # seqr = seqr_dict["import_ncbi.dna"]
        # self.assertIsInstance(seqr.source, NCBISequenceSource)
        seqr = seqr_dict["import_addgene.dna"]
        self.assertIsInstance(seqr.source, AddgeneIdSource)

        # Check warnings
        self.assertEqual(
            warning_dict["circularize_then_linearize_without_enzyme.dna"],
            ["Stopped at linearize operation without enzymes"],
        )
        del warning_dict["circularize_then_linearize_without_enzyme.dna"]
        self.assertEqual(
            warning_dict["circularize.dna"],
            ["Stopped at change topology operation"],
        )
        del warning_dict["circularize.dna"]

        # Manual editing warnings
        for f in TEST_FILES:
            basename = os.path.basename(f)
            if basename.startswith("manual_"):
                self.assertEqual(
                    warning_dict[basename],
                    ["Manual editing of sequences not supported"],
                )
                del warning_dict[basename]

        # Not other warning should remain
        for key in warning_dict:
            self.assertEqual(warning_dict[key], [])
