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
from pydna.snapgene_history_parser import parse_snapgene_history

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
        self.assertEqual(len(TEST_FILES), 42)

    def test_correctly_parsed(self):
        for file in TEST_FILES:
            try:
                parse_snapgene_history(file)
            except NotImplementedError as e:
                filename = os.path.basename(file)
                if filename not in METHOD_NOT_SUPPORTED:
                    raise AssertionError(f"File {filename} not supported") from e
            except ValueError as e:
                filename = os.path.basename(file)
                if (
                    filename not in EXPECTED_VALUE_ERROR
                    or "No product found for expected SEGUID" not in str(e)
                ):
                    raise AssertionError(f"File {filename} not supported") from e
