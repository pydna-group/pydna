#!/usr/bin/env python3
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Tests for pydna.snapgene_history_parser."""

import glob
import os
import warnings

import pytest
from pydna.opencloning_models import AddgeneIdSource
from pydna.snapgene_history_parser import (
    parse_snapgene_history,
    SnapgeneHistoryParserWarning,
)

TEST_FOLDER = os.path.join(os.path.dirname(__file__), "snapgene_history_files")

TEST_FILES = sorted(glob.glob(os.path.join(TEST_FOLDER, "*.dna")))

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


class TestSnapgeneHistoryParser:

    def test_files_exist(self):
        assert len(TEST_FILES) == 50

    @pytest.mark.parametrize("file", TEST_FILES, ids=os.path.basename)
    def test_correctly_parsed(self, file):
        filename = os.path.basename(file)
        try:
            with warnings.catch_warnings(record=True) as wlist:
                warnings.simplefilter("always")
                seqr = parse_snapgene_history(file)
                file_warnings = [
                    str(w.message)
                    for w in wlist
                    if issubclass(w.category, SnapgeneHistoryParserWarning)
                ]
        except NotImplementedError as e:
            if filename not in METHOD_NOT_SUPPORTED:
                raise AssertionError(f"File {filename} not supported") from e
            return
        except ValueError as e:
            if (
                filename not in EXPECTED_VALUE_ERROR
                or "No product found for expected SEGUID" not in str(e)
            ):
                raise AssertionError(f"File {filename} not supported") from e
            return

        # Check special cases
        if filename == "import_addgene_then_clone.dna":
            assert isinstance(seqr.source.input[0].sequence.source, AddgeneIdSource)
        # Can't do this anymore, (see _source_from_metadata)
        # if filename == "import_ncbi.dna":
        #     assert isinstance(seqr.source, NCBISequenceSource)
        if filename == "import_addgene.dna":
            assert isinstance(seqr.source, AddgeneIdSource)

        # Check warnings
        if filename == "circularize_then_linearize_without_enzyme.dna":
            expected_warnings = ["Stopped at linearize operation without enzymes"]
        elif filename == "circularize.dna":
            expected_warnings = ["Stopped at change topology operation"]
        elif filename.startswith("manual_"):
            expected_warnings = ["Manual editing of sequences not supported"]
        else:
            expected_warnings = []

        assert file_warnings == expected_warnings

    def test_parse_snapgene_history_from_bytes(self):
        example_file = os.path.join(TEST_FOLDER, "circularize.dna")
        with open(example_file, "rb") as f:
            bytes_data = f.read()
        seqr = parse_snapgene_history(bytes_data, file_name="circularize.dna")
        assert seqr.name == "circularize"

        # Can overwrite file_name, even for files
        seqr = parse_snapgene_history(example_file, file_name="overwrite.dna")
        assert seqr.name == "overwrite"

        # Otherwise takes from file name
        seqr = parse_snapgene_history(example_file)
        assert seqr.name == "circularize"
