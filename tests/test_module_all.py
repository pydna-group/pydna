#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_repr():
    pytest.importorskip("requests")

    from pydna import all

    assert all.__all__ == [
        "Anneal",
        "pcr",
        "Assembly",
        "genbank",
        "Genbank",
        "Dseqrecord",
        "Dseq",
        "read",
        "read_primer",
        "parse",
        "parse_primers",
        "primer_design",
        "assembly_fragments",
        "eq",
        "gbtext_clean",
    ]
