#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from pydna.readers import read


def test_mark_budde():
    """test mark budde"""
    a = read("pGREG505.gb")
    assert a.name == "pGREG505"
    assert a.looped().name == "pGREG505"
    assert a.id == "pGREG505"
    assert a.looped().id == "pGREG505"

    """
    >>> a = pydna.read('pGREG505.gb')
    >>> a.name
    'pGREG505'
    >>> a.looped().name
    'pGREG505'
    >>> a.annotations
    {'comment': 'XYZZYApEinfo:methylated:1', 'source': '', 'taxonomy': [], 'keywords': [''], 'accessions': ['pGREG505'], 'data_file_division': '   ', 'date': '15-DEC-2012', 'organism': '. .'}
    >>> a.id
    'pGREG505'
    >>> a.looped().id
    'pGREG505'
    >>>
    """


if __name__ == "__main__":
    pytest.cmdline.main([__file__, "-v", "-s"])
