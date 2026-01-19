#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_genbankfile():
    from pydna import genbankfile

    gbf = genbankfile.GenbankFile("aaa", path="/path/")

    assert gbf.path == "/path/"

    assert repr(gbf) == "File(id)(-3)"

    assert gbf._repr_html_() == "<a href='/path/' target='_blank'>/path/</a><br>"

    from unittest.mock import MagicMock

    pp = MagicMock()
    gbf._repr_pretty_(pp, None)
    pp.text.assert_called_with(repr(gbf))

    gbf = gbf.rc()


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
