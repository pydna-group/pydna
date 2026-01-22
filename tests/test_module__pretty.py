#!/usr/bin/env python
# -*- coding: utf-8 -*-


def test_repr_pretty_():
    from pydna._pretty import pretty_str
    from unittest.mock import MagicMock

    pp = MagicMock()
    s = pretty_str("123")
    s._repr_pretty_(pp, None)
    pp.text.assert_called_with("123")


def test_PrettyTable():
    from pydna._pretty import pretty_str
    from pydna._pretty import PrettyTable

    x = PrettyTable()

    x.field_names = ["City name", "Area", "Population", "Annual Rainfall"]
    x.add_rows(
        [
            ["Adelaide", 1295, 1158259, 600.5],
            ["Brisbane", 5905, 1857594, 1146.4],
            ["Darwin", 112, 120900, 1714.7],
            ["Hobart", 1357, 205556, 619.5],
            ["Sydney", 2058, 4336374, 1214.8],
            ["Melbourne", 1566, 3806092, 646.9],
            ["Perth", 5386, 1554769, 869.4],
        ]
    )

    assert x.lol() == [
        ["City name", "Area", "Population", "Annual Rainfall"],
        ["Adelaide", 1295, 1158259, 600.5],
        ["Brisbane", 5905, 1857594, 1146.4],
        ["Darwin", 112, 120900, 1714.7],
        ["Hobart", 1357, 205556, 619.5],
        ["Sydney", 2058, 4336374, 1214.8],
        ["Melbourne", 1566, 3806092, 646.9],
        ["Perth", 5386, 1554769, 869.4],
    ]

    assert x.get_string() == x.__repr__()
    from prettytable import TableStyle
    from copy import copy

    c = copy(x)
    c.set_style(TableStyle.MARKDOWN)
    assert x._repr_markdown_() == pretty_str(c.get_string())
