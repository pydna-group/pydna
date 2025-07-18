#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from IPython.lib.latextools import latex_to_png


class Circle(object):
    def __init__(self, radius):
        self.radius = radius

    def _repr_pretty_(self, p, cycle):
        p.text("\u25cb")

    def _repr_html_(self):
        return "<h1>Cirle: radius=%s</h1>" % self.radius

    def _repr_svg_(self):
        return """<svg>
<circle cx="100" cy="50" r="40" stroke="black" stroke-width="2" fill="red"/>
</svg>"""

    def _repr_png_(self):
        return latex_to_png("$\circle$")


class MyCircle(object):
    def _repr_html_(self):
        return "&#x25CB; (<b>html</b>)"

    def _repr_svg_(self):
        return """<svg width="100px" height="100px">
           <circle cx="50" cy="50" r="20" stroke="black" stroke-width="1" fill="blue"/>
        </svg>"""

    def _repr_latex_(self):
        return r"$\bigcirc \LaTeX$"

    def _repr_javascript_(self):
        return "alert('I am a circle!');"


if __name__ == "__main__":
    import os

    cache = os.getenv("pydna_cache")
    os.environ["pydna_cache"] = "nocache"
    import doctest

    doctest.testmod()
    os.environ["pydna_cache"] = cache
