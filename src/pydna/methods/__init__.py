# -*- coding: utf-8 -*-
"""pydna.methods — the cloning techniques pydna can simulate.

Each technique lives in its own module: the sites, enzymes or nucleases it
recognises, and the :class:`~pydna.methods._engine.Method` that declares how the
reaction runs.

A ``Method`` is a description, not code — an overlap algorithm, an assembly
shape, the provenance it records, and optional hooks. One engine executes them
all, so every technique records its history the same way and adding a technique
does not touch the shared assembly graph.
"""

from pydna.methods._engine import Method, Shape, run
from pydna.methods._registry import (
    register,
    get,
    available,
    methods,
)

__all__ = [
    "Method",
    "Shape",
    "run",
    "register",
    "get",
    "available",
    "methods",
]
