# -*- coding: utf-8 -*-
"""Backward-compatibility shim: ``pydna.legacy.dseqrecord`` re-exports ``pydna.core.dseqrecord``.

The sequence data model moved into :mod:`pydna.core`. Code written against the
pydna 5.x layout keeps working by switching imports to
``from pydna.legacy.dseqrecord import ...``.
"""

from pydna.core import dseqrecord as _mod


def __getattr__(name):
    return getattr(_mod, name)


def __dir__():
    return dir(_mod)
