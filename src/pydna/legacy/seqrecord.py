# -*- coding: utf-8 -*-
"""Backward-compatibility shim: ``pydna.legacy.seqrecord`` re-exports ``pydna.core.seqrecord``.

The sequence data model moved into :mod:`pydna.core`. Code written against the
pydna 5.x layout keeps working by switching imports to
``from pydna.legacy.seqrecord import ...``.
"""

from pydna.core import seqrecord as _mod


def __getattr__(name):
    return getattr(_mod, name)


def __dir__():
    return dir(_mod)
