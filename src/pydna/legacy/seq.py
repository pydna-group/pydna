# -*- coding: utf-8 -*-
"""Backward-compatibility shim: ``pydna.legacy.seq`` re-exports ``pydna.core.seq``.

The sequence data model moved into :mod:`pydna.core`. Code written against the
pydna 5.x layout keeps working by switching imports to
``from pydna.legacy.seq import ...``.
"""

from pydna.core import seq as _mod


def __getattr__(name):
    return getattr(_mod, name)


def __dir__():
    return dir(_mod)
