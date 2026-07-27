# -*- coding: utf-8 -*-
"""Backward-compatibility shim: ``pydna.legacy.alphabet`` re-exports ``pydna.core.alphabet``.

The sequence data model moved into :mod:`pydna.core`. Code written against the
pydna 5.x layout keeps working by switching imports to
``from pydna.legacy.alphabet import ...``.
"""

from pydna.core import alphabet as _mod


def __getattr__(name):
    return getattr(_mod, name)


def __dir__():
    return dir(_mod)
