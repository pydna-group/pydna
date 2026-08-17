# -*- coding: utf-8 -*-
"""Backward-compatibility shim: ``pydna.legacy.dseq`` re-exports ``pydna.core.dseq``.

The sequence data model moved into :mod:`pydna.core`. Code written against the
pydna 5.x layout keeps working by switching imports to
``from pydna.legacy.dseq import ...``.
"""

from pydna.core import dseq as _mod


def __getattr__(name):
    return getattr(_mod, name)


def __dir__():
    return dir(_mod)
