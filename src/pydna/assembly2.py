# -*- coding: utf-8 -*-
"""Deprecated alias for :mod:`pydna.assembly`.

``pydna.assembly2`` was renamed to ``pydna.assembly`` when the rewritten
implementation became the canonical one. This shim forwards every attribute to
``pydna.assembly`` so existing ``from pydna.assembly2 import ...`` code keeps
working. It will be removed in a future version.
"""

import warnings

from pydna import _PydnaDeprecationWarning
from pydna import assembly as _assembly

warnings.warn(
    "pydna.assembly2 was renamed to pydna.assembly; update your imports to "
    "'from pydna.assembly import ...'. This alias will be removed in a future version.",
    _PydnaDeprecationWarning,
    stacklevel=2,
)


def __getattr__(name):
    return getattr(_assembly, name)


def __dir__():
    return dir(_assembly)
