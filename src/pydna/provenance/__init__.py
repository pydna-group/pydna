# -*- coding: utf-8 -*-
"""pydna.provenance — pydna's own record of how a sequence was made.

A :class:`Step` records the method, its inputs and its parameters, which is all
:func:`replay` needs to reproduce the result. Interoperability with OpenCloning
lives in :mod:`pydna.provenance.adapters.opencloning` as a translation of that
record, so pydna can describe a technique before — or without — OpenCloning
having a matching source type.
"""

from pydna.provenance.step import Step, replay

__all__ = ["Step", "replay"]
