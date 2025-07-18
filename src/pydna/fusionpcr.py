#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""docstring."""
from pydna.common_sub_strings import terminal_overlap
from pydna.assembly import Assembly


def fuse_by_pcr(fragments, limit=15):
    """docstring."""
    asm = Assembly(fragments, limit, algorithm=terminal_overlap)
    results = asm.assemble_linear()
    return results


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
