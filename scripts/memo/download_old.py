#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Provides a class for downloading sequences from genbank."""
import pickle as _pickle
import shelve as _shelve
import os as _os
import requests as _requests
import textwrap as _textwrap

import logging as _logging

_module_logger = _logging.getLogger("pydna." + __name__)


def download_text(url):
    cached = False
    refresh = False
    cache = _shelve.open(
        _os.path.join(_os.environ["pydna_data_dir"], "web"),
        protocol=_pickle.HIGHEST_PROTOCOL,
        writeback=False,
    )
    key = str(url)

    if _os.environ["pydna_cache"] in ("compare", "cached"):
        try:
            cached = cache[key]
        except KeyError:
            if _os.environ["pydna_cache"] == "compare":
                raise Exception("no result for this key!")
            else:
                refresh = True

    if refresh or _os.environ["pydna_cache"] in ("compare", "refresh", "nocache"):
        req = _requests.get(url)
        result = req.text
    if _os.environ["pydna_cache"] == "compare":
        if result != cached:
            _module_logger.warning("download error")

    if refresh or _os.environ["pydna_cache"] == "refresh":
        cache = _shelve.open(
            _os.path.join(_os.environ["pydna_data_dir"], "genbank"),
            protocol=_pickle.HIGHEST_PROTOCOL,
            writeback=False,
        )
        cache[key] = result

    elif cached and _os.environ["pydna_cache"] not in ("nocache", "refresh"):
        result = cached

    cache.close()

    result = _textwrap.dedent(result).strip()
    result = result.replace("\r\n", "\n")
    result = result.replace("\r", "\n")
    return result


if __name__ == "__main__":
    cache = os.getenv("pydna_cache")
    os.environ["pydna_cache"] = "nocache"
    import doctest

    doctest.testmod()
    os.environ["pydna_cache"] = cache
