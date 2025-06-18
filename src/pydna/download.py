#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Provides a function for downloading online text files."""

import textwrap as _textwrap

# import os as _os
from pydna._pretty import pretty_str as _pretty_str

# from pydna.utils import memorize as _memorize
# import logging as _logging

# _module_logger = _logging.getLogger("pydna." + __name__)


# @_memorize("pydna.download.download_text")
def download_text(url):
    """docstring."""
    import requests

    #    _module_logger.info("#### DOWNLOAD TEXT ####")
    #    _module_logger.info("url = %s", url)
    req = requests.get(url)
    #    _module_logger.info("url = %s", str(req))
    result = _textwrap.dedent(req.text).strip()
    result = result.replace("\r\n", "\n").replace("\r", "\n")
    #    _module_logger.info("result[:160] = %s", result[:160])
    return _pretty_str(result)
