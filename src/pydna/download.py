#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Provides a function for downloading online text files."""

import textwrap

from pydna._pretty import pretty_str as ps


def download_text(url):
    """docstring."""
    import requests

    req = requests.get(url)

    result = textwrap.dedent(req.text).strip()
    result = result.replace("\r\n", "\n").replace("\r", "\n")

    return ps(result)
