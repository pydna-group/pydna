#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""Provides two functions, read and read_primer."""
from pydna.parsers import parse
from pydna.primer import Primer


def read(data, ds=True):
    """This function is similar the :func:`parse` function but expects one and only
    one sequence or and exception is thrown.

    Parameters
    ----------
    data : string
        see below
    ds : bool
        Double stranded or single stranded DNA, if True return
        Dseqrecord objects, else Bio.SeqRecord objects.

    Returns
    -------
    Dseqrecord
        contains the first Dseqrecord or SeqRecord object parsed.

    Notes
    -----

    The data parameter is similar to the data parameter for :func:`parse`.

    See Also
    --------
    parse

    """

    try:
        (result,) = parse(data, ds)
    except ValueError as err:
        msg = str(err)

        if "too many" in msg:
            raise ValueError(
                f"More than one sequence found in data ({str(data)[:79]})"
            ) from err
        elif "not enough" in msg:
            raise ValueError(f"No sequence found in data ({str(data)[:79]})") from err
        else:  # pragma: no cover
            raise err  # re-raises the same ValueError with original traceback
    return result


def read_primer(data):
    """Use this function to read a primer sequence from a string or a local file.
    The usage is similar to the :func:`parse_primer` function."""

    return Primer(read(data, ds=False))
