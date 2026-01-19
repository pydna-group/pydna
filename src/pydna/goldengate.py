#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

#                                          ^^
#     ^^      ..                                       ..
#             []                                       []
#           .:[]:_          ^^                       ,:[]:.
#         .: :[]: :-.                             ,-: :[]: :.
#       .: : :[]: : :`._                       ,.': : :[]: : :.
#     .: : : :[]: : : : :-._               _,-: : : : :[]: : : :.
# _..: : : : :[]: : : : : : :-._________.-: : : : : : :[]: : : : :-._
# _:_:_:_:_:_:[]:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:_:[]:_:_:_:_:_:_
# !!!!!!!!!!!![]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!![]!!!!!!!!!!!!!
# ^^^^^^^^^^^^[]^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^[]^^^^^^^^^^^^^
#             []                                       []
#             []                                       []
#             []                                       []
#  ~~^-~^_~^~/  \~^-~^~_~^-~_^~-^~_^~~-^~_~^~-~_~-^~_^/  \~^-~_~^-~~-
# ~ _~~- ~^-^~-^~~- ^~_^-^~~_ -~^_ -~_-~~^- _~~_~-^_ ~^-^~~-_^-~ ~^

"""Assembly of sequences by GoldenGate ligation assembly."""
from Bio.Restriction import BsaI, BsmBI, BbsI, FokI
from pydna.dseqrecord import Dseqrecord


BsaI, BsmBI, BbsI, FokI

DNA = Dseqrecord("gatcGAAGACtagagtctgattcg")

a, b = DNA.cut(BbsI)

assert a + b == DNA


# MoClo

# https://edinburgh-genome-foundry.github.io/GoldenHinges
