#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Agarose gel DNA ladders.

A DNA ladder is a list of FakeSeq objects that has to be initiated with
Size (bp), amount of substance (mol) and Relative mobility (Rf).

Rf is a float value between 0.000 and 1.000. These are used together with
the cubic spline interpolator in the gel module to calculate migartion
distance from fragment length. The Rf values are calculated manually from
a gel image. Exampel can be found in scripts/molecular_weight_standards.ods.
"""


from pydna.fakeseq import FakeSeq


PennStateLadder = [
    FakeSeq(int(n)) for n in (10000, 7750, 5000, 4000, 3000, 2000, 1500, 1000, 750, 500)
]


GeneRuler_1kb = [
    FakeSeq(int(n))
    for n in (
        10000,
        8000,
        6000,
        5000,
        4000,
        3500,
        3000,
        2500,
        2000,
        1500,
        1000,
        750,
        500,
        250,
    )
]

# Google spreadsheet to make the ladders below
# https://docs.google.com/spreadsheets/d/1vN0y75ibxPrG6yJQjq1uF2FXP0L-qGSn_fzInUHeTs4/edit#gid=0

GeneRuler_1kb_plus = [
    FakeSeq(ln, n=n * 1e-15, rf=rf)
    for ln, n, rf in (
        # (length, fmol, Rf )
        (20000, 1.538, 0.000),
        (10000, 3.077, 0.040),
        (7000, 4.396, 0.096),
        (5000, 23.077, 0.154),
        (4000, 7.692, 0.201),
        (3000, 10.256, 0.261),
        (2000, 15.385, 0.362),
        (1500, 82.051, 0.443),
        (1000, 38.462, 0.562),
        (700, 54.945, 0.667),
        (500, 230.769, 0.755),
        (400, 96.154, 0.808),
        (300, 128.205, 0.866),
        (200, 192.308, 0.931),
        (75, 512.821, 1.000),
    )
]


HI_LO_DNA_MARKER = [
    FakeSeq(ln, n=n * 1e-15, rf=rf)
    for ln, n, rf in (
        # (length, fmol, Rf )
        (10000, 4.545, 0.000),
        (8000, 5.682, 0.013),
        (6000, 11.364, 0.037),
        (4000, 22.727, 0.081),
        (3000, 42.929, 0.131),
        (2000, 113.636, 0.236),
        (1550, 97.752, 0.302),
        (1400, 108.225, 0.331),
        (1000, 181.818, 0.438),
        (750, 60.606, 0.522),
        (500, 181.818, 0.638),
        (400, 75.758, 0.703),
        (300, 202.02, 0.772),
        (200, 227.273, 0.856),
        (100, 303.03, 0.955),
        (50, 454.545, 1.000),
    )
]


# GeneRuler_Mix = [_FakeSeq(int(ln),) for ln,n in (
#                  10000,
#                  8000,
#                  6000,
#                  5000,
#                  4000,
#                  3500,
#                  3000,
#                  2500,
#                  2000,
#                  1500,
#                  1200,
#                  1000,
#                  900,
#                  800,
#                  700,
#                  600,
#                  500,
#                  400,
#                  300,
#                  200,
#                  100,)]


FakeGel = [
    [
        FakeSeq(1000),
        FakeSeq(2000),
    ],
    [
        FakeSeq(3000),
        FakeSeq(4000),
    ],
    [
        FakeSeq(5000),
        FakeSeq(6000),
    ],
    PennStateLadder,
]
