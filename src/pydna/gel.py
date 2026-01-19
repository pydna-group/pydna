#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# doctest: +NORMALIZE_WHITESPACE
# doctest: +SKIP
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""docstring."""

import math
from pydna.ladders import GeneRuler_1kb_plus as _mwstd


def interpolator(mwstd):
    """docstring."""
    from scipy.interpolate import CubicSpline

    interpolator = CubicSpline(
        [len(fr) for fr in mwstd[::-1]],
        [fr.rf for fr in mwstd[::-1]],
        bc_type="natural",
        extrapolate=False,
    )
    interpolator.mwstd = mwstd
    return interpolator


def gel(
    samples=None, gel_length=600, margin=50, interpolator=interpolator(mwstd=_mwstd)
):
    import numpy as np
    from PIL import Image
    from PIL import ImageDraw

    """docstring."""
    max_intensity = 256
    lane_width = 50
    lanesep = 10
    start = 10
    samples = samples or [interpolator.mwstd]
    width = int(60 + (lane_width + lanesep) * len(samples))
    lanes = np.zeros((len(samples), gel_length), dtype=int)
    image = Image.new("RGB", (width, gel_length), "#ddd")
    draw = ImageDraw.Draw(image)
    draw.rectangle((0, 0, (width, gel_length)), fill=(0, 0, 0))
    scale = (gel_length - margin) / interpolator(min(interpolator.x))

    for labelsource in samples[0]:
        peak_centre = (interpolator(len(labelsource))) * scale - 5 + start
        label = f"{len(labelsource):<5} -"
        draw.text((2, peak_centre), label, fill=(255, 255, 255))

    for lane_number, lane in enumerate(samples):
        for band in lane:
            log = math.log(len(band), 10)
            height = (band.m() / (240 * log)) * 1e10
            peak_centre = interpolator(len(band)) * scale + start
            max_spread = 10
            # if len(band) <= 50:
            #     peak_centre += 50
            #     max_spread *= 4
            #     max_intensity /= 10
            band_spread = max_spread / log
            for i in range(max_spread, 0, -1):
                y1 = peak_centre - i
                y2 = peak_centre + i
                intensity = (
                    height
                    * math.exp(
                        -float(((y1 - peak_centre) ** 2)) / (2 * (band_spread**2))
                    )
                    * max_intensity
                )
                for y in range(int(y1), int(y2)):
                    try:
                        lanes[lane_number][y] += intensity
                    except IndexError:
                        pass

    for i, lane in enumerate(lanes):
        max_intensity = np.amax(lanes[i])
        if max_intensity > 256:
            lanes[i] = np.multiply(lanes[i], 256)
            lanes[i] = np.divide(lanes[i], max_intensity)

    for i, lane in enumerate(lanes):
        x1 = 50 + i * (lane_width + lanesep)
        x2 = x1 + lane_width
        for y, intensity in enumerate(lane):
            y1 = y
            y2 = y + 1
            draw.rectangle((x1, y1, x2, y2), fill=(intensity, intensity, intensity))

    return image


# Inverting and rotating the gel
# im = gel([ GeneRuler_1kb_plus, [band, ]])
# from PIL import ImageOps
# im_invert = ImageOps.invert(im)
# im.rotate(90, expand=1)
