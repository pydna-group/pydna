# -*- coding: utf-8 -*-
import pytest


def test_gel():
    pytest.importorskip("PIL")

    import numpy as np

    from PIL import Image

    from pydna.gel import gel

    from pydna.gel import interpolator

    from pydna.ladders import PennStateLadder

    from pydna.ladders import HI_LO_DNA_MARKER

    from pydna.dseqrecord import Dseqrecord

    mygel = gel([PennStateLadder, [Dseqrecord("A" * 2000)]])

    im = Image.open("pydna_gel.png")

    frame1 = np.asarray(mygel)
    frame2 = np.asarray(im)

    assert np.array_equal(frame1, frame2)

    ip = interpolator(HI_LO_DNA_MARKER)

    x = Dseqrecord("A" * 50)

    x.n *= 10

    mygel2 = gel([[x]], interpolator=ip)

    im = Image.open("pydna_gel2.png")

    frame3 = np.asarray(mygel2)
    frame4 = np.asarray(im)

    assert np.array_equal(frame3, frame4)
