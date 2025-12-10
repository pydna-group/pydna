#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from pydna.dseq import Dseq

def dseq():

    x = Dseq( "gGGATCC",
             "  CCTAGG"[::-1], 1)

    y = Dseq(" gGGATCC",
             "  CCTAGG"[::-1], 0)

    z = Dseq("  gGGATCC",
              "  CCTAGG"[::-1], -1)

    assert x == y == z

    x = Dseq( " GGATCC",
             " gCCTAGG"[::-1], 1)

    y = Dseq("  GGATCC",
             " gCCTAGG"[::-1], 0)

    z = Dseq("   GGATCC",
              " gCCTAGG"[::-1], -1)

    assert x == y == z
