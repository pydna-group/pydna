# -*- coding: utf-8 -*-
from pydna.assembly2 import cre_lox_excision
from pydna.cre_lox import LOXP_SEQUENCE
from pydna.dseqrecord import Dseqrecord

genome = Dseqrecord(f"cccccc{LOXP_SEQUENCE}aaaa{LOXP_SEQUENCE}cccccc", name="genome")
products = cre_lox_excision(genome)

for p in products:
    p.source.validate(p)
