# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord
from pydna.assembly2 import gibson_assembly
from pydna.opencloning_models import CloningStrategy

fragments = [
    Dseqrecord("TTTTacgatAAtgctccCCCC", circular=False),
    Dseqrecord("CCCCtcatGGGG", circular=False),
    Dseqrecord("GGGGatataTTTT", circular=False),
]

products = gibson_assembly(fragments, limit=4)


cs = CloningStrategy.from_dseqrecord(products[0])
with open("dummy2.json", "w") as f:
    f.write(cs.model_dump_json())
