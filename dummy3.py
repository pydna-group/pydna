# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord
from Bio.Restriction import EcoRI, SalI
from pydna.opencloning_models import CloningStrategy

a = Dseqrecord("aaGAATTCccGTCGACaa")

c, d, e = a.cut(EcoRI, SalI)

cs = CloningStrategy.from_dseqrecord(c)

cs.add_dseqrecord(d)
cs.add_dseqrecord(e)

with open("dummy3.json", "w") as f:
    f.write(cs.model_dump_json())
