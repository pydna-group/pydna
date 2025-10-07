# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord
from Bio.Restriction import EcoRI, SalI
from pydna.opencloning_models import CloningStrategy
from pydna.assembly2 import ligation_assembly

a = Dseqrecord("aaGAATTCccGTCGACaa")

c, d, e = a.cut(EcoRI, SalI)

f = ligation_assembly([c, d, e])

cs = CloningStrategy.from_dseqrecord(f[0])

with open("dummy3.json", "w") as f:
    f.write(cs.model_dump_json())
