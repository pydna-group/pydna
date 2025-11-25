# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord
from Bio.Restriction import EcoRI

c, d = Dseqrecord("aaaaaaGAATTCtttttttt").cut(EcoRI)

print(c.source.model_dump_json())
