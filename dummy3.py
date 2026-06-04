# -*- coding: utf-8 -*-
# flake8: noqa
from pydna.all import Dseqrecord
from pydna.opencloning_models import CloningStrategy
from pydna.assembly2 import (
    homologous_recombination_integration,
    Assembly,
    common_sub_strings,
    assemble,
)
from pydna.readers import read
from Bio.SeqFeature import SimpleLocation

SimpleLocation.__repr__ = SimpleLocation.__str__


genome = Dseqrecord("CCGAGGGGAATC")
del_G = Dseqrecord("CCGAGGGAATC")
add_G = Dseqrecord("CCGAGGGGGAATC")

for repair_template in [del_G, add_G]:
    asm = Assembly([genome, repair_template], limit=4, algorithm=common_sub_strings)
    asms = asm.get_insertion_assemblies()
    asms = [a for a in asms if a[0][0] == 1]
    print(asms)
    dseqr = assemble(asm.fragments, asms[0], True)
    print(dseqr.seq)
