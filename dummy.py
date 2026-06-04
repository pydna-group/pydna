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

homology1_options = ["AAGTCCGTTCGTTTTACCTc", "AAGTCCGTTCGTTTTACCTa"]
homology2 = "ATTACAGCATGGGAAGAAAG"

for homology1 in homology1_options:
    genome = Dseqrecord("aaaaaa" + homology1 + "ccccc" + homology2 + "aaaaaa")
    repair_template = Dseqrecord(homology1 + homology2)

    # products = homologous_recombination_integration(genome, [repair_template], limit=15)
    # print(len(products))
    asm = Assembly([genome, repair_template], limit=15, algorithm=common_sub_strings)
    asms = [a for a in asm.get_insertion_assemblies() if a[0][0] == 1]

    dseqr = assemble(asm.fragments, asms[0], True)
    print(dseqr.seq)
