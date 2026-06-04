# -*- coding: utf-8 -*-
# flake8: noqa
from pydna.dseqrecord import Dseqrecord
from pydna.assembly2 import crispr_integration, homologous_recombination_integration
from Bio.SeqFeature import SimpleLocation, CompoundLocation

SimpleLocation.__repr__ = SimpleLocation.__str__
CompoundLocation.__repr__ = CompoundLocation.__str__

homology1 = "AAGTCACCTG"
homology2 = "ctgaaaaaaa"
genome = Dseqrecord(f"aa{homology1}cc{homology2}")
insert = Dseqrecord(f"{homology1}gggg{homology2}", circular=True).shifted(12)

# print(crispr_integration(genome, [insert], [guide], 20)[0])
print(homologous_recombination_integration(genome, [insert], 10)[0])
