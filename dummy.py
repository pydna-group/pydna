# -*- coding: utf-8 -*-
from pydna.assembly2 import homologous_recombination_integration
from pydna.dseqrecord import Dseqrecord

homology = "AAGTCCGTTCGTTTTACCTG"
genome = Dseqrecord(f"aaaaaa{homology}cccc", name="genome")
insert_seq = Dseqrecord(f"{homology}tttt{homology}", name="insert")
products = homologous_recombination_integration(genome, [insert_seq], 20)
for p in products:
    p.source.validate(p)
