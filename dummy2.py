# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord
from pydna.assembly2 import homologous_recombination_integration

a = Dseqrecord("1CGTACGCACAxxxxC2")
b = Dseqrecord("3CGTACGCACAyyyyCGTACGCACAT4")
res = homologous_recombination_integration(a, b, minimal_homology=10)

results = ["1CGTACGCACAyyyyCGTACGCACAxxxxC2"]

print(res[0].seq)
