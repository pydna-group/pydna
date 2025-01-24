# -*- coding: utf-8 -*-
from pydna.design import primer_design
from pydna.dseqrecord import Dseqrecord

t = Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
from pydna.primer import Primer

pf = Primer("atgactgctaacccttccttggtg")
print(pf.name)
ampl = primer_design(t, fp=pf)
print(ampl.forward_primer.name)
