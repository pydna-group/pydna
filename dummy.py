# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord
from pydna.assembly2 import recombinase_excision
from pydna.recombinase import Recombinase
from pydna.assembly2 import SingleFragmentAssembly


site1 = "ATGCCCTAAaaCT"
site2 = "AAaaTTTTTTTCCCT"

rec = Recombinase(site1, site2)
genome = Dseqrecord(f"cccccc{site1.upper()}tttt{site2.upper()}aaaaa")

products = recombinase_excision(genome, rec)

asm = SingleFragmentAssembly([genome], None, rec.overlap)
print(*asm.G.edges(data=True), sep="\n")
print(asm.get_circular_assemblies())
print()

seq2 = Dseqrecord(f"ag{site1.upper()}gtca{site1.upper()}aa")

asm = SingleFragmentAssembly([seq2], limit=13)
print(*asm.G.edges(data=True), sep="\n")
# print(asm.get_circular_assemblies())
