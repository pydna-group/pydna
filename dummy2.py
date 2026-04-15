# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord
from pydna.assembly2 import SingleFragmentAssembly
from Bio.SeqFeature import SimpleLocation
from pydna.utils import rc
from pydna.opencloning_models import CloningStrategy

SimpleLocation.__repr__ = SimpleLocation.__str__

hom = "ACAACTTTGTACAAAAAAGCAGAAG"


seq1 = Dseqrecord("ggg" + hom + "cca" + rc(hom) + "tttt")
seq1.add_feature(3, len("ggg" + hom), strand=1)
seq1.add_feature(
    len("ggg" + hom + "cca"), len("ggg" + hom + "cca" + rc(hom)), strand=-1
)
seq1.add_feature(len("ggg" + hom), len("ggg" + hom) + 3, strand=1)

asm = SingleFragmentAssembly([seq1], limit=20)


prod = asm.assemble_inversion()[0]

cs = CloningStrategy.from_dseqrecords([prod])
with open("seq2.json", "w") as f:
    f.write(cs.model_dump_json())
