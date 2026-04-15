# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord
from pydna.assembly2 import SingleFragmentAssembly
from pydna.gateway import gateway_overlap, annotate_gateway_sites
from Bio.SeqFeature import SimpleLocation
from pydna.utils import rc
from pydna.opencloning_models import CloningStrategy

SimpleLocation.__repr__ = SimpleLocation.__str__

attB1 = "ACAACTTTGTACAAAAAAGCAGAAG"
attP1 = "AAAATAATGATTTTATTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCTGAACGAGAAGCGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA"


seq1 = Dseqrecord("aaa" + attB1 + "cggt" + rc(attP1) + "ttt")
seq1.add_feature(len("aaa" + attB1), len("aaa" + attB1) + 4, strand=1)

seq1 = annotate_gateway_sites(seq1, False)


def algo(x, y, limit):
    return gateway_overlap(x, y, "BP", False)


asm = SingleFragmentAssembly([seq1], algorithm=algo)


prod = asm.assemble_inversion()[0]
annotated = annotate_gateway_sites(prod, False)
prod.features += annotated.features


cs = CloningStrategy.from_dseqrecords([prod])
with open("seq1.json", "w") as f:
    f.write(cs.model_dump_json())
