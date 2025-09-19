# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord
from pydna.assembly2 import golden_gate_assembly
from Bio.Restriction import BsaI
from pydna.opencloning_models import CloningStrategy


insert1 = Dseqrecord("GGTCTCAattaAAAAAttaaAGAGACC")
insert2 = Dseqrecord("GGTCTCAttaaCCCCCatatAGAGACC")
insert3 = Dseqrecord("GGTCTCAatatGGGGGccggAGAGACC")


vector = Dseqrecord("TTTTattaAGAGACCTTTTTGGTCTCAccggTTTT", circular=True)


assembly_output = golden_gate_assembly(
    [insert1, insert2, insert3, vector], [BsaI], circular_only=True
)

# print(assembly_output[0].source)

cs = CloningStrategy.from_dseqrecord(assembly_output[0])


with open("dummy.json", "w") as f:
    f.write(cs.model_dump_json())

CloningStrategy.model_validate(cs.model_dump())
