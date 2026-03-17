# -*- coding: utf-8 -*-
from pydna.cre_lox import LOXP_SEQUENCE
from pydna.dseqrecord import Dseqrecord
from pydna.assembly2 import cre_lox_integration, cre_lox_excision
from pydna.opencloning_models import CloningStrategy

a = Dseqrecord(f"cccccc{LOXP_SEQUENCE}aaaaa")
b = Dseqrecord(f"{LOXP_SEQUENCE}bbbbb", circular=True)
integration_product, *_ = cre_lox_integration(a, [b])
excision_product, *_ = cre_lox_excision(integration_product)

a.name = "genome"
b.name = "plasmid"
integration_product.name = "integration_product"
excision_product.name = "excision_product"

cs = CloningStrategy.from_dseqrecords([excision_product])
with open("cre_lox.json", "w") as f:
    f.write(cs.model_dump_json())


print(cs.to_dseqrecords())

print(cs.to_dseqrecords()[0].history())
