# -*- coding: utf-8 -*-
import os
from pydna.assembly2 import gateway_assembly
from pydna.parsers import parse_snapgene
from pydna.opencloning_models import CloningStrategy

test_files = "tests"

line = "pDONRtm221.dna	pcr_product-attP1_1-attP2_1.dna	entry-attP1_1-attP2_1.dna	pET-53-DESTtm.dna	expression-attP1_1-attP2_1.dna".split(
    "\t"
)
backbone, pcr_product, _, backbone_expression, _ = line
backbone = parse_snapgene(os.path.join("tests/gateway_manual_cloning/" + backbone))[0]
pcr_product = parse_snapgene(
    os.path.join("tests/gateway_manual_cloning/" + pcr_product)
)[0]
backbone_expression = parse_snapgene(
    os.path.join("tests/gateway_manual_cloning/" + backbone_expression)
)[0]

# Works with the right reaction
entry_vector = gateway_assembly(
    [backbone, pcr_product], "BP", multi_site_only=True, circular_only=True
)[0]

expression_clone, _ = gateway_assembly(
    [backbone_expression, entry_vector], "LR", multi_site_only=True, circular_only=True
)

cs = CloningStrategy.from_dseqrecord(expression_clone)
with open("dummy4.json", "w") as f:
    f.write(cs.model_dump_json())
