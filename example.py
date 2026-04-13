# -*- coding: utf-8 -*-
from pydna.opencloning_models import CloningStrategy
from pydna.assembly2 import Assembly, restriction_ligation_overlap
from Bio.Restriction import BbsI
from Bio.SeqFeature import SimpleLocation

SimpleLocation.__repr__ = SimpleLocation.__str__


def print_edges(asm):
    for d in asm.G.edges(data=True):
        print(d[2]["uid"])


with open("example.json", "r") as f:
    cs = CloningStrategy.model_validate_json(f.read())

product, *_ = cs.to_dseqrecords()

inputs = [inp.sequence for inp in product.source.input]

# prods = restriction_ligation_assembly(inputs, [BbsI])

input_2 = inputs[2]
input_2_digested = input_2.cut([BbsI])[1]

# print(repr(input_2_digested.seq))
new_inputs = inputs[:2] + [input_2_digested]

# prods = restriction_ligation_assembly(new_inputs, [BbsI])
# print(len(prods))


def algo(x, y, l):
    return restriction_ligation_overlap(x, y, [BbsI])


asm1 = Assembly(
    inputs, algorithm=algo, use_fragment_order=False, use_all_fragments=True
)
# print(asm1.get_locations_on_fragments())
# print_edges(asm1)

asm2 = Assembly(
    new_inputs, algorithm=algo, use_fragment_order=False, use_all_fragments=True
)
# print_edges(asm2)
# print(asm1.get_locations_on_fragments()[2])
# print(asm2.get_locations_on_fragments()[2])
# for d in asm.G.edges(data=True):
#     print(d[2]['uid'])
# print(asm2.get_circular_assemblies()[0])
# print(asm2.get_locations_on_fragments())
# print()
# print(*asm2.G.edges(data=True), sep='\n')

print(*asm1.G.edges(data=True), sep="\n")
print(asm1.get_locations_on_fragments()[1])
print(asm1.get_locations_on_fragments()[-1])
asm1.get_circular_assemblies(only_adjacent_edges=True)
print()
print(*asm2.G.edges(data=True), sep="\n")
print(asm2.get_locations_on_fragments()[1])
print(asm2.get_locations_on_fragments()[-1])
asm2.get_circular_assemblies(only_adjacent_edges=True)

# print(asm.get_circular_assemblies())

# print(asm2.assemble_circular()[0].seq.seguid() == product.seq.seguid())

# print()
# print(algo(new_inputs[0], new_inputs[1], None))
# print(algo(new_inputs[1], new_inputs[2], None))
# print(algo(new_inputs[2], new_inputs[0], None))
