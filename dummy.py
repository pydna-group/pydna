# -*- coding: utf-8 -*-
from pydna.dseqrecord import Dseqrecord
from pydna.assembly import Assembly
from pydna.assembly2 import Assembly as Assembly2

a = Dseqrecord("tcgatgctatactgtgCCNCCtgtgctgtgctcta")
a.add_feature(0, 10, label="a_feat")
a_feat_seq = a.features[0].extract(a)
# 12345678901234
b = Dseqrecord("tgtgctgtgctctaTTTTTTTtattctggctgtatcCCCCCC")
b.add_feature(0, 10, label="b_feat")
b_feat_seq = b.features[0].extract(b)

# 123456789012345
c = Dseqrecord(
    "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGtattctggctgtatcGGGGGtacgatgctatactgtg"
)
c.add_feature(0, 10, label="c_feat")
c_feat_seq = c.features[0].extract(c)

feature_sequences = {
    "a_feat": a_feat_seq,
    "b_feat": b_feat_seq,
    "c_feat": c_feat_seq,
}

a.name = "aaa"  # 1234567890123456
b.name = "bbb"
c.name = "ccc"
asm = Assembly((a, b, c), limit=14)
asm2 = Assembly2((a, b, c), limit=14)
x = asm.assemble_linear()[1]
x2 = asm2.assemble_linear()[0]
# print(x.features)
# print(x)
answer = "aaa|14\n    \\/\n    /\\\n    14|bbb|15\n           \\/\n           /\\\n           15|ccc"

# print(x.figure())
# print(x2.figure())

x = asm.assemble_circular()[0]
x2 = asm2.assemble_circular()[0]

print(x.figure())
print(x2.figure())
