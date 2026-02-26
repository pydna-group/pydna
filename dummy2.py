# -*- coding: utf-8 -*-
from pydna.dseq import Dseq

seq = Dseq("AGEEGaGJJJg", circular=True)

cutsite_pairs = seq.get_cutsite_pairs(seq.get_ds_meltsites(3))

shifted_cutsite_pairs = seq.shift_melt_cutsite_pairs(cutsite_pairs)

print(cutsite_pairs)
print(shifted_cutsite_pairs)
assert shifted_cutsite_pairs == [((10, 6), None), ((7, 5), None)]

expected_product = seq.apply_cut(((10, 6), None), ((7, 5), None), allow_overlap=True)
