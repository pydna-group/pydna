# -*- coding: utf-8 -*-
from pydna.dseq import Dseq

seq = Dseq("AGEEGaGJJJg", circular=True)

cutsite_pairs = seq.get_cutsite_pairs(seq.get_ds_meltsites(3))

shifted_cutsite_pairs = seq.shift_melt_cutsite_pairs(cutsite_pairs)

assert shifted_cutsite_pairs == [(((10, 6), None), ((7, 5), None))]

expected_product = seq.apply_cut(((10, 6), None), ((7, 5), None), allow_overlap=True)


print(expected_product._data)
for shift in range(len(seq)):
    new_seq = seq.shifted(shift)

    cutsite_pairs = new_seq.get_cutsite_pairs(new_seq.get_ds_meltsites(3))
    shifted_cutsite_pairs = new_seq.shift_melt_cutsite_pairs(cutsite_pairs)
    assert len(shifted_cutsite_pairs) == 1
    prod = new_seq.apply_cut(*shifted_cutsite_pairs[0], allow_overlap=True)
    print(prod._data)
    assert prod == expected_product
