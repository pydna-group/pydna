# -*- coding: utf-8 -*-
from pydna.dseq import Dseq

seq = Dseq("AGJJJGaGEEg")

cutsite_pairs = seq.get_cutsite_pairs(seq.get_ds_meltsites(2))
new_cutsite_pairs = seq.shift_melt_cutsite_pairs(cutsite_pairs)

assert new_cutsite_pairs[0] == (None, ((2, 2), None))
assert new_cutsite_pairs[1] == (((5, 5), None), ((11, 3), None))
assert new_cutsite_pairs[2] == (((11, 1), None), None)

seq = Dseq("AGEEGaGJJJg")


cutsite_pairs = seq.get_cutsite_pairs(seq.get_ds_meltsites(2))
new_cutsite_pairs = seq.shift_melt_cutsite_pairs(cutsite_pairs)

assert new_cutsite_pairs[0] == (None, ((0, -2), None))
assert new_cutsite_pairs[1] == (((0, -4), None), ((7, -4), None))
assert new_cutsite_pairs[2] == (((10, -1), None), None)

seq = Dseq("AGEEGaGJJJg", circular=True)


for shift in range(len(seq)):
    new_seq = seq.shifted(shift)
    cutsite_pairs = new_seq.get_ds_meltsites(2)
    assert len(cutsite_pairs) == 0

    print("==")
    cutsite_pairs = new_seq.get_cutsite_pairs(new_seq.get_ds_meltsites(3))

    print(repr(new_seq))
    expected_product = new_seq.apply_cut(((10, 6), None), ((4, 4), None))
    print(repr(expected_product))
    break

    print(cutsite_pairs)
    print(new_seq.shift_melt_cutsite_pairs(cutsite_pairs))
    print(repr(new_seq))
    print(new_seq.melt(2))
