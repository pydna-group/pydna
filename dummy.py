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
