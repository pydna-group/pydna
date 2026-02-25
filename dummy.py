# -*- coding: utf-8 -*-
from pydna.dseq import Dseq

seq = Dseq("AGJGaGEg")

cutsite_pairs = seq.get_cutsite_pairs(seq.get_ds_meltsites(2))
new_cutsite_pairs = seq.shift_melt_cutsite_pairs(cutsite_pairs)

assert new_cutsite_pairs[0] == (None, ((2, 2), None))
assert new_cutsite_pairs[1] == (((3, 3), None), ((8, 2), None))
assert new_cutsite_pairs[2] == (((8, 1), None), None)


seq = Dseq("AGEGaGJg")


cutsite_pairs = seq.get_cutsite_pairs(seq.get_ds_meltsites(2))
new_cutsite_pairs = seq.shift_melt_cutsite_pairs(cutsite_pairs)

seq2 = Dseq("AGGGaGGg")


assert new_cutsite_pairs[0] == (None, ((0, -2), None))
assert new_cutsite_pairs[1] == (((0, -3), None), ((6, -2), None))
assert new_cutsite_pairs[2] == (((7, -1), None), None)
