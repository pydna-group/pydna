# -*- coding: utf-8 -*-
import numpy as np
from pydivsufsort import divsufsort, kasai

string_inp = "banana$"
string_suffix_array = divsufsort(string_inp)
string_lcp_array = kasai(string_inp, string_suffix_array)
print(string_suffix_array, string_lcp_array)
# [6 5 3 1 0 4 2] [0 1 3 0 0 2 0]


from SuffixAutomaton import SuffixAutomaton, lcs1, lcs2

poet = " "
doc = poet.split()

import random

a = "".join(random.choice("ACGTURYKSWHBVDN") for i in range(100000))
b = "".join(random.choice("ACGTURYKSWHBVDN") for i in range(100000))

print(lcs1(a, b, 8))

from pydna.common_sub_strings import common_sub_strings

print(common_sub_strings(a, b, 8))
