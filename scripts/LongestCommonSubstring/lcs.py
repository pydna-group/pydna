# -*- coding: utf-8 -*-
from time import time
import numpy as np

from pydna.common_sub_strings import common_sub_strings
from pydivsufsort import common_substrings


def profile(fun, s1, s2):
    ans = float("inf")
    for _ in range(3):
        t = time()
        solutions = len(fun(s1, s2, limit=15))
        ans = min(ans, time() - t)

    print(fun, solutions, ans)


size = 100_000
s1 = "".join(map(str, np.random.randint(0, 4, size=size, dtype=np.uint8)))
s2 = "".join(map(str, np.random.randint(0, 4, size=size, dtype=np.uint8)))

profile(common_sub_strings, s1, s2)
profile(common_substrings, s1, s2)


common_sub_strings(s1, s2, limit=15)

reps = 2_000
s1 = "banana" * reps
s2 = "ananas" * reps
profile(common_sub_strings, s1, s2)
profile(common_substrings, s1, s2)
