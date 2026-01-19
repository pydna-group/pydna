#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""This module is based on the Py-rstr-max package that
was written by Romain Brixtel (rbrixtel_at_gmail_dot_com)
(https://brixtel.users.greyc.fr) and is available from
https://code.google.com/p/py-rstr-max
https://github.com/gip0/py-rstr-max
the original code was covered by an MIT licence."""


from operator import itemgetter
from typing import List, Tuple

Match = Tuple[int, int, int]  # (x_start, y_start, length)

# def _kark_sort(s, SA, n, K):
#     def radixpass(a, b, r, s, n, k):
#         c = _array("i", [0] * (k + 1))
#         for i in range(n):
#             c[r[a[i] + s]] += 1

#         somme = 0
#         for i in range(k + 1):
#             freq, c[i] = c[i], somme
#             somme += freq

#         for i in range(n):
#             b[c[r[a[i] + s]]] = a[i]
#             c[r[a[i] + s]] += 1

#     n0 = (n + 2) // 3
#     n1 = (n + 1) // 3
#     n2 = n // 3
#     n02 = n0 + n2

#     SA12 = _array("i", [0] * (n02 + 3))
#     SA0 = _array("i", [0] * n0)

#     s12 = [i for i in range(n + (n0 - n1)) if i % 3]
#     s12.extend([0] * 3)
#     s12 = _array("i", s12)

#     radixpass(s12, SA12, s, 2, n02, K)
#     radixpass(SA12, s12, s, 1, n02, K)
#     radixpass(s12, SA12, s, 0, n02, K)

#     name = 0
#     c0, c1, c2 = -1, -1, -1
#     for i in range(n02):
#         if s[SA12[i]] != c0 or s[SA12[i] + 1] != c1 or s[SA12[i] + 2] != c2:
#             name += 1
#             c0 = s[SA12[i]]
#             c1 = s[SA12[i] + 1]
#             c2 = s[SA12[i] + 2]
#         if SA12[i] % 3 == 1:
#             s12[SA12[i] // 3] = name
#         else:
#             s12[SA12[i] // 3 + n0] = name

#     if name < n02:
#         _kark_sort(s12, SA12, n02, name + 1)
#         for i in range(n02):
#             s12[SA12[i]] = i + 1
#     else:
#         for i in range(n02):
#             SA12[s12[i] - 1] = i

#     s0 = _array("i", [SA12[i] * 3 for i in range(n02) if SA12[i] < n0])
#     radixpass(s0, SA0, s, 0, n0, K)

#     p = j = k = 0
#     t = n0 - n1
#     while k < n:
#         i = SA12[t] * 3 + 1 if SA12[t] < n0 else (SA12[t] - n0) * 3 + 2
#         j = SA0[p] if p < n0 else 0

#         if SA12[t] < n0:
#             test = (s12[SA12[t] + n0] <= s12[j // 3]) if (s[i] == s[j]) else (s[i] < s[j])
#         elif s[i] == s[j]:
#             test = s12[SA12[t] - n0 + 1] <= s12[j // 3 + n0] if (s[i + 1] == s[j + 1]) else s[i + 1] < s[j + 1]
#         else:
#             test = s[i] < s[j]
#         if test:
#             SA[k] = i
#             t += 1
#             if t == n02:
#                 k += 1
#                 while p < n0:
#                     SA[k] = SA0[p]
#                     p += 1
#                     k += 1
#         else:
#             SA[k] = j
#             p += 1
#             if p == n0:
#                 k += 1
#                 while t < n02:
#                     SA[k] = (SA12[t] * 3) + 1 if SA12[t] < n0 else ((SA12[t] - n0) * 3) + 2
#                     t += 1
#                     k += 1
#         k += 1


# class _Rstr_max:
#     def __init__(self):
#         self.array_str = []

#     def add_str(self, str_unicode):
#         self.array_str.append(str_unicode)

#     def _step1_sort_suffix(self):
#         char_frontier = chr(2)

#         self.global_suffix = char_frontier.join(self.array_str)

#         nbChars = len(self.global_suffix)
#         init = [-1] * nbChars
#         self.idxString = _array("i", init)
#         self.idxPos = _array("i", init)
#         self.endAt = _array("i", init)

#         k = idx = 0
#         for mot in self.array_str:
#             last = k + len(mot)
#             for p in range(len(mot)):
#                 self.idxString[k] = idx
#                 self.idxPos[k] = p
#                 self.endAt[k] = last
#                 k += 1
#             idx += 1
#             k += 1

#         s = self.global_suffix
#         alphabet = [None] + sorted(set(s))
#         k = len(alphabet)
#         n = len(s)
#         t = dict((c, i) for i, c in enumerate(alphabet))
#         SA = _array("i", [0] * (n + 3))
#         _kark_sort(_array("i", [t[c] for c in s] + [0] * 3), SA, n, k)
#         self.res = SA[:n]

#     def _step2_lcp(self):
#         n = len(self.res)
#         init = [0] * n
#         rank = _array("i", init)
#         LCP = _array("i", init)

#         s = self.global_suffix
#         suffix__array = self.res
#         endAt = self.endAt

#         for i in range(len(self.array_str), n):
#             v = self.res[i]
#             rank[v] = i

#         l = 0
#         for j in range(n):
#             if l > 0:
#                 l -= 1
#             i = rank[j]
#             j2 = suffix__array[i - 1]
#             if i:
#                 while l + j < endAt[j] and l + j2 < endAt[j2] and s[j + l] == s[j2 + l]:
#                     l += 1
#                 LCP[i - 1] = l
#             else:
#                 l = 0
#         self.lcp = LCP

#     def _step3_rstr(self):
#         prev_len = 0
#         idx = 0
#         results = {}
#         len_lcp = len(self.lcp) - 1

#         class Stack:
#             pass

#         stack = Stack()
#         stack._top = 0
#         stack.lst_max = []

#         # if len(self.res) == 0 :
#         #  return {}

#         pos1 = self.res[0]
#         for idx in range(len_lcp):
#             current_len = self.lcp[idx]
#             pos2 = self.res[idx + 1]
#             end_ = max(pos1, pos2) + current_len
#             n = prev_len - current_len
#             if n < 0:
#                 # pushMany
#                 stack.lst_max.append([-n, idx, end_])
#                 stack._top += -n
#             elif n > 0:
#                 self.removeMany(stack, results, n, idx)
#             elif stack._top > 0 and end_ > stack.lst_max[-1][-1]:
#                 # setMax
#                 stack.lst_max[-1][-1] = end_

#             prev_len = current_len
#             pos1 = pos2

#         if stack._top > 0:
#             self.removeMany(stack, results, stack._top, idx + 1)

#         return results

#     def removeMany(self, stack, results, m, idxEnd):
#         prevStart = -1
#         while m > 0:
#             n, idxStart, maxEnd = stack.lst_max.pop()
#             if prevStart != idxStart:
#                 # idStr = self.idxString[maxEnd-1]
#                 # pos = self.idxPos[maxEnd-1]
#                 id_ = (maxEnd, idxEnd - idxStart + 1)
#                 if id_ not in results or results[id_][0] < stack._top:
#                     results[id_] = (stack._top, idxStart)
#                 prevStart = idxStart
#             m -= n
#             stack._top -= n
#         if m < 0:
#             stack.lst_max.append([-m, idxStart, maxEnd - n - m])
#             stack._top -= m

#     def go(self):
#         self._step1_sort_suffix()
#         self._step2_lcp()
#         return self._step3_rstr()


# def common_sub_strings_py(stringx: str, stringy: str, limit=25):
#     """Finds all common substrings between stringx and stringy
#     longer than limit. This function is case sensitive.
#     The substrings may overlap.

#     returns a list of tuples describing the substrings
#     The list is sorted longest -> shortest.

#     Parameters
#     ----------
#     stringx : str
#     stringy : str
#     limit : int, optional

#     Returns
#     -------
#     list of tuple
#         [(startx1, starty1, length1),(startx2, starty2, length2), ...]

#         startx1 = startposition in x, where substring 1 starts
#         starty1 = position in y where substring 1 starts
#         length1 = lenght of substring


#     Examples
#     --------

#     >>> from pydna.common_sub_strings import common_sub_strings
#     >>> common_sub_strings("gatgatttcggtagtta", "gtcagtatgtctatctatcgcg", limit=3)
#     [(1, 6, 3), (7, 17, 3), (10, 4, 3), (12, 3, 3)]

#     ::

#         Overlaps   Symbols
#         (1, 6,  3)   ---
#         (7, 17, 3)   +++
#         (10, 4, 3)   ...
#         (12, 3, 3)   ===


#                     ===
#         gatgatttcggtagtta           stringx
#          ---   +++...

#             ...
#         gtcagtatgtctatctatcgcg      stringy
#            ===---        +++

#     """
#     rstr = _Rstr_max()
#     rstr.add_str("&".join((stringx, stringy)))
#     r = rstr.go()
#     match = {}
#     for (offset_end, nb), (l, start_plage) in r.items():
#         if l < limit:
#             continue
#         startsx = []
#         startsy = []
#         for o in range(start_plage, start_plage + nb):
#             offset = rstr.idxPos[rstr.res[o]]
#             if offset > len(stringx):
#                 startsy.append(offset - len(stringx) - 1)
#             else:
#                 startsx.append(offset)

#         for a, b in _itertools.product(startsx, startsy):
#             match[(a, b)] = max(match.get((a, b)) or 0, l)

#     match = [(key[0], key[1], val) for key, val in list(match.items())]

#     match.sort()

#     match.sort(key=_itemgetter(2), reverse=True)

#     return match


def common_sub_strings(stringx: str, stringy: str, limit: int = 25) -> List[Match]:
    """
    Finds all common substrings between stringx and stringy, and returns
    them sorted by length.

    This function is case sensitive.

    Parameters
    ----------
    stringx : str
    stringy : str
    limit : int, optional

    Returns
    -------
    list of tuple
        [(startx1, starty1, length1),(startx2, starty2, length2), ...]

        startx1 = startposition in x, where substring 1 starts
        starty1 = position in y where substring 1 starts
        length1 = lenght of substring
    """
    from pydivsufsort import common_substrings

    matches = common_substrings(stringx, stringy, limit=limit)
    matches.sort()
    matches.sort(key=itemgetter(2), reverse=True)
    return matches


def terminal_overlap(stringx: str, stringy: str, limit: int = 15) -> List[Match]:
    """Finds the the flanking common substrings between stringx and stringy
    longer than limit. This means that the results only contains substrings
    that starts or ends at the the ends of stringx and stringy.

    This function is case sensitive.

    returns a list of tuples describing the substrings
    The list is sorted longest -> shortest.

    Parameters
    ----------
    stringx : str
    stringy : str
    limit : int, optional

    Returns
    -------
    list of tuple
        [(startx1,starty1,length1),(startx2,starty2,length2), ...]

        startx1 = startposition in x, where substring 1 starts
        starty1 = position in y where substring 1 starts
        length1 = lenght of substring


    Examples
    --------

    >>> from pydna.common_sub_strings import terminal_overlap
    >>> terminal_overlap("agctatgtatcttgcatcgta", "gcatcgtagtctatttgcttac", limit=8)
    [(13, 0, 8)]

    ::

                        <-- 8 ->
           <---- 13 --->
           agctatgtatcttgcatcgta                    stringx
                        gcatcgtagtctatttgcttac      stringy
                        0

    """
    return [
        m
        for m in common_sub_strings(stringx, stringy, limit)
        if (m[0] == 0 and m[1] + m[2] == len(stringy))
        or (m[1] == 0 and m[0] + m[2] == len(stringx))
    ]
