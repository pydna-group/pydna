# -*- coding: utf-8 -*-
"""Lyndon.py
Algorithms on strings and sequences based on Lyndon words.
David Eppstein, October 2011."""


def SmallestRotation(s):
    """Find the rotation of s that is smallest in lexicographic order.
    Duval 1983 describes how to modify his algorithm to do so but I think
    it's cleaner and more general to work from the ChenFoxLyndon output."""

    def ChenFoxLyndonBreakpoints(s):
        """Find starting positions of Chen-Fox-Lyndon decomposition of s.
        The decomposition is a set of Lyndon words that start at 0 and
        continue until the next position. 0 itself is not output, but
        the final breakpoint at the end of s is. The argument s must be
        of a type that can be indexed (e.g. a list, tuple, or string).
        The algorithm follows Duval, J. Algorithms 1983, but uses 0-based
        indexing rather than Duval's choice of 1-based indexing."""
        k = 0
        while k < len(s):
            i, j = k, k + 1
            while j < len(s) and s[i] <= s[j]:
                i = (s[i] == s[j]) and i + 1 or k  # Python cond?yes:no syntax
                j += 1
            while k < i + 1:
                k += j - i
                yield k

    def ChenFoxLyndon(s):
        """Decompose s into Lyndon words according to the Chen-Fox-Lyndon theorem.
        The arguments are the same as for ChenFoxLyndonBreakpoints but the
        return values are subsequences of s rather than indices of breakpoints."""
        old = 0
        for k in ChenFoxLyndonBreakpoints(s):
            yield s[old:k]
            old = k

    prev, rep = None, 0
    for w in ChenFoxLyndon(s + s):
        if w == prev:
            rep += 1
        else:
            prev, rep = w, 1
        if len(w) * rep == len(s):
            return w * rep
    raise Exception("Reached end of factorization with no shortest rotation")
