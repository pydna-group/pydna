#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""docstring."""
from operator import add
from functools import reduce
import networkx as nx
from itertools import permutations


def ligate(fragments: list):
    """docstring."""
    G = nx.DiGraph()
    G.add_nodes_from(["begin", "end"])
    fragments = fragments[:]

    fragments.extend(f.rc() for f in fragments[1:])

    for node in fragments:
        G.add_edge("begin", node)
        G.add_edge(node, "end")

    for seq1, seq2 in permutations(fragments, 2):
        try:
            seq1 + seq2
        except TypeError as err:
            if not str(err).startswith("sticky ends not compatible"):
                raise
        else:
            if seq1.seq.three_prime_end() != (
                "blunt",
                "",
            ) and seq2.seq.five_prime_end() != ("blunt", ""):
                G.add_edge(seq1, seq2)
                try:
                    G.remove_edge("begin", seq2)
                except nx.NetworkXError as err:
                    if "not in graph" not in str(err):
                        raise
                try:
                    G.remove_edge(seq1, "end")
                except nx.NetworkXError as err:
                    if "not in graph" not in str(err):
                        raise

    cpaths = [p for p in sorted(nx.simple_cycles(G), key=len) if len(p) > 1]
    csequences = [reduce(add, x).looped() for x in cpaths]
    lpaths = [
        p for p in sorted(nx.all_simple_paths(G, "begin", "end"), key=len) if len(p) > 3
    ]
    lsequences = [reduce(add, lp[1:-1]) for lp in lpaths]

    return csequences, lsequences
