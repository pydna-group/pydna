#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pydna.common_sub_strings import terminal_overlap
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.utils import rc
from itertools import combinations
from copy import copy


def fuse_by_pcr(fragments, limit=15):
    """docstring."""

    def anneal(x, y, limit=limit):
        """docstring."""
        new = None
        for a, b in [(x, y), (x, y.rc()), (x.rc(), y)]:
            try:
                ((s1, s2, ln), *r) = terminal_overlap(
                    a.seq.watson.lower(), rc(b.seq.crick.lower()), limit=limit
                )
            except ValueError as err:
                if "not enough values to unpack" not in str(err):
                    raise err
            else:
                if s2 == 0:
                    new = Dseqrecord(Dseq(a.seq.watson, b.seq.crick, ovhg=-s1))
                    new.features = copy(a.features)
                    new.features.extend([f._shift(s1) for f in b.features])
                    new.left = b  # Setting a new property dynamically
                    new.right = a  # Setting a new property dynamically
                elif s1 == 0:
                    new = Dseqrecord(Dseq(b.seq.watson, a.seq.crick, ovhg=-s2))
                    new.features = copy(b.features)
                    new.features.extend([f._shift(s2) for f in a.features])
                    new.left = a  # Setting a new property dynamically
                    new.right = b  # Setting a new property dynamically
        return new

    argument = fragments
    for arg in argument:
        arg.left = None
        arg.right = None
    newfragments = []
    while True:
        for a, b in combinations(fragments, 2):
            new = anneal(a, b, limit)
            if not new:
                continue
            new.anneal = new.seq
            new.seq = new.seq.fill_in()
            new.features = list({repr(f): f for f in new.features}.values())
            newfragments.append(new)
        if newfragments:
            fragments = newfragments
            newfragments = []
        else:
            break

    return [x for x in fragments if x not in argument]


def list_parts(fusion_pcr_fragment):
    stack = [fusion_pcr_fragment]
    processed = []

    while stack:
        obj = stack.pop()
        try:
            a, b = obj.right, obj.left
        except AttributeError:
            pass
        else:
            stack.extend((a, b))
        processed.append(obj)

    parts = [x for x in processed[::-1] if x]

    msg = "---\n\n"

    for part in parts:
        if hasattr(part, "anneal"):
            msg += repr(part.anneal) + "\n\n"
            msg += f"{part.name}\n{part.seq.watson}\n{part.seq.crick[::-1]}\n\n---\n\n"
        msg += repr(part.seq) + "\n\n"

    return msg
