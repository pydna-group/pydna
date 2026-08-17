# -*- coding: utf-8 -*-
"""pydna.core — the fundamental sequence data model.

This package holds pydna's core sequence types and their primitives:

    - :class:`pydna.core.seq.Seq`
    - :class:`pydna.core.dseq.Dseq`
    - :class:`pydna.core.seqrecord.SeqRecord`
    - :class:`pydna.core.dseqrecord.Dseqrecord`
    - :mod:`pydna.core.alphabet`  (double-stranded encoding tables)
    - :mod:`pydna.core.types`     (shared type aliases)

Import them from their submodules, e.g. ``from pydna.core.dseq import Dseq``.
Code written against the pydna 5.x layout can keep the old behaviour by importing
from the :mod:`pydna.legacy` shims instead, e.g. ``from pydna.legacy.dseq import Dseq``.
"""
