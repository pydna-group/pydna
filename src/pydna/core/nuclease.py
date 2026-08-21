#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2013-2026 Björn Johansson
# SPDX-FileCopyrightText: 2023-2026 The Project Contributors
# SPDX-License-Identifier: BSD-3-Clause

"""The nuclease interface shared by everything that cuts DNA.

A cut site may be produced by a restriction enzyme or by a programmable
nuclease, so the abstract base lives here in the core alongside the cut-site
type it appears in. Concrete nucleases are defined by the cloning method that
uses them — Cas9 in :mod:`pydna.methods.crispr`.
"""

from abc import ABC, abstractmethod
from typing import List

__all__ = ["_cas"]


class _cas(ABC):
    """
    Abstract base class for CRISPR-associated nucleases.
    pam, scaffold and cut location is set by a subclass
    such as Cas9 below.

    The meaning of size, fst5 and fst3 are the same as for the restriciton
    enzymes in the Biopython restriction module (Bio.Restriction).
    """

    scaffold: str = "ND"
    pam: str = "ND"
    size: int = 0
    fst5: int = 0
    fst3: int = 0

    def __init__(self, protospacer: str) -> None:
        """
        Initialize the nuclease with a protospacer sequence.
        The sequence is a string. Use the protospacer function
        to extract a sequence from a Dseqrecord.

        Args:
            protospacer: Protospacer sequence used to build the search pattern.
        """
        from pydna.sequence_regex import compute_regex_site

        self.protospacer: str = protospacer.upper()
        self.compsite = compute_regex_site(f"{self.protospacer}{self.pam}")

    @abstractmethod
    def search(self, dna, linear: bool = True) -> List[int]:
        """Return a list of cutting sites of the enzyme in the sequence.

        dna must be an instance of:

            - pydna.core.dseq.Dseq
            - Bio.Seq.Seq
            - Bio.Seq.MutableSeq

        pydna.core.dseqrecord.Dseqrecord or Bio.SeqRecord.SeqRecord will not work.
        This limitation is by design t omirror enzymes in the
        Biopython Bio.Restriction class

        The linear argument is laso there for compatibility with the
        Biopython Bio.Restriction class.

        An important caveat is that search ignores the circular property of
        pydna.core.dseq.Dseq.

        If linear is False, the restriction sites that span over the boundaries
        will be included.

        The positions are the first base of the 3' fragment,
        i.e. the first base after the position the enzyme will cut.
        """
        raise NotImplementedError  # pragma: no cover

    def __repr__(self) -> str:
        """
        Return a compact representation of the Cas9+gRNA nuclease instance.

        Returns:
            String representation with abbreviated protospacer.
        """
        return f"{type(self).__name__}({self.protospacer[:3]}..{self.protospacer[-3:]})"

    def __str__(self) -> str:
        """
        Return the guide RNA protospacer and scaffold as FASTA-like string.
        """
        return f">{type(self).__name__} protospacer scaffold\n{self.protospacer} {self.scaffold}"
