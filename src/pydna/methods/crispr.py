#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2013-2026 Björn Johansson
# SPDX-FileCopyrightText: 2023-2026 The Project Contributors
# SPDX-License-Identifier: BSD-3-Clause

"""
Utilities for CRISPR/Cas target searching and protospacer extraction.

"""

import re
from typing import Type
from typing import TYPE_CHECKING
from typing import List
from typing import TypeVar

if TYPE_CHECKING:  # pragma: no cover
    from pydna.core.dseqrecord import Dseqrecord
from pydna.core.nuclease import _cas

DseqrecordType = TypeVar("DseqrecordType", bound="Dseqrecord")


class cas9(_cas):
    """docstring.

    .. code-block::


            fst5              --|fst3
            |----------------
                                 PAM
       5'-NNGGAAGAGTAATACACTA-AAANGGNN-3'
          ||||||||||||||||||| ||||||||
       3'-NNCCTTCTCATTATGTGAT-TTTNCCNN-5'
            ||||||||||||||||| |||
         5'-GGAAGAGTAATACACTA AAA-g-u-a-a-g-g-3'  Scaffold (lower case)
            ---gRNA spacer------- u-a
                                  u-a
                                  u-a
                                  u-a
                                  a-u
                                  g-u-g
                                  a    a
                                  g-c-a
                                  c-g
                                  u-a
                                  a-u
                                  g   a
                                   a-a
    """

    scaffold: str = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG"
    pam: str = ".GG"
    size: int = 20
    fst5: int = 17
    fst3: int = -3
    ovhg: int = fst5 - (size + fst3)

    def search(self, dna, linear: bool = True) -> List[int]:
        """
        Search for Cas9 target sites in a DNA sequence.

        Args:
            dna: string, Bio.Seq.Seq or pydna.core.dseq.Dseq
            linear: Whether the DNA is linear or circular.

        Returns:
            A list of cut site positions.
        """
        from pydna.core.dseqrecord import Dseqrecord
        from pydna.sequence_regex import dseqrecord_finditer

        if not hasattr(dna, "_data"):
            raise TypeError
        results: List[int] = []
        query = Dseqrecord(dna, circular=(not linear))
        matches_fwd = dseqrecord_finditer(self.compsite, query)
        matches_rev = dseqrecord_finditer(self.compsite, query.reverse_complement())
        for mobj in matches_fwd:
            results.append((mobj.start() + self.fst5 + 1) % len(dna))
        for mobj in matches_rev:
            results.append((len(dna) - (mobj.start() + self.fst5) + 1) % len(dna))
        return results


def protospacer(guide_construct: DseqrecordType, cas: Type[_cas] = cas9) -> List[str]:
    """
    Extract protospacer sequences from a guide construct. This can for example
    be a plasmid containing the guide construct. This function returns a
    list since a several protospacers can be present.

    Args:
        guide_construct: Sequence construct containing protospacer and scaffold.
        cas: CRISPR nuclease class defining spacer size and scaffold.

    Returns:
        A list of protospacer sequences found in Watson and Crick orientations.
    """
    if guide_construct.circular:
        total_length = cas.size + len(cas.scaffold)
        guide_construct = guide_construct[:] + guide_construct[: total_length - 1]

    result = []

    for s in guide_construct.seq.watson.upper(), guide_construct.seq.crick.upper():
        result.extend(
            mobj.group("ps")
            for mobj in re.finditer(
                f"(?P<ps>.{{{cas.size}}})(?:{cas.scaffold})",
                s,
            )
        )

    return result
