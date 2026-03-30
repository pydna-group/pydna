#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Utilities for CRISPR/Cas target searching and protospacer extraction.

"""
import re
from abc import ABC
from abc import abstractmethod
from typing import Type
from typing import TYPE_CHECKING
from typing import List
from typing import TypeVar

if TYPE_CHECKING:  # pragma: no cover
    from pydna.dseqrecord import Dseqrecord

DseqrecordType = TypeVar("DseqrecordType", bound="Dseqrecord")


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

            - pydna.dseq.Dseq
            - Bio.Seq.Seq
            - Bio.Seq.MutableSeq

        pydna.dseqrecord.Dseqrecord or Bio.SeqRecord.SeqRecord will not work.
        This limitation is by design t omirror enzymes in the
        Biopython Bio.Restriction class

        The linear argument is laso there for compatibility with the
        Biopython Bio.Restriction class.

        An important caveat is that search ignores the circular property of
        pydna.dseq.Dseq.

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
            dna: string, Bio.Seq.Seq or pydna.dseq.Dseq
            linear: Whether the DNA is linear or circular.

        Returns:
            A list of cut site positions.
        """
        from pydna.dseqrecord import Dseqrecord
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
