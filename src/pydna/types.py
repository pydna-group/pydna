# -*- coding: utf-8 -*-
"""
Types used in the pydna package.
"""

from typing import (
    TYPE_CHECKING,
    Tuple as _Tuple,
    Union as _Union,
    TypeVar as _TypeVar,
    Iterable as _Iterable,
    Callable as _Callable,
)

if TYPE_CHECKING:
    from Bio.Restriction import AbstractCut as _AbstractCut
    from Bio.Restriction import RestrictionBatch as _RestrictionBatch
    from pydna.dseq import Dseq
    from Bio.SeqFeature import Location as _Location
    from pydna.dseqrecord import Dseqrecord as _Dseqrecord


# To represent any subclass of Dseq
DseqType = _TypeVar("DseqType", bound="Dseq")
EnzymesType = _TypeVar(
    "EnzymesType", "_RestrictionBatch", _Iterable["_AbstractCut"], "_AbstractCut"
)
CutSiteType = _Tuple[_Tuple[int, int], _Union["_AbstractCut", None]]
AssemblyEdgeType = _Tuple[int, int, "_Location | None", "_Location | None"]
AssemblySubFragmentType = _Tuple[int, "_Location | None", "_Location | None"]
EdgeRepresentationAssembly = list[AssemblyEdgeType]
SubFragmentRepresentationAssembly = list[AssemblySubFragmentType]


# Type alias that describes overlap between two sequences x and y
# the two first numbers are the positions where the overlap starts on x and y
# the third number is the length of the overlap
SequenceOverlap = _Tuple[int, int, int]
AssemblyAlgorithmType = _Callable[
    ["_Dseqrecord", "_Dseqrecord", int], list[SequenceOverlap]
]
