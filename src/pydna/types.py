# -*- coding: utf-8 -*-
"""
Types used in the pydna package.
"""

from typing import (
    TYPE_CHECKING,
    Tuple,
    Union,
    TypeVar,
    Iterable,
    Callable,
)

# Import AbstractCut at runtime for CutSiteType
from Bio.Restriction.Restriction import AbstractCut as _AbstractCut
from pydna.crispr import _cas as __cas

if TYPE_CHECKING:
    from Bio.Restriction import RestrictionBatch as _RestrictionBatch
    from pydna.dseq import Dseq
    from Bio.SeqFeature import Location as _Location
    from pydna.dseqrecord import Dseqrecord as _Dseqrecord


# To represent any subclass of Dseq
DseqType = TypeVar("DseqType", bound="Dseq")
EnzymesType = TypeVar(
    "EnzymesType", "_RestrictionBatch", Iterable["_AbstractCut"], "_AbstractCut"
)
CutSiteType = Tuple[Tuple[int, int], Union[_AbstractCut, None, __cas]]
AssemblyEdgeType = Tuple[int, int, "_Location | None", "_Location | None"]
AssemblySubFragmentType = Tuple[int, "_Location | None", "_Location | None"]
EdgeRepresentationAssembly = list[AssemblyEdgeType]
SubFragmentRepresentationAssembly = list[AssemblySubFragmentType]


# Type alias that describes overlap between two sequences x and y
# the two first numbers are the positions where the overlap starts on x and y
# the third number is the length of the overlap
SequenceOverlap = Tuple[int, int, int]
AssemblyAlgorithmType = Callable[
    ["_Dseqrecord", "_Dseqrecord", int], list[SequenceOverlap]
]
