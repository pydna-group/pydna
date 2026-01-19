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
from Bio.Restriction.Restriction import AbstractCut
from pydna.crispr import _cas

if TYPE_CHECKING:
    from Bio.Restriction import RestrictionBatch
    from pydna.dseq import Dseq
    from Bio.SeqFeature import Location
    from pydna.dseqrecord import Dseqrecord


# To represent any subclass of Dseq
DseqType = TypeVar("DseqType", bound="Dseq")
EnzymesType = TypeVar(
    "EnzymesType", "RestrictionBatch", Iterable["AbstractCut"], "AbstractCut"
)
CutSiteType = Tuple[Tuple[int, int], Union[AbstractCut, None, _cas]]
AssemblyEdgeType = Tuple[int, int, "Location | None", "Location | None"]
AssemblySubFragmentType = Tuple[int, "Location | None", "Location | None"]
EdgeRepresentationAssembly = list[AssemblyEdgeType]
SubFragmentRepresentationAssembly = list[AssemblySubFragmentType]


# Type alias that describes overlap between two sequences x and y
# the two first numbers are the positions where the overlap starts on x and y
# the third number is the length of the overlap
SequenceOverlap = Tuple[int, int, int]
AssemblyAlgorithmType = Callable[
    ["Dseqrecord", "Dseqrecord", int], list[SequenceOverlap]
]
