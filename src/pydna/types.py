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
)

if TYPE_CHECKING:
    from Bio.Restriction import AbstractCut as _AbstractCut
    from Bio.Restriction import RestrictionBatch as _RestrictionBatch
    from pydna.dseq import Dseq


# To represent any subclass of Dseq
DseqType = _TypeVar("DseqType", bound="Dseq")
EnzymesType = _TypeVar("EnzymesType", "_RestrictionBatch", _Iterable["_AbstractCut"], "_AbstractCut")
CutSiteType = _Tuple[_Tuple[int, int], _Union["_AbstractCut", None]]
