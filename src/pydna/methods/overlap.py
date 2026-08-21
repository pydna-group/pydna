# -*- coding: utf-8 -*-
"""Overlap-based assembly: Gibson, In-Fusion and overlap-extension PCR.

These techniques all join fragments through terminal sequence overlaps and
differ only in which strand end is chewed back to expose the overlap — and in
the provenance they record. Each is therefore a one-line
:class:`~pydna.methods._engine.Method` declaration over the shared engine.
"""

from __future__ import annotations

from functools import partial
from typing import Literal

from pydna.assembly import terminal_overlap
from pydna.core.dseqrecord import Dseqrecord
from pydna.methods._engine import Method, Shape
from pydna.methods._registry import register
from pydna.opencloning_models import (
    GibsonAssemblySource,
    InFusionSource,
    OverlapExtensionPCRLigationSource,
)

__all__ = [
    "gibson_assembly",
    "in_fusion_assembly",
    "fusion_pcr_assembly",
    "overlap_assembly",
]


gibson_assembly = register(
    Method(
        name="gibson",
        shape=Shape.ASSEMBLY,
        algorithm=partial(terminal_overlap, trim_ends="5'"),
        source=GibsonAssemblySource,
        params={"ends": "5'"},
        positional=("frags", "limit", "circular_only"),
        summary="Gibson assembly: join fragments by terminal overlaps exposed "
        "by 5' chew-back.",
    )
)

in_fusion_assembly = register(
    Method(
        name="in_fusion",
        shape=Shape.ASSEMBLY,
        algorithm=partial(terminal_overlap, trim_ends="3'"),
        source=InFusionSource,
        params={"ends": "3'"},
        positional=("frags", "limit", "circular_only"),
        summary="In-Fusion cloning: as Gibson, but the overlap is exposed by "
        "3' chew-back.",
    )
)

fusion_pcr_assembly = register(
    Method(
        name="fusion_pcr",
        shape=Shape.ASSEMBLY,
        algorithm=partial(terminal_overlap, trim_ends=None),
        source=OverlapExtensionPCRLigationSource,
        params={"ends": None},
        positional=("frags", "limit", "circular_only"),
        summary="Overlap-extension PCR: fragments joined by terminal overlaps "
        "with no chew-back.",
    )
)


_BY_ENDS = {
    "5'": gibson_assembly,
    "3'": in_fusion_assembly,
    None: fusion_pcr_assembly,
}


def overlap_assembly(
    frags: list[Dseqrecord],
    *,
    ends: Literal["5'", "3'"] | None = "5'",
    **params,
) -> list[Dseqrecord]:
    """Assemble fragments sharing terminal overlaps, selecting the technique.

    A convenience wrapper over the three named techniques for callers that want
    to pick one by parameter rather than by name.

    Parameters
    ----------
    frags : list[Dseqrecord]
        Fragments to assemble.
    ends : {"5'", "3'", None}, optional
        Which end is chewed back: ``"5'"`` for Gibson (default), ``"3'"`` for
        In-Fusion, ``None`` for overlap-extension PCR.
    **params
        Forwarded to the technique, e.g. ``limit`` or ``circular_only``.

    Returns
    -------
    list[Dseqrecord]
        Assembled molecules, each tagged with the matching assembly ``Source``.
    """
    try:
        technique = _BY_ENDS[ends]
    except KeyError:
        raise ValueError(f"""ends must be "5'", "3'" or None, not {ends!r}""") from None
    return technique(frags, **params)
