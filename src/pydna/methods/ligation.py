# -*- coding: utf-8 -*-
"""Ligation of pre-digested fragments, and in-vivo assembly.

Both join fragments that are already compatible: ligation matches sticky (and
optionally blunt) ends, in-vivo assembly relies on the cell's homologous
recombination between shared sequence. Neither cuts anything itself — for
restriction plus ligation in one step see
:mod:`pydna.methods.restriction_ligation`.
"""

from __future__ import annotations

from pydna.assembly import (
    blunt_overlap,
    combine_algorithms,
    common_sub_strings,
    sticky_end_sub_strings,
)
from pydna.methods._engine import Method, Shape
from pydna.methods._registry import register
from pydna.opencloning_models import InVivoAssemblySource, LigationSource

__all__ = ["ligation_assembly", "in_vivo_assembly"]


def _ligation_algorithm(params: dict):
    """Sticky-end matching, widened to blunt ends when the caller allows it."""
    allow_partial_overlap = params.pop("allow_partial_overlap", False)
    allow_blunt = params.pop("allow_blunt", False)

    def sticky_end_algorithm(x, y, _l):
        return sticky_end_sub_strings(x, y, allow_partial_overlap)

    if allow_blunt:
        return combine_algorithms(sticky_end_algorithm, blunt_overlap)
    return sticky_end_algorithm


ligation_assembly = register(
    Method(
        name="ligation",
        doc="""Returns the products for ligation assembly, as inputs pass the fragments (digested if needed) that
    will be ligated.

    For most cases, you probably should use ``restriction_ligation_assembly`` instead.

    Parameters
    ----------
    frags : list[Dseqrecord]
        List of DNA fragments to assemble
    allow_blunt : bool, optional
        If True, allow blunt end ligations, by default False
    allow_partial_overlap : bool, optional
        If True, allow partial overlaps between sticky ends, by default False
    circular_only : bool, optional
        If True, only return circular assemblies, by default False

    Returns
    -------
    list[Dseqrecord]
        List of assembled DNA molecules


    Examples
    --------
    In the example below, we plan to assemble a plasmid from a backbone and an insert,
    using the EcoRI enzyme. The insert and insertion site in the backbone are flanked by
    EcoRI sites, so there are two possible products depending on the orientation of the insert.

    >>> from pydna.methods import ligation_assembly
    >>> from pydna.core.dseqrecord import Dseqrecord
    >>> from Bio.Restriction import EcoRI
    >>> backbone = Dseqrecord("cccGAATTCaaaGAATTCccc", circular=True)
    >>> backbone_cut = backbone.cut(EcoRI)[1]
    >>> insert = Dseqrecord("ggGAATTCaggtGAATTCgg")
    >>> insert_cut = insert.cut(EcoRI)[1]
    >>> products = ligation_assembly([backbone_cut, insert_cut])
    >>> products[0].seq
    Dseq(o22)
    AATTCccccccGAATTCaggtG
    TTAAGggggggCTTAAGtccaC
    >>> products[1].seq
    Dseq(o22)
    AATTCccccccGAATTCacctG
    TTAAGggggggCTTAAGtggaC
        """,
        shape=Shape.ASSEMBLY,
        algorithm_factory=_ligation_algorithm,
        source=LigationSource,
        limit=None,
        positional=("frags", "allow_blunt", "allow_partial_overlap", "circular_only"),
        summary="Ligate fragments with compatible sticky (optionally blunt) ends.",
    )
)

in_vivo_assembly = register(
    Method(
        name="in_vivo",
        shape=Shape.ASSEMBLY,
        algorithm=common_sub_strings,
        source=InVivoAssemblySource,
        positional=("frags", "limit", "circular_only"),
        summary="In-vivo assembly (IVA): the host recombines fragments sharing "
        "terminal homology.",
    )
)
