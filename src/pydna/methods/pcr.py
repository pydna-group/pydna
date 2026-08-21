# -*- coding: utf-8 -*-
"""PCR amplification of a template between two primers.

Unlike the assembly techniques, PCR's join rule is primer annealing, which the
``PCRAssembly`` graph implements itself — so this method declares the PCR shape
rather than an overlap algorithm.
"""

from __future__ import annotations

from pydna.assembly import annotate_primer_binding_sites
from pydna.core.dseqrecord import Dseqrecord
from pydna.methods._engine import Method, Shape
from pydna.methods._registry import register
from pydna.opencloning_models import PCRSource

__all__ = ["pcr_assembly"]


def _prepare_inputs(params: dict) -> list[Dseqrecord]:
    """PCRAssembly expects the template flanked by its two primers."""
    return [
        params.pop("fwd_primer"),
        params.pop("template"),
        params.pop("rvs_primer"),
    ]


def _annotate(products, inputs, params):
    if not params.get("add_primer_features", False):
        return products
    return [annotate_primer_binding_sites(prod, inputs) for prod in products]


def _drop_redundant_orientation(products, inputs, params):
    """With identical primers both orientations are found; keep one."""
    fwd_primer, _template, rvs_primer = inputs
    if str(fwd_primer.seq).upper() == str(rvs_primer.seq).upper():
        return [p for p in products if not p.source.input[1].reverse_complemented]
    return products


pcr_assembly = register(
    Method(
        name="pcr",
        shape=Shape.PCR,
        source=PCRSource,
        limit=14,
        annotate=_annotate,
        filter_products=_drop_redundant_orientation,
        source_fields=lambda inputs, params: {
            "add_primer_features": params.get("add_primer_features", False)
        },
        prepare_inputs=_prepare_inputs,
        positional=(
            "template",
            "fwd_primer",
            "rvs_primer",
            "add_primer_features",
            "limit",
            "mismatches",
        ),
        summary="Amplify a template between a forward and a reverse primer.",
    )
)
