# -*- coding: utf-8 -*-
"""pydna.methods — cloning techniques as declarations over a shared engine.

A cloning method is described by a :class:`~pydna.methods._engine.Method`: an
overlap algorithm, an assembly shape, the provenance it records, and optional
hooks. One engine runs them all, so every technique records history the same
way and adding a technique never touches the assembly graph.

    >>> from pydna.methods import available, get
    >>> "gibson" in available()
    True
    >>> get("gibson").summary
    "Gibson assembly: join fragments by terminal overlaps exposed by 5' chew-back."

Each technique is also exposed under its conventional name, e.g.
:func:`gibson_assembly`, which is the name its provenance ``Source`` is keyed to.
"""

from pydna.methods._engine import Method, Shape, run
from pydna.methods._registry import (
    register,
    get,
    available,
    methods,
    for_source,
    entrypoint_for,
)

# Importing the names below pulls in every technique module, so all
# declarations are registered by the time this package finishes importing.

from pydna.methods.overlap import (
    gibson_assembly,
    in_fusion_assembly,
    fusion_pcr_assembly,
    overlap_assembly,
)
from pydna.methods.ligation import ligation_assembly, in_vivo_assembly
from pydna.methods.restriction_ligation import (
    restriction_ligation_assembly,
    golden_gate_assembly,
)
from pydna.methods.gateway import gateway_assembly
from pydna.methods.crispr import crispr_integration
from pydna.methods.pcr import pcr_assembly
from pydna.methods.homologous_recombination import (
    homologous_recombination_integration,
    homologous_recombination_excision_or_inversion,
    homologous_recombination_excision,
)
from pydna.methods.cre_lox import (
    cre_lox_integration,
    cre_lox_excision_or_inversion,
    cre_lox_excision,
)
from pydna.methods.recombinase import (
    recombinase_integration,
    recombinase_excision_or_inversion,
    recombinase_excision,
    recombinase_assembly,
)

__all__ = [
    # engine
    "Method",
    "Shape",
    "run",
    # registry
    "register",
    "get",
    "available",
    "methods",
    "for_source",
    "entrypoint_for",
    # techniques
    "gibson_assembly",
    "in_fusion_assembly",
    "fusion_pcr_assembly",
    "overlap_assembly",
    "ligation_assembly",
    "in_vivo_assembly",
    "restriction_ligation_assembly",
    "golden_gate_assembly",
    "gateway_assembly",
    "homologous_recombination_integration",
    "homologous_recombination_excision_or_inversion",
    "homologous_recombination_excision",
    "cre_lox_integration",
    "cre_lox_excision_or_inversion",
    "cre_lox_excision",
    "recombinase_integration",
    "recombinase_excision_or_inversion",
    "recombinase_excision",
    "recombinase_assembly",
    "crispr_integration",
    "pcr_assembly",
]
