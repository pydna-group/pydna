# -*- coding: utf-8 -*-
"""The cloning-method engine.

A cloning method is *data*, not code. Every technique in pydna is fully
described by four things — an overlap algorithm, an assembly shape, the
provenance it records, and a few flags — plus optional hooks for the handful of
techniques that need extra behaviour.

:class:`Method` holds that description; :func:`run` executes it. Adding a
technique therefore means writing an algorithm and declaring one ``Method``,
without touching the shared assembly graph.

The pipeline is fixed, which is what makes provenance uniform::

    validate -> shape -> filter_assemblies -> annotate -> filter_products
             -> record provenance -> on_empty
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable

from pydna.core.dseqrecord import Dseqrecord

__all__ = ["Shape", "Method", "run"]


class Shape(Enum):
    """How a method's inputs are turned into candidate products.

    Each value corresponds to one of the assembly shapes implemented in
    :mod:`pydna.assembly`.
    """

    #: N fragments joined into linear and/or circular products. A single
    #: fragment is circularized, spliced or inverted.
    ASSEMBLY = "assembly"

    #: One or more inserts integrated into a (linear) genome.
    INTEGRATION = "integration"

    #: A single sequence excised or inverted at two matching sites.
    EXCISION = "excision"

    #: A template amplified between two primers.
    PCR = "pcr"


@dataclass(frozen=True)
class Method:
    """A declarative description of a cloning technique.

    Parameters
    ----------
    name : str
        Registry key, e.g. ``"gibson"``.
    shape : Shape
        Which assembly shape turns inputs into candidate products.
    algorithm : Callable
        ``(seqx, seqy, limit) -> list[SequenceOverlap]``; decides where two
        sequences may join. This is the biology of the technique.
    source : type
        Provenance class recorded on every product.
    limit : int or None, optional
        Minimum overlap length, or None when the algorithm does not use one.
    circular_only : bool, optional
        Keep only circular products (``Shape.ASSEMBLY`` only).
    only_adjacent_edges : bool, optional
        Restrict assemblies to adjacent fragments (``Shape.ASSEMBLY`` only).
    summary : str, optional
        One-line description, shown by :func:`pydna.methods.available`.

    Other Parameters
    ----------------
    validate : Callable, optional
        ``(inputs, params) -> None``; raise to reject arguments before any
        work is done.
    filter_assemblies : Callable, optional
        Predicate on *edge assemblies*, applied before products are built.
    annotate : Callable, optional
        ``(products, inputs, params) -> products``; decorate the products, e.g.
        marking the sites a recombinase acted on.
    filter_products : Callable, optional
        ``(products, inputs, params) -> products``; drop unwanted products
        after provenance has been recorded.
    source_fields : Callable, optional
        ``(inputs, params) -> dict``; extra fields stored on the ``source``.
    on_empty : Callable, optional
        ``(inputs, params) -> None``; called when no product survives, so a
        method can raise a diagnostic instead of returning an empty list.
    """

    name: str
    shape: Shape
    source: type
    algorithm: Callable | None = None
    #: For techniques whose join rule depends on the call — the enzymes used,
    #: the Gateway reaction, the recombinase — a factory ``(params) -> algorithm``
    #: is declared instead of a fixed ``algorithm``.
    algorithm_factory: Callable | None = None
    limit: int | None = 25
    circular_only: bool = False
    only_adjacent_edges: bool = False
    summary: str = ""
    #: Full docstring for the generated entry point, including any ``Examples``
    #: section. Doctests here are collected exactly as they were when the
    #: technique was a hand-written function.
    doc: str = ""

    # Hooks — every technique-specific behaviour in pydna fits one of these.
    validate: Callable | None = None
    filter_assemblies: Callable | None = None
    #: For techniques where whether to filter depends on the call — Gateway's
    #: ``multi_site_only`` — a factory ``(params) -> filter or None``.
    filter_assemblies_factory: Callable | None = None
    annotate: Callable | None = None
    filter_products: Callable | None = None
    source_fields: Callable | None = None
    on_empty: Callable | None = None

    #: Parameters that define this method's identity, stored on the provenance
    #: record so a strategy can be replayed without guessing them back.
    params: dict[str, Any] = field(default_factory=dict)

    #: Names of the parameters callers may pass positionally, in order, so a
    #: declared method reads like the hand-written function it replaces —
    #: ``gibson_assembly(frags, 25, True)``.
    positional: tuple[str, ...] = ()

    #: Turns the call's arguments into the flat list of sequences the shape
    #: works on, consuming the arguments it uses. Techniques that take a genome
    #: and inserts declare one; the default takes the first positional argument.
    prepare_inputs: Callable | None = None

    def resolve_algorithm(self, params: dict) -> Callable:
        """Return the join rule for a call, from the factory if one is declared."""
        if self.algorithm_factory is not None:
            return self.algorithm_factory(params)
        if self.algorithm is None:
            raise ValueError(
                f"Method {self.name!r} declares neither algorithm nor algorithm_factory"
            )
        return self.algorithm

    def resolve_filter_assemblies(self, params: dict) -> Callable | None:
        """Return the assembly filter for a call, if the method declares one."""
        if self.filter_assemblies_factory is not None:
            return self.filter_assemblies_factory(params)
        return self.filter_assemblies

    def resolve_inputs(self, params: dict) -> list[Dseqrecord]:
        """Extract the sequences to work on, consuming them from *params*."""
        if self.prepare_inputs is not None:
            return self.prepare_inputs(params)
        return params.pop(self.positional[0])


def _build(method: Method, inputs: list[Dseqrecord], params: dict) -> list[Dseqrecord]:
    """Turn inputs into candidate products using the method's shape."""
    # Imported here because pydna.assembly imports the provenance models, which
    # in turn reach back into the methods package.
    from pydna.assembly import (
        common_function_assembly_products,
        common_function_excision_or_inversion_products,
        common_function_integration_products,
    )

    limit = params.get("limit", method.limit)

    if method.shape is Shape.PCR:
        # PCRAssembly owns its own join rule (primer annealing), so no
        # algorithm is declared; mismatches widen the required annealing.
        from pydna.assembly import PCRAssembly

        mismatches = params.get("mismatches", 0)
        asm = PCRAssembly(inputs, limit=limit + mismatches, mismatches=mismatches)
        return asm.assemble_linear()

    algorithm = method.resolve_algorithm(params)

    if method.shape is Shape.ASSEMBLY:
        return common_function_assembly_products(
            inputs,
            limit,
            algorithm,
            params.get("circular_only", method.circular_only),
            method.resolve_filter_assemblies(params),
            method.only_adjacent_edges,
        )
    if method.shape is Shape.INTEGRATION:
        return common_function_integration_products(inputs, limit, algorithm)
    if method.shape is Shape.EXCISION:
        (genome,) = inputs
        return common_function_excision_or_inversion_products(genome, limit, algorithm)
    raise ValueError(f"Unknown shape: {method.shape}")  # pragma: no cover


def run(method: Method, inputs: list[Dseqrecord], **params: Any) -> list[Dseqrecord]:
    """Execute *method* over *inputs* and return provenance-tagged products.

    The pipeline is the same for every technique, which is what keeps history
    uniform: whatever a method does, its products always carry the ``source``
    declared on the :class:`Method`.
    """
    from pydna.assembly import _recast_sources

    if method.validate is not None:
        method.validate(inputs, params)

    products = _build(method, inputs, params)

    if method.annotate is not None:
        products = method.annotate(products, inputs, params)

    extra = method.source_fields(inputs, params) if method.source_fields else {}
    products = _recast_sources(products, method.source, **extra)

    if method.filter_products is not None:
        products = method.filter_products(products, inputs, params)

    if not products and method.on_empty is not None:
        method.on_empty(inputs, params)

    return products
