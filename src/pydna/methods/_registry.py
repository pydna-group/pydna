# -*- coding: utf-8 -*-
"""Registry of cloning methods.

Every :class:`~pydna.methods._engine.Method` declared anywhere in
:mod:`pydna.methods` registers itself here, so the set of available techniques
is discoverable at runtime rather than hardcoded in a table:

    >>> from pydna.methods import available, get
    >>> "gibson" in available()
    True
    >>> get("gibson").source.__name__
    'GibsonAssemblySource'

Adding a technique means declaring a ``Method`` in its own module — nothing
central needs editing.
"""

from __future__ import annotations

import sys
from typing import Callable

from pydna.core.dseqrecord import Dseqrecord
from pydna.methods._engine import Method, run

__all__ = ["register", "get", "available", "methods", "for_source", "entrypoint_for"]

_REGISTRY: dict[str, Method] = {}
_ENTRYPOINTS: dict[str, Callable[..., list[Dseqrecord]]] = {}


def register(method: Method) -> Callable[..., list[Dseqrecord]]:
    """Register *method* and return a callable that runs it.

    The returned function is the technique's public entry point::

        gibson_assembly = register(Method(name="gibson", ...))
        products = gibson_assembly(fragments, limit=25)
    """
    if method.name in _REGISTRY:
        raise ValueError(f"A method named {method.name!r} is already registered")
    _REGISTRY[method.name] = method

    def entrypoint(*args, **params) -> list[Dseqrecord]:
        # Positional arguments are mapped onto the names the method declares,
        # so a declaration is a drop-in for the function it replaces.
        if len(args) > len(method.positional):
            raise TypeError(
                f"{method.name}() takes at most {len(method.positional)} "
                f"positional arguments but {len(args)} were given"
            )
        for name, value in zip(method.positional, args):
            if name in params:
                raise TypeError(
                    f"{method.name}() got multiple values for argument {name!r}"
                )
            params[name] = value
        inputs = method.resolve_inputs(params)
        return run(method, inputs, **params)

    entrypoint.__name__ = method.name
    entrypoint.__qualname__ = method.name
    entrypoint.__doc__ = (
        method.doc or method.summary or f"Run the {method.name} cloning method."
    )
    # Claim the declaring module so the docstring's examples are collected as
    # doctests, exactly as they were when this was a hand-written function.
    entrypoint.__module__ = sys._getframe(1).f_globals.get("__name__", __name__)
    entrypoint.method = method
    _ENTRYPOINTS[method.name] = entrypoint
    return entrypoint


def get(name: str) -> Method:
    """Return the :class:`Method` registered under *name*."""
    return _REGISTRY[name]


def for_source(source) -> Method:
    """Return the :class:`Method` that records *source* as its provenance.

    *source* may be a provenance class or an instance of one. This is the
    lookup a consumer needs to go from a recorded ``Source`` back to the
    technique that produced it, instead of maintaining its own table::

        method = for_source(GibsonAssemblySource)

    Raises
    ------
    KeyError
        If no registered method records that provenance class.
    """
    cls = source if isinstance(source, type) else type(source)
    # Prefer an exact match; fall back to a subclass relationship so that a
    # specialised source (e.g. CRISPRSource from HomologousRecombinationSource)
    # resolves to its own method rather than its parent's.
    for method in _REGISTRY.values():
        if method.source is cls:
            return method
    for method in _REGISTRY.values():
        if method.source is not None and issubclass(cls, method.source):
            return method
    raise KeyError(f"No registered method records {cls.__name__}")


def entrypoint_for(source) -> Callable[..., list[Dseqrecord]]:
    """Return the callable that runs the method recording *source*.

    Lets a caller dispatch on recorded provenance without hand-writing a
    source-to-function table::

        products = entrypoint_for(source)(fragments, limit)
    """
    return _ENTRYPOINTS[for_source(source).name]


def available() -> list[str]:
    """Sorted list of registered method names."""
    return sorted(_REGISTRY)


def methods() -> dict[str, Method]:
    """Mapping of name to :class:`Method` for every registered technique."""
    return dict(_REGISTRY)
