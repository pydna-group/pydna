# -*- coding: utf-8 -*-
"""Mapping between pydna :class:`~pydna.provenance.step.Step` records and
OpenCloning ``Source`` objects.

pydna records history in its own schema; OpenCloning support is expressed here
as a translation. Keeping it in one place means a technique with no OpenCloning
counterpart is still fully recorded by pydna, and the gap is reported honestly
rather than silently dropping the step's history.

The mapping is derived from the method registry — each
:class:`~pydna.methods._engine.Method` already declares the ``source`` it
records — so there is no second table to keep in sync.
"""

from __future__ import annotations

from pydna.provenance.step import Step

__all__ = ["source_class_for", "method_name_for", "to_source", "unsupported"]


def source_class_for(method_name: str):
    """Return the OpenCloning ``Source`` class a method records, or None.

    >>> from pydna.provenance.adapters.opencloning import source_class_for
    >>> source_class_for("gibson").__name__
    'GibsonAssemblySource'
    >>> source_class_for("not_a_method") is None
    True
    """
    from pydna.methods import methods

    method = methods().get(method_name)
    return None if method is None else method.source


def unsupported() -> list[str]:
    """Names of registered methods with no OpenCloning ``Source``.

    A method appears here when pydna can record it but OpenCloning cannot
    represent it. Exporting such a step must fail loudly — see :func:`to_source`.

    >>> from pydna.provenance.adapters.opencloning import unsupported
    >>> "gibson" in unsupported()
    False
    """
    from pydna.methods import methods

    return sorted(name for name, m in methods().items() if m.source is None)


def to_source(step: Step, **extra_fields):
    """Build the OpenCloning ``Source`` describing *step*.

    Raises
    ------
    ValueError
        If the step's method has no OpenCloning counterpart. Failing here is
        deliberate: silently returning nothing would lose provenance.
    """
    source_cls = source_class_for(step.method)
    if source_cls is None:
        raise ValueError(
            f"Method {step.method!r} has no OpenCloning source; its history "
            "cannot be exported to a CloningStrategy."
        )
    return source_cls(**extra_fields)


def method_name_for(source) -> str:
    """Name of the method that records *source* as its provenance.

    This is the half of the mapping that can be derived: the registry already
    knows which technique records each ``Source`` class. The other half — the
    parameters the technique was called with — is not in the class, so the
    ``Source`` itself supplies it (see ``AssemblySource._replay_params``).

    >>> from pydna.opencloning_models import GibsonAssemblySource
    >>> from pydna.provenance.adapters.opencloning import method_name_for
    >>> method_name_for(GibsonAssemblySource)
    'gibson'

    Raises
    ------
    NotImplementedError
        If no registered method records that ``Source`` class, so pydna cannot
        replay it.
    """
    from pydna.methods import methods

    cls = source if isinstance(source, type) else type(source)
    registered = methods().values()
    # An exact match first, so a specialised source resolves to its own method
    # rather than the one recording its parent — CRISPRSource subclasses
    # HomologousRecombinationSource but is written by crispr_integration.
    for method in registered:
        if method.source is cls:
            return method.name
    for method in registered:
        if method.source is not None and issubclass(cls, method.source):
            return method.name
    raise NotImplementedError(f"_replay_products() not implemented for {cls.__name__}")
