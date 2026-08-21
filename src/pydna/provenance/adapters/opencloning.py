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

__all__ = ["source_class_for", "to_source", "step_from_source", "unsupported"]


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


def step_from_source(source) -> Step:
    """Recover the :class:`Step` that produced *source*.

    The method is identified by the ``Source`` class it records, and the
    declared parameters of that method are restored, so the returned step can be
    replayed with :func:`pydna.provenance.replay`.
    """
    from pydna.methods import methods

    for name, method in methods().items():
        if method.source is not None and isinstance(source, method.source):
            inputs = getattr(source, "_get_input_sequences", list)()
            return Step(method=name, inputs=list(inputs), params=dict(method.params))
    raise ValueError(f"No registered method records {type(source).__name__}")
