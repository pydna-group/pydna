# -*- coding: utf-8 -*-
"""pydna's own record of a cloning step.

A :class:`Step` states plainly what happened: which method ran, on which
sequences, with which parameters. That is enough to replay the step, so pydna
does not have to guess parameters back from the geometry of the product — which
is what the ``Source`` classes in :mod:`pydna.opencloning_models` must do today
(``allow_blunt`` inferred from an overlap of zero, integration-vs-excision from
the number of inputs, and so on).

Because the schema is pydna's own, a new technique can record its provenance
immediately. Interoperability with OpenCloning is then a *mapping*, in
:mod:`pydna.provenance.adapters.opencloning`, rather than a precondition.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

__all__ = ["Step", "replay"]


@dataclass(frozen=True)
class Step:
    """One recorded cloning operation.

    Parameters
    ----------
    method : str
        Name of the method in the :mod:`pydna.methods` registry, e.g.
        ``"gibson"``.
    inputs : list
        The sequences the method consumed, in the order it received them.
    params : dict
        The parameters the method was called with. Stored explicitly so replay
        is exact rather than reconstructed.
    """

    method: str
    inputs: list = field(default_factory=list)
    params: dict[str, Any] = field(default_factory=dict)

    def describe(self) -> str:
        """One-line human-readable summary of the step.

        >>> from pydna.provenance import Step
        >>> Step("gibson", [], {"limit": 25}).describe()
        'gibson(limit=25)'
        """
        args = ", ".join(f"{k}={v!r}" for k, v in sorted(self.params.items()))
        return f"{self.method}({args})"


def replay(step: Step) -> list:
    """Re-run *step* and return its products.

    This is the whole replay implementation: because the parameters were
    recorded, every method replays through the same two lines, instead of one
    hand-written ``_replay_products`` per provenance class.
    """
    from pydna.methods import get, run

    return run(get(step.method), step.inputs, **step.params)
