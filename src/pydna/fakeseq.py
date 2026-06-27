#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2013-2026 Björn Johansson
# SPDX-FileCopyrightText: 2023-2026 The Project Contributors
# SPDX-License-Identifier: BSD-3-Clause
"""docstring."""


class FakeSeq:
    """docstring."""

    def __init__(
        self,
        length: int,
        n: float = 50e-15,  # 50 fmol = 0.05 pmol
        rf: float = 0.0,
    ) -> None:
        self._length = length
        self.n = n
        self.rf = rf

    def m(self) -> float:
        """Mass of the DNA molecule in grams."""
        # M(Da) * n (mol) = g
        return self.M() * self.n

    def M(self) -> float:
        """M grams/mol."""
        return (308.9 * self._length + 79.0) * 2

    def __len__(self) -> int:
        """docstring."""
        return self._length

    def __lt__(self, other) -> bool:
        """docstring."""
        return self._length < len(other)

    def __repr__(self) -> str:
        """docstring."""
        return f"FakeSeq({self._length:.1e})"

    def __str__(self) -> str:
        """docstring."""
        return self.__repr__()
