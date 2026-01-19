#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Classes for nicer Jupyter output.

The pretty_str class is similar to str but has a _repr_pretty_ method
for for nicer string output in the IPython shell and Jupyter notebook.
"""

from prettytable import PrettyTable as Pt
from prettytable import TableStyle
from copy import copy
from typing import List


class pretty_str(str):
    """Thanks to Min RK, UC Berkeley for this."""

    def _repr_pretty_(self, p, cycle):
        p.text(self)


class PrettyTable(Pt):
    """docstring."""

    def lol(self) -> List[list]:
        """docstring."""
        return [self._field_names] + self._rows

    def __repr__(self) -> str:
        """docstring."""
        return self.get_string()

    def _repr_markdown_(self) -> pretty_str:
        c = copy(self)
        c.set_style(TableStyle.MARKDOWN)
        return pretty_str(c.get_string())
