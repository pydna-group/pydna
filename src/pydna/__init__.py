#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.


"""
:copyright: Copyright 2013-2023 by Björn Johansson. All rights reserved.
:license:   This code is part of the pydna package, governed by the
            license in LICENSE.txt that should be included as part
            of this package.

pydna
=====
Pydna is a python package providing code for simulation of the creation of
recombinant DNA molecules using
`molecular biology <https://en.wikipedia.org/wiki/Molecular_biology>`_
techniques. Development of pydna happens in this Github `repository <https://github.com/pydna-group/pydna>`_.

Provided:
  1. PCR simulation
  2. Assembly simulation based on shared identical sequences
  3. Primer design for amplification of a given sequence
  4. Automatic design of primer tails for Gibson assembly
     or homologous recombination.
  5. Restriction digestion and cut&paste cloning
  6. Agarose gel simulation
  7. Download sequences from Genbank
  8. Parsing various sequence formats including the capacity to
     handle broken Genbank format

pydna package layout
--------------------

The most important modules and how to import functions or classes from
them are listed below. Class names starts with a capital letter,
functions with a lowercase letter:

::

      from pydna.module import function
      from pydna.module import Class

      Example: from pydna.gel import Gel

      pydna
         ├── amplify
         │         ├── Anneal
         │         └── pcr
         │
         ├── assembly
         │          └── Assembly
         │
         ├── design
         │        ├── assembly_fragments
         │        └── primer_design
         │
         ├── dseqrecord
         │            └── Dseqrecord
         ├── gel
         │     └── Gel
         │
         ├── genbank
         │         ├── genbank
         │         └── Genbank
         │
         ├── parsers
         │         ├── parse
         │         └── parse_primers
         │
         └── readers
                   ├── read
                   └── read_primers



How to use the documentation
----------------------------
Documentation is available as docstrings provided in the source code for
each module.
These docstrings can be inspected by reading the source code directly.
See further below on how to obtain the code for pydna.

In the python shell, use the built-in ``help`` function to view a
function's docstring::

  >>> from pydna import readers
  >>> help(readers.read)
  ... # doctest: +SKIP

The doctrings are also used to provide an automaticly generated reference
manual available online at
`read the docs <https://pydna-group.github.io/pydna>`_.

Docstrings can be explored using `IPython <http://ipython.org/>`_, an
advanced Python shell with
TAB-completion and introspection capabilities. To see which functions
are available in `pydna`,
type `pydna.<TAB>` (where `<TAB>` refers to the TAB key).
Use `pydna.open_config_folder?<ENTER>`to view the docstring or
`pydna.open_config_folder??<ENTER>` to view the source code.

In the `Spyder IDE <https://github.com/spyder-ide/spyder>`_ it is possible
to place the cursor immediately before the name of a module,class or
function and press ctrl+i to bring up docstrings in a separate window in Spyder

Code snippets are indicated by three greater-than signs::

    >>> x=41
    >>> x=x+1
    >>> x
    42

pydna source code
-----------------

The pydna source code is
`available on Github <https://github.com/pydna-group/pydna>`_.

How to get more help
--------------------

Please join the
`Google group <https://groups.google.com/forum/#!forum/pydna>`_
for pydna, this is the preferred location for help. If you find bugs
in pydna itself, open an issue at the
`Github repository <https://github.com/pydna-group/pydna/issues>`_.

Examples of pydna in use
------------------------

See this repository for a collection of
 `examples <https://github.com/MetabolicEngineeringGroupCBMA/pydna-examples?tab=readme-ov-file#pydna-examples>`_.

"""

from pydna._pretty import PrettyTable
from Bio.Restriction import FormattedSeq
import os

__author__ = "Björn Johansson"
__copyright__ = "Copyright 2013 - 2023 Björn Johansson"
__credits__ = ["Björn Johansson", "Mark Budde"]
__license__ = "BSD"
__maintainer__ = "Björn Johansson"
__email__ = "bjorn_johansson@bio.uminho.pt"
__status__ = "Development"  # "Production" #"Prototype"
__version__ = "0.0.0"


class _PydnaWarning(Warning):
    """Pydna warning.

    Pydna uses this warning (or subclasses of it), to make it easy to
    silence all warning messages:

    >>> import warnings
    >>> from pydna import _PydnaWarning
    >>> warnings.simplefilter('ignore', _PydnaWarning)

    Consult the warnings module documentation for more details.
    """

    pass


class _PydnaDeprecationWarning(_PydnaWarning):
    """pydna deprecation warning.

    Pydna uses this warning instead of the built in DeprecationWarning
    since those are ignored by default since Python 2.7.

    To silence all our deprecation warning messages, use:

    >>> import warnings
    >>> from pydna import _PydnaDeprecationWarning
    >>> warnings.simplefilter('ignore', _PydnaDeprecationWarning)

    Code marked as deprecated will be removed in a future version
    of Pydna. This can be discussed in the Pydna google group:
    https://groups.google.com/forum/#!forum/pydna

    """

    pass


def get_env():
    """Print a an ascii table containing all environmental variables.

    Pydna related variables have names that starts with `pydna_`
    """
    _table = PrettyTable(["Variable", "Value"])
    # _table.set_style(_prettytable.DEFAULT)
    _table.align["Variable"] = "l"  # Left align
    _table.align["Value"] = "l"  # Left align
    _table.padding_width = 1  # One space between column edges and contents
    for k, v in sorted(os.environ.items()):
        if k.lower().startswith("pydna"):
            _table.add_row([k, v])
    return _table


def logo():
    """Ascii-art logotype of pydna."""
    from pydna._pretty import pretty_str

    message = f"pydna {__version__}"
    try:
        from pyfiglet import Figlet
    except ModuleNotFoundError:
        pass
    else:
        f = Figlet()
        message = f.renderText(message)
    return pretty_str(message)


## Override Bio.Restriction.FormattedSeq._table


def _make_FormattedSeq_table() -> bytes:
    table = bytearray(256)
    upper_to_lower = ord("A") - ord("a")
    for c in b"ABCDEFGHIJKLMNOPQRSTUVWXYZ":  # Only allow IUPAC letters
        table[c] = c  # map uppercase to uppercase
        table[c - upper_to_lower] = c  # map lowercase to uppercase
    return bytes(table)


FormattedSeq._table = _make_FormattedSeq_table()
