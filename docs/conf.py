# -*- coding: utf-8 -*-
# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
from importlib.metadata import version

project = "pydna"
copyright = "2024, Björn F. Johansson"
author = "Björn F. Johansson"

sys.path.insert(0, os.path.abspath("../src/pydna"))
# contents of docs/conf.py

release = version("pydna")
# for example take major/minor
version = ".".join(release.split(".")[:3])

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autodoc.typehints",  # Automatically document type hints in function signatures
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",  # Include links to the source code in the documentation
    "sphinx.ext.doctest",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "numpydoc",
    "sphinx.ext.intersphinx",
    "sphinx_rtd_theme",
    # new:
    # "sphinx_new_tab_link",  # each link opens in a new tab
    "myst_nb",  # Markdown and Jupyter Notebook support
    # "sphinx_copybutton",  # add copy button to code blocks
]


#  https://myst-nb.readthedocs.io/en/latest/computation/execute.html
nb_execution_mode = "off"

myst_enable_extensions = ["dollarmath", "amsmath"]

# Plotly support through require javascript library
# https://myst-nb.readthedocs.io/en/latest/render/interactive.html#plotly
html_js_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js"
]

# https://myst-nb.readthedocs.io/en/latest/configuration.html
# Execution
nb_execution_raise_on_error = True
# Rendering
nb_merge_streams = True


# Add mappings https://kev.inburke.com/kevin/sphinx-interlinks
intersphinx_mapping = {
    "biopython": ("https://biopython.org/docs/latest", None),
    "python": ("http://docs.python.org/3.8", None),
}

# Settings to support markdown files
source_suffix = {
    ".rst": "restructuredtext",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store",
                    "future_release_names.md",
                    "cookbook",
                    "jupyter_execute", # local build folder for notebooks (dev)
                    "pydna_cheat_sheet",
                    "pydna_session",
                    "Makefile",
                    "make",
                    "autogen_docs",
                    "README.md",
                    "example_gallery.md",
                    ]

autodoc_member_order = "bysource"
autodoc_preserve_defaults = True

numpydoc_show_class_members = False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

texinfo_documents = [
    (
        "index",
        "Pydna",
        "Pydna Documentation",
        "Björn Johansson",
        "Pydna",
        "One line description of project.",
        "Miscellaneous",
    ),
]

# Add custom css
html_css_files = [
    "custom.css",
]

html_favicon = "_static/favicon.ico"
