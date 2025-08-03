# Documentation

Documentation is built using [Sphinx](http://www.sphinx-doc.org/) from [docstrings](https://www.python.org/dev/peps/pep-0257/)
using a GitHub [action](https://github.com/pydna-group/pydna/actions/workflows/publish-docs.yml).
The [numpy](www.numpy.org) [docstring format](https://numpy.org/doc/stable/dev/howto-docs.html#docstring-intro) is used.

Below the commands to run a local sphinx server that auto-updated when files are changed.

```bash
# Install docs dependency group
poetry install --with docs

# Generate the example gallery and the API reference
bash autogen_docs.sh

# Start the sphinx server to see docs live by default at http://127.0.0.1:8000/
sphinx-autobuild --watch src/ docs docs/_build/html
```

When you run the sphinx server live, it not always updates things if you make relevant changes (e.g. adding new pages,
changing the css, etc.). To force an update, you can run the command below which deletes the existing build folder and
then rebuilds it.

```bash
rm -rf docs/_build/html && sphinx-autobuild --watch src/ docs docs/_build/html
```

## Adding new sections to the documentation

You can add new sections (equivalent to "Getting started" or "Example gallery") by creating a new `.rst` or `.md` file in the `docs` folder, and then adding a reference to it in the `.. toctree::` directive in the `docs/index.rst` file.

## Auto-generated files

The script `autogen_docs.sh` is run in the github action before creating the docs. If you want to reproduce locally, you
should run it from the `docs` folder

```bash
cd docs
bash autogen_docs.sh
```


* It converts all notebooks in the `docs/notebooks` folder to `.md` in `docs/markdown_notebooks` (excluded from git)
* It creates all files in `docs/modules`, which are used to generate the API reference. For instance, it will create
  a `docs/modules/index.rst` file that starts like this:

  ```rst
  Modules
  =======

  .. toctree::
     :maxdepth: 1
     :caption: Modules:

     pydna
     pydna__pretty
  ```
  And then individual files for each module and submodule, e.g. `docs/modules/pydna.rst`.

  ```
  pydna
  =====

  .. automodule:: pydna
      :members:
      :undoc-members:
      :show-inheritance:
  ```

## Text imported from README.md

To avoid having to maintain the same text in multiple files, fragments of the `README.md` are imported using the directive
`include`. For instance, in the `installation.rst` file, you can find the code below. What this does is to import the text of the README.md file between the start and end markers, which are markdown comments and therefore not rendered.

```rst
.. include:: ../README.md
   :parser: myst_parser.sphinx_
   :start-after: <!-- docs/installation.rst-start -->
   :end-before: <!-- docs/installation.rst-end -->
```

## Including notebooks in the getting started and example sections

You can see the example of how to do this in the `getting_started.md` file. Note that the notebooks present in the `docs/notebooks` folder will automatically be converted to markdown in the `docs/markdown_notebooks` folder. So if you have a notebook `docs/notebooks/Example_Gibson.ipynb`, it will be converted to `docs/markdown_notebooks/Example_Gibson.md` and you can use that file path to make a link to it.

## Custom CSS

For now, I have used css to make notebook outputs that are too long scrollable, and to add a small label `python code` to the code cells and `output` to the output cells.

For further customization, you can edit the `custom.css` file.

## Misc

Other changes, such as changing the favicon, the css etc., can be made in the `conf.py` file. See the [sphinx docs](https://www.sphinx-doc.org/en/master/usage/configuration.html) and the [sphinx-rtd-theme](https://sphinx-rtd-theme.readthedocs.io/en/stable/configuring.html) docs for more information.
