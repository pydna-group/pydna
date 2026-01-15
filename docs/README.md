# Documentation

Documentation is built using [Sphinx](http://www.sphinx-doc.org/) from [docstrings](https://www.python.org/dev/peps/pep-0257/)
using a GitHub [action](https://github.com/pydna-group/pydna/actions/workflows/publish-docs.yml).
The [numpy](www.numpy.org) [docstring format](https://numpy.org/doc/stable/dev/howto-docs.html#docstring-intro) is used.

Below the commands to build the documentation locally.

```bash
# Install docs dependency group
poetry install --with docs
cd docs
bash build_docs.sh
```

Then, serve the file `docs/_build/index.html` to see the documentation locally. You can use the `python -m http.server -d _build` command to serve the file at `http://127.0.0.1:8000/`.

## Adding new sections to the documentation

You can add new sections (equivalent to "Getting started" or "Example gallery") by creating a new `.rst` or `.md` file in the `docs` folder, and then adding a reference to it in the `.. toctree::` directive in the `docs/index.rst` file.

To add sub-sections, create a toctree in child files. For instance, see [getting_started.rst](getting_started.rst) for how to add sub-sections to the getting started section.

## Text imported from README.md

To avoid having to maintain the same text in multiple files, fragments of the `README.md` are imported using the directive
`include`. For instance, in the `installation.rst` file, you can find the code below. What this does is to import the text of the README.md file between the start and end markers, which are markdown comments and therefore not rendered.

```rst
.. include:: ../README.md
   :parser: myst_parser.sphinx_
   :start-after: <!-- docs/installation.rst-start -->
   :end-before: <!-- docs/installation.rst-end -->
```

## Including notebooks in the documentation

All notebooks in the `docs/notebooks` folder will automatically be converted to markdown in the `docs/markdown_notebooks` folder. So if you have a notebook `docs/notebooks/Example_Gibson.ipynb`, it will be converted to `docs/markdown_notebooks/Example_Gibson.md` and you can use that file path to make a link to it.

You can see the example of how to do this in the `getting_started.rst` file. Note that when making a link, you should use the path with an `.html` extension, not `.ipynb` or `.md`.

## Custom CSS

For now, I have used css to make notebook outputs that are too long scrollable, and to add a small label `python code` to the code cells and `output` to the output cells.

For further customization, you can edit the `custom.css` file.

## Misc

Other changes, such as changing the favicon, the css etc., can be made in the `conf.py` file. See the [sphinx docs](https://www.sphinx-doc.org/en/master/usage/configuration.html) and the [sphinx-rtd-theme](https://sphinx-rtd-theme.readthedocs.io/en/stable/configuring.html) docs for more information.

## Build docs using Sphinx command line tools

See [build_docs.sh](build_docs.sh) for the commands to build the documentation.

## Debugging notebook execution

If notebook execution fails, it might be that jupyter is using another Kernel than the one you want to use. In my case, I had to delete existing kernels which mapped to other python environments. You can find them:

```
poetry run jupyter kernelspec list
```

Then, go to those folders and open the `kernel.json` file. You will see which env they use. I deleted all kernel folders and all envs, since I favour creating envs in the project directory.
