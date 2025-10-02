rm -rf _build
sphinx-apidoc --force --implicit-namespaces --module-first -o reference ../src/pydna
sphinx-build -n -W --keep-going -b html ./ ./_build/
