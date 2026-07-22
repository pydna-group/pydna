Guidance for AI coding agents working on pydna.

## Project Overview

pydna is a Python package for representing double-stranded DNA and for
simulating cloning and genetic assembly workflows. Its users include biologists
who may have limited Python experience, so changes should preserve readable
APIs, clear error messages, and scientifically sensible behavior.

Core package code lives under `src/pydna`. Tests live under `tests`, and doctests
may also be collected from `src/pydna`.

## Working Principles

- Read the relevant module, tests, and documentation before editing.
- Keep public APIs stable unless the task explicitly asks for an API change.
- Prefer small, behavior-preserving changes over broad refactors.
- Preserve compatibility with supported Python versions from `pyproject.toml`.
- Be careful with biological semantics: strand orientation, circular vs linear
  molecules, sticky ends, feature coordinates, reverse complements, ambiguous
  nucleotide symbols, and GenBank/FASTA formatting often matter.
- Do not silently discard sequence features, topology, annotations, or metadata.
- When changing user-visible behavior, update or add focused tests and, when
  appropriate, documentation or examples.

## Development Setup

This project uses Poetry.

Common commands:

```bash
uv sync --all-extras --group test --group dev
uv run pytest
uv run pytest tests/path/to/test_file.py
uv run pytest src/pydna/path/to/module.py
uv run ruff check . && uv run ruff format --check .
