# Contributing to Pydna

Thank you for considering a contribution to Pydna. Contributions can include
bug reports, documentation improvements, examples, tests, and code changes.

## Getting Started

1. Fork or clone the repository.
2. Create a branch for your change.
3. Install the project in a local development environment.
4. Run the test suite before opening a pull request.

## Development Setup

This project uses Python and Poetry metadata. A typical local setup is:

```bash
poetry install
poetry run pytest
```

If you do not use Poetry, install the package with the development dependencies
supported by your environment and run `pytest`.

## Making Changes

- Keep changes focused on one bug, feature, or documentation improvement.
- Add or update tests when behavior changes.
- Preserve the existing style of nearby code.
- Keep public API changes deliberate and documented.
- Use clear commit messages that explain the reason for the change.

## Headers and Licensing

Python source files should keep the SPDX license header used in the repository.
New project source files should use the BSD-3-Clause identifier unless there is
a documented reason to do otherwise.

## Pull Requests

Before opening a pull request:

1. Run the relevant tests.
2. Update documentation if user-facing behavior changes.
3. Describe what changed and why.
4. Mention any compatibility concerns or follow-up work.

Review may ask for tests, smaller changes, clearer documentation, or a different
implementation approach. That process helps keep the project maintainable.

## Reporting Issues

When reporting a bug, include:

- The Pydna version or commit.
- The Python version.
- A small reproducible example.
- The expected behavior and the actual behavior.
- Any relevant traceback or output.

Feature requests are welcome. Please describe the use case and why the feature
belongs in Pydna.
