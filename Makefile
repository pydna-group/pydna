.PHONY: setup dev build test test-all lint fmt clean

setup:
	uv sync --all-extras --group test --group dev

dev:
	uv run python -c "import pydna; print(f'pydna {pydna.__version__} ready')"

build:
	uv build

test:
	uv run pytest tests/ -q

test-all:
	uv run python run_test.py
	uv run ruff check src/pydna
	uv run ruff format --check src/pydna tests

lint:
	uv run ruff check src/pydna

fmt:
	uv run ruff format src/pydna tests

clean:
	rm -rf dist/ .pytest_cache/ .mypy_cache/ *.egg-info coverage.xml htmlcov/ prof/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
