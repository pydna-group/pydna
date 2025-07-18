# -*- coding: utf-8 -*-
import sys
import os
import logging
import tempfile
import pathlib
import pytest
import pstats

# Remove old coverage files
pathlib.Path("coverage.xml").unlink(missing_ok=True)
pathlib.Path(".coverage").unlink(missing_ok=True)


def main():
    """Run doctests and unit tests in one pytest call with coverage and profiling."""
    # Set up temp dirs and logging level for pydna
    os.environ["pydna_data_dir"] = tempfile.mkdtemp(prefix="pydna_data_dir_")
    os.environ["pydna_log_dir"] = tempfile.mkdtemp(prefix="pydna_log_dir_")
    os.environ["pydna_config_dir"] = tempfile.mkdtemp(prefix="pydna_config_dir_")
    os.environ["pydna_loglevel"] = str(logging.DEBUG)

    # Ensure profile directory exists
    pth = pathlib.Path("prof/combined.prof")
    pth.parent.mkdir(parents=True, exist_ok=True)
    pth.write_bytes(b"")

    # Single pytest call for both src and tests
    args = [
        "src",
        "tests",
        "--cov=pydna",
        "--cov-append",
        "--cov-report=html",
        "--cov-report=xml",
        "--capture=no",
        "--durations=10",
        "--nbval",
        "--current-env",
        "--doctest-modules",
        "--profile",
        "-vvv",
    ]

    return_value = pytest.main(args)

    # Process profiling stats
    stats = pstats.Stats(str(pth))
    stats.sort_stats("cumulative")
    stats.print_stats("pydna/src", 0.1)

    print("run_test.py return code:", return_value)
    return return_value


if __name__ == "__main__":
    sys.exit(main())
