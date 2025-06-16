#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import logging
import tempfile
import pathlib
import pytest

pathlib.Path("coverage.xml").unlink(missing_ok=True)
pathlib.Path(".coverage").unlink(missing_ok=True)


def main():
    """docstring."""
    # os.environ["pydna_data_dir"] = tempfile.mkdtemp(prefix="pydna_data_dir_")
    # os.environ["pydna_log_dir"] = tempfile.mkdtemp(prefix="pydna_log_dir_")
    # os.environ["pydna_config_dir"] = tempfile.mkdtemp(prefix="pydna_config_dir_")
    # os.environ["pydna_loglevel"] = str(logging.DEBUG)

    args = [
        "src",  # doctestdir
        "--cov=pydna",
        "--cov-append",
        "--cov-report=html",
        "--cov-report=xml",
        "--capture=no",
        "--durations=10",
        # "--import-mode=importlib",
        "--nbval",
        "--current-env",
        "--doctest-modules",
        "-vvv",
    ]

    return_value_doc_tests = pytest.main(args)

    args = [
        "tests",  # test suite
        "--cov=pydna",
        "--cov-append",
        "--cov-report=html",
        "--cov-report=xml",
        "--capture=no",
        "--durations=10",
        # "--import-mode=importlib",
        "--nbval",
        "--current-env",
        "--doctest-modules",
        "--capture=no",
        "-vvv",
        "--profile",  # profiling
    ]

    pth = pathlib.Path("prof/combined.prof")
    pth.parent.mkdir(parents=True, exist_ok=True)
    pth.write_bytes(b"")
    return_value_unit_tests = pytest.main(args)

    import pstats

    stats = pstats.Stats("./prof/combined.prof")
    stats.sort_stats("cumulative")
    stats.print_stats("pydna/src", 0.1)

    # Or alternatively
    # stats.print_stats("local_path", 20) # Only show 20 of the listings
    # stats.sort_stats('cumulative').print_stats('dir_name', 20) # Sort by cumulative time

    return_value = return_value_doc_tests + return_value_unit_tests

    print("run_test.py return code:", return_value)

    return return_value


if __name__ == "__main__":
    sys.exit(main())
