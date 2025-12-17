#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from pathlib import Path
import subprocess
import pstats


def run_pytest(args):
    """
    Run pytest in a fresh Python subprocess.
    Returns the pytest exit code.
    """
    cmd = [sys.executable, "-m", "pytest", *args]
    print("\n>>>", " ".join(cmd), "\n")
    return subprocess.run(cmd, check=False).returncode


def main() -> int:

    # ------------------------------------------------------------------
    # Anchor execution to project root
    # ------------------------------------------------------------------
    ROOT = Path(__file__).resolve().parent
    os.chdir(ROOT)

    # ------------------------------------------------------------------
    # Clean previous coverage artifacts
    # ------------------------------------------------------------------
    Path("coverage.xml").unlink(missing_ok=True)
    Path(".coverage").unlink(missing_ok=True)

    prof_path = ROOT / Path("prof/combined.prof")
    prof_path.parent.mkdir(parents=True, exist_ok=True)
    prof_path.write_bytes(b"")

    # ------------------------------------------------------------------
    # 1) Doctests (source code)
    # ------------------------------------------------------------------
    doctest_args = [
        "src/pydna",
        "--doctest-modules",
        "--cov=pydna",
        "--cov-report=html",
        "--cov-report=xml",
        "--cov-append",
        "--durations=10",
        "-vvv",
    ]

    rc_doctests = run_pytest(doctest_args)

    # ------------------------------------------------------------------
    # 2) Unit tests (test suite)
    # ------------------------------------------------------------------
    unit_test_args = [
        "tests",
        "--cov=pydna",
        "--cov-report=html",
        "--cov-report=xml",
        "--cov-append",
        "--durations=10",
        "-vvv",
        "--profile",
    ]

    rc_unit_tests = run_pytest(unit_test_args)

    # ------------------------------------------------------------------
    # Optional: profile summary
    # ------------------------------------------------------------------
    if prof_path.exists() and prof_path.stat().st_size > 0:
        print("\n=== Profiling summary (top cumulative) ===\n")
        stats = pstats.Stats(str(prof_path))
        SRC = (Path(ROOT) / "src" / "pydna").resolve()
        stats.stats = {
            k: v
            for k, v in stats.stats.items()
            if Path(k[0]).resolve().is_relative_to(SRC)
        }
        stats.sort_stats("tottime")
        stats.print_stats()

    # ------------------------------------------------------------------
    # Final report
    # ------------------------------------------------------------------
    print(f"\nrun_test.py return code for doctests   : {rc_doctests}")
    print(f"run_test.py return code for unit tests: {rc_unit_tests}")

    return rc_doctests + rc_unit_tests


if __name__ == "__main__":
    sys.exit(main())
