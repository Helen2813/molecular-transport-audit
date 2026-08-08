"""Run reliability unit tests and locked regression."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


SUITE_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]
CORE_DIR = ROOT / "core"

REGRESSION_RUNNER = (
    ROOT
    / "tools"
    / "run_preservation_reliability_regression.py"
)


def run_unit_tests() -> None:
    env = os.environ.copy()

    existing = env.get(
        "PYTHONPATH",
        "",
    )

    env["PYTHONPATH"] = (
        str(CORE_DIR)
        + (
            os.pathsep
            + existing
            if existing
            else ""
        )
    )

    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "unittest",
            "discover",
            "-s",
            str(
                ROOT
                / "tests"
                / "unit"
            ),
            "-p",
            (
                "test_preservation_"
                "reliability.py"
            ),
            "-v",
        ],
        cwd=ROOT,
        env=env,
        check=False,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "Reliability unit tests failed."
        )


def run_regression() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            str(
                REGRESSION_RUNNER
            ),
        ],
        cwd=ROOT,
        check=False,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "Locked reliability regression failed."
        )


def main() -> None:
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- preservation reliability validation suite"
    )
    print("=" * 80)
    print(
        f"Suite version: "
        f"{SUITE_VERSION}"
    )
    print(
        f"Python: "
        f"{sys.version.split()[0]}"
    )

    print("")
    print("=" * 80)
    print(
        "Stage: Reliability unit "
        "and RNG-contract tests"
    )
    print("=" * 80)

    run_unit_tests()

    print("")
    print(
        "Reliability unit tests: PASS"
    )

    print("")
    print("=" * 80)
    print(
        "Stage: Locked GSE239948 "
        "split-half and gene-LOO regression"
    )
    print("=" * 80)

    run_regression()

    print("")
    print("=" * 80)
    print(
        "Preservation reliability "
        "validation: PASS"
    )
    print("")
    print(
        "Validated:"
    )
    print(
        "  Frozen split-half seed schedule"
    )
    print(
        "  2000 non-overlapping gene splits"
    )
    print(
        "  Frozen signed half-scores"
    )
    print(
        "  Pearson half-score agreement"
    )
    print(
        "  Median / q05 / q95 summaries"
    )
    print(
        "  Valid-repeat accounting"
    )
    print(
        "  Per-gene leave-one-out stability"
    )
    print(
        "  Minimum and median LOO summaries"
    )
    print(
        "  No outcome data"
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
