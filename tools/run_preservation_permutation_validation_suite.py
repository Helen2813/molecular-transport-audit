"""Run preservation permutation unit tests and locked regression."""

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
    / "run_preservation_permutation_regression.py"
)


def run_unit_tests() -> None:
    env = os.environ.copy()

    existing_pythonpath = env.get(
        "PYTHONPATH",
        "",
    )

    env["PYTHONPATH"] = (
        str(CORE_DIR)
        + (
            os.pathsep
            + existing_pythonpath
            if existing_pythonpath
            else ""
        )
    )

    command = [
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
        "test_preservation_inference.py",
        "-v",
    ]

    completed = subprocess.run(
        command,
        cwd=ROOT,
        env=env,
        check=False,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "Preservation permutation "
            "unit tests failed."
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
            "Locked preservation "
            "permutation regression failed."
        )


def main() -> None:
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- preservation permutation validation suite"
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
        "Stage: Unit and RNG-contract tests"
    )
    print("=" * 80)

    run_unit_tests()

    print("")
    print(
        "Permutation unit tests: PASS"
    )

    print("")
    print("=" * 80)
    print(
        "Stage: Locked GSE239948 "
        "gene-label permutation regression"
    )
    print("=" * 80)

    run_regression()

    print("")
    print("=" * 80)
    print(
        "Preservation permutation "
        "validation: PASS"
    )
    print("")

    print(
        "Validated:"
    )
    print(
        "  Frozen module seed schedule"
    )
    print(
        "  Edge gene-label permutations"
    )
    print(
        "  Loading-label permutations"
    )
    print(
        "  Two-sided absolute-statistic null"
    )
    print(
        "  Monte-Carlo +1 correction"
    )
    print(
        "  Exact legacy extreme counts"
    )
    print(
        "  BH family across 8 direct tests"
    )
    print(
        "  No outcome data"
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
