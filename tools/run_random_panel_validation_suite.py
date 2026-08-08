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
    / "run_random_panel_validation.py"
)


def run_unit_tests() -> None:
    env = os.environ.copy()

    existing = env.get(
        "PYTHONPATH",
        "",
    )

    env[
        "PYTHONPATH"
    ] = (
        str(
            CORE_DIR
        )
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
            "test_random_controls.py",
            "-v",
        ],
        cwd=ROOT,
        env=env,
        check=False,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "Random-control unit tests failed."
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
            "Random-panel real-data "
            "validation failed."
        )


def main() -> None:
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- random-panel validation suite"
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
        "Stage: random-panel unit tests"
    )
    print("=" * 80)

    run_unit_tests()

    print("")
    print(
        "Random-panel unit tests: PASS"
    )

    print("")
    print("=" * 80)
    print(
        "Stage: legacy exact + "
        "corrected GSE239948 controls"
    )
    print("=" * 80)

    run_regression()

    print("")
    print("=" * 80)
    print(
        "Random-panel validation suite: PASS"
    )
    print("")
    print(
        "Validated:"
    )
    print(
        "  Frozen 1000-panel schedule"
    )
    print(
        "  Frozen candidate-gene exclusion"
    )
    print(
        "  Without-replacement sampling"
    )
    print(
        "  Same-bin preference"
    )
    print(
        "  Legacy fallback rule"
    )
    print(
        "  Shared RNG stream"
    )
    print(
        "  Edge-preservation null statistic"
    )
    print(
        "  Two-sided +1 empirical p"
    )
    print(
        "  Legacy classification regression"
    )
    print(
        "  Corrected pre-standardization "
        "variability sensitivity"
    )
    print(
        "  No outcome data"
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
