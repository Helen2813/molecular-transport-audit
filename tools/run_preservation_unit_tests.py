"""Run direct-preservation unit tests."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


RUNNER_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]
CORE_DIR = ROOT / "core"
TEST_DIR = ROOT / "tests" / "unit"


def main() -> None:
    if str(CORE_DIR) not in sys.path:
        sys.path.insert(
            0,
            str(CORE_DIR),
        )

    if not TEST_DIR.exists():
        raise FileNotFoundError(
            "Unit-test directory not found: "
            f"{TEST_DIR}"
        )

    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- preservation unit tests"
    )
    print("=" * 80)
    print(
        f"Runner version: "
        f"{RUNNER_VERSION}"
    )
    print(
        f"Test directory: "
        f"{TEST_DIR}"
    )
    print("")

    suite = (
        unittest.defaultTestLoader
        .discover(
            start_dir=str(
                TEST_DIR
            ),
            pattern=(
                "test_preservation.py"
            ),
        )
    )

    result = (
        unittest.TextTestRunner(
            verbosity=2
        )
        .run(
            suite
        )
    )

    print("")
    print("=" * 80)

    if result.wasSuccessful():
        print(
            "Preservation unit tests: PASS"
        )
        print("=" * 80)
        return

    print(
        "Preservation unit tests: FAIL"
    )
    print("=" * 80)

    raise SystemExit(1)


if __name__ == "__main__":
    main()
