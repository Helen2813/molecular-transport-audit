"""Run the direct-preservation validation suite."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SUITE_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]

UNIT_TEST_RUNNER = (
    ROOT
    / "tools"
    / "run_preservation_unit_tests.py"
)

REGRESSION_RUNNER = (
    ROOT
    / "tools"
    / "run_preservation_regression.py"
)


def run_stage(
    name: str,
    script: Path,
) -> None:
    if not script.exists():
        raise FileNotFoundError(
            "Validation stage is missing: "
            f"{script}"
        )

    print("")
    print("=" * 80)
    print(
        f"Stage: {name}"
    )
    print("=" * 80)

    completed = subprocess.run(
        [
            sys.executable,
            str(script),
        ],
        cwd=ROOT,
        check=False,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "Validation stage failed: "
            f"{name} "
            f"(exit code "
            f"{completed.returncode})"
        )


def main() -> None:
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- direct preservation validation suite"
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

    run_stage(
        "Unit and edge-case tests",
        UNIT_TEST_RUNNER,
    )

    run_stage(
        (
            "Locked GSE239948 "
            "direct-preservation regression"
        ),
        REGRESSION_RUNNER,
    )

    print("")
    print("=" * 80)
    print(
        "Direct preservation validation: PASS"
    )
    print("")
    print(
        "Validated:"
    )
    print(
        "  Canine gene-symbol normalization"
    )
    print(
        "  Duplicate-symbol handling"
    )
    print(
        "  Median imputation"
    )
    print(
        "  Sample-SD standardization"
    )
    print(
        "  Zero-variance filtering"
    )
    print(
        "  Pearson edge construction"
    )
    print(
        "  Upper-triangle extraction"
    )
    print(
        "  Edge Spearman preservation"
    )
    print(
        "  Frozen-score PC1 orientation"
    )
    print(
        "  Frozen-loading concordance"
    )
    print(
        "  Minimum common-gene rule"
    )
    print(
        "  GSE239948 locked direct regression"
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
