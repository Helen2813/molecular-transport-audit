"""Verify Nextflow scoring outputs against the locked fixture."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_VERSION = "0.1.0"

ABSOLUTE_TOLERANCE = 1.0e-10
RELATIVE_TOLERANCE = 1.0e-8


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--observed-scores",
        required=True,
    )
    parser.add_argument(
        "--observed-coverage",
        required=True,
    )
    parser.add_argument(
        "--expected-scores",
        required=True,
    )
    parser.add_argument(
        "--expected-coverage",
        required=True,
    )
    parser.add_argument(
        "--output",
        required=True,
    )

    return parser.parse_args()


def normalize_bool(
    value: object,
) -> bool:
    if isinstance(
        value,
        (bool, np.bool_),
    ):
        return bool(value)

    text = (
        str(value)
        .strip()
        .lower()
    )

    if text in {
        "true",
        "1",
        "yes",
    }:
        return True

    if text in {
        "false",
        "0",
        "no",
    }:
        return False

    raise ValueError(
        f"Cannot interpret boolean: "
        f"{value!r}"
    )


def normalize_text(
    value: object,
) -> str:
    if pd.isna(value):
        return ""

    return str(value)


def compare_scores(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
) -> tuple[
    bool,
    float,
    list[dict[str, object]],
]:
    if not observed.index.equals(
        expected.index
    ):
        return (
            False,
            float("inf"),
            [],
        )

    if list(
        observed.columns
    ) != list(
        expected.columns
    ):
        return (
            False,
            float("inf"),
            [],
        )

    results: list[
        dict[str, object]
    ] = []

    overall_match = True
    maximum_difference = 0.0

    for column in expected.columns:
        observed_values = (
            pd.to_numeric(
                observed[column],
                errors="coerce",
            )
            .to_numpy(
                dtype=float
            )
        )

        expected_values = (
            pd.to_numeric(
                expected[column],
                errors="coerce",
            )
            .to_numpy(
                dtype=float
            )
        )

        finite = (
            np.isfinite(
                observed_values
            )
            & np.isfinite(
                expected_values
            )
        )

        if np.any(finite):
            difference = float(
                np.max(
                    np.abs(
                        observed_values[
                            finite
                        ]
                        - expected_values[
                            finite
                        ]
                    )
                )
            )
        else:
            difference = 0.0

        matched = bool(
            np.allclose(
                observed_values,
                expected_values,
                rtol=(
                    RELATIVE_TOLERANCE
                ),
                atol=(
                    ABSOLUTE_TOLERANCE
                ),
                equal_nan=True,
            )
        )

        maximum_difference = max(
            maximum_difference,
            difference,
        )

        overall_match = (
            overall_match
            and matched
        )

        results.append(
            {
                "score_column":
                    column,
                "matched":
                    matched,
                "max_absolute_difference":
                    difference,
            }
        )

    return (
        overall_match,
        maximum_difference,
        results,
    )


def compare_coverage(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
) -> tuple[
    bool,
    list[str],
]:
    key_columns = [
        "module_label",
        "mapping",
    ]

    compare_columns = [
        "n_frozen_genes",
        "n_available_genes",
        "coverage_fraction",
        "minimum_rule_passed",
        "available_genes",
        "missing_genes",
    ]

    required = (
        key_columns
        + compare_columns
    )

    for name, table in [
        ("observed", observed),
        ("expected", expected),
    ]:
        missing = [
            column
            for column in required
            if column
            not in table.columns
        ]

        if missing:
            raise ValueError(
                f"{name} coverage table "
                "is missing columns: "
                + ", ".join(
                    missing
                )
            )

    left = (
        expected[
            required
        ]
        .sort_values(
            key_columns
        )
        .reset_index(
            drop=True
        )
    )

    right = (
        observed[
            required
        ]
        .sort_values(
            key_columns
        )
        .reset_index(
            drop=True
        )
    )

    mismatches: list[str] = []

    if (
        left[
            key_columns
        ].to_dict(
            "records"
        )
        != right[
            key_columns
        ].to_dict(
            "records"
        )
    ):
        return (
            False,
            ["coverage_key_family"],
        )

    for column in [
        "n_frozen_genes",
        "n_available_genes",
    ]:
        if not np.array_equal(
            pd.to_numeric(
                left[column],
                errors="coerce",
            ).to_numpy(),
            pd.to_numeric(
                right[column],
                errors="coerce",
            ).to_numpy(),
        ):
            mismatches.append(
                column
            )

    if not np.allclose(
        pd.to_numeric(
            left[
                "coverage_fraction"
            ],
            errors="coerce",
        ).to_numpy(
            dtype=float
        ),
        pd.to_numeric(
            right[
                "coverage_fraction"
            ],
            errors="coerce",
        ).to_numpy(
            dtype=float
        ),
        rtol=RELATIVE_TOLERANCE,
        atol=ABSOLUTE_TOLERANCE,
        equal_nan=True,
    ):
        mismatches.append(
            "coverage_fraction"
        )

    left_boolean = [
        normalize_bool(value)
        for value
        in left[
            "minimum_rule_passed"
        ]
    ]

    right_boolean = [
        normalize_bool(value)
        for value
        in right[
            "minimum_rule_passed"
        ]
    ]

    if left_boolean != right_boolean:
        mismatches.append(
            "minimum_rule_passed"
        )

    for column in [
        "available_genes",
        "missing_genes",
    ]:
        left_values = [
            normalize_text(value)
            for value
            in left[column]
        ]

        right_values = [
            normalize_text(value)
            for value
            in right[column]
        ]

        if left_values != right_values:
            mismatches.append(
                column
            )

    return (
        len(mismatches) == 0,
        mismatches,
    )


def main() -> None:
    args = parse_args()

    observed_scores_path = Path(
        args.observed_scores
    )
    observed_coverage_path = Path(
        args.observed_coverage
    )
    expected_scores_path = Path(
        args.expected_scores
    )
    expected_coverage_path = Path(
        args.expected_coverage
    )
    output_path = Path(
        args.output
    )

    observed_scores = pd.read_csv(
        observed_scores_path,
        index_col=0,
    )

    expected_scores = pd.read_csv(
        expected_scores_path,
        index_col=0,
    )

    observed_coverage = pd.read_csv(
        observed_coverage_path
    )

    expected_coverage = pd.read_csv(
        expected_coverage_path
    )

    observed_scores.index = (
        observed_scores.index
        .astype(str)
    )

    expected_scores.index = (
        expected_scores.index
        .astype(str)
    )

    (
        scores_match,
        maximum_difference,
        score_columns,
    ) = compare_scores(
        observed_scores,
        expected_scores,
    )

    (
        coverage_match,
        coverage_mismatches,
    ) = compare_coverage(
        observed_coverage,
        expected_coverage,
    )

    passed = bool(
        scores_match
        and coverage_match
    )

    payload = {
        "script_version":
            SCRIPT_VERSION,
        "created_at_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),
        "passed":
            passed,
        "scores_match":
            scores_match,
        "coverage_match":
            coverage_match,
        "n_samples":
            int(
                expected_scores.shape[0]
            ),
        "n_score_columns":
            int(
                expected_scores.shape[1]
            ),
        "max_absolute_score_difference":
            maximum_difference,
        "absolute_tolerance":
            ABSOLUTE_TOLERANCE,
        "relative_tolerance":
            RELATIVE_TOLERANCE,
        "coverage_mismatches":
            coverage_mismatches,
        "score_columns":
            score_columns,
        "artifacts": {
            "observed_scores_sha256":
                sha256_file(
                    observed_scores_path
                ),
            "observed_coverage_sha256":
                sha256_file(
                    observed_coverage_path
                ),
            "expected_scores_sha256":
                sha256_file(
                    expected_scores_path
                ),
            "expected_coverage_sha256":
                sha256_file(
                    expected_coverage_path
                ),
        },
    }

    output_path.write_text(
        json.dumps(
            payload,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("=" * 80)
    print(
        "Nextflow scoring regression verifier"
    )
    print("=" * 80)
    print(
        "Score regression: "
        + (
            "PASS"
            if scores_match
            else "FAIL"
        )
    )
    print(
        "Coverage regression: "
        + (
            "PASS"
            if coverage_match
            else "FAIL"
        )
    )
    print(
        "Maximum absolute difference: "
        f"{maximum_difference:.3e}"
    )
    print("=" * 80)

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
