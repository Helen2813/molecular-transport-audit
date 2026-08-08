"""Verify Nextflow preservation outputs against the locked fixture."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SCRIPT_VERSION = "0.1.0"

ABSOLUTE_TOLERANCE = 1.0e-10
RELATIVE_TOLERANCE = 1.0e-8

MODULE_ORDER = [
    "M34",
    "M11",
    "M24",
    "M40",
]

DIRECT_METRICS = [
    "edge_spearman",
    "loading_spearman",
    "external_pc1_variance_explained",
    (
        "pc1_orientation_correlation_"
        "with_frozen_score"
    ),
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--observed-structure",
        required=True,
    )
    parser.add_argument(
        "--observed-coverage",
        required=True,
    )
    parser.add_argument(
        "--expected-structure",
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


def normalize_bool(value: Any) -> bool:
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
        f"Cannot interpret boolean "
        f"value: {value!r}"
    )


def normalize_text(value: Any) -> str:
    if pd.isna(value):
        return ""

    return str(value)


def order_modules(
    table: pd.DataFrame,
) -> pd.DataFrame:
    if "module_label" not in table.columns:
        raise ValueError(
            "Table is missing module_label."
        )

    result = table.copy()

    result["module_label"] = (
        result["module_label"]
        .astype(str)
    )

    if result[
        "module_label"
    ].duplicated().any():
        raise ValueError(
            "Duplicate module labels "
            "were found."
        )

    result = (
        result.set_index(
            "module_label"
        )
        .reindex(
            MODULE_ORDER
        )
    )

    if result.isna().all(
        axis=1
    ).any():
        missing = (
            result.index[
                result.isna().all(
                    axis=1
                )
            ]
            .tolist()
        )

        raise ValueError(
            "Missing modules: "
            + ", ".join(missing)
        )

    return result.reset_index()


def compare_structure(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
) -> tuple[
    bool,
    float,
    list[dict[str, Any]],
    list[str],
]:
    observed = order_modules(
        observed
    )

    expected = order_modules(
        expected
    )

    required = [
        "module_label",
        "n_common_genes",
        *DIRECT_METRICS,
        "estimable",
        "nonestimable_reason",
    ]

    for label, table in [
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
                f"{label} structure table "
                "is missing: "
                + ", ".join(missing)
            )

    mismatches: list[str] = []
    summaries: list[
        dict[str, Any]
    ] = []

    maximum_difference = 0.0

    if not np.array_equal(
        pd.to_numeric(
            observed[
                "n_common_genes"
            ],
            errors="coerce",
        ).to_numpy(),
        pd.to_numeric(
            expected[
                "n_common_genes"
            ],
            errors="coerce",
        ).to_numpy(),
    ):
        mismatches.append(
            "n_common_genes"
        )

    observed_estimable = [
        normalize_bool(value)
        for value
        in observed["estimable"]
    ]

    expected_estimable = [
        normalize_bool(value)
        for value
        in expected["estimable"]
    ]

    if (
        observed_estimable
        != expected_estimable
    ):
        mismatches.append(
            "estimable"
        )

    observed_reasons = [
        normalize_text(value)
        for value
        in observed[
            "nonestimable_reason"
        ]
    ]

    expected_reasons = [
        normalize_text(value)
        for value
        in expected[
            "nonestimable_reason"
        ]
    ]

    if (
        observed_reasons
        != expected_reasons
    ):
        mismatches.append(
            "nonestimable_reason"
        )

    for row_index, module in enumerate(
        MODULE_ORDER
    ):
        for metric in DIRECT_METRICS:
            observed_value = float(
                pd.to_numeric(
                    pd.Series(
                        [
                            observed.loc[
                                row_index,
                                metric,
                            ]
                        ]
                    ),
                    errors="coerce",
                ).iloc[0]
            )

            expected_value = float(
                pd.to_numeric(
                    pd.Series(
                        [
                            expected.loc[
                                row_index,
                                metric,
                            ]
                        ]
                    ),
                    errors="coerce",
                ).iloc[0]
            )

            if (
                np.isfinite(
                    observed_value
                )
                and np.isfinite(
                    expected_value
                )
            ):
                difference = abs(
                    observed_value
                    - expected_value
                )
            elif (
                np.isnan(
                    observed_value
                )
                and np.isnan(
                    expected_value
                )
            ):
                difference = 0.0
            else:
                difference = float(
                    "inf"
                )

            matched = bool(
                np.isclose(
                    observed_value,
                    expected_value,
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

            if not matched:
                mismatches.append(
                    f"{module}:{metric}"
                )

            summaries.append(
                {
                    "module_label":
                        module,
                    "metric":
                        metric,
                    "observed":
                        observed_value,
                    "expected":
                        expected_value,
                    "absolute_difference":
                        difference,
                    "matched":
                        matched,
                }
            )

    return (
        len(mismatches) == 0,
        maximum_difference,
        summaries,
        mismatches,
    )


def compare_coverage(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
) -> tuple[
    bool,
    list[str],
]:
    observed = order_modules(
        observed
    )

    expected = order_modules(
        expected
    )

    integer_columns = [
        "n_frozen_genes",
        "n_common_genes",
    ]

    float_columns = [
        "coverage_fraction",
    ]

    text_columns = [
        "common_genes",
        "missing_external_genes",
    ]

    required = [
        "module_label",
        *integer_columns,
        *float_columns,
        *text_columns,
    ]

    for label, table in [
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
                f"{label} coverage table "
                "is missing: "
                + ", ".join(missing)
            )

    mismatches: list[str] = []

    for column in integer_columns:
        if not np.array_equal(
            pd.to_numeric(
                observed[column],
                errors="coerce",
            ).to_numpy(),
            pd.to_numeric(
                expected[column],
                errors="coerce",
            ).to_numpy(),
        ):
            mismatches.append(column)

    for column in float_columns:
        if not np.allclose(
            pd.to_numeric(
                observed[column],
                errors="coerce",
            ).to_numpy(
                dtype=float
            ),
            pd.to_numeric(
                expected[column],
                errors="coerce",
            ).to_numpy(
                dtype=float
            ),
            rtol=RELATIVE_TOLERANCE,
            atol=ABSOLUTE_TOLERANCE,
            equal_nan=True,
        ):
            mismatches.append(column)

    for column in text_columns:
        observed_values = [
            normalize_text(value)
            for value
            in observed[column]
        ]

        expected_values = [
            normalize_text(value)
            for value
            in expected[column]
        ]

        if (
            observed_values
            != expected_values
        ):
            mismatches.append(column)

    return (
        len(mismatches) == 0,
        mismatches,
    )


def main() -> None:
    args = parse_args()

    observed_structure_path = Path(
        args.observed_structure
    )

    observed_coverage_path = Path(
        args.observed_coverage
    )

    expected_structure_path = Path(
        args.expected_structure
    )

    expected_coverage_path = Path(
        args.expected_coverage
    )

    output_path = Path(args.output)

    observed_structure = pd.read_csv(
        observed_structure_path
    )

    observed_coverage = pd.read_csv(
        observed_coverage_path
    )

    expected_structure = pd.read_csv(
        expected_structure_path
    )

    expected_coverage = pd.read_csv(
        expected_coverage_path
    )

    (
        structure_match,
        maximum_difference,
        metric_summary,
        structure_mismatches,
    ) = compare_structure(
        observed_structure,
        expected_structure,
    )

    (
        coverage_match,
        coverage_mismatches,
    ) = compare_coverage(
        observed_coverage,
        expected_coverage,
    )

    passed = bool(
        structure_match
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
        "structure_match":
            structure_match,
        "coverage_match":
            coverage_match,
        "modules":
            MODULE_ORDER,
        "max_absolute_metric_difference":
            maximum_difference,
        "absolute_tolerance":
            ABSOLUTE_TOLERANCE,
        "relative_tolerance":
            RELATIVE_TOLERANCE,
        "structure_mismatches":
            structure_mismatches,
        "coverage_mismatches":
            coverage_mismatches,
        "metrics":
            metric_summary,
        "artifacts": {
            "observed_structure_sha256":
                sha256_file(
                    observed_structure_path
                ),
            "observed_coverage_sha256":
                sha256_file(
                    observed_coverage_path
                ),
            "expected_structure_sha256":
                sha256_file(
                    expected_structure_path
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
        "Nextflow preservation verifier"
    )
    print("=" * 80)

    print(
        "Direct preservation regression: "
        + (
            "PASS"
            if structure_match
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
        "Maximum absolute metric "
        "difference: "
        f"{maximum_difference:.3e}"
    )

    print("=" * 80)

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
