"""Verify reliability outputs against the locked fixture."""

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

MODULES = [
    "M34",
    "M11",
    "M24",
    "M40",
]

EXPECTED_SEEDS = {
    "M34": 44,
    "M11": 144,
    "M24": 244,
    "M40": 344,
}

ABSOLUTE_TOLERANCE = 1.0e-12
RELATIVE_TOLERANCE = 1.0e-10

SUMMARY_METRICS = [
    "split_half_median",
    "split_half_q05",
    "split_half_q95",
    "minimum_gene_loo_correlation",
    "median_gene_loo_correlation",
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
        "--observed-reliability",
        required=True,
    )
    parser.add_argument(
        "--observed-loo",
        required=True,
    )
    parser.add_argument(
        "--expected-reliability",
        required=True,
    )
    parser.add_argument(
        "--expected-loo",
        required=True,
    )
    parser.add_argument(
        "--output",
        required=True,
    )

    return parser.parse_args()


def order_reliability(
    table: pd.DataFrame,
) -> pd.DataFrame:
    result = table.copy()

    if (
        "module_label"
        not in result.columns
    ):
        raise ValueError(
            "Reliability table is missing "
            "module_label."
        )

    result[
        "module_label"
    ] = (
        result[
            "module_label"
        ]
        .astype(str)
    )

    if result[
        "module_label"
    ].duplicated().any():
        raise ValueError(
            "Duplicate module labels."
        )

    result = (
        result.set_index(
            "module_label"
        )
        .reindex(
            MODULES
        )
        .reset_index()
    )

    if result.isna().all(
        axis=1
    ).any():
        raise ValueError(
            "One or more required "
            "modules are missing."
        )

    return result


def normalize_loo(
    table: pd.DataFrame,
) -> pd.DataFrame:
    result = table.copy()

    required = {
        "module_label",
        "left_out_gene",
        "n_genes_remaining",
        "correlation_with_full_score",
    }

    missing = sorted(
        required.difference(
            result.columns
        )
    )

    if missing:
        raise ValueError(
            "LOO table is missing: "
            + ", ".join(missing)
        )

    result[
        "module_label"
    ] = (
        result[
            "module_label"
        ]
        .astype(str)
    )

    result[
        "left_out_gene"
    ] = (
        result[
            "left_out_gene"
        ]
        .astype(str)
    )

    return (
        result.sort_values(
            [
                "module_label",
                "left_out_gene",
            ]
        )
        .reset_index(
            drop=True
        )
    )


def main() -> None:
    args = parse_args()

    observed_reliability_path = Path(
        args.observed_reliability
    )

    observed_loo_path = Path(
        args.observed_loo
    )

    expected_reliability_path = Path(
        args.expected_reliability
    )

    expected_loo_path = Path(
        args.expected_loo
    )

    output_path = Path(
        args.output
    )

    observed_reliability = (
        order_reliability(
            pd.read_csv(
                observed_reliability_path
            )
        )
    )

    expected_reliability = (
        order_reliability(
            pd.read_csv(
                expected_reliability_path
            )
        )
    )

    observed_loo = normalize_loo(
        pd.read_csv(
            observed_loo_path
        )
    )

    expected_loo = normalize_loo(
        pd.read_csv(
            expected_loo_path
        )
    )

    metric_checks: list[
        dict[str, Any]
    ] = []

    summary_match = True
    maximum_summary_difference = 0.0

    for index, module in enumerate(
        MODULES
    ):
        for metric in SUMMARY_METRICS:
            observed_value = float(
                observed_reliability.loc[
                    index,
                    metric,
                ]
            )

            expected_value = float(
                expected_reliability.loc[
                    index,
                    metric,
                ]
            )

            if (
                np.isnan(
                    observed_value
                )
                and np.isnan(
                    expected_value
                )
            ):
                difference = 0.0
            else:
                difference = abs(
                    observed_value
                    - expected_value
                )

            matched = bool(
                np.isclose(
                    observed_value,
                    expected_value,
                    rtol=RELATIVE_TOLERANCE,
                    atol=ABSOLUTE_TOLERANCE,
                    equal_nan=True,
                )
            )

            summary_match = (
                summary_match
                and matched
            )

            maximum_summary_difference = max(
                maximum_summary_difference,
                difference,
            )

            metric_checks.append(
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

    valid_repeats_match = bool(
        np.array_equal(
            pd.to_numeric(
                observed_reliability[
                    "split_half_valid_repeats"
                ],
                errors="coerce",
            ).to_numpy(),
            pd.to_numeric(
                expected_reliability[
                    "split_half_valid_repeats"
                ],
                errors="coerce",
            ).to_numpy(),
        )
    )

    common_genes_match = bool(
        np.array_equal(
            pd.to_numeric(
                observed_reliability[
                    "n_common_genes"
                ],
                errors="coerce",
            ).to_numpy(),
            pd.to_numeric(
                expected_reliability[
                    "n_common_genes"
                ],
                errors="coerce",
            ).to_numpy(),
        )
    )

    seeds_match = True

    for index, module in enumerate(
        MODULES
    ):
        observed_seed = int(
            observed_reliability.loc[
                index,
                "split_half_seed",
            ]
        )

        if (
            observed_seed
            != EXPECTED_SEEDS[module]
        ):
            seeds_match = False

    observed_keys = (
        observed_loo[
            [
                "module_label",
                "left_out_gene",
                "n_genes_remaining",
            ]
        ]
    )

    expected_keys = (
        expected_loo[
            [
                "module_label",
                "left_out_gene",
                "n_genes_remaining",
            ]
        ]
    )

    loo_keys_match = bool(
        observed_keys.equals(
            expected_keys
        )
    )

    if (
        observed_loo.shape[0]
        != expected_loo.shape[0]
    ):
        loo_values_match = False
        maximum_loo_difference = (
            float("inf")
        )
    else:
        observed_values = pd.to_numeric(
            observed_loo[
                "correlation_with_full_score"
            ],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        expected_values = pd.to_numeric(
            expected_loo[
                "correlation_with_full_score"
            ],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        differences = np.abs(
            observed_values
            - expected_values
        )

        maximum_loo_difference = float(
            np.nanmax(
                differences
            )
        )

        loo_values_match = bool(
            np.allclose(
                observed_values,
                expected_values,
                rtol=RELATIVE_TOLERANCE,
                atol=ABSOLUTE_TOLERANCE,
                equal_nan=True,
            )
        )

    passed = bool(
        summary_match
        and valid_repeats_match
        and common_genes_match
        and seeds_match
        and loo_keys_match
        and loo_values_match
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
        "modules":
            MODULES,
        "split_half_repeats":
            2000,
        "checks": {
            "summary_metrics_match":
                summary_match,
            "valid_repeat_counts_match":
                valid_repeats_match,
            "common_gene_counts_match":
                common_genes_match,
            "seed_schedule_match":
                seeds_match,
            "loo_keys_match":
                loo_keys_match,
            "loo_correlations_match":
                loo_values_match,
        },
        "maximum_summary_metric_difference":
            maximum_summary_difference,
        "maximum_loo_correlation_difference":
            maximum_loo_difference,
        "observed_loo_rows":
            int(
                observed_loo.shape[0]
            ),
        "expected_loo_rows":
            int(
                expected_loo.shape[0]
            ),
        "metric_checks":
            metric_checks,
        "artifacts": {
            "observed_reliability_sha256":
                sha256_file(
                    observed_reliability_path
                ),
            "observed_loo_sha256":
                sha256_file(
                    observed_loo_path
                ),
            "expected_reliability_sha256":
                sha256_file(
                    expected_reliability_path
                ),
            "expected_loo_sha256":
                sha256_file(
                    expected_loo_path
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
        "Nextflow preservation "
        "reliability verifier"
    )
    print("=" * 80)

    print(
        "Reliability summary: "
        + (
            "PASS"
            if summary_match
            else "FAIL"
        )
    )

    print(
        "Valid repeat counts: "
        + (
            "PASS"
            if valid_repeats_match
            else "FAIL"
        )
    )

    print(
        "Common gene counts: "
        + (
            "PASS"
            if common_genes_match
            else "FAIL"
        )
    )

    print(
        "Seed schedule: "
        + (
            "PASS"
            if seeds_match
            else "FAIL"
        )
    )

    print(
        "LOO keys: "
        + (
            "PASS"
            if loo_keys_match
            else "FAIL"
        )
    )

    print(
        "LOO correlations: "
        + (
            "PASS"
            if loo_values_match
            else "FAIL"
        )
    )

    print(
        "Maximum summary difference: "
        f"{maximum_summary_difference:.3e}"
    )

    print(
        "Maximum LOO difference: "
        f"{maximum_loo_difference:.3e}"
    )

    print("=" * 80)

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
