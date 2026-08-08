"""Independent smoke verifier for the unified GSE239948 audit."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SCRIPT_VERSION = "0.1.0"

ATOL = 1.0e-10
RTOL = 1.0e-8

DIRECT_COLUMNS = [
    "edge_spearman",
    "loading_spearman",
    "external_pc1_variance_explained",
    (
        "pc1_orientation_correlation_"
        "with_frozen_score"
    ),
]

RELIABILITY_COLUMNS = [
    "split_half_median",
    "split_half_q05",
    "split_half_q95",
    "minimum_gene_loo_correlation",
    "median_gene_loo_correlation",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--module-summary",
        required=True,
    )
    parser.add_argument(
        "--summary",
        required=True,
    )
    parser.add_argument(
        "--expected-direct",
        required=True,
    )
    parser.add_argument(
        "--expected-coverage",
        required=True,
    )
    parser.add_argument(
        "--expected-reliability",
        required=True,
    )
    parser.add_argument(
        "--config",
        required=True,
    )
    parser.add_argument(
        "--output",
        required=True,
    )

    return parser.parse_args()


def ordered(
    table: pd.DataFrame,
    modules: list[str],
) -> pd.DataFrame:
    result = table.copy()

    result[
        "module_label"
    ] = (
        result[
            "module_label"
        ]
        .astype(str)
    )

    return (
        result.set_index(
            "module_label"
        )
        .reindex(
            modules
        )
        .reset_index()
    )


def compare_numeric(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
    columns: list[str],
) -> tuple[
    bool,
    float,
    list[str],
]:
    mismatches: list[str] = []
    maximum_difference = 0.0

    for column in columns:
        left = pd.to_numeric(
            observed[column],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        right = pd.to_numeric(
            expected[column],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        differences = np.abs(
            left - right
        )

        if np.isfinite(
            differences
        ).any():
            maximum_difference = max(
                maximum_difference,
                float(
                    np.nanmax(
                        differences
                    )
                ),
            )

        if not np.allclose(
            left,
            right,
            rtol=RTOL,
            atol=ATOL,
            equal_nan=True,
        ):
            mismatches.append(
                column
            )

    return (
        len(mismatches) == 0,
        maximum_difference,
        mismatches,
    )


def main() -> None:
    args = parse_args()

    config = json.loads(
        Path(
            args.config
        ).read_text(
            encoding="utf-8"
        )
    )

    modules = [
        str(module)
        for module
        in config[
            "module_order"
        ]
    ]

    observed = ordered(
        pd.read_csv(
            args.module_summary
        ),
        modules,
    )

    expected_direct = ordered(
        pd.read_csv(
            args.expected_direct
        ),
        modules,
    )

    expected_coverage = ordered(
        pd.read_csv(
            args.expected_coverage
        ),
        modules,
    )

    expected_reliability = ordered(
        pd.read_csv(
            args.expected_reliability
        ),
        modules,
    )

    summary = json.loads(
        Path(
            args.summary
        ).read_text(
            encoding="utf-8"
        )
    )

    (
        direct_match,
        direct_difference,
        direct_mismatches,
    ) = compare_numeric(
        observed,
        expected_direct,
        DIRECT_COLUMNS,
    )

    (
        permutation_match,
        permutation_difference,
        permutation_mismatches,
    ) = compare_numeric(
        observed,
        expected_direct,
        [
            "edge_permutation_p",
            "loading_permutation_p",
        ],
    )

    (
        reliability_match,
        reliability_difference,
        reliability_mismatches,
    ) = compare_numeric(
        observed,
        expected_reliability,
        RELIABILITY_COLUMNS,
    )

    coverage_counts_match = bool(
        np.array_equal(
            pd.to_numeric(
                observed[
                    "n_frozen_genes"
                ],
                errors="coerce",
            ).to_numpy(),
            pd.to_numeric(
                expected_coverage[
                    "n_frozen_genes"
                ],
                errors="coerce",
            ).to_numpy(),
        )
        and np.array_equal(
            pd.to_numeric(
                observed[
                    "n_common_genes"
                ],
                errors="coerce",
            ).to_numpy(),
            pd.to_numeric(
                expected_coverage[
                    "n_common_genes"
                ],
                errors="coerce",
            ).to_numpy(),
        )
    )

    coverage_fraction_match = bool(
        np.allclose(
            pd.to_numeric(
                observed[
                    "coverage_fraction"
                ],
                errors="coerce",
            ).to_numpy(
                dtype=float
            ),
            pd.to_numeric(
                expected_coverage[
                    "coverage_fraction"
                ],
                errors="coerce",
            ).to_numpy(
                dtype=float
            ),
            rtol=RTOL,
            atol=ATOL,
            equal_nan=True,
        )
    )

    expected_extreme = (
        config[
            (
                "expected_corrected_"
                "random_extreme_counts"
            )
        ]
    )

    random_p_match = True

    for index, module in enumerate(
        modules
    ):
        expected_p = (
            int(
                expected_extreme[
                    module
                ]
            )
            + 1
        ) / (
            int(
                config[
                    "random_panels"
                ]
            )
            + 1
        )

        observed_p = float(
            observed.loc[
                index,
                "random_panel_empirical_p",
            ]
        )

        if not np.isclose(
            observed_p,
            expected_p,
            rtol=0.0,
            atol=1.0e-15,
        ):
            random_p_match = False

    random_panel_counts_match = bool(
        (
            pd.to_numeric(
                observed[
                    "n_random_panels"
                ],
                errors="coerce",
            )
            == int(
                config[
                    "random_panels"
                ]
            )
        ).all()
    )

    no_fallbacks = bool(
        (
            pd.to_numeric(
                observed[
                    "fallback_draws"
                ],
                errors="coerce",
            )
            == 0
        ).all()
    )

    expected_classes = (
        config[
            "expected_classification"
        ]
    )

    classification_match = True

    for index, module in enumerate(
        modules
    ):
        observed_class = str(
            observed.loc[
                index,
                (
                    "external_canine_"
                    "representation_class"
                ),
            ]
        )

        if (
            observed_class
            != expected_classes[
                module
            ]
        ):
            classification_match = False

    summary_modules = {
        str(item["module"]):
            item
        for item
        in summary[
            "modules"
        ]
    }

    summary_contract_match = bool(
        summary[
            "schema_version"
        ]
        == "0.1.0"
        and summary[
            "run"
        ][
            "status"
        ]
        == "completed"
        and summary[
            "run"
        ][
            "outcome_loaded"
        ]
        is False
        and summary[
            "parameters"
        ][
            "random_panel_matching_mode"
        ]
        == (
            "prestandardization_"
            "variance"
        )
        and summary[
            "diagnostics"
        ][
            "gene_leave_one_out_rows"
        ]
        == 272
        and set(
            summary_modules
        )
        == set(
            modules
        )
    )

    summary_classification_match = (
        summary_contract_match
    )

    if summary_contract_match:
        for module in modules:
            if (
                summary_modules[
                    module
                ][
                    "classification"
                ][
                    "full"
                ]
                != expected_classes[
                    module
                ]
            ):
                summary_classification_match = (
                    False
                )

    maximum_difference = max(
        direct_difference,
        permutation_difference,
        reliability_difference,
    )

    passed = bool(
        direct_match
        and permutation_match
        and reliability_match
        and coverage_counts_match
        and coverage_fraction_match
        and random_p_match
        and random_panel_counts_match
        and no_fallbacks
        and classification_match
        and summary_contract_match
        and summary_classification_match
    )

    payload: dict[
        str,
        Any,
    ] = {
        "script_version":
            SCRIPT_VERSION,
        "generated_at_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),
        "passed":
            passed,
        "checks": {
            "direct_metrics_match":
                direct_match,
            "permutation_p_values_match":
                permutation_match,
            "reliability_metrics_match":
                reliability_match,
            "coverage_counts_match":
                coverage_counts_match,
            "coverage_fraction_match":
                coverage_fraction_match,
            "corrected_random_p_values_match":
                random_p_match,
            "random_panel_counts_match":
                random_panel_counts_match,
            "no_random_panel_fallbacks":
                no_fallbacks,
            "classification_match":
                classification_match,
            "summary_contract_match":
                summary_contract_match,
            "summary_classification_match":
                summary_classification_match,
        },
        "maximum_reference_metric_difference":
            maximum_difference,
        "mismatches": {
            "direct":
                direct_mismatches,
            "permutation":
                permutation_mismatches,
            "reliability":
                reliability_mismatches,
        },
    }

    Path(
        args.output
    ).write_text(
        json.dumps(
            payload,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("=" * 80)
    print(
        "Unified molecular transport "
        "audit verifier"
    )
    print("=" * 80)

    for name, value in payload[
        "checks"
    ].items():
        print(
            f"{name}: "
            + (
                "PASS"
                if value
                else "FAIL"
            )
        )

    print(
        "Maximum reference metric "
        "difference: "
        f"{maximum_difference:.3e}"
    )

    print("=" * 80)

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
