"""Assemble the stable UI-facing molecular transport audit summary."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from transport_audit.random_controls import (
    classify_external_representation,
)


SCRIPT_VERSION = "0.1.0"
SCHEMA_VERSION = "0.1.0"

DIRECT_METRICS = [
    "edge_spearman",
    "loading_spearman",
    "external_pc1_variance_explained",
    (
        "pc1_orientation_correlation_"
        "with_frozen_score"
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--direct",
        required=True,
    )
    parser.add_argument(
        "--coverage",
        required=True,
    )
    parser.add_argument(
        "--inference",
        required=True,
    )
    parser.add_argument(
        "--reliability",
        required=True,
    )
    parser.add_argument(
        "--loo",
        required=True,
    )
    parser.add_argument(
        "--random-controls",
        required=True,
    )
    parser.add_argument(
        "--config",
        required=True,
    )
    parser.add_argument(
        "--module-summary",
        required=True,
    )
    parser.add_argument(
        "--summary",
        required=True,
    )

    return parser.parse_args()


def order_modules(
    table: pd.DataFrame,
    modules: list[str],
) -> pd.DataFrame:
    result = table.copy()

    if "module_label" not in result.columns:
        raise ValueError(
            "Table is missing module_label."
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
            modules
        )
        .reset_index()
    )

    if result.isna().all(
        axis=1
    ).any():
        raise ValueError(
            "One or more modules are missing."
        )

    return result


def assert_direct_consistency(
    direct: pd.DataFrame,
    inference: pd.DataFrame,
) -> None:
    for metric in DIRECT_METRICS:
        left = pd.to_numeric(
            direct[metric],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        right = pd.to_numeric(
            inference[metric],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        if not np.allclose(
            left,
            right,
            rtol=1.0e-8,
            atol=1.0e-10,
            equal_nan=True,
        ):
            raise RuntimeError(
                "Direct and permutation "
                "calculations disagree for "
                f"{metric}."
            )


def json_value(
    value: Any,
) -> Any:
    if value is None:
        return None

    if isinstance(
        value,
        (np.integer,),
    ):
        return int(value)

    if isinstance(
        value,
        (np.floating,),
    ):
        if np.isnan(value):
            return None
        return float(value)

    if isinstance(
        value,
        float,
    ):
        if np.isnan(value):
            return None
        return value

    if pd.isna(value):
        return None

    if isinstance(
        value,
        (np.bool_,),
    ):
        return bool(value)

    return value


def short_classification(
    value: str,
) -> str:
    mapping = {
        (
            "strong_external_canine_"
            "representation_preservation"
        ):
            "strong",
        (
            "partial_external_canine_"
            "representation_preservation"
        ):
            "partial",
        (
            "limited_specific_"
            "external_canine_signal"
        ):
            "limited",
        (
            "no_clear_external_canine_"
            "representation_preservation"
        ):
            "no_clear",
    }

    if value not in mapping:
        raise ValueError(
            "Unknown classification: "
            f"{value}"
        )

    return mapping[value]


def main() -> None:
    args = parse_args()

    config = json.loads(
        Path(
            args.config
        ).read_text(
            encoding="utf-8"
        )
    )

    if bool(
        config.get(
            "outcome_loaded",
            True,
        )
    ):
        raise RuntimeError(
            "Unified audit requires "
            "outcome_loaded=false."
        )

    modules = [
        str(module)
        for module
        in config[
            "module_order"
        ]
    ]

    direct = order_modules(
        pd.read_csv(
            args.direct
        ),
        modules,
    )

    coverage = order_modules(
        pd.read_csv(
            args.coverage
        ),
        modules,
    )

    inference = order_modules(
        pd.read_csv(
            args.inference
        ),
        modules,
    )

    reliability = order_modules(
        pd.read_csv(
            args.reliability
        ),
        modules,
    )

    random_controls = order_modules(
        pd.read_csv(
            args.random_controls
        ),
        modules,
    )

    loo = pd.read_csv(
        args.loo
    )

    assert_direct_consistency(
        direct,
        inference,
    )

    classification_input = (
        inference.merge(
            reliability[
                [
                    "module_label",
                    "split_half_median",
                ]
            ],
            on="module_label",
            how="left",
            validate="one_to_one",
        )
    )

    classified = (
        classify_external_representation(
            classification_input,
            random_controls,
        )
    )

    classification_columns = (
        classified[
            [
                "module_label",
                "edge_q_bh_8",
                "loading_q_bh_8",
                (
                    "external_canine_"
                    "representation_class"
                ),
            ]
        ]
    )

    module_summary = (
        direct[
            [
                "module_label",
                "n_common_genes",
                *DIRECT_METRICS,
            ]
        ]
        .merge(
            coverage[
                [
                    "module_label",
                    "n_frozen_genes",
                    "coverage_fraction",
                ]
            ],
            on="module_label",
            how="left",
            validate="one_to_one",
        )
        .merge(
            inference[
                [
                    "module_label",
                    "edge_permutation_p",
                    "edge_extreme_count",
                    "edge_seed",
                    "loading_permutation_p",
                    "loading_extreme_count",
                    "loading_seed",
                    "n_permutations",
                ]
            ],
            on="module_label",
            how="left",
            validate="one_to_one",
        )
        .merge(
            reliability[
                [
                    "module_label",
                    "split_half_seed",
                    "split_half_median",
                    "split_half_q05",
                    "split_half_q95",
                    "split_half_valid_repeats",
                    (
                        "minimum_gene_"
                        "loo_correlation"
                    ),
                    (
                        "median_gene_"
                        "loo_correlation"
                    ),
                ]
            ],
            on="module_label",
            how="left",
            validate="one_to_one",
        )
        .merge(
            random_controls[
                [
                    "module_label",
                    "matching_mode",
                    "n_random_panels",
                    "random_edge_median",
                    "random_edge_q05",
                    "random_edge_q95",
                    "random_panel_empirical_p",
                    "fallback_draws",
                    "candidate_genes",
                ]
            ],
            on="module_label",
            how="left",
            validate="one_to_one",
        )
        .merge(
            classification_columns,
            on="module_label",
            how="left",
            validate="one_to_one",
        )
    )

    module_summary[
        "classification_short"
    ] = (
        module_summary[
            (
                "external_canine_"
                "representation_class"
            )
        ]
        .astype(str)
        .map(
            short_classification
        )
    )

    module_summary_path = Path(
        args.module_summary
    )

    summary_path = Path(
        args.summary
    )

    module_summary.to_csv(
        module_summary_path,
        index=False,
    )

    module_payload: list[
        dict[str, Any]
    ] = []

    for row in module_summary.to_dict(
        orient="records"
    ):
        module_payload.append(
            {
                "module":
                    row[
                        "module_label"
                    ],
                "coverage": {
                    "n_frozen_genes":
                        json_value(
                            row[
                                "n_frozen_genes"
                            ]
                        ),
                    "n_common_genes":
                        json_value(
                            row[
                                "n_common_genes"
                            ]
                        ),
                    "fraction":
                        json_value(
                            row[
                                "coverage_fraction"
                            ]
                        ),
                },
                "direct_preservation": {
                    "edge_spearman":
                        json_value(
                            row[
                                "edge_spearman"
                            ]
                        ),
                    "loading_spearman":
                        json_value(
                            row[
                                "loading_spearman"
                            ]
                        ),
                    "external_pc1_variance_explained":
                        json_value(
                            row[
                                (
                                    "external_pc1_"
                                    "variance_explained"
                                )
                            ]
                        ),
                    "pc1_orientation_correlation":
                        json_value(
                            row[
                                (
                                    "pc1_orientation_"
                                    "correlation_with_"
                                    "frozen_score"
                                )
                            ]
                        ),
                },
                "inference": {
                    "edge_permutation_p":
                        json_value(
                            row[
                                "edge_permutation_p"
                            ]
                        ),
                    "edge_q_bh_8":
                        json_value(
                            row[
                                "edge_q_bh_8"
                            ]
                        ),
                    "edge_extreme_count":
                        json_value(
                            row[
                                "edge_extreme_count"
                            ]
                        ),
                    "loading_permutation_p":
                        json_value(
                            row[
                                "loading_permutation_p"
                            ]
                        ),
                    "loading_q_bh_8":
                        json_value(
                            row[
                                "loading_q_bh_8"
                            ]
                        ),
                    "loading_extreme_count":
                        json_value(
                            row[
                                "loading_extreme_count"
                            ]
                        ),
                },
                "reliability": {
                    "split_half_median":
                        json_value(
                            row[
                                "split_half_median"
                            ]
                        ),
                    "split_half_q05":
                        json_value(
                            row[
                                "split_half_q05"
                            ]
                        ),
                    "split_half_q95":
                        json_value(
                            row[
                                "split_half_q95"
                            ]
                        ),
                    "valid_repeats":
                        json_value(
                            row[
                                (
                                    "split_half_"
                                    "valid_repeats"
                                )
                            ]
                        ),
                    "minimum_gene_loo_correlation":
                        json_value(
                            row[
                                (
                                    "minimum_gene_"
                                    "loo_correlation"
                                )
                            ]
                        ),
                    "median_gene_loo_correlation":
                        json_value(
                            row[
                                (
                                    "median_gene_"
                                    "loo_correlation"
                                )
                            ]
                        ),
                },
                "random_control": {
                    "matching_mode":
                        row[
                            "matching_mode"
                        ],
                    "n_random_panels":
                        json_value(
                            row[
                                "n_random_panels"
                            ]
                        ),
                    "null_median":
                        json_value(
                            row[
                                "random_edge_median"
                            ]
                        ),
                    "null_q05":
                        json_value(
                            row[
                                "random_edge_q05"
                            ]
                        ),
                    "null_q95":
                        json_value(
                            row[
                                "random_edge_q95"
                            ]
                        ),
                    "empirical_p":
                        json_value(
                            row[
                                (
                                    "random_panel_"
                                    "empirical_p"
                                )
                            ]
                        ),
                    "fallback_draws":
                        json_value(
                            row[
                                "fallback_draws"
                            ]
                        ),
                },
                "classification": {
                    "short":
                        row[
                            "classification_short"
                        ],
                    "full":
                        row[
                            (
                                "external_canine_"
                                "representation_class"
                            )
                        ],
                },
            }
        )

    summary = {
        "schema_version":
            SCHEMA_VERSION,
        "assembler_version":
            SCRIPT_VERSION,
        "generated_at_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),
        "run": {
            "analysis_id":
                config[
                    "analysis_id"
                ],
            "status":
                "completed",
            "reference_cohort":
                config[
                    "reference_cohort"
                ],
            "external_cohort":
                config[
                    "external_cohort"
                ],
            "outcome_loaded":
                False,
        },
        "parameters": {
            "permutations":
                int(
                    config[
                        "permutations"
                    ]
                ),
            "split_half_repeats":
                int(
                    config[
                        "split_half_repeats"
                    ]
                ),
            "random_panels":
                int(
                    config[
                        "random_panels"
                    ]
                ),
            "random_panel_matching_mode":
                config[
                    (
                        "random_panel_"
                        "matching_mode"
                    )
                ],
        },
        "diagnostics": {
            "gene_leave_one_out_rows":
                int(
                    loo.shape[0]
                ),
            "module_count":
                int(
                    module_summary.shape[0]
                ),
        },
        "modules":
            module_payload,
    }

    summary_path.write_text(
        json.dumps(
            summary,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("=" * 80)
    print(
        "Molecular transport audit summary"
    )
    print("=" * 80)

    print(
        module_summary[
            [
                "module_label",
                "edge_spearman",
                "loading_spearman",
                "split_half_median",
                "random_panel_empirical_p",
                "classification_short",
            ]
        ].to_string(
            index=False
        )
    )

    print("")
    print(
        "LOO rows: "
        f"{loo.shape[0]}"
    )
    print(
        "Outcome loaded: False"
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
