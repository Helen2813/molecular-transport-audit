"""Reproduce locked GSE239948 split-half and gene-LOO reliability."""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


RUNNER_VERSION = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]

CORE_DIR = ROOT / "core"

if str(CORE_DIR) not in sys.path:
    sys.path.insert(
        0,
        str(CORE_DIR),
    )

from transport_audit.preservation_reliability import (
    analyze_reliability,
)


FIXTURE_DIR = (
    ROOT
    / "reference_results"
    / "osteosarcoma_locked"
    / "preservation_fixture"
)

BASE_FIXTURE_MANIFEST = (
    FIXTURE_DIR
    / "fixture_manifest.json"
)

RELIABILITY_FIXTURE_MANIFEST = (
    FIXTURE_DIR
    / "reliability_fixture_manifest.json"
)

EXTERNAL_FILE = (
    FIXTURE_DIR
    / "GSE239948_external_module_expression.csv"
)

WEIGHTS_FILE = (
    FIXTURE_DIR
    / "primary_canine_program_weights.csv"
)

EXPECTED_RELIABILITY_FILE = (
    FIXTURE_DIR
    / "expected_split_half_reliability.csv"
)

EXPECTED_LOO_FILE = (
    FIXTURE_DIR
    / "expected_gene_leave_one_out.csv"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "preservation_reliability_regression"
)

OBSERVED_RELIABILITY_FILE = (
    OUTPUT_DIR
    / "GSE239948_observed_reliability.csv"
)

OBSERVED_LOO_FILE = (
    OUTPUT_DIR
    / "GSE239948_observed_gene_leave_one_out.csv"
)

VERIFICATION_FILE = (
    OUTPUT_DIR
    / "GSE239948_reliability_verification.json"
)

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

RELIABILITY_METRICS = [
    "split_half_median",
    "split_half_q05",
    "split_half_q95",
    "minimum_gene_loo_correlation",
    "median_gene_loo_correlation",
]


def sha256_file(
    path: Path,
) -> str:
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


def load_json(
    path: Path,
) -> dict[str, Any]:
    payload = json.loads(
        path.read_text(
            encoding="utf-8"
        )
    )

    if not isinstance(
        payload,
        dict,
    ):
        raise ValueError(
            f"Expected JSON object: {path}"
        )

    return payload


def verify_fixture_hashes(
    manifest: dict[str, Any],
) -> dict[str, bool]:
    files = manifest.get(
        "fixture_files",
        {},
    )

    checks: dict[
        str,
        bool,
    ] = {}

    for filename, metadata in (
        files.items()
    ):
        path = (
            FIXTURE_DIR
            / filename
        )

        if not path.exists():
            raise FileNotFoundError(
                f"Fixture file missing: {path}"
            )

        expected = str(
            metadata.get(
                "sha256",
                "",
            )
        )

        observed = (
            sha256_file(
                path
            )
        )

        checks[
            filename
        ] = bool(
            expected
            and expected
            == observed
        )

    return checks


def ordered_reliability(
    table: pd.DataFrame,
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
            MODULES
        )
        .reset_index()
    )


def normalize_loo(
    table: pd.DataFrame,
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
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- GSE239948 reliability regression"
    )
    print("=" * 80)

    for path in [
        BASE_FIXTURE_MANIFEST,
        RELIABILITY_FIXTURE_MANIFEST,
        EXTERNAL_FILE,
        WEIGHTS_FILE,
        EXPECTED_RELIABILITY_FILE,
        EXPECTED_LOO_FILE,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Required fixture missing:\n{path}"
            )

    reliability_manifest = (
        load_json(
            RELIABILITY_FIXTURE_MANIFEST
        )
    )

    hash_checks = (
        verify_fixture_hashes(
            reliability_manifest
        )
    )

    if not all(
        hash_checks.values()
    ):
        failed = [
            name
            for name, passed
            in hash_checks.items()
            if not passed
        ]

        raise RuntimeError(
            "Reliability fixture hash "
            "verification failed: "
            + ", ".join(failed)
        )

    contract = reliability_manifest[
        "contract"
    ]

    if int(
        contract[
            "split_half_repeats"
        ]
    ) != 2000:
        raise RuntimeError(
            "Fixture does not specify "
            "2000 split-half repeats."
        )

    if bool(
        reliability_manifest.get(
            "outcome_loaded",
            True,
        )
    ):
        raise RuntimeError(
            "Reliability fixture does not "
            "confirm outcome_loaded=false."
        )

    external = pd.read_csv(
        EXTERNAL_FILE,
        index_col=0,
    )

    weights = pd.read_csv(
        WEIGHTS_FILE
    )

    external.index = (
        external.index.astype(str)
    )

    expected_reliability = (
        ordered_reliability(
            pd.read_csv(
                EXPECTED_RELIABILITY_FILE
            )
        )
    )

    expected_loo = normalize_loo(
        pd.read_csv(
            EXPECTED_LOO_FILE
        )
    )

    (
        observed_reliability,
        observed_loo,
    ) = analyze_reliability(
        external_expression=external,
        weights=weights,
        module_order=MODULES,
        minimum_genes=3,
        n_split_half_repeats=2000,
        base_seed=42,
        module_seed_stride=100,
        split_half_seed_offset=2,
    )

    observed_reliability = (
        ordered_reliability(
            observed_reliability
        )
    )

    observed_loo = (
        normalize_loo(
            observed_loo
        )
    )

    metric_checks: list[
        dict[str, Any]
    ] = []

    metrics_match = True
    maximum_difference = 0.0

    for row_index, module in enumerate(
        MODULES
    ):
        for metric in (
            RELIABILITY_METRICS
        ):
            observed_value = float(
                observed_reliability.loc[
                    row_index,
                    metric,
                ]
            )

            expected_value = float(
                expected_reliability.loc[
                    row_index,
                    metric,
                ]
            )

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

            maximum_difference = max(
                maximum_difference,
                difference,
            )

            metrics_match = (
                metrics_match
                and matched
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

    common_gene_counts_match = bool(
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

    seeds_match = all(
        int(
            observed_reliability.loc[
                index,
                "split_half_seed",
            ]
        )
        == EXPECTED_SEEDS[module]
        for index, module
        in enumerate(MODULES)
    )

    loo_keys_match = bool(
        observed_loo[
            [
                "module_label",
                "left_out_gene",
                "n_genes_remaining",
            ]
        ].equals(
            expected_loo[
                [
                    "module_label",
                    "left_out_gene",
                    "n_genes_remaining",
                ]
            ]
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
        all(hash_checks.values())
        and metrics_match
        and valid_repeats_match
        and common_gene_counts_match
        and seeds_match
        and loo_keys_match
        and loo_values_match
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    observed_reliability.to_csv(
        OBSERVED_RELIABILITY_FILE,
        index=False,
    )

    observed_loo.to_csv(
        OBSERVED_LOO_FILE,
        index=False,
    )

    verification = {
        "runner_version":
            RUNNER_VERSION,
        "generated_at_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),
        "passed":
            passed,
        "outcome_loaded":
            False,
        "fixture_hashes_match":
            bool(
                all(
                    hash_checks.values()
                )
            ),
        "split_half_repeats":
            2000,
        "checks": {
            "summary_metrics_match":
                metrics_match,
            "valid_repeat_counts_match":
                valid_repeats_match,
            "common_gene_counts_match":
                common_gene_counts_match,
            "seed_schedule_match":
                seeds_match,
            "loo_keys_match":
                loo_keys_match,
            "loo_correlations_match":
                loo_values_match,
        },
        "maximum_summary_metric_difference":
            maximum_difference,
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
        "outputs": {
            "reliability": {
                "path":
                    str(
                        OBSERVED_RELIABILITY_FILE
                        .relative_to(ROOT)
                    ),
                "sha256":
                    sha256_file(
                        OBSERVED_RELIABILITY_FILE
                    ),
            },
            "gene_loo": {
                "path":
                    str(
                        OBSERVED_LOO_FILE
                        .relative_to(ROOT)
                    ),
                "sha256":
                    sha256_file(
                        OBSERVED_LOO_FILE
                    ),
            },
        },
    }

    VERIFICATION_FILE.write_text(
        json.dumps(
            verification,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        f"External samples: "
        f"{external.shape[0]}"
    )
    print(
        "Split-half repeats per module: "
        "2000"
    )
    print(
        "Seeds: "
        "M34=44, M11=144, "
        "M24=244, M40=344"
    )
    print(
        "Fixture hashes: PASS"
    )
    print(
        "Reliability summary metrics: "
        + (
            "PASS"
            if metrics_match
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
            if common_gene_counts_match
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
        "Maximum summary metric "
        "difference: "
        f"{maximum_difference:.3e}"
    )
    print(
        "Maximum LOO correlation "
        "difference: "
        f"{maximum_loo_difference:.3e}"
    )

    print("")
    print(
        observed_reliability[
            [
                "module_label",
                "n_common_genes",
                "split_half_seed",
                "split_half_median",
                "split_half_q05",
                "split_half_q95",
                "split_half_valid_repeats",
                "minimum_gene_loo_correlation",
                "median_gene_loo_correlation",
            ]
        ].to_string(
            index=False
        )
    )

    print("")
    print(
        "Verification: "
        + str(
            VERIFICATION_FILE
            .relative_to(ROOT)
        )
    )
    print("=" * 80)

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
