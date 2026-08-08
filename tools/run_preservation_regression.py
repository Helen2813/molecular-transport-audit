"""Run the locked GSE239948 direct-preservation regression."""

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

from transport_audit.preservation import (
    analyze_preservation,
)


FIXTURE_DIR = (
    ROOT
    / "reference_results"
    / "osteosarcoma_locked"
    / "preservation_fixture"
)

MANIFEST_FILE = (
    FIXTURE_DIR
    / "fixture_manifest.json"
)

REFERENCE_FILE = (
    FIXTURE_DIR
    / "DOG2_reference_module_expression.csv"
)

EXTERNAL_FILE = (
    FIXTURE_DIR
    / "GSE239948_external_module_expression.csv"
)

WEIGHTS_FILE = (
    FIXTURE_DIR
    / "primary_canine_program_weights.csv"
)

EXPECTED_COVERAGE_FILE = (
    FIXTURE_DIR
    / "expected_preservation_coverage.csv"
)

EXPECTED_DIRECT_FILE = (
    FIXTURE_DIR
    / "expected_direct_preservation.csv"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "preservation_regression"
)

OBSERVED_DIRECT_FILE = (
    OUTPUT_DIR
    / "GSE239948_observed_direct_preservation.csv"
)

OBSERVED_COVERAGE_FILE = (
    OUTPUT_DIR
    / "GSE239948_observed_preservation_coverage.csv"
)

VERIFICATION_FILE = (
    OUTPUT_DIR
    / "GSE239948_preservation_verification.json"
)

ABSOLUTE_TOLERANCE = 1.0e-10
RELATIVE_TOLERANCE = 1.0e-8

DIRECT_METRICS = [
    "edge_spearman",
    "loading_spearman",
    "external_pc1_variance_explained",
    (
        "pc1_orientation_correlation_"
        "with_frozen_score"
    ),
]


def sha256_file(
    path: Path,
) -> str:
    digest = hashlib.sha256()

    with path.open(
        "rb"
    ) as handle:
        for chunk in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):
            digest.update(
                chunk
            )

    return digest.hexdigest()


def load_manifest() -> dict[str, Any]:
    if not MANIFEST_FILE.exists():
        raise FileNotFoundError(
            "Preservation fixture manifest "
            "not found:\n"
            f"{MANIFEST_FILE}"
        )

    payload = json.loads(
        MANIFEST_FILE.read_text(
            encoding="utf-8"
        )
    )

    if not isinstance(
        payload,
        dict,
    ):
        raise ValueError(
            "Fixture manifest must contain "
            "a JSON object."
        )

    return payload


def verify_fixture_hashes(
    manifest: dict[str, Any],
) -> dict[str, bool]:
    fixture_files = manifest.get(
        "fixture_files",
        {},
    )

    if not isinstance(
        fixture_files,
        dict,
    ):
        raise ValueError(
            "fixture_files is missing from "
            "the preservation manifest."
        )

    checks: dict[str, bool] = {}

    for filename, metadata in (
        fixture_files.items()
    ):
        path = (
            FIXTURE_DIR
            / filename
        )

        if not path.exists():
            raise FileNotFoundError(
                f"Fixture file not found: "
                f"{path}"
            )

        if not isinstance(
            metadata,
            dict,
        ):
            raise ValueError(
                "Invalid fixture metadata "
                f"for {filename}."
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


def normalize_bool(
    value: Any,
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
        "Cannot interpret boolean "
        f"value: {value!r}"
    )


def normalize_text(
    value: Any,
) -> str:
    if pd.isna(value):
        return ""

    return str(value)


def ordered(
    table: pd.DataFrame,
    modules: list[str],
) -> pd.DataFrame:
    if "module_label" not in table.columns:
        raise ValueError(
            "Table is missing "
            "module_label."
        )

    indexed = table.copy()

    indexed[
        "module_label"
    ] = (
        indexed[
            "module_label"
        ]
        .astype(str)
    )

    indexed = (
        indexed.set_index(
            "module_label"
        )
        .reindex(
            modules
        )
    )

    if indexed.index.has_duplicates:
        raise ValueError(
            "Duplicate module labels "
            "were found."
        )

    if indexed.isna().all(
        axis=1
    ).any():
        missing = (
            indexed.index[
                indexed.isna().all(
                    axis=1
                )
            ]
            .tolist()
        )

        raise ValueError(
            "Missing module rows: "
            + ", ".join(
                missing
            )
        )

    return indexed.reset_index()


def compare_direct_metrics(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
    modules: list[str],
) -> tuple[
    bool,
    float,
    list[dict[str, Any]],
    list[str],
]:
    observed = ordered(
        observed,
        modules,
    )

    expected = ordered(
        expected,
        modules,
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
                f"{label} direct table is "
                "missing columns: "
                + ", ".join(
                    missing
                )
            )

    mismatches: list[str] = []

    summaries: list[
        dict[str, Any]
    ] = []

    max_difference = 0.0

    observed_counts = pd.to_numeric(
        observed[
            "n_common_genes"
        ],
        errors="coerce",
    ).to_numpy()

    expected_counts = pd.to_numeric(
        expected[
            "n_common_genes"
        ],
        errors="coerce",
    ).to_numpy()

    if not np.array_equal(
        observed_counts,
        expected_counts,
    ):
        mismatches.append(
            "n_common_genes"
        )

    observed_estimable = [
        normalize_bool(value)
        for value
        in observed[
            "estimable"
        ]
    ]

    expected_estimable = [
        normalize_bool(value)
        for value
        in expected[
            "estimable"
        ]
    ]

    if (
        observed_estimable
        != expected_estimable
    ):
        mismatches.append(
            "estimable"
        )

    observed_reason = [
        normalize_text(value)
        for value
        in observed[
            "nonestimable_reason"
        ]
    ]

    expected_reason = [
        normalize_text(value)
        for value
        in expected[
            "nonestimable_reason"
        ]
    ]

    if (
        observed_reason
        != expected_reason
    ):
        mismatches.append(
            "nonestimable_reason"
        )

    for module_index, module in enumerate(
        modules
    ):
        for metric in DIRECT_METRICS:
            observed_value = float(
                pd.to_numeric(
                    pd.Series(
                        [
                            observed.loc[
                                module_index,
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
                                module_index,
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

            max_difference = max(
                max_difference,
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
        max_difference,
        summaries,
        mismatches,
    )


def compare_coverage(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
    modules: list[str],
) -> tuple[
    bool,
    list[str],
]:
    observed = ordered(
        observed,
        modules,
    )

    expected = ordered(
        expected,
        modules,
    )

    numeric_integer = [
        "n_frozen_genes",
        "n_common_genes",
    ]

    numeric_float = [
        "coverage_fraction",
    ]

    text_columns = [
        "common_genes",
        "missing_external_genes",
    ]

    required = [
        "module_label",
        *numeric_integer,
        *numeric_float,
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
                f"{label} coverage table is "
                "missing columns: "
                + ", ".join(
                    missing
                )
            )

    mismatches: list[str] = []

    for column in numeric_integer:
        left = pd.to_numeric(
            observed[
                column
            ],
            errors="coerce",
        ).to_numpy()

        right = pd.to_numeric(
            expected[
                column
            ],
            errors="coerce",
        ).to_numpy()

        if not np.array_equal(
            left,
            right,
        ):
            mismatches.append(
                column
            )

    for column in numeric_float:
        left = pd.to_numeric(
            observed[
                column
            ],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        right = pd.to_numeric(
            expected[
                column
            ],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        if not np.allclose(
            left,
            right,
            rtol=(
                RELATIVE_TOLERANCE
            ),
            atol=(
                ABSOLUTE_TOLERANCE
            ),
            equal_nan=True,
        ):
            mismatches.append(
                column
            )

    for column in text_columns:
        left = [
            normalize_text(value)
            for value
            in observed[
                column
            ]
        ]

        right = [
            normalize_text(value)
            for value
            in expected[
                column
            ]
        ]

        if left != right:
            mismatches.append(
                column
            )

    return (
        len(mismatches) == 0,
        mismatches,
    )


def main() -> None:
    manifest = load_manifest()

    hash_checks = (
        verify_fixture_hashes(
            manifest
        )
    )

    if not all(
        hash_checks.values()
    ):
        failed = [
            filename
            for filename, passed
            in hash_checks.items()
            if not passed
        ]

        raise RuntimeError(
            "Preservation fixture hash "
            "verification failed: "
            + ", ".join(
                failed
            )
        )

    contract = manifest.get(
        "fixture_contract",
        {},
    )

    modules = [
        str(module)
        for module
        in contract.get(
            "primary_modules",
            [
                "M34",
                "M11",
                "M24",
                "M40",
            ],
        )
    ]

    minimum_genes = int(
        contract.get(
            "minimum_module_genes",
            3,
        )
    )

    if bool(
        contract.get(
            "outcome_loaded",
            True,
        )
    ):
        raise RuntimeError(
            "Fixture contract does not "
            "confirm outcome_loaded=false."
        )

    reference = pd.read_csv(
        REFERENCE_FILE,
        index_col=0,
    )

    external = pd.read_csv(
        EXTERNAL_FILE,
        index_col=0,
    )

    weights = pd.read_csv(
        WEIGHTS_FILE
    )

    expected_coverage = pd.read_csv(
        EXPECTED_COVERAGE_FILE
    )

    expected_direct = pd.read_csv(
        EXPECTED_DIRECT_FILE
    )

    reference.index = (
        reference.index
        .astype(str)
    )

    external.index = (
        external.index
        .astype(str)
    )

    (
        observed_direct,
        observed_coverage,
    ) = analyze_preservation(
        reference_expression=reference,
        external_expression=external,
        weights=weights,
        module_order=modules,
        minimum_genes=minimum_genes,
        pca_random_state=42,
    )

    (
        direct_match,
        max_difference,
        metric_summary,
        direct_mismatches,
    ) = compare_direct_metrics(
        observed=observed_direct,
        expected=expected_direct,
        modules=modules,
    )

    (
        coverage_match,
        coverage_mismatches,
    ) = compare_coverage(
        observed=observed_coverage,
        expected=expected_coverage,
        modules=modules,
    )

    passed = bool(
        all(
            hash_checks.values()
        )
        and direct_match
        and coverage_match
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    observed_direct.to_csv(
        OBSERVED_DIRECT_FILE,
        index=False,
    )

    observed_coverage.to_csv(
        OBSERVED_COVERAGE_FILE,
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
        "fixture_hashes_match":
            bool(
                all(
                    hash_checks.values()
                )
            ),
        "hash_checks":
            hash_checks,
        "contract": {
            "modules":
                modules,
            "minimum_genes":
                minimum_genes,
            "outcome_loaded":
                False,
            "scope":
                (
                    "deterministic_direct_"
                    "preservation"
                ),
        },
        "dimensions": {
            "reference_samples":
                int(
                    reference.shape[0]
                ),
            "external_samples":
                int(
                    external.shape[0]
                ),
            "fixture_genes":
                int(
                    reference.shape[1]
                ),
            "weight_rows":
                int(
                    weights.shape[0]
                ),
        },
        "comparison": {
            "direct_metrics_match":
                direct_match,
            "coverage_match":
                coverage_match,
            (
                "max_absolute_metric_"
                "difference"
            ):
                max_difference,
            "absolute_tolerance":
                ABSOLUTE_TOLERANCE,
            "relative_tolerance":
                RELATIVE_TOLERANCE,
            "direct_mismatches":
                direct_mismatches,
            "coverage_mismatches":
                coverage_mismatches,
            "metrics":
                metric_summary,
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

    display = (
        observed_direct[
            [
                "module_label",
                "n_common_genes",
                "edge_spearman",
                "loading_spearman",
                (
                    "external_pc1_"
                    "variance_explained"
                ),
                (
                    "pc1_orientation_"
                    "correlation_with_"
                    "frozen_score"
                ),
            ]
        ]
    )

    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- GSE239948 direct preservation regression"
    )
    print("=" * 80)

    print(
        "Reference samples: "
        f"{reference.shape[0]}"
    )

    print(
        "External samples: "
        f"{external.shape[0]}"
    )

    print(
        "Fixture genes: "
        f"{reference.shape[1]}"
    )

    print(
        "Frozen weight rows: "
        f"{weights.shape[0]}"
    )

    print(
        "Fixture hashes: "
        + (
            "PASS"
            if all(
                hash_checks.values()
            )
            else "FAIL"
        )
    )

    print(
        "Direct preservation regression: "
        + (
            "PASS"
            if direct_match
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
        f"{max_difference:.3e}"
    )

    print("")
    print(
        display.to_string(
            index=False
        )
    )

    if direct_mismatches:
        print("")
        print(
            "Direct mismatches: "
            + ", ".join(
                direct_mismatches
            )
        )

    if coverage_mismatches:
        print("")
        print(
            "Coverage mismatches: "
            + ", ".join(
                coverage_mismatches
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
