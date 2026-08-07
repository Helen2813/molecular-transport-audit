"""Run the locked TARGET scoring regression."""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

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

from transport_audit.scoring import (
    score_programs,
)
from transport_audit.schemas import (
    ScoringRule,
)


FIXTURE_DIR = (
    ROOT
    / "reference_results"
    / "osteosarcoma_locked"
    / "scoring_fixture"
)

MANIFEST_FILE = (
    FIXTURE_DIR
    / "fixture_manifest.json"
)

EXPRESSION_FILE = (
    FIXTURE_DIR
    / "TARGET_OS_primary_strict_expression.csv"
)

WEIGHTS_FILE = (
    FIXTURE_DIR
    / "primary_strict_gene_weights.csv"
)

EXPECTED_SCORES_FILE = (
    FIXTURE_DIR
    / "TARGET_OS_expected_primary_scores.csv"
)

EXPECTED_COVERAGE_FILE = (
    FIXTURE_DIR
    / "TARGET_OS_expected_strict_coverage.csv"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "scoring_regression"
)

OBSERVED_SCORES_FILE = (
    OUTPUT_DIR
    / "TARGET_OS_observed_primary_scores.csv"
)

OBSERVED_COVERAGE_FILE = (
    OUTPUT_DIR
    / "TARGET_OS_observed_strict_coverage.csv"
)

VERIFICATION_FILE = (
    OUTPUT_DIR
    / "TARGET_OS_scoring_verification.json"
)

ABSOLUTE_TOLERANCE = 1.0e-10
RELATIVE_TOLERANCE = 1.0e-8


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
            digest.update(chunk)

    return digest.hexdigest()


def load_manifest() -> dict:
    if not MANIFEST_FILE.exists():
        raise FileNotFoundError(
            "Fixture manifest not found: "
            f"{MANIFEST_FILE}"
        )

    return json.loads(
        MANIFEST_FILE.read_text(
            encoding="utf-8"
        )
    )


def verify_fixture_hashes(
    manifest: dict,
) -> dict[str, bool]:
    file_map = {
        EXPRESSION_FILE.name:
            EXPRESSION_FILE,
        WEIGHTS_FILE.name:
            WEIGHTS_FILE,
        EXPECTED_SCORES_FILE.name:
            EXPECTED_SCORES_FILE,
        EXPECTED_COVERAGE_FILE.name:
            EXPECTED_COVERAGE_FILE,
    }

    expected = manifest.get(
        "fixture_files",
        {},
    )

    checks: dict[
        str,
        bool,
    ] = {}

    for name, path in file_map.items():
        if not path.exists():
            raise FileNotFoundError(
                f"Fixture file not found: "
                f"{path}"
            )

        expected_sha = str(
            expected.get(
                name,
                {},
            ).get(
                "sha256",
                "",
            )
        )

        observed_sha = (
            sha256_file(path)
        )

        checks[name] = bool(
            expected_sha
            and observed_sha
            == expected_sha
        )

    return checks


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
        "Cannot interpret boolean "
        f"value: {value!r}"
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

    summaries = []

    all_match = True
    max_absolute = 0.0

    for column in expected.columns:
        observed_values = (
            pd.to_numeric(
                observed[column],
                errors="coerce",
            )
            .to_numpy(dtype=float)
        )

        expected_values = (
            pd.to_numeric(
                expected[column],
                errors="coerce",
            )
            .to_numpy(dtype=float)
        )

        same_nan = np.array_equal(
            np.isnan(
                observed_values
            ),
            np.isnan(
                expected_values
            ),
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
            column_max = float(
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
            column_max = 0.0

        column_match = (
            same_nan
            and np.allclose(
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

        max_absolute = max(
            max_absolute,
            column_max,
        )

        all_match = (
            all_match
            and column_match
        )

        summaries.append(
            {
                "score_column":
                    column,
                "matched":
                    bool(column_match),
                "max_absolute_difference":
                    column_max,
            }
        )

    return (
        all_match,
        max_absolute,
        summaries,
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

    missing_expected = [
        column
        for column
        in key_columns
        + compare_columns
        if column
        not in expected.columns
    ]

    if missing_expected:
        raise ValueError(
            "Expected coverage table is "
            "missing columns: "
            + ", ".join(
                missing_expected
            )
        )

    missing_observed = [
        column
        for column
        in key_columns
        + compare_columns
        if column
        not in observed.columns
    ]

    if missing_observed:
        raise ValueError(
            "Observed coverage table is "
            "missing columns: "
            + ", ".join(
                missing_observed
            )
        )

    left = expected[
        key_columns
        + compare_columns
    ].copy()

    right = observed[
        key_columns
        + compare_columns
    ].copy()

    left = (
        left.sort_values(
            key_columns
        )
        .reset_index(drop=True)
    )

    right = (
        right.sort_values(
            key_columns
        )
        .reset_index(drop=True)
    )

    mismatches: list[
        str
    ] = []

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
        mismatches.append(
            "coverage_key_family"
        )

        return (
            False,
            mismatches,
        )

    for column in [
        "n_frozen_genes",
        "n_available_genes",
    ]:
        left_values = (
            pd.to_numeric(
                left[column],
                errors="coerce",
            )
            .to_numpy()
        )

        right_values = (
            pd.to_numeric(
                right[column],
                errors="coerce",
            )
            .to_numpy()
        )

        if not np.array_equal(
            left_values,
            right_values,
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

    left_bool = [
        normalize_bool(value)
        for value
        in left[
            "minimum_rule_passed"
        ]
    ]

    right_bool = [
        normalize_bool(value)
        for value
        in right[
            "minimum_rule_passed"
        ]
    ]

    if left_bool != right_bool:
        mismatches.append(
            "minimum_rule_passed"
        )

    for column in [
        "available_genes",
        "missing_genes",
    ]:
        left_text = [
            normalize_text(value)
            for value
            in left[column]
        ]

        right_text = [
            normalize_text(value)
            for value
            in right[column]
        ]

        if left_text != right_text:
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
        raise RuntimeError(
            "One or more scoring fixture "
            "hashes do not match the "
            "locked fixture manifest."
        )

    contract = manifest[
        "scoring_contract"
    ]

    cohort = str(
        manifest["cohort"]
    )

    modules = [
        str(module)
        for module
        in manifest[
            "primary_modules"
        ]
    ]

    expression = pd.read_csv(
        EXPRESSION_FILE,
        index_col=0,
    )

    weights = pd.read_csv(
        WEIGHTS_FILE,
    )

    expected_scores = pd.read_csv(
        EXPECTED_SCORES_FILE,
        index_col=0,
    )

    expected_coverage = (
        pd.read_csv(
            EXPECTED_COVERAGE_FILE
        )
    )

    expression.index = (
        expression.index
        .astype(str)
    )

    expected_scores.index = (
        expected_scores.index
        .astype(str)
    )

    rule = ScoringRule(
        mapping_name="strict",
        minimum_genes=int(
            contract[
                "minimum_genes"
            ]
        ),
        minimum_fraction=float(
            contract[
                "minimum_fraction"
            ]
        ),
    )

    (
        observed_scores,
        observed_coverage,
    ) = score_programs(
        expression=expression,
        weights=weights,
        cohort=cohort,
        module_order=modules,
        rule=rule,
    )

    observed_scores = (
        observed_scores.reindex(
            index=(
                expected_scores.index
            ),
            columns=(
                expected_scores.columns
            ),
        )
    )

    (
        score_match,
        max_diff,
        score_summary,
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
        score_match
        and coverage_match
        and all(
            hash_checks.values()
        )
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    observed_scores.to_csv(
        OBSERVED_SCORES_FILE,
        index=True,
        index_label=(
            observed_scores.index.name
            or "sample_id"
        ),
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
        "scoring_contract":
            contract,
        "comparison": {
            "cohort":
                cohort,
            "modules":
                modules,
            "n_samples":
                int(
                    expected_scores
                    .shape[0]
                ),
            "n_score_columns":
                int(
                    expected_scores
                    .shape[1]
                ),
            "scores_match":
                bool(
                    score_match
                ),
            "coverage_match":
                bool(
                    coverage_match
                ),
            "max_absolute_score_difference":
                max_diff,
            "absolute_tolerance":
                ABSOLUTE_TOLERANCE,
            "relative_tolerance":
                RELATIVE_TOLERANCE,
            "coverage_mismatches":
                coverage_mismatches,
            "score_columns":
                score_summary,
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

    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- TARGET scoring regression"
    )
    print("=" * 80)

    print(
        f"Cohort: {cohort}"
    )

    print(
        "Samples: "
        f"{expected_scores.shape[0]}"
    )

    print(
        "Score columns: "
        f"{expected_scores.shape[1]}"
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
        "Score regression: "
        + (
            "PASS"
            if score_match
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
        "Maximum absolute score "
        "difference: "
        f"{max_diff:.3e}"
    )

    print("")

    print(
        pd.DataFrame(
            score_summary
        ).to_string(
            index=False
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
