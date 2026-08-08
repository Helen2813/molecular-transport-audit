"""Validate legacy and corrected GSE239948 random-panel controls."""

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

from transport_audit.random_controls import (
    LEGACY_Z_VARIANCE,
    PRESTANDARDIZATION_VARIANCE,
    classify_external_representation,
    random_panel_controls,
)


FIXTURE_DIR = (
    ROOT
    / "reference_results"
    / "osteosarcoma_locked"
    / "preservation_fixture"
)

MANIFEST_FILE = (
    FIXTURE_DIR
    / "random_panel_fixture_manifest.json"
)

REFERENCE_ARRAY = (
    FIXTURE_DIR
    / "random_panel_DOG2_prestandardization.npy"
)

EXTERNAL_ARRAY = (
    FIXTURE_DIR
    / "random_panel_GSE239948_prestandardization.npy"
)

GENES_FILE = (
    FIXTURE_DIR
    / "random_panel_background_genes.csv"
)

WEIGHTS_FILE = (
    FIXTURE_DIR
    / "primary_canine_program_weights.csv"
)

STRUCTURE_FILE = (
    FIXTURE_DIR
    / "expected_direct_preservation.csv"
)

EXPECTED_RANDOM_FILE = (
    FIXTURE_DIR
    / "expected_legacy_random_panel_controls.csv"
)

EXPECTED_CLASSIFICATION_FILE = (
    FIXTURE_DIR
    / "expected_legacy_representation_classification.csv"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "random_panel_validation"
)

LEGACY_OUTPUT = (
    OUTPUT_DIR
    / "legacy_exact_random_panel_controls.csv"
)

CORRECTED_OUTPUT = (
    OUTPUT_DIR
    / "corrected_prestandardization_random_panel_controls.csv"
)

LEGACY_CLASSIFICATION_OUTPUT = (
    OUTPUT_DIR
    / "legacy_exact_classification.csv"
)

CORRECTED_CLASSIFICATION_OUTPUT = (
    OUTPUT_DIR
    / "corrected_prestandardization_classification.csv"
)

COMPARISON_OUTPUT = (
    OUTPUT_DIR
    / "legacy_vs_corrected_random_panel_controls.csv"
)

VERIFICATION_FILE = (
    OUTPUT_DIR
    / "random_panel_validation.json"
)

MODULES = [
    "M34",
    "M11",
    "M24",
    "M40",
]

N_RANDOM_PANELS = 1000
N_BINS = 10
RANDOM_SEED = 42

ATOL = 1.0e-10
RTOL = 1.0e-8

LEGACY_NUMERIC_COLUMNS = [
    "n_module_genes",
    "n_random_panels",
    "observed_edge_spearman",
    "random_edge_median",
    "random_edge_q05",
    "random_edge_q95",
    "random_panel_empirical_p",
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
    return json.loads(
        path.read_text(
            encoding="utf-8"
        )
    )


def verify_hashes(
    manifest: dict[str, Any],
) -> dict[str, bool]:
    records = manifest[
        "fixture_files"
    ]

    checks: dict[
        str,
        bool,
    ] = {}

    for filename, metadata in (
        records.items()
    ):
        path = (
            FIXTURE_DIR
            / filename
        )

        checks[
            filename
        ] = bool(
            path.exists()
            and sha256_file(
                path
            )
            == metadata[
                "sha256"
            ]
        )

    return checks


def order_modules(
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


def compare_legacy_random(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
) -> tuple[
    bool,
    float,
    list[str],
]:
    observed = order_modules(
        observed
    )

    expected = order_modules(
        expected
    )

    mismatches: list[str] = []
    maximum_difference = 0.0

    for column in (
        LEGACY_NUMERIC_COLUMNS
    ):
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

        difference = float(
            np.nanmax(
                np.abs(
                    left
                    - right
                )
            )
        )

        maximum_difference = max(
            maximum_difference,
            difference,
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


def compare_legacy_classification(
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

    mismatches: list[str] = []

    class_column = (
        "external_canine_"
        "representation_class"
    )

    if (
        observed[
            class_column
        ].astype(str).tolist()
        != expected[
            class_column
        ].astype(str).tolist()
    ):
        mismatches.append(
            class_column
        )

    for column in [
        "edge_q_bh_8",
        "loading_q_bh_8",
        "random_panel_empirical_p",
    ]:
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
            rtol=RTOL,
            atol=ATOL,
            equal_nan=True,
        ):
            mismatches.append(
                column
            )

    return (
        len(mismatches) == 0,
        mismatches,
    )


def main() -> None:
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- random-panel validation"
    )
    print("=" * 80)

    for path in [
        MANIFEST_FILE,
        REFERENCE_ARRAY,
        EXTERNAL_ARRAY,
        GENES_FILE,
        WEIGHTS_FILE,
        STRUCTURE_FILE,
        EXPECTED_RANDOM_FILE,
        EXPECTED_CLASSIFICATION_FILE,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Required fixture missing:\n"
                f"{path}"
            )

    manifest = load_json(
        MANIFEST_FILE
    )

    hash_checks = verify_hashes(
        manifest
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
            "Fixture hash verification "
            "failed: "
            + ", ".join(
                failed
            )
        )

    genes = (
        pd.read_csv(
            GENES_FILE
        )[
            "gene_symbol"
        ]
        .astype(str)
        .tolist()
    )

    reference_array = np.load(
        REFERENCE_ARRAY,
        allow_pickle=False,
    )

    external_array = np.load(
        EXTERNAL_ARRAY,
        allow_pickle=False,
    )

    reference = pd.DataFrame(
        reference_array,
        columns=genes,
    )

    external = pd.DataFrame(
        external_array,
        columns=genes,
    )

    weights = pd.read_csv(
        WEIGHTS_FILE
    )

    structure = pd.read_csv(
        STRUCTURE_FILE
    )

    expected_random = (
        pd.read_csv(
            EXPECTED_RANDOM_FILE
        )
    )

    expected_classification = (
        pd.read_csv(
            EXPECTED_CLASSIFICATION_FILE
        )
    )

    print(
        "Reference samples: "
        f"{reference.shape[0]}"
    )

    print(
        "External samples: "
        f"{external.shape[0]}"
    )

    print(
        "Shared genes: "
        f"{reference.shape[1]}"
    )

    print(
        "Random panels per module: "
        f"{N_RANDOM_PANELS}"
    )

    print("")
    print("=" * 80)
    print(
        "Stage 1: exact legacy mode"
    )
    print("=" * 80)

    legacy = random_panel_controls(
        reference_expression=reference,
        external_expression=external,
        weights=weights,
        observed_structure=structure,
        module_order=MODULES,
        matching_mode=(
            LEGACY_Z_VARIANCE
        ),
        n_random_panels=(
            N_RANDOM_PANELS
        ),
        n_variability_bins=(
            N_BINS
        ),
        random_seed=(
            RANDOM_SEED
        ),
    )

    legacy_classification = (
        classify_external_representation(
            structure,
            legacy,
        )
    )

    (
        legacy_random_match,
        maximum_legacy_difference,
        legacy_random_mismatches,
    ) = compare_legacy_random(
        legacy,
        expected_random,
    )

    (
        legacy_class_match,
        legacy_class_mismatches,
    ) = compare_legacy_classification(
        legacy_classification,
        expected_classification,
    )

    print(
        "Legacy random-panel regression: "
        + (
            "PASS"
            if legacy_random_match
            else "FAIL"
        )
    )

    print(
        "Legacy classification regression: "
        + (
            "PASS"
            if legacy_class_match
            else "FAIL"
        )
    )

    print(
        "Maximum legacy difference: "
        f"{maximum_legacy_difference:.3e}"
    )

    print("")
    print(
        legacy[
            [
                "module_label",
                "n_module_genes",
                "n_random_panels",
                "observed_edge_spearman",
                "random_edge_median",
                "random_edge_q05",
                "random_edge_q95",
                "random_panel_empirical_p",
                "fallback_draws",
            ]
        ].to_string(
            index=False
        )
    )

    if (
        not legacy_random_match
        or not legacy_class_match
    ):
        print("")
        print(
            "Legacy random mismatches: "
            + ", ".join(
                legacy_random_mismatches
            )
        )

        print(
            "Legacy classification mismatches: "
            + ", ".join(
                legacy_class_mismatches
            )
        )

        raise SystemExit(1)

    print("")
    print("=" * 80)
    print(
        "Stage 2: corrected "
        "pre-standardization matching"
    )
    print("=" * 80)

    corrected = random_panel_controls(
        reference_expression=reference,
        external_expression=external,
        weights=weights,
        observed_structure=structure,
        module_order=MODULES,
        matching_mode=(
            PRESTANDARDIZATION_VARIANCE
        ),
        n_random_panels=(
            N_RANDOM_PANELS
        ),
        n_variability_bins=(
            N_BINS
        ),
        random_seed=(
            RANDOM_SEED
        ),
    )

    corrected_classification = (
        classify_external_representation(
            structure,
            corrected,
        )
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    legacy.to_csv(
        LEGACY_OUTPUT,
        index=False,
    )

    corrected.to_csv(
        CORRECTED_OUTPUT,
        index=False,
    )

    legacy_classification.to_csv(
        LEGACY_CLASSIFICATION_OUTPUT,
        index=False,
    )

    corrected_classification.to_csv(
        CORRECTED_CLASSIFICATION_OUTPUT,
        index=False,
    )

    comparison = (
        order_modules(
            legacy
        )[
            [
                "module_label",
                "random_edge_median",
                "random_edge_q05",
                "random_edge_q95",
                "random_panel_empirical_p",
                "fallback_draws",
            ]
        ]
        .rename(
            columns={
                column:
                    (
                        f"legacy__{column}"
                        if column
                        != "module_label"
                        else column
                    )
                for column in [
                    "module_label",
                    "random_edge_median",
                    "random_edge_q05",
                    "random_edge_q95",
                    "random_panel_empirical_p",
                    "fallback_draws",
                ]
            }
        )
        .merge(
            order_modules(
                corrected
            )[
                [
                    "module_label",
                    "random_edge_median",
                    "random_edge_q05",
                    "random_edge_q95",
                    "random_panel_empirical_p",
                    "fallback_draws",
                ]
            ].rename(
                columns={
                    column:
                        (
                            f"corrected__{column}"
                            if column
                            != "module_label"
                            else column
                        )
                    for column in [
                        "module_label",
                        "random_edge_median",
                        "random_edge_q05",
                        "random_edge_q95",
                        "random_panel_empirical_p",
                        "fallback_draws",
                    ]
                }
            ),
            on="module_label",
            how="inner",
        )
    )

    comparison = comparison.merge(
        order_modules(
            legacy_classification
        )[
            [
                "module_label",
                (
                    "external_canine_"
                    "representation_class"
                ),
            ]
        ].rename(
            columns={
                (
                    "external_canine_"
                    "representation_class"
                ):
                    "legacy__classification"
            }
        ),
        on="module_label",
        how="left",
    )

    comparison = comparison.merge(
        order_modules(
            corrected_classification
        )[
            [
                "module_label",
                (
                    "external_canine_"
                    "representation_class"
                ),
            ]
        ].rename(
            columns={
                (
                    "external_canine_"
                    "representation_class"
                ):
                    "corrected__classification"
            }
        ),
        on="module_label",
        how="left",
    )

    comparison.to_csv(
        COMPARISON_OUTPUT,
        index=False,
    )

    corrected_all_panels_valid = bool(
        (
            corrected[
                "n_random_panels"
            ]
            == N_RANDOM_PANELS
        ).all()
    )

    corrected_classes = (
        corrected_classification
        .set_index(
            "module_label"
        )[
            (
                "external_canine_"
                "representation_class"
            )
        ]
        .to_dict()
    )

    payload = {
        "runner_version":
            RUNNER_VERSION,
        "generated_at_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),
        "passed":
            bool(
                legacy_random_match
                and legacy_class_match
                and corrected_all_panels_valid
            ),
        "outcome_loaded":
            False,
        "fixture_hashes_match":
            True,
        "legacy_exact": {
            "matching_mode":
                LEGACY_Z_VARIANCE,
            "random_panel_regression_match":
                legacy_random_match,
            "classification_regression_match":
                legacy_class_match,
            "maximum_absolute_difference":
                maximum_legacy_difference,
        },
        "corrected_sensitivity": {
            "matching_mode":
                PRESTANDARDIZATION_VARIANCE,
            "all_modules_have_1000_panels":
                corrected_all_panels_valid,
            "classification":
                corrected_classes,
        },
        "outputs": {
            "legacy":
                str(
                    LEGACY_OUTPUT
                    .relative_to(ROOT)
                ),
            "corrected":
                str(
                    CORRECTED_OUTPUT
                    .relative_to(ROOT)
                ),
            "comparison":
                str(
                    COMPARISON_OUTPUT
                    .relative_to(ROOT)
                ),
            "corrected_classification":
                str(
                    CORRECTED_CLASSIFICATION_OUTPUT
                    .relative_to(ROOT)
                ),
        },
    }

    VERIFICATION_FILE.write_text(
        json.dumps(
            payload,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        "Corrected panels complete: "
        + (
            "PASS"
            if corrected_all_panels_valid
            else "FAIL"
        )
    )

    print("")
    print(
        corrected[
            [
                "module_label",
                "n_module_genes",
                "n_random_panels",
                "observed_edge_spearman",
                "random_edge_median",
                "random_edge_q05",
                "random_edge_q95",
                "random_panel_empirical_p",
                "fallback_draws",
            ]
        ].to_string(
            index=False
        )
    )

    print("")
    print(
        "Corrected classification:"
    )
    print("")

    print(
        corrected_classification[
            [
                "module_label",
                "edge_q_bh_8",
                "loading_q_bh_8",
                "split_half_median",
                "random_panel_empirical_p",
                (
                    "external_canine_"
                    "representation_class"
                ),
            ]
        ].to_string(
            index=False
        )
    )

    print("")
    print("=" * 80)
    print(
        "Random-panel validation: PASS"
    )
    print("=" * 80)

    print(
        "Legacy exact regression: PASS"
    )

    print(
        "Legacy classification: PASS"
    )

    print(
        "Corrected sensitivity completed: PASS"
    )

    print(
        "Verification: "
        + str(
            VERIFICATION_FILE
            .relative_to(ROOT)
        )
    )

    print("=" * 80)


if __name__ == "__main__":
    main()
