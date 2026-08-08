"""Reproduce locked GSE239948 edge/loading permutation inference."""

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

from transport_audit.multiplicity import (
    benjamini_hochberg,
)
from transport_audit.preservation_inference import (
    analyze_permutation_preservation,
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

EXPECTED_FILE = (
    FIXTURE_DIR
    / "expected_direct_preservation.csv"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "preservation_permutation_regression"
)

OBSERVED_FILE = (
    OUTPUT_DIR
    / "GSE239948_observed_permutation_preservation.csv"
)

VERIFICATION_FILE = (
    OUTPUT_DIR
    / "GSE239948_permutation_verification.json"
)

MODULES = [
    "M34",
    "M11",
    "M24",
    "M40",
]

N_PERMUTATIONS = 5000
BASE_SEED = 42
MODULE_SEED_STRIDE = 100

NUMERIC_TOLERANCE = 1.0e-10
P_VALUE_TOLERANCE = 1.0e-15


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


def load_manifest() -> dict[str, Any]:
    return json.loads(
        MANIFEST_FILE.read_text(
            encoding="utf-8"
        )
    )


def verify_fixture_hashes(
    manifest: dict[str, Any],
) -> dict[str, bool]:
    records = manifest.get(
        "fixture_files",
        {},
    )

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

        if not path.exists():
            raise FileNotFoundError(
                f"Missing fixture file: "
                f"{path}"
            )

        expected_sha = str(
            metadata.get(
                "sha256",
                "",
            )
        )

        observed_sha = (
            sha256_file(
                path
            )
        )

        checks[
            filename
        ] = bool(
            expected_sha
            and expected_sha
            == observed_sha
        )

    return checks


def ordered(
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

    result = (
        result.set_index(
            "module_label"
        )
        .reindex(
            MODULES
        )
        .reset_index()
    )

    return result


def expected_extreme_count(
    p_value: float,
) -> int:
    raw = (
        float(p_value)
        * (N_PERMUTATIONS + 1)
        - 1.0
    )

    rounded = int(
        round(raw)
    )

    if abs(
        raw - rounded
    ) > 1.0e-8:
        raise RuntimeError(
            "Locked permutation p-value "
            "cannot be represented by the "
            "expected Monte-Carlo count: "
            f"p={p_value}, raw_count={raw}"
        )

    return rounded


def main() -> None:
    print("=" * 80)
    print(
        "Molecular Transport Audit "
        "- GSE239948 permutation preservation regression"
    )
    print("=" * 80)

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
            "Preservation fixture hash "
            "verification failed."
        )

    contract = manifest[
        "fixture_contract"
    ]

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

    expected = ordered(
        pd.read_csv(
            EXPECTED_FILE
        )
    )

    reference.index = (
        reference.index.astype(str)
    )

    external.index = (
        external.index.astype(str)
    )

    observed = ordered(
        analyze_permutation_preservation(
            reference_expression=reference,
            external_expression=external,
            weights=weights,
            module_order=MODULES,
            minimum_genes=3,
            n_permutations=(
                N_PERMUTATIONS
            ),
            base_seed=BASE_SEED,
            module_seed_stride=(
                MODULE_SEED_STRIDE
            ),
            pca_random_state=42,
        )
    )

    direct_metrics = [
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

    metric_checks: list[
        dict[str, Any]
    ] = []

    maximum_metric_difference = 0.0
    metrics_match = True

    for row_index, module in enumerate(
        MODULES
    ):
        for metric in direct_metrics:
            observed_value = float(
                observed.loc[
                    row_index,
                    metric,
                ]
            )

            expected_value = float(
                expected.loc[
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
                    rtol=1.0e-8,
                    atol=(
                        NUMERIC_TOLERANCE
                    ),
                )
            )

            maximum_metric_difference = max(
                maximum_metric_difference,
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

    count_checks: list[
        dict[str, Any]
    ] = []

    counts_match = True
    p_values_match = True

    for row_index, module in enumerate(
        MODULES
    ):
        for family in [
            "edge",
            "loading",
        ]:
            p_column = (
                f"{family}_permutation_p"
            )

            count_column = (
                f"{family}_extreme_count"
            )

            observed_p = float(
                observed.loc[
                    row_index,
                    p_column,
                ]
            )

            expected_p = float(
                expected.loc[
                    row_index,
                    p_column,
                ]
            )

            observed_count = int(
                observed.loc[
                    row_index,
                    count_column,
                ]
            )

            locked_count = (
                expected_extreme_count(
                    expected_p
                )
            )

            count_match = (
                observed_count
                == locked_count
            )

            p_match = bool(
                np.isclose(
                    observed_p,
                    expected_p,
                    rtol=0.0,
                    atol=(
                        P_VALUE_TOLERANCE
                    ),
                )
            )

            counts_match = (
                counts_match
                and count_match
            )

            p_values_match = (
                p_values_match
                and p_match
            )

            count_checks.append(
                {
                    "module_label":
                        module,
                    "test":
                        family,
                    "observed_count":
                        observed_count,
                    "expected_count":
                        locked_count,
                    "observed_p":
                        observed_p,
                    "expected_p":
                        expected_p,
                    "count_match":
                        count_match,
                    "p_match":
                        p_match,
                }
            )

    expected_family = np.concatenate(
        [
            pd.to_numeric(
                expected[
                    "edge_permutation_p"
                ],
                errors="coerce",
            ).to_numpy(
                dtype=float
            ),
            pd.to_numeric(
                expected[
                    "loading_permutation_p"
                ],
                errors="coerce",
            ).to_numpy(
                dtype=float
            ),
        ]
    )

    expected_q = (
        benjamini_hochberg(
            expected_family
        )
    )

    observed_q = np.concatenate(
        [
            observed[
                "edge_q_bh_8"
            ].to_numpy(
                dtype=float
            ),
            observed[
                "loading_q_bh_8"
            ].to_numpy(
                dtype=float
            ),
        ]
    )

    q_values_match = bool(
        np.allclose(
            observed_q,
            expected_q,
            rtol=0.0,
            atol=1.0e-15,
            equal_nan=True,
        )
    )

    seed_expected = [
        (42, 43),
        (142, 143),
        (242, 243),
        (342, 343),
    ]

    seeds_match = True

    for row_index, (
        expected_edge_seed,
        expected_loading_seed,
    ) in enumerate(
        seed_expected
    ):
        seeds_match = (
            seeds_match
            and int(
                observed.loc[
                    row_index,
                    "edge_seed",
                ]
            )
            == expected_edge_seed
            and int(
                observed.loc[
                    row_index,
                    "loading_seed",
                ]
            )
            == expected_loading_seed
        )

    passed = bool(
        all(
            hash_checks.values()
        )
        and metrics_match
        and counts_match
        and p_values_match
        and q_values_match
        and seeds_match
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    observed.to_csv(
        OBSERVED_FILE,
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
        "n_permutations":
            N_PERMUTATIONS,
        "base_seed":
            BASE_SEED,
        "module_seed_stride":
            MODULE_SEED_STRIDE,
        "outcome_loaded":
            False,
        "checks": {
            "direct_metrics_match":
                metrics_match,
            "permutation_counts_match":
                counts_match,
            "permutation_p_values_match":
                p_values_match,
            "bh_q_values_match":
                q_values_match,
            "seed_schedule_match":
                seeds_match,
        },
        "maximum_direct_metric_difference":
            maximum_metric_difference,
        "metric_checks":
            metric_checks,
        "permutation_checks":
            count_checks,
        "expected_q_family":
            expected_q.tolist(),
        "observed_q_family":
            observed_q.tolist(),
        "observed_artifact": {
            "path":
                str(
                    OBSERVED_FILE
                    .relative_to(ROOT)
                ),
            "sha256":
                sha256_file(
                    OBSERVED_FILE
                ),
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
        f"Reference samples: "
        f"{reference.shape[0]}"
    )

    print(
        f"External samples: "
        f"{external.shape[0]}"
    )

    print(
        f"Permutations per test: "
        f"{N_PERMUTATIONS}"
    )

    print(
        "Fixture hashes: PASS"
    )

    print(
        "Direct metrics: "
        + (
            "PASS"
            if metrics_match
            else "FAIL"
        )
    )

    print(
        "Exact permutation counts: "
        + (
            "PASS"
            if counts_match
            else "FAIL"
        )
    )

    print(
        "Permutation p-values: "
        + (
            "PASS"
            if p_values_match
            else "FAIL"
        )
    )

    print(
        "BH across 8 direct tests: "
        + (
            "PASS"
            if q_values_match
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
        "Maximum direct metric "
        "difference: "
        f"{maximum_metric_difference:.3e}"
    )

    print("")

    display_columns = [
        "module_label",
        "edge_spearman",
        "edge_permutation_p",
        "edge_q_bh_8",
        "edge_extreme_count",
        "edge_seed",
        "loading_spearman",
        "loading_permutation_p",
        "loading_q_bh_8",
        "loading_extreme_count",
        "loading_seed",
    ]

    print(
        observed[
            display_columns
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
