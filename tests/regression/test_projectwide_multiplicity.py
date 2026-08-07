from __future__ import annotations

import numpy as np
import pandas as pd

from transport_audit.multiplicity import (
    benjamini_hochberg,
    significance_mask,
)
from transport_audit.reference import (
    get_reference_artifact_path,
    load_tolerance_policy,
)


REFERENCE_ID = (
    "projectwide_primary_multiplicity"
)


def test_projectwide_bh_reproduces_locked_result() -> None:
    reference_path = (
        get_reference_artifact_path(
            REFERENCE_ID
        )
    )

    table = pd.read_csv(
        reference_path
    )

    assert table.shape[0] == 12

    required_columns = {
        "cohort",
        "endpoint",
        "module_label",
        "primary_p",
        "projectwide_q_12",
        "projectwide_fdr_supported",
    }

    missing = (
        required_columns
        - set(table.columns)
    )

    assert not missing, (
        "Locked reference table is missing "
        f"required columns: {sorted(missing)}"
    )

    observed_q = (
        benjamini_hochberg(
            table["primary_p"].to_numpy(
                dtype=float
            )
        )
    )

    policy = load_tolerance_policy()

    deterministic = policy[
        "classes"
    ]["deterministic_table"]

    atol = float(
        deterministic[
            "absolute_tolerance"
        ]
    )

    rtol = float(
        deterministic[
            "relative_tolerance"
        ]
    )

    expected_q = table[
        "projectwide_q_12"
    ].to_numpy(
        dtype=float
    )

    np.testing.assert_allclose(
        observed_q,
        expected_q,
        atol=atol,
        rtol=rtol,
    )


def test_projectwide_fdr_classification_is_unchanged() -> None:
    reference_path = (
        get_reference_artifact_path(
            REFERENCE_ID
        )
    )

    table = pd.read_csv(
        reference_path
    )

    observed_q = (
        benjamini_hochberg(
            table["primary_p"].to_numpy(
                dtype=float
            )
        )
    )

    observed_supported = (
        significance_mask(
            observed_q,
            alpha=0.05,
        )
    )

    expected_supported = table[
        "projectwide_fdr_supported"
    ].astype(bool).to_numpy()

    np.testing.assert_array_equal(
        observed_supported,
        expected_supported,
    )


def test_locked_family_contains_expected_program_cohort_pairs() -> None:
    reference_path = (
        get_reference_artifact_path(
            REFERENCE_ID
        )
    )

    table = pd.read_csv(
        reference_path
    )

    observed_pairs = set(
        zip(
            table["cohort"],
            table["module_label"],
            strict=True,
        )
    )

    expected_cohorts = {
        "TARGET_OS",
        "GSE21257",
        "GSE39055",
    }

    expected_modules = {
        "M34",
        "M11",
        "M24",
        "M40",
    }

    expected_pairs = {
        (
            cohort,
            module,
        )
        for cohort in expected_cohorts
        for module in expected_modules
    }

    assert observed_pairs == expected_pairs


def test_only_locked_m34_gse21257_result_survives_fdr() -> None:
    reference_path = (
        get_reference_artifact_path(
            REFERENCE_ID
        )
    )

    table = pd.read_csv(
        reference_path
    )

    observed_q = (
        benjamini_hochberg(
            table["primary_p"].to_numpy(
                dtype=float
            )
        )
    )

    supported = table.loc[
        observed_q < 0.05,
        [
            "cohort",
            "module_label",
        ],
    ]

    observed = {
        tuple(row)
        for row in supported.to_numpy()
    }

    assert observed == {
        (
            "GSE21257",
            "M34",
        )
    }
