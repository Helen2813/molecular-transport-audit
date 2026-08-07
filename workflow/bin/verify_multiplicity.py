from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


VERIFIER_VERSION = "0.1.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Verify a multiplicity result against "
            "a locked scientific reference."
        )
    )

    parser.add_argument(
        "--observed",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--reference",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--tolerances",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--keys",
        default=(
            "cohort,endpoint,module_label"
        ),
    )

    parser.add_argument(
        "--q-column",
        default="projectwide_q_12",
    )

    parser.add_argument(
        "--support-column",
        default=(
            "projectwide_fdr_supported"
        ),
    )

    return parser.parse_args()


def parse_boolean_series(
    values: pd.Series,
) -> pd.Series:
    if pd.api.types.is_bool_dtype(
        values
    ):
        return values.astype(bool)

    mapping = {
        "true": True,
        "false": False,
        "1": True,
        "0": False,
    }

    normalized = (
        values.astype(str)
        .str.strip()
        .str.lower()
    )

    unknown = (
        set(normalized.unique())
        - set(mapping)
    )

    if unknown:
        raise ValueError(
            "Unable to parse boolean values: "
            f"{sorted(unknown)}"
        )

    return normalized.map(
        mapping
    ).astype(bool)


def main() -> None:
    args = parse_args()

    observed = pd.read_csv(
        args.observed
    )

    reference = pd.read_csv(
        args.reference
    )

    policy = yaml.safe_load(
        args.tolerances.read_text(
            encoding="utf-8"
        )
    )

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

    keys = [
        value.strip()
        for value in args.keys.split(",")
        if value.strip()
    ]

    required = set(
        keys
        + [
            args.q_column,
            args.support_column,
        ]
    )

    missing_observed = (
        required
        - set(observed.columns)
    )

    missing_reference = (
        required
        - set(reference.columns)
    )

    if missing_observed:
        raise KeyError(
            "Observed table missing columns: "
            f"{sorted(missing_observed)}"
        )

    if missing_reference:
        raise KeyError(
            "Reference table missing columns: "
            f"{sorted(missing_reference)}"
        )

    observed = observed[
        list(required)
    ].copy()

    reference = reference[
        list(required)
    ].copy()

    merged = observed.merge(
        reference,
        on=keys,
        how="outer",
        suffixes=(
            "_observed",
            "_reference",
        ),
        indicator=True,
        validate="one_to_one",
    )

    pair_complete = bool(
        merged["_merge"]
        .eq("both")
        .all()
    )

    observed_q = merged[
        f"{args.q_column}_observed"
    ].to_numpy(
        dtype=float
    )

    reference_q = merged[
        f"{args.q_column}_reference"
    ].to_numpy(
        dtype=float
    )

    q_match = (
        pair_complete
        and np.allclose(
            observed_q,
            reference_q,
            atol=atol,
            rtol=rtol,
            equal_nan=True,
        )
    )

    if pair_complete:
        max_abs_difference = float(
            np.nanmax(
                np.abs(
                    observed_q
                    - reference_q
                )
            )
        )

        observed_supported = (
            parse_boolean_series(
                merged[
                    f"{args.support_column}_observed"
                ]
            )
        )

        reference_supported = (
            parse_boolean_series(
                merged[
                    f"{args.support_column}_reference"
                ]
            )
        )

        support_match = bool(
            np.array_equal(
                observed_supported.to_numpy(),
                reference_supported.to_numpy(),
            )
        )

        support_mismatch_count = int(
            (
                observed_supported
                != reference_supported
            ).sum()
        )

    else:
        max_abs_difference = None
        support_match = False
        support_mismatch_count = None

    passed = bool(
        pair_complete
        and q_match
        and support_match
    )

    report = {
        "verifier_version": (
            VERIFIER_VERSION
        ),
        "generated_at_utc": (
            datetime.now(
                timezone.utc
            ).isoformat()
        ),
        "passed": passed,
        "checks": {
            "same_key_family": (
                pair_complete
            ),
            "q_values_match": bool(
                q_match
            ),
            "support_classification_match": (
                support_match
            ),
        },
        "comparison": {
            "key_columns": keys,
            "q_column": args.q_column,
            "support_column": (
                args.support_column
            ),
            "absolute_tolerance": (
                atol
            ),
            "relative_tolerance": (
                rtol
            ),
            "max_absolute_q_difference": (
                max_abs_difference
            ),
            "support_mismatch_count": (
                support_mismatch_count
            ),
            "observed_rows": int(
                observed.shape[0]
            ),
            "reference_rows": int(
                reference.shape[0]
            ),
        },
    }

    args.output.write_text(
        json.dumps(
            report,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        f"Same key family: {pair_complete}"
    )

    print(
        f"Q-values match: {q_match}"
    )

    print(
        "Scientific classification match: "
        f"{support_match}"
    )

    print(
        f"Verification passed: {passed}"
    )

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
