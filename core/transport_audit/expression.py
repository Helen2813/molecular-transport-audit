"""Expression preprocessing for frozen-program scoring."""

from __future__ import annotations

import numpy as np
import pandas as pd


def normalize_expression_columns(
    expression: pd.DataFrame,
) -> pd.DataFrame:
    """Upper-case genes and keep first duplicates."""
    result = expression.copy()

    result.columns = (
        result.columns
        .astype(str)
        .str.strip()
        .str.upper()
    )

    result = result.loc[
        :,
        ~pd.Index(
            result.columns
        ).duplicated(
            keep="first"
        ),
    ].copy()

    return result


def standardize_genes(
    expression: pd.DataFrame,
    genes: list[str] | tuple[str, ...],
) -> pd.DataFrame:
    """
    Standardize genes within the target cohort.

    Missing values are median-imputed per gene.
    Sample standard deviation uses ddof=1.
    """
    selected = list(genes)

    if not selected:
        return pd.DataFrame(
            index=expression.index
        )

    x = expression.loc[
        :,
        selected,
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )

    x = x.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    x = x.fillna(
        x.median(axis=0)
    )

    means = x.mean(axis=0)

    stds = (
        x.std(
            axis=0,
            ddof=1,
        )
        .replace(
            0.0,
            np.nan,
        )
    )

    z = (x - means) / stds

    valid = z.columns[
        z.notna().all(axis=0)
    ]

    return z.loc[
        :,
        valid,
    ].copy()


def zscore_series(
    values: pd.Series,
) -> pd.Series:
    """Standardize a score using sample SD."""
    numeric = pd.to_numeric(
        values,
        errors="coerce",
    )

    std = numeric.std(
        ddof=1
    )

    if (
        not np.isfinite(std)
        or std == 0
    ):
        return pd.Series(
            np.nan,
            index=numeric.index,
            dtype=float,
        )

    return (
        numeric
        - numeric.mean()
    ) / std
