"""Reliability diagnostics for frozen molecular programs."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from transport_audit.preservation import (
    clean_gene_symbol,
    frozen_signed_score,
    standardize_expression,
)


DEFAULT_N_SPLIT_HALF_REPEATS = 2000
DEFAULT_RANDOM_SEED = 42
DEFAULT_MODULE_SEED_STRIDE = 100
DEFAULT_SPLIT_HALF_SEED_OFFSET = 2
DEFAULT_MINIMUM_GENES = 3


def reliability_seed(
    module_index: int,
    *,
    base_seed: int = DEFAULT_RANDOM_SEED,
    module_seed_stride: int = DEFAULT_MODULE_SEED_STRIDE,
    seed_offset: int = DEFAULT_SPLIT_HALF_SEED_OFFSET,
) -> int:
    """Return the frozen split-half seed for one module."""
    if module_index < 0:
        raise ValueError(
            "module_index must be non-negative."
        )

    return int(
        base_seed
        + module_index
        * module_seed_stride
        + seed_offset
    )


def split_half_reliability(
    expression_z: pd.DataFrame,
    genes: Sequence[str],
    loadings: pd.Series,
    *,
    seed: int,
    n_repeats: int = DEFAULT_N_SPLIT_HALF_REPEATS,
) -> dict[str, float | int]:
    """Reproduce the frozen non-overlapping split-half reliability."""
    selected = [
        str(gene)
        for gene in genes
    ]

    if n_repeats < 1:
        raise ValueError(
            "n_repeats must be at least 1."
        )

    if len(selected) < 4:
        return {
            "split_half_median":
                np.nan,
            "split_half_q05":
                np.nan,
            "split_half_q95":
                np.nan,
            "split_half_valid_repeats":
                0,
        }

    rng = np.random.default_rng(
        seed
    )

    correlations: list[
        float
    ] = []

    for _ in range(
        n_repeats
    ):
        shuffled = np.asarray(
            selected,
            dtype=object,
        )

        shuffled = rng.permutation(
            shuffled
        )

        midpoint = (
            len(shuffled)
            // 2
        )

        first = (
            shuffled[
                :midpoint
            ].tolist()
        )

        second = (
            shuffled[
                midpoint:
            ].tolist()
        )

        first_score = (
            frozen_signed_score(
                expression_z,
                first,
                loadings,
            )
        )

        second_score = (
            frozen_signed_score(
                expression_z,
                second,
                loadings,
            )
        )

        correlation = (
            first_score.corr(
                second_score
            )
        )

        if np.isfinite(
            correlation
        ):
            correlations.append(
                float(
                    correlation
                )
            )

    if not correlations:
        return {
            "split_half_median":
                np.nan,
            "split_half_q05":
                np.nan,
            "split_half_q95":
                np.nan,
            "split_half_valid_repeats":
                0,
        }

    values = np.asarray(
        correlations,
        dtype=float,
    )

    return {
        "split_half_median":
            float(
                np.median(
                    values
                )
            ),
        "split_half_q05":
            float(
                np.quantile(
                    values,
                    0.05,
                )
            ),
        "split_half_q95":
            float(
                np.quantile(
                    values,
                    0.95,
                )
            ),
        "split_half_valid_repeats":
            int(
                len(values)
            ),
    }


def gene_leave_one_out(
    expression_z: pd.DataFrame,
    module_label: str,
    genes: Sequence[str],
    loadings: pd.Series,
    *,
    minimum_genes: int = DEFAULT_MINIMUM_GENES,
) -> pd.DataFrame:
    """Reproduce frozen single-gene leave-one-out stability."""
    selected = [
        str(gene)
        for gene in genes
    ]

    if len(selected) <= minimum_genes:
        return pd.DataFrame(
            columns=[
                "module_label",
                "left_out_gene",
                "n_genes_remaining",
                "correlation_with_full_score",
            ]
        )

    full_score = frozen_signed_score(
        expression_z,
        selected,
        loadings,
    )

    rows: list[
        dict[str, Any]
    ] = []

    for gene in selected:
        subset = [
            item
            for item in selected
            if item != gene
        ]

        score = frozen_signed_score(
            expression_z,
            subset,
            loadings,
        )

        rows.append(
            {
                "module_label":
                    str(
                        module_label
                    ),
                "left_out_gene":
                    gene,
                "n_genes_remaining":
                    int(
                        len(subset)
                    ),
                "correlation_with_full_score":
                    float(
                        full_score.corr(
                            score
                        )
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


def _module_weights(
    weights: pd.DataFrame,
    module_label: str,
) -> pd.DataFrame:
    required = {
        "module_label",
        "canine_gene_symbol",
        "risk_oriented_loading",
    }

    missing = sorted(
        required.difference(
            weights.columns
        )
    )

    if missing:
        raise ValueError(
            "Frozen weights are missing: "
            + ", ".join(
                missing
            )
        )

    part = weights[
        weights[
            "module_label"
        ]
        .astype(str)
        .eq(
            str(module_label)
        )
    ].copy()

    if part.empty:
        raise ValueError(
            f"No frozen genes for "
            f"{module_label}."
        )

    part[
        "canine_gene_symbol"
    ] = (
        part[
            "canine_gene_symbol"
        ]
        .map(
            clean_gene_symbol
        )
    )

    part = part[
        part[
            "canine_gene_symbol"
        ].ne("")
    ].copy()

    part = (
        part.drop_duplicates(
            "canine_gene_symbol",
            keep="first",
        )
        .set_index(
            "canine_gene_symbol"
        )
    )

    part[
        "risk_oriented_loading"
    ] = pd.to_numeric(
        part[
            "risk_oriented_loading"
        ],
        errors="coerce",
    ).fillna(0.0)

    return part


def analyze_reliability(
    external_expression: pd.DataFrame,
    weights: pd.DataFrame,
    *,
    module_order: Sequence[str],
    minimum_genes: int = DEFAULT_MINIMUM_GENES,
    n_split_half_repeats: int = DEFAULT_N_SPLIT_HALF_REPEATS,
    base_seed: int = DEFAULT_RANDOM_SEED,
    module_seed_stride: int = DEFAULT_MODULE_SEED_STRIDE,
    split_half_seed_offset: int = DEFAULT_SPLIT_HALF_SEED_OFFSET,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
]:
    """Run split-half and per-gene LOO reliability for frozen programs."""
    external_z = standardize_expression(
        external_expression
    )

    reliability_rows: list[
        dict[str, Any]
    ] = []

    loo_tables: list[
        pd.DataFrame
    ] = []

    for module_index, module in enumerate(
        module_order
    ):
        part = _module_weights(
            weights,
            str(module),
        )

        frozen_genes = [
            str(gene)
            for gene in part.index
        ]

        genes = [
            gene
            for gene in frozen_genes
            if gene
            in external_z.columns
        ]

        seed = reliability_seed(
            module_index,
            base_seed=base_seed,
            module_seed_stride=(
                module_seed_stride
            ),
            seed_offset=(
                split_half_seed_offset
            ),
        )

        if len(genes) < minimum_genes:
            reliability_rows.append(
                {
                    "module_label":
                        str(module),
                    "n_common_genes":
                        int(
                            len(genes)
                        ),
                    "split_half_median":
                        np.nan,
                    "split_half_q05":
                        np.nan,
                    "split_half_q95":
                        np.nan,
                    "split_half_valid_repeats":
                        0,
                    "minimum_gene_loo_correlation":
                        np.nan,
                    "median_gene_loo_correlation":
                        np.nan,
                    "split_half_seed":
                        int(seed),
                }
            )

            continue

        loadings = (
            part.loc[
                genes,
                "risk_oriented_loading",
            ]
        )

        reliability = (
            split_half_reliability(
                expression_z=external_z,
                genes=genes,
                loadings=loadings,
                seed=seed,
                n_repeats=(
                    n_split_half_repeats
                ),
            )
        )

        loo = gene_leave_one_out(
            expression_z=external_z,
            module_label=str(
                module
            ),
            genes=genes,
            loadings=loadings,
            minimum_genes=(
                minimum_genes
            ),
        )

        if not loo.empty:
            loo_tables.append(
                loo
            )

        reliability_rows.append(
            {
                "module_label":
                    str(module),
                "n_common_genes":
                    int(
                        len(genes)
                    ),
                **reliability,
                "minimum_gene_loo_correlation":
                    (
                        float(
                            loo[
                                "correlation_with_full_score"
                            ].min()
                        )
                        if not loo.empty
                        else np.nan
                    ),
                "median_gene_loo_correlation":
                    (
                        float(
                            loo[
                                "correlation_with_full_score"
                            ].median()
                        )
                        if not loo.empty
                        else np.nan
                    ),
                "split_half_seed":
                    int(seed),
            }
        )

    reliability_table = (
        pd.DataFrame(
            reliability_rows
        )
    )

    gene_loo_table = (
        pd.concat(
            loo_tables,
            ignore_index=True,
        )
        if loo_tables
        else pd.DataFrame(
            columns=[
                "module_label",
                "left_out_gene",
                "n_genes_remaining",
                "correlation_with_full_score",
            ]
        )
    )

    return (
        reliability_table,
        gene_loo_table,
    )
