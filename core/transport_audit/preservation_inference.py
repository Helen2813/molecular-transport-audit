"""Permutation inference for molecular-program preservation."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from transport_audit.multiplicity import (
    benjamini_hochberg,
)
from transport_audit.preservation import (
    clean_gene_symbol,
    correlation_matrix,
    edge_spearman_from_matrices,
    fit_oriented_external_pc1,
    loading_preservation,
    standardize_expression,
    upper_triangle,
)


DEFAULT_N_PERMUTATIONS = 5000
DEFAULT_RANDOM_SEED = 42
DEFAULT_MODULE_SEED_STRIDE = 100
DEFAULT_MINIMUM_GENES = 3
DEFAULT_PCA_RANDOM_STATE = 42


@dataclass(frozen=True)
class PermutationTestResult:
    """Result of one two-sided permutation test."""

    p_value: float
    extreme_count: int
    n_permutations: int
    seed: int


def _validate_n_permutations(
    n_permutations: int,
) -> None:
    if n_permutations < 1:
        raise ValueError(
            "n_permutations must be at least 1."
        )


def empirical_two_sided_p(
    extreme_count: int,
    n_permutations: int,
) -> float:
    """Apply the frozen +1 Monte-Carlo correction."""
    _validate_n_permutations(
        n_permutations
    )

    if (
        extreme_count < 0
        or extreme_count
        > n_permutations
    ):
        raise ValueError(
            "extreme_count must lie between "
            "0 and n_permutations."
        )

    return float(
        (extreme_count + 1)
        / (n_permutations + 1)
    )


def permutation_seeds(
    module_index: int,
    *,
    base_seed: int = DEFAULT_RANDOM_SEED,
    module_seed_stride: int = DEFAULT_MODULE_SEED_STRIDE,
) -> tuple[int, int]:
    """Return frozen edge and loading seeds for one module."""
    if module_index < 0:
        raise ValueError(
            "module_index must be non-negative."
        )

    module_seed = (
        int(base_seed)
        + int(module_index)
        * int(module_seed_stride)
    )

    return (
        module_seed,
        module_seed + 1,
    )


def edge_gene_label_permutation(
    reference_edges: np.ndarray,
    external_correlation: np.ndarray,
    observed: float,
    *,
    n_permutations: int = DEFAULT_N_PERMUTATIONS,
    seed: int,
) -> PermutationTestResult:
    """Run the frozen edge gene-label permutation test."""
    _validate_n_permutations(
        n_permutations
    )

    reference = np.asarray(
        reference_edges,
        dtype=float,
    )

    external = np.asarray(
        external_correlation,
        dtype=float,
    )

    if (
        external.ndim != 2
        or external.shape[0]
        != external.shape[1]
    ):
        raise ValueError(
            "external_correlation must be square."
        )

    n_genes = int(
        external.shape[0]
    )

    expected_edges = (
        n_genes
        * (n_genes - 1)
        // 2
    )

    if reference.size != expected_edges:
        raise ValueError(
            "reference_edges has an "
            "unexpected length."
        )

    if not np.isfinite(
        observed
    ):
        return PermutationTestResult(
            p_value=np.nan,
            extreme_count=0,
            n_permutations=n_permutations,
            seed=int(seed),
        )

    rng = np.random.default_rng(
        seed
    )

    extreme_count = 0

    for _ in range(
        n_permutations
    ):
        permutation = rng.permutation(
            n_genes
        )

        permuted = external[
            np.ix_(
                permutation,
                permutation,
            )
        ]

        value = stats.spearmanr(
            reference,
            upper_triangle(
                permuted
            ),
        ).statistic

        if (
            np.isfinite(value)
            and abs(value)
            >= abs(observed)
        ):
            extreme_count += 1

    return PermutationTestResult(
        p_value=empirical_two_sided_p(
            extreme_count,
            n_permutations,
        ),
        extreme_count=extreme_count,
        n_permutations=n_permutations,
        seed=int(seed),
    )


def loading_gene_label_permutation(
    frozen_loadings: np.ndarray,
    external_loadings: np.ndarray,
    observed: float,
    *,
    n_permutations: int = DEFAULT_N_PERMUTATIONS,
    seed: int,
) -> PermutationTestResult:
    """Run the frozen loading-label permutation test."""
    _validate_n_permutations(
        n_permutations
    )

    frozen = np.asarray(
        frozen_loadings,
        dtype=float,
    )

    external = np.asarray(
        external_loadings,
        dtype=float,
    )

    if frozen.ndim != 1:
        raise ValueError(
            "frozen_loadings must be one-dimensional."
        )

    if external.ndim != 1:
        raise ValueError(
            "external_loadings must be one-dimensional."
        )

    if frozen.shape != external.shape:
        raise ValueError(
            "Frozen and external loading vectors "
            "must have equal length."
        )

    if not np.isfinite(
        observed
    ):
        return PermutationTestResult(
            p_value=np.nan,
            extreme_count=0,
            n_permutations=n_permutations,
            seed=int(seed),
        )

    rng = np.random.default_rng(
        seed
    )

    extreme_count = 0

    for _ in range(
        n_permutations
    ):
        permuted = rng.permutation(
            external
        )

        value = stats.spearmanr(
            frozen,
            permuted,
        ).statistic

        if (
            np.isfinite(value)
            and abs(value)
            >= abs(observed)
        ):
            extreme_count += 1

    return PermutationTestResult(
        p_value=empirical_two_sided_p(
            extreme_count,
            n_permutations,
        ),
        extreme_count=extreme_count,
        n_permutations=n_permutations,
        seed=int(seed),
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
            "Frozen weights are missing columns: "
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
            f"No frozen genes found for "
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


def adjust_direct_permutation_family(
    table: pd.DataFrame,
) -> pd.DataFrame:
    """Apply BH jointly across edge and loading tests."""
    required = {
        "edge_permutation_p",
        "loading_permutation_p",
    }

    missing = sorted(
        required.difference(
            table.columns
        )
    )

    if missing:
        raise ValueError(
            "Permutation table is missing: "
            + ", ".join(
                missing
            )
        )

    result = table.copy()

    edge_p = pd.to_numeric(
        result[
            "edge_permutation_p"
        ],
        errors="coerce",
    ).to_numpy(
        dtype=float
    )

    loading_p = pd.to_numeric(
        result[
            "loading_permutation_p"
        ],
        errors="coerce",
    ).to_numpy(
        dtype=float
    )

    family = np.concatenate(
        [
            edge_p,
            loading_p,
        ]
    )

    q_values = (
        benjamini_hochberg(
            family
        )
    )

    n_modules = int(
        result.shape[0]
    )

    result[
        "edge_q_bh_8"
    ] = q_values[
        :n_modules
    ]

    result[
        "loading_q_bh_8"
    ] = q_values[
        n_modules:
    ]

    result[
        "direct_permutation_family_size"
    ] = int(
        np.isfinite(
            family
        ).sum()
    )

    return result


def analyze_permutation_preservation(
    reference_expression: pd.DataFrame,
    external_expression: pd.DataFrame,
    weights: pd.DataFrame,
    *,
    module_order: Sequence[str],
    minimum_genes: int = DEFAULT_MINIMUM_GENES,
    n_permutations: int = DEFAULT_N_PERMUTATIONS,
    base_seed: int = DEFAULT_RANDOM_SEED,
    module_seed_stride: int = DEFAULT_MODULE_SEED_STRIDE,
    pca_random_state: int = DEFAULT_PCA_RANDOM_STATE,
) -> pd.DataFrame:
    """Run frozen edge/loading permutation inference."""
    _validate_n_permutations(
        n_permutations
    )

    reference_z = (
        standardize_expression(
            reference_expression
        )
    )

    external_z = (
        standardize_expression(
            external_expression
        )
    )

    rows: list[
        dict[str, Any]
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
            if (
                gene
                in reference_z.columns
                and gene
                in external_z.columns
            )
        ]

        (
            edge_seed,
            loading_seed,
        ) = permutation_seeds(
            module_index,
            base_seed=base_seed,
            module_seed_stride=(
                module_seed_stride
            ),
        )

        if len(genes) < minimum_genes:
            rows.append(
                {
                    "module_label":
                        str(module),
                    "n_common_genes":
                        int(
                            len(genes)
                        ),
                    "edge_spearman":
                        np.nan,
                    "edge_permutation_p":
                        np.nan,
                    "edge_extreme_count":
                        0,
                    "edge_seed":
                        edge_seed,
                    "loading_spearman":
                        np.nan,
                    "loading_permutation_p":
                        np.nan,
                    "loading_extreme_count":
                        0,
                    "loading_seed":
                        loading_seed,
                    (
                        "external_pc1_"
                        "variance_explained"
                    ):
                        np.nan,
                    (
                        "pc1_orientation_"
                        "correlation_with_"
                        "frozen_score"
                    ):
                        np.nan,
                    "n_permutations":
                        n_permutations,
                    "estimable":
                        False,
                    "nonestimable_reason":
                        (
                            "fewer_than_"
                            f"{minimum_genes}_"
                            "common_genes"
                        ),
                }
            )

            continue

        reference_correlation = (
            correlation_matrix(
                reference_z,
                genes,
            )
        )

        external_correlation = (
            correlation_matrix(
                external_z,
                genes,
            )
        )

        reference_edges = (
            upper_triangle(
                reference_correlation
            )
        )

        edge_spearman = (
            edge_spearman_from_matrices(
                reference_correlation,
                external_correlation,
            )
        )

        edge_test = (
            edge_gene_label_permutation(
                reference_edges=(
                    reference_edges
                ),
                external_correlation=(
                    external_correlation
                ),
                observed=(
                    edge_spearman
                ),
                n_permutations=(
                    n_permutations
                ),
                seed=edge_seed,
            )
        )

        frozen_loadings = (
            part.loc[
                genes,
                "risk_oriented_loading",
            ]
        )

        external_pc1 = (
            fit_oriented_external_pc1(
                external_z=external_z,
                genes=genes,
                frozen_loadings=(
                    frozen_loadings
                ),
                random_state=(
                    pca_random_state
                ),
            )
        )

        loading_spearman = (
            loading_preservation(
                frozen_loadings=(
                    frozen_loadings
                ),
                external_pc1_loadings=(
                    external_pc1.loadings
                ),
                genes=genes,
            )
        )

        loading_test = (
            loading_gene_label_permutation(
                frozen_loadings=(
                    frozen_loadings
                    .to_numpy(
                        dtype=float
                    )
                ),
                external_loadings=(
                    external_pc1.loadings
                ),
                observed=(
                    loading_spearman
                ),
                n_permutations=(
                    n_permutations
                ),
                seed=loading_seed,
            )
        )

        rows.append(
            {
                "module_label":
                    str(module),
                "n_common_genes":
                    int(
                        len(genes)
                    ),
                "edge_spearman":
                    float(
                        edge_spearman
                    ),
                "edge_permutation_p":
                    float(
                        edge_test
                        .p_value
                    ),
                "edge_extreme_count":
                    int(
                        edge_test
                        .extreme_count
                    ),
                "edge_seed":
                    int(
                        edge_seed
                    ),
                "loading_spearman":
                    float(
                        loading_spearman
                    ),
                "loading_permutation_p":
                    float(
                        loading_test
                        .p_value
                    ),
                "loading_extreme_count":
                    int(
                        loading_test
                        .extreme_count
                    ),
                "loading_seed":
                    int(
                        loading_seed
                    ),
                (
                    "external_pc1_"
                    "variance_explained"
                ):
                    float(
                        external_pc1
                        .explained_variance_ratio
                    ),
                (
                    "pc1_orientation_"
                    "correlation_with_"
                    "frozen_score"
                ):
                    float(
                        external_pc1
                        .orientation_correlation
                    ),
                "n_permutations":
                    int(
                        n_permutations
                    ),
                "estimable":
                    True,
                "nonestimable_reason":
                    "",
            }
        )

    result = pd.DataFrame(
        rows
    )

    return (
        adjust_direct_permutation_family(
            result
        )
    )
