"""Outcome-blind random-panel specificity controls."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from transport_audit.multiplicity import (
    benjamini_hochberg,
)
from transport_audit.preservation import (
    clean_gene_symbol,
    normalize_expression_columns,
    standardize_expression,
    upper_triangle,
)


LEGACY_Z_VARIANCE = "legacy_z_variance"

PRESTANDARDIZATION_VARIANCE = (
    "prestandardization_variance"
)

VALID_MATCHING_MODES = {
    LEGACY_Z_VARIANCE,
    PRESTANDARDIZATION_VARIANCE,
}

DEFAULT_N_RANDOM_PANELS = 1000
DEFAULT_N_VARIABILITY_BINS = 10
DEFAULT_RANDOM_SEED = 42
DEFAULT_MINIMUM_GENES = 3


def prepare_prestandardization_expression(
    expression: pd.DataFrame,
) -> pd.DataFrame:
    """Prepare expression before gene-wise z scaling."""
    prepared = normalize_expression_columns(
        expression
    )

    numeric = prepared.apply(
        pd.to_numeric,
        errors="coerce",
    )

    numeric = numeric.fillna(
        numeric.median(axis=0)
    )

    return numeric


def variability_bins(
    reference_expression: pd.DataFrame,
    external_expression: pd.DataFrame,
    *,
    matching_mode: str,
    n_bins: int = DEFAULT_N_VARIABILITY_BINS,
) -> pd.Series:
    """Build legacy or corrected outcome-blind variability bins."""
    if matching_mode not in (
        VALID_MATCHING_MODES
    ):
        raise ValueError(
            "Unknown matching_mode: "
            f"{matching_mode}"
        )

    if n_bins < 1:
        raise ValueError(
            "n_bins must be at least 1."
        )

    reference_pre = (
        prepare_prestandardization_expression(
            reference_expression
        )
    )

    external_pre = (
        prepare_prestandardization_expression(
            external_expression
        )
    )

    reference_z = (
        standardize_expression(
            reference_pre
        )
    )

    external_z = (
        standardize_expression(
            external_pre
        )
    )

    common = (
        reference_z.columns
        .intersection(
            external_z.columns
        )
    )

    if matching_mode == LEGACY_Z_VARIANCE:
        reference_basis = (
            reference_z[
                common
            ]
        )

        external_basis = (
            external_z[
                common
            ]
        )
    else:
        reference_basis = (
            reference_pre[
                common
            ]
        )

        external_basis = (
            external_pre[
                common
            ]
        )

    reference_variance = (
        reference_basis.var(
            axis=0,
            ddof=1,
        )
    )

    external_variance = (
        external_basis.var(
            axis=0,
            ddof=1,
        )
    )

    combined_rank = (
        reference_variance.rank(
            pct=True
        )
        + external_variance.rank(
            pct=True
        )
    ) / 2.0

    return pd.qcut(
        combined_rank.rank(
            method="average"
        ),
        q=min(
            n_bins,
            len(
                combined_rank
            ),
        ),
        labels=False,
        duplicates="drop",
    )


def random_panel_empirical_p(
    null_values: np.ndarray,
    observed: float,
) -> float:
    """Frozen two-sided empirical p with +1 correction."""
    null_array = np.asarray(
        null_values,
        dtype=float,
    )

    null_array = null_array[
        np.isfinite(
            null_array
        )
    ]

    if (
        not np.isfinite(
            observed
        )
        or null_array.size == 0
    ):
        return np.nan

    count = int(
        np.sum(
            np.abs(
                null_array
            )
            >= abs(
                observed
            )
        )
    )

    return float(
        (count + 1)
        / (
            null_array.size
            + 1
        )
    )


def draw_matched_panel(
    candidate_genes: Sequence[str],
    bins: pd.Series,
    target_bins: Sequence[Any],
    rng: np.random.Generator,
) -> tuple[
    list[str],
    int,
]:
    """
    Draw one panel using the frozen without-replacement rule.

    A same-bin candidate is preferred. If none remains, the
    legacy fallback samples from all unused candidates.
    """
    candidates = [
        str(gene)
        for gene in candidate_genes
    ]

    available_by_bin: dict[
        Any,
        list[str],
    ] = {}

    for gene in candidates:
        bin_value = bins.get(
            gene,
            np.nan,
        )

        if pd.isna(
            bin_value
        ):
            continue

        available_by_bin.setdefault(
            bin_value,
            [],
        ).append(
            gene
        )

    selected: list[str] = []
    used: set[str] = set()

    fallback_count = 0

    for target_bin in target_bins:
        pool = (
            available_by_bin.get(
                target_bin,
                [],
            )
        )

        if pool:
            gene = str(
                rng.choice(
                    pool
                )
            )

            pool.remove(
                gene
            )
        else:
            fallback_pool = [
                gene
                for gene in candidates
                if gene not in used
            ]

            if not fallback_pool:
                return (
                    [],
                    fallback_count,
                )

            gene = str(
                rng.choice(
                    fallback_pool
                )
            )

            fallback_count += 1

            gene_bin = bins.get(
                gene,
                np.nan,
            )

            if (
                not pd.isna(
                    gene_bin
                )
                and gene_bin
                in available_by_bin
                and gene
                in available_by_bin[
                    gene_bin
                ]
            ):
                available_by_bin[
                    gene_bin
                ].remove(
                    gene
                )

        selected.append(
            gene
        )

        used.add(
            gene
        )

    return (
        selected,
        fallback_count,
    )


def _clean_weights(
    weights: pd.DataFrame,
) -> pd.DataFrame:
    required = {
        "module_label",
        "canine_gene_symbol",
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

    result = weights.copy()

    result[
        "canine_gene_symbol"
    ] = (
        result[
            "canine_gene_symbol"
        ]
        .map(
            clean_gene_symbol
        )
    )

    result = result[
        result[
            "canine_gene_symbol"
        ].ne("")
    ].copy()

    return result


def random_panel_controls(
    reference_expression: pd.DataFrame,
    external_expression: pd.DataFrame,
    weights: pd.DataFrame,
    observed_structure: pd.DataFrame,
    *,
    module_order: Sequence[str],
    matching_mode: str,
    n_random_panels: int = DEFAULT_N_RANDOM_PANELS,
    n_variability_bins: int = DEFAULT_N_VARIABILITY_BINS,
    random_seed: int = DEFAULT_RANDOM_SEED,
    minimum_genes: int = DEFAULT_MINIMUM_GENES,
) -> pd.DataFrame:
    """Run legacy-compatible or corrected random-panel controls."""
    if matching_mode not in (
        VALID_MATCHING_MODES
    ):
        raise ValueError(
            "Unknown matching_mode: "
            f"{matching_mode}"
        )

    if n_random_panels < 1:
        raise ValueError(
            "n_random_panels must be at least 1."
        )

    reference_pre = (
        prepare_prestandardization_expression(
            reference_expression
        )
    )

    external_pre = (
        prepare_prestandardization_expression(
            external_expression
        )
    )

    reference_z = (
        standardize_expression(
            reference_pre
        )
    )

    external_z = (
        standardize_expression(
            external_pre
        )
    )

    common = (
        reference_z.columns
        .intersection(
            external_z.columns
        )
    )

    bins = variability_bins(
        reference_pre[
            common
        ],
        external_pre[
            common
        ],
        matching_mode=(
            matching_mode
        ),
        n_bins=(
            n_variability_bins
        ),
    )

    cleaned_weights = (
        _clean_weights(
            weights
        )
    )

    frozen_union = set(
        cleaned_weights.loc[
            cleaned_weights[
                "module_label"
            ]
            .astype(str)
            .isin(
                [
                    str(module)
                    for module
                    in module_order
                ]
            ),
            "canine_gene_symbol",
        ]
    )

    candidate_genes = [
        str(gene)
        for gene in common
        if gene
        not in frozen_union
    ]

    observed = (
        observed_structure.copy()
    )

    observed[
        "module_label"
    ] = (
        observed[
            "module_label"
        ]
        .astype(str)
    )

    rng = np.random.default_rng(
        random_seed
    )

    rows: list[
        dict[str, Any]
    ] = []

    for module in module_order:
        module = str(
            module
        )

        part = cleaned_weights[
            cleaned_weights[
                "module_label"
            ]
            .astype(str)
            .eq(
                module
            )
        ].copy()

        part = part.drop_duplicates(
            "canine_gene_symbol",
            keep="first",
        )

        target_genes = [
            gene
            for gene
            in part[
                "canine_gene_symbol"
            ]
            if gene in common
        ]

        if (
            len(target_genes)
            < minimum_genes
        ):
            continue

        target_bins = (
            bins.reindex(
                target_genes
            )
            .to_numpy()
        )

        null_values: list[
            float
        ] = []

        total_fallbacks = 0

        for _ in range(
            n_random_panels
        ):
            (
                selected,
                fallback_count,
            ) = draw_matched_panel(
                candidate_genes=(
                    candidate_genes
                ),
                bins=bins,
                target_bins=(
                    target_bins
                ),
                rng=rng,
            )

            total_fallbacks += (
                fallback_count
            )

            if (
                len(selected)
                != len(
                    target_genes
                )
            ):
                continue

            reference_correlation = (
                np.corrcoef(
                    reference_z[
                        selected
                    ].to_numpy(
                        dtype=float
                    ),
                    rowvar=False,
                )
            )

            external_correlation = (
                np.corrcoef(
                    external_z[
                        selected
                    ].to_numpy(
                        dtype=float
                    ),
                    rowvar=False,
                )
            )

            edge = stats.spearmanr(
                upper_triangle(
                    reference_correlation
                ),
                upper_triangle(
                    external_correlation
                ),
            ).statistic

            if np.isfinite(
                edge
            ):
                null_values.append(
                    float(edge)
                )

        module_observed = observed[
            observed[
                "module_label"
            ].eq(
                module
            )
        ]

        if (
            module_observed.shape[0]
            != 1
        ):
            raise ValueError(
                "Observed structure must "
                "contain exactly one row "
                f"for {module}."
            )

        observed_edge = float(
            module_observed[
                "edge_spearman"
            ].iloc[0]
        )

        null_array = np.asarray(
            null_values,
            dtype=float,
        )

        rows.append(
            {
                "module_label":
                    module,
                "matching_mode":
                    matching_mode,
                "n_module_genes":
                    int(
                        len(
                            target_genes
                        )
                    ),
                "n_random_panels":
                    int(
                        len(
                            null_array
                        )
                    ),
                "observed_edge_spearman":
                    observed_edge,
                "random_edge_median":
                    (
                        float(
                            np.median(
                                null_array
                            )
                        )
                        if len(
                            null_array
                        )
                        else np.nan
                    ),
                "random_edge_q05":
                    (
                        float(
                            np.quantile(
                                null_array,
                                0.05,
                            )
                        )
                        if len(
                            null_array
                        )
                        else np.nan
                    ),
                "random_edge_q95":
                    (
                        float(
                            np.quantile(
                                null_array,
                                0.95,
                            )
                        )
                        if len(
                            null_array
                        )
                        else np.nan
                    ),
                "random_panel_empirical_p":
                    random_panel_empirical_p(
                        null_array,
                        observed_edge,
                    ),
                "fallback_draws":
                    int(
                        total_fallbacks
                    ),
                "candidate_genes":
                    int(
                        len(
                            candidate_genes
                        )
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


def classify_external_representation(
    structure: pd.DataFrame,
    random_controls: pd.DataFrame,
) -> pd.DataFrame:
    """Reproduce the frozen script-46 classification rule."""
    result = structure.copy()

    result[
        "module_label"
    ] = (
        result[
            "module_label"
        ]
        .astype(str)
    )

    random_table = (
        random_controls[
            [
                "module_label",
                "random_panel_empirical_p",
            ]
        ]
        .copy()
    )

    random_table[
        "module_label"
    ] = (
        random_table[
            "module_label"
        ]
        .astype(str)
    )

    result = result.merge(
        random_table,
        on="module_label",
        how="left",
    )

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

    classes: list[str] = []

    for row in result.itertuples(
        index=False
    ):
        edge_supported = bool(
            np.isfinite(
                row.edge_q_bh_8
            )
            and row.edge_q_bh_8
            < 0.05
        )

        loading_supported = bool(
            np.isfinite(
                row.loading_q_bh_8
            )
            and row.loading_q_bh_8
            < 0.05
        )

        random_supported = bool(
            np.isfinite(
                row.random_panel_empirical_p
            )
            and row.random_panel_empirical_p
            < 0.05
        )

        reliable = bool(
            np.isfinite(
                row.split_half_median
            )
            and row.split_half_median
            >= 0.60
        )

        if (
            edge_supported
            and loading_supported
            and reliable
        ):
            label = (
                "strong_external_canine_"
                "representation_preservation"
            )

        elif (
            (
                edge_supported
                or loading_supported
            )
            and reliable
        ):
            label = (
                "partial_external_canine_"
                "representation_preservation"
            )

        elif (
            random_supported
            or edge_supported
            or loading_supported
        ):
            label = (
                "limited_specific_"
                "external_canine_signal"
            )

        else:
            label = (
                "no_clear_external_canine_"
                "representation_preservation"
            )

        classes.append(
            label
        )

    result[
        "external_canine_representation_class"
    ] = classes

    return result
