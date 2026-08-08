"""Molecular-program representation preservation."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA


DEFAULT_MINIMUM_GENES = 3
DEFAULT_PCA_RANDOM_STATE = 42


@dataclass(frozen=True)
class ExternalPC1:
    """Outcome-blind external PC1 oriented to a frozen score."""

    scores: pd.Series
    loadings: np.ndarray
    explained_variance_ratio: float
    orientation_correlation: float


def clean_gene_symbol(
    value: Any,
) -> str:
    """Normalize a canine gene symbol using the frozen legacy rule."""
    text = str(value).strip().upper()

    if text in {
        "",
        "NAN",
        "NONE",
        "NA",
    }:
        return ""

    if "|" in text:
        parts = [
            part.strip()
            for part in text.split("|")
            if part.strip()
        ]

        symbol_like = [
            part
            for part in parts
            if not part.startswith(
                "ENSCAFG"
            )
        ]

        if parts:
            text = (
                symbol_like[-1]
                if symbol_like
                else parts[-1]
            )

    text = re.sub(
        r"\.\d+$",
        "",
        text,
    )

    return text


def normalize_expression_columns(
    expression: pd.DataFrame,
) -> pd.DataFrame:
    """Clean symbols and keep the first duplicate column."""
    result = expression.copy()

    result.columns = [
        clean_gene_symbol(column)
        for column in result.columns
    ]

    nonempty = np.asarray(
        [
            bool(column)
            for column in result.columns
        ],
        dtype=bool,
    )

    result = result.loc[
        :,
        nonempty,
    ].copy()

    result = result.loc[
        :,
        ~pd.Index(
            result.columns
        ).duplicated(
            keep="first"
        ),
    ].copy()

    return result


def standardize_expression(
    expression: pd.DataFrame,
) -> pd.DataFrame:
    """
    Standardize genes within one cohort.

    Missing values are median-imputed.
    Sample standard deviation uses ddof=1.
    Non-estimable genes are removed.
    """
    if expression.shape[0] < 2:
        raise ValueError(
            "Expression must contain at least "
            "two samples."
        )

    prepared = normalize_expression_columns(
        expression
    )

    x = prepared.apply(
        pd.to_numeric,
        errors="coerce",
    )

    x = x.fillna(
        x.median(axis=0)
    )

    standard_deviation = (
        x.std(
            axis=0,
            ddof=1,
        )
        .replace(
            0.0,
            np.nan,
        )
    )

    z = (
        x
        - x.mean(axis=0)
    ) / standard_deviation

    return z.loc[
        :,
        z.notna().all(axis=0),
    ].copy()


def upper_triangle(
    matrix: np.ndarray,
) -> np.ndarray:
    """Return the strict upper triangle of a square matrix."""
    array = np.asarray(
        matrix,
        dtype=float,
    )

    if (
        array.ndim != 2
        or array.shape[0]
        != array.shape[1]
    ):
        raise ValueError(
            "matrix must be square."
        )

    indices = np.triu_indices(
        array.shape[0],
        k=1,
    )

    return array[indices]


def frozen_signed_score(
    expression_z: pd.DataFrame,
    genes: Sequence[str],
    loadings: pd.Series,
) -> pd.Series:
    """Compute the unstandardized frozen signed-mean score."""
    selected = [
        str(gene)
        for gene in genes
    ]

    if not selected:
        return pd.Series(
            np.nan,
            index=expression_z.index,
            dtype=float,
        )

    aligned = pd.to_numeric(
        loadings.reindex(
            selected
        ),
        errors="coerce",
    ).fillna(0.0)

    signs = np.sign(
        aligned
    ).replace(
        0.0,
        1.0,
    )

    return (
        expression_z[
            selected
        ]
        .mul(
            signs,
            axis=1,
        )
        .mean(axis=1)
    )


def correlation_matrix(
    expression_z: pd.DataFrame,
    genes: Sequence[str],
) -> np.ndarray:
    """Compute the gene-gene Pearson correlation matrix."""
    selected = [
        str(gene)
        for gene in genes
    ]

    if len(selected) < 2:
        raise ValueError(
            "At least two genes are required "
            "for a correlation matrix."
        )

    return np.corrcoef(
        expression_z[
            selected
        ].to_numpy(
            dtype=float
        ),
        rowvar=False,
    )


def edge_spearman_from_matrices(
    reference_correlation: np.ndarray,
    target_correlation: np.ndarray,
) -> float:
    """Compare correlation-matrix edges using Spearman correlation."""
    reference = np.asarray(
        reference_correlation,
        dtype=float,
    )

    target = np.asarray(
        target_correlation,
        dtype=float,
    )

    if reference.shape != target.shape:
        raise ValueError(
            "Reference and target correlation "
            "matrices must have the same shape."
        )

    reference_edges = upper_triangle(
        reference
    )

    target_edges = upper_triangle(
        target
    )

    return float(
        stats.spearmanr(
            reference_edges,
            target_edges,
        ).statistic
    )


def edge_preservation(
    reference_z: pd.DataFrame,
    target_z: pd.DataFrame,
    genes: Sequence[str],
) -> float:
    """Compute direct within-program edge preservation."""
    reference_correlation = (
        correlation_matrix(
            reference_z,
            genes,
        )
    )

    target_correlation = (
        correlation_matrix(
            target_z,
            genes,
        )
    )

    return edge_spearman_from_matrices(
        reference_correlation,
        target_correlation,
    )


def fit_oriented_external_pc1(
    external_z: pd.DataFrame,
    genes: Sequence[str],
    frozen_loadings: pd.Series,
    *,
    random_state: int = DEFAULT_PCA_RANDOM_STATE,
) -> ExternalPC1:
    """
    Fit external PC1 and orient it using only the frozen signed score.
    """
    selected = [
        str(gene)
        for gene in genes
    ]

    if not selected:
        raise ValueError(
            "At least one gene is required "
            "for external PC1."
        )

    signed_score = frozen_signed_score(
        external_z,
        selected,
        frozen_loadings,
    )

    pca = PCA(
        n_components=1,
        random_state=random_state,
    )

    pc_score = pd.Series(
        pca.fit_transform(
            external_z[
                selected
            ].to_numpy(
                dtype=float
            )
        ).ravel(),
        index=external_z.index,
        dtype=float,
    )

    pc_loadings = (
        pca.components_[0]
        .copy()
    )

    orientation_correlation = (
        pc_score.corr(
            signed_score
        )
    )

    if (
        np.isfinite(
            orientation_correlation
        )
        and orientation_correlation
        < 0.0
    ):
        pc_score = -pc_score
        pc_loadings = (
            -pc_loadings
        )
        orientation_correlation = (
            -orientation_correlation
        )

    return ExternalPC1(
        scores=pc_score,
        loadings=np.asarray(
            pc_loadings,
            dtype=float,
        ),
        explained_variance_ratio=float(
            pca.explained_variance_ratio_[
                0
            ]
        ),
        orientation_correlation=(
            float(
                orientation_correlation
            )
            if np.isfinite(
                orientation_correlation
            )
            else np.nan
        ),
    )


def loading_preservation(
    frozen_loadings: pd.Series,
    external_pc1_loadings: np.ndarray,
    genes: Sequence[str],
) -> float:
    """Compare frozen loadings with oriented external PC1 loadings."""
    selected = [
        str(gene)
        for gene in genes
    ]

    aligned = pd.to_numeric(
        frozen_loadings.reindex(
            selected
        ),
        errors="coerce",
    ).fillna(0.0)

    external = np.asarray(
        external_pc1_loadings,
        dtype=float,
    )

    if len(aligned) != len(
        external
    ):
        raise ValueError(
            "Frozen and external loading "
            "vectors must have equal length."
        )

    return float(
        stats.spearmanr(
            aligned.to_numpy(
                dtype=float
            ),
            external,
        ).statistic
    )


def _prepare_module_weights(
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
            "Frozen weight table is missing "
            "required columns: "
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
            "No frozen genes found for "
            f"module {module_label}."
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

    if part.empty:
        raise ValueError(
            f"Module {module_label} has no "
            "usable frozen genes."
        )

    return part


def analyze_standardized_program(
    reference_z: pd.DataFrame,
    external_z: pd.DataFrame,
    weights: pd.DataFrame,
    module_label: str,
    *,
    minimum_genes: int = DEFAULT_MINIMUM_GENES,
    pca_random_state: int = DEFAULT_PCA_RANDOM_STATE,
) -> tuple[
    dict[str, Any],
    dict[str, Any],
]:
    """Analyze one frozen program using standardized cohorts."""
    if minimum_genes < 1:
        raise ValueError(
            "minimum_genes must be at least 1."
        )

    part = _prepare_module_weights(
        weights,
        module_label,
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

    coverage = {
        "module_label":
            str(module_label),
        "n_frozen_genes":
            int(
                len(
                    frozen_genes
                )
            ),
        "n_common_genes":
            int(
                len(
                    genes
                )
            ),
        "coverage_fraction":
            (
                float(
                    len(genes)
                    / len(
                        frozen_genes
                    )
                )
                if frozen_genes
                else np.nan
            ),
        "common_genes":
            ";".join(
                genes
            ),
        "missing_external_genes":
            ";".join(
                gene
                for gene
                in frozen_genes
                if gene
                not in external_z.columns
            ),
    }

    if len(genes) < minimum_genes:
        result = {
            "module_label":
                str(module_label),
            "n_common_genes":
                int(
                    len(
                        genes
                    )
                ),
            "edge_spearman":
                np.nan,
            "loading_spearman":
                np.nan,
            "external_pc1_variance_explained":
                np.nan,
            (
                "pc1_orientation_correlation_"
                "with_frozen_score"
            ):
                np.nan,
            "estimable":
                False,
            "nonestimable_reason":
                (
                    "fewer_than_"
                    f"{minimum_genes}_"
                    "common_genes"
                ),
        }

        return (
            result,
            coverage,
        )

    frozen_loadings = (
        part.loc[
            genes,
            "risk_oriented_loading",
        ]
    )

    edge_spearman = (
        edge_preservation(
            reference_z,
            external_z,
            genes,
        )
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

    result = {
        "module_label":
            str(module_label),
        "n_common_genes":
            int(
                len(
                    genes
                )
            ),
        "edge_spearman":
            float(
                edge_spearman
            ),
        "loading_spearman":
            float(
                loading_spearman
            ),
        "external_pc1_variance_explained":
            float(
                external_pc1
                .explained_variance_ratio
            ),
        (
            "pc1_orientation_correlation_"
            "with_frozen_score"
        ):
            float(
                external_pc1
                .orientation_correlation
            ),
        "estimable":
            True,
        "nonestimable_reason":
            "",
    }

    return (
        result,
        coverage,
    )


def analyze_program_preservation(
    reference_expression: pd.DataFrame,
    external_expression: pd.DataFrame,
    weights: pd.DataFrame,
    module_label: str,
    *,
    minimum_genes: int = DEFAULT_MINIMUM_GENES,
    pca_random_state: int = DEFAULT_PCA_RANDOM_STATE,
) -> tuple[
    dict[str, Any],
    dict[str, Any],
]:
    """Analyze one program from raw cohort expression matrices."""
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

    return analyze_standardized_program(
        reference_z=reference_z,
        external_z=external_z,
        weights=weights,
        module_label=module_label,
        minimum_genes=minimum_genes,
        pca_random_state=(
            pca_random_state
        ),
    )


def analyze_preservation(
    reference_expression: pd.DataFrame,
    external_expression: pd.DataFrame,
    weights: pd.DataFrame,
    *,
    module_order: Sequence[str] | None = None,
    minimum_genes: int = DEFAULT_MINIMUM_GENES,
    pca_random_state: int = DEFAULT_PCA_RANDOM_STATE,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
]:
    """Analyze direct preservation for multiple frozen programs."""
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

    if module_order is None:
        modules = (
            weights[
                "module_label"
            ]
            .dropna()
            .astype(str)
            .drop_duplicates()
            .tolist()
        )
    else:
        modules = [
            str(module)
            for module
            in module_order
        ]

    structure_rows: list[
        dict[str, Any]
    ] = []

    coverage_rows: list[
        dict[str, Any]
    ] = []

    for module in modules:
        (
            result,
            coverage,
        ) = analyze_standardized_program(
            reference_z=reference_z,
            external_z=external_z,
            weights=weights,
            module_label=module,
            minimum_genes=(
                minimum_genes
            ),
            pca_random_state=(
                pca_random_state
            ),
        )

        structure_rows.append(
            result
        )

        coverage_rows.append(
            coverage
        )

    return (
        pd.DataFrame(
            structure_rows
        ),
        pd.DataFrame(
            coverage_rows
        ),
    )
