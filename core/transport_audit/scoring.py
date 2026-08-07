"""Frozen molecular-program scoring."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from transport_audit.expression import (
    normalize_expression_columns,
    standardize_genes,
    zscore_series,
)
from transport_audit.schemas import (
    FrozenProgram,
    ProgramScoreResult,
    ScoreCoverage,
    ScoringRule,
)


def score_signed_mean(
    standardized_expression: pd.DataFrame,
    loadings: pd.Series,
) -> pd.Series:
    """Compute the frozen signed-mean score."""
    aligned = pd.to_numeric(
        loadings.reindex(
            standardized_expression.columns
        ),
        errors="coerce",
    ).fillna(0.0)

    signs = (
        np.sign(aligned)
        .replace(
            0.0,
            1.0,
        )
    )

    raw_score = (
        standardized_expression
        .mul(
            signs,
            axis=1,
        )
        .mean(axis=1)
    )

    return zscore_series(
        raw_score
    )


def score_canine_pca_weighted(
    standardized_expression: pd.DataFrame,
    loadings: pd.Series,
) -> pd.Series:
    """Compute the frozen PCA-weighted score."""
    aligned = pd.to_numeric(
        loadings.reindex(
            standardized_expression.columns
        ),
        errors="coerce",
    ).fillna(0.0)

    denominator = float(
        aligned.abs().sum()
    )

    if (
        not np.isfinite(
            denominator
        )
        or denominator <= 0.0
    ):
        return pd.Series(
            np.nan,
            index=standardized_expression.index,
            dtype=float,
        )

    normalized = (
        aligned
        / denominator
    )

    raw_score = (
        standardized_expression
        .mul(
            normalized,
            axis=1,
        )
        .sum(axis=1)
    )

    return zscore_series(
        raw_score
    )


def score_frozen_program(
    expression: pd.DataFrame,
    program: FrozenProgram,
    cohort: str,
    rule: ScoringRule | None = None,
) -> ProgramScoreResult:
    """Score one frozen program without outcome tuning."""
    active_rule = (
        rule
        or ScoringRule()
    )

    prepared = (
        normalize_expression_columns(
            expression
        )
    )

    requested = list(
        program.gene_symbols
    )

    available = [
        gene
        for gene in requested
        if gene in prepared.columns
    ]

    missing = [
        gene
        for gene in requested
        if gene not in prepared.columns
    ]

    n_frozen = len(
        requested
    )

    n_available = len(
        available
    )

    coverage_fraction = (
        n_available / n_frozen
        if n_frozen
        else 0.0
    )

    minimum_rule_passed = (
        n_available
        >= active_rule.minimum_genes
        and coverage_fraction
        >= active_rule.minimum_fraction
    )

    if minimum_rule_passed:
        standardized = (
            standardize_genes(
                prepared,
                available,
            )
        )
    else:
        standardized = pd.DataFrame(
            index=prepared.index
        )

    n_scored = int(
        standardized.shape[1]
    )

    coverage = ScoreCoverage(
        cohort=str(cohort),
        module_label=program.module_label,
        mapping=active_rule.mapping_name,
        n_frozen_genes=n_frozen,
        n_available_genes=n_available,
        n_scored_genes=n_scored,
        coverage_fraction=float(
            coverage_fraction
        ),
        minimum_rule_passed=bool(
            minimum_rule_passed
        ),
        available_genes=tuple(
            available
        ),
        missing_genes=tuple(
            missing
        ),
    )

    if (
        not minimum_rule_passed
        or n_scored
        < active_rule.minimum_genes
    ):
        return ProgramScoreResult(
            module_label=program.module_label,
            coverage=coverage,
            signed_mean_z=None,
            canine_pca_weighted_z=None,
        )

    loadings = (
        program.loading_series()
    )

    signed_mean = (
        score_signed_mean(
            standardized,
            loadings,
        )
    )

    weighted = (
        score_canine_pca_weighted(
            standardized,
            loadings,
        )
    )

    return ProgramScoreResult(
        module_label=program.module_label,
        coverage=coverage,
        signed_mean_z=signed_mean,
        canine_pca_weighted_z=weighted,
    )


def score_programs(
    expression: pd.DataFrame,
    weights: pd.DataFrame,
    cohort: str,
    module_order: Sequence[str] | None = None,
    rule: ScoringRule | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score multiple frozen programs."""
    active_rule = (
        rule
        or ScoringRule()
    )

    prepared = (
        normalize_expression_columns(
            expression
        )
    )

    if module_order is None:
        modules = (
            weights["module_label"]
            .dropna()
            .astype(str)
            .drop_duplicates()
            .tolist()
        )
    else:
        modules = [
            str(module)
            for module in module_order
        ]

    scores = pd.DataFrame(
        index=prepared.index
    )

    coverage_rows: list[
        dict[str, object]
    ] = []

    for module in modules:
        program = (
            FrozenProgram.from_frame(
                weights,
                module,
            )
        )

        result = (
            score_frozen_program(
                expression=prepared,
                program=program,
                cohort=cohort,
                rule=active_rule,
            )
        )

        coverage_rows.append(
            result.coverage
            .to_legacy_record()
        )

        if not result.estimable:
            continue

        prefix = (
            f"{module}"
            f"__{active_rule.mapping_name}"
        )

        scores[
            f"{prefix}__signed_mean_z"
        ] = result.signed_mean_z

        scores[
            f"{prefix}"
            "__canine_pca_weighted_z"
        ] = (
            result.canine_pca_weighted_z
        )

    coverage = pd.DataFrame(
        coverage_rows
    )

    return scores, coverage
