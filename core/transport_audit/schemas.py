"""Data models for frozen-program scoring."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class ScoringRule:
    """Frozen scoring rule."""

    mapping_name: str = "strict"
    minimum_genes: int = 3
    minimum_fraction: float = 0.50

    def __post_init__(self) -> None:
        if self.minimum_genes < 1:
            raise ValueError("minimum_genes must be at least 1.")

        if not 0.0 <= self.minimum_fraction <= 1.0:
            raise ValueError(
                "minimum_fraction must lie within [0, 1]."
            )


@dataclass(frozen=True)
class FrozenGeneWeight:
    """One frozen gene and its risk-oriented loading."""

    gene_symbol: str
    risk_oriented_loading: float


@dataclass(frozen=True)
class FrozenProgram:
    """Frozen molecular program."""

    module_label: str
    genes: tuple[FrozenGeneWeight, ...]

    @classmethod
    def from_frame(
        cls,
        weights: pd.DataFrame,
        module_label: str,
    ) -> "FrozenProgram":
        required = {
            "module_label",
            "human_gene_symbol",
            "risk_oriented_loading",
        }
        missing = sorted(
            required.difference(weights.columns)
        )

        if missing:
            raise ValueError(
                "Weight table is missing required columns: "
                + ", ".join(missing)
            )

        part = weights[
            weights["module_label"]
            .astype(str)
            .eq(str(module_label))
        ].copy()

        if part.empty:
            raise ValueError(
                f"No frozen weights found for module "
                f"{module_label}."
            )

        part["human_gene_symbol"] = (
            part["human_gene_symbol"]
            .astype(str)
            .str.strip()
            .str.upper()
        )

        part = part[
            ~part["human_gene_symbol"].isin(
                ["", "NAN", "NONE"]
            )
        ].copy()

        part = part.drop_duplicates(
            "human_gene_symbol",
            keep="first",
        )

        part["risk_oriented_loading"] = (
            pd.to_numeric(
                part["risk_oriented_loading"],
                errors="coerce",
            )
            .fillna(0.0)
        )

        genes = tuple(
            FrozenGeneWeight(
                gene_symbol=str(
                    row.human_gene_symbol
                ),
                risk_oriented_loading=float(
                    row.risk_oriented_loading
                ),
            )
            for row in part.itertuples(
                index=False
            )
        )

        if not genes:
            raise ValueError(
                f"Module {module_label} has no usable "
                "frozen genes."
            )

        return cls(
            module_label=str(module_label),
            genes=genes,
        )

    @property
    def gene_symbols(self) -> tuple[str, ...]:
        return tuple(
            gene.gene_symbol
            for gene in self.genes
        )

    def loading_series(self) -> pd.Series:
        return pd.Series(
            {
                gene.gene_symbol:
                gene.risk_oriented_loading
                for gene in self.genes
            },
            dtype=float,
        )


@dataclass(frozen=True)
class ScoreCoverage:
    """Coverage of one frozen program."""

    cohort: str
    module_label: str
    mapping: str
    n_frozen_genes: int
    n_available_genes: int
    n_scored_genes: int
    coverage_fraction: float
    minimum_rule_passed: bool
    available_genes: tuple[str, ...]
    missing_genes: tuple[str, ...]

    def to_legacy_record(
        self,
    ) -> dict[str, Any]:
        return {
            "cohort": self.cohort,
            "module_label": self.module_label,
            "mapping": self.mapping,
            "n_frozen_genes":
                self.n_frozen_genes,
            "n_available_genes":
                self.n_available_genes,
            "coverage_fraction":
                self.coverage_fraction,
            "minimum_rule_passed":
                self.minimum_rule_passed,
            "available_genes":
                ";".join(self.available_genes),
            "missing_genes":
                ";".join(self.missing_genes),
        }


@dataclass(frozen=True)
class ProgramScoreResult:
    """Scores and coverage for one frozen program."""

    module_label: str
    coverage: ScoreCoverage
    signed_mean_z: pd.Series | None
    canine_pca_weighted_z: pd.Series | None

    @property
    def estimable(self) -> bool:
        return (
            self.coverage.minimum_rule_passed
            and self.signed_mean_z is not None
            and self.canine_pca_weighted_z
            is not None
        )
