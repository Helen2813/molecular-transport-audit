"""Unit tests for frozen molecular-program scoring."""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from transport_audit.expression import (
    normalize_expression_columns,
    standardize_genes,
    zscore_series,
)
from transport_audit.schemas import (
    FrozenProgram,
    ScoringRule,
)
from transport_audit.scoring import (
    score_canine_pca_weighted,
    score_frozen_program,
    score_signed_mean,
)


class ExpressionPreprocessingTests(unittest.TestCase):
    def test_gene_symbols_are_uppercased_and_first_duplicate_is_kept(
        self,
    ) -> None:
        expression = pd.DataFrame(
            {
                "genea": [1.0, 2.0, 3.0],
                "GENEA": [100.0, 200.0, 300.0],
                " geneb ": [4.0, 5.0, 6.0],
            },
            index=["S1", "S2", "S3"],
        )

        observed = normalize_expression_columns(
            expression
        )

        self.assertEqual(
            list(observed.columns),
            ["GENEA", "GENEB"],
        )

        np.testing.assert_allclose(
            observed["GENEA"].to_numpy(),
            np.array(
                [1.0, 2.0, 3.0]
            ),
        )

    def test_median_imputation_and_sample_standardization(
        self,
    ) -> None:
        expression = pd.DataFrame(
            {
                "A": [
                    1.0,
                    2.0,
                    np.nan,
                    4.0,
                ],
                "B": [
                    10.0,
                    20.0,
                    30.0,
                    40.0,
                ],
            },
            index=[
                "S1",
                "S2",
                "S3",
                "S4",
            ],
        )

        observed = standardize_genes(
            expression,
            ["A", "B"],
        )

        expected_a_raw = pd.Series(
            [1.0, 2.0, 2.0, 4.0],
            index=expression.index,
        )

        expected_a = (
            expected_a_raw
            - expected_a_raw.mean()
        ) / expected_a_raw.std(
            ddof=1
        )

        np.testing.assert_allclose(
            observed["A"].to_numpy(),
            expected_a.to_numpy(),
            rtol=1.0e-12,
            atol=1.0e-12,
        )

        self.assertAlmostEqual(
            float(
                observed["A"].mean()
            ),
            0.0,
            places=12,
        )

        self.assertAlmostEqual(
            float(
                observed["A"].std(
                    ddof=1
                )
            ),
            1.0,
            places=12,
        )

    def test_zero_variance_gene_is_removed(
        self,
    ) -> None:
        expression = pd.DataFrame(
            {
                "A": [1.0, 2.0, 3.0, 4.0],
                "B": [5.0, 5.0, 5.0, 5.0],
                "C": [4.0, 3.0, 2.0, 1.0],
            }
        )

        observed = standardize_genes(
            expression,
            ["A", "B", "C"],
        )

        self.assertEqual(
            list(observed.columns),
            ["A", "C"],
        )

    def test_constant_score_returns_nan_zscore(
        self,
    ) -> None:
        score = pd.Series(
            [3.0, 3.0, 3.0, 3.0]
        )

        observed = zscore_series(
            score
        )

        self.assertTrue(
            observed.isna().all()
        )


class FrozenProgramTests(unittest.TestCase):
    def test_duplicate_frozen_gene_keeps_first_loading(
        self,
    ) -> None:
        weights = pd.DataFrame(
            {
                "module_label": [
                    "MTEST",
                    "MTEST",
                    "MTEST",
                ],
                "human_gene_symbol": [
                    "A",
                    "A",
                    "B",
                ],
                "risk_oriented_loading": [
                    1.0,
                    -10.0,
                    -2.0,
                ],
            }
        )

        program = FrozenProgram.from_frame(
            weights,
            "MTEST",
        )

        self.assertEqual(
            program.gene_symbols,
            ("A", "B"),
        )

        loadings = (
            program.loading_series()
        )

        self.assertEqual(
            float(loadings["A"]),
            1.0,
        )

        self.assertEqual(
            float(loadings["B"]),
            -2.0,
        )


class ScoreFunctionTests(unittest.TestCase):
    def test_zero_loading_is_positive_for_signed_mean(
        self,
    ) -> None:
        z = pd.DataFrame(
            {
                "A": [-1.0, 0.0, 1.0],
                "B": [-1.0, 0.0, 1.0],
                "C": [1.0, 0.0, -1.0],
            },
            index=["S1", "S2", "S3"],
        )

        loadings = pd.Series(
            {
                "A": 0.0,
                "B": -2.0,
                "C": 3.0,
            }
        )

        observed = score_signed_mean(
            z,
            loadings,
        )

        signs = pd.Series(
            {
                "A": 1.0,
                "B": -1.0,
                "C": 1.0,
            }
        )

        expected_raw = (
            z.mul(
                signs,
                axis=1,
            )
            .mean(axis=1)
        )

        expected = zscore_series(
            expected_raw
        )

        np.testing.assert_allclose(
            observed.to_numpy(),
            expected.to_numpy(),
            rtol=1.0e-12,
            atol=1.0e-12,
        )

    def test_all_zero_loadings_return_nan_weighted_score(
        self,
    ) -> None:
        z = pd.DataFrame(
            {
                "A": [-1.0, 0.0, 1.0],
                "B": [1.0, 0.0, -1.0],
                "C": [-1.0, 1.0, 0.0],
            }
        )

        loadings = pd.Series(
            {
                "A": 0.0,
                "B": 0.0,
                "C": 0.0,
            }
        )

        observed = (
            score_canine_pca_weighted(
                z,
                loadings,
            )
        )

        self.assertTrue(
            observed.isna().all()
        )

    def test_weighted_score_uses_only_available_loadings(
        self,
    ) -> None:
        expression = pd.DataFrame(
            {
                "A": [1.0, 2.0, 3.0, 4.0],
                "B": [4.0, 3.0, 2.0, 1.0],
                "D": [2.0, 4.0, 1.0, 3.0],
            },
            index=[
                "S1",
                "S2",
                "S3",
                "S4",
            ],
        )

        weights = pd.DataFrame(
            {
                "module_label": [
                    "MTEST",
                    "MTEST",
                    "MTEST",
                    "MTEST",
                ],
                "human_gene_symbol": [
                    "A",
                    "B",
                    "C",
                    "D",
                ],
                "risk_oriented_loading": [
                    1.0,
                    2.0,
                    100.0,
                    -1.0,
                ],
            }
        )

        program = FrozenProgram.from_frame(
            weights,
            "MTEST",
        )

        rule = ScoringRule(
            minimum_genes=3,
            minimum_fraction=0.50,
        )

        result = score_frozen_program(
            expression=expression,
            program=program,
            cohort="TEST",
            rule=rule,
        )

        self.assertTrue(
            result.estimable
        )

        self.assertEqual(
            result.coverage.n_frozen_genes,
            4,
        )

        self.assertEqual(
            result.coverage.n_available_genes,
            3,
        )

        self.assertEqual(
            result.coverage.missing_genes,
            ("C",),
        )

        z = standardize_genes(
            expression,
            ["A", "B", "D"],
        )

        expected_weights = pd.Series(
            {
                "A": 1.0 / 4.0,
                "B": 2.0 / 4.0,
                "D": -1.0 / 4.0,
            }
        )

        expected_raw = (
            z.mul(
                expected_weights,
                axis=1,
            )
            .sum(axis=1)
        )

        expected = zscore_series(
            expected_raw
        )

        np.testing.assert_allclose(
            result.canine_pca_weighted_z
            .to_numpy(),
            expected.to_numpy(),
            rtol=1.0e-12,
            atol=1.0e-12,
        )


class CoverageRuleTests(unittest.TestCase):
    @staticmethod
    def build_program(
        n_genes: int,
    ) -> FrozenProgram:
        weights = pd.DataFrame(
            {
                "module_label": [
                    "MTEST"
                ]
                * n_genes,
                "human_gene_symbol": [
                    f"G{i}"
                    for i
                    in range(
                        1,
                        n_genes + 1,
                    )
                ],
                "risk_oriented_loading": [
                    1.0
                    if i % 2
                    else -1.0
                    for i
                    in range(
                        1,
                        n_genes + 1,
                    )
                ],
            }
        )

        return FrozenProgram.from_frame(
            weights,
            "MTEST",
        )

    def test_exact_fifty_percent_coverage_passes(
        self,
    ) -> None:
        program = self.build_program(
            4
        )

        expression = pd.DataFrame(
            {
                "G1": [
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                ],
                "G2": [
                    4.0,
                    3.0,
                    2.0,
                    1.0,
                ],
            }
        )

        rule = ScoringRule(
            minimum_genes=2,
            minimum_fraction=0.50,
        )

        result = score_frozen_program(
            expression=expression,
            program=program,
            cohort="TEST",
            rule=rule,
        )

        self.assertTrue(
            result.coverage
            .minimum_rule_passed
        )

        self.assertTrue(
            result.estimable
        )

        self.assertAlmostEqual(
            result.coverage
            .coverage_fraction,
            0.50,
            places=12,
        )

    def test_gene_count_can_pass_while_fraction_fails(
        self,
    ) -> None:
        program = self.build_program(
            7
        )

        expression = pd.DataFrame(
            {
                "G1": [
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                ],
                "G2": [
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                ],
                "G3": [
                    5.0,
                    4.0,
                    3.0,
                    2.0,
                ],
            }
        )

        result = score_frozen_program(
            expression=expression,
            program=program,
            cohort="TEST",
            rule=ScoringRule(),
        )

        self.assertEqual(
            result.coverage
            .n_available_genes,
            3,
        )

        self.assertLess(
            result.coverage
            .coverage_fraction,
            0.50,
        )

        self.assertFalse(
            result.coverage
            .minimum_rule_passed
        )

        self.assertFalse(
            result.estimable
        )

    def test_zero_variance_filter_can_make_score_nonestimable(
        self,
    ) -> None:
        program = self.build_program(
            3
        )

        expression = pd.DataFrame(
            {
                "G1": [
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                ],
                "G2": [
                    5.0,
                    5.0,
                    5.0,
                    5.0,
                ],
                "G3": [
                    4.0,
                    3.0,
                    2.0,
                    1.0,
                ],
            }
        )

        result = score_frozen_program(
            expression=expression,
            program=program,
            cohort="TEST",
            rule=ScoringRule(),
        )

        self.assertTrue(
            result.coverage
            .minimum_rule_passed
        )

        self.assertEqual(
            result.coverage
            .n_available_genes,
            3,
        )

        self.assertEqual(
            result.coverage
            .n_scored_genes,
            2,
        )

        self.assertFalse(
            result.estimable
        )


if __name__ == "__main__":
    unittest.main()
