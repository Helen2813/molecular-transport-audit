"""Unit tests for molecular-program preservation."""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from transport_audit.preservation import (
    analyze_program_preservation,
    clean_gene_symbol,
    edge_spearman_from_matrices,
    fit_oriented_external_pc1,
    standardize_expression,
    upper_triangle,
)


class GeneSymbolTests(
    unittest.TestCase
):
    def test_clean_gene_symbol(
        self,
    ) -> None:
        self.assertEqual(
            clean_gene_symbol(
                " tp53.2 "
            ),
            "TP53",
        )

        self.assertEqual(
            clean_gene_symbol(
                (
                    "ENSCAFG000001"
                    "|tp53.3"
                )
            ),
            "TP53",
        )

        self.assertEqual(
            clean_gene_symbol(
                "nan"
            ),
            "",
        )


class StandardizationTests(
    unittest.TestCase
):
    def test_standardization_uses_sample_sd(
        self,
    ) -> None:
        expression = pd.DataFrame(
            {
                "A": [
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                ],
                "B": [
                    4.0,
                    3.0,
                    2.0,
                    1.0,
                ],
            }
        )

        observed = (
            standardize_expression(
                expression
            )
        )

        self.assertAlmostEqual(
            float(
                observed[
                    "A"
                ].mean()
            ),
            0.0,
            places=12,
        )

        self.assertAlmostEqual(
            float(
                observed[
                    "A"
                ].std(
                    ddof=1
                )
            ),
            1.0,
            places=12,
        )

    def test_median_imputation(
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
            }
        )

        observed = (
            standardize_expression(
                expression
            )
        )

        raw = pd.Series(
            [
                1.0,
                2.0,
                2.0,
                4.0,
            ]
        )

        expected = (
            raw
            - raw.mean()
        ) / raw.std(
            ddof=1
        )

        np.testing.assert_allclose(
            observed[
                "A"
            ].to_numpy(),
            expected.to_numpy(),
            rtol=1.0e-12,
            atol=1.0e-12,
        )

    def test_zero_variance_gene_is_removed(
        self,
    ) -> None:
        expression = pd.DataFrame(
            {
                "A": [
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                ],
                "B": [
                    5.0,
                    5.0,
                    5.0,
                    5.0,
                ],
            }
        )

        observed = (
            standardize_expression(
                expression
            )
        )

        self.assertEqual(
            list(
                observed.columns
            ),
            ["A"],
        )

    def test_duplicate_symbol_keeps_first_column(
        self,
    ) -> None:
        expression = pd.DataFrame(
            np.array(
                [
                    [1.0, 10.0, 4.0],
                    [2.0, 20.0, 3.0],
                    [3.0, 30.0, 2.0],
                    [4.0, 40.0, 1.0],
                ]
            ),
            columns=[
                "a",
                "A",
                "B",
            ],
        )

        observed = (
            standardize_expression(
                expression
            )
        )

        expected_a = pd.Series(
            [
                1.0,
                2.0,
                3.0,
                4.0,
            ]
        )

        expected_a = (
            expected_a
            - expected_a.mean()
        ) / expected_a.std(
            ddof=1
        )

        np.testing.assert_allclose(
            observed[
                "A"
            ].to_numpy(),
            expected_a.to_numpy(),
            rtol=1.0e-12,
            atol=1.0e-12,
        )


class EdgePreservationTests(
    unittest.TestCase
):
    def test_upper_triangle_order(
        self,
    ) -> None:
        matrix = np.array(
            [
                [
                    1.0,
                    0.1,
                    0.8,
                ],
                [
                    0.1,
                    1.0,
                    -0.4,
                ],
                [
                    0.8,
                    -0.4,
                    1.0,
                ],
            ]
        )

        observed = (
            upper_triangle(
                matrix
            )
        )

        np.testing.assert_allclose(
            observed,
            np.array(
                [
                    0.1,
                    0.8,
                    -0.4,
                ]
            ),
        )

    def test_monotonic_edge_order_gives_one(
        self,
    ) -> None:
        reference = np.array(
            [
                [
                    1.0,
                    0.1,
                    0.8,
                ],
                [
                    0.1,
                    1.0,
                    -0.4,
                ],
                [
                    0.8,
                    -0.4,
                    1.0,
                ],
            ]
        )

        target = np.array(
            [
                [
                    1.0,
                    0.2,
                    0.9,
                ],
                [
                    0.2,
                    1.0,
                    -0.2,
                ],
                [
                    0.9,
                    -0.2,
                    1.0,
                ],
            ]
        )

        observed = (
            edge_spearman_from_matrices(
                reference,
                target,
            )
        )

        self.assertAlmostEqual(
            observed,
            1.0,
            places=12,
        )


class PC1OrientationTests(
    unittest.TestCase
):
    def test_pc1_is_oriented_to_frozen_score(
        self,
    ) -> None:
        expression = pd.DataFrame(
            {
                "A": [
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                ],
                "B": [
                    1.1,
                    2.2,
                    2.7,
                    4.3,
                    4.8,
                    6.2,
                ],
                "C": [
                    6.0,
                    5.0,
                    4.0,
                    3.0,
                    2.0,
                    1.0,
                ],
            }
        )

        z = standardize_expression(
            expression
        )

        loadings = pd.Series(
            {
                "A": 3.0,
                "B": 2.0,
                "C": -1.0,
            }
        )

        observed = (
            fit_oriented_external_pc1(
                external_z=z,
                genes=[
                    "A",
                    "B",
                    "C",
                ],
                frozen_loadings=loadings,
            )
        )

        self.assertGreater(
            observed.orientation_correlation,
            0.99,
        )

        self.assertGreater(
            observed.explained_variance_ratio,
            0.90,
        )

        self.assertEqual(
            list(
                observed.scores.index
            ),
            list(
                expression.index
            ),
        )


class ProgramPreservationTests(
    unittest.TestCase
):
    @staticmethod
    def weights() -> pd.DataFrame:
        return pd.DataFrame(
            {
                "module_label": [
                    "MTEST",
                    "MTEST",
                    "MTEST",
                ],
                "canine_gene_symbol": [
                    "A",
                    "B",
                    "C",
                ],
                "risk_oriented_loading": [
                    2.0,
                    1.0,
                    -1.0,
                ],
            }
        )

    def test_exact_three_common_genes_is_estimable(
        self,
    ) -> None:
        reference = pd.DataFrame(
            {
                "A": [
                    1.0,
                    2.0,
                    3.0,
                    5.0,
                    4.0,
                ],
                "B": [
                    2.0,
                    1.0,
                    4.0,
                    3.0,
                    5.0,
                ],
                "C": [
                    5.0,
                    4.0,
                    2.0,
                    1.0,
                    3.0,
                ],
            }
        )

        external = pd.DataFrame(
            {
                "A": [
                    1.0,
                    3.0,
                    2.0,
                    5.0,
                    4.0,
                ],
                "B": [
                    2.0,
                    1.0,
                    5.0,
                    3.0,
                    4.0,
                ],
                "C": [
                    5.0,
                    3.0,
                    4.0,
                    1.0,
                    2.0,
                ],
            }
        )

        (
            result,
            coverage,
        ) = analyze_program_preservation(
            reference_expression=reference,
            external_expression=external,
            weights=self.weights(),
            module_label="MTEST",
        )

        self.assertTrue(
            result[
                "estimable"
            ]
        )

        self.assertEqual(
            result[
                "n_common_genes"
            ],
            3,
        )

        self.assertEqual(
            coverage[
                "n_frozen_genes"
            ],
            3,
        )

        self.assertEqual(
            coverage[
                "n_common_genes"
            ],
            3,
        )

        self.assertAlmostEqual(
            coverage[
                "coverage_fraction"
            ],
            1.0,
            places=12,
        )

    def test_zero_variance_can_make_program_nonestimable(
        self,
    ) -> None:
        reference = pd.DataFrame(
            {
                "A": [
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                ],
                "B": [
                    2.0,
                    4.0,
                    1.0,
                    3.0,
                ],
                "C": [
                    4.0,
                    3.0,
                    2.0,
                    1.0,
                ],
            }
        )

        external = pd.DataFrame(
            {
                "A": [
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                ],
                "B": [
                    4.0,
                    3.0,
                    2.0,
                    1.0,
                ],
                "C": [
                    7.0,
                    7.0,
                    7.0,
                    7.0,
                ],
            }
        )

        (
            result,
            coverage,
        ) = analyze_program_preservation(
            reference_expression=reference,
            external_expression=external,
            weights=self.weights(),
            module_label="MTEST",
        )

        self.assertFalse(
            result[
                "estimable"
            ]
        )

        self.assertEqual(
            coverage[
                "n_common_genes"
            ],
            2,
        )

    def test_duplicate_frozen_symbol_keeps_first(
        self,
    ) -> None:
        weights = pd.DataFrame(
            {
                "module_label": [
                    "MTEST",
                    "MTEST",
                    "MTEST",
                    "MTEST",
                ],
                "canine_gene_symbol": [
                    "A",
                    "A",
                    "B",
                    "C",
                ],
                "risk_oriented_loading": [
                    2.0,
                    -100.0,
                    1.0,
                    -1.0,
                ],
            }
        )

        reference = pd.DataFrame(
            {
                "A": [
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                ],
                "B": [
                    2.0,
                    4.0,
                    1.0,
                    3.0,
                ],
                "C": [
                    4.0,
                    3.0,
                    2.0,
                    1.0,
                ],
            }
        )

        external = reference.copy()

        (
            result,
            coverage,
        ) = analyze_program_preservation(
            reference_expression=reference,
            external_expression=external,
            weights=weights,
            module_label="MTEST",
        )

        self.assertTrue(
            result[
                "estimable"
            ]
        )

        self.assertEqual(
            coverage[
                "n_frozen_genes"
            ],
            3,
        )


if __name__ == "__main__":
    unittest.main()
