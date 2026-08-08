"""Unit tests for frozen-program reliability."""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from transport_audit.preservation import (
    frozen_signed_score,
)
from transport_audit.preservation_reliability import (
    gene_leave_one_out,
    reliability_seed,
    split_half_reliability,
)


class ReliabilitySeedTests(
    unittest.TestCase
):
    def test_frozen_seed_schedule(
        self,
    ) -> None:
        self.assertEqual(
            reliability_seed(0),
            44,
        )
        self.assertEqual(
            reliability_seed(1),
            144,
        )
        self.assertEqual(
            reliability_seed(2),
            244,
        )
        self.assertEqual(
            reliability_seed(3),
            344,
        )


class SplitHalfTests(
    unittest.TestCase
):
    def test_fewer_than_four_genes_is_not_estimable(
        self,
    ) -> None:
        expression = pd.DataFrame(
            {
                "A": [1.0, 2.0, 3.0],
                "B": [2.0, 1.0, 3.0],
                "C": [3.0, 2.0, 1.0],
            }
        )

        loadings = pd.Series(
            {
                "A": 1.0,
                "B": -1.0,
                "C": 1.0,
            }
        )

        result = split_half_reliability(
            expression,
            ["A", "B", "C"],
            loadings,
            seed=44,
            n_repeats=10,
        )

        self.assertEqual(
            result[
                "split_half_valid_repeats"
            ],
            0,
        )

        self.assertTrue(
            np.isnan(
                result[
                    "split_half_median"
                ]
            )
        )

    def test_split_half_matches_manual_legacy_loop(
        self,
    ) -> None:
        expression = pd.DataFrame(
            {
                "A": [-1.0, -0.5, 0.2, 0.7, 1.1],
                "B": [-0.8, -0.4, 0.1, 0.9, 1.0],
                "C": [1.0, 0.6, 0.2, -0.4, -0.9],
                "D": [0.9, 0.5, 0.3, -0.5, -0.8],
                "E": [-0.3, 0.8, -0.5, 0.9, -0.1],
                "F": [0.2, -0.7, 0.8, -0.4, 0.1],
            }
        )

        genes = [
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
        ]

        loadings = pd.Series(
            {
                "A": 1.0,
                "B": 2.0,
                "C": -1.0,
                "D": -2.0,
                "E": 1.0,
                "F": -1.0,
            }
        )

        seed = 77
        repeats = 37

        rng = np.random.default_rng(
            seed
        )

        manual: list[float] = []

        for _ in range(
            repeats
        ):
            shuffled = np.asarray(
                genes,
                dtype=object,
            )

            shuffled = rng.permutation(
                shuffled
            )

            midpoint = (
                len(shuffled)
                // 2
            )

            first = shuffled[
                :midpoint
            ].tolist()

            second = shuffled[
                midpoint:
            ].tolist()

            first_score = (
                frozen_signed_score(
                    expression,
                    first,
                    loadings,
                )
            )

            second_score = (
                frozen_signed_score(
                    expression,
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
                manual.append(
                    float(
                        correlation
                    )
                )

        observed = (
            split_half_reliability(
                expression,
                genes,
                loadings,
                seed=seed,
                n_repeats=repeats,
            )
        )

        self.assertEqual(
            observed[
                "split_half_valid_repeats"
            ],
            len(manual),
        )

        self.assertAlmostEqual(
            observed[
                "split_half_median"
            ],
            float(
                np.median(
                    manual
                )
            ),
            places=15,
        )

        self.assertAlmostEqual(
            observed[
                "split_half_q05"
            ],
            float(
                np.quantile(
                    manual,
                    0.05,
                )
            ),
            places=15,
        )

        self.assertAlmostEqual(
            observed[
                "split_half_q95"
            ],
            float(
                np.quantile(
                    manual,
                    0.95,
                )
            ),
            places=15,
        )


class GeneLeaveOneOutTests(
    unittest.TestCase
):
    def test_loo_matches_manual_full_score_correlations(
        self,
    ) -> None:
        expression = pd.DataFrame(
            {
                "A": [-1.0, -0.5, 0.2, 0.7, 1.1],
                "B": [-0.8, -0.4, 0.1, 0.9, 1.0],
                "C": [1.0, 0.6, 0.2, -0.4, -0.9],
                "D": [0.3, -0.7, 0.8, -0.2, 0.1],
            }
        )

        genes = [
            "A",
            "B",
            "C",
            "D",
        ]

        loadings = pd.Series(
            {
                "A": 2.0,
                "B": 1.0,
                "C": -1.0,
                "D": 1.0,
            }
        )

        full_score = frozen_signed_score(
            expression,
            genes,
            loadings,
        )

        observed = gene_leave_one_out(
            expression_z=expression,
            module_label="MTEST",
            genes=genes,
            loadings=loadings,
        )

        self.assertEqual(
            observed.shape[0],
            4,
        )

        for row in observed.itertuples(
            index=False
        ):
            subset = [
                gene
                for gene in genes
                if gene
                != row.left_out_gene
            ]

            manual_score = (
                frozen_signed_score(
                    expression,
                    subset,
                    loadings,
                )
            )

            expected = float(
                full_score.corr(
                    manual_score
                )
            )

            self.assertAlmostEqual(
                row.correlation_with_full_score,
                expected,
                places=15,
            )

            self.assertEqual(
                row.n_genes_remaining,
                3,
            )

    def test_three_gene_module_has_no_loo(
        self,
    ) -> None:
        expression = pd.DataFrame(
            {
                "A": [1.0, 2.0, 3.0, 4.0],
                "B": [4.0, 3.0, 2.0, 1.0],
                "C": [1.0, 3.0, 2.0, 4.0],
            }
        )

        loadings = pd.Series(
            {
                "A": 1.0,
                "B": -1.0,
                "C": 1.0,
            }
        )

        observed = gene_leave_one_out(
            expression_z=expression,
            module_label="MTEST",
            genes=[
                "A",
                "B",
                "C",
            ],
            loadings=loadings,
            minimum_genes=3,
        )

        self.assertTrue(
            observed.empty
        )


if __name__ == "__main__":
    unittest.main()
