"""Unit tests for preservation permutation inference."""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd
from scipy import stats

from transport_audit.multiplicity import (
    benjamini_hochberg,
)
from transport_audit.preservation import (
    upper_triangle,
)
from transport_audit.preservation_inference import (
    adjust_direct_permutation_family,
    edge_gene_label_permutation,
    empirical_two_sided_p,
    loading_gene_label_permutation,
    permutation_seeds,
)


class EmpiricalPValueTests(
    unittest.TestCase
):
    def test_plus_one_correction(
        self,
    ) -> None:
        observed = (
            empirical_two_sided_p(
                extreme_count=0,
                n_permutations=5000,
            )
        )

        self.assertAlmostEqual(
            observed,
            1.0 / 5001.0,
            places=15,
        )

    def test_invalid_extreme_count_raises(
        self,
    ) -> None:
        with self.assertRaises(
            ValueError
        ):
            empirical_two_sided_p(
                extreme_count=11,
                n_permutations=10,
            )


class SeedScheduleTests(
    unittest.TestCase
):
    def test_frozen_module_seed_schedule(
        self,
    ) -> None:
        self.assertEqual(
            permutation_seeds(0),
            (42, 43),
        )

        self.assertEqual(
            permutation_seeds(1),
            (142, 143),
        )

        self.assertEqual(
            permutation_seeds(2),
            (242, 243),
        )

        self.assertEqual(
            permutation_seeds(3),
            (342, 343),
        )


class EdgePermutationTests(
    unittest.TestCase
):
    def test_edge_permutation_matches_manual_legacy_loop(
        self,
    ) -> None:
        reference = np.array(
            [
                [1.0, 0.2, 0.7, -0.4],
                [0.2, 1.0, -0.1, 0.5],
                [0.7, -0.1, 1.0, 0.3],
                [-0.4, 0.5, 0.3, 1.0],
            ],
            dtype=float,
        )

        external = np.array(
            [
                [1.0, 0.4, 0.8, -0.2],
                [0.4, 1.0, 0.1, 0.6],
                [0.8, 0.1, 1.0, 0.2],
                [-0.2, 0.6, 0.2, 1.0],
            ],
            dtype=float,
        )

        reference_edges = (
            upper_triangle(
                reference
            )
        )

        external_edges = (
            upper_triangle(
                external
            )
        )

        observed = float(
            stats.spearmanr(
                reference_edges,
                external_edges,
            ).statistic
        )

        n_permutations = 31
        seed = 7

        rng = np.random.default_rng(
            seed
        )

        manual_count = 0

        for _ in range(
            n_permutations
        ):
            permutation = (
                rng.permutation(4)
            )

            permuted = external[
                np.ix_(
                    permutation,
                    permutation,
                )
            ]

            value = stats.spearmanr(
                reference_edges,
                upper_triangle(
                    permuted
                ),
            ).statistic

            if (
                np.isfinite(value)
                and abs(value)
                >= abs(observed)
            ):
                manual_count += 1

        result = (
            edge_gene_label_permutation(
                reference_edges,
                external,
                observed,
                n_permutations=(
                    n_permutations
                ),
                seed=seed,
            )
        )

        self.assertEqual(
            result.extreme_count,
            manual_count,
        )

        self.assertEqual(
            result.p_value,
            empirical_two_sided_p(
                manual_count,
                n_permutations,
            ),
        )


class LoadingPermutationTests(
    unittest.TestCase
):
    def test_loading_permutation_matches_manual_legacy_loop(
        self,
    ) -> None:
        frozen = np.array(
            [
                0.8,
                -0.4,
                0.2,
                0.6,
                -0.1,
            ]
        )

        external = np.array(
            [
                0.7,
                -0.2,
                0.4,
                0.5,
                -0.3,
            ]
        )

        observed = float(
            stats.spearmanr(
                frozen,
                external,
            ).statistic
        )

        n_permutations = 37
        seed = 11

        rng = np.random.default_rng(
            seed
        )

        manual_count = 0

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
                manual_count += 1

        result = (
            loading_gene_label_permutation(
                frozen,
                external,
                observed,
                n_permutations=(
                    n_permutations
                ),
                seed=seed,
            )
        )

        self.assertEqual(
            result.extreme_count,
            manual_count,
        )

        self.assertEqual(
            result.p_value,
            empirical_two_sided_p(
                manual_count,
                n_permutations,
            ),
        )


class DirectFamilyTests(
    unittest.TestCase
):
    def test_bh_is_joint_across_eight_tests(
        self,
    ) -> None:
        table = pd.DataFrame(
            {
                "module_label": [
                    "M34",
                    "M11",
                    "M24",
                    "M40",
                ],
                "edge_permutation_p": [
                    0.001,
                    0.20,
                    0.30,
                    0.002,
                ],
                "loading_permutation_p": [
                    0.003,
                    0.70,
                    0.25,
                    0.004,
                ],
            }
        )

        observed = (
            adjust_direct_permutation_family(
                table
            )
        )

        family = np.array(
            [
                0.001,
                0.20,
                0.30,
                0.002,
                0.003,
                0.70,
                0.25,
                0.004,
            ]
        )

        expected = (
            benjamini_hochberg(
                family
            )
        )

        np.testing.assert_allclose(
            observed[
                "edge_q_bh_8"
            ].to_numpy(
                dtype=float
            ),
            expected[:4],
            rtol=0.0,
            atol=1.0e-15,
        )

        np.testing.assert_allclose(
            observed[
                "loading_q_bh_8"
            ].to_numpy(
                dtype=float
            ),
            expected[4:],
            rtol=0.0,
            atol=1.0e-15,
        )

        self.assertTrue(
            (
                observed[
                    "direct_permutation_family_size"
                ]
                == 8
            ).all()
        )


if __name__ == "__main__":
    unittest.main()
