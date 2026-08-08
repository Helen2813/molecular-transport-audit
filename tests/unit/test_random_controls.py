"""Unit tests for random-panel specificity controls."""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from transport_audit.random_controls import (
    LEGACY_Z_VARIANCE,
    PRESTANDARDIZATION_VARIANCE,
    classify_external_representation,
    draw_matched_panel,
    random_panel_empirical_p,
    variability_bins,
)


class EmpiricalPTests(
    unittest.TestCase
):
    def test_plus_one_correction(
        self,
    ) -> None:
        null = np.array(
            [
                0.1,
                0.2,
                0.3,
            ]
        )

        observed = 0.5

        self.assertEqual(
            random_panel_empirical_p(
                null,
                observed,
            ),
            0.25,
        )

    def test_two_sided_absolute_statistic(
        self,
    ) -> None:
        null = np.array(
            [
                -0.9,
                0.2,
                0.4,
            ]
        )

        observed = 0.8

        self.assertEqual(
            random_panel_empirical_p(
                null,
                observed,
            ),
            0.5,
        )


class MatchedPanelTests(
    unittest.TestCase
):
    def test_same_bin_and_without_replacement(
        self,
    ) -> None:
        candidates = [
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
        ]

        bins = pd.Series(
            {
                "A": 0,
                "B": 0,
                "C": 0,
                "D": 1,
                "E": 1,
                "F": 1,
            }
        )

        rng = np.random.default_rng(
            42
        )

        panel, fallbacks = (
            draw_matched_panel(
                candidate_genes=(
                    candidates
                ),
                bins=bins,
                target_bins=[
                    0,
                    0,
                    1,
                    1,
                ],
                rng=rng,
            )
        )

        self.assertEqual(
            len(panel),
            4,
        )

        self.assertEqual(
            len(
                set(panel)
            ),
            4,
        )

        self.assertEqual(
            fallbacks,
            0,
        )

        self.assertEqual(
            [
                bins.loc[gene]
                for gene
                in panel
            ],
            [
                0,
                0,
                1,
                1,
            ],
        )

    def test_fallback_is_counted(
        self,
    ) -> None:
        candidates = [
            "A",
            "B",
            "C",
        ]

        bins = pd.Series(
            {
                "A": 0,
                "B": 1,
                "C": 1,
            }
        )

        rng = np.random.default_rng(
            5
        )

        panel, fallbacks = (
            draw_matched_panel(
                candidate_genes=(
                    candidates
                ),
                bins=bins,
                target_bins=[
                    0,
                    0,
                ],
                rng=rng,
            )
        )

        self.assertEqual(
            len(panel),
            2,
        )

        self.assertEqual(
            fallbacks,
            1,
        )


class VariabilityModeTests(
    unittest.TestCase
):
    def test_corrected_mode_uses_prestandardization_variance(
        self,
    ) -> None:
        reference = pd.DataFrame(
            {
                "A": [
                    0.0,
                    1.0,
                    0.0,
                    1.0,
                ],
                "B": [
                    0.0,
                    10.0,
                    0.0,
                    10.0,
                ],
                "C": [
                    0.0,
                    100.0,
                    0.0,
                    100.0,
                ],
                "D": [
                    0.0,
                    1000.0,
                    0.0,
                    1000.0,
                ],
            }
        )

        external = reference.copy()

        corrected = variability_bins(
            reference,
            external,
            matching_mode=(
                PRESTANDARDIZATION_VARIANCE
            ),
            n_bins=2,
        )

        self.assertEqual(
            corrected.loc[
                "A"
            ],
            0,
        )

        self.assertEqual(
            corrected.loc[
                "D"
            ],
            1,
        )

        legacy = variability_bins(
            reference,
            external,
            matching_mode=(
                LEGACY_Z_VARIANCE
            ),
            n_bins=2,
        )

        self.assertEqual(
            legacy.shape[0],
            4,
        )


class ClassificationTests(
    unittest.TestCase
):
    def test_strong_class_does_not_require_random_panel_support(
        self,
    ) -> None:
        structure = pd.DataFrame(
            {
                "module_label": [
                    "MTEST"
                ],
                "edge_permutation_p": [
                    0.001
                ],
                "loading_permutation_p": [
                    0.001
                ],
                "split_half_median": [
                    0.90
                ],
            }
        )

        random_controls = pd.DataFrame(
            {
                "module_label": [
                    "MTEST"
                ],
                "random_panel_empirical_p": [
                    0.80
                ],
            }
        )

        observed = (
            classify_external_representation(
                structure,
                random_controls,
            )
        )

        self.assertEqual(
            observed.loc[
                0,
                (
                    "external_canine_"
                    "representation_class"
                ),
            ],
            (
                "strong_external_canine_"
                "representation_preservation"
            ),
        )


if __name__ == "__main__":
    unittest.main()
