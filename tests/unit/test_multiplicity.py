from __future__ import annotations

import numpy as np
import pytest

from transport_audit.multiplicity import (
    benjamini_hochberg,
    significance_mask,
)


def test_benjamini_hochberg_known_example() -> None:
    pvalues = np.array(
        [
            0.01,
            0.04,
            0.03,
            0.002,
        ],
        dtype=float,
    )

    observed = benjamini_hochberg(
        pvalues
    )

    expected = np.array(
        [
            0.02,
            0.04,
            0.04,
            0.008,
        ],
        dtype=float,
    )

    np.testing.assert_allclose(
        observed,
        expected,
        rtol=0.0,
        atol=1e-12,
    )


def test_benjamini_hochberg_restores_original_order() -> None:
    sorted_pvalues = np.array(
        [
            0.001,
            0.01,
            0.04,
            0.20,
        ]
    )

    permutation = np.array(
        [
            2,
            0,
            3,
            1,
        ]
    )

    unsorted_pvalues = sorted_pvalues[
        permutation
    ]

    sorted_qvalues = benjamini_hochberg(
        sorted_pvalues
    )

    unsorted_qvalues = benjamini_hochberg(
        unsorted_pvalues
    )

    np.testing.assert_allclose(
        unsorted_qvalues,
        sorted_qvalues[permutation],
        rtol=0.0,
        atol=1e-12,
    )


def test_benjamini_hochberg_preserves_nan_positions() -> None:
    pvalues = np.array(
        [
            0.01,
            np.nan,
            0.20,
            0.04,
        ]
    )

    observed = benjamini_hochberg(
        pvalues
    )

    assert np.isnan(
        observed[1]
    )

    valid_observed = observed[
        [0, 2, 3]
    ]

    expected = benjamini_hochberg(
        np.array(
            [
                0.01,
                0.20,
                0.04,
            ]
        )
    )

    np.testing.assert_allclose(
        valid_observed,
        expected,
        rtol=0.0,
        atol=1e-12,
    )


@pytest.mark.parametrize(
    "invalid_value",
    [
        -0.01,
        1.01,
        np.inf,
        -np.inf,
    ],
)
def test_invalid_pvalues_raise(
    invalid_value: float,
) -> None:
    pvalues = np.array(
        [
            0.01,
            invalid_value,
            0.20,
        ]
    )

    if np.isinf(invalid_value):
        # Infinite values are treated as non-finite/missing.
        result = benjamini_hochberg(
            pvalues
        )

        assert np.isnan(
            result[1]
        )

    else:
        with pytest.raises(
            ValueError
        ):
            benjamini_hochberg(
                pvalues
            )


def test_multidimensional_input_raises() -> None:
    pvalues = np.array(
        [
            [0.01, 0.02],
            [0.03, 0.04],
        ]
    )

    with pytest.raises(
        ValueError
    ):
        benjamini_hochberg(
            pvalues
        )


def test_significance_mask() -> None:
    qvalues = np.array(
        [
            0.01,
            0.049,
            0.05,
            0.20,
            np.nan,
        ]
    )

    observed = significance_mask(
        qvalues,
        alpha=0.05,
    )

    expected = np.array(
        [
            True,
            True,
            False,
            False,
            False,
        ]
    )

    np.testing.assert_array_equal(
        observed,
        expected,
    )


@pytest.mark.parametrize(
    "alpha",
    [
        0.0,
        -0.1,
        1.1,
    ],
)
def test_invalid_alpha_raises(
    alpha: float,
) -> None:
    with pytest.raises(
        ValueError
    ):
        significance_mask(
            [0.01, 0.20],
            alpha=alpha,
        )
