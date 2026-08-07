"""Multiple-testing correction utilities."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np


def _as_float_array(
    values: Sequence[float] | np.ndarray,
) -> np.ndarray:
    """Convert p-values to a one-dimensional float array."""

    array = np.asarray(values, dtype=float)

    if array.ndim != 1:
        raise ValueError(
            "P-values must be a one-dimensional sequence."
        )

    return array


def validate_pvalues(
    values: Sequence[float] | np.ndarray,
) -> np.ndarray:
    """
    Validate p-values while allowing missing values.

    Finite values must lie in the closed interval [0, 1].
    NaN values are preserved and excluded from correction.
    """

    pvalues = _as_float_array(values)

    finite = np.isfinite(pvalues)

    invalid = finite & (
        (pvalues < 0.0)
        | (pvalues > 1.0)
    )

    if np.any(invalid):
        bad_values = pvalues[invalid]

        raise ValueError(
            "Finite p-values must lie within [0, 1]. "
            f"Invalid values: {bad_values.tolist()}"
        )

    return pvalues


def benjamini_hochberg(
    values: Sequence[float] | np.ndarray,
) -> np.ndarray:
    """
    Adjust p-values using the Benjamini-Hochberg procedure.

    Missing values are excluded from the correction and restored
    as NaN in the returned array.

    Parameters
    ----------
    values:
        One-dimensional sequence of raw p-values.

    Returns
    -------
    numpy.ndarray
        BH-adjusted q-values in the same order as the input.

    Notes
    -----
    The implementation reproduces the procedure used in the
    frozen osteosarcoma analysis:

    1. remove non-finite values from the multiplicity family;
    2. sort valid p-values in ascending order;
    3. multiply by m / rank;
    4. enforce monotonicity from largest to smallest rank;
    5. clip adjusted values to [0, 1];
    6. restore the original input order.
    """

    pvalues = validate_pvalues(values)

    adjusted_full = np.full(
        pvalues.shape,
        np.nan,
        dtype=float,
    )

    valid = np.isfinite(pvalues)

    if not np.any(valid):
        return adjusted_full

    valid_pvalues = pvalues[valid]

    order = np.argsort(
        valid_pvalues,
        kind="stable",
    )

    ranked = valid_pvalues[order]

    n_tests = ranked.size

    ranks = np.arange(
        1,
        n_tests + 1,
        dtype=float,
    )

    adjusted_ranked = (
        ranked
        * n_tests
        / ranks
    )

    adjusted_ranked = np.minimum.accumulate(
        adjusted_ranked[::-1]
    )[::-1]

    adjusted_ranked = np.clip(
        adjusted_ranked,
        0.0,
        1.0,
    )

    restored = np.empty(
        n_tests,
        dtype=float,
    )

    restored[order] = adjusted_ranked

    adjusted_full[valid] = restored

    return adjusted_full


def significance_mask(
    adjusted_pvalues: Sequence[float] | np.ndarray,
    alpha: float = 0.05,
) -> np.ndarray:
    """
    Return a boolean FDR-support mask.

    Missing values are always classified as False.
    """

    if not 0.0 < alpha <= 1.0:
        raise ValueError(
            "alpha must satisfy 0 < alpha <= 1."
        )

    qvalues = validate_pvalues(
        adjusted_pvalues
    )

    return (
        np.isfinite(qvalues)
        & (qvalues < alpha)
    )
