"""Restricted hard and soft decoders for noisy Steane syndromes.

Week 3 intentionally considers only the 22 equally likely hypotheses

    {I, X_1, ..., X_7, Y_1, ..., Y_7, Z_1, ..., Z_7}.

This small hypothesis space makes the hard-versus-soft information question
transparent. It is not a general multi-qubit or degeneracy-aware decoder.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .analog_syndrome import (
    hard_threshold_syndrome,
    syndrome_bits_to_means,
)
from .css import css_syndrome
from .gf2 import BinaryArray
from .steane import (
    N_PHYSICAL,
    STEANE_H_X,
    STEANE_H_Z,
)


def single_pauli_hypotheses(
) -> tuple[tuple[str, ...], BinaryArray, BinaryArray]:
    """Return labels and binary X/Z components for the 22 hypotheses."""
    labels = (
        "I",
        *(
            f"X{qubit}"
            for qubit in range(1, N_PHYSICAL + 1)
        ),
        *(
            f"Y{qubit}"
            for qubit in range(1, N_PHYSICAL + 1)
        ),
        *(
            f"Z{qubit}"
            for qubit in range(1, N_PHYSICAL + 1)
        ),
    )

    error_x = np.zeros(
        (len(labels), N_PHYSICAL),
        dtype=np.uint8,
    )

    error_z = np.zeros_like(error_x)

    for qubit in range(N_PHYSICAL):
        # X_i has only an X component.
        error_x[1 + qubit, qubit] = 1

        # Y_i has both X and Z components.
        error_x[
            1 + N_PHYSICAL + qubit,
            qubit,
        ] = 1

        error_z[
            1 + N_PHYSICAL + qubit,
            qubit,
        ] = 1

        # Z_i has only a Z component.
        error_z[
            1 + 2 * N_PHYSICAL + qubit,
            qubit,
        ] = 1

    return labels, error_x, error_z


def restricted_syndrome_table() -> BinaryArray:
    """Return the unique six-bit syndrome of every restricted hypothesis.

    Syndrome order is

        [X-check responses | Z-check responses].

    X checks detect Z components.
    Z checks detect X components.
    """
    _, errors_x, errors_z = single_pauli_hypotheses()

    rows = []

    for error_x, error_z in zip(
        errors_x,
        errors_z,
        strict=True,
    ):
        syndrome_x, syndrome_z = css_syndrome(
            error_x,
            error_z,
            STEANE_H_X,
            STEANE_H_Z,
        )

        rows.append(
            np.concatenate(
                [syndrome_x, syndrome_z]
            )
        )

    return np.asarray(
        rows,
        dtype=np.uint8,
    )


def ideal_analog_signatures() -> NDArray[np.float64]:
    """Return the 22 ideal +1/-1 analog syndrome signatures."""
    return syndrome_bits_to_means(
        restricted_syndrome_table()
    )


def _as_observation_matrix(
    observations: ArrayLike,
) -> tuple[NDArray[np.float64], bool]:
    """Validate observations and convert one sample into a row matrix."""
    values = np.asarray(
        observations,
        dtype=np.float64,
    )

    squeeze = values.ndim == 1

    if squeeze:
        values = values[None, :]

    if values.ndim != 2 or values.shape[1] != 6:
        raise ValueError(
            "observations must have shape (6,) or (trials,6)"
        )

    if not np.all(np.isfinite(values)):
        raise ValueError(
            "observations must be finite"
        )

    return values, squeeze


def hard_nearest_syndrome_decode(
    observations: ArrayLike,
) -> int | NDArray[np.int64]:
    """Threshold first and choose the nearest valid binary syndrome.

    The distance between the thresholded observation and every valid
    syndrome signature is measured using Hamming distance.

    numpy.argmin supplies a fixed first-index tie break.

    The returned indices refer to the ordering produced by
    single_pauli_hypotheses().
    """
    values, squeeze = _as_observation_matrix(
        observations
    )

    hard_bits = hard_threshold_syndrome(
        values
    )

    candidates = restricted_syndrome_table()

    distances = np.count_nonzero(
        hard_bits[:, None, :]
        != candidates[None, :, :],
        axis=2,
    )

    decoded = np.argmin(
        distances,
        axis=1,
    )

    if squeeze:
        return int(decoded[0])

    return decoded


def soft_ml_syndrome_decode(
    observations: ArrayLike,
) -> int | NDArray[np.int64]:
    """Perform Gaussian ML on the continuous syndrome observations.

    With equal hypothesis priors and equal independent Gaussian variance,
    maximizing the likelihood is equivalent to minimizing the squared
    Euclidean distance to an ideal analog syndrome signature.
    """
    values, squeeze = _as_observation_matrix(
        observations
    )

    candidates = ideal_analog_signatures()

    squared_distances = np.sum(
        (
            values[:, None, :]
            - candidates[None, :, :]
        )
        ** 2,
        axis=2,
    )

    decoded = np.argmin(
        squared_distances,
        axis=1,
    )

    if squeeze:
        return int(decoded[0])

    return decoded


def corrections_from_indices(
    indices: ArrayLike,
) -> tuple[BinaryArray, BinaryArray]:
    """Return X/Z correction components selected by hypothesis indices."""
    values = np.asarray(indices)

    if not np.issubdtype(
        values.dtype,
        np.integer,
    ):
        raise ValueError(
            "indices must be integers"
        )

    labels, corrections_x, corrections_z = (
        single_pauli_hypotheses()
    )

    if (
        np.any(values < 0)
        or np.any(values >= len(labels))
    ):
        raise ValueError(
            "hypothesis index out of range"
        )

    return (
        corrections_x[values].copy(),
        corrections_z[values].copy(),
    )


def restricted_failure_count(
    actual: ArrayLike,
    decoded: ArrayLike,
) -> int:
    """Count correction failures inside the restricted hypothesis model.

    Any two distinct allowed hypotheses differ by a Pauli operator of
    weight at most two. Every non-identity Steane stabilizer has weight
    four.

    Therefore, within this restricted set, correction success modulo
    stabilizers is equivalent to exact hypothesis identification.
    """
    reference = np.asarray(actual)
    estimate = np.asarray(decoded)

    if reference.shape != estimate.shape:
        raise ValueError(
            f"shape mismatch: "
            f"{reference.shape} != {estimate.shape}"
        )

    return int(
        np.count_nonzero(
            reference != estimate
        )
    )
