"""Analog measurement model for the six Steane stabilizer syndromes.

An ideal binary syndrome bit is represented by the stabilizer eigenvalue

    0 -> +1
    1 -> -1.

The phenomenological readout model adds independent Gaussian noise to each
of the six values. It isolates the loss caused by hard thresholding without
claiming to model a particular quantum-computing platform.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .gf2 import BinaryArray, as_binary_array

N_STEANE_CHECKS = 6


def syndrome_bits_to_means(
    syndrome_bits: ArrayLike,
) -> NDArray[np.float64]:
    """Map final-axis syndrome bits to ideal analog means in {+1, -1}."""
    bits = as_binary_array(
        syndrome_bits,
        name="syndrome_bits",
    )

    if bits.ndim == 0 or bits.shape[-1] != N_STEANE_CHECKS:
        raise ValueError(
            "syndrome_bits must have final dimension 6"
        )

    return 1.0 - 2.0 * bits.astype(np.float64)


def sample_analog_syndrome(
    syndrome_bits: ArrayLike,
    sigma_m: float,
    rng: np.random.Generator,
) -> NDArray[np.float64]:
    """Sample y = (-1)^s + epsilon with Gaussian readout noise."""
    if not np.isfinite(sigma_m) or sigma_m < 0.0:
        raise ValueError(
            "sigma_m must be finite and non-negative"
        )

    means = syndrome_bits_to_means(
        syndrome_bits
    )

    noise = rng.normal(
        loc=0.0,
        scale=sigma_m,
        size=means.shape,
    )

    return means + noise


def hard_threshold_syndrome(
    observations: ArrayLike,
) -> BinaryArray:
    """Threshold analog stabilizer measurements at zero.

    Non-negative values map to syndrome bit zero.
    Negative values map to syndrome bit one.
    """
    values = np.asarray(
        observations,
        dtype=np.float64,
    )

    if values.ndim == 0 or values.shape[-1] != N_STEANE_CHECKS:
        raise ValueError(
            "observations must have final dimension 6"
        )

    if not np.all(np.isfinite(values)):
        raise ValueError(
            "observations must be finite"
        )

    return (values < 0.0).astype(np.uint8)
