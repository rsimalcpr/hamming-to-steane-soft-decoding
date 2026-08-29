"""Classical channels and BPSK modulation helpers."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .gf2 import BinaryArray, as_binary_array


def binary_symmetric_channel(
    bits: ArrayLike,
    flip_probability: float,
    rng: np.random.Generator,
) -> BinaryArray:
    """Flip each bit independently with probability p."""
    if not 0.0 <= flip_probability <= 1.0:
        raise ValueError("flip_probability must lie in [0,1]")
    binary = as_binary_array(bits, name="bits")
    flips = rng.random(binary.shape) < flip_probability
    return np.bitwise_xor(binary, flips.astype(np.uint8))


def bpsk_modulate(bits: ArrayLike) -> NDArray[np.float64]:
    """Map 0 to +1 and 1 to -1."""
    binary = as_binary_array(bits, name="bits")
    return 1.0 - 2.0 * binary.astype(np.float64)


def hard_demodulate(observations: ArrayLike) -> BinaryArray:
    """Threshold BPSK observations at zero; ties map to bit zero."""
    values = np.asarray(observations, dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise ValueError("observations must be finite")
    return (values < 0.0).astype(np.uint8)


def noise_sigma_from_ebn0_db(ebn0_db: float, rate: float) -> float:
    """Return real-AWGN sigma for unit symbol energy and information-bit Eb/N0."""
    if not 0.0 < rate <= 1.0:
        raise ValueError("rate must lie in (0,1]")
    gamma_b = 10.0 ** (float(ebn0_db) / 10.0)
    return float(np.sqrt(1.0 / (2.0 * rate * gamma_b)))


def transmit_bpsk_awgn(
    bits: ArrayLike,
    ebn0_db: float,
    rate: float,
    rng: np.random.Generator,
) -> NDArray[np.float64]:
    """Modulate bits and add independent zero-mean Gaussian noise."""
    symbols = bpsk_modulate(bits)
    sigma = noise_sigma_from_ebn0_db(ebn0_db, rate)
    return symbols + rng.normal(loc=0.0, scale=sigma, size=symbols.shape)
