"""Hard syndrome and exact soft-ML decoders for Hamming [7,4,3]."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from .classical_channels import bpsk_modulate
from .classical_codes import HammingCode74
from .gf2 import BinaryArray, as_binary_array


def syndrome_correction_table(code: HammingCode74) -> BinaryArray:
    """Map each integer-valued three-bit syndrome to a correction vector."""
    table = np.zeros((2 ** (code.n - code.k), code.n), dtype=np.uint8)
    powers = 2 ** np.arange(code.n - code.k - 1, -1, -1)
    for bit_index, column in enumerate(code.parity_check.T):
        syndrome_index = int(column @ powers)
        table[syndrome_index, bit_index] = 1
    return table


def syndrome_to_index(syndromes: ArrayLike) -> np.ndarray:
    """Convert final-axis binary syndrome vectors to integer lookup indices."""
    values = as_binary_array(syndromes, name="syndromes")
    if values.ndim == 0:
        raise ValueError("syndromes must have at least one dimension")
    powers = 2 ** np.arange(values.shape[-1] - 1, -1, -1)
    return values.astype(np.int64) @ powers


def hard_syndrome_decode(
    received: ArrayLike,
    code: HammingCode74,
) -> tuple[BinaryArray, BinaryArray]:
    """Apply minimum-weight single-error syndrome lookup decoding."""
    words = as_binary_array(received, name="received")
    if words.ndim == 0 or words.shape[-1] != code.n:
        raise ValueError("received words must have final dimension 7")
    corrections = syndrome_correction_table(code)[syndrome_to_index(code.syndrome(words))]
    corrected = np.bitwise_xor(words, corrections)
    return corrected, corrections


def exact_soft_ml_decode(
    observations: ArrayLike,
    code: HammingCode74,
    *,
    chunk_size: int = 50_000,
) -> BinaryArray:
    """Decode AWGN observations by exact ML search over all 16 codewords."""
    values = np.asarray(observations, dtype=np.float64)
    if values.ndim == 1:
        values = values[None, :]
        squeeze = True
    else:
        squeeze = False
    if values.ndim != 2 or values.shape[1] != code.n:
        raise ValueError("observations must have shape (7,) or (frames,7)")
    if not np.all(np.isfinite(values)):
        raise ValueError("observations must be finite")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")

    candidates = code.codewords()
    candidate_symbols = bpsk_modulate(candidates)
    decoded = np.empty((values.shape[0], code.n), dtype=np.uint8)
    for start in range(0, values.shape[0], chunk_size):
        stop = min(start + chunk_size, values.shape[0])
        scores = values[start:stop] @ candidate_symbols.T
        decoded[start:stop] = candidates[np.argmax(scores, axis=1)]
    return decoded[0] if squeeze else decoded
