"""Classical linear block codes used by the project."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import gf2

BinaryArray = NDArray[np.uint8]


DEFAULT_G = np.array(
    [
        [1, 0, 0, 0, 1, 1, 0],
        [0, 1, 0, 0, 1, 0, 1],
        [0, 0, 1, 0, 0, 1, 1],
        [0, 0, 0, 1, 1, 1, 1],
    ],
    dtype=np.uint8,
)

DEFAULT_H = np.array(
    [
        [1, 1, 0, 1, 1, 0, 0],
        [1, 0, 1, 1, 0, 1, 0],
        [0, 1, 1, 1, 0, 0, 1],
    ],
    dtype=np.uint8,
)


@dataclass(frozen=True)
class HammingCode74:
    """A systematic Hamming [7,4,3] code with a fixed coordinate order."""

    generator: BinaryArray = field(default_factory=lambda: DEFAULT_G.copy())
    parity_check: BinaryArray = field(default_factory=lambda: DEFAULT_H.copy())

    def __post_init__(self) -> None:
        generator = gf2.as_binary_array(self.generator, name="generator").copy()
        parity_check = gf2.as_binary_array(self.parity_check, name="parity_check").copy()
        gf2.validate_generator_parity_check(generator, parity_check)
        if generator.shape != (4, 7) or parity_check.shape != (3, 7):
            raise ValueError("HammingCode74 requires G shape (4,7) and H shape (3,7)")
        if not np.array_equal(generator[:, :4], np.eye(4, dtype=np.uint8)):
            raise ValueError("generator must be systematic with I4 in its first four columns")
        if len({tuple(column) for column in parity_check.T}) != 7:
            raise ValueError("parity-check columns must be distinct")
        if np.any(np.sum(parity_check, axis=0) == 0):
            raise ValueError("parity-check columns must be nonzero")
        object.__setattr__(self, "generator", generator)
        object.__setattr__(self, "parity_check", parity_check)

    @property
    def n(self) -> int:
        return 7

    @property
    def k(self) -> int:
        return 4

    @property
    def rate(self) -> float:
        return self.k / self.n

    def encode(self, messages: ArrayLike) -> BinaryArray:
        """Encode one message or a batch whose final dimension is four."""
        message_array = gf2.as_binary_array(messages, name="messages")
        if message_array.ndim == 0 or message_array.shape[-1] != self.k:
            raise ValueError("messages must have final dimension 4")
        return gf2.matmul(message_array, self.generator)

    def syndrome(self, received: ArrayLike) -> BinaryArray:
        """Return H r^T, represented with syndrome bits on the final axis."""
        words = gf2.as_binary_array(received, name="received")
        if words.ndim == 0 or words.shape[-1] != self.n:
            raise ValueError("received words must have final dimension 7")
        return gf2.matmul(words, self.parity_check.T)

    def messages(self) -> BinaryArray:
        return gf2.binary_vectors(self.k)

    def codewords(self) -> BinaryArray:
        return self.encode(self.messages())

    def extract_messages(self, codewords: ArrayLike) -> BinaryArray:
        words = gf2.as_binary_array(codewords, name="codewords")
        if words.ndim == 0 or words.shape[-1] != self.n:
            raise ValueError("codewords must have final dimension 7")
        return words[..., : self.k].copy()

    def minimum_distance(self) -> int:
        weights = gf2.hamming_weight(self.codewords())
        return int(np.min(weights[weights > 0]))
