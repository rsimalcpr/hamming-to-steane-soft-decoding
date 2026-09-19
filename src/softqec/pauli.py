"""Pauli operators in binary symplectic form.

This module provides the minimal algebra needed for stabilizer and CSS codes.

An n-qubit Pauli operator is represented by two binary vectors

    (x | z)

where, for each qubit,

    I -> (0, 0)
    X -> (1, 0)
    Z -> (0, 1)
    Y -> (1, 1)

Global phases (+1, -1, +i, -i) are intentionally ignored because they do not
affect commutation relations or syndrome calculations used in this project.
"""


from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .gf2 import as_binary_array

BinaryArray = NDArray[np.uint8]


def _as_binary_vector(values: ArrayLike, *, name: str) -> BinaryArray:
    """Return a validated one-dimensional GF(2) vector."""
    vector = as_binary_array(values, name=name)

    if vector.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")

    return vector


def pauli_string_to_symplectic(
    pauli: str,
) -> tuple[BinaryArray, BinaryArray]:
    """Convert a Pauli string to binary symplectic form.

    The mapping is

        I -> (0, 0)
        X -> (1, 0)
        Z -> (0, 1)
        Y -> (1, 1)

    Global phases are ignored.
    """
    pauli = pauli.upper()

    if not pauli:
        raise ValueError("Pauli string cannot be empty")

    allowed = {"I", "X", "Y", "Z"}

    if any(operator not in allowed for operator in pauli):
        raise ValueError("Pauli string may contain only I, X, Y and Z")

    x = np.zeros(len(pauli), dtype=np.uint8)
    z = np.zeros(len(pauli), dtype=np.uint8)

    for index, operator in enumerate(pauli):
        if operator == "X":
            x[index] = 1

        elif operator == "Z":
            z[index] = 1

        elif operator == "Y":
            x[index] = 1
            z[index] = 1

    return x, z


def symplectic_product(
    x1: ArrayLike,
    z1: ArrayLike,
    x2: ArrayLike,
    z2: ArrayLike,
) -> int:
    """Return the binary symplectic inner product.

    For P1 = (x1 | z1) and P2 = (x2 | z2),

        <P1, P2> = x1·z2 + z1·x2  (mod 2)

    0 means that the Pauli operators commute.
    1 means that they anticommute.
    """
    x1 = _as_binary_vector(x1, name="x1")
    z1 = _as_binary_vector(z1, name="z1")
    x2 = _as_binary_vector(x2, name="x2")
    z2 = _as_binary_vector(z2, name="z2")

    if not (x1.shape == z1.shape == x2.shape == z2.shape):
        raise ValueError("all Pauli vectors must have the same shape")

    value = (
        int(np.dot(x1.astype(np.int64), z2.astype(np.int64)))
        + int(np.dot(z1.astype(np.int64), x2.astype(np.int64)))
    ) % 2

    return value


def commutes(
    x1: ArrayLike,
    z1: ArrayLike,
    x2: ArrayLike,
    z2: ArrayLike,
) -> bool:
    """Return True if two Pauli operators commute."""
    return symplectic_product(x1, z1, x2, z2) == 0


def pauli_weight(
    x: ArrayLike,
    z: ArrayLike,
) -> int:
    """Return the number of non-identity qubit positions."""
    x = _as_binary_vector(x, name="x")
    z = _as_binary_vector(z, name="z")

    if x.shape != z.shape:
        raise ValueError("x and z must have the same shape")

    return int(np.count_nonzero(np.bitwise_or(x, z)))
