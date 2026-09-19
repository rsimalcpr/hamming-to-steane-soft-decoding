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


def _as_binary_vector(v: np.ndarray | list[int]) -> np.ndarray:
    """Return a validated one-dimensional GF(2) vector."""

    arr = np.asarray(v, dtype=np.uint8)

    if arr.ndim != 1:
        raise ValueError("Pauli vectors must be one-dimensional.")

    if not np.all((arr == 0) | (arr == 1)):
        raise ValueError("Pauli vectors must contain only 0 and 1.")

    return arr


def pauli_string_to_symplectic(pauli: str) -> tuple[np.ndarray, np.ndarray]:
    """Convert a Pauli string to binary symplectic form.

    Parameters
    ----------
    pauli:
        String containing only I, X, Y and Z.

    Returns
    -------
    x, z:
        Binary vectors describing the X and Z components.

    Examples
    --------
    XZI -> x = [1, 0, 0], z = [0, 1, 0]
    YII -> x = [1, 0, 0], z = [1, 0, 0]
    """

    pauli = pauli.upper()

    if not pauli:
        raise ValueError("Pauli string cannot be empty.")

    allowed = {"I", "X", "Y", "Z"}

    if any(p not in allowed for p in pauli):
        raise ValueError("Pauli string may contain only I, X, Y and Z.")

    x = np.zeros(len(pauli), dtype=np.uint8)
    z = np.zeros(len(pauli), dtype=np.uint8)

    for i, operator in enumerate(pauli):
        if operator == "X":
            x[i] = 1
        elif operator == "Z":
            z[i] = 1
        elif operator == "Y":
            x[i] = 1
            z[i] = 1

    return x, z


def symplectic_product(
    x1: np.ndarray | list[int],
    z1: np.ndarray | list[int],
    x2: np.ndarray | list[int],
    z2: np.ndarray | list[int],
) -> int:
    """Compute the binary symplectic inner product.

    For Pauli operators P1=(x1|z1) and P2=(x2|z2),

        <P1, P2> = x1·z2 + z1·x2  (mod 2)

    A result of 0 means that the operators commute.
    A result of 1 means that they anticommute.
    """

    x1 = _as_binary_vector(x1)
    z1 = _as_binary_vector(z1)
    x2 = _as_binary_vector(x2)
    z2 = _as_binary_vector(z2)

    if not (len(x1) == len(z1) == len(x2) == len(z2)):
        raise ValueError("All Pauli vectors must have the same length.")

    value = (np.dot(x1, z2) + np.dot(z1, x2)) % 2

    return int(value)


def commutes(
    x1: np.ndarray | list[int],
    z1: np.ndarray | list[int],
    x2: np.ndarray | list[int],
    z2: np.ndarray | list[int],
) -> bool:
    """Return True if two Pauli operators commute."""

    return symplectic_product(x1, z1, x2, z2) == 0


def pauli_weight(
    x: np.ndarray | list[int],
    z: np.ndarray | list[int],
) -> int:
    """Return the number of qubits on which the Pauli is non-identity."""

    x = _as_binary_vector(x)
    z = _as_binary_vector(z)

    if len(x) != len(z):
        raise ValueError("X and Z vectors must have the same length.")

    return int(np.count_nonzero(x | z))
