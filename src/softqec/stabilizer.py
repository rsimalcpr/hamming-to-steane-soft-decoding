"""Basic stabilizer-code utilities.

This module builds on the binary symplectic representation from ``pauli.py``.

A stabilizer generator is represented by two binary vectors

    (x | z)

and a collection of m generators acting on n qubits is represented by two
binary matrices of shape (m, n):

    stabilizer_x
    stabilizer_z

For an error E=(e_x | e_z), the syndrome bit associated with stabilizer S_j
is the binary symplectic product

    s_j = <S_j, E>_sp

A syndrome bit is therefore

    0  -> error commutes with the stabilizer
    1  -> error anticommutes with the stabilizer

Equivalently, ideal stabilizer measurement returns

    +1 -> syndrome bit 0
    -1 -> syndrome bit 1
"""

from __future__ import annotations

import numpy as np

from .pauli import symplectic_product


def _as_binary_matrix(matrix: np.ndarray | list[list[int]]) -> np.ndarray:
    """Return a validated two-dimensional binary matrix."""

    arr = np.asarray(matrix, dtype=np.uint8)

    if arr.ndim != 2:
        raise ValueError("Stabilizer representation must be a 2D matrix.")

    if not np.all((arr == 0) | (arr == 1)):
        raise ValueError("Stabilizer matrices must contain only 0 and 1.")

    return arr


def validate_stabilizers(
    stabilizer_x: np.ndarray | list[list[int]],
    stabilizer_z: np.ndarray | list[list[int]],
) -> tuple[np.ndarray, np.ndarray]:
    """Validate X and Z parts of a stabilizer-generator matrix."""

    sx = _as_binary_matrix(stabilizer_x)
    sz = _as_binary_matrix(stabilizer_z)

    if sx.shape != sz.shape:
        raise ValueError(
            "The X and Z stabilizer matrices must have the same shape."
        )

    return sx, sz


def commutation_matrix(
    stabilizer_x: np.ndarray | list[list[int]],
    stabilizer_z: np.ndarray | list[list[int]],
) -> np.ndarray:
    """Return pairwise symplectic products of stabilizer generators.

    Entry (i, j) is

        0 if stabilizers i and j commute
        1 if stabilizers i and j anticommute

    A valid stabilizer group requires this matrix to be all zeros.
    """

    sx, sz = validate_stabilizers(stabilizer_x, stabilizer_z)

    return (sx @ sz.T + sz @ sx.T) % 2


def stabilizers_commute(
    stabilizer_x: np.ndarray | list[list[int]],
    stabilizer_z: np.ndarray | list[list[int]],
) -> bool:
    """Return True if every pair of stabilizer generators commutes."""

    return bool(
        np.all(
            commutation_matrix(stabilizer_x, stabilizer_z) == 0
        )
    )


def syndrome(
    error_x: np.ndarray | list[int],
    error_z: np.ndarray | list[int],
    stabilizer_x: np.ndarray | list[list[int]],
    stabilizer_z: np.ndarray | list[list[int]],
) -> np.ndarray:
    """Calculate the ideal stabilizer syndrome of a Pauli error.

    Parameters
    ----------
    error_x, error_z:
        Binary symplectic representation of the error.

    stabilizer_x, stabilizer_z:
        Binary symplectic representation of the stabilizer generators.

    Returns
    -------
    numpy.ndarray
        Binary syndrome vector. A 1 means that the error anticommutes
        with the corresponding stabilizer generator.
    """

    sx, sz = validate_stabilizers(stabilizer_x, stabilizer_z)

    error_x = np.asarray(error_x, dtype=np.uint8)
    error_z = np.asarray(error_z, dtype=np.uint8)

    if error_x.ndim != 1 or error_z.ndim != 1:
        raise ValueError("Error X and Z parts must be one-dimensional.")

    if len(error_x) != len(error_z):
        raise ValueError("Error X and Z vectors must have equal length.")

    if len(error_x) != sx.shape[1]:
        raise ValueError(
            "Error length must match the number of physical qubits."
        )

    if not np.all((error_x == 0) | (error_x == 1)):
        raise ValueError("Error vectors must contain only 0 and 1.")

    if not np.all((error_z == 0) | (error_z == 1)):
        raise ValueError("Error vectors must contain only 0 and 1.")

    result = np.zeros(sx.shape[0], dtype=np.uint8)

    for j in range(sx.shape[0]):
        result[j] = symplectic_product(
            sx[j],
            sz[j],
            error_x,
            error_z,
        )

    return result
