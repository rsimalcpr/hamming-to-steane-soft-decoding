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

from numpy.typing import ArrayLike

from .gf2 import BinaryArray, add, as_binary_array, matmul


def _as_binary_matrix(
    values: ArrayLike,
    *,
    name: str,
) -> BinaryArray:
    """Return a validated two-dimensional GF(2) matrix."""
    matrix = as_binary_array(values, name=name)

    if matrix.ndim != 2:
        raise ValueError(f"{name} must be two-dimensional")

    return matrix


def _as_binary_vector(
    values: ArrayLike,
    *,
    name: str,
) -> BinaryArray:
    """Return a validated one-dimensional GF(2) vector."""
    vector = as_binary_array(values, name=name)

    if vector.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")

    return vector


def validate_stabilizers(
    stabilizer_x: ArrayLike,
    stabilizer_z: ArrayLike,
) -> tuple[BinaryArray, BinaryArray]:
    """Validate the X and Z parts of stabilizer generators."""
    sx = _as_binary_matrix(
        stabilizer_x,
        name="stabilizer_x",
    )
    sz = _as_binary_matrix(
        stabilizer_z,
        name="stabilizer_z",
    )

    if sx.shape != sz.shape:
        raise ValueError(
            "stabilizer_x and stabilizer_z must have the same shape"
        )

    return sx, sz


def commutation_matrix(
    stabilizer_x: ArrayLike,
    stabilizer_z: ArrayLike,
) -> BinaryArray:
    """Return all pairwise stabilizer symplectic products.

    For stabilizer generators represented by matrices S_X and S_Z,

        C = S_X S_Z^T + S_Z S_X^T  (mod 2)

    Entry C[i, j] equals

        0 -> generators i and j commute
        1 -> generators i and j anticommute.
    """
    sx, sz = validate_stabilizers(
        stabilizer_x,
        stabilizer_z,
    )

    first_term = matmul(sx, sz.T)
    second_term = matmul(sz, sx.T)

    return add(first_term, second_term)


def stabilizers_commute(
    stabilizer_x: ArrayLike,
    stabilizer_z: ArrayLike,
) -> bool:
    """Return True if every pair of stabilizer generators commutes."""
    matrix = commutation_matrix(
        stabilizer_x,
        stabilizer_z,
    )

    return not matrix.any()


def syndrome(
    error_x: ArrayLike,
    error_z: ArrayLike,
    stabilizer_x: ArrayLike,
    stabilizer_z: ArrayLike,
) -> BinaryArray:
    """Calculate the ideal syndrome of a Pauli error.

    For stabilizers S = (S_X | S_Z) and error E = (e_X | e_Z),

        s = S_X e_Z + S_Z e_X  (mod 2).

    A syndrome bit is

        0 -> error commutes with the corresponding stabilizer
        1 -> error anticommutes with the corresponding stabilizer.
    """
    sx, sz = validate_stabilizers(
        stabilizer_x,
        stabilizer_z,
    )

    ex = _as_binary_vector(
        error_x,
        name="error_x",
    )
    ez = _as_binary_vector(
        error_z,
        name="error_z",
    )

    if ex.shape != ez.shape:
        raise ValueError(
            "error_x and error_z must have the same shape"
        )

    if ex.shape[0] != sx.shape[1]:
        raise ValueError(
            "error length must match the number of physical qubits"
        )

    x_z_term = matmul(sx, ez)
    z_x_term = matmul(sz, ex)

    return add(x_z_term, z_x_term)
