"""Utilities for Calderbank-Shor-Steane (CSS) stabilizer codes.

A CSS code separates stabilizer generators into two sets:

    X-type checks defined by H_X
    Z-type checks defined by H_Z

where H_X and H_Z are binary matrices.

The essential CSS commutation condition is

    H_X H_Z^T = 0  (mod 2)

which guarantees that every X-type stabilizer commutes with every
Z-type stabilizer.

For Pauli errors written as

    E = (e_x | e_z),

X errors are detected by Z-type checks:

    s_Z = H_Z e_x  (mod 2)

and Z errors are detected by X-type checks:

    s_X = H_X e_z  (mod 2).
"""


from __future__ import annotations

from numpy.typing import ArrayLike

from .gf2 import BinaryArray, as_binary_array, matmul, rank


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


def validate_css_checks(
    h_x: ArrayLike,
    h_z: ArrayLike,
) -> tuple[BinaryArray, BinaryArray]:
    """Validate CSS X-check and Z-check matrices."""
    hx = _as_binary_matrix(
        h_x,
        name="h_x",
    )
    hz = _as_binary_matrix(
        h_z,
        name="h_z",
    )

    if hx.shape[1] != hz.shape[1]:
        raise ValueError(
            "h_x and h_z must act on the same number of qubits"
        )

    return hx, hz


def css_commutation_matrix(
    h_x: ArrayLike,
    h_z: ArrayLike,
) -> BinaryArray:
    """Return H_X H_Z^T over GF(2).

    Entry (i, j) equals 1 exactly when X-type stabilizer i
    anticommutes with Z-type stabilizer j.
    """
    hx, hz = validate_css_checks(
        h_x,
        h_z,
    )

    return matmul(hx, hz.T)


def css_checks_commute(
    h_x: ArrayLike,
    h_z: ArrayLike,
) -> bool:
    """Return True if H_X H_Z^T = 0 over GF(2)."""
    return not css_commutation_matrix(
        h_x,
        h_z,
    ).any()


def number_logical_qubits(
    h_x: ArrayLike,
    h_z: ArrayLike,
) -> int:
    """Return the number of logical qubits encoded by a CSS code.

    For independent CSS stabilizer checks,

        k = n - rank(H_X) - rank(H_Z).
    """
    hx, hz = validate_css_checks(
        h_x,
        h_z,
    )

    if not css_checks_commute(hx, hz):
        raise ValueError(
            "h_x and h_z do not satisfy the CSS commutation condition"
        )

    n = hx.shape[1]

    k = n - rank(hx) - rank(hz)

    if k < 0:
        raise ValueError("invalid CSS dimensions produced k < 0")

    return k


def x_error_syndrome(
    error_x: ArrayLike,
    h_z: ArrayLike,
) -> BinaryArray:
    """Return the Z-check syndrome produced by an X error.

    s_Z = H_Z e_X  (mod 2).
    """
    error = _as_binary_vector(
        error_x,
        name="error_x",
    )

    hz = _as_binary_matrix(
        h_z,
        name="h_z",
    )

    if error.shape[0] != hz.shape[1]:
        raise ValueError(
            "error length must match the number of physical qubits"
        )

    return matmul(hz, error)


def z_error_syndrome(
    error_z: ArrayLike,
    h_x: ArrayLike,
) -> BinaryArray:
    """Return the X-check syndrome produced by a Z error.

    s_X = H_X e_Z  (mod 2).
    """
    error = _as_binary_vector(
        error_z,
        name="error_z",
    )

    hx = _as_binary_matrix(
        h_x,
        name="h_x",
    )

    if error.shape[0] != hx.shape[1]:
        raise ValueError(
            "error length must match the number of physical qubits"
        )

    return matmul(hx, error)


def css_syndrome(
    error_x: ArrayLike,
    error_z: ArrayLike,
    h_x: ArrayLike,
    h_z: ArrayLike,
) -> tuple[BinaryArray, BinaryArray]:
    """Return both components of a CSS syndrome.

    Returns
    -------
    syndrome_x_checks:
        X-stabilizer outcomes. These detect Z errors.

    syndrome_z_checks:
        Z-stabilizer outcomes. These detect X errors.
    """
    hx, hz = validate_css_checks(
        h_x,
        h_z,
    )

    syndrome_x_checks = z_error_syndrome(
        error_z,
        hx,
    )

    syndrome_z_checks = x_error_syndrome(
        error_x,
        hz,
    )

    return syndrome_x_checks, syndrome_z_checks
