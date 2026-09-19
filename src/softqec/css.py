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

import numpy as np


def _as_binary_matrix(
    matrix: np.ndarray | list[list[int]],
) -> np.ndarray:
    """Return a validated two-dimensional GF(2) matrix."""

    arr = np.asarray(matrix, dtype=np.uint8)

    if arr.ndim != 2:
        raise ValueError("CSS check matrices must be two-dimensional.")

    if not np.all((arr == 0) | (arr == 1)):
        raise ValueError("CSS check matrices must contain only 0 and 1.")

    return arr


def _as_binary_vector(
    vector: np.ndarray | list[int],
) -> np.ndarray:
    """Return a validated one-dimensional GF(2) vector."""

    arr = np.asarray(vector, dtype=np.uint8)

    if arr.ndim != 1:
        raise ValueError("Error vectors must be one-dimensional.")

    if not np.all((arr == 0) | (arr == 1)):
        raise ValueError("Error vectors must contain only 0 and 1.")

    return arr


def gf2_rank(
    matrix: np.ndarray | list[list[int]],
) -> int:
    """Compute matrix rank over GF(2) using Gaussian elimination."""

    a = _as_binary_matrix(matrix).copy()

    n_rows, n_cols = a.shape
    rank = 0

    for col in range(n_cols):
        pivot_candidates = np.flatnonzero(a[rank:, col])

        if len(pivot_candidates) == 0:
            continue

        pivot = rank + pivot_candidates[0]

        if pivot != rank:
            a[[rank, pivot]] = a[[pivot, rank]]

        for row in range(n_rows):
            if row != rank and a[row, col] == 1:
                a[row] ^= a[rank]

        rank += 1

        if rank == n_rows:
            break

    return rank


def css_commutation_matrix(
    h_x: np.ndarray | list[list[int]],
    h_z: np.ndarray | list[list[int]],
) -> np.ndarray:
    """Return H_X H_Z^T over GF(2).

    Entry (i, j) is 1 exactly when X-type stabilizer i
    anticommutes with Z-type stabilizer j.
    """

    hx = _as_binary_matrix(h_x)
    hz = _as_binary_matrix(h_z)

    if hx.shape[1] != hz.shape[1]:
        raise ValueError(
            "H_X and H_Z must act on the same number of qubits."
        )

    return (hx @ hz.T) % 2


def css_checks_commute(
    h_x: np.ndarray | list[list[int]],
    h_z: np.ndarray | list[list[int]],
) -> bool:
    """Return True when the CSS commutation condition is satisfied."""

    return bool(
        np.all(
            css_commutation_matrix(h_x, h_z) == 0
        )
    )


def number_logical_qubits(
    h_x: np.ndarray | list[list[int]],
    h_z: np.ndarray | list[list[int]],
) -> int:
    """Return the number k of encoded logical qubits.

    For a valid CSS stabilizer code,

        k = n - rank(H_X) - rank(H_Z).
    """

    hx = _as_binary_matrix(h_x)
    hz = _as_binary_matrix(h_z)

    if hx.shape[1] != hz.shape[1]:
        raise ValueError(
            "H_X and H_Z must act on the same number of qubits."
        )

    if not css_checks_commute(hx, hz):
        raise ValueError(
            "H_X and H_Z do not satisfy the CSS commutation condition."
        )

    n = hx.shape[1]

    k = n - gf2_rank(hx) - gf2_rank(hz)

    if k < 0:
        raise ValueError("Invalid CSS dimensions produced k < 0.")

    return int(k)


def x_error_syndrome(
    error_x: np.ndarray | list[int],
    h_z: np.ndarray | list[list[int]],
) -> np.ndarray:
    """Return the syndrome of an X-type error using Z checks."""

    error = _as_binary_vector(error_x)
    hz = _as_binary_matrix(h_z)

    if len(error) != hz.shape[1]:
        raise ValueError(
            "Error length must equal the number of physical qubits."
        )

    return (hz @ error) % 2


def z_error_syndrome(
    error_z: np.ndarray | list[int],
    h_x: np.ndarray | list[list[int]],
) -> np.ndarray:
    """Return the syndrome of a Z-type error using X checks."""

    error = _as_binary_vector(error_z)
    hx = _as_binary_matrix(h_x)

    if len(error) != hx.shape[1]:
        raise ValueError(
            "Error length must equal the number of physical qubits."
        )

    return (hx @ error) % 2


def css_syndrome(
    error_x: np.ndarray | list[int],
    error_z: np.ndarray | list[int],
    h_x: np.ndarray | list[list[int]],
    h_z: np.ndarray | list[list[int]],
) -> tuple[np.ndarray, np.ndarray]:
    """Return both CSS syndrome components.

    Returns
    -------
    syndrome_x_checks:
        Outcomes of X-type stabilizers. These detect Z errors.

    syndrome_z_checks:
        Outcomes of Z-type stabilizers. These detect X errors.
    """

    syndrome_x_checks = z_error_syndrome(error_z, h_x)
    syndrome_z_checks = x_error_syndrome(error_x, h_z)

    return syndrome_x_checks, syndrome_z_checks
