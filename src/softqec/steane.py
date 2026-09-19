"""The Steane [[7,1,3]] CSS quantum error-correcting code.

The Steane code is constructed from the parity-check matrix of the
classical Hamming [7,4,3] code.

For this CSS construction,

    H_X = H_Z = H,

with H chosen so that

    H H^T = 0  (mod 2).

The code therefore has

    n = 7
    k = 7 - rank(H_X) - rank(H_Z) = 1
    d = 3

and can correct any single-qubit Pauli error.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from .css import (
    css_checks_commute,
    css_syndrome,
    number_logical_qubits,
)
from .gf2 import BinaryArray, add, as_binary_array, rank


# ---------------------------------------------------------------------------
# Steane CSS construction
# ---------------------------------------------------------------------------

STEANE_H = np.array(
    [
        [1, 0, 1, 0, 1, 0, 1],
        [0, 1, 1, 0, 0, 1, 1],
        [0, 0, 0, 1, 1, 1, 1],
    ],
    dtype=np.uint8,
)

STEANE_H_X = STEANE_H.copy()
STEANE_H_Z = STEANE_H.copy()

N_PHYSICAL = 7
N_LOGICAL = 1
DISTANCE = 3


# Logical Pauli operators.
#
# X_bar = X X X X X X X
# Z_bar = Z Z Z Z Z Z Z

LOGICAL_X = np.ones(N_PHYSICAL, dtype=np.uint8)
LOGICAL_Z = np.ones(N_PHYSICAL, dtype=np.uint8)


def steane_stabilizers() -> tuple[BinaryArray, BinaryArray]:
    """Return the six Steane stabilizer generators in symplectic form.

    The first three rows are X-type stabilizers.
    The final three rows are Z-type stabilizers.

    Returns
    -------
    stabilizer_x, stabilizer_z
        Two 6 x 7 binary matrices.
    """
    zeros = np.zeros_like(STEANE_H)

    stabilizer_x = np.vstack(
        [
            STEANE_H_X,
            zeros,
        ]
    )

    stabilizer_z = np.vstack(
        [
            zeros,
            STEANE_H_Z,
        ]
    )

    return stabilizer_x, stabilizer_z


def steane_parameters() -> tuple[int, int, int]:
    """Return the Steane code parameters (n, k, d)."""
    return (
        N_PHYSICAL,
        number_logical_qubits(
            STEANE_H_X,
            STEANE_H_Z,
        ),
        DISTANCE,
    )


def _as_error_vector(
    values: ArrayLike,
    *,
    name: str,
) -> BinaryArray:
    """Return a validated seven-qubit binary error vector."""
    vector = as_binary_array(
        values,
        name=name,
    )

    if vector.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")

    if vector.shape[0] != N_PHYSICAL:
        raise ValueError(
            f"{name} must contain exactly {N_PHYSICAL} entries"
        )

    return vector


def single_qubit_lookup() -> dict[tuple[int, ...], BinaryArray]:
    """Build the ideal single-qubit syndrome lookup table.

    Each nonzero syndrome of the Hamming parity-check matrix uniquely
    identifies one of the seven physical qubits.
    """
    lookup: dict[tuple[int, ...], BinaryArray] = {}

    for qubit in range(N_PHYSICAL):
        syndrome = tuple(
            int(value)
            for value in STEANE_H[:, qubit]
        )

        error = np.zeros(
            N_PHYSICAL,
            dtype=np.uint8,
        )
        error[qubit] = 1

        lookup[syndrome] = error

    return lookup


def decode_single_qubit_syndrome(
    syndrome: ArrayLike,
) -> BinaryArray:
    """Return the minimum-weight correction for an ideal syndrome.

    The zero syndrome returns the identity correction.

    A nonzero syndrome is matched to the corresponding single-qubit
    correction.
    """
    s = as_binary_array(
        syndrome,
        name="syndrome",
    )

    if s.ndim != 1 or s.shape[0] != 3:
        raise ValueError(
            "Steane syndrome must contain exactly three bits"
        )

    if not s.any():
        return np.zeros(
            N_PHYSICAL,
            dtype=np.uint8,
        )

    lookup = single_qubit_lookup()

    key = tuple(int(value) for value in s)

    try:
        return lookup[key].copy()
    except KeyError as exc:
        raise ValueError(
            f"unknown Steane syndrome: {key}"
        ) from exc


def decode_ideal_css(
    error_x: ArrayLike,
    error_z: ArrayLike,
) -> tuple[BinaryArray, BinaryArray]:
    """Decode an ideal Steane syndrome.

    X errors are detected by Z-type stabilizers.
    Z errors are detected by X-type stabilizers.

    Returns
    -------
    correction_x, correction_z
        Binary correction vectors.
    """
    ex = _as_error_vector(
        error_x,
        name="error_x",
    )

    ez = _as_error_vector(
        error_z,
        name="error_z",
    )

    syndrome_x_checks, syndrome_z_checks = css_syndrome(
        ex,
        ez,
        STEANE_H_X,
        STEANE_H_Z,
    )

    correction_x = decode_single_qubit_syndrome(
        syndrome_z_checks
    )

    correction_z = decode_single_qubit_syndrome(
        syndrome_x_checks
    )

    return correction_x, correction_z


def in_row_space(
    vector: ArrayLike,
    matrix: ArrayLike,
) -> bool:
    """Return True if a binary vector belongs to a matrix row space."""
    v = as_binary_array(
        vector,
        name="vector",
    )

    m = as_binary_array(
        matrix,
        name="matrix",
    )

    if v.ndim != 1:
        raise ValueError(
            "vector must be one-dimensional"
        )

    if m.ndim != 2:
        raise ValueError(
            "matrix must be two-dimensional"
        )

    if v.shape[0] != m.shape[1]:
        raise ValueError(
            "vector length must equal matrix column count"
        )

    augmented = np.vstack(
        [
            m,
            v,
        ]
    )

    return rank(augmented) == rank(m)


def is_stabilizer(
    pauli_x: ArrayLike,
    pauli_z: ArrayLike,
) -> bool:
    """Return True if a Pauli belongs to the Steane stabilizer group.

    For a CSS stabilizer group, the X component must belong to the
    row space of H_X and the Z component must belong to the row space
    of H_Z.
    """
    x = _as_error_vector(
        pauli_x,
        name="pauli_x",
    )

    z = _as_error_vector(
        pauli_z,
        name="pauli_z",
    )

    return (
        in_row_space(x, STEANE_H_X)
        and in_row_space(z, STEANE_H_Z)
    )


def correction_succeeds(
    error_x: ArrayLike,
    error_z: ArrayLike,
    correction_x: ArrayLike,
    correction_z: ArrayLike,
) -> bool:
    """Return True if error and correction differ only by a stabilizer.

    Success does not require the residual Pauli to be the identity.

    It is sufficient that

        E C ∈ S,

    because stabilizers act trivially on the logical state.
    """
    ex = _as_error_vector(
        error_x,
        name="error_x",
    )

    ez = _as_error_vector(
        error_z,
        name="error_z",
    )

    cx = _as_error_vector(
        correction_x,
        name="correction_x",
    )

    cz = _as_error_vector(
        correction_z,
        name="correction_z",
    )

    residual_x = add(ex, cx)
    residual_z = add(ez, cz)

    return is_stabilizer(
        residual_x,
        residual_z,
    )


def validate_steane_construction() -> None:
    """Raise ValueError if basic Steane-code identities fail."""

    if not css_checks_commute(
        STEANE_H_X,
        STEANE_H_Z,
    ):
        raise ValueError(
            "Steane X and Z checks do not commute"
        )

    if rank(STEANE_H_X) != 3:
        raise ValueError(
            "Steane H_X must have rank 3"
        )

    if rank(STEANE_H_Z) != 3:
        raise ValueError(
            "Steane H_Z must have rank 3"
        )

    if number_logical_qubits(
        STEANE_H_X,
        STEANE_H_Z,
    ) != 1:
        raise ValueError(
            "Steane code must encode exactly one logical qubit"
        )
