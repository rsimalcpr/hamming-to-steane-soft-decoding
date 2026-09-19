import numpy as np

from softqec.pauli import pauli_string_to_symplectic
from softqec.stabilizer import (
    commutation_matrix,
    stabilizers_commute,
    syndrome,
)


def build_three_qubit_repetition_stabilizers():
    """Return stabilizers ZZI and IZZ in symplectic form."""

    stabilizers = ["ZZI", "IZZ"]

    x_rows = []
    z_rows = []

    for stabilizer in stabilizers:
        x, z = pauli_string_to_symplectic(stabilizer)
        x_rows.append(x)
        z_rows.append(z)

    return np.array(x_rows), np.array(z_rows)


def test_repetition_stabilizers_commute():
    sx, sz = build_three_qubit_repetition_stabilizers()

    matrix = commutation_matrix(sx, sz)

    np.testing.assert_array_equal(
        matrix,
        np.zeros((2, 2), dtype=np.uint8),
    )

    assert stabilizers_commute(sx, sz)


def test_no_error_has_zero_syndrome():
    sx, sz = build_three_qubit_repetition_stabilizers()

    error_x, error_z = pauli_string_to_symplectic("III")

    result = syndrome(
        error_x,
        error_z,
        sx,
        sz,
    )

    np.testing.assert_array_equal(result, [0, 0])


def test_x_error_on_first_qubit():
    sx, sz = build_three_qubit_repetition_stabilizers()

    error_x, error_z = pauli_string_to_symplectic("XII")

    result = syndrome(
        error_x,
        error_z,
        sx,
        sz,
    )

    np.testing.assert_array_equal(result, [1, 0])


def test_x_error_on_second_qubit():
    sx, sz = build_three_qubit_repetition_stabilizers()

    error_x, error_z = pauli_string_to_symplectic("IXI")

    result = syndrome(
        error_x,
        error_z,
        sx,
        sz,
    )

    np.testing.assert_array_equal(result, [1, 1])


def test_x_error_on_third_qubit():
    sx, sz = build_three_qubit_repetition_stabilizers()

    error_x, error_z = pauli_string_to_symplectic("IIX")

    result = syndrome(
        error_x,
        error_z,
        sx,
        sz,
    )

    np.testing.assert_array_equal(result, [0, 1])


def test_z_error_is_not_detected_by_z_checks():
    sx, sz = build_three_qubit_repetition_stabilizers()

    error_x, error_z = pauli_string_to_symplectic("ZII")

    result = syndrome(
        error_x,
        error_z,
        sx,
        sz,
    )

    np.testing.assert_array_equal(result, [0, 0])
