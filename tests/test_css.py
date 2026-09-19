import numpy as np
import pytest

from softqec.css import (
    css_checks_commute,
    css_commutation_matrix,
    css_syndrome,
    number_logical_qubits,
    x_error_syndrome,
    z_error_syndrome,
)


# A small CSS example.
#
# X stabilizers:
#
#   XXII
#   IIXX
#
# Z stabilizer:
#
#   ZZZZ
#
# Every X stabilizer overlaps ZZZZ on exactly two qubits,
# so all generators commute.
H_X = np.array(
    [
        [1, 1, 0, 0],
        [0, 0, 1, 1],
    ],
    dtype=np.uint8,
)

H_Z = np.array(
    [
        [1, 1, 1, 1],
    ],
    dtype=np.uint8,
)


def test_css_commutation_condition():
    result = css_commutation_matrix(H_X, H_Z)

    np.testing.assert_array_equal(
        result,
        np.zeros((2, 1), dtype=np.uint8),
    )

    assert css_checks_commute(H_X, H_Z)


def test_invalid_css_checks_are_detected():
    bad_h_z = np.array(
        [
            [1, 0, 0, 0],
        ],
        dtype=np.uint8,
    )

    result = css_commutation_matrix(H_X, bad_h_z)

    assert np.any(result == 1)
    assert not css_checks_commute(H_X, bad_h_z)



def test_number_of_logical_qubits():
    # n = 4
    # rank(H_X) = 2
    # rank(H_Z) = 1
    #
    # k = 4 - 2 - 1 = 1

    assert number_logical_qubits(H_X, H_Z) == 1


def test_x_error_is_detected_by_z_checks():
    error_x = np.array([1, 0, 0, 0], dtype=np.uint8)

    syndrome = x_error_syndrome(
        error_x,
        H_Z,
    )

    np.testing.assert_array_equal(syndrome, [1])


def test_z_error_is_detected_by_x_checks():
    error_z = np.array([1, 0, 0, 0], dtype=np.uint8)

    syndrome = z_error_syndrome(
        error_z,
        H_X,
    )

    np.testing.assert_array_equal(syndrome, [1, 0])


def test_y_error_has_both_syndrome_components():
    # Y = XZ, so both binary components are 1
    # on the first physical qubit.

    error_x = np.array([1, 0, 0, 0], dtype=np.uint8)
    error_z = np.array([1, 0, 0, 0], dtype=np.uint8)

    syndrome_x_checks, syndrome_z_checks = css_syndrome(
        error_x,
        error_z,
        H_X,
        H_Z,
    )

    np.testing.assert_array_equal(
        syndrome_x_checks,
        [1, 0],
    )

    np.testing.assert_array_equal(
        syndrome_z_checks,
        [1],
    )


def test_noncommuting_css_code_has_no_valid_k():
    bad_h_z = np.array(
        [
            [1, 0, 0, 0],
        ],
        dtype=np.uint8,
    )

    with pytest.raises(ValueError):
        number_logical_qubits(H_X, bad_h_z)
