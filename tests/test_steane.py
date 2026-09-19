import numpy as np

from softqec.css import (
    css_checks_commute,
    css_syndrome,
    number_logical_qubits,
)
from softqec.gf2 import rank
from softqec.pauli import (
    commutes,
    pauli_string_to_symplectic,
)
from softqec.steane import (
    DISTANCE,
    LOGICAL_X,
    LOGICAL_Z,
    N_LOGICAL,
    N_PHYSICAL,
    STEANE_H,
    STEANE_H_X,
    STEANE_H_Z,
    correction_succeeds,
    decode_ideal_css,
    in_row_space,
    single_qubit_lookup,
    steane_parameters,
    steane_stabilizers,
    validate_steane_construction,
)


def test_steane_css_checks_commute():
    assert css_checks_commute(
        STEANE_H_X,
        STEANE_H_Z,
    )


def test_h_h_transpose_is_zero():
    result = (STEANE_H @ STEANE_H.T) % 2

    np.testing.assert_array_equal(
        result,
        np.zeros((3, 3), dtype=np.uint8),
    )


def test_steane_check_ranks():
    assert rank(STEANE_H_X) == 3
    assert rank(STEANE_H_Z) == 3


def test_steane_encodes_one_logical_qubit():
    assert number_logical_qubits(
        STEANE_H_X,
        STEANE_H_Z,
    ) == 1

    assert steane_parameters() == (7, 1, 3)


def test_six_stabilizer_generators():
    sx, sz = steane_stabilizers()

    assert sx.shape == (6, 7)
    assert sz.shape == (6, 7)


def test_logical_operators_are_not_stabilizers():
    assert not in_row_space(
        LOGICAL_X,
        STEANE_H_X,
    )

    assert not in_row_space(
        LOGICAL_Z,
        STEANE_H_Z,
    )


def test_logical_x_commutes_with_stabilizers():
    sx, sz = steane_stabilizers()

    logical_x_symplectic = LOGICAL_X
    logical_z_zero = np.zeros(
        N_PHYSICAL,
        dtype=np.uint8,
    )

    for row in range(6):
        assert commutes(
            logical_x_symplectic,
            logical_z_zero,
            sx[row],
            sz[row],
        )


def test_logical_z_commutes_with_stabilizers():
    sx, sz = steane_stabilizers()

    logical_x_zero = np.zeros(
        N_PHYSICAL,
        dtype=np.uint8,
    )
    logical_z_symplectic = LOGICAL_Z

    for row in range(6):
        assert commutes(
            logical_x_zero,
            logical_z_symplectic,
            sx[row],
            sz[row],
        )


def test_logical_x_and_z_anticommute():
    zero = np.zeros(
        N_PHYSICAL,
        dtype=np.uint8,
    )

    assert not commutes(
        LOGICAL_X,
        zero,
        zero,
        LOGICAL_Z,
    )


def test_all_seven_nonzero_syndromes_are_present():
    lookup = single_qubit_lookup()

    assert len(lookup) == 7

    assert (0, 0, 0) not in lookup


def test_every_single_qubit_x_error_is_corrected():
    for qubit in range(N_PHYSICAL):
        error_x = np.zeros(
            N_PHYSICAL,
            dtype=np.uint8,
        )
        error_z = np.zeros(
            N_PHYSICAL,
            dtype=np.uint8,
        )

        error_x[qubit] = 1

        correction_x, correction_z = decode_ideal_css(
            error_x,
            error_z,
        )

        assert correction_succeeds(
            error_x,
            error_z,
            correction_x,
            correction_z,
        )


def test_every_single_qubit_z_error_is_corrected():
    for qubit in range(N_PHYSICAL):
        error_x = np.zeros(
            N_PHYSICAL,
            dtype=np.uint8,
        )
        error_z = np.zeros(
            N_PHYSICAL,
            dtype=np.uint8,
        )

        error_z[qubit] = 1

        correction_x, correction_z = decode_ideal_css(
            error_x,
            error_z,
        )

        assert correction_succeeds(
            error_x,
            error_z,
            correction_x,
            correction_z,
        )


def test_every_single_qubit_y_error_is_corrected():
    for qubit in range(N_PHYSICAL):
        error_x = np.zeros(
            N_PHYSICAL,
            dtype=np.uint8,
        )
        error_z = np.zeros(
            N_PHYSICAL,
            dtype=np.uint8,
        )

        # Y = XZ, ignoring the global phase.
        error_x[qubit] = 1
        error_z[qubit] = 1

        correction_x, correction_z = decode_ideal_css(
            error_x,
            error_z,
        )

        assert correction_succeeds(
            error_x,
            error_z,
            correction_x,
            correction_z,
        )


def test_validate_steane_construction():
    validate_steane_construction()
