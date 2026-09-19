import numpy as np
import pytest

from softqec.pauli import (
    commutes,
    pauli_string_to_symplectic,
    pauli_weight,
    symplectic_product,
)


def test_single_qubit_pauli_mapping():
    x, z = pauli_string_to_symplectic("I")
    np.testing.assert_array_equal(x, [0])
    np.testing.assert_array_equal(z, [0])

    x, z = pauli_string_to_symplectic("X")
    np.testing.assert_array_equal(x, [1])
    np.testing.assert_array_equal(z, [0])

    x, z = pauli_string_to_symplectic("Z")
    np.testing.assert_array_equal(x, [0])
    np.testing.assert_array_equal(z, [1])

    x, z = pauli_string_to_symplectic("Y")
    np.testing.assert_array_equal(x, [1])
    np.testing.assert_array_equal(z, [1])


def test_x_and_z_anticommute():
    x1, z1 = pauli_string_to_symplectic("X")
    x2, z2 = pauli_string_to_symplectic("Z")

    assert symplectic_product(x1, z1, x2, z2) == 1
    assert not commutes(x1, z1, x2, z2)


def test_x_and_x_commute():
    x1, z1 = pauli_string_to_symplectic("X")
    x2, z2 = pauli_string_to_symplectic("X")

    assert commutes(x1, z1, x2, z2)


def test_two_anticommutations_cancel():
    x1, z1 = pauli_string_to_symplectic("XX")
    x2, z2 = pauli_string_to_symplectic("ZZ")

    assert commutes(x1, z1, x2, z2)


def test_disjoint_paulis_commute():
    x1, z1 = pauli_string_to_symplectic("XI")
    x2, z2 = pauli_string_to_symplectic("IZ")

    assert commutes(x1, z1, x2, z2)


def test_pauli_weight():
    x, z = pauli_string_to_symplectic("IXYZI")

    assert pauli_weight(x, z) == 3


def test_invalid_pauli_character():
    with pytest.raises(ValueError):
        pauli_string_to_symplectic("XA")
