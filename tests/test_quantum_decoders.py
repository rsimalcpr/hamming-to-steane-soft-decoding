import numpy as np
import pytest

from softqec.quantum_decoders import (
    corrections_from_indices,
    hard_nearest_syndrome_decode,
    ideal_analog_signatures,
    restricted_failure_count,
    restricted_syndrome_table,
    single_pauli_hypotheses,
    soft_ml_syndrome_decode,
)
from softqec.steane import (
    N_PHYSICAL,
    STEANE_H,
    correction_succeeds,
)


def test_restricted_hypothesis_set_contains_identity_and_21_paulis():
    labels, errors_x, errors_z = (
        single_pauli_hypotheses()
    )

    assert len(labels) == 22
    assert labels[0] == "I"

    assert labels[1:8] == tuple(
        f"X{i}" for i in range(1, 8)
    )

    assert labels[8:15] == tuple(
        f"Y{i}" for i in range(1, 8)
    )

    assert labels[15:22] == tuple(
        f"Z{i}" for i in range(1, 8)
    )

    assert errors_x.shape == (
        22,
        N_PHYSICAL,
    )

    assert errors_z.shape == (
        22,
        N_PHYSICAL,
    )


def test_all_22_restricted_syndromes_are_unique():
    table = restricted_syndrome_table()

    assert table.shape == (22, 6)

    unique_syndromes = {
        tuple(row)
        for row in table
    }

    assert len(unique_syndromes) == 22


def test_single_pauli_syndrome_structure():
    table = restricted_syndrome_table()

    zero = np.zeros(
        3,
        dtype=np.uint8,
    )

    for qubit in range(N_PHYSICAL):
        column = STEANE_H[:, qubit]

        # X_i -> [000 | H_i]
        np.testing.assert_array_equal(
            table[1 + qubit],
            np.r_[zero, column],
        )

        # Y_i -> [H_i | H_i]
        np.testing.assert_array_equal(
            table[8 + qubit],
            np.r_[column, column],
        )

        # Z_i -> [H_i | 000]
        np.testing.assert_array_equal(
            table[15 + qubit],
            np.r_[column, zero],
        )


@pytest.mark.parametrize(
    "index",
    range(22),
)
def test_both_decoders_recover_every_ideal_signature(
    index,
):
    observation = ideal_analog_signatures()[
        index
    ]

    assert (
        hard_nearest_syndrome_decode(
            observation
        )
        == index
    )

    assert (
        soft_ml_syndrome_decode(
            observation
        )
        == index
    )


def test_soft_decoder_can_use_reliability_discarded_by_thresholding():
    rng = np.random.default_rng(19)

    signatures = ideal_analog_signatures()

    for _ in range(10_000):
        actual = int(
            rng.integers(0, 22)
        )

        observation = (
            signatures[actual]
            + rng.normal(
                0.0,
                0.9,
                size=6,
            )
        )

        hard = hard_nearest_syndrome_decode(
            observation
        )

        soft = soft_ml_syndrome_decode(
            observation
        )

        if (
            hard != soft
            and soft == actual
        ):
            break

    else:
        pytest.fail(
            "no deterministic hard/soft "
            "reliability witness found"
        )


def test_restricted_success_matches_stabilizer_success_criterion():
    labels, errors_x, errors_z = (
        single_pauli_hypotheses()
    )

    for actual in range(len(labels)):
        for decoded in range(len(labels)):
            correction_x, correction_z = (
                corrections_from_indices(
                    decoded
                )
            )

            succeeds = correction_succeeds(
                errors_x[actual],
                errors_z[actual],
                correction_x,
                correction_z,
            )

            assert succeeds == (
                actual == decoded
            )


def test_restricted_failure_count_validates_shapes():
    failures = restricted_failure_count(
        [0, 1, 2],
        [0, 1, 3],
    )

    assert failures == 1

    with pytest.raises(ValueError):
        restricted_failure_count(
            [0, 1],
            [0],
        )
