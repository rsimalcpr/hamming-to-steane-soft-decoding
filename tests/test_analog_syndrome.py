import numpy as np
import pytest

from softqec.analog_syndrome import (
    hard_threshold_syndrome,
    sample_analog_syndrome,
    syndrome_bits_to_means,
)


def test_bit_to_eigenvalue_mapping_and_threshold_inverse():
    bits = np.array(
        [0, 1, 0, 1, 1, 0],
        dtype=np.uint8,
    )

    means = syndrome_bits_to_means(bits)

    np.testing.assert_array_equal(
        means,
        [1, -1, 1, -1, -1, 1],
    )

    np.testing.assert_array_equal(
        hard_threshold_syndrome(means),
        bits,
    )


def test_zero_noise_returns_ideal_means():
    bits = np.array(
        [[0, 0, 1, 1, 0, 1]],
        dtype=np.uint8,
    )

    observed = sample_analog_syndrome(
        bits,
        sigma_m=0.0,
        rng=np.random.default_rng(7),
    )

    np.testing.assert_array_equal(
        observed,
        syndrome_bits_to_means(bits),
    )


def test_sampling_is_reproducible_with_fixed_seed():
    bits = np.zeros(
        (4, 6),
        dtype=np.uint8,
    )

    first = sample_analog_syndrome(
        bits,
        sigma_m=0.5,
        rng=np.random.default_rng(12),
    )

    second = sample_analog_syndrome(
        bits,
        sigma_m=0.5,
        rng=np.random.default_rng(12),
    )

    np.testing.assert_array_equal(
        first,
        second,
    )


def test_invalid_sigma_is_rejected():
    rng = np.random.default_rng(1)

    bits = np.zeros(
        6,
        dtype=np.uint8,
    )

    with pytest.raises(ValueError):
        sample_analog_syndrome(
            bits,
            sigma_m=-0.1,
            rng=rng,
        )


def test_nonfinite_observation_is_rejected():
    with pytest.raises(ValueError):
        hard_threshold_syndrome(
            [1, 1, 1, 1, 1, np.nan]
        )
