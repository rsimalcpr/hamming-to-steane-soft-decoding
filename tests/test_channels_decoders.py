import unittest

import numpy as np

from softqec.classical_channels import (
    binary_symmetric_channel,
    bpsk_modulate,
    hard_demodulate,
    noise_sigma_from_ebn0_db,
)
from softqec.classical_codes import HammingCode74
from softqec.classical_decoders import exact_soft_ml_decode, hard_syndrome_decode


class TestChannelsAndDecoders(unittest.TestCase):
    def setUp(self) -> None:
        self.code = HammingCode74()

    def test_bsc_limits(self) -> None:
        bits = np.array([[0, 1, 0], [1, 1, 0]], dtype=np.uint8)
        np.testing.assert_array_equal(binary_symmetric_channel(bits, 0.0, np.random.default_rng(1)), bits)
        np.testing.assert_array_equal(binary_symmetric_channel(bits, 1.0, np.random.default_rng(1)), 1 - bits)

    def test_bpsk_and_threshold_are_inverse_without_noise(self) -> None:
        bits = np.array([[0, 1, 1, 0]], dtype=np.uint8)
        np.testing.assert_array_equal(hard_demodulate(bpsk_modulate(bits)), bits)

    def test_ebn0_sigma_convention(self) -> None:
        self.assertAlmostEqual(noise_sigma_from_ebn0_db(0.0, 1.0), np.sqrt(0.5))
        self.assertAlmostEqual(noise_sigma_from_ebn0_db(0.0, 4 / 7), np.sqrt(7 / 8))

    def test_noiseless_hard_and_soft_decode_every_codeword(self) -> None:
        codewords = self.code.codewords()
        hard_words, _ = hard_syndrome_decode(hard_demodulate(bpsk_modulate(codewords)), self.code)
        soft_words = exact_soft_ml_decode(bpsk_modulate(codewords), self.code)
        np.testing.assert_array_equal(hard_words, codewords)
        np.testing.assert_array_equal(soft_words, codewords)

    def test_soft_decoder_uses_reliability(self) -> None:
        transmitted = self.code.codewords()[0]
        observation = bpsk_modulate(transmitted)
        observation[[0, 1]] = [-0.05, -0.05]
        decoded = exact_soft_ml_decode(observation, self.code)
        np.testing.assert_array_equal(decoded, transmitted)


if __name__ == "__main__":
    unittest.main()
