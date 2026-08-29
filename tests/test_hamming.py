import unittest

import numpy as np

from softqec import gf2
from softqec.classical_codes import HammingCode74
from softqec.classical_decoders import hard_syndrome_decode


class TestHammingCode74(unittest.TestCase):
    def setUp(self) -> None:
        self.code = HammingCode74()

    def test_generator_and_parity_check_are_compatible(self) -> None:
        np.testing.assert_array_equal(
            gf2.matmul(self.code.parity_check, self.code.generator.T),
            np.zeros((3, 4), dtype=np.uint8),
        )

    def test_codebook_has_16_unique_words_and_distance_three(self) -> None:
        codewords = self.code.codewords()
        self.assertEqual(len(np.unique(codewords, axis=0)), 16)
        self.assertEqual(self.code.minimum_distance(), 3)

    def test_all_codewords_have_zero_syndrome(self) -> None:
        np.testing.assert_array_equal(
            self.code.syndrome(self.code.codewords()),
            np.zeros((16, 3), dtype=np.uint8),
        )

    def test_all_112_single_bit_cases_are_corrected(self) -> None:
        codewords = self.code.codewords()
        for transmitted in codewords:
            for bit_index in range(self.code.n):
                received = transmitted.copy()
                received[bit_index] ^= 1
                corrected, correction = hard_syndrome_decode(received, self.code)
                np.testing.assert_array_equal(corrected, transmitted)
                self.assertEqual(int(np.sum(correction)), 1)
                self.assertEqual(int(correction[bit_index]), 1)


if __name__ == "__main__":
    unittest.main()
