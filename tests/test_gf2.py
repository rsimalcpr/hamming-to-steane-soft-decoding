import unittest

import numpy as np

from softqec import gf2


class TestGF2(unittest.TestCase):
    def test_addition_is_xor(self) -> None:
        left = np.array([0, 1, 1], dtype=np.uint8)
        right = np.array([1, 1, 0], dtype=np.uint8)
        np.testing.assert_array_equal(gf2.add(left, right), [1, 0, 1])

    def test_rref_and_rank(self) -> None:
        matrix = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]], dtype=np.uint8)
        reduced, pivots = gf2.rref(matrix)
        self.assertEqual(pivots, (0, 1))
        self.assertEqual(gf2.rank(matrix), 2)
        np.testing.assert_array_equal(reduced, [[1, 0, 1], [0, 1, 1], [0, 0, 0]])

    def test_null_space(self) -> None:
        matrix = np.array([[1, 1, 0], [0, 1, 1]], dtype=np.uint8)
        basis = gf2.null_space(matrix)
        self.assertEqual(basis.shape, (1, 3))
        np.testing.assert_array_equal(gf2.matmul(matrix, basis.T), np.zeros((2, 1), dtype=np.uint8))

    def test_rejects_nonbinary_input(self) -> None:
        with self.assertRaises(ValueError):
            gf2.rank([[0, 2]])

    def test_rejects_incompatible_multiplication(self) -> None:
        with self.assertRaises(ValueError):
            gf2.matmul(np.zeros((2, 3), dtype=np.uint8), np.zeros((4, 2), dtype=np.uint8))


if __name__ == "__main__":
    unittest.main()
