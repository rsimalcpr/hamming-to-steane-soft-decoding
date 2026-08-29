import unittest

import numpy as np

from softqec.metrics import binomial_estimate, frame_error_count, wilson_interval


class TestMetrics(unittest.TestCase):
    def test_frame_error_count(self) -> None:
        reference = np.zeros((3, 4), dtype=np.uint8)
        estimate = reference.copy()
        estimate[0, 0] = 1
        estimate[2, 1:3] = 1
        self.assertEqual(frame_error_count(reference, estimate), 2)

    def test_wilson_interval_contains_observed_rate(self) -> None:
        low, high = wilson_interval(10, 100)
        self.assertLess(low, 0.1)
        self.assertGreater(high, 0.1)

    def test_zero_error_interval_is_not_zero_width(self) -> None:
        estimate = binomial_estimate(0, 100)
        self.assertEqual(estimate.rate, 0.0)
        self.assertEqual(estimate.ci_low, 0.0)
        self.assertGreater(estimate.ci_high, 0.0)


if __name__ == "__main__":
    unittest.main()
