"""Error-rate metrics and binomial confidence intervals."""

from __future__ import annotations

from dataclasses import dataclass
from statistics import NormalDist

import numpy as np
from numpy.typing import ArrayLike


@dataclass(frozen=True)
class BinomialEstimate:
    errors: int
    trials: int
    rate: float
    ci_low: float
    ci_high: float


def frame_error_count(reference: ArrayLike, estimate: ArrayLike) -> int:
    """Count rows containing at least one mismatch."""
    expected = np.asarray(reference)
    actual = np.asarray(estimate)
    if expected.shape != actual.shape:
        raise ValueError(f"shape mismatch: {expected.shape} != {actual.shape}")
    if expected.ndim == 1:
        return int(np.any(expected != actual))
    if expected.ndim != 2:
        raise ValueError("frame arrays must be one- or two-dimensional")
    return int(np.count_nonzero(np.any(expected != actual, axis=1)))


def wilson_interval(errors: int, trials: int, confidence: float = 0.95) -> tuple[float, float]:
    """Compute a two-sided Wilson score interval for a binomial proportion."""
    if trials <= 0:
        raise ValueError("trials must be positive")
    if not 0 <= errors <= trials:
        raise ValueError("errors must lie between zero and trials")
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must lie in (0,1)")
    z = NormalDist().inv_cdf(0.5 + confidence / 2.0)
    proportion = errors / trials
    denominator = 1.0 + z * z / trials
    center = (proportion + z * z / (2.0 * trials)) / denominator
    margin = z * np.sqrt(
        proportion * (1.0 - proportion) / trials + z * z / (4.0 * trials * trials)
    ) / denominator
    return max(0.0, float(center - margin)), min(1.0, float(center + margin))


def binomial_estimate(errors: int, trials: int, confidence: float = 0.95) -> BinomialEstimate:
    low, high = wilson_interval(errors, trials, confidence)
    return BinomialEstimate(errors, trials, errors / trials, low, high)
