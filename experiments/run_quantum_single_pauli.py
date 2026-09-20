"""Run the 3rd Section restricted Steane hard-versus-soft ML experiment."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import time
import tomllib
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from softqec.analog_syndrome import (
    sample_analog_syndrome,
)
from softqec.metrics import (
    binomial_estimate,
)
from softqec.quantum_decoders import (
    hard_nearest_syndrome_decode,
    restricted_failure_count,
    restricted_syndrome_table,
    single_pauli_hypotheses,
    soft_ml_syndrome_decode,
)

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def load_config(
    path: Path,
) -> dict[str, Any]:
    """Load a TOML experiment configuration."""
    with path.open("rb") as stream:
        return tomllib.load(stream)


def simulate_point(
    sigma_m: float,
    trials: int,
    seed: int,
) -> dict[str, float | int]:
    """Simulate one measurement-noise value with paired decoder inputs."""
    if trials <= 0:
        raise ValueError(
            "trials must be positive"
        )

    if (
        not np.isfinite(sigma_m)
        or sigma_m < 0.0
    ):
        raise ValueError(
            "sigma_m must be finite and non-negative"
        )

    labels, _, _ = single_pauli_hypotheses()

    syndrome_table = restricted_syndrome_table()

    rng = np.random.default_rng(seed)

    actual = rng.integers(
        0,
        len(labels),
        size=trials,
    )

    observations = sample_analog_syndrome(
        syndrome_table[actual],
        sigma_m,
        rng,
    )

    hard_started = time.perf_counter()

    hard = hard_nearest_syndrome_decode(
        observations
    )

    hard_runtime_ms = (
        1_000.0
        * (
            time.perf_counter()
            - hard_started
        )
    )

    soft_started = time.perf_counter()

    soft = soft_ml_syndrome_decode(
        observations
    )

    soft_runtime_ms = (
        1_000.0
        * (
            time.perf_counter()
            - soft_started
        )
    )

    estimates = {
        "hard": binomial_estimate(
            restricted_failure_count(
                actual,
                hard,
            ),
            trials,
        ),
        "soft": binomial_estimate(
            restricted_failure_count(
                actual,
                soft,
            ),
            trials,
        ),
    }

    row: dict[str, float | int] = {
        "sigma_m": float(sigma_m),
        "trials": int(trials),
        "hard_runtime_ms": hard_runtime_ms,
        "soft_runtime_ms": soft_runtime_ms,
    }

    for name, estimate in estimates.items():
        row.update(
            {
                f"{name}_{key}": value
                for key, value
                in asdict(estimate).items()
                if key != "trials"
            }
        )

    return row


def run_experiment(
    config: dict[str, Any],
    *,
    trials_override: int | None = None,
) -> tuple[
    list[dict[str, Any]],
    dict[str, Any],
]:
    """Run every configured syndrome-noise point."""
    experiment = config["experiment"]

    sigmas = [
        float(value)
        for value in experiment["sigma_m"]
    ]

    trials = int(
        experiment["trials_per_point"]
        if trials_override is None
        else trials_override
    )

    seed = int(
        experiment["seed"]
    )

    if trials <= 0:
        raise ValueError(
            "trials_per_point must be positive"
        )

    seed_sequence = np.random.SeedSequence(
        seed
    )

    children = seed_sequence.spawn(
        len(sigmas)
    )

    point_seeds = [
        int(
            child.generate_state(
                1,
                dtype=np.uint64,
            )[0]
        )
        for child in children
    ]

    started = time.perf_counter()

    rows = [
        simulate_point(
            sigma,
            trials,
            point_seed,
        )
        for sigma, point_seed in zip(
            sigmas,
            point_seeds,
            strict=True,
        )
    ]

    elapsed = (
        time.perf_counter()
        - started
    )

    metadata = {
        "experiment": (
            "quantum_single_pauli_noisy_syndrome"
        ),
        "scope": (
            "equally likely I, Xi, Yi, Zi hypotheses"
        ),
        "top_level_seed": seed,
        "point_seeds": point_seeds,
        "sigma_m": sigmas,
        "trials_per_point": trials,
        "elapsed_seconds": elapsed,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "conventions": {
            "syndrome_order": (
                "X-check responses followed by "
                "Z-check responses"
            ),
            "analog_mapping": (
                "0 -> +1, 1 -> -1"
            ),
            "measurement_model": (
                "independent N(0, sigma_m^2) noise"
            ),
            "hard_decoder": (
                "threshold then nearest valid "
                "syndrome in Hamming distance"
            ),
            "hard_tie_break": (
                "first hypothesis index"
            ),
            "soft_decoder": (
                "minimum squared Euclidean "
                "distance (Gaussian ML)"
            ),
            "success": (
                "residual Pauli belongs to the "
                "Steane stabilizer group"
            ),
        },
        "excluded": [
            "multi-qubit errors",
            "nonuniform priors and MAP decoding",
            (
                "stabilizer-coset probability "
                "aggregation"
            ),
            "repeated syndrome rounds",
            "circuit-level noise",
        ],
    }

    return rows, metadata


def write_csv(
    rows: list[dict[str, Any]],
    path: Path,
) -> None:
    """Write experiment rows to CSV."""
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=list(rows[0]),
        )

        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    """Parse arguments, run the experiment, and save its outputs."""
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--config",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--trials",
        type=int,
        help=(
            "Override trials_per_point "
            "for a smoke run"
        ),
    )

    parser.add_argument(
        "--no-figure",
        action="store_true",
        help=(
            "Save data without rebuilding "
            "the plot"
        ),
    )

    args = parser.parse_args()

    config_path = (
        args.config
        if args.config.is_absolute()
        else REPOSITORY_ROOT / args.config
    )

    config = load_config(
        config_path
    )

    rows, metadata = run_experiment(
        config,
        trials_override=args.trials,
    )

    output = config["output"]

    csv_path = (
        REPOSITORY_ROOT
        / output["results_csv"]
    )

    metadata_path = (
        REPOSITORY_ROOT
        / output["metadata_json"]
    )

    figure_path = (
        REPOSITORY_ROOT
        / output["figure_path"]
    )

    write_csv(
        rows,
        csv_path,
    )

    metadata_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    metadata_path.write_text(
        json.dumps(
            metadata,
            indent=2,
        ),
        encoding="utf-8",
    )

    if not args.no_figure:
        from experiments.make_figures import (
            plot_quantum_single_pauli_failure,
        )

        plot_quantum_single_pauli_failure(
            csv_path,
            figure_path,
        )

    print(
        "Saved "
        f"{csv_path.relative_to(REPOSITORY_ROOT)}"
    )

    print(
        "Saved "
        f"{metadata_path.relative_to(REPOSITORY_ROOT)}"
    )

    if not args.no_figure:
        print(
            "Saved "
            f"{figure_path.relative_to(REPOSITORY_ROOT)}"
        )

    print(
        f"Runtime: "
        f"{metadata['elapsed_seconds']:.3f} s"
    )


if __name__ == "__main__":
    main()
