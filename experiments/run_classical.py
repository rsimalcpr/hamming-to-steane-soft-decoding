"""Run the Week 1 classical Hamming hard-versus-soft experiment."""

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

from softqec.classical_channels import hard_demodulate, transmit_bpsk_awgn
from softqec.classical_codes import HammingCode74
from softqec.classical_decoders import exact_soft_ml_decode, hard_syndrome_decode
from softqec.metrics import binomial_estimate, frame_error_count

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def load_config(path: Path) -> dict[str, Any]:
    with path.open("rb") as stream:
        return tomllib.load(stream)


def simulate_point(
    ebn0_db: float,
    trials: int,
    seed: int,
    code: HammingCode74,
    *,
    chunk_size: int,
) -> dict[str, float | int]:
    """Simulate all three curves at one Eb/N0 value."""
    rng = np.random.default_rng(seed)
    messages = rng.integers(0, 2, size=(trials, code.k), dtype=np.uint8)

    uncoded_observations = transmit_bpsk_awgn(messages, ebn0_db, rate=1.0, rng=rng)
    uncoded_estimates = hard_demodulate(uncoded_observations)

    transmitted_codewords = code.encode(messages)
    coded_observations = transmit_bpsk_awgn(
        transmitted_codewords,
        ebn0_db,
        rate=code.rate,
        rng=rng,
    )

    hard_words, _ = hard_syndrome_decode(hard_demodulate(coded_observations), code)
    hard_messages = code.extract_messages(hard_words)

    soft_words = exact_soft_ml_decode(coded_observations, code, chunk_size=chunk_size)
    soft_messages = code.extract_messages(soft_words)

    estimates = {
        "uncoded": binomial_estimate(frame_error_count(messages, uncoded_estimates), trials),
        "hard": binomial_estimate(frame_error_count(messages, hard_messages), trials),
        "soft": binomial_estimate(frame_error_count(messages, soft_messages), trials),
    }

    row: dict[str, float | int] = {"ebn0_db": float(ebn0_db), "trials": int(trials)}
    for label, estimate in estimates.items():
        row.update({f"{label}_{key}": value for key, value in asdict(estimate).items() if key != "trials"})
    return row


def run_experiment(config: dict[str, Any], *, trials_override: int | None = None) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    experiment = config["experiment"]
    ebn0_values = [float(value) for value in experiment["ebn0_db"]]
    trials = int(experiment["trials_per_point"] if trials_override is None else trials_override)
    seed = int(experiment["seed"])
    chunk_size = int(experiment.get("chunk_size", 50_000))
    if trials <= 0:
        raise ValueError("trials_per_point must be positive")

    code = HammingCode74()
    seed_sequence = np.random.SeedSequence(seed)
    child_sequences = seed_sequence.spawn(len(ebn0_values))
    child_seeds = [int(child.generate_state(1, dtype=np.uint64)[0]) for child in child_sequences]

    started = time.perf_counter()
    rows = [
        simulate_point(value, trials, point_seed, code, chunk_size=chunk_size)
        for value, point_seed in zip(ebn0_values, child_seeds, strict=True)
    ]
    elapsed = time.perf_counter() - started

    metadata = {
        "experiment": "classical_hamming_fer",
        "top_level_seed": seed,
        "point_seeds": child_seeds,
        "ebn0_db": ebn0_values,
        "trials_per_point": trials,
        "chunk_size": chunk_size,
        "elapsed_seconds": elapsed,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "conventions": {
            "bpsk": "0 -> +1, 1 -> -1",
            "noise_variance": "1 / (2 * R * 10**(EbN0_dB/10))",
            "frame": "four information bits",
        },
    }
    return rows, metadata


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--trials", type=int, help="Override trials_per_point for a smoke run")
    parser.add_argument("--no-figure", action="store_true", help="Save data without rebuilding the plot")
    args = parser.parse_args()

    config_path = args.config if args.config.is_absolute() else REPOSITORY_ROOT / args.config
    config = load_config(config_path)
    rows, metadata = run_experiment(config, trials_override=args.trials)

    output = config["output"]
    csv_path = REPOSITORY_ROOT / output["results_csv"]
    metadata_path = REPOSITORY_ROOT / output["metadata_json"]
    figure_path = REPOSITORY_ROOT / output["figure_path"]
    write_csv(rows, csv_path)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    if not args.no_figure:
        from experiments.make_figures import plot_classical_fer

        plot_classical_fer(csv_path, figure_path)

    print(f"Saved {csv_path.relative_to(REPOSITORY_ROOT)}")
    print(f"Saved {metadata_path.relative_to(REPOSITORY_ROOT)}")
    if not args.no_figure:
        print(f"Saved {figure_path.relative_to(REPOSITORY_ROOT)}")
    print(f"Runtime: {metadata['elapsed_seconds']:.3f} s")


if __name__ == "__main__":
    main()
