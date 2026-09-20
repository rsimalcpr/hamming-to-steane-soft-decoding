"""Build publication-ready figures from saved experiment tables."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


def read_rows(
    path: Path,
) -> list[dict[str, float]]:
    """Read a numeric CSV table."""
    with path.open(
        newline="",
        encoding="utf-8",
    ) as stream:
        return [
            {
                key: float(value)
                for key, value in row.items()
            }
            for row in csv.DictReader(stream)
        ]


def plot_classical_fer(
    input_path: Path,
    output_path: Path,
) -> None:
    """Plot the Week 1 classical frame-error rates."""
    rows = read_rows(input_path)

    ebn0 = np.array(
        [
            row["ebn0_db"]
            for row in rows
        ]
    )

    trials = np.array(
        [
            row["trials"]
            for row in rows
        ]
    )

    styles = {
        "uncoded": {
            "label": (
                "Uncoded hard decision "
                "(4-bit frame)"
            ),
            "marker": "o",
            "color": "#777777",
        },
        "hard": {
            "label": (
                "Hamming hard syndrome"
            ),
            "marker": "s",
            "color": "#d95f02",
        },
        "soft": {
            "label": (
                "Hamming exact soft ML"
            ),
            "marker": "^",
            "color": "#007c83",
        },
    }

    figure, axis = plt.subplots(
        figsize=(7.2, 4.8),
        constrained_layout=True,
    )

    for name, style in styles.items():
        rate = np.array(
            [
                row[f"{name}_rate"]
                for row in rows
            ]
        )

        low = np.array(
            [
                row[f"{name}_ci_low"]
                for row in rows
            ]
        )

        high = np.array(
            [
                row[f"{name}_ci_high"]
                for row in rows
            ]
        )

        display_floor = 0.5 / trials

        axis.semilogy(
            ebn0,
            np.maximum(
                rate,
                display_floor,
            ),
            linewidth=2.0,
            markersize=5.5,
            **style,
        )

        axis.fill_between(
            ebn0,
            np.maximum(
                low,
                display_floor,
            ),
            np.maximum(
                high,
                display_floor,
            ),
            color=style["color"],
            alpha=0.14,
            linewidth=0,
        )

    axis.set_xlabel(
        r"$E_b/N_0$ (dB)"
    )

    axis.set_ylabel(
        "Frame-error rate"
    )

    axis.set_title(
        "Hamming [7,4,3]: "
        "hard versus soft decoding"
    )

    axis.grid(
        True,
        which="both",
        alpha=0.28,
    )

    axis.legend(
        frameon=False
    )

    axis.set_ylim(
        bottom=max(
            1e-6,
            0.4 / float(
                np.max(trials)
            ),
        ),
        top=1.0,
    )

    axis.text(
        0.012,
        0.018,
        (
            "Zero-count points are displayed "
            "at 0.5/N on the log axis; "
            "bands are 95% Wilson intervals."
        ),
        transform=axis.transAxes,
        fontsize=7.5,
        color="#555555",
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure.savefig(
        output_path,
        dpi=220,
    )

    plt.close(figure)


def plot_quantum_single_pauli_failure(
    input_path: Path,
    output_path: Path,
) -> None:
    """Plot Week 3 hard and soft decoding failure rates."""
    rows = read_rows(input_path)

    sigma_m = np.array(
        [
            row["sigma_m"]
            for row in rows
        ]
    )

    styles = {
        "hard": {
            "label": (
                "Hard threshold "
                "+ nearest syndrome"
            ),
            "marker": "o",
            "color": "#d95f02",
        },
        "soft": {
            "label": (
                "Soft Gaussian ML"
            ),
            "marker": "s",
            "color": "#007c83",
        },
    }

    figure, axis = plt.subplots(
        figsize=(7.2, 4.8),
        constrained_layout=True,
    )

    for name, style in styles.items():
        rate = np.array(
            [
                row[f"{name}_rate"]
                for row in rows
            ]
        )

        low = np.array(
            [
                row[f"{name}_ci_low"]
                for row in rows
            ]
        )

        high = np.array(
            [
                row[f"{name}_ci_high"]
                for row in rows
            ]
        )

        axis.plot(
            sigma_m,
            rate,
            linewidth=2.0,
            markersize=5.5,
            **style,
        )

        axis.fill_between(
            sigma_m,
            low,
            high,
            color=style["color"],
            alpha=0.14,
            linewidth=0,
        )

    axis.set_xlabel(
        r"Syndrome measurement noise $\sigma_m$"
    )

    axis.set_ylabel(
        "Decoding failure rate"
    )

    axis.set_title(
        "Steane [[7,1,3]] "
        "single-Pauli decoding"
    )

    axis.grid(
        True,
        alpha=0.28,
    )

    axis.legend(
        frameon=False
    )

    axis.set_ylim(
        bottom=0.0,
        top=1.0,
    )

    axis.text(
        0.012,
        0.018,
        (
            "Restricted equally likely "
            "{I, Xi, Yi, Zi} model; "
            "bands are 95% Wilson intervals."
        ),
        transform=axis.transAxes,
        fontsize=7.5,
        color="#555555",
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure.savefig(
        output_path,
        dpi=220,
    )

    plt.close(figure)


def main() -> None:
    """Rebuild one figure from an existing CSV table."""
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--output",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--kind",
        choices=(
            "classical",
            "quantum-single-pauli",
        ),
        default="classical",
    )

    args = parser.parse_args()

    if args.kind == "classical":
        plot_classical_fer(
            args.input,
            args.output,
        )
    else:
        plot_quantum_single_pauli_failure(
            args.input,
            args.output,
        )

    print(
        f"Saved {args.output}"
    )


if __name__ == "__main__":
    main()
