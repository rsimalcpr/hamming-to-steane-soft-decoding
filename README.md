# From Hamming to Steane: Hard and Soft Decoding

This repository is a transparent classical-to-quantum error-correction study. The current Week 1 milestone implements the classical Hamming [7,4,3] foundation and compares hard syndrome decoding with exact soft maximum-likelihood (ML) decoding over a BPSK/AWGN channel.

## Research question

How much frame-error-rate improvement is obtained when a Hamming decoder retains the continuous BPSK/AWGN observations instead of thresholding them immediately?

## Week 1 pipeline

1. Encode four information bits with a systematic Hamming [7,4,3] generator matrix.
2. Map code bits to BPSK symbols using `0 -> +1` and `1 -> -1`.
3. Add Gaussian noise using a documented information-bit `Eb/N0` convention.
4. Decode the same coded observations in two ways:
   - hard path: threshold, compute syndrome, and apply a single-bit lookup correction;
   - soft path: compare the observation with all 16 valid BPSK codewords and choose exact ML.
5. Estimate frame-error rates and 95% Wilson confidence intervals.

The figure also includes an uncoded four-bit hard-decision baseline at the same `Eb/N0` definition.

## Current reproducible result

With the committed seed and 30,000 frames per point, at `Eb/N0 = 4 dB` the observed Hamming FER is `0.0355` for hard syndrome decoding and `0.01097` for exact soft ML. That is an observed relative reduction of about 69%. The corresponding 95% Wilson intervals are `[0.03346, 0.03765]` and `[0.00985, 0.01221]`. This is a result for the stated finite simulation, not a universal coding-gain claim.

## Quick start

Python 3.11 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Validate the mathematics and implementation:

```bash
make test
```

Run the configured Monte Carlo experiment and rebuild the figure:

```bash
make classical
```

The main outputs are:

- `results/data/classical_fer.csv`: error counts, FER estimates, and confidence intervals;
- `results/data/classical_run_metadata.json`: seed, configuration, and environment metadata;
- `results/figures/classical_fer_vs_ebn0.png`: the Week 1 comparison plot.

To make a fast smoke run while editing:

```bash
python -m experiments.run_classical --config configs/classical.toml --trials 1000
```

## Mathematical conventions

The matrices are systematic:

```text
G = [I4 | P]
H = [P^T | I3]
```

All arithmetic in encoding, parity checks, and syndromes is over GF(2). The implementation verifies

```text
H G^T = 0  (mod 2).
```

For unit-energy BPSK symbols and code rate `R = k/n`, the real AWGN variance is

```text
sigma^2 = 1 / (2 R 10^(Eb/N0_dB / 10)).
```

This means `Eb` is energy per information bit. The uncoded baseline uses `R=1`; the Hamming paths use `R=4/7`.

## Validation included

- GF(2) row reduction, rank, null space, and multiplication checks;
- `H G^T = 0`;
- 16 unique codewords and minimum distance 3;
- exhaustive correction of every single-bit error for every codeword;
- noiseless hard and exact-soft decoder checks;
- channel and confidence-interval sanity checks.

## Scope and limitations

This Week 1 code is a classical baseline, not yet a quantum decoder. It uses independent Gaussian channel noise, exact enumeration of only 16 codewords, and software Monte Carlo timing. Week 2 will reuse the same parity-check structure to construct and validate the Steane [[7,1,3]] CSS code. No employer data, internal specifications, or proprietary implementation is used.

## Repository guide

See [WEEK1_GUIDE.md](WEEK1_GUIDE.md) for the concepts, intended outputs, daily milestones, and interpretation rules.

## Licensing status

No open-source license is granted yet. 
