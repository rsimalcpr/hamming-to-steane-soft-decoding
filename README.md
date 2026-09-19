# From Hamming to Steane: Hard and Soft Decoding

This repository develops a transparent classical-to-quantum error-correction study.

The completed Week 1 milestone implements the classical Hamming [7,4,3] foundation and compares hard syndrome decoding with exact soft maximum-likelihood (ML) decoding over a BPSK/AWGN channel.

The Week 2 milestone extends the same binary coding ideas into quantum error correction through Pauli symplectic algebra, stabilizer syndromes, CSS check matrices, and the Steane [[7,1,3]] code. The implementation validates logical operators and exhaustively tests all 21 single-qubit Pauli errors.

The next milestone will introduce noisy analog syndrome measurements and compare hard-threshold decoding with soft-information decoding.

---

## Research questions

The repository is organized around a connected sequence of questions.

### Classical question

How much frame-error-rate improvement is obtained when a Hamming decoder retains the continuous BPSK/AWGN observations instead of thresholding them immediately?

### Classical-to-quantum bridge

How do parity checks, binary syndromes, and GF(2) linear algebra from the Hamming [7,4,3] code extend into stabilizer and CSS quantum error correction?

### Planned soft-QEC question

How much useful information is lost when noisy continuous stabilizer measurements are thresholded into binary syndrome bits before decoding?

---

# Week 1: classical Hamming decoding

## Week 1 pipeline

1. Encode four information bits with a systematic Hamming [7,4,3] generator matrix.
2. Map code bits to BPSK symbols using `0 -> +1` and `1 -> -1`.
3. Add Gaussian noise using a documented information-bit `Eb/N0` convention.
4. Decode the same coded observations in two ways:
   - hard path: threshold, compute syndrome, and apply a single-bit lookup correction;
   - soft path: compare the observation with all 16 valid BPSK codewords and choose exact ML.
5. Estimate frame-error rates and 95% Wilson confidence intervals.

The classical figure also includes an uncoded four-bit hard-decision baseline at the same `Eb/N0` definition.

## Current classical result

With the committed seed and 30,000 frames per point, at `Eb/N0 = 4 dB` the observed Hamming FER is

```text
hard syndrome decoding: 0.0355
exact soft ML:          0.01097
```

This corresponds to an observed relative FER reduction of about 69%.

The corresponding 95% Wilson intervals are

```text
hard syndrome decoding: [0.03346, 0.03765]
exact soft ML:          [0.00985, 0.01221]
```

This is a finite Monte Carlo result for the documented channel, code-rate convention, seed, and simulation size.

It is not presented as a universal coding-gain claim.

---

# Week 2: Hamming to Steane

Week 2 builds the quantum-coding layer through the following progression:

```text
Pauli operators
        ↓
binary symplectic representation
        ↓
stabilizers and quantum syndromes
        ↓
CSS check matrices
        ↓
Steane [[7,1,3]] code
        ↓
logical error correction
```

## Binary symplectic Pauli representation

An n-qubit Pauli operator is represented by two binary vectors

```text
(x | z)
```

using the mapping

```text
I -> (0,0)
X -> (1,0)
Z -> (0,1)
Y -> (1,1).
```

Two Pauli operators

```text
P1 = (x1 | z1)
P2 = (x2 | z2)
```

commute when their symplectic product is zero:

```text
x1 · z2 + z1 · x2 = 0 mod 2.
```

A result of one means that they anticommute.

This allows Pauli commutation to be handled using binary linear algebra rather than exponentially large operator matrices.

---

## Stabilizer syndromes

An encoded state satisfies

```text
S_j |psi_L> = |psi_L>
```

for every stabilizer generator `S_j`.

If a Pauli error commutes with a stabilizer, its syndrome bit is zero.

If it anticommutes, its syndrome bit is one.

For

```text
S = (S_X | S_Z)
```

and

```text
E = (e_X | e_Z),
```

the syndrome relation is

```text
s = S_X e_Z + S_Z e_X mod 2.
```

The Week 2 implementation reuses the GF(2) utilities developed during Week 1.

---

## CSS construction

CSS codes separate stabilizer checks into two binary matrices:

```text
H_X
H_Z.
```

The required commutation condition is

```text
H_X H_Z^T = 0 mod 2.
```

This condition means that every X-type stabilizer overlaps every Z-type stabilizer on an even number of qubits and therefore commutes with it.

X errors are detected by Z-type checks:

```text
s_Z = H_Z e_X.
```

Z errors are detected by X-type checks:

```text
s_X = H_X e_Z.
```

A Y error contains both X and Z components and can therefore activate both syndrome families.

---

## Steane [[7,1,3]] construction

Week 2 uses a parity-check representation of the classical Hamming [7,4,3] code:

```text
H =
[1 0 1 0 1 0 1
 0 1 1 0 0 1 1
 0 0 0 1 1 1 1].
```

For the Steane CSS construction,

```text
H_X = H_Z = H.
```

The implementation verifies

```text
H_X H_Z^T = 0 mod 2
rank(H_X) = 3
rank(H_Z) = 3.
```

The number of logical qubits is therefore

```text
k
= n - rank(H_X) - rank(H_Z)
= 7 - 3 - 3
= 1.
```

Together with distance three, the code has parameters

```text
[[7,1,3]].
```

The implementation also constructs six stabilizer generators:

```text
3 X-type generators
3 Z-type generators.
```

---

## Logical operators

A convenient choice of logical Pauli operators is

```text
X_bar = X X X X X X X
Z_bar = Z Z Z Z Z Z Z.
```

The implementation verifies that both logical operators:

```text
commute with every stabilizer
```

while neither belongs to the stabilizer group.

It also verifies

```text
X_bar Z_bar = - Z_bar X_bar,
```

as required for the logical Pauli operators of the encoded qubit.

---

## Exhaustive single-qubit validation

A distance-three quantum code must correct every weight-one Pauli error.

For seven physical qubits this gives

```text
7 X errors
7 Y errors
7 Z errors
---------
21 total.
```

The automated Week 2 tests generate every one of these errors and verify successful logical correction.

Success is defined modulo the stabilizer group.

If the physical error is `E` and the decoder correction is `C`, decoding succeeds whenever

```text
E C
```

belongs to the stabilizer group.

The correction therefore does not need to reproduce the exact physical error pattern.

This stabilizer equivalence will become important for the Week 3 soft coset-MAP decoder.

---

# Classical-to-quantum dictionary

| Classical coding | Stabilizer / CSS quantum coding |
|---|---|
| GF(2) binary algebra | GF(2) symplectic algebra |
| parity-check matrix `H` | CSS check matrices `H_X`, `H_Z` |
| binary error vector `e` | Pauli components `(e_X \| e_Z)` |
| syndrome `H e` | `s_Z = H_Z e_X`, `s_X = H_X e_Z` |
| parity constraints | stabilizer measurements |
| single-bit error location | single-qubit Pauli error location |
| decoding to the transmitted codeword | restoring the encoded logical state |
| exact physical error | stabilizer-equivalent error class |

---

# Quick start

Python 3.11 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Run the full automated test suite:

```bash
make test
```

or

```bash
pytest -v
```

Run the configured Week 1 classical Monte Carlo experiment and rebuild its figure:

```bash
make classical
```

For a fast classical smoke run while editing:

```bash
python -m experiments.run_classical \
    --config configs/classical.toml \
    --trials 1000
```

---

# Mathematical conventions

All coding and stabilizer arithmetic is performed over GF(2).

## Classical Hamming code

The Week 1 systematic convention is

```text
G = [I4 | P]
H = [P^T | I3].
```

The implementation verifies

```text
H G^T = 0 mod 2.
```

For unit-energy BPSK symbols and code rate `R = k/n`, the real AWGN variance is

```text
sigma^2 = 1 / (2 R 10^(Eb/N0_dB / 10)).
```

Thus `Eb` refers to energy per information bit.

The uncoded baseline uses

```text
R = 1,
```

while the Hamming code uses

```text
R = 4/7.
```

## Stabilizer / CSS code

For CSS checks,

```text
H_X H_Z^T = 0 mod 2.
```

For a Pauli error

```text
(e_X | e_Z),
```

the two syndrome components are

```text
s_X = H_X e_Z
s_Z = H_Z e_X.
```

All operations above are evaluated modulo two.

---

# Validation included

## Week 1 validation

The classical test suite covers:

```text
GF(2) row reduction
GF(2) rank
GF(2) null space
GF(2) matrix multiplication
H G^T = 0
16 unique Hamming codewords
minimum distance d = 3
all 112 codeword / single-bit-error combinations
noiseless hard decoding
noiseless exact-soft decoding
channel sanity checks
confidence-interval sanity checks.
```

## Week 2 validation

The quantum test suite covers:

```text
Pauli-string to symplectic conversion
known Pauli commuting cases
known Pauli anticommuting cases
Pauli weight
stabilizer commutation
repetition-code syndromes
CSS commutation
CSS logical-qubit count
X-error syndrome calculation
Z-error syndrome calculation
Steane H H^T = 0
rank(H_X) = 3
rank(H_Z) = 3
k = 1
six Steane stabilizer generators
logical X validity
logical Z validity
logical X/Z anticommutation
seven nonzero single-qubit syndromes
all seven X errors
all seven Y errors
all seven Z errors.
```

GitHub Actions runs the automated test workflow on committed changes.

---

# Repository structure

```text
hamming-to-steane-soft-decoding/
│
├── README.md
├── WEEK1_GUIDE.md
├── WEEK2_GUIDE.md
├── pyproject.toml
├── Makefile
│
├── src/
│   └── softqec/
│       ├── __init__.py
│       ├── gf2.py
│       ├── classical_codes.py
│       ├── classical_channels.py
│       ├── classical_decoders.py
│       ├── metrics.py
│       ├── pauli.py
│       ├── stabilizer.py
│       ├── css.py
│       └── steane.py
│
├── tests/
│   ├── test_gf2.py
│   ├── test_hamming.py
│   ├── test_channels_decoders.py
│   ├── test_metrics.py
│   ├── test_pauli.py
│   ├── test_stabilizer.py
│   ├── test_css.py
│   └── test_steane.py
│
├── notebooks/
│   ├── 01_classical_hamming.ipynb
│   └── 02_hamming_to_steane.ipynb
│
├── configs/
├── experiments/
└── results/
```

---

# Main outputs

## Week 1

```text
results/data/classical_fer.csv
results/data/classical_run_metadata.json
results/figures/classical_fer_vs_ebn0.png
WEEK1_GUIDE.md
```

These provide the classical hard-versus-soft baseline.

## Week 2

```text
src/softqec/pauli.py
src/softqec/stabilizer.py
src/softqec/css.py
src/softqec/steane.py

tests/test_pauli.py
tests/test_stabilizer.py
tests/test_css.py
tests/test_steane.py

notebooks/02_hamming_to_steane.ipynb
WEEK2_GUIDE.md
```

These provide the validated Hamming-to-Steane quantum-coding bridge.

---

# Scope and limitations

## Week 1

The classical experiment uses:

```text
independent Gaussian noise
BPSK modulation
Hamming [7,4,3]
exact enumeration over 16 codewords
finite Monte Carlo simulation.
```

Exact soft ML is practical here because the code is very small.

It is not intended as a scalable long-code decoder.

## Week 2

The Steane implementation assumes:

```text
ideal syndrome extraction
perfect binary syndrome information
independent algebraic Pauli-error reasoning
single-shot decoding.
```

It does not yet include:

```text
measurement noise
repeated syndrome rounds
gate faults
ancilla faults
circuit-level error propagation
leakage
correlated errors
hardware-calibrated readout distributions
large-code decoding
hardware decoder latency.
```

The Steane code is used as a transparent classical-to-quantum bridge.

It is not presented as a scalable QEC performance benchmark.

No employer data, internal specifications, proprietary source code, or proprietary implementation is used.

---

# Next milestone: soft syndrome decoding

Week 3 will replace the ideal binary syndrome assumption with a continuous noisy measurement model.

For syndrome bit `s_j`,

```text
y_j = (-1)^s_j + epsilon_j

epsilon_j ~ Normal(0, sigma_m^2).
```

The hard path will threshold `y_j` into a binary syndrome before lookup decoding.

The soft path will retain its analog reliability information.

For the small seven-qubit code, the planned exact decoder will enumerate physical error patterns and aggregate posterior probability over stabilizer-equivalent classes.

The central comparison will therefore become

```text
hard threshold
    +
minimum-weight syndrome lookup

versus

analog syndrome likelihood
    +
exact coset-MAP decoding.
```

The experiment will study the trade-off between logical failure rate and computational cost under the stated phenomenological noise model.

---

# Repository guides

Detailed explanations are available in:

- [WEEK1_GUIDE.md](WEEK1_GUIDE.md) — GF(2), Hamming [7,4,3], BPSK/AWGN, hard decoding, and exact soft ML.
- [WEEK2_GUIDE.md](WEEK2_GUIDE.md) — Pauli algebra, stabilizers, CSS construction, Steane [[7,1,3]], logical operators, and exhaustive single-qubit validation.

The explanatory notebook

```text
notebooks/02_hamming_to_steane.ipynb
```

provides a compact executable walkthrough of the Week 2 classical-to-quantum bridge.

---

# Licensing status

No open-source license is granted yet.
