# Week 3 guide: hard and soft ML decoding of noisy Steane syndromes

## 1. What Week 3 is designed to prove

Week 3 transfers the hard-versus-soft comparison from the classical Hamming experiment to stabilizer measurements of the Steane `[[7,1,3]]` code.

The milestone is designed to leave visible evidence of four abilities:

1. **Quantum measurement modeling:** convert ideal binary stabilizer outcomes into continuous noisy measurements.
2. **Decoder design:** distinguish a decoder that thresholds first from one that retains analog reliability information.
3. **Statistical reasoning:** derive a Gaussian maximum-likelihood rule under an explicit equal-prior hypothesis model.
4. **Scope control:** demonstrate one clear idea without prematurely adding multi-qubit priors, MAP inference, or stabilizer-coset aggregation.

The three-week progression is

```text
Week 1: noisy Hamming codeword
        -> hard decision or soft ML

Week 2: Pauli algebra
        -> stabilizers
        -> ideal Steane syndromes

Week 3: noisy Steane syndrome
        -> hard decision or soft ML
```

---

## 2. Deliberately restricted error model

The candidate set is

```text
{I, X1, ..., X7, Y1, ..., Y7, Z1, ..., Z7}.
```

It contains:

```text
1 identity hypothesis
7 single-qubit X errors
7 single-qubit Y errors
7 single-qubit Z errors
-------------------------
22 hypotheses in total.
```

All 22 hypotheses are assumed to be equally likely.

This restriction is important. The experiment asks only:

> Does retaining continuous stabilizer-measurement information improve the identification and correction of a no-error or single-qubit Pauli event compared with thresholding the same measurements first?

Because the hypotheses have equal priors,

```text
MAP = ML.
```

No physical error-rate prior is required for this milestone.

The experiment deliberately excludes arbitrary multi-qubit errors and therefore does not require coset-MAP decoding.

---

## 3. The six-bit Steane syndrome table

The repository uses the Hamming parity-check matrix

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

The complete Steane syndrome contains six bits.

This project orders them as

```text
[X-check responses | Z-check responses].
```

The first three bits are produced by the X-type stabilizer checks.

The final three bits are produced by the Z-type stabilizer checks.

X-type checks detect Z components:

```text
s_X = H_X e_Z.
```

Z-type checks detect X components:

```text
s_Z = H_Z e_X.
```

Therefore, if `H_i` is column `i` of `H`, the restricted error signatures are

```text
I   -> [000 | 000]

X_i -> [000 | H_i]

Z_i -> [H_i | 000]

Y_i -> [H_i | H_i].
```

A Y error has both X and Z components, so both syndrome families respond.

The seven columns of `H` are all distinct and nonzero. Consequently, all 22 restricted hypotheses have unique six-bit syndromes.

The hypothesis set and syndrome table are implemented in

```text
src/softqec/quantum_decoders.py.
```

---

## 4. Analog syndrome model

An ideal stabilizer measurement has eigenvalue `+1` or `-1`.

The binary syndrome convention is mapped to these eigenvalues as

```text
s_j = 0 -> mu_j = +1

s_j = 1 -> mu_j = -1.
```

Equivalently,

```text
mu_j = (-1)^s_j.
```

The observed continuous value is

```text
y_j = (-1)^s_j + epsilon_j,
```

where

```text
epsilon_j ~ Normal(0, sigma_m^2).
```

The parameter `sigma_m` controls syndrome-measurement uncertainty.

### Small measurement noise

When `sigma_m` is small, observations remain close to `+1` and `-1`.

The two conditional Gaussian distributions overlap very little.

The hard and soft decoders should therefore both recover almost every syndrome correctly.

### Intermediate measurement noise

As `sigma_m` increases, some observations approach zero or cross the threshold.

The sign alone becomes less reliable.

However, the analog magnitude still tells the decoder which measurements are confident and which are uncertain.

This is the main region in which soft information can help.

### Large measurement noise

When `sigma_m` is large, the Gaussian distributions overlap strongly.

Both decoders become unreliable.

The analog magnitudes can still contain useful information, but the measurement itself becomes increasingly ambiguous.

The channel is implemented in

```text
src/softqec/analog_syndrome.py.
```

This is a phenomenological one-shot model.

It is not calibrated to a specific superconducting, trapped-ion, neutral-atom, or photonic readout system.

---

## 5. Hard syndrome decoder

The hard decoder first thresholds every continuous measurement:

```text
y_j >= 0 -> s_hat_j = 0

y_j < 0 -> s_hat_j = 1.
```

It then compares the resulting six-bit string with the 22 valid syndrome signatures.

The selected error hypothesis is

```text
E_hat_hard
    =
argmin_E d_H(s_hat, s(E)),
```

where `d_H` is Hamming distance.

If the thresholded string exactly matches one of the 22 valid syndromes, that candidate has distance zero.

If noise produces one of the other

```text
64 - 22 = 42
```

possible six-bit strings, the decoder chooses the nearest valid signature.

When multiple candidates have the same minimum distance, the implementation uses a fixed first-index tie rule.

This rule makes the simulation reproducible, but it is also an explicit limitation of the hard baseline.

### Information lost during thresholding

Consider two observations:

```text
y = +0.01

y = +3.00.
```

Both become syndrome bit zero.

However, they do not provide the same amount of evidence:

- `+0.01` is very uncertain;
- `+3.00` is strong evidence for syndrome bit zero.

After thresholding, this difference is permanently lost.

---

## 6. Soft Gaussian maximum-likelihood decoder

The soft decoder does not threshold the observation.

For every candidate error `E`, it compares the measured vector `y` with the ideal analog syndrome signature `mu_E`.

Under independent Gaussian measurement noise,

```text
p(y | E)
    =
product over j of
p(y_j | s_j(E)).
```

Ignoring normalization factors that are identical for every candidate,

```text
p(y | E)
    proportional to
exp(
    -||y - mu_E||^2
    /
    (2 sigma_m^2)
).
```

Maximizing this likelihood is therefore equivalent to minimizing squared Euclidean distance:

```text
E_hat_soft
    =
argmin_E ||y - mu_E||^2.
```

Because all 22 error hypotheses have equal prior probability,

```text
MAP = ML.
```

The soft decoder therefore does not need a separate prior term.

This is the direct syndrome-level counterpart of the Week 1 soft Hamming decoder:

```text
Week 1:
continuous noisy codeword
    -> nearest valid BPSK codeword

Week 3:
continuous noisy syndrome
    -> nearest valid analog syndrome signature.
```

Both hard and soft decoders are implemented in

```text
src/softqec/quantum_decoders.py.
```

---

## 7. Correction success

After estimating an error, the decoder applies the estimated Pauli as the correction.

In general stabilizer quantum error correction, success means

```text
E C belongs to the stabilizer group,
```

where:

- `E` is the physical error;
- `C` is the selected correction.

The correction does not always need to reproduce the exact physical error.

Two errors differing by a stabilizer have the same effect on the logical state.

However, the Week 3 hypothesis set contains only:

```text
identity
or
one single-qubit Pauli error.
```

Two distinct hypotheses from this set differ by a Pauli operator of weight at most two.

Every non-identity Steane stabilizer has weight four.

Therefore, inside this restricted hypothesis set,

```text
success modulo stabilizers
    is equivalent to
correct hypothesis identification.
```

The automated tests verify this equivalence exhaustively for all

```text
22 x 22 = 484
```

possible actual-error and selected-correction pairs.

This restricted equivalence should not be generalized to arbitrary multi-qubit errors.

---

## 8. Reproducible experiment

The committed experiment configuration is stored in

```text
configs/quantum_single_pauli.toml.
```

It specifies

```text
seed = 20260920

trials per sigma_m = 30000

sigma_m =
0.10
0.20
0.35
0.50
0.75
1.00
1.25
1.50
2.00.
```

For every trial:

1. Select one of the 22 hypotheses uniformly.
2. Calculate its ideal six-bit Steane syndrome.
3. Convert the syndrome to a `+1/-1` analog signature.
4. Add independent Gaussian measurement noise.
5. Give the same observation to both decoders.
6. Apply hard threshold plus nearest-syndrome decoding.
7. Apply continuous Gaussian ML decoding.
8. Count correction success or failure.

The same observation is used for both decoders.

This paired comparison ensures that a difference in performance comes from decoder information use rather than different random samples.

Run the complete experiment with

```bash
make quantum
```

or

```bash
python -m experiments.run_quantum_single_pauli \
    --config configs/quantum_single_pauli.toml
```

For a fast smoke test:

```bash
make smoke-quantum
```

---

## 9. Measured quantities

The main metric is decoding failure rate:

```text
decoding failure rate
    =
number of failed corrections
/
number of trials.
```

The experiment records this separately for

```text
hard threshold + nearest syndrome

soft Gaussian ML.
```

For each rate, it also calculates a 95% Wilson confidence interval.

The Wilson interval is preferred to the simplest normal approximation because it remains meaningful near zero and one.

The experiment additionally records software runtimes for both decoder calls.

These Python runtimes describe only this implementation.

They are not hardware-decoder latency measurements.

---

## 10. Main Week 3 result

With the committed random seed and 30,000 trials per noise point, at

```text
sigma_m = 1.0
```

the observed decoding failure rates are

```text
hard threshold + nearest syndrome: 0.52933

soft Gaussian ML:                  0.45193.
```

The absolute reduction is approximately

```text
0.52933 - 0.45193 = 0.0774.
```

The observed relative reduction is approximately

```text
(0.52933 - 0.45193) / 0.52933
    = 0.146
    = 14.6%.
```

The corresponding 95% Wilson intervals are

```text
hard:
[0.52368, 0.53498]

soft:
[0.44631, 0.45757].
```

At very small `sigma_m`, both decoders recover the ideal signatures.

At intermediate noise, the curves separate because the soft decoder retains the analog magnitudes discarded by thresholding.

At very high noise, both decoders become unreliable, although soft ML continues to perform better under the stated model.

The result demonstrates the value of preserving analog syndrome reliability information.

It is not a universal quantum-coding gain claim.

---

## 11. Experiment outputs

The full experiment produces

```text
results/data/quantum_single_pauli_failure.csv
```

containing the failure counts, failure rates, confidence intervals, and decoder runtimes.

It also produces

```text
results/data/quantum_single_pauli_metadata.json
```

containing:

- the random seeds;
- the noise grid;
- the number of trials;
- the Python and NumPy versions;
- syndrome ordering;
- analog mapping;
- decoder definitions;
- correction-success convention;
- deliberately excluded features.

The final figure is

```text
results/figures/quantum_single_pauli_failure_vs_sigma.png.
```

---

## 12. Automated validation

The Week 3 test suite verifies:

```text
0 -> +1 and 1 -> -1 analog mapping

thresholding of ideal analog values

zero-noise channel behavior

fixed-seed channel reproducibility

invalid sigma detection

nonfinite observation detection

22 correctly ordered hypotheses

22 unique six-bit syndromes

X, Y, and Z syndrome structure

noiseless hard decoding of all 22 hypotheses

noiseless soft decoding of all 22 hypotheses

a deterministic case where analog reliability helps

all 484 restricted error-correction combinations.
```

Run all repository tests with

```bash
pytest -v
```

or

```bash
make test
```

The expected total after Week 3 is

```text
84 passed.
```

GitHub Actions also performs small smoke runs of both the Week 1 and Week 3 experiments.

---

## 13. Honest limitations

This milestone intentionally excludes

```text
arbitrary multi-qubit errors

physical data-qubit error probability p

weight-dependent error priors

nonuniform priors

MAP decoding

stabilizer-coset probability aggregation

repeated syndrome measurements

faulty syndrome-extraction circuits

gate errors

ancilla errors

leakage

correlated noise

hardware-calibrated readout distributions.
```

The hard decoder also uses a fixed tie-breaking rule when multiple valid signatures are equally close to a thresholded syndrome.

The experiment demonstrates the value of analog information inside a small controlled model.

It is not:

- a fault-tolerance threshold study;
- a scalable quantum decoder;
- a circuit-level QEC simulation;
- a hardware-real-time benchmark;
- a complete treatment of quantum-code degeneracy.

---

## 14. Future research-grade extension

The natural next extension is to:

1. enumerate broader physical error patterns;
2. assign a physical-noise prior;
3. group stabilizer-equivalent errors into cosets;
4. aggregate posterior probability across each coset;
5. select the most probable logical correction class.

That extension would introduce degeneracy-aware coset-MAP decoding.
