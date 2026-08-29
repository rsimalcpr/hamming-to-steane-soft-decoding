# Week 1 guide: ingredients, outputs, and interpretation

## 1. What Week 1 is designed to prove

Week 1 should leave visible evidence of four abilities:

1. **Algebra:** you can construct and validate a binary linear code rather than treating it as a black-box library object.
2. **Probability:** you understand the difference between a hard bit decision and a likelihood-based decision.
3. **Experiment design:** the two coded decoders receive the same messages and noisy observations, and `Eb/N0` is normalized correctly for code rate.
4. **Research engineering:** the result is reproducible from a saved configuration and is protected by deterministic tests.

The target is not simply a pretty FER curve. The real artifact is the chain of reasoning connecting matrices, syndromes, channel evidence, decoders, error counting, and uncertainty.

## 2. Algebraic ingredients

### 2.1 Binary field GF(2)

GF(2) contains only 0 and 1. Addition is XOR and multiplication is ordinary bit multiplication:

| a | b | a + b mod 2 | a b mod 2 |
|---|---|----------------|-----------|
| 0 | 0 | 0 | 0 |
| 0 | 1 | 1 | 0 |
| 1 | 0 | 1 | 0 |
| 1 | 1 | 0 | 1 |

Because `1 + 1 = 0`, subtracting a row is identical to adding/XORing it. The `gf2.py` module supplies multiplication, reduced row-echelon form, rank, null space, and matrix-pair validation.

### 2.2 Linear block-code parameters

The Hamming code has parameters `[n,k,d] = [7,4,3]`:

- `k=4`: four independent information bits;
- `n=7`: seven transmitted code bits;
- `n-k=3`: three redundant parity constraints;
- `R=k/n=4/7`: the code rate;
- `d=3`: any two valid codewords differ in at least three positions.

A minimum distance of three means that radius-one Hamming balls around codewords do not overlap. Therefore every pattern containing at most one flipped bit has one unambiguous nearest codeword.

### 2.3 Generator and parity-check matrices

The repository fixes a systematic convention:

```text
G = [I4 | P],       H = [P^T | I3].
```

An information row vector `u` becomes

```text
c = u G mod 2.
```

Every valid codeword must satisfy

```text
H c^T = 0 mod 2.
```

The compatibility condition `H G^T = 0` proves that every encoded message satisfies all parity checks. Since `G` is systematic, the first four bits of a valid codeword equal the information vector.

The columns of `H` are the seven distinct nonzero three-bit vectors. This gives every single-bit position a unique syndrome.

### 2.4 Hamming weight, distance, and enumeration

The Hamming weight `wt(c)` is the number of ones in `c`. The distance between two binary words is `wt(c1 + c2)`, where addition is XOR. For a linear code, the minimum distance equals the smallest nonzero codeword weight. The implementation enumerates all `2^4=16` messages and codewords and verifies that this minimum is three.

## 3. Hard decoding ingredients

### 3.1 Binary symmetric channel

The BSC is the simplest discrete model: each bit flips independently with probability `p`. It is useful for testing syndrome logic separately from modulation and Gaussian noise.

If `r = c + e mod 2`, then

```text
s = H r^T = H(c+e)^T = H e^T,
```

because `Hc^T=0`. The syndrome therefore depends on the error pattern, not on which codeword was sent.

For a one-bit error at position `i`, `s` is exactly column `i` of `H`. The decoder looks up that column, flips bit `i`, and returns a codeword.

### 3.2 What the exhaustive test proves

There are 16 codewords and 7 possible one-bit locations, so the central test checks all 112 codeword/error pairs. Each corrected result must equal the transmitted codeword. This is much stronger than showing one example and is the main Week 1 correctness certificate.

The decoder is not guaranteed to correct two errors. With distance three, a two-bit pattern can be miscorrected to a different valid codeword.

## 4. Soft-decoding ingredients

### 4.1 BPSK and AWGN

Each bit becomes a real symbol:

```text
x_i = (-1)^c_i,
```

so `0 -> +1` and `1 -> -1`. The receiver observes

```text
y_i = x_i + n_i,       n_i ~ Normal(0, sigma^2).
```

For unit symbol energy and energy per information bit `Eb`, the code-rate-aware variance is

```text
sigma^2 = 1 / (2 R gamma_b),
gamma_b = 10^(Eb/N0_dB / 10).
```

Documenting `R` is essential: using the same noise variance for coded and uncoded transmissions while labeling the axis `Eb/N0` would compare unequal information-bit energies.

### 4.2 Hard AWGN path

Hard demodulation keeps only the sign:

```text
r_i = 0 if y_i >= 0 else 1.
```

The magnitude is thrown away. A sample at `y=+0.01` and one at `y=+3.0` both become bit zero even though the second is much more reliable. The thresholded seven-bit word then goes through the syndrome lookup decoder.

### 4.3 Exact soft ML path

There are only 16 valid codewords, so the decoder can compare the observation with all of them. With equal message priors and common Gaussian variance,

```text
c_hat = argmin_c ||y - x(c)||^2,
```

where `x(c)` is the BPSK image of codeword `c`. Since all BPSK codewords have equal energy, this is also

```text
c_hat = argmax_c y dot x(c).
```

This decoder uses sign and magnitude together. It is exact for the stated code, channel, and equal-prior model. It is not a scalable strategy for long codes because the candidate list grows as `2^k`.

## 5. What is measured

### 5.1 Frame-error rate

A frame is wrong when at least one of the four recovered information bits differs from the transmitted message:

```text
FER = number of wrong frames / number of transmitted frames.
```

The plot contains:

- uncoded hard decision on four-bit frames;
- Hamming hard syndrome decoding;
- Hamming exact soft ML decoding.

Bit-error rate is deliberately not the primary metric because the code operates on frames and the Week 1 objective is to compare complete decoder outcomes.

### 5.2 Confidence intervals

A Monte Carlo FER is an estimate. The code reports a 95% Wilson binomial interval instead of pretending the plotted point is exact. At very low FER, increase the trial count; zero observed failures does not prove a zero true failure probability.

### 5.3 Reproducibility

The configuration records the `Eb/N0` grid, number of frames, top-level random seed, and output paths. A deterministic child seed is derived for every grid point. Hard and soft Hamming decoding use the same coded observations at each point, preventing channel-sample differences from being mistaken for decoder differences.

## 6. Concrete Week 1 outputs

| Output | What it demonstrates | Acceptance check |
|---|---|---|
| `src/softqec/gf2.py` | Binary linear algebra | Rank/null-space tests pass |
| `src/softqec/classical_codes.py` | Hamming construction | `HG^T=0`, 16 words, `d=3` |
| `src/softqec/classical_channels.py` | BSC and rate-aware BPSK/AWGN | Channel sanity tests pass |
| `src/softqec/classical_decoders.py` | Hard lookup and exact soft ML | Noiseless and exhaustive tests pass |
| `experiments/run_classical.py` | Reproducible paired simulation | Same config and seed reproduce counts |
| `results/data/classical_fer.csv` | Raw evidence | Counts, FER, and intervals retained |
| `results/figures/classical_fer_vs_ebn0.png` | Main classical result | Three labeled curves, readable axes |
| Notebook | Explanatory walkthrough | Calls package functions; no duplicate hidden implementation |
| README | Reviewer-facing entry point | Clean setup, one-command test and experiment |

## 7. Daily implementation map

### Monday: definitions and scope

Write down `[n,k,d]`, rate, weight, distance, `G`, `H`, syndrome, hard information, and soft information. Freeze the Week 1 claim and do not add quantum code yet.

### Tuesday: GF(2)

Implement row reduction, rank, null space, multiplication, and validation. Tests should include both correct examples and rejected non-binary or shape-incompatible inputs.

### Wednesday: Hamming code

Fix `G` and `H`, encode all messages, enumerate codewords, verify `HG^T=0`, and compute the minimum distance.

### Thursday: hard decoding

Implement BSC sampling, syndrome lookup, exhaustive one-bit correction, and a smoke run.

### Saturday: soft decoding

Implement BPSK, rate-aware AWGN, thresholding, exact ML over 16 codewords, and vectorized Monte Carlo. Run a small pilot first and inspect whether the ordering of curves is physically plausible.

### Sunday: release

Run all tests, generate the saved CSV and figure, execute the notebook, check setup from a fresh environment, update the README with evidence-based observations, and mark the Week 1 milestone.

## 8. How to interpret the expected figure

The likely ordering over a useful middle-to-high `Eb/N0` region is:

```text
Hamming soft ML FER < Hamming hard FER,
```

because the soft decoder has access to all information available to the hard decoder plus reliability magnitude. Coding gain relative to uncoded transmission may not appear uniformly at very low `Eb/N0`; the extra transmitted symbols and frequent multi-bit corruption can create crossings. Report the observed regime rather than forcing a universal claim.

If soft ML looks worse than hard decoding with enough trials, first check the bit-to-symbol convention, code-rate factor in `sigma`, codeword ordering, and whether both paths used the same observations.

## 9. Week 1 completion gate

Week 1 is complete when a reviewer can install the package in a clean environment, run one documented command, reproduce the classical FER figure, and run tests proving correction of all 112 single-bit cases. The README must state the channel model, `Eb/N0` convention, FER definition, seed policy, and limitations.
