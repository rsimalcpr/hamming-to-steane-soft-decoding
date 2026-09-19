# Week 2 guide: stabilizers, CSS codes, and the Steane [[7,1,3]] code

## 1. What Week 2 is designed to prove

Week 2 extends the classical coding foundation developed in Week 1 into quantum error correction.

The milestone is designed to leave visible evidence of five abilities:

1. **Pauli algebra:** represent multi-qubit Pauli operators efficiently and determine commutation without constructing exponentially large matrices.
2. **Stabilizer reasoning:** understand how commutation and anticommutation generate quantum error syndromes.
3. **CSS construction:** translate binary parity-check matrices into compatible X-type and Z-type quantum checks.
4. **Logical reasoning:** distinguish physical errors, stabilizers, logical operators, and logical failure.
5. **Verification:** exhaustively validate every single-qubit Pauli error for the Steane [[7,1,3]] code.

The main objective is not merely to reproduce a known seven-qubit code.

The important result is the conceptual and computational bridge

```text
classical binary linear algebra
        ↓
Pauli symplectic representation
        ↓
stabilizer syndromes
        ↓
CSS construction
        ↓
Steane [[7,1,3]] code
        ↓
logical error correction
```

Week 1 and Week 2 therefore form one continuous project rather than two independent exercises.

---

## 2. Reusing the Week 1 GF(2) foundation

Week 1 introduced binary linear algebra over GF(2) through

```text
src/softqec/gf2.py
```

including:

- binary-array validation;
- addition as XOR;
- matrix multiplication modulo two;
- reduced row-echelon form;
- rank;
- null space;
- binary-vector enumeration.

Week 2 deliberately reuses these utilities.

This is important because the same mathematical operations that appear in classical coding,

```text
H G^T = 0 mod 2
```

and

```text
s = H e
```

also appear in quantum coding as

```text
H_X H_Z^T = 0 mod 2
```

and

```text
s_Z = H_Z e_X
s_X = H_X e_Z.
```

The quantum layer therefore builds on the classical algebra instead of introducing a second independent mathematical implementation.

---

## 3. Pauli operators

### 3.1 Single-qubit Pauli operators

The four single-qubit Pauli operators are

```text
I, X, Y, Z.
```

For stabilizer calculations, an overall global phase does not affect commutation or syndrome extraction.

A Pauli operator can therefore be represented by two binary values:

| Pauli | x | z |
|---|---:|---:|
| I | 0 | 0 |
| X | 1 | 0 |
| Z | 0 | 1 |
| Y | 1 | 1 |

The Y operator contains both an X and a Z component because, up to phase,

```text
Y = i X Z.
```

Ignoring the phase gives the binary correspondence

```text
Y ↔ XZ ↔ (1,1).
```

### 3.2 Multi-qubit symplectic representation

An n-qubit Pauli operator is represented as

```text
(x | z),
```

where `x` and `z` are binary vectors of length `n`.

For example,

```text
X Z I
```

is represented by

```text
x = [1,0,0]
z = [0,1,0].
```

This is much more efficient than explicitly constructing a `2^n × 2^n` operator matrix.

The implementation is contained in

```text
src/softqec/pauli.py
```

and the corresponding tests are in

```text
tests/test_pauli.py.
```

---

## 4. Binary symplectic product

For two Pauli operators

```text
P1 = (x1 | z1)
P2 = (x2 | z2),
```

their binary symplectic product is

```text
<P1,P2>_sp
    = x1 · z2 + z1 · x2 mod 2.
```

The result determines their commutation relation:

```text
0 -> commute
1 -> anticommute.
```

### Example: X and Z

For

```text
X = (1|0)
Z = (0|1),
```

the symplectic product is

```text
1*1 + 0*0 = 1 mod 2.
```

Therefore X and Z anticommute.

### Example: XX and ZZ

Consider

```text
XX
ZZ.
```

There are two positions at which X meets Z.

Each position individually contributes an anticommutation, but

```text
(-1)^2 = +1.
```

The two-qubit operators therefore commute.

The symplectic product counts these X/Z overlaps modulo two.

An even number of anticommuting overlaps produces zero.

An odd number produces one.

---

## 5. Stabilizers

### 5.1 Stabilizer definition

A state `|psi>` is stabilized by an operator `S` if

```text
S |psi> = |psi>.
```

A stabilizer code is defined as the simultaneous `+1` eigenspace of a collection of mutually commuting Pauli operators.

These operators are called stabilizer generators.

For an encoded state `|psi_L>`,

```text
S_j |psi_L> = |psi_L>
```

for every stabilizer generator `S_j`.

---

## 6. How an error creates a syndrome

Suppose an error `E` occurs.

The state becomes

```text
E |psi_L>.
```

There are two possible relationships between the error and a stabilizer generator.

### Commuting case

If

```text
S E = E S,
```

then

```text
S E|psi_L>
    = E S|psi_L>
    = E|psi_L>.
```

The corrupted state remains a `+1` eigenstate of the stabilizer.

The syndrome bit is therefore encoded as

```text
0.
```

### Anticommuting case

If

```text
S E = -E S,
```

then

```text
S E|psi_L>
    = -E S|psi_L>
    = -E|psi_L>.
```

The corrupted state is now a `-1` eigenstate.

The syndrome bit is encoded as

```text
1.
```

Therefore a stabilizer syndrome records which stabilizer generators anticommute with the error.

---

## 7. Stabilizer syndrome in symplectic form

Let a stabilizer generator be

```text
S = (S_X | S_Z)
```

and the error be

```text
E = (e_X | e_Z).
```

The corresponding syndrome bit is

```text
s = S_X e_Z + S_Z e_X mod 2.
```

For many stabilizer generators, this becomes a vector operation.

The implementation in

```text
src/softqec/stabilizer.py
```

reuses the Week 1 GF(2) matrix multiplication and addition functions.

This makes the classical-to-quantum connection explicit.

---

## 8. Repetition-code sanity check

Before constructing the Steane code, the stabilizer implementation is tested using the three-qubit bit-flip repetition code.

Its stabilizers are

```text
S1 = Z Z I
S2 = I Z Z.
```

Consider the three possible single-qubit X errors.

For

```text
X I I
```

the first stabilizer anticommutes with the error while the second commutes:

```text
XII -> 10.
```

Similarly,

```text
IXI -> 11
IIX -> 01.
```

The error location can therefore be identified from the syndrome.

However,

```text
ZII -> 00.
```

A Z error commutes with these Z-type stabilizers.

The repetition code therefore illustrates why one family of checks is insufficient for general quantum errors.

This motivates CSS codes.

---

## 9. CSS codes

CSS stands for Calderbank-Shor-Steane.

A CSS stabilizer code separates stabilizer checks into two families:

```text
H_X
```

and

```text
H_Z.
```

Rows of `H_X` describe X-type stabilizer generators.

Rows of `H_Z` describe Z-type stabilizer generators.

The implementation is contained in

```text
src/softqec/css.py.
```

---

## 10. Why CSS codes require H_X H_Z^T = 0

All stabilizer generators must commute.

Two X-type generators automatically commute.

Two Z-type generators also automatically commute.

The nontrivial condition is therefore between an X-type generator and a Z-type generator.

Suppose a row of `H_X` specifies the positions of X operators and a row of `H_Z` specifies the positions of Z operators.

Their binary dot product counts the number of positions at which an X and Z overlap.

Each such overlap contributes

```text
XZ = -ZX.
```

If the number of overlaps is even,

```text
(-1)^even = +1,
```

so the two multi-qubit operators commute.

If the number is odd,

```text
(-1)^odd = -1,
```

so they anticommute.

Therefore the CSS commutation condition is

```text
H_X H_Z^T = 0 mod 2.
```

This condition is not an arbitrary matrix identity.

It is the binary-linear-algebra statement that every X-type stabilizer commutes with every Z-type stabilizer.

---

## 11. Which checks detect which errors?

Suppose an X error occurs.

An X operator commutes with another X operator, so X-type stabilizers do not detect it directly.

However,

```text
XZ = -ZX.
```

Therefore X errors are detected by Z-type stabilizers:

```text
s_Z = H_Z e_X.
```

Likewise, Z errors are detected by X-type stabilizers:

```text
s_X = H_X e_Z.
```

A Y error contains both components because

```text
Y ~ XZ.
```

Therefore Y errors can activate both syndrome families.

A useful memory rule is

```text
X error <-> Z checks
Z error <-> X checks.
```

---

## 12. Number of logical qubits

An unconstrained n-qubit Hilbert space has dimension

```text
2^n.
```

Every independent stabilizer condition halves the allowed subspace.

With `r` independent stabilizers,

```text
dimension(code space) = 2^(n-r).
```

If the code encodes `k` logical qubits, its code space has dimension

```text
2^k.
```

Therefore

```text
k = n - r.
```

For a CSS code,

```text
r = rank(H_X) + rank(H_Z),
```

assuming the check rows are independent.

Hence

```text
k = n - rank(H_X) - rank(H_Z).
```

The Week 2 implementation uses the same GF(2) `rank()` routine developed in Week 1.

---

## 13. From the Hamming code to the Steane code

### 13.1 Classical starting point

The classical Hamming code has parameters

```text
[7,4,3].
```

This means

```text
n = 7
k = 4
d = 3.
```

A three-row parity-check matrix distinguishes the seven possible single-bit error locations through its seven distinct nonzero columns.

Week 2 uses a parity-check representation of the Hamming [7,4,3] code suitable for the CSS construction:

```text
H =
[1 0 1 0 1 0 1
 0 1 1 0 0 1 1
 0 0 0 1 1 1 1].
```

This Hamming representation is equivalent to other valid Hamming parity-check conventions up to coordinate ordering.

It is therefore not required to be literally identical to the systematic matrix representation used for the Week 1 classical encoder.

---

## 14. Steane CSS construction

For the Steane code,

```text
H_X = H_Z = H.
```

The repository verifies

```text
H H^T = 0 mod 2.
```

Therefore the corresponding X-type and Z-type stabilizer generators commute.

The ranks are

```text
rank(H_X) = 3
rank(H_Z) = 3.
```

Since

```text
n = 7,
```

the number of encoded logical qubits is

```text
k
= 7 - 3 - 3
= 1.
```

The Steane code has distance

```text
d = 3.
```

Its parameters are therefore

```text
[[7,1,3]].
```

The notation means:

```text
7 physical qubits
1 logical qubit
distance 3.
```

Distance three allows correction of arbitrary single-qubit Pauli errors.

---

## 15. Steane stabilizer generators

The first row of `H`,

```text
1 0 1 0 1 0 1,
```

defines the X-type generator

```text
X I X I X I X.
```

The second row defines

```text
I X X I I X X.
```

The third defines

```text
I I I X X X X.
```

The corresponding Z-type generators are

```text
Z I Z I Z I Z
I Z Z I I Z Z
I I I Z Z Z Z.
```

The code therefore contains six independent stabilizer generators:

```text
3 X-type
3 Z-type.
```

These six independent constraints reduce the seven-qubit physical Hilbert space to a two-dimensional logical code space.

---

## 16. Logical operators

A logical Pauli operator must satisfy two important conditions:

1. it commutes with every stabilizer generator;
2. it is not itself a member of the stabilizer group.

A convenient choice for the Steane code is

```text
X_bar = X X X X X X X
Z_bar = Z Z Z Z Z Z Z.
```

Both operators commute with all six stabilizers.

However, they overlap as X and Z on all seven qubits.

Therefore

```text
(-1)^7 = -1,
```

so

```text
X_bar Z_bar = - Z_bar X_bar.
```

This reproduces the ordinary Pauli anticommutation relation at the logical-qubit level.

The repository tests all of these properties.

---

## 17. Ideal Steane syndrome decoding

Every nonzero three-bit syndrome corresponds to one of the seven columns of `H`.

For an X error,

```text
s_Z = H_Z e_X
```

identifies the error location.

For a Z error,

```text
s_X = H_X e_Z
```

does the same.

For a Y error, both components are present:

```text
e_X != 0
e_Z != 0.
```

Both syndrome families are therefore used.

The ideal Week 2 decoder uses these syndromes to apply minimum-weight single-qubit corrections.

Syndrome extraction is assumed to be perfect at this stage.

Measurement uncertainty is intentionally postponed until Week 3.

---

## 18. Physical error versus logical error

An important quantum-specific idea is that successful decoding does not require the correction to exactly reproduce the physical error.

Suppose the actual error is

```text
E
```

and the decoder applies

```text
C.
```

The residual operator is

```text
E C.
```

If

```text
E C = I,
```

the correction obviously succeeds.

But suppose instead

```text
E C = S,
```

where `S` is a stabilizer.

For an encoded state,

```text
S |psi_L> = |psi_L>.
```

Therefore

```text
E C |psi_L>
    = S |psi_L>
    = |psi_L>.
```

The logical information has still been restored.

The correct success criterion is therefore

```text
E C belongs to the stabilizer group.
```

---

## 19. Stabilizer equivalence and degeneracy

Two physical Pauli errors may differ by a stabilizer but have exactly the same effect on the encoded logical state.

If

```text
E2 = E1 S,
```

for some stabilizer `S`, then `E1` and `E2` belong to the same logical equivalence class.

This phenomenon is known as quantum-code degeneracy.

It becomes especially important in Week 3.

A probabilistic decoder should not necessarily ask

```text
Which single physical error pattern is most probable?
```

Instead, for logical decoding it is often more meaningful to ask

```text
Which stabilizer-equivalence class has the largest total probability?
```

This motivates the later exact coset-MAP decoder.

---

## 20. Exhaustive single-qubit validation

The Steane code has distance three.

It must therefore correct every weight-one Pauli error.

There are seven physical qubits and three possible non-identity single-qubit Paulis:

```text
X
Y
Z.
```

Therefore the test set contains

```text
7 X errors
7 Y errors
7 Z errors
---------
21 total.
```

The test suite explicitly generates every one of these errors.

For each case it:

1. constructs the binary X and Z error components;
2. calculates the ideal CSS syndrome;
3. obtains the corresponding correction;
4. calculates the residual error;
5. verifies that the residual belongs to the stabilizer group.

This exhaustive validation is stronger than demonstrating only a few selected examples.

---

## 21. Week 2 implementation files

| File | Purpose |
|---|---|
| `src/softqec/gf2.py` | Shared GF(2) linear algebra from Week 1 |
| `src/softqec/pauli.py` | Pauli binary symplectic representation |
| `src/softqec/stabilizer.py` | Stabilizer commutation and syndrome calculation |
| `src/softqec/css.py` | CSS checks, commutation, and X/Z syndromes |
| `src/softqec/steane.py` | Steane [[7,1,3]] construction and ideal decoder |
| `tests/test_pauli.py` | Pauli commutation and representation tests |
| `tests/test_stabilizer.py` | Stabilizer and repetition-code syndrome tests |
| `tests/test_css.py` | CSS compatibility and syndrome tests |
| `tests/test_steane.py` | Logical operators and exhaustive Steane validation |
| `notebooks/02_hamming_to_steane.ipynb` | Reviewer-facing classical-to-quantum walkthrough |

---

## 22. Dependency structure

The implementation follows this conceptual dependency chain:

```text
gf2.py
   |
   +------> pauli.py
   |           |
   |           +------> stabilizer.py
   |
   +------> css.py
                 |
                 +------> steane.py
```

The notebook does not contain a second hidden implementation.

It imports and demonstrates the same package functions that are covered by automated tests.

---

## 23. Automated validation

The repository test suite checks:

```text
Pauli string conversion
Pauli weight
known commuting Pauli pairs
known anticommuting Pauli pairs
stabilizer commutation
repetition-code syndromes
CSS commutation
CSS logical-qubit count
X-error syndromes
Z-error syndromes
Steane H H^T = 0
Steane check ranks
Steane k = 1
logical X validity
logical Z validity
logical X/Z anticommutation
all seven nonzero single-qubit syndromes
all 7 single-qubit X errors
all 7 single-qubit Y errors
all 7 single-qubit Z errors.
```

GitHub Actions runs the automated test workflow on committed changes.

---

## 24. What Week 2 does not claim

The Week 2 implementation is an algebraic and decoding baseline.

It does not yet model:

```text
noisy syndrome measurements
repeated syndrome rounds
ancilla faults
gate faults
circuit-level propagation
leakage
correlated errors
hardware-calibrated readout distributions
large-code decoding
real-time hardware latency.
```

The Steane code is used because it is small enough for transparent verification while providing a direct bridge from classical Hamming coding to quantum stabilizer coding.

It is not presented as a scalable QEC performance benchmark.

---

## 25. Connection to Week 1

Week 1 studied

```text
information bits
    ↓
Hamming encoding
    ↓
channel noise
    ↓
syndrome / likelihood decoding.
```

Week 2 studies

```text
Pauli errors
    ↓
stabilizer checks
    ↓
CSS syndromes
    ↓
logical correction.
```

The common mathematical language is GF(2).

In Week 1,

```text
s = H e.
```

In Week 2,

```text
s_Z = H_Z e_X
s_X = H_X e_Z.
```

This is the central classical-to-quantum bridge of the repository.

---

## 26. Connection to Week 3

Week 2 assumes an exact binary syndrome.

Week 3 will instead model each syndrome bit as a hidden binary variable observed through a noisy continuous measurement.

The planned phenomenological model is

```text
y_j = (-1)^s_j + epsilon_j

epsilon_j ~ Normal(0, sigma_m^2).
```

A hard decoder will first threshold

```text
y_j
```

into a binary syndrome.

A soft decoder will retain the analog reliability information.

The Week 2 concept of stabilizer equivalence will then become essential.

Rather than selecting only the most likely physical error pattern, the planned soft decoder will aggregate posterior probability over stabilizer-equivalent errors.

This leads to the Week 3 comparison

```text
hard threshold
    +
syndrome lookup

versus

analog syndrome likelihood
    +
exact coset-MAP decoding.
```
