"""Linear-algebra utilities over the binary field GF(2)."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

BinaryArray = NDArray[np.uint8]


def as_binary_array(values: ArrayLike, *, name: str = "array") -> BinaryArray:
    """Return a uint8 array after checking that every entry is 0 or 1."""
    array = np.asarray(values)
    if array.size and not np.all((array == 0) | (array == 1)):
        raise ValueError(f"{name} must contain only 0 and 1")
    return array.astype(np.uint8, copy=False)


def add(left: ArrayLike, right: ArrayLike) -> BinaryArray:
    """Add two arrays over GF(2), which is elementwise XOR."""
    a = as_binary_array(left, name="left")
    b = as_binary_array(right, name="right")
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch for GF(2) addition: {a.shape} != {b.shape}")
    return np.bitwise_xor(a, b)


def matmul(left: ArrayLike, right: ArrayLike) -> BinaryArray:
    """Multiply arrays using NumPy matmul semantics, then reduce modulo 2."""
    a = as_binary_array(left, name="left")
    b = as_binary_array(right, name="right")
    try:
        product = np.matmul(a.astype(np.int64), b.astype(np.int64))
    except ValueError as exc:
        raise ValueError(f"incompatible shapes for GF(2) multiplication: {a.shape}, {b.shape}") from exc
    return np.asarray(product % 2, dtype=np.uint8)


def rref(matrix: ArrayLike) -> tuple[BinaryArray, tuple[int, ...]]:
    """Compute reduced row-echelon form over GF(2)."""
    reduced = as_binary_array(matrix, name="matrix").copy()
    if reduced.ndim != 2:
        raise ValueError("matrix must be two-dimensional")

    rows, columns = reduced.shape
    pivot_row = 0
    pivots: list[int] = []

    for column in range(columns):
        candidates = np.flatnonzero(reduced[pivot_row:, column])
        if candidates.size == 0:
            continue

        selected = pivot_row + int(candidates[0])
        if selected != pivot_row:
            reduced[[pivot_row, selected]] = reduced[[selected, pivot_row]]

        for row in range(rows):
            if row != pivot_row and reduced[row, column]:
                reduced[row] ^= reduced[pivot_row]

        pivots.append(column)
        pivot_row += 1
        if pivot_row == rows:
            break

    return reduced, tuple(pivots)


def rank(matrix: ArrayLike) -> int:
    """Return matrix rank over GF(2)."""
    return len(rref(matrix)[1])


def null_space(matrix: ArrayLike) -> BinaryArray:
    """Return a row-basis for the right null space of a binary matrix."""
    reduced, pivots = rref(matrix)
    columns = reduced.shape[1]
    free_columns = [column for column in range(columns) if column not in pivots]
    basis = np.zeros((len(free_columns), columns), dtype=np.uint8)

    for basis_row, free_column in enumerate(free_columns):
        basis[basis_row, free_column] = 1
        for row, pivot_column in enumerate(pivots):
            basis[basis_row, pivot_column] = reduced[row, free_column]
    return basis


def validate_generator_parity_check(generator: ArrayLike, parity_check: ArrayLike) -> None:
    """Raise ValueError unless G and H are compatible full-rank code matrices."""
    g = as_binary_array(generator, name="generator")
    h = as_binary_array(parity_check, name="parity_check")
    if g.ndim != 2 or h.ndim != 2:
        raise ValueError("generator and parity_check must be two-dimensional")
    if g.shape[1] != h.shape[1]:
        raise ValueError("generator and parity_check must have the same block length")
    if rank(g) != g.shape[0]:
        raise ValueError("generator rows must be linearly independent")
    if rank(h) != h.shape[0]:
        raise ValueError("parity-check rows must be linearly independent")
    if g.shape[0] + h.shape[0] != g.shape[1]:
        raise ValueError("matrix dimensions do not define complementary code spaces")
    if np.any(matmul(h, g.T)):
        raise ValueError("H G^T must equal zero over GF(2)")


def binary_vectors(length: int) -> BinaryArray:
    """Enumerate all binary row vectors of a requested length."""
    if length < 0:
        raise ValueError("length must be non-negative")
    integers = np.arange(2**length, dtype=np.uint64)[:, None]
    shifts = np.arange(length - 1, -1, -1, dtype=np.uint64)
    return ((integers >> shifts) & 1).astype(np.uint8)


def hamming_weight(vectors: ArrayLike, *, axis: int = -1) -> NDArray[np.int64]:
    """Count nonzero binary entries along an axis."""
    values = as_binary_array(vectors, name="vectors")
    return np.sum(values, axis=axis, dtype=np.int64)
