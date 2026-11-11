#!/usr/bin/env python3
"""Expand affine difference-matrix and orthogonal-array certificates into MOLS.

All arithmetic is componentwise modulo the stated prime. Certificates contain ordinary
rectangular matrices and integer-encoded vectors. The common affine construction
is described in manuscript/mols_consolidated.tex.
"""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

def prime(n: int) -> bool:
    """Trial-division primality, for the small arithmetic bases used here."""
    return isinstance(n, int) and n >= 2 and all(n % k for k in range(2, math.isqrt(n)+1))


def digits(value: int, radices: list[int]) -> list[int]:
    """Decode a nonnegative integer into little-endian mixed-radix digits."""
    if value < 0:
        raise ValueError("Negative symbol encoding")
    result = []
    for radix in radices:
        result.append(value % radix)
        value //= radix
    if value:
        raise ValueError("Integer outside alphabet")
    return result


def encode(vector: list[int], radices: list[int]) -> int:
    """Encode a vector; its entries must already be reduced modulo the radices."""
    if len(vector) != len(radices):
        raise ValueError("Wrong vector length")
    result, weight = 0, 1
    for entry, radix in zip(vector, radices):
        if not 0 <= entry < radix:
            raise ValueError("Unreduced vector entry")
        result += entry * weight
        weight *= radix
    return result


def group_table(radices: list[int], sign: int = 1) -> list[list[int]]:
    """Build componentwise addition (sign=1) or subtraction (sign=-1).

    These small lookup tables accelerate expansion but do not change the
    arithmetic. In particular, this is NOT addition modulo the group order.
    """
    vectors = [digits(x, radices) for x in range(math.prod(radices))]
    table = []
    for left in vectors:
        row = []
        for right in vectors:
            value = [(x + sign * y) % p for x, y, p in zip(left, right, radices)]
            row.append(encode(value, radices))
        table.append(row)
    return table


def load(order: int, filename: Path) -> dict:
    """Read a certificate and validate its dimensions and symbol ranges."""
    cert = json.loads(filename.read_text())
    if cert["order"] != order:
        raise ValueError("Order mismatch")
    if cert["kind"] not in ("difference_matrix", "translated_array"):
        raise ValueError("Unknown construction type")
    if not cert["coarse_radices"] or any(not prime(p) for p in cert["coarse_radices"]):
        raise ValueError("Invalid coarse group")

    s = math.prod(cert["coarse_radices"])
    p, d, a = cert["p"], cert["d"], cert["a"]
    H, maps, offsets = cert["H"], cert["maps"], cert["offsets"]
    r = len(maps)
    if not prime(p) or not 0 < a < d or s * p**d != order:
        raise ValueError("Invalid lift parameters")
    if not H or len(offsets) != len(H) or any(len(row) != r for row in H + offsets):
        raise ValueError("Ragged construction data")
    if any(not 0 <= x < s for row in H for x in row):
        raise ValueError("Coarse symbol outside alphabet")
    if any(not 0 <= x < p**d for row in offsets for x in row):
        raise ValueError("Offset outside alphabet")
    for matrix in maps:
        if len(matrix) != d or any(len(row) != a for row in matrix):
            raise ValueError("Wrong matrix shape")
        if any(not 0 <= x < p for row in matrix for x in row):
            raise ValueError("Invalid matrix entry")
    extra = 1 if cert["kind"] == "difference_matrix" else 2
    if r != cert["mols"] + extra:
        raise ValueError("Wrong number of output columns")
    return cert


def affine_values(cert: dict) -> list[list[list[int]]]:
    """Return E[t][u][i] = encode(M_i * digits(u) + offset[t][i]).

    Matrix multiplication uses scalar arithmetic modulo p. The input index u
    ranges over all a-digit vectors; each output has d digits.
    """
    p, d, a = cert["p"], cert["d"], cert["a"]
    r = len(cert["maps"])
    radices = [p] * d
    inputs = [digits(u, [p] * a) for u in range(p**a)]
    images = []
    for matrix in cert["maps"]:
        column_images = []
        for vector in inputs:
            image = [sum(x * y for x, y in zip(row, vector)) % p for row in matrix]
            column_images.append(image)
        images.append(column_images)

    result = []
    for offsets in cert["offsets"]:
        vectors = [digits(value, radices) for value in offsets]
        block = []
        for u in range(p**a):
            row = []
            for i in range(r):
                value = [(x + y) % p for x, y in zip(images[i][u], vectors[i])]
                row.append(encode(value, radices))
            block.append(row)
        result.append(block)
    return result


def verify_compact(cert: dict) -> None:
    """Exhaustively check the paper's small difference-partition conditions."""
    p, d, a = cert["p"], cert["d"], cert["a"]
    q, H, r = p**d, cert["H"], len(cert["maps"])
    s = math.prod(cert["coarse_radices"])
    E = affine_values(cert)
    fine_subtract = group_table([p] * d, -1)
    coarse_subtract = group_table(cert["coarse_radices"], -1)
    is_dm = cert["kind"] == "difference_matrix"
    for i, j in itertools.combinations(range(r), 2):
        groups = {}
        for t, row in enumerate(H):
            key = coarse_subtract[row[i]][row[j]] if is_dm else (row[i], row[j])
            groups.setdefault(key, []).append(t)
        if len(groups) != (s if is_dm else s * s):
            raise ValueError("Missing coarse pair group")
        for key, source_rows in groups.items():
            if len(source_rows) != p**(d - a):
                raise ValueError("Wrong coarse multiplicity")
            seen = bytearray(q)
            for t in source_rows:
                for u in range(p**a):
                    value = fine_subtract[E[t][u][i]][E[t][u][j]]
                    if seen[value]:
                        raise ValueError(f"Repeated fine difference: columns {i},{j}, group {key}")
                    seen[value] = 1
            if not all(seen):
                raise ValueError("Incomplete fine difference partition")


def difference_matrix(cert: dict) -> list[list[int]]:
    """Generate the normalized difference matrix; its first column is zero."""
    if cert["kind"] != "difference_matrix":
        raise ValueError("This instance uses translated OA blocks, not an ordinary DM")
    s = math.prod(cert["coarse_radices"])
    E = affine_values(cert)
    subtract = group_table(cert["coarse_radices"] + [cert["p"]] * cert["d"], -1)
    result = []
    for t, coarse_row in enumerate(cert["H"]):
        for values in E[t]:
            raw = [h + s * v for h, v in zip(coarse_row, values)]
            result.append([subtract[value][raw[0]] for value in raw])
    if len(result) != cert["order"]:
        raise ValueError("Wrong difference-matrix size")
    return result


def squares(cert: dict) -> tuple[list[list[list[int]]], list[list[int]] | None]:
    """Expand a certificate into squares, and also return its DM when applicable."""
    n, k = cert["order"], cert["mols"]
    if cert["kind"] == "difference_matrix":
        D = difference_matrix(cert)
        add = group_table(cert["coarse_radices"] + [cert["p"]] * cert["d"])
        # The row number x is an index, not a group element used in arithmetic.
        result = [
            [[add[D[x][j + 1]][y] for y in range(n)] for x in range(n)]
            for j in range(k)
        ]
        return result, D

    # The other two instances generate an OA directly. Its first two outputs
    # identify a grid cell; all the other outputs are its square entries.
    s = math.prod(cert["coarse_radices"])
    q = cert["p"] ** cert["d"]
    E = affine_values(cert)
    add = group_table([cert["p"]] * cert["d"])
    grid = [None] * (n * n)
    for t, coarse_row in enumerate(cert["H"]):
        for values in E[t]:
            for z in range(q):
                row = [h + s * add[v][z] for h, v in zip(coarse_row, values)]
                position = n * row[0] + row[1]
                if grid[position] is not None:
                    raise ValueError("Repeated grid position")
                grid[position] = row[2:]
    if any(row is None for row in grid):
        raise ValueError("Unfilled grid cell")
    result = [
        [[grid[n * x + y][j] for y in range(n)] for x in range(n)]
        for j in range(k)
    ]
    return result, None
