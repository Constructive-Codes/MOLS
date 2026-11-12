"""Shared development algorithms. All construction arithmetic is exact.

Group elements use little-endian mixed-radix integer encodings. Quasi-
difference matrices use integers m,...,m+h-1 for the infinity symbols.
"""

from itertools import combinations
from math import prod


class Group:
    """An additive product of cyclic groups with mixed-radix labels."""

    def __init__(self, moduli=None):
        self.moduli = list(moduli or [])
        self.order = prod(self.moduli)
        self.table = [[self._add(a, b) for b in range(self.order)]
                      for a in range(self.order)]
        self.inverses = [self.table[a].index(0) for a in range(self.order)]

    def _add(self, a, b):
        value, place = 0, 1
        for radix in self.moduli:
            value += ((a % radix + b % radix) % radix) * place
            a //= radix
            b //= radix
            place *= radix
        return value

    def add(self, a, b):
        return self.table[a][b]

    def difference(self, b, a):
        """Additive difference b-a, invariant under common translation."""
        return self.table[b][self.inverses[a]]


def group_from_spec(spec):
    return Group(spec["moduli"])


def check_difference_matrix(matrix, group, hole=0):
    """Exhaustively check the precise DM/QDM conditions used in the paper."""
    m = group.order
    if len(matrix) != m + 2 * hole or not matrix:
        raise ValueError("Wrong number of difference-matrix rows")
    width = len(matrix[0])
    if width < 2 or any(len(row) != width for row in matrix):
        raise ValueError("Ragged or too narrow difference matrix")
    if any(type(v) is not int or not 0 <= v < m + hole
           for row in matrix for v in row):
        raise ValueError("Difference-matrix entry outside its alphabet")
    if any(sum(v >= m for v in row) > 1 for row in matrix):
        raise ValueError("Two infinity entries in one row")
    for infinity in range(m, m + hole):
        for j in range(width):
            if sum(row[j] == infinity for row in matrix) != 1:
                raise ValueError("Infinity must occur once in each column")
    for i, j in combinations(range(width), 2):
        differences = [group.difference(row[j], row[i]) for row in matrix
                       if row[i] < m and row[j] < m]
        if len(differences) != m or set(differences) != set(range(m)):
            raise ValueError(f"Difference condition fails for columns {i},{j}")
    return {"matrix_rows": len(matrix), "matrix_columns": width,
            "difference_pairs": width * (width - 1) // 2,
            "group_order": m, "hole": hole}


def dm_to_squares(matrix, group):
    """Normalize column zero, then L_j(r,g)=Dbar[r,j]+g."""
    normalized = [[group.difference(v, row[0]) for v in row] for row in matrix]
    n = group.order
    squares = [[[group.add(row[j], g) for g in range(n)]
                for row in normalized] for j in range(1, len(matrix[0]))]
    return squares, normalized


def field_product(a, b, q):
    """Multiplication needed by the small hole arrays, q=7, q=8, or q=9."""
    if q == 7:
        return a * b % 7
    if q == 8:
        # F_8=F_2[z]/(z^3+z+1), with binary coefficient labels.
        value = 0
        while b:
            if b & 1:
                value ^= a
            b >>= 1
            a <<= 1
            if a & 8:
                a ^= 11
        return value
    if q == 9:
        # F_9=F_3[z]/(z^2+1), enc(a0+a1*z)=a0+3*a1.
        a0, a1, b0, b1 = a % 3, a // 3, b % 3, b // 3
        return (a0*b0 - a1*b1) % 3 + 3*((a0*b1 + a1*b0) % 3)
    raise ValueError("Unsupported hole field")


def hole_array(hole, width):
    """Field arrays for the one-, seven-, eight-, and nine-symbol holes."""
    if hole == 1:
        return [[0] * width]
    if hole == 8 and width == 8:
        return [[y ^ field_product(m, x, 8) for m in range(8)]
                for x in range(8) for y in range(8)]
    if hole not in (7, 9) or width != hole + 1:
        raise ValueError("Unsupported hole size or number of coordinates")
    additive = Group([7] if hole == 7 else [3, 3])
    return [[u, v] + [additive.add(u, field_product(a, v, hole))
                       for a in range(1, width - 1)]
            for u in range(hole) for v in range(hole)]


def qdm_to_array(matrix, group, hole):
    """Add every group translation to each seed row and fill the hole."""
    m, width = group.order, len(matrix[0])
    rows = [[group.add(v, g) if v < m else v for v in row]
            for row in matrix for g in range(m)]
    rows.extend([[m + v for v in row] for row in hole_array(hole, width)])
    return rows


def check_array(rows, n):
    if len(rows) != n*n or not rows:
        raise ValueError("An index-one OA must have n^2 rows")
    width = len(rows[0])
    if width < 3 or any(len(row) != width for row in rows):
        raise ValueError("Ragged or too narrow orthogonal array")
    if any(type(v) is not int or not 0 <= v < n for row in rows for v in row):
        raise ValueError("OA symbol out of range")
    for i, j in combinations(range(width), 2):
        if len({n*row[i]+row[j] for row in rows}) != n*n:
            raise ValueError(f"OA condition fails for coordinates {i},{j}")
    return {"oa_rows": len(rows), "oa_columns": width,
            "oa_pairs": width*(width-1)//2}


def array_to_squares(rows, n):
    squares = [[[-1] * n for _ in range(n)] for _ in range(len(rows[0])-2)]
    occupied = set()
    for row in rows:
        r, c = row[:2]
        if (r, c) in occupied:
            raise ValueError("Repeated grid cell in OA")
        occupied.add((r, c))
        for j, square in enumerate(squares):
            square[r][c] = row[j+2]
    if len(occupied) != n*n:
        raise ValueError("Missing grid cell in OA")
    return squares
