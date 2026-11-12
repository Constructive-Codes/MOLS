"""Nine MOLS of order 225 from a difference matrix over Z_15 x Z_15.

The two tables and ten coordinate formulas are the fixed data in the paper.
Group elements (u, v) have label u + 15*v; row (x, y) has index x + 15*y.
The shared development normalizes column zero and constructs the squares.
"""

B = (0, 4, 6, 8, 13, 3, 5, 14, 1, 10, 12, 2, 7, 9, 11)
C = (0, 12, 10, 13, 7, 6, 4, 1, 14, 11, 9, 8, 2, 5, 3)


def build():
    """Return the 225-by-10 difference matrix in increasing row-label order."""
    matrix = []
    for y in range(15):
        for x in range(15):
            entries = (
                (B[x] + B[y], x + y),
                (C[x] + y, x + B[y]),
                (B[y], B[x] + B[y]),
                (y, C[x] + y),
                (x + C[y], x + C[y]),
                (C[y], x),
                (B[x] + y, B[x]),
                (B[x] + C[y], C[x] + B[y]),
                (C[x] + C[y], B[x] + y),
                (C[x] + B[y], C[x]),
            )
            matrix.append([u % 15 + 15 * (v % 15) for u, v in entries])
    return {"kind": "dm", "moduli": [15, 15], "matrix": matrix}
