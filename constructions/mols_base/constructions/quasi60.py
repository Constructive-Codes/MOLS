"""Six MOLS of order 60 from a 68 by 8 quasi-difference matrix.

The group is Z_26 x Z_2, with label t+26*e for (t,e).
Labels 52,...,59 are the eight infinity symbols.
"""


TAILS = (
    (40, 45, 26, 28, 38, 39),
    (23, 39, 29, 5, 26, 32),
    (13, 5, 30, 13, 22, 41),
    (33, 43, 22, 14, 33, 8),
    (47, 50, 21, 25, 1, 3),
    (16, 22, 11, 33, 37, 49),
    (32, 34, 31, 37, 2, 20),
    (24, 10, 20, 4, 21, 25),
)
FULL = (0, 0, 14, 14, 16, 16, 9, 9)
PERMUTATIONS = (
    (1, 0, 3, 2, 5, 4, 7, 6),
    (2, 3, 0, 1, 6, 7, 4, 5),
    (4, 5, 6, 7, 0, 1, 2, 3),
)
OFFSETS = (
    (0, 0, 0, 0, 0, 0, 0, 0),
    (26, 39, 0, 13, 26, 39, 0, 13),
    (39, 13, 39, 13, 0, 26, 0, 26),
)


def add(a, b):
    """Add in the product group, fixing any infinity symbol in a."""
    if a >= 52:
        return a
    return (a % 26 + b % 26) % 26 + 26*((a // 26 + b // 26) % 2)


def transform(row, mask):
    """Apply operations 2, 1, then 0, selected by bits of mask."""
    for j in (2, 1, 0):
        if (mask >> j) & 1:
            row = tuple(add(row[p], c)
                        for p, c in zip(PERMUTATIONS[j], OFFSETS[j]))
    return list(row)


def build():
    """Construct the QDM; shared code develops it and fills the eight-symbol hole."""
    matrix = [transform((52+h, 0, *tail), mask)
              for h, tail in enumerate(TAILS) for mask in range(8)]
    matrix.extend(transform(FULL, mask) for mask in (0, 2, 4, 6))
    return {"kind": "qdm", "moduli": [26, 2], "hole": 8,
            "matrix": matrix}
