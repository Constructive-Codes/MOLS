"""Seven MOLS of orders 78 and 92 by one relative-difference algorithm.

Construction needs only Python's standard library and the compact JSON
certificates.  There are fourteen seed rows over F_7 x F_p, p=11 or 13.
Multiplication of their fine coordinates by the nonzero squares of F_p
gives a relative difference matrix.  Translate it, then fill every
seven-symbol hole together with one common infinity using the field F_8.
"""

import json
from pathlib import Path


CERTIFICATES = Path(__file__).resolve().parents[1] / "certificates"
ORDERS = (78, 92)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def certificate(n):
    """Read the data; expand the two quadratic seed sequences at order 78."""
    require(n in ORDERS, "This family supports only orders 78 and 92")
    data = json.loads((CERTIFICATES / f"rdm{n}.json").read_text())
    require(data["n"] == n and n == 7 * data["p"] + 1,
            "Wrong relative-difference certificate parameters")
    if n == 78:
        H, Y = [], []
        for c, b, U in zip(data["c"], data["b"], data["U"]):
            for a in range(7):
                H.append([0, (a - pow(4*c, -1, 7)) % 7]
                         + [c * (a+t)**2 % 7 for t in range(7)])
                Y.append([0, b] + [U[(a+t) % 7] for t in range(7)])
        data["H"], data["Y"] = H, Y
    return data


def relative_matrix(H, Y, p):
    """Return rows in the encoding (a,b) -> a+7*b."""
    squares = sorted({x*x % p for x in range(1, p)})
    return [[h + 7 * (m*y % p) for h, y in zip(hrow, yrow)]
            for hrow, yrow in zip(H, Y) for m in squares]


def multiply8(x, y):
    """Polynomial bit encodings in F_2[z]/(z^3+z+1)."""
    answer = 0
    while y:
        if y & 1:
            answer ^= x
        y >>= 1
        x <<= 1
        if x & 8:
            x ^= 11  # binary 1011 is z^3+z+1
    return answer


def develop_and_fill(matrix, p):
    """Translate the relative matrix and fill all p holes with one infinity.

    The nine F_8 forms, in order, are x, y, x+y, x+2*y, ..., x+7*y,
    where multiplication uses the field, not integer arithmetic.
    The finite label (a,b) is a+7*b and infinity is 7*p.
    """
    infinity = 7*p
    rows = [[(value % 7 + a) % 7 + 7*((value // 7 + b) % p)
             for value in row]
            for row in matrix for a in range(7) for b in range(p)]
    for b in range(p):
        for x in range(8):
            for y in range(8):
                if x == y == 0:
                    continue
                labels = [x, y] + [x ^ multiply8(a, y) for a in range(1, 8)]
                rows.append([infinity if e == 0 else e-1 + 7*b for e in labels])
    rows.append([infinity] * 9)
    return rows


def build_from_seeds(H, Y, p):
    """Construct an OA from the fourteen coarse and fine seed rows."""
    matrix = relative_matrix(H, Y, p)
    return {"kind": "oa", "n": 7*p+1,
            "rows": develop_and_fill(matrix, p)}


def build(n):
    data = certificate(n)
    return build_from_seeds(data["H"], data["Y"], data["p"])
