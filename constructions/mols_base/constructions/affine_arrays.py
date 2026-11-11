"""Affine lifts of binary arrays, giving seven MOLS(50), nine MOLS(54),
and seven MOLS(98).

Only the Python standard library is needed.  The finite data are in
certificates/oa50.json, oa54.json and oa98.json; all indices start at zero.

For each binary row t and x in F_p^d, output the coordinate
    H[t][i] + 2 * encode(M[i][H[t][i]] x + C[t][i]).
Here encode(v) = sum(v[j] * p**j).  At order 98, the matrices are built
from their two slopes W, so that M(u,x,y) = (x,y) + u W.
"""

from itertools import product
import json
from pathlib import Path


CERTIFICATES = Path(__file__).resolve().parents[1] / "certificates"
PARAMETERS = {50: (5, 3, 2, 9), 54: (3, 5, 3, 11), 98: (7, 3, 2, 9)}


def require(condition, message):
    """Validation must remain active even when Python is run with -O."""
    if not condition:
        raise ValueError(message)


def load(n):
    """Return p, input dimension, output dimension, H, M, C."""
    if n not in PARAMETERS:
        raise ValueError("This family supports only orders 50, 54 and 98")
    p, d, r, columns = PARAMETERS[n]
    data = json.loads((CERTIFICATES / f"oa{n}.json").read_text())
    require(data["n"] == n, "Wrong certificate order")
    if n == 50:
        H, matrices, C = expand50(data)
        return p, d, r, H, matrices, C
    if n == 54:
        matrices = data["M"]
    else:
        slopes = data["W"]
        require(len(slopes) == columns, "Wrong number of slopes")
        require(all(len(pair) == 2 and all(len(w) == 2 for w in pair)
                    for pair in slopes), "Wrong slope shape")
        matrices = [[[[w[0], 1, 0], [w[1], 0, 1]] for w in pair]
                    for pair in slopes]
    return p, d, r, data["H"], matrices, data["C"]


def rotate_triples(values):
    """Rotate each consecutive triple (a,b,c) to (c,a,b)."""
    require(len(values) % 3 == 0, "Rotation needs complete triples")
    return [values[b + j] for b in range(0, len(values), 3) for j in (2, 0, 1)]


def expand50(data):
    """Expand the order-50 numerical seeds.

    The binary symbol does not affect the fine matrix in this case, so
    copy each matrix into both slots of the shared family representation.
    """
    require(len(data["matrix_seeds"]) == 3, "Wrong matrix-seed count")
    matrices = []
    for seed in data["matrix_seeds"]:
        require(len(seed) == 2 and all(len(row) == 3 for row in seed),
                "Wrong matrix-seed shape")
        matrix = [list(row) for row in seed]
        for _ in range(3):
            matrices.append([[list(row) for row in matrix] for _ in range(2)])
            matrix = [rotate_triples(row) for row in matrix]
    H, C = [], []
    for seed in data["row_seeds"]:
        binary, packed = list(map(int, seed["H"])), list(seed["D"])
        require(seed["length"] in (1, 3), "Wrong coarse-orbit length")
        require(len(binary) == len(packed) == 9, "Wrong row-seed length")
        require(all(isinstance(v, int) and 0 <= v < 25 for v in packed),
                "Invalid packed offset")
        if seed["length"] == 1:
            require(rotate_triples(binary) == binary and
                    rotate_triples(packed) == packed,
                    "Length-one row is not fixed by rotation")
        for _ in range(seed["length"]):
            H.append(binary)
            C.append([[v % 5, v // 5] for v in packed])
            binary, packed = rotate_triples(binary), rotate_triples(packed)
    return H, matrices, C


def build(n):
    """Expand the certificate into an OA; coordinates 0 and 1 index the grid."""
    p, d, r, H, M, C = load(n)
    rows = []
    for t, binary_row in enumerate(H):
        for x in product(range(p), repeat=d):
            row = []
            for i, a in enumerate(binary_row):
                fine = [(sum(coefficient * value for coefficient, value in zip(equation, x))
                         + offset) % p
                        for equation, offset in zip(M[i][a], C[t][i])]
                encoded = sum(value * p ** j for j, value in enumerate(fine))
                row.append(a + 2 * encoded)
            rows.append(row)
    return {"kind": "oa", "n": n, "rows": rows}
