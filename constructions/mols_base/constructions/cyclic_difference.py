"""Seven MOLS of order 85 from three cycles over the additive group Z/85Z.

The certificate is the 21 entries printed in the paper.  Expansion is purely
deterministic; common.dm_to_squares performs the shared DM development.
Only the Python standard library is used.
"""

import json
from pathlib import Path


CERTIFICATE = Path(__file__).resolve().parents[1] / "certificates" / "dm85.json"


def read_certificate():
    """Read the three cycles and check their dimensions."""
    data = json.loads(CERTIFICATE.read_text())
    if data.get("n") != 85 or data.get("beta") != 47:
        raise ValueError("Expected order 85 and multiplier 47")
    cycles = data.get("seedcycles")
    if not isinstance(cycles, list) or len(cycles) != 3:
        raise ValueError("Expected three seed cycles")
    if any(not isinstance(cycle, list) or len(cycle) != 7 or
           any(type(x) is not int or not 0 <= x < 85 for x in cycle)
           for cycle in cycles):
        raise ValueError("Each seed cycle must contain seven residues modulo 85")
    return data


def build():
    """Return the 85-by-8 difference matrix, with an all-zero first column.

    Row zero is all zero.  The remaining rows are ordered by seed cycle,
    then multiplier exponent 0,...,3, then cyclic shift 0,...,6.  This is
    exactly the ordering of the supplied order-85 squares.
    """
    data = read_certificate()
    n, beta = data["n"], data["beta"]
    matrix = [[0] * 8]
    for cycle in data["seedcycles"]:
        for phase in range(4):
            multiplier = pow(beta, phase, n)
            for shift in range(7):
                matrix.append([0] + [
                    multiplier * cycle[(j + shift) % 7] % n
                    for j in range(7)
                ])
    return {"kind": "dm", "moduli": [n], "matrix": matrix}
