"""Affine difference-matrix lifts for orders 45, 72, 80, and 96.

The compact seeds are expanded for the shared affine construction engine.
A vector is encoded by its base-p digits, least significant digit first.
Matrix M[i] is stored by its encoded columns, so it represents a map
F_p^a -> F_p^d.  The JSON files are exactly the data printed in the paper.
"""
import json
from pathlib import Path

CERTIFICATES = Path(__file__).resolve().parents[1] / 'certificates'
ORDERS = (45, 72, 80, 96)


def digits(value, p, length):
    """Decode an integer as a vector of the requested length."""
    return [(value // p**j) % p for j in range(length)]


def field16_multiply(x, y):
    """Multiply polynomial encodings modulo X^4 + X + 1."""
    result = 0
    while y:
        if y & 1:
            result ^= x
        y >>= 1
        x <<= 1
        if x & 16:
            x ^= 19  # X^4 + X + 1
    return result


def certificate(n):
    """Read the compact certificate and expand its short generation rules."""
    if n not in ORDERS:
        raise ValueError(f'Unsupported affine-lift order: {n}')
    data = json.loads((CERTIFICATES / f'dm{n}.json').read_text())
    if n == 72:
        def power_T(x, exponent):
            for _ in range(exponent % 7):
                x = ((x << 1) & 7) ^ (3 if x & 4 else 0)
            return x

        data['H'] = [[0, x] + [power_T(f[power_T(x, -j)], j)
                               for j in range(7)]
                     for f in data['permutations'] for x in range(8)]
        C = [[None] * 9 for _ in range(24)]
        C[:8] = [row[:] for row in data['seed_offsets']]
        S, tau, P = (data[key] for key in
                     ('row_permutation', 'layer_permutation', 'column_permutation'))
        for t in (0, 2):  # recurrence first fills layer 2, then layer 1
            for x in range(8):
                for i in range(9):
                    u, v = digits(C[8*t + x][i], 3, 2)
                    delta = i if i < 3 else 0
                    C[8*tau[t] + S[x]][P[i]] = (u + v) % 3 + 3*((v + delta) % 3)
        data['C'] = C
    elif n == 80:
        data['H'] = [[0] + [A[(x+j) % 5] for j in range(5)]
                      + [B[(x+j) % 5] for j in range(5)]
                     for A, B in zip(data['A'], data['B']) for x in range(5)]
        data['M'] = [[a, field16_multiply(a, 6)] for a in data['slopes']]
    return data
