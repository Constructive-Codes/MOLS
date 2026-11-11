"""Orbit constructions over Z_21 and Z_67 for orders 22 and 74."""

CYCLIC22_BASES = [
    [21, 0, 0, 0, 1, 4], [0, 1, 4, 21, 14, 14],
    [0, 6, 16, 3, 11, 8], [0, 4, 6, 16, 12, 3],
    [0, 7, 11, 14, 3, 17], [0, 10, 3, 0, 7, 11],
]
CYCLIC22_FIXED = [
    [0, 11, 19, 7, 18, 5], [0, 18, 12, 6, 9, 15],
    [0, 3, 9, 8, 5, 20], [0, 2, 13, 17, 1, 18],
    [0, 5, 1, 18, 20, 10],
]
QUATERNION74_BASES = [
    [67, 1, 29, 37, 0, 35, 10, 22],
    [68, 29, 37, 1, 0, 50, 43, 41],
    [69, 39, 59, 36, 0, 19, 15, 33],
    [70, 59, 36, 39, 0, 52, 34, 48],
    [71, 10, 22, 35, 0, 40, 21, 6],
    [72, 43, 41, 50, 0, 7, 2, 58],
    [73, 25, 55, 54, 0, 3, 20, 44],
    [0, 33, 37, 60, 5, 21, 23, 43],
    [0, 65, 19, 1, 11, 41, 6, 64],
    [0, 29, 9, 15, 51, 47, 50, 40],
]


def cyclic_transform(row):
    """T(a,b,c,d,e,f)=16(c,a,b,f,d,e), fixing infinity=21."""
    return [21 if row[j] == 21 else 16*row[j] % 21
            for j in [2, 0, 1, 5, 3, 4]]


def quaternion_product(a, b):
    """Labels 0,...,7 mean 1,i,j,k,-1,-i,-j,-k."""
    positive_products = [[0, 1, 2, 3], [1, 4, 3, 6],
                         [2, 7, 4, 1], [3, 2, 5, 4]]
    return positive_products[a % 4][b % 4] ^ (4*((a//4) ^ (b//4)))


def build(order):
    if order == 22:
        rows = []
        for base in CYCLIC22_BASES:
            row = base[:]
            for _ in range(3):
                rows.append(row)
                row = cyclic_transform(row)
        rows.extend(row[:] for row in CYCLIC22_FIXED)
        return {"kind": "qdm", "moduli": [21], "hole": 1, "matrix": rows}
    if order == 74:
        rows = [[base[quaternion_product(g, u)] for u in range(8)]
                for base in QUATERNION74_BASES for g in range(8)]
        rows.append([0]*8)
        return {"kind": "qdm", "moduli": [67], "hole": 7, "matrix": rows}
    raise ValueError(f"Unknown small quasi construction: {order}")
