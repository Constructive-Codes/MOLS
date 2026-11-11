"""Eight MOLS of order 90 from a 99 by 10 quasi-difference matrix.

F = F_3[u]/(u^2+1), K = F[v]/(v^2-(1+u)).  Encode
a0+a1*u+b0*v+b1*u*v as a0+3*a1+9*b0+27*b1.  All code in this
module uses the Python standard library.
"""


BETA = 31
ALPHA = 62
PLANE = (0, 1, 2, 9, 10, 11, 18, 19, 20)
INTERCEPTS = (0, 54, 59, 47, 51, 0, 69, 13, 4, 40)


def add(a, b):
    """Add four base-three coordinates; labels are not integers mod 81."""
    return sum(((a // 3**i + b // 3**i) % 3) * 3**i for i in range(4))


def multiply9(a, b):
    """Multiply a0+a1*u and b0+b1*u, using u^2=-1."""
    a0, a1, b0, b1 = a % 3, a // 3, b % 3, b // 3
    return (a0*b0-a1*b1) % 3 + 3*((a0*b1+a1*b0) % 3)


def multiply(a, b):
    """Multiply (a0+a1*v)(b0+b1*v), using v^2=1+u (label 4)."""
    a0, a1, b0, b1 = a % 9, a // 9, b % 9, b // 9
    constant = add(multiply9(a0, b0), multiply9(4, multiply9(a1, b1)))
    linear = add(multiply9(a0, b1), multiply9(a1, b0))
    return constant + 9*linear


def power(a, exponent):
    value = 1
    for _ in range(exponent):
        value = multiply(value, a)
    return value


def build():
    """Construct the QDM; shared code develops it and fills the nine-symbol hole."""
    intercepts = INTERCEPTS
    if len(intercepts) != 10 or any(type(z) is not int or not 0 <= z < 81
                                    for z in intercepts):
        raise ValueError("Expected ten field labels in 0,...,80")
    beta_powers = [power(BETA, c) for c in range(10)]
    alpha_powers = [power(ALPHA, t) for t in range(10)]
    matrix = [[multiply(w, beta_powers[c]) for c in range(10)] for w in PLANE]
    for h, w in enumerate(PLANE):
        for t in range(10):
            row = []
            for c in range(10):
                if c == t:
                    row.append(81+h)
                else:
                    seed = add(multiply(w, beta_powers[(c-t) % 10]),
                               intercepts[(c-t) % 10])
                    row.append(multiply(alpha_powers[t], seed))
            matrix.append(row)
    return {"kind": "qdm", "moduli": [3, 3, 3, 3], "hole": 9,
            "matrix": matrix}
