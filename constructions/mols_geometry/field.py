"""Small finite fields with the exact encodings used in the manuscript.

No external computer algebra system is required.  An integer label is a
base-p coefficient vector, NOT an integer residue modulo p**d.
"""
from __future__ import annotations
from functools import lru_cache
from math import isqrt

# X**d + sum(low[i] * X**i).  Coefficients are listed constant-first.
EXTENSIONS = {
    8: (2, (1, 1, 0)),          # X^3 + X + 1
    9: (3, (1, 0)),             # X^2 + 1
    16: (2, (1, 1, 0, 0)),      # X^4 + X + 1
    25: (5, (2, 0)),            # X^2 + 2
    27: (3, (2, 2, 0)),         # X^3 + 2X + 2
    32: (2, (1, 0, 1, 0, 0)),   # X^5 + X^2 + 1
}


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


class Field:
    """Arithmetic tables for one of the small fields used by the recipes."""
    def __init__(self, order: int):
        self.order = order
        if order in EXTENSIONS:
            self.p, self.low = EXTENSIONS[order]
        elif is_prime(order):
            self.p, self.low = order, (0,)  # prime field; degree one
        else:
            raise ValueError(f"No specified field representation for {order}")
        self.degree = len(self.low)
        p, d = self.p, self.degree
        digits = [tuple((x // p**i) % p for i in range(d))
                  for x in range(order)]

        def encode(coefficients):
            return sum((a % p) * p**i for i, a in enumerate(coefficients))

        def multiply(a, b):
            coefficients = [0] * (2*d - 1)
            for i, x in enumerate(a):
                for j, y in enumerate(b):
                    coefficients[i+j] += x*y
            # Polynomial long division by the specified monic modulus.
            for power in range(2*d - 2, d - 1, -1):
                lead = coefficients[power] % p
                for j, coefficient in enumerate(self.low):
                    coefficients[power-d+j] -= lead * coefficient
            return encode(coefficients[:d])

        self.add_table = [[encode(a+b for a, b in zip(x, y)) for y in digits]
                          for x in digits]
        self.neg_table = [encode(-a for a in x) for x in digits]
        self.mul_table = [[multiply(x, y) for y in digits] for x in digits]
        self.inverse_table = [0]
        for x in range(1, order):
            if sorted(self.mul_table[x]) != list(range(order)):
                raise ValueError("The supplied polynomial does not define a field")
            self.inverse_table.append(self.mul_table[x].index(1))
        self.nonzero_squares = frozenset(self.mul(x, x) for x in range(1, order))

    def add(self, x: int, y: int) -> int:
        return self.add_table[x][y]

    def neg(self, x: int) -> int:
        return self.neg_table[x]

    def sub(self, x: int, y: int) -> int:
        return self.add(x, self.neg(y))

    def mul(self, x: int, y: int) -> int:
        return self.mul_table[x][y]

    def inv(self, x: int) -> int:
        if x == 0:
            raise ZeroDivisionError("Zero has no multiplicative inverse")
        return self.inverse_table[x]

    def div(self, x: int, y: int) -> int:
        return self.mul(x, self.inv(y))

    def pow(self, x: int, exponent: int) -> int:
        if exponent < 0:
            return self.pow(self.inv(x), -exponent)
        result = 1
        while exponent:
            if exponent & 1:
                result = self.mul(result, x)
            x = self.mul(x, x)
            exponent //= 2
        return result

    def character(self, x: int) -> int:
        """Quadratic character in an ODD field, with character(0) = 0."""
        if self.p == 2:
            raise ValueError("Quadratic-character recipes require odd characteristic")
        return 0 if x == 0 else (1 if x in self.nonzero_squares else -1)


@lru_cache(maxsize=None)
def field(order: int) -> Field:
    return Field(order)
