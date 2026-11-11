"""Explicit point selections for all twelve geometric parameter gains.

Constants and integer field labels use the encodings in the manuscript.
"""
from __future__ import annotations
from dataclasses import dataclass
from .field import field


@dataclass(frozen=True)
class Recipe:
    order: int
    squares: int
    plane_order: int
    plane_kind: str
    selected: frozenset[int]
    histogram: dict[int, int]     # includes empty/singleton ambient intersections
    exceptional_size: int | None = None
    center_ambient_id: int | None = None


ORDERS = (175, 183, 207, 328, 335, 427, 430, 442, 456, 469, 482, 494)


def ray_points(q, positive, negative, full, infinity, include_origin=False):
    """Select quadratic half-rays, full rays, and designated infinity points."""
    f = field(q)
    positive, negative, full = set(positive), set(negative), set(full)
    if positive & negative or full & (positive | negative):
        raise ValueError("Ray types overlap")
    points = {0} if include_origin else set()
    for direction in range(q+1):
        sign = 1 if direction in positive else -1 if direction in negative else 0
        for r in range(1, q):
            if direction not in full and (not sign or f.character(r) != sign):
                continue
            x, y = (r, f.mul(r, direction)) if direction < q else (0, r)
            points.add(q*x+y)
    points.update(q*q+d for d in infinity)
    return frozenset(points)


def norm25(z):
    """Norm(a+b*alpha)=a^2+2*b^2 for alpha^2=3 over F5."""
    a, b = z % 5, z // 5
    return (a*a+2*b*b) % 5


def arc_points(coefficient, subgroup, translate_y=0):
    f = field(16)
    points = set()
    for x in range(16):
        for y in range(16):
            z = f.add(y, translate_y)
            value = f.add(f.add(f.mul(x, x), f.mul(x, z)),
                          f.mul(coefficient, f.mul(z, z)))
            if value in subgroup:
                points.add(16*x+y)
    return points


def make_recipe(n: int) -> Recipe:
    if n == 175:
        f = field(16)
        arc_a = arc_points(8, {0, 1, 2, 3})
        arc_b = set()
        # T[X:Y:Z]=[Z:alpha*X:Y], with alpha labelled 2.
        for point in arc_a:
            x, y = divmod(point, 16)
            if y:
                arc_b.add(16*f.inv(y)+f.div(f.mul(2, x), y))
            else:
                arc_b.add(256+f.mul(2, x))
        points = set(range(273)) - (arc_a ^ arc_b)
        return Recipe(n, 7, 16, "ordinary", frozenset(points),
                      {9: 124, 11: 45, 13: 101, 17: 3})

    if n == 207:
        arc_a = arc_points(8, {0, 1, 2, 3})
        arc_b = arc_points(13, {0, 6, 10, 12}, translate_y=1)
        points = set(range(273)) - (arc_a ^ arc_b)
        return Recipe(n, 9, 16, "ordinary", frozenset(points),
                      {11: 77, 13: 165, 17: 31})

    if n == 183:
        q, f = 16, field(16)
        allowed_products = set(range(10))  # 0 and labels 1,...,9
        points = {q*x+y for x in range(q) for y in range(q)
                  if f.mul(x, y) in allowed_products}
        points.update(range(q*q, q*q+q+1))
        return Recipe(n, 7, q, "ordinary", frozenset(points),
                      {9: 45, 11: 165, 13: 45, 17: 18})

    if n == 328:
        half = {0, 1, 2, 3, 4, 25}  # subfield directions, plus vertical
        full = set(range(5, 15))
        points = ray_points(25, [], half, full, half | full)
        return Recipe(n, 9, 25, "ordinary", points,
                      {0: 10, 11: 72, 13: 486, 16: 73, 25: 10})

    if n == 335:
        circle = {z for z in range(25) if norm25(z) == 2}
        signs = [1 if z // 5 <= 2 else -1 for z in range(25)] + [1]
        positive = {z for z in range(26) if z not in circle and signs[z] == 1}
        negative = {z for z in range(26) if z not in circle and signs[z] == -1}
        full = {z for z in circle if signs[z] == 1}
        points = ray_points(25, positive, negative, full, positive | negative | full)
        return Recipe(n, 9, 25, "ordinary", points,
                      {0: 3, 11: 240, 13: 164, 16: 240, 23: 1, 25: 3})

    if n == 442:
        positive = {5, 14, 16, 21}
        negative = {6, 17, 22, 24, 27}
        full = {0, 2, 3, 4, 7, 8, 9, 10, 11, 12, 13, 15}
        points = ray_points(27, positive, negative, full, full | {1})
        return Recipe(n, 11, 27, "ordinary", points,
                      {0: 6, 1: 1, 13: 88, 16: 325, 17: 234, 19: 91, 27: 12})

    if n == 494:
        positive = {2, 4, 9, 12, 17, 23, 25, 26, 28}
        negative = {0, 1, 3}
        full = {5, 10, 11, 13, 14, 15, 19, 21, 22, 31}
        infinity = positive | negative | {31}
        points = ray_points(31, positive, negative, full, infinity, True)
        return Recipe(n, 11, 31, "ordinary", points,
                      {1: 10, 13: 226, 16: 330, 17: 327, 19: 90, 31: 9, 32: 1})

    if n in (430, 469):
        q, f = 27, field(27)
        if n == 430:
            allowed = {2,3,4,5,6,9,10,11,13,14,15,16,18,19,20,21}
        else:
            allowed = {1,5,6,7,8,9,10,11,13,14,15,16,18,19,21,22,0}
        points = {q*x+y for x in range(q) for y in range(q)
                  if f.mul(x, y) in allowed}
        if n == 430:
            points.add(0)
            points.update(q*q+m for m in range(1, q) if f.character(m) == -1)
            histogram = {1: 2, 13: 183, 16: 455, 19: 91, 25: 26}
        else:
            histogram = {0: 1, 13: 78, 16: 260, 17: 78, 19: 312, 25: 26, 27: 2}
        return Recipe(n, 11, q, "ordinary", frozenset(points), histogram)

    if n == 427:
        q, f = 25, field(25)
        excluded = {z for z in range(q) if norm25(z) == 3} | {5, 20}
        points = {q*x+f.mul(c, f.pow(x, 4)) for x in range(1, q)
                  for c in range(q) if c not in excluded}
        points.update(y for y in range(q) if norm25(y) != 2)
        return Recipe(n, 12, q, "ordinary", frozenset(points),
                      {0: 1, 1: 6, 13: 8, 16: 246, 17: 216, 19: 161, 25: 13},
                      exceptional_size=13, center_ambient_id=0)

    if n == 456:
        q, f = 27, field(27)
        def polynomial(c):
            value = f.pow(c, 5)
            for term in (f.pow(c, 3), f.mul(c, c), c):
                value = f.sub(value, term)
            return f.add(value, 1)
        coefficients = {c for c in range(1, q) if f.character(polynomial(c)) == 1}
        points = {q*x+f.mul(c, f.pow(x, 23)) for x in range(1, q)
                  for c in coefficients}
        points.update(range(q))  # complete vertical axis
        points.update(q*q+m for m in f.nonzero_squares)
        return Recipe(n, 11, q, "twisted27", frozenset(points),
                      {1: 1, 13: 157, 16: 117, 17: 208, 19: 247, 23: 26, 27: 1})

    if n == 482:
        q, f = 27, field(27)
        beta = 3
        special = {f.neg(1), f.neg(f.mul(beta, beta)),
                   f.neg(f.pow(f.add(1, beta), 2))}
        coefficients = set(f.nonzero_squares) | special
        points = {q*x for x in range(q)} | set(range(q))  # both axes
        points.update(q*x+f.mul(c, f.pow(x, 9)) for x in range(1, q)
                      for c in coefficients)
        points.update(q*q+m for m in range(1, q) if f.character(m) == -1)
        return Recipe(n, 7, q, "switched27", frozenset(points),
                      {8: 13, 13: 1, 16: 78, 17: 390, 19: 234, 25: 26, 27: 15},
                      exceptional_size=8, center_ambient_id=0)

    raise ValueError(f"Supported geometric orders are {ORDERS}")
