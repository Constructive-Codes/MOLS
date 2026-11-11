"""Projective incidence structures as ordinary lists of point IDs.

Affine (x,y) has ID q*x+y.  Infinity point for slope m has ID q*q+m.
The vertical infinity point has ID q*q+q.  This fixes reproducible labels.
"""
from __future__ import annotations
from itertools import combinations
from .field import field


def direction_value(q: int, kind: str, m: int, x: int) -> int:
    """Evaluate one of the three direction families used by the constructions."""
    f = field(q)
    if kind == "ordinary":
        return f.mul(m, x)
    if kind == "switched27" and q == 27:
        exponent = 1 if m == 0 or m in f.nonzero_squares else 3
        return f.mul(m, f.pow(x, exponent))
    if kind == "twisted27" and q == 27:
        return f.add(f.mul(m, x), f.mul(f.pow(m, 3), f.pow(x, 9)))
    raise ValueError(f"Unknown plane type {kind!r} at order {q}")


def plane_lines(q: int, kind: str = "ordinary") -> list[list[int]]:
    """Enumerate all q^2+q+1 lines, each of size q+1."""
    f = field(q)
    lines = []
    for m in range(q):
        values = [direction_value(q, kind, m, x) for x in range(q)]
        for intercept in range(q):
            line = [q*x + f.add(values[x], intercept) for x in range(q)]
            lines.append(line + [q*q + m])
    for x in range(q):
        lines.append([q*x+y for y in range(q)] + [q*q+q])
    lines.append(list(range(q*q, q*q+q+1)))
    return lines


def check_plane(q: int, lines: list[list[int]]) -> None:
    """Exhaustive pair-incidence check; inexpensive for q <= 31."""
    size = q*q+q+1
    if len(lines) != size:
        raise ValueError("Wrong number of lines")
    seen = bytearray(size*size)
    for line in lines:
        if len(line) != q+1 or len(set(line)) != q+1:
            raise ValueError("Malformed line")
        if not all(0 <= x < size for x in line):
            raise ValueError("Point outside ambient set")
        for x, y in combinations(sorted(line), 2):
            index = size*x+y
            if seen[index]:
                raise ValueError("Ambient pair covered more than once")
            seen[index] = 1
    if sum(seen) != size*(size-1)//2:
        raise ValueError("An ambient pair is missing")
