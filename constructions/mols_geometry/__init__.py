"""Geometric constructions of MOLS from selected finite-plane points.

construct(n) returns a list of ordinary n-by-n integer matrices.
"""
from .assembly import construct as construct_flat


def construct(order):
    flat = construct_flat(order)
    return [[list(square[row*order:(row+1)*order]) for row in range(order)]
            for square in flat]
