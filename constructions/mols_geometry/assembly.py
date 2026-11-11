"""Turn a selected incidence structure into complete mutually orthogonal squares."""
from __future__ import annotations
from array import array
from collections import Counter
from itertools import combinations
from .field import field
from .geometry import plane_lines, check_plane
from .recipes import make_recipe


def build_design(order: int) -> dict:
    recipe = make_recipe(order)
    points = sorted(recipe.selected)
    if len(points) != order:
        raise ValueError("Recipe has the wrong cardinality")
    lines = plane_lines(recipe.plane_order, recipe.plane_kind)
    check_plane(recipe.plane_order, lines)
    index = {point: i for i, point in enumerate(points)}
    intersections = [sorted(index[x] for x in line if x in index) for line in lines]
    histogram = dict(sorted(Counter(map(len, intersections)).items()))
    if histogram != recipe.histogram:
        raise ValueError(f"Unexpected line histogram: {histogram}")
    blocks = [block for block in intersections if len(block) > 1]
    seen = bytearray(order*order)
    for block in blocks:
        for x, y in combinations(block, 2):
            position = order*x+y
            if seen[position]:
                raise ValueError("Selected pair covered more than once")
            seen[position] = 1
    if sum(seen) != order*(order-1)//2:
        raise ValueError("A selected pair is missing")

    exceptional = []
    center = None
    if recipe.exceptional_size is not None:
        center = index[recipe.center_ambient_id]
        exceptional = [i for i, block in enumerate(blocks)
                       if len(block) == recipe.exceptional_size]
        for i in exceptional:
            if center not in blocks[i]:
                raise ValueError("An exceptional block misses its common center")
            # Give the common center local field label zero.
            blocks[i] = [center] + [x for x in blocks[i] if x != center]
        for i, j in combinations(exceptional, 2):
            if set(blocks[i]) & set(blocks[j]) != {center}:
                raise ValueError("Exceptional blocks have an extra common point")

    for i, block in enumerate(blocks):
        minimum = recipe.squares+1 if i in exceptional else recipe.squares+2
        if len(block) < minimum:
            raise ValueError("Insufficient local field capacity")
        field(len(block))  # rejects an unspecified or nonfield local order

    return {
        "mols": recipe.squares, "blocks": blocks,
        "exceptional_block_indices": exceptional,
    }


def construct(order: int) -> list[array]:
    """Return k flattened n-by-n arrays.

    Ordinary blocks: (1-t)*x+t*y, t=2,...,k+1.
    Exceptional blocks: x+t*y, t=1,...,k, with their center labelled zero.
    """
    design = build_design(order)
    k, n = design["mols"], order
    squares = [array("H", [0])*(n*n) for _ in range(k)]
    for x in range(n):
        for square in squares:
            square[n*x+x] = x

    exceptional = set(design["exceptional_block_indices"])
    local_tables = {}
    for block_index, block in enumerate(design["blocks"]):
        size = len(block)
        special = block_index in exceptional
        key = (size, special)
        if key not in local_tables:
            f = field(size)
            tables = []
            for j in range(k):
                t = j+1 if special else j+2
                a = 1 if special else f.sub(1, t)
                tables.append([[f.add(f.mul(a, x), f.mul(t, y))
                                for y in range(size)] for x in range(size)])
            local_tables[key] = tables
        for j, table in enumerate(local_tables[key]):
            square = squares[j]
            for x, global_x in enumerate(block):
                for y, global_y in enumerate(block):
                    if x != y or special:
                        square[n*global_x+global_y] = block[table[x][y]]
    return squares
