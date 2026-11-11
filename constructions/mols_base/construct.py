"""Construct the small-order families from the manuscript's fixed data.

Affine difference matrices use the same engine as the larger finite
constructions; the adapter below expands compact orbit
seeds into that engine's ordinary matrix representation.
"""
from .common import (array_to_squares, check_array, check_difference_matrix,
                     dm_to_squares, group_from_spec, qdm_to_array)
from .constructions import (affine_arrays, cyclic90, cyclic_difference,
                           difference_lifts, relative_difference,
                           small_quasi)

CASES = (22, 45, 50, 54, 72, 74, 78, 80, 85, 90, 92, 96, 98)


def affine_certificate(order):
    """Translate the printed DM seeds into the shared affine-engine schema."""
    data = difference_lifts.certificate(order)
    p, d, a = (data[key] for key in ('p', 'd', 'a'))
    maps = [[[column // p**digit % p for column in columns]
             for digit in range(d)] for columns in data['M']]
    return dict(order=order, mols=len(maps)-1, kind='difference_matrix',
                coarse_radices=data['coarse_moduli'], p=p, d=d, a=a,
                H=data['H'], maps=maps, offsets=data['C'])


def build(order):
    """Return a cyclic DM, QDM, or directly expanded OA certificate."""
    if order in (22, 74):
        return small_quasi.build(order)
    if order in (50, 54, 98):
        return affine_arrays.build(order)
    if order in relative_difference.ORDERS:
        return relative_difference.build(order)
    if order == 85:
        return cyclic_difference.build()
    if order == 90:
        return cyclic90.build()
    raise ValueError(f'Unsupported direct-construction order: {order}')


def construct(order):
    """Return k squares as k lists of n rows of n integer symbols."""
    if order not in CASES:
        raise ValueError(f'Unsupported order: {order}')
    if order in difference_lifts.ORDERS:
        from construction.affine import verify_compact, squares
        certificate = affine_certificate(order)
        verify_compact(certificate)
        return squares(certificate)[0]
    spec = build(order)
    if spec['kind'] == 'oa':
        rows, n = spec['rows'], spec['n']
    else:
        group = group_from_spec(spec)
        hole = spec.get('hole', 0)
        check_difference_matrix(spec['matrix'], group, hole)
        if spec['kind'] == 'dm':
            return dm_to_squares(spec['matrix'], group)[0]
        n = group.order + hole
        rows = qdm_to_array(spec['matrix'], group, hole)
    check_array(rows, n)
    return array_to_squares(rows, n)
