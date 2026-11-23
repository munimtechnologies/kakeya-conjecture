# ag33_lines.py -- how many line directions can a small point set in AG(3,3) carry?
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0. Standard library only.
#
# fits(k, cap) is True exactly when some k lines of AG(3,3) = F_3^3, in k pairwise distinct
# directions, have a union of at most cap points. The search is exhaustive. It may assume the
# configuration contains one fixed line `start` (through the origin, in the first direction),
# because the affine group of F_3^3 is transitive on lines and maps directions to directions
# bijectively; the other lines are then chosen in increasing direction order.
#
# D(m), the largest number of directions of lines inside an m-point subset of AG(3,3), is the
# largest k with fits(k, m). moments/N24.md proves by hand D(5) = 2, D(6) = 3, D(7) = D(8) = 4;
# test_n24.py checks those values with this search.

import time
from itertools import product


class Timeout(Exception):
    pass


def _popcount(x):
    return bin(x).count("1")


def geometry():
    """(points, index, directions, lines) of AG(3,3). lines[d] is the sorted list of the 9 lines
    (as bitmasks over the 27 point indices) with direction directions[d]."""
    pts = list(product(range(3), repeat=3))
    idx = {p: i for i, p in enumerate(pts)}
    dirs = [v for v in pts if any(v) and next(c for c in v if c) == 1]
    lines = []
    for v in dirs:
        masks = set()
        for p in pts:
            m = 0
            for t in range(3):
                q = tuple((p[i] + t * v[i]) % 3 for i in range(3))
                m |= 1 << idx[q]
            masks.add(m)
        lines.append(sorted(masks))
    return pts, idx, dirs, lines


def mask_of(points):
    _, idx, _, _ = geometry()
    m = 0
    for p in points:
        m |= 1 << idx[tuple(c % 3 for c in p)]
    return m


def directions_inside(mask):
    """Number of directions d such that some line of direction d lies inside the point set."""
    _, _, _, lines = geometry()
    return sum(1 for fam in lines if any(m & ~mask == 0 for m in fam))


def fits(k, cap, deadline=None):
    if k <= 0:
        return True
    if cap < 3:
        return False
    _, _, _, lines = geometry()
    nd = len(lines)
    start = next(m for m in lines[0] if m & 1)

    def rec(mask, count, d):
        if count == k:
            return True
        if count + (nd - d) < k:
            return False
        if deadline is not None and time.monotonic() > deadline:
            raise Timeout()
        for dd in range(d, nd):
            for m in lines[dd]:
                nm = mask | m
                if _popcount(nm) <= cap and rec(nm, count + 1, dd + 1):
                    return True
        return False

    return rec(start, 1, 1)


def max_directions(cap, deadline=None):
    k = 0
    while k < 13 and fits(k + 1, cap, deadline):
        k += 1
    return k
