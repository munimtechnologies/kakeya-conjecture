# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Stand-alone Kakeya check in F_3^n from the definition (shares no code with any search).

A point is stored as the integer sum_k x_k 3^(n-1-k) (first coordinate most significant).
K is Kakeya iff for every direction d (nonzero, first nonzero entry 1) some x in K has
x + d and x + 2d in K."""
import itertools


def decode(n, i):
    return tuple((i // 3 ** (n - 1 - k)) % 3 for k in range(n))


def directions(n):
    for d in itertools.product(range(3), repeat=n):
        nz = [c for c in d if c]
        if nz and nz[0] == 1:
            yield d


def uncovered(n, idx, stop_after=None):
    pts = [decode(n, i) for i in idx]
    S = set(pts)
    bad = []
    for d in directions(n):
        if not any(tuple((a + b) % 3 for a, b in zip(x, d)) in S and
                   tuple((a + 2 * b) % 3 for a, b in zip(x, d)) in S for x in pts):
            bad.append(d)
            if stop_after is not None and len(bad) >= stop_after:
                break
    return bad


def is_kakeya(n, idx):
    if len(set(idx)) != len(idx) or not all(0 <= i < 3 ** n for i in idx):
        return False
    return not uncovered(n, idx, stop_after=1)
