# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Stand-alone Kakeya check on coordinate tuples.  Uses nothing from gf.py, ksym.py or bb.py:
only the definition (a set K in F_3^n is Kakeya if for every nonzero direction d it contains
a whole line {x, x + d, x + 2d})."""
import itertools


def to_tuples(n, ints):
    return [tuple((x // 3 ** i) % 3 for i in range(n)) for x in ints]


def is_kakeya(n, K):
    S = set(tuple(x) for x in K)
    for d in itertools.product(range(3), repeat=n):
        nz = [c for c in d if c]
        if not nz or nz[0] != 1:
            continue
        ok = False
        for x in S:
            if (tuple((a + b) % 3 for a, b in zip(x, d)) in S and
                    tuple((a + 2 * b) % 3 for a, b in zip(x, d)) in S):
                ok = True
                break
        if not ok:
            return False
    return True


def classical(n):
    """{0,1}^n u {0,2}^n, the classical Kakeya set with 2^(n+1) - 1 points."""
    S = set()
    for pat in itertools.product((0, 1), repeat=n):
        S.add(pat)
        S.add(tuple(2 * c for c in pat))
    return sorted(S)
