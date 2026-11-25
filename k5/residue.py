#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""The residue lemma for hyperplane sections in F_3^n (RESIDUE.md), and three explicit
points showing where the profile method stops.

Residue lemma.  For a multiset K of N points in F_3^n and a linear form xi, let
n_j = #{x in K : xi.x = j}, l = sum_x xi.x and R = sum_x (xi.x)^2 (both in F_3).  Since
[t=0] = 1 - t^2, [t=1] = 2t^2 + 2t and [t=2] = 2t^2 + t over F_3,

    n_0 = N - R,   n_1 = -R - l,   n_2 = -R + l      (mod 3).

l is a linear form and R a quadratic form in xi.  The quadratic form of QUAD45.md is
Q_K = N R - l^2, so for N = 0 (mod 3) Q_K = -l^2 has rank at most 1.

The profile LP (moment identities, condition (T), SD2 and all residue information) has the
integer points FEASIBLE below.  Each is a list of section profiles, one per hyperplane class,
with an integer line count L >= H, that meets every constraint of the method, together with
a linear form l and a quadratic form R whose predicted residues match the profiles class by
class.  So the method by itself cannot exclude 26 points in F_3^4, 45 points in F_3^5 or
78 points in F_3^6.
"""
from collections import Counter
from math import comb

try:
    from . import sd2, quad
except ImportError:
    import sd2
    import quad


def moments(K, xi):
    """(l, R) in F_3 for the class xi: l = sum xi.x, R = sum (xi.x)^2."""
    l = 0
    R = 0
    for x in K:
        t = sum(a * b for a, b in zip(xi, x)) % 3
        l += t
        R += t * t
    return l % 3, R % 3


def predicted_sections(N, l, R):
    """(n_0, n_1, n_2) mod 3 as given by the residue lemma."""
    return ((N - R) % 3, (-R - l) % 3, (-R + l) % 3)


def residue_multiset(p):
    return tuple(sorted(t % 3 for t in p))


def group_of(p, N):
    """'P' when the three sizes have distinct residues (needs N = 0 mod 3, l != 0); else the
    value r of R on the class, read off from sizes = {N - r, -r, -r} (mod 3); None if neither."""
    m = residue_multiset(p)
    for r in range(3):
        if residue_multiset(((N - r) % 3, (-r) % 3, (-r) % 3)) == m:
            return r
    if m == (0, 1, 2) and N % 3 == 0:
        return "P"
    return None


# ---- explicit points of the profile LP ------------------------------------------------
def _zero(xi):
    return 0


def _first(xi):
    return xi[0] % 3


def _second_sq(xi):
    return (xi[1] * xi[1]) % 3


FEASIBLE = {
    # F_3^4, N = 26 (Lev's 27-point set shows k_4 <= 27): s = 0, R = 0.
    "f3_4_n26": {"n": 4, "N": 26, "L": 53, "l": _zero, "R": _zero,
                 "counts": {(6, 6, 14): 10, (6, 8, 12): 1, (6, 9, 11): 1, (8, 9, 9): 28}},
    # F_3^5, N = 45: s = 0, R = 0 (every section size divisible by 3).
    "f3_5_n45": {"n": 5, "N": 45, "L": 125, "l": _zero, "R": _zero,
                 "counts": {(9, 12, 24): 23, (12, 12, 21): 1, (12, 15, 18): 1,
                            (15, 15, 15): 96}},
    # F_3^6, N = 78: l = xi_1 (s = e_1), R = xi_2^2.
    "f3_6_n78": {"n": 6, "N": 78, "L": 364, "l": _first, "R": _second_sq,
                 "counts": {(18, 18, 42): 39, (18, 21, 39): 1, (17, 19, 42): 3,
                            (25, 26, 27): 234, (24, 26, 28): 3, (24, 25, 29): 3,
                            (26, 26, 26): 81}},
}


def check_point(pt):
    """Return a list of failed checks (empty if the point satisfies every constraint)."""
    n, N, L, counts = pt["n"], pt["N"], pt["L"], pt["counts"]
    T, H, a1, a2 = sd2.constants(n)
    bad = []
    if sum(counts.values()) != H:
        bad.append("class count %d != %d" % (sum(counts.values()), H))
    for p in counts:
        if sum(p) != N:
            bad.append("profile %s does not sum to N" % (p,))
        if not sd2.allowed(tuple(sorted(p)), T, "sd2"):
            bad.append("profile %s not allowed by (T)/SD2" % (p,))
    S2 = sum(c * sd2.s2(p) for p, c in counts.items())
    S3 = sum(c * sd2.s3(p) for p, c in counts.items())
    if S2 != a1 * comb(N, 2):
        bad.append("pair identity: %d != %d" % (S2, a1 * comb(N, 2)))
    if S3 != a2 * comb(N, 3) + (a1 - a2) * L:
        bad.append("triple identity: %d != %d" % (S3, a2 * comb(N, 3) + (a1 - a2) * L))
    if L < H:
        bad.append("L < H")
    # residues: the multiset of predicted residue triples over all classes must equal the
    # multiset of residue triples of the profiles (then a class-by-class assignment exists)
    pred = Counter()
    for xi in quad.proj_points(n):
        pred[residue_multiset(predicted_sections(N, pt["l"](xi), pt["R"](xi)))] += 1
    have = Counter()
    for p, c in counts.items():
        have[residue_multiset(p)] += c
    if pred != have:
        bad.append("residues: predicted %s, profiles %s" % (dict(pred), dict(have)))
    return bad


def q_type(pt):
    """(z0, z1, z2) of Q_K = N R - l^2 over the classes, for the structure of a point."""
    z = [0, 0, 0]
    for xi in quad.proj_points(pt["n"]):
        l, R = pt["l"](xi), pt["R"](xi)
        z[(pt["N"] * R - l * l) % 3] += 1
    return z


if __name__ == "__main__":
    for name, pt in FEASIBLE.items():
        print(name, check_point(pt) or "all checks pass", "Q type", q_type(pt))
