# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Refined sum-difference bound (SD2) and the two-moment profile LP for Kakeya sets in F_3^n.

A Kakeya set K in F_3^n with |K| = N cuts every one of the H = (3^n - 1)/2 parallel
classes of hyperplanes into three sections of sizes (n1, n2, n3).  With T = 3^(n-1):

  sum over classes of C(n1,2)+C(n2,2)+C(n3,2) = a1 C(N,2)                 (pairs)
  sum over classes of C(n1,3)+C(n2,3)+C(n3,3) = a2 C(N,3) + (a1 - a2) L   (triples)

a1 = (3^(n-1)-1)/2, a2 = (3^(n-2)-1)/2, L = number of lines in K >= H.
(T): any two sections multiply to at least T.
SD  (moments/SD.md):  third section >= T^2 / Emax(x, y).
SD2 (this note):      third section >= T^2 / (Emax(x, y) - (x y - T)).
Proofs: SD2.md.  Exact integer / Fraction arithmetic throughout.
"""
from fractions import Fraction
from math import comb


def emax(x, y, T, cap=None):
    """Largest sum of r_v^2 over T integers with 1 <= r_v <= min(x, y, cap) and sum x*y.

    Returns None when no such integers exist (x*y < T, or too few columns)."""
    m = min(x, y) if cap is None else min(cap, x, y)
    W = x * y - T
    if W < 0:
        return None
    if W == 0:
        return T
    if m <= 1:
        return None
    q, rem = divmod(W, m - 1)
    cols = q + (1 if rem else 0)
    if cols > T:
        return None
    return q * m * m + ((1 + rem) ** 2 if rem else 0) + (T - cols)


def min_third(x, y, T, refined=True):
    """Least size of the third section next to sections of sizes x and y.

    None when x*y < T, i.e. condition (T) already fails."""
    E = emax(x, y, T)
    if E is None:
        return None
    den = E - (x * y - T) if refined else E
    return -(-(T * T) // den)


def allowed(p, T, rule="sd2"):
    """Is the profile p = (n1, n2, n3) allowed by (T) and, for rule sd1/sd2, by SD/SD2?"""
    a, b, c = p
    for (i, j, k) in ((a, b, c), (a, c, b), (b, c, a)):
        if i * j < T:
            return False
        if rule in ("sd1", "sd2"):
            if k < min_third(i, j, T, refined=(rule == "sd2")):
                return False
    return True


def profiles(N, T, rule="sd2"):
    out = []
    for a in range(N // 3 + 1):
        for b in range(a, (N - a) // 2 + 1):
            c = N - a - b
            if allowed((a, b, c), T, rule):
                out.append((a, b, c))
    return out


def s2(p):
    return sum(comb(t, 2) for t in p)


def s3(p):
    return sum(comb(t, 3) for t in p)


def constants(n):
    T = 3 ** (n - 1)
    H = (3 ** n - 1) // 2
    a1 = (3 ** (n - 1) - 1) // 2
    a2 = (3 ** (n - 2) - 1) // 2
    return T, H, a1, a2


def upper_hull_value(points, m):
    """Max of y over convex combinations of points (x, y) whose mean x equals m.

    None if m is outside the x-range."""
    best = {}
    for x, y in points:
        if x not in best or y > best[x]:
            best[x] = y
    pts = sorted(best.items())
    if not pts or m < pts[0][0] or m > pts[-1][0]:
        return None
    hull = []
    for p in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (y2 - y1) * (p[0] - x1) <= (p[1] - y1) * (x2 - x1):
                hull.pop()
            else:
                break
        hull.append(p)
    if len(hull) == 1:
        return Fraction(hull[0][1]) if hull[0][0] == m else None
    for (x1, y1), (x2, y2) in zip(hull, hull[1:]):
        if x1 <= m <= x2:
            return Fraction(y1) + Fraction(y2 - y1) * (m - x1) / (x2 - x1)
    return None


def lp_feasible(n, N, rule="sd2", exclude=()):
    """Can the two moment identities, L >= H and the allowed profiles hold at once (LP relaxation)?"""
    T, H, a1, a2 = constants(n)
    pts = [(s2(p), s3(p)) for p in profiles(N, T, rule) if p not in exclude]
    m = Fraction(a1 * comb(N, 2), H)
    need = Fraction(a2 * comb(N, 3) + (a1 - a2) * H, H)
    v = upper_hull_value(pts, m)
    return v is not None and v >= need


def scan(n, lo, hi, rule="sd2"):
    """List of N in [lo, hi] that the LP excludes.  k_n >= max(excluded) + 1."""
    return [N for N in range(lo, hi + 1) if not lp_feasible(n, N, rule)]
