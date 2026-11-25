#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Pencil-level linear program for Kakeya sets in F_3^n (k5/PENCIL.md).

A *pencil* is a codimension-2 subspace V of F_3^n.  Its 9 cosets ("cells") carry the counts
m(a, b) = #{x in K : xi1.x = a, xi2.x = b} of a set K, where xi1, xi2 span the annihilator of
V.  The 4 hyperplane classes through V (xi1, xi2, xi1+xi2, xi1+2 xi2) cut the 3 x 3 array
into its 4 parallel classes of lines; their line sums are the section profiles of those
classes.  Facts used (each proved in PENCIL.md):

  * each class lies in a1 = (3^(n-1)-1)/2 pencils; there are P = [n choose 2]_3 pencils;
    any two distinct classes lie in exactly one common pencil;
  * sum over pencils of sum over cells C(m,3) = B2 C(N,3) + (B1 - B2) L, where
    B1 = [n-1 choose 2]_3, B2 = [n-2 choose 2]_3, L = number of lines in K;
  * pencil Kakeya filter: for each of the 4 classes, the 3^(n-2) directions inside the class
    but not inside V need lines inside the class's sections; a section made of cells with
    a, b, c points holds lines in at most min(ab, bc, ca) directions; the (3^(n-2)-1)/2
    directions inside V need lines inside single cells;
  * line counts: every line of K lies inside a cell for B1 pencils and crosses three cells
    for P - B1 pencils (inside a section for a1 classes, crossing for 3^(n-1) classes).

The LP: variables x_p (classes with profile p), y_f (pencils with feature f) and
L' = L - H >= 0, all >= 0; the constraints of build_lp.  An infeasible LP, backed by an exact
Farkas certificate (check_certificate), proves that no Kakeya set of N points exists, and by
monotonicity (adding points keeps a set Kakeya) none of fewer points.

Exact integer / Fraction arithmetic throughout this module.  The solver lives in pencil_lp.py.
"""
import itertools
import random
from collections import Counter
from fractions import Fraction
from math import comb

try:
    from . import sd2
except ImportError:
    import sd2

DIRS = ((1, 0), (0, 1), (1, 1), (1, 2))
LINES = [[[a * 3 + b for a in range(3) for b in range(3) if (al * a + be * b) % 3 == c]
          for c in range(3)] for (al, be) in DIRS]
KMIN = {0: 1, 1: 3, 2: 7, 3: 13}       # smallest Kakeya sets in F_3^d (d = 3: Lev 2009)
MAXLINES2 = [0, 0, 0, 1, 1, 2, 3, 5, 8, 12]   # max lines in m points of AG(2,3)
MAXDIRS2 = [0, 0, 0, 1, 1, 2, 3, 4, 4, 4]     # max directions covered by m points of AG(2,3)


# ---- constants -------------------------------------------------------------------------
def gauss2(m):
    """Gaussian binomial [m choose 2]_3 (number of 2-dim subspaces of F_3^m)."""
    if m < 2:
        return 0
    return (3 ** m - 1) * (3 ** (m - 1) - 1) // 16


def consts(n):
    T, H, a1, a2 = sd2.constants(n)
    return {"T": T, "H": H, "a1": a1, "a2": a2, "P": gauss2(n), "B1": gauss2(n - 1),
            "B2": gauss2(n - 2), "D": 3 ** (n - 2), "DV": (3 ** (n - 2) - 1) // 2}


def max_lines(m, d):
    """Upper bound on the number of lines inside m points of AG(d, 3) (exact for d <= 2)."""
    if d == 1:
        return 1 if m >= 3 else 0
    if d == 2:
        return MAXLINES2[m]
    return m * (m - 1) // 6


def max_dirs(m, d):
    """Upper bound on the number of directions in which m points of AG(d, 3) contain a line."""
    if d == 1:
        return 1 if m >= 3 else 0
    if d == 2:
        return MAXDIRS2[m]
    dv = (3 ** d - 1) // 2
    cap = dv if m >= KMIN.get(d, 0) else dv - 1
    return min(m * (m - 1) // 6, cap)


# ---- profiles ----------------------------------------------------------------------------
def profile_list(n, N, rule="sd2"):
    T = 3 ** (n - 1)
    return [p for p in sd2.profiles(N, T, rule) if p[2] <= T]


def residue(p, N):
    """Residue group of a profile (the value of the quadratic form Q_K of QUAD45.md)."""
    return (comb(N, 2) - sd2.s2(p)) % 3


def variance(p, N):
    """3 sum p_j^2 - N^2 = 2 |Fourier coefficient|^2 of the class."""
    return 3 * sum(t * t for t in p) - N * N


def inside_class(p, n):
    return sum(max_lines(t, n - 1) for t in p)


def cross_class(p):
    a, b, c = p
    return min(a * b, b * c, a * c)


# ---- quadric types (copied from quad.py, QUAD45.md) -------------------------------------
def proj_points(n):
    pts = []
    for v in itertools.product(range(3), repeat=n):
        for c in v:
            if c:
                if c == 1:
                    pts.append(v)
                break
    return pts


def quadric_types(n):
    """(z0, z1, z2) for the 2n + 1 equivalence classes of quadratic forms over F_3."""
    pts = proj_points(n)
    forms = [[0] * n]
    for r in range(1, n + 1):
        for a in (1, 2):
            forms.append([1] * (r - 1) + [a] + [0] * (n - r))
    out = []
    for diag in forms:
        z = [0, 0, 0]
        for v in pts:
            z[sum(c * x * x for c, x in zip(diag, v)) % 3] += 1
        out.append(tuple(z))
    return out


# ---- arrays and features -----------------------------------------------------------------
def array_features(m, n, allowed):
    """Feature of a 3x3 pencil array m (tuple of 9, index a*3+b), or None if it fails a filter.

    Feature = (sorted tuple of the 4 sorted line-sum profiles, sum C(m,3), lines-inside bound,
    lines-crossing bound)."""
    k = consts(n)
    D, d = k["D"], n - 2
    if max(m) > D:
        return None
    profs = []
    cross = 0
    for L in LINES:
        s = tuple(sorted(sum(m[i] for i in l) for l in L))
        if s not in allowed:
            return None
        cov = 0
        for l in L:
            x, y, z = m[l[0]], m[l[1]], m[l[2]]
            mn = min(x * y, y * z, x * z)
            cross += mn
            cov += min(mn, D)
        if cov < D:
            return None
        profs.append(s)
    if sum(max_dirs(v, d) for v in m) < k["DV"]:
        return None
    inside = sum(max_lines(v, d) for v in m)
    tsame = sum(comb(v, 3) for v in m)
    return (tuple(sorted(profs)), tsame, inside, cross)


def enumerate_features(n, N, rule="sd2"):
    """All features of arrays that pass the filters, each with one witness array.

    Symmetry: AGL(2,3) permutes the 4 directions as S4 and permutes rows / columns freely, so
    every array is equivalent to one whose rows carry the lexicographically least profile,
    whose columns carry the least of the other three, both sorted ascending.  Features are
    invariant, so it suffices to enumerate those arrays."""
    k = consts(n)
    D = k["D"]
    profs = profile_list(n, N, rule)
    allowed = set(profs)
    feats = {}
    L2, L3 = LINES[2], LINES[3]
    for i, r in enumerate(profs):
        r0, r1, r2 = r
        for c in profs[i:]:
            c0, c1, c2 = c
            for m00 in range(min(r0, c0, D) + 1):
                for m01 in range(min(r0 - m00, c1, D) + 1):
                    m02 = r0 - m00 - m01
                    if m02 > c2 or m02 > D:
                        continue
                    for m10 in range(min(r1, c0 - m00, D) + 1):
                        m20 = c0 - m00 - m10
                        if m20 > r2 or m20 > D:
                            continue
                        for m11 in range(min(r1 - m10, c1 - m01, D) + 1):
                            m12 = r1 - m10 - m11
                            if m12 > c2 - m02 or m12 > D:
                                continue
                            m21 = c1 - m01 - m11
                            m22 = c2 - m02 - m12
                            if m21 > D or m22 > D or m20 + m21 + m22 != r2:
                                continue
                            m = (m00, m01, m02, m10, m11, m12, m20, m21, m22)
                            d2 = tuple(sorted(m[l[0]] + m[l[1]] + m[l[2]] for l in L2))
                            if d2 < c or d2 not in allowed:
                                continue
                            d3 = tuple(sorted(m[l[0]] + m[l[1]] + m[l[2]] for l in L3))
                            if d3 < c or d3 not in allowed:
                                continue
                            f = array_features(m, n, allowed)
                            if f is not None and f not in feats:
                                feats[f] = m
    return feats


def prune(feats):
    """Drop dominated features: among features with the same profile multiset and the same
    sum C(m,3), keep only the Pareto-maximal (inside, cross) pairs.  Valid because inside and
    cross appear only in the <= rows I1, I2 with negative coefficients: moving the weight of a
    dominated column to its dominator keeps every row satisfied.  So if the pruned LP is
    infeasible, so is the full one."""
    groups = {}
    for f in feats:
        groups.setdefault((f[0], f[1]), []).append((f[2], f[3]))
    out = []
    for key, vals in groups.items():
        vals = sorted(set(vals), reverse=True)
        best_cross = -1
        for ins, cr in vals:             # inside descending; keep if cross beats all before
            if cr > best_cross:
                out.append((key[0], key[1], ins, cr))
                best_cross = cr
    return sorted(out)


# ---- the LP ------------------------------------------------------------------------------
def unbalanced_set(profs, N, level=1):
    """Profiles whose variance exceeds the level-th smallest variance value among profs."""
    vals = sorted(set(variance(p, N) for p in profs))
    if len(vals) <= level - 1:
        return set()
    t = vals[level - 1]
    return set(p for p in profs if variance(p, N) > t)


def build_lp(n, N, profs, feats, z=None, U=None, u=None):
    """Rows of the LP.  Returns (eq, ub, bounds): eq / ub are lists of (dict col->int, rhs)
    meaning sum = rhs / sum <= rhs; bounds[j] is a valid upper bound for column j.
    Columns: x_p (p in profs), then y_f (f in feats), then L' = L - H."""
    k = consts(n)
    T, H, a1, a2 = k["T"], k["H"], k["a1"], k["a2"]
    P, B1, B2 = k["P"], k["B1"], k["B2"]
    feats = list(feats)
    npf, F = len(profs), len(feats)
    yo, Lc = npf, npf + len(feats)
    eq, ub = [], []
    eq.append(({i: 1 for i in range(npf)}, H))
    eq.append(({i: sd2.s2(p) for i, p in enumerate(profs) if sd2.s2(p)}, a1 * comb(N, 2)))
    r = {i: sd2.s3(p) for i, p in enumerate(profs) if sd2.s3(p)}
    r[Lc] = -(a1 - a2)
    eq.append((r, a2 * comb(N, 3) + (a1 - a2) * H))
    eq.append(({yo + j: 1 for j in range(F)}, P))
    pidx = {p: i for i, p in enumerate(profs)}
    link = [{i: -a1} for i in range(npf)]
    for j, f in enumerate(feats):
        for p, c in Counter(f[0]).items():
            link[pidx[p]][yo + j] = c
    for i in range(npf):
        eq.append((link[i], 0))
    r = {yo + j: f[1] for j, f in enumerate(feats) if f[1]}
    if B1 - B2:
        r[Lc] = -(B1 - B2)
    eq.append((r, B2 * comb(N, 3) + (B1 - B2) * H))
    if z is not None:
        for c in range(3):
            eq.append(({i: 1 for i, p in enumerate(profs) if residue(p, N) == c}, z[c]))
    if U is not None:
        eq.append(({i: 1 for i, p in enumerate(profs) if p in U}, u))
        r = {}
        for j, f in enumerate(feats):
            kk = sum(1 for p in f[0] if p in U)
            if kk >= 2:
                r[yo + j] = comb(kk, 2)
        eq.append((r, comb(u, 2)))
    r = {yo + j: -f[2] for j, f in enumerate(feats) if f[2]}
    r[Lc] = B1
    ub.append((r, -B1 * H))
    r = {yo + j: -f[3] for j, f in enumerate(feats) if f[3]}
    r[Lc] = P - B1
    ub.append((r, -(P - B1) * H))
    r = {i: -inside_class(p, n) for i, p in enumerate(profs) if inside_class(p, n)}
    r[Lc] = a1
    ub.append((r, -a1 * H))
    r = {i: -cross_class(p) for i, p in enumerate(profs) if cross_class(p)}
    r[Lc] = T
    ub.append((r, -T * H))
    bounds = [H] * npf + [P] * F + [comb(N, 2) // 3]
    return eq, ub, bounds


def check_point(eq, ub, v):
    """Does the integer vector v satisfy every row exactly?  Returns list of failed rows."""
    bad = []
    for k, (r, rhs) in enumerate(eq):
        if sum(c * v[j] for j, c in r.items()) != rhs:
            bad.append(("eq", k))
    for k, (r, rhs) in enumerate(ub):
        if sum(c * v[j] for j, c in r.items()) > rhs:
            bad.append(("ub", k))
    return bad


def check_certificate(eq, ub, bounds, lam, mu):
    """Exact Farkas check.  lam: {eq row: value}, mu: {ub row: value >= 0} (Fractions or
    strings).  With r_j = sum_i lam_i A_ij + sum_k mu_k G_kj, any v with A v = b, G v <= h,
    0 <= v <= bounds has sum_j r_j v_j <= b.lam + h.mu, while sum_j r_j v_j >=
    sum_{r_j < 0} r_j bounds_j.  So S = b.lam + h.mu - sum_{r_j<0} r_j bounds_j < 0 proves
    that no such v exists.  Returns (proved, S)."""
    lam = {int(i): Fraction(x) for i, x in lam.items()}
    mu = {int(i): Fraction(x) for i, x in mu.items()}
    if any(x < 0 for x in mu.values()):
        return False, None
    red = Counter()
    S = Fraction(0)
    for i, x in lam.items():
        r, rhs = eq[i]
        S += x * rhs
        for j, c in r.items():
            red[j] += x * c
    for i, x in mu.items():
        r, rhs = ub[i]
        S += x * rhs
        for j, c in r.items():
            red[j] += x * c
    for j, val in red.items():
        if val < 0:
            S -= val * bounds[j]
    return S < 0, S


# ---- real sets (for controls) ------------------------------------------------------------
def dot(u, v):
    return sum(a * b for a, b in zip(u, v)) % 3


def norm(v):
    for c in v:
        if c:
            return tuple(v) if c == 1 else tuple((2 * a) % 3 for a in v)
    return tuple(v)


def pencil_bases(n):
    """One basis (xi1, xi2) for each of the [n choose 2]_3 pencils."""
    pts = proj_points(n)
    seen, out = set(), []
    for i, x in enumerate(pts):
        for y in pts[i + 1:]:
            key = frozenset(norm(tuple((a * s + b * t) % 3 for s, t in zip(x, y)))
                            for a, b in DIRS)
            if key not in seen:
                seen.add(key)
                out.append((x, y))
    return out


def is_kakeya(K, n):
    S = set(tuple(x) for x in K)
    for d in proj_points(n):
        if not any(tuple((a + b) % 3 for a, b in zip(x, d)) in S and
                   tuple((a + 2 * b) % 3 for a, b in zip(x, d)) in S for x in S):
            return False
    return True


def count_lines(K):
    S = set(tuple(x) for x in K)
    Kl = list(S)
    c = 0
    for i in range(len(Kl)):
        for j in range(i + 1, len(Kl)):
            t = tuple((-(a + b)) % 3 for a, b in zip(Kl[i], Kl[j]))
            if t in S:
                c += 1
    return c // 3


def real_data(K, n, rule="sd2"):
    """Profiles per class, pencil arrays, features and L of an actual set K."""
    N = len(K)
    allowed = set(profile_list(n, N, rule))
    classes = Counter()
    for xi in proj_points(n):
        s = [0, 0, 0]
        for x in K:
            s[dot(xi, x)] += 1
        classes[tuple(sorted(s))] += 1
    arrays, feats = [], Counter()
    for x1, x2 in pencil_bases(n):
        m = [0] * 9
        for x in K:
            m[dot(x1, x) * 3 + dot(x2, x)] += 1
        m = tuple(m)
        arrays.append(m)
        feats[array_features(m, n, allowed)] += 1
    return {"N": N, "classes": classes, "arrays": arrays, "feats": feats, "L": count_lines(K)}


def real_vector(n, data, profs, featlist):
    """The LP vector (x, y, L') of an actual set, columns as in build_lp."""
    H = consts(n)["H"]
    v = [data["classes"].get(p, 0) for p in profs]
    v += [data["feats"].get(f, 0) for f in featlist]
    v.append(data["L"] - H)
    return v


def greedy_kakeya(n, seed, reduce=True):
    """A random Kakeya set: add one line per uncovered direction (reusing points), then
    optionally remove points while the set stays Kakeya."""
    rng = random.Random(seed)
    dirs = proj_points(n)
    rng.shuffle(dirs)
    allp = list(itertools.product(range(3), repeat=n))
    K = set()
    for d in dirs:
        if any(tuple((a + b) % 3 for a, b in zip(x, d)) in K and
               tuple((a + 2 * b) % 3 for a, b in zip(x, d)) in K for x in K):
            continue
        best, bc = [], -1
        for p in allp:
            line = [p, tuple((a + b) % 3 for a, b in zip(p, d)),
                    tuple((a + 2 * b) % 3 for a, b in zip(p, d))]
            c = sum(1 for q in line if q in K)
            if c > bc:
                best, bc = [line], c
            elif c == bc:
                best.append(line)
        K.update(rng.choice(best))
    if reduce:
        order = sorted(K)
        rng.shuffle(order)
        for p in order:
            K.discard(p)
            if not is_kakeya(K, n):
                K.add(p)
    return sorted(K)
