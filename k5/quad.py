#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""The quadric constraint on hyperplane sections of Kakeya sets in F_3^n (QUAD45.md).

For K in F_3^n with |K| = N and a hyperplane class xi (a nonzero linear form up to sign),
s2(xi) = sum over the three sections of C(size, 2) counts the pairs {x, y} of K with
xi(x - y) = 0.  Over F_3, [t = 0] = 1 - t^2, so

    s2(xi) = C(N, 2) - Q_K(xi)   (mod 3),   Q_K(xi) = sum over pairs of xi(x - y)^2,

and Q_K is a quadratic form in xi.  Every quadratic form over F_3 is equivalent to 0 or to
diag(1, ..., 1, a) of some rank r with a in {1, 2}, so the numbers (z0, z1, z2) of classes on
which Q_K takes the values 0, 1, 2 are one of 2n + 1 vectors.  A class whose profile p has
(C(N,2) - s2(p)) mod 3 = c must lie where Q_K = c.

This gives a sharper profile LP: for a quadric type z, the profile counts x_p >= 0 must put
exactly z_c classes into residue group c, satisfy the pair identity exactly and the triple
identity with L >= H.  Weak duality: for any real mu,

    sum_p x_p s3(p) <= phi(mu) = mu * S2 + sum_c z_c * max_{p in group c} (s3(p) - mu s2(p)),

so phi(mu) < (needed s3 total) excludes the type.  If every type is excluded, no Kakeya set of
size N exists in F_3^n, and (adding points keeps a set Kakeya) none of any smaller size.
"""
import argparse
import itertools
import json
import time
from fractions import Fraction
from math import comb

try:
    from . import sd2
except ImportError:
    import sd2


# ---- projective points and quadratic forms over F_3 --------------------------------
def proj_points(n):
    """One vector for each of the (3^n - 1)/2 points of PG(n-1, 3): first nonzero entry 1."""
    pts = []
    for v in itertools.product(range(3), repeat=n):
        for c in v:
            if c:
                if c == 1:
                    pts.append(v)
                break
    return pts


def monomials(n):
    return [(i, j) for i in range(n) for j in range(i, n)]


def value_counts(coefs, pts, mons):
    """(z0, z1, z2): number of projective points where the form sum c_ij x_i x_j is 0, 1, 2."""
    z = [0, 0, 0]
    for v in pts:
        s = 0
        for c, (i, j) in zip(coefs, mons):
            if c:
                s += c * v[i] * v[j]
        z[s % 3] += 1
    return tuple(z)


def type_reps(n):
    """Representatives of all equivalence classes of quadratic forms in n variables over F_3:
    the zero form and diag(1, ..., 1, a) of rank r = 1..n, a in {1, 2}."""
    mons = monomials(n)
    reps = [("rank0", [0] * len(mons))]
    for r in range(1, n + 1):
        for a in (1, 2):
            coefs = []
            for (i, j) in mons:
                if i == j and i < r - 1:
                    coefs.append(1)
                elif i == j == r - 1:
                    coefs.append(a)
                else:
                    coefs.append(0)
            reps.append(("rank%d-a%d" % (r, a), coefs))
    return reps


def quadric_types(n):
    """List of (name, (z0, z1, z2)) for the 2n + 1 types."""
    pts = proj_points(n)
    mons = monomials(n)
    return [(name, value_counts(c, pts, mons)) for name, c in type_reps(n)]


# ---- point sets: sections and the form Q_K ------------------------------------------
def sections(K, xi):
    z = [0, 0, 0]
    for x in K:
        z[sum(a * b for a, b in zip(xi, x)) % 3] += 1
    return z


def s2_of(K, xi):
    return sum(comb(t, 2) for t in sections(K, xi))


def QK(K, xi):
    """Q_K(xi) = sum over unordered pairs {x, y} of K of xi(x - y)^2, in F_3."""
    s = 0
    for i in range(len(K)):
        for j in range(i + 1, len(K)):
            t = sum(a * (b - c) for a, b, c in zip(xi, K[i], K[j])) % 3
            s += t * t
    return s % 3


def is_kakeya(K, n):
    Ks = set(tuple(x) for x in K)
    for d in proj_points(n):
        ok = False
        for x in Ks:
            y = tuple((a + b) % 3 for a, b in zip(x, d))
            w = tuple((a + 2 * b) % 3 for a, b in zip(x, d))
            if y in Ks and w in Ks:
                ok = True
                break
        if not ok:
            return False
    return True


# ---- the quadric-refined profile LP ---------------------------------------------------
def residue(p, N):
    """Value of Q_K on a class whose sections have sizes p."""
    return (comb(N, 2) - sd2.s2(p)) % 3


def groups(n, N, rule="sd2"):
    T = 3 ** (n - 1)
    g = {0: [], 1: [], 2: []}
    for p in sd2.profiles(N, T, rule):
        g[residue(p, N)].append(p)
    return g


def totals(n, N):
    T, H, a1, a2 = sd2.constants(n)
    S2 = a1 * comb(N, 2)
    need = a2 * comb(N, 3) + (a1 - a2) * H
    return S2, need


def phi(grp, z, mu, S2):
    """Weak-duality bound on sum x_p s3(p); None when some z_c > 0 has an empty group."""
    tot = mu * S2
    for c in range(3):
        if z[c] == 0:
            continue
        if not grp[c]:
            return None
        tot += z[c] * max(sd2.s3(p) - mu * sd2.s2(p) for p in grp[c])
    return tot


def upper_hull(points):
    best = {}
    for x, y in points:
        if x not in best or y > best[x]:
            best[x] = y
    pts = sorted(best.items())
    hull = []
    for p in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1) >= 0:
                hull.pop()
            else:
                break
        hull.append(p)
    return hull


def candidate_mus(pts):
    """Every breakpoint of mu -> max (s3 - mu s2) over each group (upper-hull slopes), plus
    far points on both sides (to detect an unbounded LP dual)."""
    mus = set()
    for c in range(3):
        if not pts[c]:
            continue
        h = upper_hull(pts[c])
        for (x1, y1), (x2, y2) in zip(h, h[1:]):
            mus.add(Fraction(y2 - y1, x2 - x1))
    if not mus:
        mus.add(Fraction(0))
    lo, hi = min(mus), max(mus)
    for j in range(0, 48):
        mus.add(hi + 2 ** j)
        mus.add(lo - 2 ** j)
    return sorted(mus)


class Dual:
    """Precomputed group maxima M_c(mu) = max over group c of (s3 - mu s2) at all candidate mu."""

    def __init__(self, n, N, rule="sd2"):
        self.n, self.N = n, N
        self.grp = groups(n, N, rule)
        self.S2, self.need = totals(n, N)
        self.pts = {c: [(sd2.s2(p), sd2.s3(p)) for p in self.grp[c]] for c in range(3)}
        self.mus = candidate_mus(self.pts)
        self.M = []
        for mu in self.mus:
            self.M.append({c: max(y - mu * x for x, y in self.pts[c])
                           for c in range(3) if self.pts[c]})

    def certify(self, z):
        for c in range(3):
            if z[c] > 0 and not self.pts[c]:
                return {"excluded": True, "reason": "empty residue group %d" % c, "mu": None,
                        "phi": None, "need": self.need}
        best_mu, best_val = None, None
        for mu, M in zip(self.mus, self.M):
            v = mu * self.S2 + sum(z[c] * M[c] for c in range(3) if z[c])
            if best_val is None or v < best_val:
                best_mu, best_val = mu, v
        return {"excluded": best_val < self.need, "reason": "lp", "mu": best_mu,
                "phi": best_val, "need": self.need}


def certify(n, N, z, rule="sd2"):
    """Return a dict: excluded (bool), mu (certificate), phi, need."""
    return Dual(n, N, rule).certify(z)


def check_certificate(n, N, z, mu, rule="sd2"):
    """Independent re-check of one certificate by direct evaluation over all profiles."""
    T = 3 ** (n - 1)
    S2, need = totals(n, N)
    tot = Fraction(mu) * S2
    for c in range(3):
        if z[c] == 0:
            continue
        vals = [sd2.s3(p) - Fraction(mu) * sd2.s2(p)
                for p in sd2.profiles(N, T, rule) if (comb(N, 2) - sd2.s2(p)) % 3 == c]
        if not vals:
            return True
        tot += z[c] * max(vals)
    return tot < need


def excluded(n, N, rule="sd2", types=None):
    if types is None:
        types = quadric_types(n)
    dual = Dual(n, N, rule)
    per = {}
    for name, z in types:
        per[name] = dict(dual.certify(z), z=list(z))
    return all(v["excluded"] for v in per.values()), per


def scan(n, lo, hi, rule="sd2"):
    types = quadric_types(n)
    out = {}
    for N in range(lo, hi + 1):
        ok, per = excluded(n, N, rule, types)
        out[N] = {"excluded": ok,
                  "types": {k: {"z": v["z"], "excluded": v["excluded"], "reason": v["reason"],
                                "mu": str(v["mu"]), "phi": str(v["phi"]), "need": v["need"]}
                            for k, v in per.items()}}
    ex = [N for N in out if out[N]["excluded"]]
    return out, (max(ex) + 1 if ex else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="quad.json")
    args = ap.parse_args()
    t0 = time.time()
    res = {"contract": "k5quad-v1", "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "types": {}, "scan": {}, "bounds": {}, "old_lp_bounds": {}}
    for n, lo, hi in ((3, 7, 13), (4, 18, 27), (5, 36, 53), (6, 66, 100), (7, 120, 160)):
        res["types"][n] = quadric_types(n)
        out, bound = scan(n, lo, hi)
        res["scan"][n] = out
        res["bounds"][n] = bound
        old = sd2.scan(n, lo, hi, "sd2")
        res["old_lp_bounds"][n] = (max(old) + 1) if old else None
        print("n=%d quadric-LP bound %s (moment LP + SD2: %s) [%.1fs]" % (
            n, bound, res["old_lp_bounds"][n], time.time() - t0), flush=True)
        with open(args.out, "w") as f:
            json.dump(res, f, indent=1, default=str)
    res["seconds"] = round(time.time() - t0, 1)
    with open(args.out, "w") as f:
        json.dump(res, f, indent=1, default=str)


if __name__ == "__main__":
    main()
