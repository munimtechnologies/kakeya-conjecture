#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""k5cut: can three sections of sizes (x, y, z) in F_3^d carry all 3^d crossing lines?

Three parallel hyperplanes of F_3^(d+1) are copies of G = F_3^d holding sections A, B, C.
The line with direction (v, 1) crosses them in a, a+v, a+2v, so the sections carry all
crossing directions iff

    for every v in G there is a in A with a + v in B and a + 2v in C.          (*)

Columns r(v) = |A cap (B - v)| count pairs (a, b) with b - a = v; rows r+(s) count pairs
with a + b = s; the third point of the line through a and b is -(a + b).  So (*) holds with
|C| <= z iff B - A = G and some z rows (sums) cover every column.

Necessary conditions (proved in SD2.md), used only to prune:
  (E)  E = sum_v r(v)^2 >= T^2/z + x*y - T, T = 3^d;
  (K)  k1 = max_v r(v) is at least the least cap whose Emax reaches (E);
  (D)  translate B so that a largest column is v = 0 and put I = A cap B (|I| = k1),
       A = I + XA, B = I + Y.  Then B - A lies in (I-I) + (I-XA) + (Y-I) + (Y-XA), and
       |I - I| = 1 + 2 D(I) with D(I) = number of directions spanned by I, so
       D(I) >= ceil((T - 1 - k1 |XA| - k1 |Y| - |XA||Y|) / 2).
I runs over AGL(d,3)-orbit representatives (all of them with the needed excess);
XA and Y are enumerated in full; every (A, B) passing (E) is decided exactly by a set-cover
search over rows.  Any (A, B, C) found is re-checked against (*) directly.
"""
import argparse
import itertools
import json
import os
import sys
import time
from fractions import Fraction
from math import comb

try:
    from . import sd2
except ImportError:
    import sd2


def popcount(m):
    return bin(m).count("1")


def gl_order(k):
    r = 1
    for i in range(k):
        r *= 3 ** k - 3 ** i
    return r


def agl_order(d):
    return 3 ** d * gl_order(d)


class Space:
    """The affine space AG(d, 3); points are integers 0 .. 3^d - 1 (base-3 digits)."""

    def __init__(self, d):
        self.d = d
        n = 3 ** d
        self.n = n
        self.coords = [tuple((p // 3 ** i) % 3 for i in range(d)) for p in range(n)]

        def enc(c):
            return sum((c[i] % 3) * 3 ** i for i in range(d))

        self.add = [[enc([a + b for a, b in zip(self.coords[p], self.coords[q])])
                     for q in range(n)] for p in range(n)]
        self.neg = [enc([-a for a in self.coords[p]]) for p in range(n)]
        self.sub = [[self.add[p][self.neg[q]] for q in range(n)] for p in range(n)]
        self.dirof = [0] + [min(v, self.neg[v]) for v in range(1, n)]
        self.full = (1 << n) - 1

    # ---- directions and excess -------------------------------------------------
    def ndirections(self, S):
        sub, dirof = self.sub, self.dirof
        dirs = set()
        for i in range(len(S)):
            si = S[i]
            for j in range(i + 1, len(S)):
                dirs.add(dirof[sub[S[j]][si]])
        return len(dirs)

    def excess(self, S):
        """Pairs minus distinct directions; never increases when points are removed."""
        return comb(len(S), 2) - self.ndirections(S)

    # ---- canonical form under AGL(d, 3) ---------------------------------------
    def canon(self, S):
        """Return (key, aut, k): key identifies the AGL(d,3)-orbit of S, k is the affine
        dimension of S and aut the order of the stabiliser of S inside AGL(k,3) acting
        on the affine span (counted as the number of allowed frames giving the key)."""
        S = list(S)
        m = len(S)
        add, sub, neg, dirof = self.add, self.sub, self.neg, self.dirof
        if m == 1:
            return (0, (0,)), 1, 0
        Sset = set(S)
        mult = {}
        for i in range(m):
            for j in range(i + 1, m):
                dd = dirof[sub[S[j]][S[i]]]
                mult[dd] = mult.get(dd, 0) + 1
        inv = {}
        for p in S:
            lines = 0
            shared = 0
            for q in S:
                if q == p:
                    continue
                if neg[add[p][q]] in Sset:
                    lines += 1
                if mult[dirof[sub[q][p]]] >= 2:
                    shared += 1
            inv[p] = (lines, shared)
        top = max(inv.values())
        starts = [p for p in S if inv[p] == top]
        # affine dimension
        base = S[0]
        span = {0}
        k = 0
        for s in S:
            v = sub[s][base]
            if v not in span:
                v2 = add[v][v]
                span = span | {add[u][v] for u in span} | {add[u][v2] for u in span}
                k += 1
        pw = [3 ** i for i in range(k)]
        state = {"best": None, "count": 0}

        def offer(img):
            b = state["best"]
            if b is None or img < b:
                state["best"] = img
                state["count"] = 1
            elif img == b:
                state["count"] += 1

        for q0 in starts:
            vecs = [sub[s][q0] for s in S]
            # first basis point: maximise an AGL-invariant of the pair (q0, q1)
            phi = {}
            for q in S:
                if q != q0:
                    phi[q] = (mult[dirof[sub[q][q0]]], inv[q])
            ptop = max(phi.values())
            first = [sub[q][q0] for q in S if q != q0 and phi[q] == ptop]

            def finish(lookup, b, p):
                img = []
                for v in vecs:
                    t = lookup.get(v)
                    if t is not None:
                        img.append(t)
                        continue
                    t = lookup.get(sub[v][b])
                    if t is not None:
                        img.append(t + p)
                    else:
                        img.append(lookup[add[v][b]] + 2 * p)
                img.sort()
                return tuple(img)

            def rec(depth, lookup, choices):
                for v in choices:
                    if v in lookup:
                        continue
                    if depth == k - 1:
                        offer(finish(lookup, v, pw[depth]))
                    else:
                        p = pw[depth]
                        v2 = add[v][v]
                        newl = dict(lookup)
                        for u, im in lookup.items():
                            newl[add[u][v]] = im + p
                            newl[add[u][v2]] = im + 2 * p
                        rec(depth + 1, newl, vecs)

            rec(0, {0: 0}, first)
        return (k, state["best"]), state["count"], k

    def orbit_size(self, aut, k):
        d = self.d
        den = aut * 3 ** (k * (d - k)) * gl_order(d - k)
        g = agl_order(d)
        if g % den:
            raise AssertionError("stabiliser order does not divide |AGL|: %r" % ((aut, k),))
        return g // den

    def orbit_levels(self, maxk, xbound, log=None):
        """Orbit representatives of k-subsets (1 <= k <= maxk) with excess <= xbound[k].

        xbound must be non-increasing in k for completeness (excess is monotone)."""
        levels = {1: [((0,), self.orbit_size(1, 0))]}
        for k in range(2, maxk + 1):
            t0 = time.time()
            seen = {}
            for (P, _) in levels[k - 1]:
                Pset = set(P)
                for p in range(self.n):
                    if p in Pset:
                        continue
                    S = P + (p,)
                    if self.excess(S) > xbound[k]:
                        continue
                    key, aut, kk = self.canon(S)
                    if key not in seen:
                        seen[key] = (tuple(sorted(S)), self.orbit_size(aut, kk))
            levels[k] = list(seen.values())
            if log:
                log("orbits d=%d k=%d xbound=%s reps=%d (%.1fs)" % (
                    self.d, k, xbound[k], len(levels[k]), time.time() - t0))
        return levels


# ---- exact row cover --------------------------------------------------------------
def cover(rows, full, z):
    """rows: dict sum -> column mask.  Return a list of at most z sums whose masks cover
    `full`, or None if there is none."""
    items = [(s, r) for s, r in rows.items() if r]

    def dfs(unc, left, chosen):
        if unc == 0:
            return list(chosen)
        if left == 0:
            return None
        bestr = 0
        for _, r in items:
            c = popcount(r & unc)
            if c > bestr:
                bestr = c
        if bestr * left < popcount(unc):
            return None
        bestlist = None
        u = unc
        while u:
            low = u & -u
            lst = [(s, r) for s, r in items if r & low]
            if bestlist is None or len(lst) < len(bestlist):
                bestlist = lst
                if len(lst) <= 1:
                    break
            u ^= low
        for s, r in bestlist:
            chosen.append(s)
            got = dfs(unc & ~r, left - 1, chosen)
            chosen.pop()
            if got is not None:
                return got
        return None

    return dfs(full, z, [])


def rows_and_columns(sp, A, B):
    add, sub = sp.add, sp.sub
    r = [0] * sp.n
    rows = {}
    for a in A:
        for b in B:
            v = sub[b][a]
            r[v] += 1
            s = add[a][b]
            rows[s] = rows.get(s, 0) | (1 << v)
    return r, rows


def min_cover_size(sp, A, B):
    r, rows = rows_and_columns(sp, A, B)
    if min(r) == 0:
        return None
    for z in range(1, sp.n + 1):
        if cover(rows, sp.full, z) is not None:
            return z
    return None


def check_star(sp, A, B, C):
    """Direct check of (*)."""
    Bs, Cs = set(B), set(C)
    add = sp.add
    for v in range(sp.n):
        v2 = add[v][v]
        if not any(add[a][v] in Bs and add[a][v2] in Cs for a in A):
            return False
    return True


# ---- the pruned search ------------------------------------------------------------
def plan(sp, x, y, z):
    """Necessary conditions (E), (K), (D): returns Emin and {k1: Xmax}."""
    T = sp.n
    Ef = Fraction(T * T, z) + x * y - T
    Emin = -(-Ef.numerator // Ef.denominator)
    need = {}
    if x * y < T:
        return Emin, need
    for k1 in range(1, min(x, y) + 1):
        e = sd2.emax(x, y, T, cap=k1)
        if e is None or e < Emin:
            continue
        nx, ny = x - k1, y - k1
        Dmin = max(0, -(-(T - 1 - k1 * nx - k1 * ny - nx * ny) // 2))
        Xmax = comb(k1, 2) - Dmin
        if Xmax >= 0:
            need[k1] = Xmax
    return Emin, need


def xbounds(need, maxk):
    xb = {}
    for k in range(1, maxk + 1):
        later = [need[j] for j in need if j >= k]
        xb[k] = max(later) if later else -1
    return xb


def search_rep(sp, x, y, z, k1, I, Emin, max_solutions=3, deadline=None):
    """All (A, B) with A cap B = I (a fixed orbit representative), largest column k1 at 0."""
    T, full, add, sub = sp.n, sp.full, sp.add, sp.sub
    nx, ny = x - k1, y - k1
    stats = {"XA": 0, "XA_kept": 0, "B": 0, "maxcol_skip": 0, "E_pass": 0, "solutions": 0}
    sols = []
    Iset = set(I)
    II = 0
    for a in I:
        for b in I:
            II |= 1 << sub[b][a]
    outside = [p for p in range(T) if p not in Iset]
    for XA in itertools.combinations(outside, nx):
        stats["XA"] += 1
        if deadline and time.time() > deadline:
            stats["timeout"] = True
            break
        A = tuple(I) + XA
        Aset = set(A)
        cov = II
        for b in I:
            for a in XA:
                cov |= 1 << sub[b][a]
        U = full & ~cov
        if popcount(U) > ny * x:
            continue
        stats["XA_kept"] += 1
        pool = [p for p in range(T) if p not in Aset]
        M = {}
        for yv in pool:
            mm = 0
            for a in A:
                mm |= 1 << sub[yv][a]
            M[yv] = mm
        Ys = set()

        def dfs(Ysel, Uc):
            if len(Ysel) == ny:
                if Uc == 0:
                    Ys.add(frozenset(Ysel))
                return
            if Uc == 0:
                rest = [p for p in pool if p not in Ysel]
                for extra in itertools.combinations(rest, ny - len(Ysel)):
                    Ys.add(frozenset(Ysel + extra))
                return
            low = Uc & -Uc
            u = low.bit_length() - 1
            for a in A:
                yv = add[u][a]
                if yv in Aset or yv in Ysel:
                    continue
                dfs(Ysel + (yv,), Uc & ~M[yv])

        dfs((), U)
        for Y in Ys:
            stats["B"] += 1
            B = tuple(I) + tuple(sorted(Y))
            r, rows = rows_and_columns(sp, A, B)
            if min(r) == 0:
                raise AssertionError("coverage bookkeeping is wrong")
            if max(r) > k1:
                stats["maxcol_skip"] += 1
                continue
            if sum(t * t for t in r) < Emin:
                continue
            stats["E_pass"] += 1
            got = cover(rows, full, z)
            if got is not None:
                C = tuple(sorted(sp.neg[s] for s in got))
                if not check_star(sp, A, B, C):
                    raise AssertionError("set cover returned a non-solution")
                stats["solutions"] += 1
                if len(sols) < max_solutions:
                    sols.append({"A": list(A), "B": list(B), "C": list(C)})
    return stats, sols


def search_profile(sp, x, y, z, levels=None, log=None, max_solutions=3, stop_at_first=False,
                   workers=1, deadline=None):
    Emin, need = plan(sp, x, y, z)
    res = {"d": sp.d, "x": x, "y": y, "z": z, "Emin": Emin, "need": need,
           "per_k1": {}, "solutions": [], "n_solutions": 0}
    if not need:
        res["verdict"] = "CLOSED-by-bounds"
        return res
    maxk = max(need)
    if levels is None:
        levels = sp.orbit_levels(maxk, xbounds(need, maxk), log=log)
    tasks = []
    for k1 in sorted(need):
        reps = [I for (I, _) in levels[k1] if sp.excess(I) <= need[k1]]
        res["per_k1"][k1] = {"Xmax": need[k1], "reps": len(reps)}
        for I in reps:
            tasks.append((k1, I))
    results = []
    if workers > 1 and len(tasks) > 1:
        import multiprocessing as mp
        with mp.Pool(workers, initializer=_init_worker, initargs=(sp.d,)) as pool:
            for out in pool.imap_unordered(_work, [(x, y, z, k1, I, Emin, max_solutions, deadline)
                                                   for (k1, I) in tasks]):
                results.append(out)
                if log:
                    log("  k1=%d I=%s %s" % (out[0], list(out[1]), out[2]))
    else:
        for (k1, I) in tasks:
            st, so = search_rep(sp, x, y, z, k1, I, Emin, max_solutions, deadline)
            results.append((k1, I, st, so))
            if log:
                log("  k1=%d I=%s %s" % (k1, list(I), st))
            if stop_at_first and so:
                break
    for k1, I, st, so in results:
        agg = res["per_k1"][k1]
        for key, val in st.items():
            if isinstance(val, bool):
                agg[key] = agg.get(key, False) or val
            else:
                agg[key] = agg.get(key, 0) + val
        res["n_solutions"] += st["solutions"]
        for s in so:
            if len(res["solutions"]) < max_solutions:
                res["solutions"].append(s)
    timed_out = any(v.get("timeout") for v in res["per_k1"].values())
    if timed_out:
        res["verdict"] = "INCOMPLETE-timeout"
    elif res["n_solutions"] or res["solutions"]:
        res["verdict"] = "FEASIBLE"
    else:
        res["verdict"] = "CLOSED"
    return res


_W = {}


def _init_worker(d):
    _W["sp"] = Space(d)


def _work(args):
    x, y, z, k1, I, Emin, max_solutions, deadline = args
    st, so = search_rep(_W["sp"], x, y, z, k1, I, Emin, max_solutions, deadline)
    return (k1, I, st, so)


# ---- the run ----------------------------------------------------------------------
LEV27_CUT_6_8_13 = [  # Lev's 27-point Kakeya set in F_3^4 (kcutc control), last coordinate = level
    [1, 2, 0, 0], [1, 0, 1, 0], [0, 1, 1, 0], [1, 2, 1, 0], [0, 0, 2, 0], [2, 1, 2, 0],
    [0, 0, 0, 1], [1, 0, 0, 1], [2, 0, 0, 1], [0, 1, 0, 1], [2, 0, 1, 1], [2, 1, 1, 1],
    [1, 2, 1, 1], [0, 0, 2, 1],
    [0, 0, 0, 2], [2, 0, 0, 2], [1, 1, 0, 2], [2, 1, 0, 2], [0, 2, 0, 2], [2, 2, 0, 2],
    [0, 0, 1, 2], [1, 1, 1, 2], [0, 2, 1, 2], [0, 0, 2, 2], [2, 0, 2, 2], [0, 1, 2, 2],
    [1, 1, 2, 2]]


def lev_sections(sp3):
    secs = {0: [], 1: [], 2: []}
    for c in LEV27_CUT_6_8_13:
        secs[c[3]].append(sum(c[i] * 3 ** i for i in range(3)))
    return secs[0], secs[1], secs[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default="k5cut.json")
    ap.add_argument("--minutes", type=float, default=330.0)
    args = ap.parse_args()
    t_start = time.time()
    deadline = t_start + 60 * args.minutes
    out = {"contract": "k5cut-v1", "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

    def log(msg):
        print("[%8.1fs] %s" % (time.time() - t_start, msg), flush=True)

    def save():
        tmp = args.out + ".tmp"
        with open(tmp, "w") as f:
            json.dump(out, f, indent=1, default=str)
        os.replace(tmp, args.out)

    # 1. moment LP with SD and SD2 for n = 4..7 (cheap)
    lp = {}
    for n, lo, hi in ((4, 18, 27), (5, 36, 53), (6, 66, 100), (7, 120, 160)):
        lp[n] = {}
        for rule in ("t", "sd1", "sd2"):
            ex = sd2.scan(n, lo, hi, rule)
            lp[n][rule] = {"excluded": ex, "bound": (max(ex) + 1) if ex else None,
                           "range": [lo, hi]}
        log("LP n=%d %s" % (n, {r: lp[n][r]["bound"] for r in lp[n]}))
    out["lp"] = lp
    save()

    # 2. control: Lev's cut (6,8,13) in F_3^3 satisfies (*) and all pruning conditions
    sp3 = Space(3)
    A, B, C = lev_sections(sp3)
    ctrl = {"sizes": [len(A), len(B), len(C)], "star": check_star(sp3, A, B, C)}
    Emin3, need3 = plan(sp3, len(A), len(B), len(C))
    r, rows = rows_and_columns(sp3, A, B)
    ctrl["E"] = sum(t * t for t in r)
    ctrl["Emin"] = Emin3
    ctrl["maxcol"] = max(r)
    ctrl["maxcol_allowed"] = max(r) in need3
    v0 = r.index(max(r))
    Bt = [sp3.sub[b][v0] for b in B]
    I = sorted(set(A) & set(Bt))
    ctrl["I_excess"] = sp3.excess(I)
    ctrl["I_excess_ok"] = max(r) in need3 and sp3.excess(I) <= need3[max(r)]
    ctrl["min_cover"] = min_cover_size(sp3, A, B)
    out["control_lev"] = ctrl
    log("Lev control %s" % ctrl)
    if not (ctrl["star"] and ctrl["E"] >= Emin3 and ctrl["I_excess_ok"] and ctrl["min_cover"] <= 13):
        out["verdict"] = "CONTROL-FAILED"
        save()
        sys.exit(1)
    save()

    # 3. the main search: (10, 11, 23) in F_3^4
    sp4 = Space(4)
    x, y, z = 10, 11, 23
    Emin, need = plan(sp4, x, y, z)
    log("plan (10,11,23): Emin=%d need=%s" % (Emin, need))
    maxk = max(need)
    levels = sp4.orbit_levels(maxk, xbounds(need, maxk), log=log)
    out["orbits"] = {k: len(v) for k, v in levels.items()}
    # orbit-stabiliser sanity check on the levels that are complete (bound >= max excess)
    full_ok = {}
    for k, reps in levels.items():
        if xbounds(need, maxk)[k] >= comb(k, 2):
            full_ok[k] = (sum(s for _, s in reps) == comb(81, k))
    out["orbit_sum_checks"] = full_ok
    log("orbit sums complete levels: %s" % full_ok)
    if not all(full_ok.values()):
        out["verdict"] = "ORBIT-CHECK-FAILED"
        save()
        sys.exit(1)
    save()
    res = search_profile(sp4, x, y, z, levels=levels, log=log, workers=args.workers,
                         deadline=deadline)
    out["main"] = res
    out["verdict"] = res["verdict"]
    out["seconds"] = round(time.time() - t_start, 1)
    save()
    log("verdict: %s" % res["verdict"])


if __name__ == "__main__":
    main()
