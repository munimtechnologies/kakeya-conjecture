# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Shared code for AAB pairs over F_3 (see aab/AAB.md).

An AAB pair of F_3^m is (A, B) with A n B non-empty such that every direction d has a line
{p, p + d, p + 2d} with p, p + d in A and p + 2d in B.  Equivalently: for every direction
there is a pair {x, y} of points of A in that direction whose third collinear point
-(x + y) (its "centre") lies in B.

Points of F_3^m are indexed by sum x_k 3^(m-1-k) (first coordinate most significant), the
order of itertools.product(range(3), repeat=m)."""
import itertools
import os
import random
import sys
import time

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)


def points(m):
    return list(itertools.product(range(3), repeat=m))


def norm(v):
    for c in v:
        if c:
            return tuple(v) if c == 1 else tuple((2 * x) % 3 for x in v)
    return None


def to_idx(p):
    i = 0
    for c in p:
        i = 3 * i + int(c) % 3
    return i


def to_pt(i, m):
    return tuple((i // 3 ** (m - 1 - k)) % 3 for k in range(m))


class Space:
    """PD[i][j]: index (in D) of the direction of j - i.  PC[i][j]: index of -(i + j), the
    third point of the line through i and j."""

    def __init__(self, m):
        self.m = m
        self.P = points(m)
        self.N = len(self.P)
        self.D = [d for d in self.P if norm(d) == d]
        self.ND = len(self.D)
        didx = {d: i for i, d in enumerate(self.D)}
        self.PD = [[-1] * self.N for _ in range(self.N)]
        self.PC = [[-1] * self.N for _ in range(self.N)]
        for i, p in enumerate(self.P):
            rowd, rowc = self.PD[i], self.PC[i]
            for j, q in enumerate(self.P):
                if i != j:
                    rowd[j] = didx[norm(tuple((b - a) % 3 for a, b in zip(p, q)))]
                    rowc[j] = to_idx(tuple((-a - b) % 3 for a, b in zip(p, q)))


# ---------------------------------------------------------------- independent checkers

def check_pair(m, A, B):
    """From the definition, with coordinate tuples: A n B non-empty and every direction d
    has p with p, p + d in A and p + 2d in B."""
    P = points(m)
    At = set(P[i] for i in A)
    Bt = set(P[i] for i in B)
    if not (At & Bt):
        return False
    for d in P:
        if norm(d) != d:
            continue
        ok = False
        for p in At:
            q = tuple((a + b) % 3 for a, b in zip(p, d))
            r = tuple((a + 2 * b) % 3 for a, b in zip(p, d))
            if q in At and r in Bt:
                ok = True
                break
        if not ok:
            return False
    return True


def kakeya_by_pairs(n, idx):
    """Every direction of F_3^n is the direction of a line inside the set."""
    P = sorted(set(to_pt(i, n) for i in idx))
    S = set(P)
    cov = set()
    for i, x in enumerate(P):
        for y in P[i + 1:]:
            if tuple((-a - b) % 3 for a, b in zip(x, y)) in S:
                cov.add(norm(tuple((b - a) % 3 for a, b in zip(x, y))))
    return len(cov) == (3 ** n - 1) // 2


def product(m1, A1, B1, m2, A2, B2):
    N2 = 3 ** m2
    return (m1 + m2, [x * N2 + y for x in A1 for y in A2],
            [x * N2 + y for x in B1 for y in B2])


# ---------------------------------------------------------------- exact search, given A

def centres(S, A):
    """cand[d]: bitmask of the centres -(x + y) of the pairs {x, y} of A in direction d;
    mask: bitmask of the directions of the pairs of A."""
    A = list(A)
    cand = [0] * S.ND
    mask = 0
    for x in range(len(A)):
        px, cx = S.PD[A[x]], S.PC[A[x]]
        for y in range(x + 1, len(A)):
            d = px[A[y]]
            mask |= 1 << d
            cand[d] |= 1 << cx[A[y]]
    return cand, mask


def min_out(cand, amask, E, limit):
    """Least f such that f centres outside A (amask) together with at most E centres inside
    A meet every cand[d].  Exact branch and bound.  Returns (f, chosen-bitmask) with
    f <= limit, or (None, None)."""
    nd = len(cand)
    cov = {}
    for d in range(nd):
        c = cand[d]
        while c:
            low = c & -c
            j = low.bit_length() - 1
            cov[j] = cov.get(j, 0) | (1 << d)
            c ^= low
    vals = list(cov.values())
    best = [limit + 1, None]

    def rec(unc, chosen, f, ein):
        if unc == 0:
            if f < best[0]:
                best[0], best[1] = f, chosen
            return
        if f >= best[0]:
            return
        cs = sorted(((v & unc).bit_count() for v in vals), reverse=True)
        if not cs or cs[0] == 0:
            return
        left, k = unc.bit_count(), 0
        for c in cs:
            if left <= 0 or c == 0:
                break
            left -= c
            k += 1
        if left > 0:
            return
        if f + max(0, k - (E - ein)) >= best[0]:
            return
        dsel, nb, u = -1, 1 << 30, unc
        while u:
            low = u & -u
            d = low.bit_length() - 1
            u ^= low
            nc = cand[d].bit_count()
            if nc < nb:
                nb, dsel = nc, d
        c = cand[dsel]
        while c:
            low = c & -c
            c ^= low
            j = low.bit_length() - 1
            if amask & low:
                if ein < E:
                    rec(unc & ~cov[j], chosen | low, f, ein + 1)
            else:
                rec(unc & ~cov[j], chosen | low, f + 1, ein)

    rec((1 << nd) - 1, 0, 0, 0)
    if best[1] is None:
        return None, None
    return best[0], best[1]


def analyse(S, A, E, limit):
    """For a fixed A: the least B (at most E points of B inside A, at most `limit` outside)
    making (A, B) an AAB pair, padded with min(A) if B misses A.  None if there is none."""
    cand, mask = centres(S, A)
    if mask != (1 << S.ND) - 1:
        return None
    amask = 0
    for x in A:
        amask |= 1 << x
    f, ch = min_out(cand, amask, E, limit)
    if f is None:
        return None
    B = [j for j in range(S.N) if ch >> j & 1]
    if not (ch & amask):
        B = sorted(B + [min(A)])
    return f, B


def scan3(a, share, nshares, Es, limit, deadline=None):
    """Exhaustive over the a-point sets A of F_3^3 with 0 in A (every pair can be translated
    so), taking every nshares-th one from `share`.  For each E in Es: the least number f of
    points of B outside A (f <= limit) over all A, with at most E points of B inside A."""
    S = Space(3)
    full = (1 << S.ND) - 1
    PD, PC = S.PD, S.PC
    res = {E: {"best": None, "count": 0, "example": None} for E in Es}
    seen, done, k = 0, True, 0
    for i, rest in enumerate(itertools.combinations(range(1, S.N), a - 1)):
        if i % nshares != share:
            continue
        k += 1
        if deadline is not None and (k & 1023) == 0 and time.time() > deadline:
            done = False
            break
        A = (0,) + rest
        cand = [0] * S.ND
        mask = 0
        for x in range(a):
            px, cx = PD[A[x]], PC[A[x]]
            for y in range(x + 1, a):
                d = px[A[y]]
                mask |= 1 << d
                cand[d] |= 1 << cx[A[y]]
        if mask != full:
            continue
        seen += 1
        amask = 0
        for x in A:
            amask |= 1 << x
        for E in Es:
            r = res[E]
            lim = limit if r["best"] is None else r["best"]
            f, ch = min_out(cand, amask, E, lim)
            if f is None:
                continue
            if r["best"] is None or f < r["best"]:
                B = [j for j in range(S.N) if ch >> j & 1]
                if not (ch & amask):
                    B = sorted(B + [0])
                r["best"], r["count"], r["example"] = f, 0, (list(A), B)
            r["count"] += 1
    return {"a": a, "share": share, "nshares": nshares, "limit": limit, "seen": seen,
            "done": done, "res": {str(E): res[E] for E in Es}}


def scan_task(t):
    return scan3(*t)


def run_pool(fn, tasks, workers=8):
    import multiprocessing as mp
    pool = mp.get_context("spawn").Pool(max(1, min(workers, len(tasks))))
    try:
        return pool.map(fn, tasks, chunksize=1)
    finally:
        pool.close()
        pool.join()


# ---------------------------------------------------------------- searches inside a set K

def lines_in(S, K):
    """Lines inside K as (direction index, (x, y, z)) with x < y < z."""
    Ks = set(K)
    Kl = sorted(Ks)
    out = []
    for i, x in enumerate(Kl):
        px, cx = S.PD[x], S.PC[x]
        for y in Kl[i + 1:]:
            z = cx[y]
            if z > y and z in Ks:
                out.append((px[y], (x, y, z)))
    return out


def label_search(S, K, bound, deadline, t=2):
    """Exact branch and bound for an AAB pair (A, B) with A, B inside K minimising
    a^t + b^t - e^t below `bound` (the cost only grows when a point joins A or B)."""
    Kl = sorted(set(K))
    pos = {x: i for i, x in enumerate(Kl)}
    opts = [[] for _ in range(S.ND)]
    for d, (x, y, z) in lines_in(S, Kl):
        bx, by, bz = 1 << pos[x], 1 << pos[y], 1 << pos[z]
        opts[d].extend([(by | bz, bx), (bx | bz, by), (bx | by, bz)])
    if not all(opts):
        return {"kakeya": False}
    order = sorted(range(S.ND), key=lambda d: len(opts[d]))

    def cost(a, b, e):
        return a ** t + b ** t - e ** t

    st = {"best": bound, "am": None, "bm": None, "nodes": 0, "timeout": False}

    def rec(am, bm):
        st["nodes"] += 1
        if (st["nodes"] & 511) == 0 and time.time() > deadline:
            st["timeout"] = True
        if st["timeout"]:
            return
        a, b, e = am.bit_count(), bm.bit_count(), (am & bm).bit_count()
        if cost(a, b, e) >= st["best"]:
            return
        dsel = -1
        for d in order:
            for an, bn in opts[d]:
                if an & am == an and bn & bm == bn:
                    break
            else:
                dsel = d
                break
        if dsel < 0:
            if e == 0:
                c1, c2 = cost(a, b + 1, 1), cost(a + 1, b, 1)
                if min(c1, c2) < st["best"]:
                    if c1 <= c2:
                        bm |= am & -am
                    else:
                        am |= bm & -bm
                    st["best"], st["am"], st["bm"] = min(c1, c2), am, bm
            else:
                st["best"], st["am"], st["bm"] = cost(a, b, e), am, bm
            return
        ch = []
        for an, bn in opts[dsel]:
            na, nb = am | an, bm | bn
            ch.append((cost(na.bit_count(), nb.bit_count(), (na & nb).bit_count()), na, nb))
        ch.sort()
        for c, na, nb in ch:
            if c >= st["best"]:
                break
            rec(na, nb)

    rec(0, 0)
    out = {"kakeya": True, "complete": not st["timeout"], "nodes": st["nodes"],
           "best": None}
    if st["am"] is not None:
        A = [Kl[i] for i in range(len(Kl)) if st["am"] >> i & 1]
        B = [Kl[i] for i in range(len(Kl)) if st["bm"] >> i & 1]
        out.update(best=st["best"], A=A, B=B)
    return out


class KWalk:
    """WalkSAT-style local search for an AAB pair inside K with |A| <= ta, |B| <= tb."""

    def __init__(self, S, K, ta, tb, seed, noise=0.1):
        self.rng = random.Random(seed)
        self.K = sorted(set(K))
        n = len(self.K)
        pos = {x: i for i, x in enumerate(self.K)}
        self.ND = S.ND
        self.opts = []
        for d, (x, y, z) in lines_in(S, self.K):
            px, py, pz = pos[x], pos[y], pos[z]
            self.opts.extend([(d, py, pz, px), (d, px, pz, py), (d, px, py, pz)])
        self.byd = [[] for _ in range(S.ND)]
        self.oA = [[] for _ in range(n)]
        self.oB = [[] for _ in range(n)]
        for o, (d, a1, a2, b) in enumerate(self.opts):
            self.byd[d].append(o)
            self.oA[a1].append(o)
            self.oA[a2].append(o)
            self.oB[b].append(o)
        self.n, self.ta, self.tb, self.noise = n, ta, tb, noise
        self.inA = [False] * n
        self.inB = [False] * n
        self.sat = [False] * len(self.opts)
        self.cnt = [0] * S.ND
        self.na = 0
        self.nb = 0
        idx = list(range(n))
        self.rng.shuffle(idx)
        for x in idx[:ta]:
            self.setA(x, True)
        self.rng.shuffle(idx)
        for x in idx[:tb]:
            self.setB(x, True)

    def _upd(self, olist):
        inA, inB, sat, cnt, opts = self.inA, self.inB, self.sat, self.cnt, self.opts
        for o in olist:
            d, a1, a2, b = opts[o]
            s = inA[a1] and inA[a2] and inB[b]
            if s != sat[o]:
                sat[o] = s
                cnt[d] += 1 if s else -1

    def setA(self, x, v):
        if self.inA[x] != v:
            self.inA[x] = v
            self.na += 1 if v else -1
            self._upd(self.oA[x])

    def setB(self, x, v):
        if self.inB[x] != v:
            self.inB[x] = v
            self.nb += 1 if v else -1
            self._upd(self.oB[x])

    def brk(self, olist):
        loss = {}
        for o in olist:
            if self.sat[o]:
                d = self.opts[o][0]
                loss[d] = loss.get(d, 0) + 1
        return sum(1 for d, l in loss.items() if l == self.cnt[d])

    def run(self, deadline, maxsteps):
        rng = self.rng
        if not all(self.byd):
            return {"unc": -1, "steps": 0}
        best = self.ND + 1
        for step in range(maxsteps):
            if (step & 255) == 0 and time.time() > deadline:
                break
            unc = [d for d in range(self.ND) if self.cnt[d] == 0]
            if len(unc) < best:
                best = len(unc)
            if not unc:
                A = [self.K[i] for i in range(self.n) if self.inA[i]]
                B = [self.K[i] for i in range(self.n) if self.inB[i]]
                if not (set(A) & set(B)):
                    B = sorted(B + [A[0]])
                return {"unc": 0, "A": A, "B": B, "steps": step}
            d = rng.choice(unc)
            olist = self.byd[d]
            if rng.random() < self.noise:
                o = rng.choice(olist)
            else:
                sc = []
                for o in olist:
                    _, a1, a2, b = self.opts[o]
                    sc.append(((not self.inA[a1]) + (not self.inA[a2]) + (not self.inB[b]),
                               rng.random(), o))
                o = min(sc)[2]
            _, a1, a2, b = self.opts[o]
            self.setA(a1, True)
            self.setA(a2, True)
            self.setB(b, True)
            keep = (a1, a2, b)
            while self.na > self.ta:
                c = [x for x in range(self.n) if self.inA[x] and x not in keep]
                if not c:
                    break
                if rng.random() < self.noise:
                    x = rng.choice(c)
                else:
                    x = min((self.brk(self.oA[y]), rng.random(), y) for y in c)[2]
                self.setA(x, False)
            while self.nb > self.tb:
                c = [x for x in range(self.n) if self.inB[x] and x not in keep]
                if not c:
                    break
                if rng.random() < self.noise:
                    x = rng.choice(c)
                else:
                    x = min((self.brk(self.oB[y]), rng.random(), y) for y in c)[2]
                self.setB(x, False)
        return {"unc": best, "steps": maxsteps}
