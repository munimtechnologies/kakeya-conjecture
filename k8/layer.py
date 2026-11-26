# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Layered search for small Kakeya sets in F_3^n (standard library only).

Write a point of F_3^n as (c, x), c in F_3, x in F_3^m, m = n - 1, and fix a Kakeya set L0
of F_3^m in the layer c = 0.  L0 covers every direction (0, d).  A direction (1, v),
v in F_3^m, is covered iff some line {(0, x), (1, x + v), (2, x + 2v)} lies in the set, i.e.
iff some x in L0 has x + v in L1 and x + 2v in L2 (every line in such a direction meets each
layer once).  So

    K = {0} x L0  u  {1} x L1  u  {2} x L2

is a Kakeya set of F_3^n iff the pairs cover all 3^m vectors v.  This module searches for
small L1 u L2 for a given L0: 3^m conditions with |L0| candidate lines each, instead of the
(3^n - 1)/2 directions and 3^(n-1) lines per direction of the full problem.

In the 211-point Kakeya set of F_3^7 (k6/sets/kakeya7_211.json) the layer L0 is the
105-point set of F_3^6 and |L1| + |L2| = 106 = 2 * 53.  Searches with a frozen 53-point
layer in F_3^6 stopped at 53 + 54 = 107.  The question here: is the extra cost about
2 k_(n-2) in general, and can it be beaten?

Point encoding: an integer i in [0, 3^m) is x = digits of i in base 3, first coordinate most
significant.  Internal points p in [0, 2 * 3^m): p < N is y = p in L1, p >= N is z = p - N in
L2; in F_3^n the point p has index N + p (and x in L0 has index x).

Search: WalkSAT style as k6/fastwalk.py.  For every point we keep cost[p] = the number of
vectors v whose only covering line uses p.  Because all lines of one v use different points
x of L0, the index sum xs[v] of the covering lines identifies the unique line when cnt[v]=1.
"""
import array
import itertools
import os
import random
import sys
import time
import traceback

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:  # so that k5 and k8 import in spawned workers too
    sys.path.insert(0, _ROOT)


# ---- arithmetic in F_3^m -----------------------------------------------------------------
def _tab(k):
    base = 3 ** k
    digs = list(itertools.product(range(3), repeat=k))
    t = [0] * (base * base)
    for x, dx in enumerate(digs):
        for y, dy in enumerate(digs):
            z = 0
            for a, b in zip(dx, dy):
                z = 3 * z + (a + b) % 3
            t[x * base + y] = z
    return t


def arith(m):
    """(add, neg) on integers 0..3^m-1 encoding F_3^m."""
    s = m // 2
    lo, hi = 3 ** s, 3 ** (m - s)
    AH, AL = _tab(m - s), _tab(s)

    def add(x, y):
        xh, xl = divmod(x, lo)
        yh, yl = divmod(y, lo)
        return AH[xh * hi + yh] * lo + AL[xl * lo + yl]

    N = 3 ** m
    NEG = [0] * N
    for i, d in enumerate(itertools.product(range(3), repeat=m)):
        z = 0
        for a in d:
            z = 3 * z + (-a) % 3
        NEG[i] = z
    return add, NEG


class Layer:
    """Incremental state for a fixed L0 in F_3^m."""

    def __init__(self, m, L0):
        self.m = m
        self.N = N = 3 ** m
        self.X = X = sorted(set(L0))
        self.L = L = len(X)
        add, NEG = arith(m)
        YV = array.array("H", [0]) * (N * L)   # YV[v*L+i] = x_i + v   (point of L1)
        ZV = array.array("H", [0]) * (N * L)   # ZV[v*L+i] = x_i - v = x_i + 2v (point of L2)
        for v in range(N):
            nv, b = NEG[v], v * L
            for i, x in enumerate(X):
                YV[b + i] = add(x, v)
                ZV[b + i] = add(x, nv)
        V1 = array.array("H", [0]) * (N * L)   # V1[y*L+i] = v with x_i + v = y
        V2 = array.array("H", [0]) * (N * L)   # V2[z*L+i] = v with x_i + 2v = z
        for v in range(N):
            b = v * L
            for i in range(L):
                V1[YV[b + i] * L + i] = v
                V2[ZV[b + i] * L + i] = v
        self.YV, self.ZV, self.V1, self.V2 = YV, ZV, V1, V2
        self.reset()

    # -- state --
    def reset(self):
        N = self.N
        self.inS = bytearray(2 * N)
        self.cnt = [0] * N
        self.xs = [0] * N
        self.cost = [0] * (2 * N)
        self.unc = list(range(N))
        self.upos = list(range(N))
        self.S = []
        self.spos = [-1] * (2 * N)

    def _unc_del(self, v):
        i = self.upos[v]
        last = self.unc.pop()
        if last != v:
            self.unc[i] = last
            self.upos[last] = i
        self.upos[v] = -1

    def _unc_add(self, v):
        self.upos[v] = len(self.unc)
        self.unc.append(v)

    def _inc(self, v, i, py, pz):
        c = self.cnt[v]
        if c == 0:
            self._unc_del(v)
            self.cost[py] += 1
            self.cost[pz] += 1
        elif c == 1:
            i0, b = self.xs[v], v * self.L
            self.cost[self.YV[b + i0]] -= 1
            self.cost[self.N + self.ZV[b + i0]] -= 1
        self.cnt[v] = c + 1
        self.xs[v] += i

    def _dec(self, v, i, py, pz):
        c = self.cnt[v] - 1
        self.cnt[v] = c
        self.xs[v] -= i
        if c == 0:
            self._unc_add(v)
            self.cost[py] -= 1
            self.cost[pz] -= 1
        elif c == 1:
            i0, b = self.xs[v], v * self.L
            self.cost[self.YV[b + i0]] += 1
            self.cost[self.N + self.ZV[b + i0]] += 1

    def _pairs(self, p):
        """(v, i, py, pz) for every line through p whose other free point is in the set."""
        N, L, inS = self.N, self.L, self.inS
        out = []
        if p < N:
            b = p * L
            V1, ZV = self.V1, self.ZV
            for i in range(L):
                v = V1[b + i]
                pz = N + ZV[v * L + i]
                if inS[pz]:
                    out.append((v, i, p, pz))
        else:
            b = (p - N) * L
            V2, YV = self.V2, self.YV
            for i in range(L):
                v = V2[b + i]
                py = YV[v * L + i]
                if inS[py]:
                    out.append((v, i, py, p))
        return out

    def add(self, p):
        if self.inS[p]:
            return
        for t in self._pairs(p):
            self._inc(*t)
        self.inS[p] = 1
        self.spos[p] = len(self.S)
        self.S.append(p)

    def remove(self, p):
        if not self.inS[p]:
            return
        self.inS[p] = 0
        for t in self._pairs(p):
            self._dec(*t)
        i = self.spos[p]
        last = self.S.pop()
        if last != p:
            self.S[i] = last
            self.spos[last] = i
        self.spos[p] = -1

    def load(self, pts):
        self.reset()
        for p in pts:
            self.add(p)

    def line_missing(self, v, i):
        b = v * self.L + i
        py, pz = self.YV[b], self.N + self.ZV[b]
        return [q for q in (py, pz) if not self.inS[q]]

    def greedy(self, rng):
        L = self.L
        while self.unc:
            v = rng.choice(self.unc)
            best, opts = 3, []
            for i in range(L):
                k = len(self.line_missing(v, i))
                if k < best:
                    best, opts = k, [i]
                elif k == best:
                    opts.append(i)
            for q in self.line_missing(v, rng.choice(opts)):
                self.add(q)

    def prune(self, rng):
        while True:
            free = [p for p in self.S if self.cost[p] == 0]
            if not free:
                return
            self.remove(rng.choice(free))

    def full_indices(self, S=None):
        """Indices in F_3^n (n = m + 1, first coordinate most significant)."""
        S = self.S if S is None else S
        return sorted(list(self.X) + [self.N + p for p in S])


# ---- search ------------------------------------------------------------------------------
def search(W, seed=1, budget=60.0, noise=0.1, tenure=5, stall=5000, ncand=8, start=None,
           max_steps=None, stop_at=0):
    """Returns {'best': |L1|+|L2|, 'points': internal points, 'reached': [(extra, s)], ...}."""
    t0 = time.time()
    rng = random.Random(seed)
    W.reset()
    if start:
        for p in start:
            W.add(p)
    W.greedy(rng)
    W.prune(rng)
    best = sorted(W.S)
    reached = [(len(best), round(time.time() - t0, 1))]
    traj = list(best)
    P, L = 2 * W.N, W.L
    added_at = [-10 ** 9] * P
    st = {"it": 0, "target": len(W.S) - 1}

    def shrink(protect):
        it = st["it"]
        cost = W.cost
        while len(W.S) > st["target"]:
            cand = [x for x in W.S if x not in protect and it - added_at[x] > tenure]
            if not cand:
                cand = [x for x in W.S if x not in protect]
                if not cand:
                    return
            if rng.random() < noise:
                x = rng.choice(cand)
            else:
                mc = min(cost[x] for x in cand)
                x = rng.choice([x for x in cand if cost[x] == mc])
            W.remove(x)

    shrink(())
    last, bestu, restarts = 0, 10 ** 9, 0
    while time.time() - t0 < budget and (max_steps is None or st["it"] < max_steps):
        st["it"] += 1
        it = st["it"]
        u = len(W.unc)
        if u == 0:
            if len(W.S) < len(traj):
                traj = sorted(W.S)
            if len(W.S) < len(best):
                best = sorted(W.S)
                reached.append((len(best), round(time.time() - t0, 1)))
                if len(best) <= stop_at:
                    break
            st["target"] = len(W.S) - 1
            shrink(())
            last, bestu = it, len(W.unc)
            continue
        if u < bestu:
            bestu, last = u, it
        if it - last > stall:
            restarts += 1
            W.load(traj)
            if W.S:
                W.remove(rng.choice(W.S))
            st["target"] = len(W.S)
            for i in range(P):
                added_at[i] = -10 ** 9
            last, bestu = it, len(W.unc)
            continue
        v = rng.choice(W.unc)
        inS, YV, ZV, N = W.inS, W.YV, W.ZV, W.N
        b = v * L
        mn, opts, near = 3, [], []
        for i in range(L):
            k = (not inS[YV[b + i]]) + (not inS[N + ZV[b + i]])
            if k < mn:
                near = opts if k == mn - 1 else []
                mn, opts = k, [i]
            elif k == mn:
                opts.append(i)
            elif k == mn + 1:
                near.append(i)
        if rng.random() < noise:
            i = rng.choice(opts + near)
        else:
            sample = opts if len(opts) <= ncand else rng.sample(opts, ncand)
            bu, bl = None, []
            for i in sample:
                miss = W.line_missing(v, i)
                for q in miss:
                    W.add(q)
                uu = len(W.unc)
                for q in miss:
                    W.remove(q)
                if bu is None or uu < bu:
                    bu, bl = uu, [i]
                elif uu == bu:
                    bl.append(i)
            i = rng.choice(bl)
        pick = W.line_missing(v, i)
        for q in pick:
            W.add(q)
            added_at[q] = it
        shrink(set(pick))
    return {"best": len(best), "points": best, "reached": reached, "steps": st["it"],
            "restarts": restarts, "seconds": round(time.time() - t0, 1)}


# ---- inputs ------------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def load_repo_set(name):
    import json
    with open(os.path.join(ROOT, "k6", "sets", name), encoding="utf-8") as f:
        return json.load(f)["points"]


def classical(m):
    S = set()
    for pat in itertools.product((0, 1), repeat=m):
        for v in (1, 2):
            S.add(sum(v * c * 3 ** (m - 1 - k) for k, c in enumerate(pat)))
    return sorted(S)


K7_PLANE = [0, 1, 2, 3, 5, 6, 8]  # a 7-point Kakeya set of F_3^2


def downloaded(m, size):
    """Published set of Chung, Du and Wesley in F_3^m (fetched at run time, not copied)."""
    from k5 import cdw
    pts, note = cdw.load_set(m, size)
    return sorted(sum(c * 3 ** (m - 1 - k) for k, c in enumerate(p)) for p in pts), note


# ---- parallel probe ----------------------------------------------------------------------
def task(t):
    label, m, L0, seed, budget, noise = t
    t0 = time.time()
    try:
        from k8 import check
        W = Layer(m, L0)
        tset = round(time.time() - t0, 1)
        r = search(W, seed=seed, budget=max(5.0, budget - (time.time() - t0)), noise=noise)
        full = W.full_indices(r["points"])
        r["verified"] = check.is_kakeya(m + 1, full)
        r["full"] = full
        r.update(label=label, m=m, L0=len(W.X), seed=seed, noise=noise, setup_s=tset)
        return r
    except Exception:  # noqa: BLE001 - the report is the product
        return {"label": label, "error": traceback.format_exc()[-600:]}


def run_pool(tasks, workers):
    try:
        import multiprocessing as mp
        ctx = mp.get_context("spawn")
        with ctx.Pool(workers) as p:
            return p.map(task, tasks, chunksize=1), "pool%d" % workers
    except Exception as e:  # noqa: BLE001
        return [task(t) for t in tasks[:1]], "serial(%s:%s)" % (type(e).__name__, str(e)[:80])


RECORD = {4: 27, 5: 53, 6: 105, 7: 211, 8: 511}   # 8: the classical construction 2^9 - 1
PATTERN = {4: 27, 5: 53, 6: 107, 7: 213, 8: 427}


def probe(seconds=450, workers=5):
    t0 = time.time()
    out = []
    try:
        L211 = load_repo_set("kakeya7_211.json")
        L211_0 = [i for i in L211 if i < 3 ** 6]       # its layer x0 = 0, the 105-point set
        L105 = load_repo_set("kakeya6_105.json")
        notes = ["L211 %d L105 %d layer0(211)==105:%s" % (len(L211), len(L105),
                                                          sorted(L211_0) == sorted(L105))]
        sets = {}
        for m, size in ((4, 27), (5, 53)):
            try:
                sets[m], note = downloaded(m, size)
                notes.append("d%d_%d ok %d" % (m, size, len(sets[m])))
            except Exception as e:  # noqa: BLE001
                notes.append("d%d_%d failed %s" % (m, size, str(e)[:60]))
        out.append("inputs " + "; ".join(notes))
        budget = seconds - (time.time() - t0) - 30
        tasks = [("n8-L211", 7, L211, 1, budget, 0.10),
                 ("n8-L211", 7, L211, 2, budget, 0.20),
                 ("n7-L105", 6, L105, 3, budget, 0.10)]
        tasks.append(("n6-L53", 5, sets[5], 4, budget, 0.10) if 5 in sets
                     else ("n6-Lcl63", 5, classical(5), 4, budget, 0.10))
        tasks.append(("n5-L27", 4, sets[4], 5, budget, 0.10) if 4 in sets
                     else ("n5-Lcl31", 4, classical(4), 5, budget, 0.10))
        res, mode = run_pool(tasks[:workers], workers)
        out.append("mode=%s cpu=%s budget=%.0fs" % (mode, os.cpu_count(), budget))
        news = []
        for r in res:
            if "error" in r:
                out.append("%s ERROR %s" % (r["label"], r["error"]))
                continue
            n = r["m"] + 1
            tot = r["L0"] + r["best"]
            rr = r["reached"]
            short = rr if len(rr) <= 7 else rr[:2] + ["..."] + rr[-4:]
            out.append("%s |L0|=%d extra=%d total=%d (record %d, pattern %d) verified=%s "
                       "setup=%ss reached %s; %d steps, %d restarts" % (
                           r["label"], r["L0"], r["best"], tot, RECORD[n], PATTERN[n],
                           r["verified"], r["setup_s"], short, r["steps"], r["restarts"]))
            if r["verified"] and tot < RECORD[n]:
                news.append((n, tot, r))
        best_by_n = {}
        for n, tot, r in news:
            if n not in best_by_n or tot < best_by_n[n][0]:
                best_by_n[n] = (tot, r)
        for n in sorted(best_by_n):
            tot, r = best_by_n[n]
            out.append("NEW n=%d size=%d (record was %d) label=%s POINTS %s" % (
                n, tot, RECORD[n], r["label"], ",".join(map(str, r["full"]))))
        if not news:
            out.append("no set below the records")
    except Exception:  # noqa: BLE001
        out.append("ERROR " + traceback.format_exc()[-700:])
    out.append("total %.0fs" % (time.time() - t0))
    return " || ".join(out)
