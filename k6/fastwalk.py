# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Fast local search for small Kakeya sets in F_3^n, n = 6 and 7 (standard library only).

Known sizes of the smallest Kakeya set in F_3^n: k_1 = 3, k_2 = 7, k_3 = 13 (Lev 2009),
k_4 = 27 (upper bound Lev/Polymath; lower bound by our computer proof, PROOF26.md), and
45 <= k_5 <= 53 (53-point set of Chung, Du and Wesley, arXiv:2608.23691; 45 is k5/QUAD45.md).
V. F. Lev guessed in 2009 that k_n = k_(n-1) + 2 k_(n-2) continues, which predicts
k_6 = 107 and k_7 = 213.  The best published constructions we found for n = 6, 7 are the
classical set {0,1}^n u {0,2}^n with 2^(n+1) - 1 points (127 and 255).

This module searches for smaller sets.  Compared with k5/walk.py it keeps, for every point p
of the set, cost[p] = the number of directions whose only line inside the set passes
through p, updated incrementally in O(number of directions) per added or removed point.
(Two lines in the same direction are parallel, hence disjoint, so a point can only be the
last support of a direction whose line count is 1.)

Search (WalkSAT style at a size budget m): pick an uncovered direction at random, complete
one of its lines that needs the fewest new points (best of a few by how many uncovered
directions it covers), then drop points of least cost (no recently added ones) until the
size is m again.  When every direction is covered the set is recorded and m drops by one.
Points listed as `frozen` are never removed.  Every set reported is re-checked by
k5.walk.is_kakeya, which uses none of these tables.
"""
import array
import itertools
import os
import random
import sys
import time
import traceback

try:
    from k5 import walk as k5walk
except ImportError:  # pragma: no cover - run from inside k6/
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from k5 import walk as k5walk

LEV = {1: 3, 2: 7, 3: 13, 4: 27, 5: 53, 6: 107, 7: 213}
CLASSICAL = {n: 2 ** (n + 1) - 1 for n in range(1, 9)}


# ---- tables ------------------------------------------------------------------------------
def _addtab(k):
    base = 3 ** k
    digs = list(itertools.product(range(3), repeat=k))
    t = array.array("H", [0]) * (base * base)
    for x, dx in enumerate(digs):
        for y, dy in enumerate(digs):
            z = 0
            for a, b in zip(dx, dy):
                z = 3 * z + (a + b) % 3
            t[x * base + y] = z
    return t


def setup(n):
    """Points are indexed row-major in base 3 (first coordinate most significant), as in
    k5/walk.py.  Directions: nonzero vectors whose first nonzero entry is 1.
    LA[p*D + j] = p + d_j, LB[p*D + j] = p + 2 d_j, anchors[j] = smallest point of each
    line in direction j."""
    P = 3 ** n
    pts = list(itertools.product(range(3), repeat=n))
    dirs = [i for i, v in enumerate(pts) if any(v) and [c for c in v if c][0] == 1]
    D = len(dirs)
    s = n // 2
    m, h = 3 ** s, 3 ** (n - s)
    AH, AL = _addtab(n - s), _addtab(s)
    dsplit = [divmod(dj, m) for dj in dirs]
    LA = array.array("H", [0]) * (P * D)
    for i in range(P):
        hi, lo = divmod(i, m)
        bh, bl, b = hi * h, lo * m, i * D
        for j, (dh, dl) in enumerate(dsplit):
            LA[b + j] = AH[bh + dh] * m + AL[bl + dl]
    LB = array.array("H", [0]) * (P * D)
    for i in range(P):
        b = i * D
        for j in range(D):
            LB[b + j] = LA[LA[b + j] * D + j]
    anchors = []
    for j in range(D):
        anchors.append(array.array("H", [i for i in range(P)
                                         if i < LA[i * D + j] and i < LB[i * D + j]]))
    return {"n": n, "P": P, "D": D, "dirs": dirs, "pts": pts, "LA": LA, "LB": LB,
            "anchors": anchors}


def is_kakeya_idx(T, K):
    inK = bytearray(T["P"])
    for p in K:
        inK[p] = 1
    D, LA, LB = T["D"], T["LA"], T["LB"]
    for j in range(D):
        if not any(inK[x] and inK[LA[x * D + j]] and inK[LB[x * D + j]] for x in K):
            return False
    return True


def independent_check(n, K):
    """k5.walk.is_kakeya on coordinate tuples (no shared tables)."""
    pts = [tuple((p // 3 ** (n - 1 - k)) % 3 for k in range(n)) for p in K]
    return len(set(pts)) == len(K) and k5walk.is_kakeya(n, pts)


# ---- constructions -----------------------------------------------------------------------
def classical(n):
    """{0,1}^n u {0,2}^n, 2^(n+1) - 1 points: the line through (d_i^2) in direction d."""
    S = set()
    for pat in itertools.product((0, 1), repeat=n):
        for v in (1, 2):
            S.add(sum(v * c * 3 ** (n - 1 - k) for k, c in enumerate(pat)))
    return sorted(S)


def product(K1, n1, K2, n2):
    """Cartesian product (first factor = leading coordinates); a product of Kakeya sets is
    Kakeya."""
    return sorted(a * 3 ** n2 + b for a in K1 for b in K2)


def full(n):
    return list(range(3 ** n))


# ---- incremental state -------------------------------------------------------------------
class FW:
    def __init__(self, T, frozen=()):
        self.T = T
        self.P, self.D = T["P"], T["D"]
        self.LA, self.LB = T["LA"], T["LB"]
        self.fz = bytearray(self.P)
        for p in frozen:
            self.fz[p] = 1
        self.frozen = sorted(set(frozen))
        self.reset()

    def reset(self):
        P, D = self.P, self.D
        self.inK = bytearray(P)
        self.lines = [set() for _ in range(D)]
        self.cost = [0] * P
        self.unc = list(range(D))
        self.upos = list(range(D))
        self.K = []
        self.kpos = [-1] * P
        for p in self.frozen:
            self.add(p)

    def _unc_del(self, j):
        i = self.upos[j]
        last = self.unc.pop()
        if last != j:
            self.unc[i] = last
            self.upos[last] = i
        self.upos[j] = -1

    def _unc_add(self, j):
        self.upos[j] = len(self.unc)
        self.unc.append(j)

    def add(self, q):
        inK = self.inK
        if inK[q]:
            return
        D, LA, LB, cost, lines = self.D, self.LA, self.LB, self.cost, self.lines
        b = q * D
        for j in range(D):
            a = LA[b + j]
            if inK[a]:
                c = LB[b + j]
                if inK[c]:
                    S = lines[j]
                    k = len(S)
                    if k == 0:
                        self._unc_del(j)
                        cost[q] += 1
                        cost[a] += 1
                        cost[c] += 1
                    elif k == 1:
                        for x in S:
                            pass
                        bx = x * D + j
                        cost[x] -= 1
                        cost[LA[bx]] -= 1
                        cost[LB[bx]] -= 1
                    S.add(min(q, a, c))
        inK[q] = 1
        self.kpos[q] = len(self.K)
        self.K.append(q)

    def remove(self, p):
        inK = self.inK
        if not inK[p]:
            return
        D, LA, LB, cost, lines = self.D, self.LA, self.LB, self.cost, self.lines
        inK[p] = 0
        b = p * D
        for j in range(D):
            a = LA[b + j]
            if inK[a]:
                c = LB[b + j]
                if inK[c]:
                    S = lines[j]
                    S.discard(min(p, a, c))
                    k = len(S)
                    if k == 0:
                        self._unc_add(j)
                        cost[p] -= 1
                        cost[a] -= 1
                        cost[c] -= 1
                    elif k == 1:
                        for x in S:
                            pass
                        bx = x * D + j
                        cost[x] += 1
                        cost[LA[bx]] += 1
                        cost[LB[bx]] += 1
        i = self.kpos[p]
        last = self.K.pop()
        if last != p:
            self.K[i] = last
            self.kpos[last] = i
        self.kpos[p] = -1

    def load(self, pts):
        self.reset()
        for p in pts:
            self.add(p)

    def missing(self, j, x):
        D, inK = self.D, self.inK
        a, c = self.LA[x * D + j], self.LB[x * D + j]
        return [y for y in (x, a, c) if not inK[y]]

    def greedy(self, rng):
        while self.unc:
            j = rng.choice(self.unc)
            best, opts = 4, []
            for x in self.T["anchors"][j]:
                m = len(self.missing(j, x))
                if m < best:
                    best, opts = m, [x]
                elif m == best:
                    opts.append(x)
            for y in self.missing(j, rng.choice(opts)):
                self.add(y)

    def prune(self, rng):
        while True:
            free = [p for p in self.K if self.cost[p] == 0 and not self.fz[p]]
            if not free:
                return
            self.remove(rng.choice(free))


# ---- search ------------------------------------------------------------------------------
def search(T, seed=1, budget=60.0, start=None, frozen=(), noise=0.1, tenure=5, stall=5000,
           ncand=8, fresh=0, stop_at=0):
    t0 = time.time()
    rng = random.Random(seed)
    W = FW(T, frozen)
    if start:
        for p in start:
            W.add(p)
    W.greedy(rng)
    W.prune(rng)
    best = sorted(W.K)
    reached = [(len(best), round(time.time() - t0, 1))]
    traj = list(best)
    P, D, LA, LB, fz = T["P"], T["D"], T["LA"], T["LB"], W.fz
    anchors = T["anchors"]
    added_at = [-10 ** 9] * P
    st = {"it": 0, "target": len(W.K) - 1}

    def shrink(protect):
        it = st["it"]
        cost = W.cost  # W.reset() replaces the list, so read it here
        while len(W.K) > st["target"]:
            cand = [x for x in W.K if not fz[x] and x not in protect and it - added_at[x] > tenure]
            if not cand:
                cand = [x for x in W.K if not fz[x] and x not in protect]
                if not cand:
                    return
            if rng.random() < noise:
                x = rng.choice(cand)
            else:
                mc = min(cost[x] for x in cand)
                x = rng.choice([x for x in cand if cost[x] == mc])
            W.remove(x)

    shrink(())
    last, bestu, restarts, fresh_starts = 0, 10 ** 9, 0, 0
    while time.time() - t0 < budget:
        st["it"] += 1
        it = st["it"]
        u = len(W.unc)
        if u == 0:
            if len(W.K) < len(traj):
                traj = sorted(W.K)
            if len(W.K) < len(best):
                best = sorted(W.K)
                reached.append((len(best), round(time.time() - t0, 1)))
                if len(best) <= stop_at:
                    break
            st["target"] = len(W.K) - 1
            shrink(())
            last, bestu = it, len(W.unc)
            continue
        if u < bestu:
            bestu, last = u, it
        if it - last > stall:
            restarts += 1
            if fresh and restarts % fresh == 0:
                W.reset()
                W.greedy(rng)
                W.prune(rng)
                fresh_starts += 1
                traj = sorted(W.K)
                if len(traj) < len(best):
                    best = list(traj)
                    reached.append((len(best), round(time.time() - t0, 1)))
            else:
                W.load(traj)
            free = [x for x in W.K if not fz[x]]
            if free:
                W.remove(rng.choice(free))
            st["target"] = len(W.K)
            for i in range(P):
                added_at[i] = -10 ** 9
            last, bestu = it, len(W.unc)
            continue
        j = rng.choice(W.unc)
        mn, opts, near = 4, [], []
        inK = W.inK
        for x in anchors[j]:
            bx = x * D + j
            m = (not inK[x]) + (not inK[LA[bx]]) + (not inK[LB[bx]])
            if m < mn:
                near = opts if m == mn - 1 else []
                mn, opts = m, [x]
            elif m == mn:
                opts.append(x)
            elif m == mn + 1:
                near.append(x)
        if rng.random() < noise:
            x = rng.choice(opts + near)
        else:
            sample = opts if len(opts) <= ncand else rng.sample(opts, ncand)
            bu, bl = None, []
            for x in sample:
                miss = W.missing(j, x)
                for y in miss:
                    W.add(y)
                uu = len(W.unc)
                for y in miss:
                    W.remove(y)
                if bu is None or uu < bu:
                    bu, bl = uu, [x]
                elif uu == bu:
                    bl.append(x)
            x = rng.choice(bl)
        pick = W.missing(j, x)
        for y in pick:
            W.add(y)
            added_at[y] = it
        shrink(set(pick))
    return {"best": len(best), "points": best, "reached": reached, "steps": st["it"],
            "restarts": restarts, "fresh_starts": fresh_starts,
            "seconds": round(time.time() - t0, 1)}


# ---- parallel probe ----------------------------------------------------------------------
_TAB = {}


def _tables(n):
    if n not in _TAB:
        _TAB.clear()
        _TAB[n] = setup(n)
    return _TAB[n]


def task(t):
    """t = (n, label, start, frozen, seed, budget, noise, fresh)."""
    n, label, start, frozen, seed, budget, noise, fresh = t
    t0 = time.time()
    try:
        T = _tables(n)
        tset = round(time.time() - t0, 1)
        r = search(T, seed=seed, budget=max(5.0, budget - (time.time() - t0)), start=start,
                   frozen=frozen, noise=noise, fresh=fresh)
        r["verified"] = independent_check(n, r["points"])
        r["label"], r["seed"], r["noise"], r["setup_s"] = label, seed, noise, tset
        r["start_size"] = len(set(start or [])) if start else None
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


def downloads():
    """The published sets of Chung, Du and Wesley (dualverse-ai/station_data_v2), fetched at
    run time and re-checked; returns {n: indices} for what could be fetched."""
    from k5 import cdw
    out, notes = {}, []
    for n, size in ((3, 13), (4, 27), (5, 53)):
        try:
            pts, note = cdw.load_set(n, size)
            ok = k5walk.is_kakeya(n, pts) and len(pts) == size
            idx = sorted(sum(c * 3 ** (n - 1 - k) for k, c in enumerate(p)) for p in pts)
            notes.append("d%d_%d %s" % (n, size, "ok" if ok else "BAD"))
            if ok:
                out[n] = idx
        except Exception as e:  # noqa: BLE001
            notes.append("d%d_%d failed %s" % (n, size, str(e)[:80]))
    return out, "; ".join(notes)


K7_PLANE = [0, 1, 2, 3, 5, 6, 8]  # {(0,0),(0,1),(0,2),(1,0),(1,2),(2,0),(2,2)}: a 7-point Kakeya set in F_3^2


def plan(n, seconds, sets):
    """Start sets for dimension n (6 or 7) from the downloaded sets."""
    full1 = full(1)
    starts = []
    if n == 6:
        starts.append(("greedy", None, (), 1, 0.10, 5))
        if 5 in sets:
            starts.append(("layer53-frozen", None, sets[5], 2, 0.10, 5))
            starts.append(("prod53x3", product(full1, 1, sets[5], 5), (), 3, 0.10, 0))
        else:
            starts.append(("greedy", None, (), 2, 0.15, 5))
            starts.append(("greedy", None, (), 3, 0.05, 5))
        starts.append(("classical127", classical(6), (), 4, 0.10, 0))
        if 4 in sets:
            starts.append(("prod27x7", product(sets[4], 4, K7_PLANE, 2), (), 5, 0.15, 0))
        else:
            starts.append(("greedy", None, (), 5, 0.20, 5))
    else:
        starts.append(("greedy", None, (), 1, 0.10, 5))
        starts.append(("classical255", classical(7), (), 2, 0.10, 0))
        if 5 in sets:
            starts.append(("prod53x7", product(sets[5], 5, K7_PLANE, 2), (), 3, 0.10, 0))
        else:
            starts.append(("greedy", None, (), 3, 0.20, 5))
        if 4 in sets and 3 in sets:
            starts.append(("prod27x13", product(sets[4], 4, sets[3], 3), (), 4, 0.10, 0))
        else:
            starts.append(("greedy", None, (), 4, 0.05, 5))
        starts.append(("greedy", None, (), 5, 0.20, 5))
    return [(n, lab, st, fr, seed, seconds, noise, fresh)
            for (lab, st, fr, seed, noise, fresh) in starts]


def probe(n, seconds=450, workers=5):
    t0 = time.time()
    out = []
    try:
        sets, note = downloads()
        out.append("fetch " + note)
        tasks = plan(n, seconds - (time.time() - t0) - 15, sets)
        res, mode = run_pool(tasks[:workers], workers)
        out.append("mode=%s cpu=%s" % (mode, os.cpu_count()))
        best = None
        for r in res:
            if "error" in r:
                out.append("%s ERROR %s" % (r["label"], r["error"]))
                continue
            rr = r["reached"]
            short = rr if len(rr) <= 8 else rr[:2] + ["..."] + rr[-5:]
            out.append("%s seed=%d noise=%.2f start=%s setup=%ss: best %d verified=%s reached %s; "
                       "%d steps, %d restarts, %d fresh" % (
                           r["label"], r["seed"], r["noise"], r["start_size"], r["setup_s"],
                           r["best"], r["verified"], short, r["steps"], r["restarts"],
                           r["fresh_starts"]))
            if r["verified"] and (best is None or r["best"] < best["best"]):
                best = r
        if best is not None:
            out.append("BEST n=%d size=%d (classical %d, Lev guess %d) label=%s POINTS %s" % (
                n, best["best"], CLASSICAL[n], LEV[n], best["label"],
                ",".join(map(str, best["points"]))))
    except Exception:  # noqa: BLE001
        out.append("ERROR " + traceback.format_exc()[-700:])
    out.append("total %.0fs" % (time.time() - t0))
    return " || ".join(out)
