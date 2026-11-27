# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Exhaustive scan for rainbow triples of F_3^m (rainbow/SCAN.md, rainbow/RAINBOW3.md).

A triple (A, M, C) of subsets of F_3^m is rainbow if every d in F_3^m (d = 0 included) has
a witness a in A with a + d in M and a + 2d in C.  tau(m) is the least |A| |M| |C|.  By the
tripling lemma of TRIPLES.md, a rainbow triple gives a Kakeya set of 3 |A| |M| |C| points in
F_3^(3m), and products of rainbow triples are rainbow.

Normal form (each step is proved in SCAN.md).  Every rainbow triple with product P can be
moved, keeping the three sizes as a multiset and keeping it rainbow, to one with
  N1  |A| <= |M| <= |C|                         (the six role permutations);
  N2  A the canonical representative of its affine class   (affine maps);
  N3  0 in M, and M the least (as a bit mask) of its translates M - m, m in M
                                                (translations (A, M + t, C + 2t)).
Then |M| >= max(|A|, ceil(3^m / |A|)) (pair bound) and |A| |M|^2 <= P.

Method.  For a fixed A and a fixed size s of M, every M of size s in normal form N3 whose
set M - A is the whole space is enumerated (a pruned walk over sorted index sets); for each,
the least C is found exactly as a minimum hitting set of the sets
  P_d = {a + 2d : a in A, a + d in M}   (d in F_3^m),
since (A, M, C) is rainbow iff C meets every P_d.

Points of F_3^m are indexed by sum x_k 3^(m-1-k); the zero vector has index 0."""
import itertools
import os
import sys
import time

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)


# ---------------------------------------------------------------- the space

def points(m):
    return list(itertools.product(range(3), repeat=m))


def to_idx(p):
    i = 0
    for c in p:
        i = 3 * i + int(c) % 3
    return i


def to_pt(i, m):
    return tuple((i // 3 ** (m - 1 - k)) % 3 for k in range(m))


_TAB = {}


def tables(m):
    """(N, add, neg): add[i][j] = index of p_i + p_j, neg[i] = index of -p_i (= 2 p_i)."""
    if m not in _TAB:
        P = points(m)
        add = [[to_idx(tuple((a + b) % 3 for a, b in zip(p, q))) for q in P] for p in P]
        neg = [to_idx(tuple((-a) % 3 for a in p)) for p in P]
        _TAB[m] = (len(P), add, neg)
    return _TAB[m]


# ---------------------------------------------------------------- independent checker

def check_triple(m, A, M, C):
    """Straight from the definition, with coordinate tuples (A, M, C: index lists)."""
    At = [to_pt(i, m) for i in A]
    Ms = set(to_pt(i, m) for i in M)
    Cs = set(to_pt(i, m) for i in C)
    for d in points(m):
        ok = False
        for a in At:
            if (tuple((x + y) % 3 for x, y in zip(a, d)) in Ms
                    and tuple((x + 2 * y) % 3 for x, y in zip(a, d)) in Cs):
                ok = True
                break
        if not ok:
            return False
    return True


def product_triple(m1, T1, m2, T2):
    """(A1 x A2, M1 x M2, C1 x C2) as index lists of F_3^(m1 + m2)."""
    N2 = 3 ** m2
    return tuple([x * N2 + y for x in X1 for y in X2] for X1, X2 in zip(T1, T2))


# ---------------------------------------------------------------- affine classes (N2)

def affine_dim(m, S):
    """Dimension of the affine hull of the point set S (index list)."""
    N, add, neg = tables(m)
    S = list(S)
    o = neg[S[0]]
    span, k = {0}, 0
    for p in S[1:]:
        q = add[p][o]
        if q not in span:
            span = {add[s][t] for s in span for t in (0, q, neg[q])}
            k += 1
    return k


def canon(m, S):
    """(k, form, stab): k = affine dimension of S; form = the least sorted index tuple of
    the images of S under the affine maps that send an ordered affine frame
    (p_0, p_1, ..., p_k) of points of S to (0, e_1, ..., e_k) (coordinates on the hull);
    stab = the number of frames giving that least image, which is the order of the
    stabiliser of S in AGL(k, 3) acting on its hull.  Two sets are affinely equivalent iff
    they have the same form."""
    N, add, neg = tables(m)
    S = list(S)
    k = affine_dim(m, S)
    pw = [3 ** (m - 1 - i) for i in range(k)]
    best, cnt = None, 0
    for fr in itertools.permutations(S, k + 1):
        no = neg[fr[0]]
        cur = [(0, 0)]
        for i in range(k):
            v = add[fr[i + 1]][no]
            v2 = neg[v]
            g = pw[i]
            cur = [(add[u][w], h + c * g) for (u, h) in cur
                   for (w, c) in ((0, 0), (v, 1), (v2, 2))]
        tab = dict(cur)
        if len(tab) < len(cur):  # the frame is affinely dependent
            continue
        img = tuple(sorted(tab[add[p][no]] for p in S))
        if best is None or img < best:
            best, cnt = img, 1
        elif img == best:
            cnt += 1
    return k, best, cnt


_CL = {}


def classes(m, x):
    """One canonical representative of every affine class of x-subsets of F_3^m, as a
    sorted list of (k, form, stab).  Built incrementally: every x-set contains an
    (x-1)-set, which some affine map sends to a class representative R; the same map
    sends the x-set to R + {p} for a point p outside R."""
    key = (m, x)
    if key not in _CL:
        N = tables(m)[0]
        out = {}
        if x == 1:
            c = canon(m, [0])
            out[c[1]] = c
        else:
            for (_, R, _) in classes(m, x - 1):
                Rs = set(R)
                for p in range(N):
                    if p not in Rs:
                        c = canon(m, sorted(R + (p,)))
                        if c[1] not in out:
                            out[c[1]] = c
        _CL[key] = sorted(out.values())
    return _CL[key]


def agl_order(k):
    o = 3 ** k
    for i in range(k):
        o *= 3 ** k - 3 ** i
    return o


# ---------------------------------------------------------------- minimum hitting set

def _hit(sets, j):
    """A bit mask of at most j points meeting every mask in `sets`, or None.  Pruned by
    the greedy packing bound: pairwise disjoint sets need distinct points."""
    if not sets:
        return 0
    if j <= 0:
        return None
    sets = sorted(sets, key=int.bit_count)
    used, cnt = 0, 0
    for s in sets:
        if not s & used:
            used |= s
            cnt += 1
            if cnt > j:
                return None
    b = sets[0]
    while b:
        low = b & -b
        b ^= low
        r = _hit([s for s in sets if not s & low], j - 1)
        if r is not None:
            return r | low
    return None


def min_hit(sets, k):
    """A least hitting set (bit mask) of the masks in `sets` if its size is <= k, else
    None.  Exact: _hit is exhaustive, and the size is lowered until _hit fails."""
    r = _hit(sets, k)
    if r is None:
        return None
    while True:
        r2 = _hit(sets, r.bit_count() - 1)
        if r2 is None:
            return r
        r = r2


def brute_min_hit(sets, n):
    """Reference: the least size of a subset of range(n) meeting every mask (small n)."""
    sets = list(sets)
    for k in range(n + 1):
        for comb in itertools.combinations(range(n), k):
            msk = 0
            for p in comb:
                msk |= 1 << p
            if all(s & msk for s in sets):
                return k
    return None


# ---------------------------------------------------------------- the scan for one A

_PREP = {}


def prep(m, A):
    """cov[p]: mask of the d that the point p can serve as M-point a + d (a in A);
    opt[d]: the pairs (index of a + d, bit of a + 2d) over a in A; sub[p][q] = p - q."""
    key = (m, tuple(A))
    if key not in _PREP:
        N, add, neg = tables(m)
        cov = [0] * N
        opt = []
        for d in range(N):
            row = []
            for a in A:
                p = add[a][d]
                cov[p] |= 1 << d
                row.append((p, 1 << add[a][neg[d]]))
            opt.append(tuple(row))
        sub = [[add[p][neg[q]] for q in range(N)] for p in range(N)]
        suf = [0] * (N + 1)
        for p in range(N - 1, -1, -1):
            suf[p] = suf[p + 1] | cov[p]
        if len(_PREP) > 64:
            _PREP.clear()
        _PREP[key] = (N, cov, opt, sub, suf)
    return _PREP[key]


def translation_least(Ml, mask, sub):
    """True iff mask(M) <= mask(M - m) for every m in M (M = Ml, Ml[0] = 0)."""
    for m in Ml[1:]:
        t = 0
        for p in Ml:
            t |= 1 << sub[p][m]
        if t < mask:
            return False
    return True


def prefixes(N, s):
    """The starts (p_1, ..., p_j) of the sorted M = {0 < p_1 < ...} of size s used to split
    the scan into tasks; every such M has exactly one of them."""
    depth = min(s - 1, 1 if s <= 7 else 2)
    out = []
    for pre in itertools.combinations(range(1, N), depth):
        if not pre or pre[-1] <= N - 1 - (s - 1 - depth):
            out.append(pre)
    return out


def scan_task(task):
    """task = (m, A, s, prefix, bound, deadline, find).  Every M of size s in normal form
    N3 starting {0} + prefix with M - A = F_3^m; for each, the least C.  Records the least
    product |A| s |C| < bound (a triple found) and, with find, stops at the first."""
    m, A, s, prefix, bound, deadline, find = task
    A = list(A)
    N, cov, opt, sub, suf = prep(m, A)
    x = len(A)
    FULL = (1 << N) - 1
    cmax = (bound - 1) // (x * s)
    st = {"nodes": 0, "covers": 0, "least": 0, "best": None, "M": None, "C": None,
          "timeout": False}

    def leaf(Ml):
        st["covers"] += 1
        mask = 0
        for p in Ml:
            mask |= 1 << p
        if not translation_least(Ml, mask, sub):
            return
        st["least"] += 1
        sets = set()
        for o in opt:
            c = 0
            for mi, cb in o:
                if mask >> mi & 1:
                    c |= cb
            if not c:  # impossible when M - A is everything; kept as a guard
                return
            sets.add(c)
        r = min_hit(sets, cmax)
        if r is not None:
            prod = x * s * r.bit_count()
            if st["best"] is None or prod < st["best"]:
                st["best"] = prod
                st["M"] = sorted(Ml)
                st["C"] = [i for i in range(N) if r >> i & 1]

    def rec(start, need, covered, Ml):
        st["nodes"] += 1
        if (st["nodes"] & 2047) == 0 and time.time() > deadline:
            st["timeout"] = True
        if st["timeout"]:
            return
        if need == 0:
            if covered == FULL:
                leaf(Ml)
            return
        unc = FULL & ~covered
        if unc & ~suf[start]:
            return
        if unc.bit_count() > need * x:
            return
        for p in range(start, N - need + 1):
            Ml.append(p)
            rec(p + 1, need - 1, covered | cov[p], Ml)
            Ml.pop()
            if st["timeout"] or (find and st["best"] is not None):
                return

    t0 = time.time()
    Ml = [0] + list(prefix)
    covered = 0
    for p in Ml:
        covered |= cov[p]
    start = Ml[-1] + 1 if prefix else 1
    if time.time() <= deadline:
        rec(start, s - len(Ml), covered, Ml)
    else:
        st["timeout"] = True
    return {"m": m, "x": x, "A": A, "s": s, "prefix": list(prefix), "bound": bound,
            "complete": not st["timeout"], "best": st["best"], "M": st["M"], "C": st["C"],
            "nodes": st["nodes"], "covers": st["covers"], "least": st["least"],
            "secs": round(time.time() - t0, 2)}


# ---------------------------------------------------------------- sweeps

def sizes_for(m, x, bound):
    """The sizes s of M allowed by the normal form below `bound`: s >= max(x, ceil(3^m/x))
    and x s^2 < bound."""
    N = 3 ** m
    s = max(x, -(-N // x))
    out = []
    while x * s * s < bound:
        out.append(s)
        s += 1
    return out


def make_tasks(m, bound, xs, deadline, full_only, find=False):
    N = 3 ** m
    tasks = []
    for x in xs:
        for (k, A, _) in classes(m, x):
            if full_only and k < m:
                continue
            for s in sizes_for(m, x, bound):
                for pre in prefixes(N, s):
                    tasks.append((m, tuple(A), s, pre, bound, deadline, find))
    # big tasks first (large s, small first index), for load balance
    tasks.sort(key=lambda t: (-t[2], t[3]))
    return tasks


def sweep(m, bound, xs, deadline, workers=8, full_only=False, find=False):
    """Every class representative A of every size in xs (with full_only, only the
    full-dimensional ones, by lemma H), every allowed size of M, below `bound`."""
    tasks = make_tasks(m, bound, xs, deadline, full_only, find)
    if workers <= 1 or len(tasks) <= 1:
        return [scan_task(t) for t in tasks]
    import multiprocessing as mp
    pool = mp.get_context("spawn").Pool(max(1, min(workers, len(tasks))))
    try:
        return list(pool.imap_unordered(scan_task, tasks, chunksize=1))
    finally:
        pool.close()
        pool.join()


def summarise(res):
    by = {}
    for r in res:
        g = by.setdefault(r["x"], {"tasks": 0, "reps": set(), "complete": 0, "best": None,
                                   "covers": 0, "least": 0, "nodes": 0, "maxsecs": 0.0,
                                   "ex": None})
        g["tasks"] += 1
        g["reps"].add(tuple(r["A"]))
        g["complete"] += 1 if r["complete"] else 0
        g["nodes"] += r["nodes"]
        g["covers"] += r["covers"]
        g["least"] += r["least"]
        g["maxsecs"] = max(g["maxsecs"], r["secs"])
        if r["best"] is not None and (g["best"] is None or r["best"] < g["best"]):
            g["best"] = r["best"]
            g["ex"] = (r["A"], r["M"], r["C"])
    for g in by.values():
        g["reps"] = len(g["reps"])
    return by


def describe(by):
    return " | ".join("x%d reps=%d tasks=%d complete=%d covers=%d least=%d best=%s maxs=%.1f"
                      % (x, g["reps"], g["tasks"], g["complete"], g["covers"], g["least"],
                         g["best"], g["maxsecs"]) for x, g in sorted(by.items()))
