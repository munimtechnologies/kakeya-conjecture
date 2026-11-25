# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Is there a 52-point Kakeya set in F_3^5 near the known 53-point one?

The smallest Kakeya set known in F_3^5 has 53 points.  It was found by Chung, Du and Wesley
(arXiv:2608.23691) and is published as kakeya_F3_d5_53.npy in
github.com/dualverse-ai/station_data_v2 (artifacts/finite_kakeya).  This module downloads
that file at run time (it is not copied into this repository), checks it independently, and
asks two questions about it:

1. Exact neighbourhood: delete r points and add r - 1 points (r = 1, 2, 3, every choice).
   Is any resulting set Kakeya?  Every such set has at most 52 points, so a single hit would
   beat the 53-point record and break V. Lev's 2009 guess k_5 = 53.  The search is complete:
   for each deletion R it enumerates every way to restore the directions R destroyed with at
   most r - 1 new points (branching on the lines of one destroyed direction at a time).
2. Local search (walk.py) started from the 53-point set with a 52-point budget.

Controls: the 53-point set plus one extra point must lead back to the 53-point set (r <= 2),
the same with two extra points and deletions containing both (exercises the 2-point branch),
and the published 27-point set in F_3^4 must have no Kakeya neighbour of 26 points (k_4 = 27).
Every set reported as Kakeya is re-checked by walk.is_kakeya.  Standard library only.
"""
import ast
import io
import itertools
import math
import os
import struct
import time
import traceback
import urllib.request
import zipfile

try:
    from . import walk
except ImportError:  # pragma: no cover - flat layout
    import walk

SRC = "dualverse-ai/station_data_v2/main/artifacts/finite_kakeya/"
URLS = ["https://raw.githubusercontent.com/" + SRC,
        "https://media.githubusercontent.com/media/" + SRC]


# ---- download and parse ------------------------------------------------------------------
def fetch(name, timeout=40):
    last = None
    for base in URLS:
        try:
            req = urllib.request.Request(base + name, headers={"User-Agent": "kakeya-check"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = r.read()
            if data.startswith(b"version https://git-lfs"):
                last = "git-lfs pointer at " + base
                continue
            return data
        except Exception as e:  # noqa: BLE001 - report and try the next mirror
            last = "%s: %s" % (type(e).__name__, str(e)[:80])
    raise RuntimeError("fetch %s failed: %s" % (name, last))


_CODES = {"b1": "?", "u1": "B", "i1": "b", "u2": "H", "i2": "h", "u4": "I", "i4": "i",
          "u8": "Q", "i8": "q", "f4": "f", "f8": "d"}


def parse_npy(data):
    """Return (shape, flat values in C order, descr) of a numeric .npy file."""
    if data[:6] != b"\x93NUMPY":
        raise ValueError("not an npy file")
    if data[6] == 1:
        hl, off = struct.unpack("<H", data[8:10])[0], 10
    else:
        hl, off = struct.unpack("<I", data[8:12])[0], 12
    hdr = ast.literal_eval(data[off:off + hl].decode("latin1"))
    descr, shape = hdr["descr"], tuple(hdr["shape"])
    if not isinstance(descr, str) or descr[1:] not in _CODES:
        raise ValueError("unsupported dtype %r" % (descr,))
    endian = ">" if descr[0] == ">" else "<"
    count = 1
    for s in shape:
        count *= s
    vals = struct.unpack_from(endian + str(count) + _CODES[descr[1:]], data, off + hl)
    vals = [int(round(v)) if isinstance(v, float) else int(v) for v in vals]
    if hdr.get("fortran_order") and len(shape) > 1:
        strides = [1]
        for s in shape[:-1]:
            strides.append(strides[-1] * s)
        vals = [vals[sum(i * st for i, st in zip(ix, strides))]
                for ix in itertools.product(*[range(s) for s in shape])]
    return shape, vals, descr


def to_points(shape, vals, n):
    """Interpret an array as a point set of F_3^n (list of n-tuples)."""
    def digits(i):
        return tuple((i // 3 ** (n - 1 - k)) % 3 for k in range(n))
    if len(shape) == 2 and shape[1] == n:
        pts = [tuple(v % 3 for v in vals[r * n:(r + 1) * n]) for r in range(shape[0])]
    elif len(shape) == 2 and shape[0] == n and shape[1] != n:
        m = shape[1]
        pts = [tuple(vals[k * m + r] % 3 for k in range(n)) for r in range(m)]
    elif shape == (3,) * n or shape == (3 ** n,):
        pts = [digits(i) for i, v in enumerate(vals) if v]
    elif len(shape) == 1:
        pts = [digits(v) for v in vals]
    else:
        raise ValueError("cannot read shape %r as points of F_3^%d" % (shape, n))
    return sorted(set(pts))


def load_set(n, size):
    """Download the published size-point Kakeya set in F_3^n; returns (points, note)."""
    name = "kakeya_F3_d%d_%d.npy" % (n, size)
    try:
        shape, vals, descr = parse_npy(fetch(name))
        return to_points(shape, vals, n), "%s shape=%s dtype=%s" % (name, shape, descr)
    except Exception as e:  # noqa: BLE001 - fall back to the certificate bundle
        first = "%s: %s" % (type(e).__name__, str(e)[:100])
    z = zipfile.ZipFile(io.BytesIO(fetch("station_point_certificates.npz")))
    key = "station_d%d_p3_points.npy" % n
    shape, vals, descr = parse_npy(z.read(key))
    return to_points(shape, vals, n), "npz:%s shape=%s dtype=%s (npy failed: %s)" % (
        key, shape, descr, first)


# ---- tables and simple invariants --------------------------------------------------------
_T = {}


def tables(n):
    if n not in _T:
        _T[n] = walk.setup(n)
    return _T[n]


def indices(T, pts):
    return sorted(T["idx"][tuple(p)] for p in pts)


def lines_inside(T, inK):
    return [[ln for ln in T["lines"][j] if inK[ln[0]] and inK[ln[1]] and inK[ln[2]]]
            for j in range(T["D"])]


def invariants(T, K):
    """Lines per direction, point degrees, hyperplane profiles, central symmetries."""
    n, pts = T["n"], T["pts"]
    inK = bytearray(T["P"])
    for p in K:
        inK[p] = 1
    KL = lines_inside(T, inK)
    per_dir = {}
    for L in KL:
        per_dir[len(L)] = per_dir.get(len(L), 0) + 1
    deg = {p: 0 for p in K}
    for L in KL:
        for ln in L:
            for x in ln:
                deg[x] += 1
    deg_hist = {}
    for v in deg.values():
        deg_hist[v] = deg_hist.get(v, 0) + 1
    prof = {}
    for d in T["dirs"]:
        c = [0, 0, 0]
        for p in K:
            c[sum(a * b for a, b in zip(d, pts[p])) % 3] += 1
        key = tuple(sorted(c))
        prof[key] = prof.get(key, 0) + 1
    S = set(pts[p] for p in K)
    centres = [c for c in pts
               if all(tuple((2 * a - b) % 3 for a, b in zip(c, x)) in S for x in S)]
    removable = sum(1 for p in K if walk.is_kakeya(n, [pts[q] for q in K if q != p]))
    return {"lines": sum(len(L) for L in KL), "per_dir": dict(sorted(per_dir.items())),
            "deg": dict(sorted(deg_hist.items())),
            "profiles": dict(sorted(prof.items(), key=lambda kv: -kv[1])),
            "centres": len(centres), "removable": removable}


# ---- exact neighbourhood -----------------------------------------------------------------
def _covered_by_new(T, base, j, added):
    D, LA, LB = T["D"], T["LA"], T["LB"]
    for q in added:
        b = q * D + j
        if base[LA[b]] and base[LB[b]]:
            return True
    return False


def complete(T, base, U, s, out, prefix):
    """Append to `out` every tuple A (|A| <= s, built line by line) such that adding A to the
    set `base` (bytearray, modified and restored) gives a line in every direction of U.  Each
    direction of U needs a line whose missing points all lie in A, so branching over the
    lines of one direction of U with at most s missing points is exhaustive."""
    if not U:
        out.append(tuple(prefix))
        return
    if s == 0:
        return
    lines = T["lines"]
    if s == 1:
        j0 = U[0]
        for (a, b, c) in lines[j0]:
            m = (not base[a]) + (not base[b]) + (not base[c])
            if m != 1:
                continue
            x = a if not base[a] else (b if not base[b] else c)
            base[x] = 1
            ok = True
            for j in U[1:]:
                if not _covered_by_new(T, base, j, (x,)):
                    ok = False
                    break
            base[x] = 0
            if ok:
                out.append(tuple(prefix) + (x,))
        return
    best = None
    for j in U:
        opts = []
        cheap = 0
        for ln in lines[j]:
            miss = [x for x in ln if not base[x]]
            if len(miss) <= s:
                opts.append(miss)
                if len(miss) < s:
                    cheap += 1
        if not opts:
            return  # direction j cannot be restored with s new points
        key = (cheap, len(opts))
        if best is None or key < best[0]:
            best = (key, opts)
    for miss in best[1]:
        for x in miss:
            base[x] = 1
        prefix.extend(miss)
        U2 = [j for j in U if not _covered_by_new(T, base, j, miss)]
        complete(T, base, U2, s - len(miss), out, prefix)
        del prefix[len(prefix) - len(miss):]
        for x in miss:
            base[x] = 0


class Nbhd:
    def __init__(self, T, K):
        self.T = T
        self.K = sorted(K)
        inK = bytearray(T["P"])
        for p in self.K:
            inK[p] = 1
        self.inK = inK
        self.KL = lines_inside(T, inK)
        self.aff = {p: set() for p in self.K}
        for j, L in enumerate(self.KL):
            for ln in L:
                for x in ln:
                    self.aff[x].add(j)
        self.uncov0 = [j for j in range(T["D"]) if not self.KL[j]]

    def uncovered(self, R):
        Rs = set(R)
        cand = set()
        for p in R:
            cand |= self.aff[p]
        U = list(self.uncov0)
        for j in sorted(cand):
            if all(ln[0] in Rs or ln[1] in Rs or ln[2] in Rs for ln in self.KL[j]):
                U.append(j)
        return U

    def search(self, R, s):
        """All sets (K minus R) plus A, |A| <= s, that are Kakeya (as frozensets)."""
        base = bytearray(self.inK)
        for p in R:
            base[p] = 0
        out = []
        complete(self.T, base, self.uncovered(R), s, out, [])
        Rs = set(R)
        keep = [p for p in self.K if p not in Rs]
        return [frozenset(keep + list(A)) for A in out]


_NB = {}


def _nb(n, K):
    key = (n, tuple(K))
    if key not in _NB:
        _NB.clear()
        _NB[key] = Nbhd(tables(n), K)
    return _NB[key]


def nb_task(task):
    """Deletions R of size r whose smallest member is K[i]; add r - 1 points."""
    n, K, r, i, deadline = task
    NB = _nb(n, K)
    m = len(K)
    total = math.comb(m - 1 - i, r - 1)
    done = 0
    found = set()
    for rest in itertools.combinations(range(i + 1, m), r - 1):
        if time.time() > deadline:
            break
        R = (K[i],) + tuple(K[t] for t in rest)
        for f in NB.search(R, r - 1):
            found.add(tuple(sorted(f)))
        done += 1
    return {"r": r, "i": i, "done": done, "total": total, "found": sorted(found)[:50],
            "nfound": len(found)}


def run_pool(fn, tasks, workers, ordered=False):
    if workers > 1:
        try:
            import multiprocessing as mp
            ctx = mp.get_context("spawn")
            with ctx.Pool(workers) as p:
                if ordered:
                    return p.map(fn, tasks, chunksize=1), "pool%d" % workers
                return list(p.imap_unordered(fn, tasks)), "pool%d" % workers
        except Exception as e:  # noqa: BLE001 - fall back to this process
            note = "serial(%s:%s)" % (type(e).__name__, str(e)[:60])
    else:
        note = "serial"
    return [fn(t) for t in tasks], note


def neighbourhood(n, K, rmax, deadline, workers):
    tasks = [(n, list(K), r, i, deadline) for r in range(1, rmax + 1) for i in range(len(K))]
    res, mode = run_pool(nb_task, tasks, workers)
    by_r = {}
    for x in res:
        a = by_r.setdefault(x["r"], {"done": 0, "total": 0, "found": set()})
        a["done"] += x["done"]
        a["total"] += x["total"]
        a["found"].update(x["found"])
    return by_r, mode


def check_found(n, T, sets):
    return [(len(f), walk.is_kakeya(n, [T["pts"][p] for p in f])) for f in sets]


def workers_default(cap=5):
    return max(1, min(cap, (os.cpu_count() or 1) - 2))


# ---- probe A: fetch, invariants, controls, exact neighbourhood ---------------------------
def probe_a(seconds=470, rmax=4):
    t0 = time.time()
    deadline = t0 + seconds
    out = []
    try:
        T5, T4 = tables(5), tables(4)
        P53, note5 = load_set(5, 53)
        out.append("fetch5 %s n=%d kakeya=%s" % (note5, len(P53), walk.is_kakeya(5, P53)))
        K53 = indices(T5, P53)
        out.append("POINTS53 " + ",".join(map(str, K53)))
        inv = invariants(T5, K53)
        out.append("inv53 %r" % (inv,))
        try:
            P27, note4 = load_set(4, 27)
            K27 = indices(T4, P27)
            out.append("fetch4 %s n=%d kakeya=%s inv27 %r" % (
                note4, len(P27), walk.is_kakeya(4, P27), invariants(T4, K27)))
        except Exception as e:  # noqa: BLE001
            K27 = None
            out.append("fetch4 failed %s" % e)
        # positive controls
        extra = [p for p in range(T5["P"]) if p not in set(K53)][:2]
        K54 = sorted(K53 + extra[:1])
        NB54 = Nbhd(T5, K54)
        hit1 = sum(1 for R in itertools.combinations(K54, 1)
                   for f in NB54.search(R, 0) if f == frozenset(K53))
        hit2 = sum(1 for R in itertools.combinations(K54, 2)
                   for f in NB54.search(R, 1) if f == frozenset(K53))
        K55 = sorted(K53 + extra)
        NB55 = Nbhd(T5, K55)
        hit3 = sum(1 for x in K53 for f in NB55.search(tuple(extra) + (x,), 2)
                   if f == frozenset(K53))
        out.append("controls+ K53 recovered: r1 %d (expect 1), r2 %d (expect 53), "
                   "r3 with both extras %d (expect 53; both 53s assume removable=0)"
                   % (hit1, hit2, hit3))
        if K27 is not None:
            by4, _ = neighbourhood(4, K27, 3, deadline, 1)
            out.append("control- F3^4 27-set, r<=3: " + "; ".join(
                "r%d %d/%d found %d" % (r, a["done"], a["total"], len(a["found"]))
                for r, a in sorted(by4.items())))
        out.append("t=%.0fs" % (time.time() - t0))
        w = workers_default()
        by5, mode = neighbourhood(5, K53, rmax, deadline - 10, w)
        out.append("mode=%s cpu=%s" % (mode, os.cpu_count()))
        for r, a in sorted(by5.items()):
            s = "MAIN r%d (delete %d, add <= %d): %d/%d deletions done, Kakeya sets found %d" % (
                r, r, r - 1, a["done"], a["total"], len(a["found"]))
            if a["found"]:
                fl = sorted(a["found"], key=len)[:3]
                s += " checks %r SETS %s" % (check_found(5, T5, fl),
                                             " | ".join(",".join(map(str, f)) for f in fl))
            out.append(s)
    except Exception:  # noqa: BLE001 - the report is the product
        out.append("ERROR " + traceback.format_exc()[-700:])
    out.append("total %.0fs" % (time.time() - t0))
    return " || ".join(out)


# ---- probe B: local search from the 53-point set -----------------------------------------
def walk_task(task):
    n, start, seed, budget, noise, fresh = task
    return walk.run(n, "trivial", seed, budget, noise=noise, fresh=fresh, T=tables(n),
                    start=start)


def probe_b(seconds=470, seeds=(11, 12, 13, 14, 15, 16),
            noises=(0.05, 0.1, 0.15, 0.2, 0.1, 0.3)):
    t0 = time.time()
    out = []
    try:
        T5 = tables(5)
        P53, note5 = load_set(5, 53)
        K53 = indices(T5, P53)
        out.append("fetch5 %s n=%d kakeya=%s" % (note5, len(P53), walk.is_kakeya(5, P53)))
        w = workers_default()
        budget = max(30.0, seconds - (time.time() - t0) - 25)
        # even workers never restart from a fresh greedy set (they stay near the 53-set)
        tasks = [(5, K53, seeds[k % len(seeds)], budget, noises[k % len(noises)],
                  0 if k % 2 == 0 else 4) for k in range(w)]
        res, mode = run_pool(walk_task, tasks, w, ordered=True)
        out.append("mode=%s budget=%.0fs" % (mode, budget))
        for t, r in zip(tasks, res):
            s = ("seed=%d noise=%.2f fresh=%d: best %d (verified %s) reached %s; %d steps, "
                 "%d restarts" % (r["seed"], t[4], t[5], r["best"], r["verified"], r["reached"],
                                  r["iterations"], r["restarts"]))
            if r["best"] <= 52:
                s += " POINTS " + ",".join(map(str, r["points"]))
            out.append(s)
    except Exception:  # noqa: BLE001
        out.append("ERROR " + traceback.format_exc()[-700:])
    out.append("total %.0fs" % (time.time() - t0))
    return " || ".join(out)
