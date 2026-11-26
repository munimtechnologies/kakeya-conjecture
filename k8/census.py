# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Census of layered Kakeya sets in F_3^n (standard library only).

A layered set is K = {0} x L0  u  {1} x L1  u  {2} x L2 where L0 is itself a Kakeya set of
F_3^(n-1); k8/layer.py states the criterion (every v in F_3^(n-1) needs some x in L0 with
x + v in L1 and x - v in L2) and searches L1, L2 for one fixed L0.  The 211-point set of
F_3^7 is layered, with L0 = the 105-point set of F_3^6.  Here L0 varies:

* exact_min_extra(m, L0, cap): least |L1| + |L2| by exhaustive enumeration (tiny m only).
  k8/CENSUS.md proves by hand that a layered Kakeya set of F_3^3 has at least 14 points
  (k_3 = 13), so no 13-point Kakeya set of F_3^3 has a Kakeya plane section;
  k8/test_census.py checks this exhaustively.
* census_phase: for many different Kakeya sets L0 of F_3^m (from k5/walk.py with fresh seeds,
  plus random extra points, plus the classical set) the least |L1| + |L2| that the layered
  search finds in a few seconds.  For m = 4 a total of at most 52 would beat the 53-point
  record of Chung, Du and Wesley (arXiv:2608.23691) in F_3^5.
* recursive_phase: layered sets of F_3^5 built this way (and the published 53-point set,
  downloaded at run time) used as L0 for F_3^6 (record 105).

Searches give upper bounds only; a census that never goes below 53 is evidence, not proof.
"""
import itertools
import os
import random
import sys
import time
import traceback

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:  # so that k5 and k8 import in spawned workers too
    sys.path.insert(0, _ROOT)

from k8 import check, layer  # noqa: E402


# ---- exact (tiny cases) -------------------------------------------------------------------
def exact_min_extra(m, L0, cap):
    """Least |L1| + |L2| <= cap making {0}xL0 u {1}xL1 u {2}xL2 Kakeya, or None."""
    N = 3 ** m
    add, NEG = layer.arith(m)
    X = sorted(set(L0))
    lines = [[(add(x, v), add(x, NEG[v])) for x in X] for v in range(N)]
    for s in range(cap + 1):
        for a in range(s + 1):
            b = s - a
            if a * b < N:   # v -> (x + v, x - v) is injective, so |L1| |L2| >= 3^m
                continue
            for L1 in itertools.combinations(range(N), a):
                s1 = set(L1)
                cand = [[z for (y, z) in lines[v] if y in s1] for v in range(N)]
                if any(not c for c in cand):
                    continue
                for L2 in itertools.combinations(range(N), b):
                    s2 = set(L2)
                    if all(any(z in s2 for z in c) for c in cand):
                        return s
    return None


# ---- census ------------------------------------------------------------------------------
def _walk_set(T, m, seed, budget):
    from k5 import walk
    r = walk.run(m, seed=seed, budget=budget, T=T)
    return sorted(r["points"]), r["verified"]


def census_phase(m, seconds, per, max_add, record, seed0, walk_budget=1.0):
    """Layered searches over many Kakeya sets L0 of F_3^m.  record = best known size in
    F_3^(m+1); a verified total below it is returned in 'new'."""
    from k5 import walk
    t0 = time.time()
    T = walk.setup(m)
    rng = random.Random(seed0)
    stats = {}
    seen = set()
    best = None
    news = []
    fails = unver = k = 0
    queue = [("classical", layer.classical(m))]
    while time.time() - t0 < seconds - per - walk_budget - 1:
        k += 1
        if queue:
            src, L0 = queue.pop()
        else:
            pts, ok = _walk_set(T, m, seed0 * 100003 + k, walk_budget)
            if not ok:
                fails += 1
                continue
            add = min(rng.choice([0, 0, 0, 1, 1, 2, 3, 4, 5]), max_add)
            sp = set(pts)
            free = [p for p in range(3 ** m) if p not in sp]
            L0 = sorted(set(pts) | set(rng.sample(free, add)))
            src = "walk%d+%d" % (len(pts), add)
        key = tuple(L0)
        if key in seen:
            continue
        seen.add(key)
        if not check.is_kakeya(m, L0):
            fails += 1
            continue
        W = layer.Layer(m, L0)
        r = layer.search(W, seed=seed0 * 7919 + k, budget=per,
                         stop_at=record - 1 - len(L0))
        full = W.full_indices(r["points"])
        ver = check.is_kakeya(m + 1, full)
        if not ver:
            unver += 1
            continue
        e, tot = r["best"], len(full)
        st = stats.setdefault(len(L0), {"tried": 0, "min_extra": 10 ** 9, "hist": {}})
        st["tried"] += 1
        st["min_extra"] = min(st["min_extra"], e)
        st["hist"][e] = st["hist"].get(e, 0) + 1
        if best is None or tot < best[0]:
            best = (tot, len(L0), e, src, full)
        if tot < record:
            news.append(full)
    return {"m": m, "stats": stats, "distinct_L0": len(seen), "fails": fails,
            "unverified": unver, "best": best, "new": news, "record": record,
            "seconds": round(time.time() - t0, 1)}


def recursive_phase(seconds, per5, per6, seed0):
    """L0 in F_3^5 = the published 53-point set, then layered sets of F_3^5 built from fresh
    Kakeya sets of F_3^4; each is the frozen layer of a layered search in F_3^6."""
    from k5 import walk
    t0 = time.time()
    T4 = walk.setup(4)
    rows, news = [], []
    queue = []
    try:
        pts, _ = layer.downloaded(5, 53)
        queue.append(("CDW53", pts))
    except Exception as e:  # noqa: BLE001
        rows.append(("CDW53 download failed %s" % str(e)[:50],))
    k = 0
    while time.time() - t0 < seconds - per6 - per5 - 4:
        k += 1
        if queue:
            src, L05 = queue.pop()
        else:
            p4, ok = _walk_set(T4, 4, seed0 * 100003 + k, 1.0)
            if not ok:
                continue
            W4 = layer.Layer(4, p4)
            r4 = layer.search(W4, seed=seed0 * 31 + k, budget=per5, stop_at=52 - len(p4))
            L05 = W4.full_indices(r4["points"])
            if not check.is_kakeya(5, L05):
                continue
            src = "layered%d+%d" % (len(p4), r4["best"])
        W5 = layer.Layer(5, L05)
        r5 = layer.search(W5, seed=seed0 * 7 + k, budget=per6, stop_at=104 - len(L05))
        full = W5.full_indices(r5["points"])
        ver = check.is_kakeya(6, full)
        rows.append((src, len(L05), r5["best"], len(full), ver, r5["reached"][-1]))
        if ver and len(full) < 105:
            news.append(full)
    return {"rows": rows, "new": news, "seconds": round(time.time() - t0, 1)}


# ---- parallel probe ----------------------------------------------------------------------
def task(t):
    label, phases, seed, seconds = t
    try:
        if phases is None:
            return {"label": label, "rec": recursive_phase(seconds, 8.0, 90.0, seed)}
        out, t0 = [], time.time()
        for (m, sec, per, max_add, record) in phases:
            sec = min(sec, seconds - (time.time() - t0))
            out.append(census_phase(m, sec, per, max_add, record, seed * 10 + m))
        return {"label": label, "phases": out}
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


def _fmt_stats(stats):
    parts = []
    for size in sorted(stats):
        st = stats[size]
        h = " ".join("%dx%d" % (e, c) for e, c in sorted(st["hist"].items()))
        parts.append("|L0|=%d: %d tried, min extra %d [%s]" % (size, st["tried"],
                                                               st["min_extra"], h))
    return "; ".join(parts)


def probe(seconds=430, workers=5):
    t0 = time.time()
    out = []
    try:
        ex = exact_min_extra(2, layer.K7_PLANE, 8)
        out.append("exact F_3^3: |L0|=7 min extra %s -> least layered size %s" % (
            ex, None if ex is None else 7 + ex))
        S = seconds - (time.time() - t0) - 20
        tasks = [("c5a", [(4, S, 4.0, 3, 53)], 11, S),
                 ("c5b", [(4, S, 4.0, 3, 53)], 12, S),
                 ("c5c", [(4, S, 2.5, 5, 53)], 13, S),
                 ("c4+c5", [(3, 80, 1.0, 2, 27), (4, S, 4.0, 3, 53)], 14, S),
                 ("rec6", None, 15, S)]
        res, mode = run_pool(tasks[:workers], workers)
        out.append("mode=%s cpu=%s S=%.0fs" % (mode, os.cpu_count(), S))
        agg = {3: {}, 4: {}}
        best = {3: None, 4: None}
        news = []
        for r in res:
            if "error" in r:
                out.append("%s ERROR %s" % (r["label"], r["error"]))
                continue
            if "rec" in r:
                rec = r["rec"]
                for row in rec["rows"]:
                    if len(row) == 1:
                        out.append("rec6 " + row[0])
                    else:
                        out.append("rec6 %s |L0|=%d extra=%d total=%d verified=%s last=%s" % row)
                news += [(6, f) for f in rec["new"]]
                continue
            for ph in r["phases"]:
                m = ph["m"]
                out.append("%s n=%d: %d distinct L0, %d fails, %d unverified, %.0fs; %s" % (
                    r["label"], m + 1, ph["distinct_L0"], ph["fails"], ph["unverified"],
                    ph["seconds"], _fmt_stats(ph["stats"])))
                for size, st in ph["stats"].items():
                    a = agg[m].setdefault(size, {"tried": 0, "min_extra": 10 ** 9, "hist": {}})
                    a["tried"] += st["tried"]
                    a["min_extra"] = min(a["min_extra"], st["min_extra"])
                    for e, c in st["hist"].items():
                        a["hist"][e] = a["hist"].get(e, 0) + c
                if ph["best"] and (best[m] is None or ph["best"][0] < best[m][0]):
                    best[m] = ph["best"]
                news += [(m + 1, f) for f in ph["new"]]
        for m in (3, 4):
            if agg[m]:
                b = best[m]
                out.append("SUMMARY n=%d layered: %s; least total %s (|L0| %s + extra %s, %s)" % (
                    m + 1, _fmt_stats(agg[m]), b[0], b[1], b[2], b[3]))
        if news:
            seen = set()
            for n, f in sorted(news, key=lambda t: (t[0], len(t[1]))):
                if (n, len(f)) in seen:
                    continue
                seen.add((n, len(f)))
                out.append("NEW n=%d size=%d POINTS %s" % (n, len(f), ",".join(map(str, f))))
        else:
            out.append("no layered set below the records (F_3^5: 53, F_3^6: 105, F_3^4: 27)")
    except Exception:  # noqa: BLE001
        out.append("ERROR " + traceback.format_exc()[-700:])
    out.append("total %.0fs" % (time.time() - t0))
    return " || ".join(out)
