# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Report-only probe: is the doubling k_(n+1) <= 2 k_n of the layered records a 'split'?

Background.  K = {0} x L0 u {1} x L1 u {2} x L2 with L0 a Kakeya set of F_3^m is a Kakeya set
of F_3^(m+1) iff every v in F_3^m has some x in L0 with x + v in L1 and x - v in L2.  The
starts L1 = L0 + t, L2 = L0 - t always work (k8/probe10.py), at cost 2 |L0|, and the records
211 = 105 + 106, 422 = 211 + 211, 844 = 422 + 422, 1695 = 844 + 851 have |L1| + |L2| close
to |L0|.  If the free layers were always about half of L0 + t and half of L0 - t, a proof of
k_(m+1) <= 2 k_m (+ small) might follow from a rule for the split.

1. OV (main process): for each stored layered record, the largest overlaps
   max_t |L1 cap (L0 + t)| and max_t |L2 cap (L0 + t)| with their t, and whether the two
   best t are negatives of each other (as in the shift start).
2. SPLIT searches: the layered search of k8/layer.py restricted to points of the start
   L1 = L0 + t, L2 = L0 - t (t = 0 is the product start).  Same WalkSAT moves, but a line may
   only be completed with allowed points.  n = 7, 8 (sequential worker), 9 and 10.
3. More layered searches in F_3^10 on the 844-point layer: continued from the 1695-point
   set's own free layers ('own'), and fresh 'prod'/'shift' starts.

Point lists are delta-coded as in k8/probe8.py; an internal point p of the free layers of
F_3^n stands for the point with index 3^(n-1) + p.  Key lines come last in the report."""
import os
import random
import sys
import time
import traceback

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from k8 import check, layer  # noqa: E402
from k8.probe8 import encode, decode  # noqa: E402
from k8.probe10 import load_any  # noqa: E402

BASE = {7: "kakeya6_105.json", 8: "kakeya7_211.json", 9: "kakeya8_422.json",
        10: "kakeya9_844.json"}
OWN = {9: "kakeya9_844.json", 10: "kakeya10_1695.json"}
RECORD = {7: 211, 8: 422, 9: 844, 10: 1695}
CLASSICAL = {7: 255, 8: 511, 9: 1023, 10: 2047}
LAYERED = [("211", 7, "kakeya7_211.json"), ("422", 8, "kakeya8_422.json"),
           ("844", 9, "kakeya9_844.json"), ("846", 9, "kakeya9_846.json"),
           ("1702", 10, "kakeya10_1702.json"), ("1695", 10, "kakeya10_1695.json")]


# ---- 1. translate overlaps ---------------------------------------------------------------
def overlaps(n, name):
    idx = load_any(name)
    m = n - 1
    N = 3 ** m
    L = [sorted(i % N for i in idx if i // N == c) for c in range(3)]
    add, NEG = layer.arith(m)
    res = []
    for c in (1, 2):
        cnt = {}
        for y in L[c]:
            for x in L[0]:
                t = add(y, NEG[x])           # y = x + t
                cnt[t] = cnt.get(t, 0) + 1
        t, best = max(cnt.items(), key=lambda kv: (kv[1], -kv[0]))
        res.append((best, t, cnt.get(0, 0)))
    (a, t1, a0), (b, t2, b0) = res
    neg = NEG[t1] == t2
    # overlap of L1 with L0 + t1 and of L2 with L0 - t1 for the SAME t (the shift pattern)
    s1 = set(add(x, t1) for x in L[0])
    s2 = set(add(x, NEG[t1]) for x in L[0])
    pair = (sum(1 for y in L[1] if y in s1), sum(1 for z in L[2] if z in s2))
    return "%s:%d/%d/%d t1=%d max1=%d t2=%d max2=%d neg=%s at0=%d,%d pair=%d,%d" % (
        name.split("_")[1].split(".")[0], len(L[0]), len(L[1]), len(L[2]), t1, a, t2, b,
        neg, a0, b0, pair[0], pair[1])


# ---- 2. restricted search ----------------------------------------------------------------
def search_r(W, allowed, seed=1, budget=60.0, noise=0.1, tenure=5, stall=5000, ncand=8,
             start=None):
    """layer.search with every added point taken from `allowed` (bytearray over 2N).
    The start must already cover every vector (it does for the prod/shift starts)."""
    t0 = time.time()
    rng = random.Random(seed)
    W.reset()
    for p in start:
        W.add(p)
    if W.unc:
        raise ValueError("start does not cover %d vectors" % len(W.unc))
    W.prune(rng)
    best = sorted(W.S)
    reached = [(len(best), round(time.time() - t0, 1))]
    traj = list(best)
    P, L, N = 2 * W.N, W.L, W.N
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
    while time.time() - t0 < budget:
        st["it"] += 1
        it = st["it"]
        u = len(W.unc)
        if u == 0:
            if len(W.S) < len(traj):
                traj = sorted(W.S)
            if len(W.S) < len(best):
                best = sorted(W.S)
                reached.append((len(best), round(time.time() - t0, 1)))
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
        inS, YV, ZV = W.inS, W.YV, W.ZV
        b = v * L
        mn, opts, near = 3, [], []
        for i in range(L):
            py, pz = YV[b + i], N + ZV[b + i]
            if not (allowed[py] and allowed[pz]):
                continue
            k = (not inS[py]) + (not inS[pz])
            if k < mn:
                near = opts if k == mn - 1 else []
                mn, opts = k, [i]
            elif k == mn:
                opts.append(i)
            elif k == mn + 1:
                near.append(i)
        if not opts:
            raise ValueError("vector %d has no allowed line" % v)
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


# ---- searches ----------------------------------------------------------------------------
def start_of(W, K, mode, rng):
    N = W.N
    if mode == "prod":
        t = 0
    elif mode == "shift":
        t = rng.randrange(1, N)
    else:
        raise ValueError(mode)
    add, NEG = layer.arith(W.m)
    L1 = set(add(k, t) for k in K)
    L2 = set(add(k, NEG[t]) for k in K)
    return sorted(L1) + [N + z for z in sorted(L2)]


def lsearch(n, seed, noise, mode, deadline, margin, restrict=False):
    t0 = time.time()
    m = n - 1
    K = [i for i in load_any(BASE[n])]
    W = layer.Layer(m, K)
    setup = round(time.time() - t0, 1)
    rng = random.Random(1000 + seed)
    if mode == "own":
        N = 3 ** m
        start = [f - N for f in load_any(OWN[n]) if f >= N]
    else:
        start = start_of(W, K, mode, rng)
    S0 = set(start)
    budget = max(5.0, deadline - time.time() - margin)
    if restrict:
        allowed = bytearray(2 * W.N)
        for p in start:
            allowed[p] = 1
        r = search_r(W, allowed, seed=seed, budget=budget, noise=noise, start=start)
    else:
        r = layer.search(W, seed=seed, budget=budget, noise=noise, start=start)
    full = W.full_indices(r["points"])
    ok = check.is_kakeya(n, full) if len(full) < CLASSICAL[n] else None
    l1 = sum(1 for p in r["points"] if p < W.N)
    kept = sum(1 for p in r["points"] if p in S0)
    return {"kind": ("S%d" if restrict else "L%d") % n, "n": n, "seed": seed, "noise": noise,
            "start": mode, "best": r["best"], "total": len(full), "split": (l1, r["best"] - l1),
            "kept": kept, "reached": r["reached"][-2:], "steps": r["steps"], "verified": ok,
            "code": encode(r["points"]), "setup": setup}


def seq_small(deadline):
    """Split searches in F_3^7 (on the 105-point layer) and F_3^8 (on the 211-point layer)."""
    mid = time.time() + (deadline - time.time()) * 0.35
    a = lsearch(7, 81, 0.10, "prod", mid, 5, restrict=True)
    b = lsearch(8, 82, 0.10, "prod", deadline, 15, restrict=True)
    return {"kind": "SEQ", "parts": [a, b]}


def run_task(t):
    kind, args = t[0], t[1:]
    try:
        if kind == "SEQ":
            return seq_small(*args)
        return lsearch(*args)
    except Exception:  # noqa: BLE001 - the report is the product
        return {"kind": kind, "error": traceback.format_exc()[-300:]}


def _short(r):
    if "error" in r:
        return "%s ERROR %s" % (r.get("kind"), r["error"])
    return "%s s%d %.2f %s x%d t%d %d+%d k%d %s st%d ok=%s" % (
        r["kind"], r["seed"], r["noise"], r["start"], r["best"], r["total"], r["split"][0],
        r["split"][1], r["kept"], r["reached"], r["steps"], r["verified"])


def probe(seconds=560, workers=8):
    t0 = time.time()
    deadline = t0 + seconds - 20
    out, tail = [], []
    try:
        sample = [0, 1, 2, 40, 41, 1000, 1031, 1032, 39000]
        assert decode(encode(sample)) == sample, "codec"
        tasks = [("L10", 10, 71, 0.05, "own", deadline, 85),
                 ("L10", 10, 72, 0.07, "own", deadline, 85),
                 ("L10", 10, 73, 0.07, "prod", deadline, 85),
                 ("L10", 10, 74, 0.06, "prod", deadline, 85),
                 ("L10", 10, 75, 0.07, "shift", deadline, 85),
                 ("S10", 10, 76, 0.07, "prod", deadline, 85, True),
                 ("S9", 9, 77, 0.07, "prod", deadline, 40, True),
                 ("SEQ", deadline)][:workers]
        import multiprocessing as mp
        pool = mp.get_context("spawn").Pool(workers)
        try:
            pending = pool.map_async(run_task, tasks, chunksize=1)
            ov = []
            for label, n, name in LAYERED:   # main process works while the pool searches
                try:
                    ov.append(overlaps(n, name))
                except Exception as e:  # noqa: BLE001
                    ov.append("%s ERR %s" % (label, str(e)[:60]))
            res = pending.get()
        finally:
            pool.close()
            pool.join()
        out.append("pool%d cpu=%s" % (workers, os.cpu_count()))
        flat = []
        for r in res:
            if r.get("kind") == "SEQ" and "parts" in r:
                flat.extend(r["parts"])
            else:
                flat.append(r)
        best = {}
        for r in flat:
            out.append(_short(r))
            if "error" in r or not r["verified"]:
                continue
            n = r["n"]
            if n not in best or r["total"] < best[n]["total"]:
                best[n] = r
        tail.append("OV " + " | ".join(ov))
        totals = sorted((r["total"], r["kind"]) for r in flat if "error" not in r)
        tail.append("TOTALS %s" % totals)
        for n in (7, 8, 9, 10):
            b = best.get(n)
            if b is not None and b["total"] < RECORD[n]:
                tail.append("NEW N%d total=%d (%s + %d = %d+%d; record %d) %s seed %d CODE %s" % (
                    n, b["total"], BASE[n], b["best"], b["split"][0], b["split"][1], RECORD[n],
                    b["kind"], b["seed"], b["code"]))
            else:
                tail.append("N%d best %s, nothing below %d" % (
                    n, b["total"] if b else None, RECORD[n]))
    except Exception:  # noqa: BLE001
        tail.append("ERROR " + traceback.format_exc()[-700:])
    tail.append("total %.0fs" % (time.time() - t0))
    return " || ".join(out + tail)
