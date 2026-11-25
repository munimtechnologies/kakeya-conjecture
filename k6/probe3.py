# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Third round of searches for small Kakeya sets in F_3^6 and F_3^7.

Records going in: 105 points in F_3^6 (k6/sets/kakeya6_105.json) and 213 in F_3^7.

F_3^6: more walks from the start that produced 105 (F_3 x the 53-point set of Chung, Du and
Wesley) with new seeds and noise levels, and two walks continuing from the 105-point set.

F_3^7: the 'layer' idea.  Put a Kakeya set A of F_3^6 in the layer x_0 = 0; it covers every
direction with d_0 = 0.  A direction (1, d') needs points (0, u), (1, u + d'), (2, u + 2d').
In F_3^6 the same idea with the 53-point set frozen in one layer stopped at
107 = 53 + 2 * 27 points.  If the extra cost in F_3^7 is again about twice the size two
dimensions down (2 * 53), a frozen 105-point layer would give about 105 + 106 = 211 < 213.
Phase 1 keeps the layer frozen, phase 2 frees it and keeps walking from the phase-1 best.
"""
import json
import os
import time
import traceback

try:
    from . import fastwalk as fw
    from . import probe2
except ImportError:  # pragma: no cover
    import fastwalk as fw
    import probe2

HERE = os.path.dirname(os.path.abspath(__file__))
RECORD = {6: "kakeya6_105.json", 7: "kakeya7_213.json"}


def load(name):
    with open(os.path.join(HERE, "sets", name), encoding="utf-8") as f:
        return json.load(f)["points"]


def task(t):
    """t = (n, label, start, frozen, seed, budget, noise, fresh, split).
    split > 0 and frozen: phase 1 (frozen) gets split * budget, phase 2 frees everything."""
    n, label, start, frozen, seed, budget, noise, fresh, split = t
    t0 = time.time()
    try:
        T = fw._tables(n)
        tset = round(time.time() - t0, 1)
        phase1 = None
        if frozen and split:
            left = budget - (time.time() - t0)
            r1 = fw.search(T, seed=seed, budget=max(5.0, left * split), start=start,
                           frozen=frozen, noise=noise, fresh=fresh)
            phase1 = (r1["best"], r1["reached"][-4:], r1["steps"])
            left = budget - (time.time() - t0)
            r = fw.search(T, seed=seed + 1000, budget=max(5.0, left), start=r1["points"],
                          frozen=(), noise=noise, fresh=0)
            if r1["best"] < r["best"]:  # cannot happen (phase 2 starts from r1), kept safe
                r["best"], r["points"] = r1["best"], r1["points"]
        else:
            r = fw.search(T, seed=seed, budget=max(5.0, budget - (time.time() - t0)),
                          start=start, frozen=frozen, noise=noise, fresh=fresh)
        r["phase1"] = phase1
        r["verified"] = fw.independent_check(n, r["points"])
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


def plan(n, seconds, sets):
    full1 = fw.full(1)
    if n == 6:
        rec = load(RECORD[6])
        p53 = fw.product(full1, 1, sets[5], 5) if 5 in sets else fw.classical(6)
        st = [("from105", rec, (), 41, 0.15, 0, 0), ("from105", rec, (), 42, 0.25, 0, 0),
              ("prod53x3", p53, (), 43, 0.15, 0, 0), ("prod53x3", p53, (), 44, 0.18, 0, 0),
              ("prod53x3", p53, (), 45, 0.12, 0, 0)]
    else:
        a105, a106 = load("kakeya6_105.json"), load("kakeya6_106.json")
        st = [("layer105-frozen", None, a105, 51, 0.10, 5, 0.6),
              ("layer105-frozen", None, a105, 52, 0.15, 0, 0.6),
              ("layer106-frozen", None, a106, 53, 0.10, 5, 0.6),
              ("prod105x3", fw.product(full1, 1, a105, 6), (), 54, 0.12, 0, 0),
              ("from213", load(RECORD[7]), (), 55, 0.20, 0, 0)]
    return [(n, lab, start, fr, seed, seconds, noise, fresh, split)
            for (lab, start, fr, seed, noise, fresh, split) in st]


def probe(n, seconds=440, workers=5):
    t0 = time.time()
    out = []
    try:
        rec = load(RECORD[n])
        sets, note = fw.downloads()
        out.append("fetch " + note)
        T = fw.setup(n)
        out.append("RECORD %d" % len(rec))
        tasks = plan(n, seconds - (time.time() - t0) - 25, sets)
        res, mode = run_pool(tasks[:workers], workers)
        out.append("mode=%s cpu=%s" % (mode, os.cpu_count()))
        best = None
        for r in res:
            if "error" in r:
                out.append("%s ERROR %s" % (r["label"], r["error"]))
                continue
            rr = r["reached"]
            short = rr if len(rr) <= 8 else rr[:2] + ["..."] + rr[-5:]
            out.append("%s seed=%d noise=%.2f start=%s phase1=%s: best %d verified=%s reached %s; "
                       "%d steps, %d restarts" % (
                           r["label"], r["seed"], r["noise"], r["start_size"], r["phase1"],
                           r["best"], r["verified"], short, r["steps"], r["restarts"]))
            if r["verified"] and (best is None or r["best"] < best["best"]):
                best = r
        if best is not None:
            out.append("BEST-STATS " + probe2.stats(T, best["points"]))
            if best["best"] < len(rec):
                out.append("NEW n=%d size=%d (record was %d) label=%s POINTS %s" % (
                    n, best["best"], len(rec), best["label"],
                    ",".join(map(str, best["points"]))))
            else:
                out.append("no set below %d (best %d)" % (len(rec), best["best"]))
    except Exception:  # noqa: BLE001
        out.append("ERROR " + traceback.format_exc()[-700:])
    out.append("total %.0fs" % (time.time() - t0))
    return " || ".join(out)
