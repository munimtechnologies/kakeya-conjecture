# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Second round of searches below the record sets in k6/sets/ (106 points in F_3^6,
213 in F_3^7), plus structure statistics of the record sets.

Runs five k6.fastwalk searches in parallel: two continue from the record set with new
seeds, three start from other structured sets (products and the classical set) with new
seeds.  Reports, for the input record and for the best set found: number of lines, number
of directions covered by exactly one line, number of points whose removal keeps the set
Kakeya (0 means minimal under deletion), and the smallest and largest hyperplane sections.
"""
import json
import os
import time
import traceback

try:
    from . import fastwalk as fw
except ImportError:  # pragma: no cover
    import fastwalk as fw

HERE = os.path.dirname(os.path.abspath(__file__))
RECORD = {6: "kakeya6_106.json", 7: "kakeya7_213.json"}


def load(n):
    with open(os.path.join(HERE, "sets", RECORD[n]), encoding="utf-8") as f:
        return json.load(f)["points"]


def stats(T, K):
    """Structure numbers of a set K (indices) using the fastwalk tables."""
    n, D = T["n"], T["D"]
    W = fw.FW(T)
    W.load(K)
    lines = sum(len(S) for S in W.lines)
    single = sum(1 for S in W.lines if len(S) == 1)
    removable = sum(1 for p in W.K if W.cost[p] == 0)
    pts = [T["pts"][p] for p in K]
    secmin, secmax, prof = 10 ** 9, -1, {}
    for j in range(D):
        xi = T["pts"][T["dirs"][j]]
        c = [0, 0, 0]
        for x in pts:
            c[sum(a * b for a, b in zip(x, xi)) % 3] += 1
        c.sort()
        secmin, secmax = min(secmin, c[0]), max(secmax, c[2])
        prof[tuple(c)] = prof.get(tuple(c), 0) + 1
    neg = set(tuple((-a) % 3 for a in x) for x in pts)
    top = sorted(prof.items(), key=lambda kv: -kv[1])[:4]
    return ("size %d uncovered %d lines %d single-line-dirs %d removable %d sections min %d "
            "max %d profiles %d top %s |K cap -K| %d" % (
                len(K), len(W.unc), lines, single, removable, secmin, secmax, len(prof),
                top, len(neg & set(pts))))


def plan(n, seconds, sets, rec):
    full1 = fw.full(1)
    if n == 6:
        p13 = fw.product(sets[3], 3, sets[3], 3) if 3 in sets else None
        p53 = fw.product(full1, 1, sets[5], 5) if 5 in sets else fw.classical(6)
        st = [("from106", rec, 21, 0.10), ("from106", rec, 22, 0.20),
              ("prod13x13", p13, 23, 0.10), ("classical127", fw.classical(6), 24, 0.05),
              ("prod53x3", p53, 25, 0.15)]
    else:
        p106 = fw.product(full1, 1, load(6), 6)
        st = [("from213", rec, 31, 0.10), ("from213", rec, 32, 0.05),
              ("prod106x3", p106, 33, 0.10), ("prod106x3", p106, 34, 0.15),
              ("classical255", fw.classical(7), 35, 0.10)]
    return [(n, lab, start, (), seed, seconds, noise, 0) for (lab, start, seed, noise) in st]


def probe(n, seconds=440, workers=5):
    t0 = time.time()
    out = []
    try:
        rec = load(n)
        sets, note = fw.downloads()
        out.append("fetch " + note)
        T = fw.setup(n)
        out.append("RECORD " + stats(T, rec))
        tasks = plan(n, seconds - (time.time() - t0) - 25, sets, rec)
        res, mode = fw.run_pool(tasks[:workers], workers)
        out.append("mode=%s cpu=%s" % (mode, os.cpu_count()))
        best = None
        for r in res:
            if "error" in r:
                out.append("%s ERROR %s" % (r["label"], r["error"]))
                continue
            rr = r["reached"]
            short = rr if len(rr) <= 8 else rr[:2] + ["..."] + rr[-5:]
            out.append("%s seed=%d noise=%.2f start=%s: best %d verified=%s reached %s; "
                       "%d steps, %d restarts" % (
                           r["label"], r["seed"], r["noise"], r["start_size"], r["best"],
                           r["verified"], short, r["steps"], r["restarts"]))
            if r["verified"] and (best is None or r["best"] < best["best"]):
                best = r
        if best is not None:
            out.append("BEST-STATS " + stats(T, best["points"]))
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
