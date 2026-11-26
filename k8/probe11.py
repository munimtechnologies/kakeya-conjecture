# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Report-only probe around the 844-point Kakeya set of F_3^9 (k6/sets/kakeya9_844.json).

1. REM: layer sizes, full lines, single-line directions and deletable points of the
   844-point set, then one more layered search in F_3^9.
2. One more layered search in F_3^9 with the 422-point set K8 fixed in the layer x0 = 0.
3. Layered searches in F_3^10 with the 844-point set K9 fixed in the layer x0 = 0, given
   more time than in k8/probe10.py (record 1702 on the 846-point layer; classical 2047,
   Lev's pattern 1707, 2 * 844 = 1688).

Starts as in k8/probe10.py ("prod": F_3 x K, "shift": K + t, K - t; both proved Kakeya).
Point lists are delta-coded as in k8/probe8.py; for F_3^n with fixed layer in F_3^(n-1)
an internal point p stands for the point with index 3^(n-1) + p.  The report puts its key
lines last (the gate keeps only the end of a long report)."""
import os
import random
import sys
import time
import traceback

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from k8 import check, layer  # noqa: E402
from k8.probe8 import encode, decode, removable_of  # noqa: E402
from k8.probe10 import load_any, start_points  # noqa: E402

CLASSICAL = {9: 1023, 10: 2047}
BASE = {9: "kakeya8_422.json", 10: "kakeya9_844.json"}
RECORD = {9: 844, 10: 1702}


def lsearch(n, seed, noise, mode, deadline, margin):
    t0 = time.time()
    m = n - 1
    K = load_any(BASE[n])
    W = layer.Layer(m, K)
    setup = round(time.time() - t0, 1)
    rng = random.Random(1000 + seed)
    start = start_points(W, K, mode, rng)
    W.load(start)
    start_unc = len(W.unc)
    budget = max(5.0, deadline - time.time() - margin)
    r = layer.search(W, seed=seed, budget=budget, noise=noise, start=start)
    full = W.full_indices(r["points"])
    ok = check.is_kakeya(n, full) if len(full) < CLASSICAL[n] else None
    l1 = sum(1 for p in r["points"] if p < W.N)
    return {"kind": "L%d" % n, "n": n, "seed": seed, "noise": noise, "layer0": len(K),
            "start": "%s/u%d" % (mode, start_unc), "best": r["best"], "total": len(full),
            "split": (l1, r["best"] - l1), "reached": r["reached"][-3:], "first": r["reached"][:2],
            "steps": r["steps"], "verified": ok, "code": encode(r["points"]), "setup": setup}


def rem844(deadline):
    out = {"kind": "REM"}
    idx = load_any("kakeya9_844.json")
    N = 3 ** 8
    out["layers"] = [sum(1 for i in idx if i // N == c) for c in range(3)]
    out["ok844"] = check.is_kakeya(9, idx)
    lines, single, rem, unc = removable_of(9, idx)
    out.update(lines=lines, single=single, removable=rem, unc=unc)
    try:
        out["search"] = lsearch(9, 50, 0.10, "prod", deadline, 40)
    except Exception:  # noqa: BLE001
        out["search"] = {"kind": "L9", "error": traceback.format_exc()[-300:]}
    return out


def run_task(t):
    kind, args = t[0], t[1:]
    try:
        if kind == "REM":
            return rem844(*args)
        return lsearch(*args)
    except Exception:  # noqa: BLE001 - the report is the product
        return {"kind": kind, "error": traceback.format_exc()[-400:]}


def _short(r):
    if "error" in r:
        return "%s ERROR %s" % (r.get("kind"), r["error"])
    return "%s s%d nz%.2f st%s: x%d t%d (%d+%d) last %s st%d su%s ok=%s" % (
        r["kind"], r["seed"], r["noise"], r["start"], r["best"], r["total"], r["split"][0],
        r["split"][1], r["reached"][-2:], r["steps"], r["setup"], r["verified"])


def probe(seconds=550, workers=8):
    t0 = time.time()
    deadline = t0 + seconds - 20
    out, tail = [], []
    try:
        sample = [0, 1, 2, 40, 41, 1000, 1031, 1032, 39000]
        assert decode(encode(sample)) == sample, "codec"
        tasks = [("REM", deadline),
                 ("L9", 9, 51, 0.07, "shift", deadline, 40),
                 ("L10", 10, 61, 0.07, "prod", deadline, 100),
                 ("L10", 10, 62, 0.05, "prod", deadline, 100),
                 ("L10", 10, 63, 0.07, "shift", deadline, 100),
                 ("L10", 10, 64, 0.10, "prod", deadline, 100),
                 ("L10", 10, 65, 0.05, "shift", deadline, 100),
                 ("L10", 10, 66, 0.07, "prod", deadline, 100)][:workers]
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(workers) as pool:
            res = pool.map(run_task, tasks, chunksize=1)
        out.append("pool%d cpu=%s" % (workers, os.cpu_count()))
        best = {9: None, 10: None}
        rem = None
        for r in res:
            if r.get("kind") == "REM" and "error" not in r:
                rem = r
                res.append(r["search"])   # handled when the loop reaches it
                continue
            out.append(_short(r))
            if "error" in r or not r["verified"]:
                continue
            n = r["n"]
            if best[n] is None or r["total"] < best[n]["total"]:
                best[n] = r
        if rem is not None:
            tail.append("REM844 layers=%s ok=%s lines=%d single=%d unc=%d removable=%s" % (
                rem["layers"], rem["ok844"], rem["lines"], rem["single"], rem["unc"],
                rem["removable"][:20]))
        else:
            tail.append("REM844 failed")
        totals = sorted((r["total"], r["n"]) for r in res if "error" not in r and "total" in r)
        tail.append("TOTALS %s" % totals)
        b9 = best[9]
        if b9 is not None and b9["total"] < RECORD[9]:
            tail.append("NEW N9 total=%d (422 + %d = %d+%d) seed %d CODE %s" % (
                b9["total"], b9["best"], b9["split"][0], b9["split"][1], b9["seed"], b9["code"]))
        else:
            tail.append("N9 best %s, nothing below 844" % (b9["total"] if b9 else None))
        b10 = best[10]
        if b10 is not None and b10["total"] < RECORD[10]:
            tail.append("NEW N10 total=%d (layer kakeya9_844 + %d = %d+%d; record 1702, pattern "
                        "1707, 2*844 = 1688) seed %d CODE %s" % (
                            b10["total"], b10["best"], b10["split"][0], b10["split"][1],
                            b10["seed"], b10["code"]))
        else:
            tail.append("N10 best %s on the 844 layer, nothing below 1702" % (
                b10["total"] if b10 else None))
    except Exception:  # noqa: BLE001
        tail.append("ERROR " + traceback.format_exc()[-700:])
    tail.append("total %.0fs" % (time.time() - t0))
    return " || ".join(out + tail)
