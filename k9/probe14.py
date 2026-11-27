# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Report-only probe: layered searches started from a doubled set.

The excess-1 split of the 105-point set (k9/splits/kakeya6_105_split.json), lifted by
Lemma 2 of k9/SPLIT.md through 211 -> 422 -> 844, gives a split (A, B) of the 844-point
Kakeya set K9 of F_3^9 with |A| + |B| = 845.  Lemma 1 puts A and B in the two free layers:
a Kakeya set of 1689 points in F_3^10.  Here the layered search of k8/layer.py keeps K9
fixed in the layer x0 = 0 and starts from exactly these free layers.  A result with fewer
than 845 free points is a Kakeya set of F_3^10 below 1689, and by Lemma 2 it inherits a split
of excess 1, so it lowers the whole doubling family.  One worker does the same one level
down (F_3^9 on the 422-point layer, start 423 free points; the record there has 422) as a
calibration.  Each result is re-checked from the layer criterion with coordinate tuples.
Key lines come last in the report."""
import json
import os
import sys
import time
import traceback

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from k8 import layer  # noqa: E402
from k9 import split  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
T32 = "abcdefghijklmnopqrstuvwxyz012345"
C32 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ6789+/"
CHAIN = {8: ("kakeya7_211.json", "kakeya8_422.json"),
         9: ("kakeya7_211.json", "kakeya8_422.json", "kakeya9_844.json")}


def encode(pts):
    """Gap code of k6/sets (sorted values, gaps from -1 in base 32)."""
    out, prev = [], -1
    for p in sorted(pts):
        d, prev = p - prev, p
        digs = []
        while True:
            digs.append(d % 32)
            d //= 32
            if d == 0:
                break
        digs.reverse()
        out.append("".join(C32[x] for x in digs[:-1]) + T32[digs[-1]])
    return "".join(out)


def decode(s):
    out, prev, val = [], -1, 0
    for ch in s:
        i = C32.find(ch)
        if i >= 0:
            val = val * 32 + i
            continue
        val = val * 32 + T32.index(ch)
        prev += val
        out.append(prev)
        val = 0
    return out


def lifted_split(level):
    """(K, A, B): the record of F_3^level and the split lifted from the 105-point set."""
    with open(os.path.join(HERE, "splits", "kakeya6_105_split.json"), encoding="utf-8") as f:
        d = json.load(f)
    m, K = split.load_set(d["set_file"])
    A, B = split.split_from_states(K, split.decode_split(d["states_code"], len(K)))
    for name in CHAIN[level]:
        n, Kn = split.load_set(name)
        A, B = split.lift(m, Kn, A, B)
        K, m = Kn, n
    assert m == level
    return K, A, B


def layer_check(m, K, pts):
    """Layer criterion with tuples: every v of F_3^m has x in K, x + v in L1, x + 2v in L2."""
    N = 3 ** m
    L0 = [split.digits(x, m) for x in K]
    L1 = [split.digits(p, m) for p in pts if p < N]
    L2 = set(split.digits(p - N, m) for p in pts if p >= N)
    cov = set()
    for x in L0:
        for y in L1:
            if tuple((-a - b) % 3 for a, b in zip(x, y)) in L2:
                cov.add(tuple((b - a) % 3 for a, b in zip(x, y)))
    return len(cov) == N


def run_task(t):
    level, seed, noise, deadline, margin = t
    try:
        t0 = time.time()
        K, A, B = lifted_split(level)
        N = 3 ** level
        start = sorted(A) + [N + b for b in sorted(B)]
        start_ok = layer_check(level, K, start)
        W = layer.Layer(level, K)
        setup = round(time.time() - t0, 1)
        budget = max(5.0, deadline - time.time() - margin)
        r = layer.search(W, seed=seed, budget=budget, noise=noise, start=start)
        pts = r["points"]
        ok = layer_check(level, K, pts)
        l1 = sum(1 for p in pts if p < N)
        kept = len(set(pts) & set(start))
        return {"n": level + 1, "seed": seed, "noise": noise, "layer": len(K),
                "start": len(start), "start_ok": start_ok, "best": r["best"],
                "split": (l1, r["best"] - l1), "kept": kept, "reached": r["reached"][-2:],
                "steps": r["steps"], "restarts": r["restarts"], "ok": ok,
                "code": encode(pts), "setup": setup}
    except Exception:  # noqa: BLE001 - the report is the product
        return {"n": level + 1, "seed": seed, "error": traceback.format_exc()[-300:]}


def _short(r):
    if "error" in r:
        return "N%d s%d ERROR %s" % (r["n"], r["seed"], r["error"])
    return "N%d s%d %.2f L%d st%d/%s x%d %d+%d k%d %s steps%d r%d ok=%s su%.0f" % (
        r["n"], r["seed"], r["noise"], r["layer"], r["start"], r["start_ok"], r["best"],
        r["split"][0], r["split"][1], r["kept"], r["reached"], r["steps"], r["restarts"],
        r["ok"], r["setup"])


def probe(seconds=540, workers=8):
    t0 = time.time()
    deadline = t0 + seconds - 20
    out, tail = [], []
    try:
        sample = [0, 1, 2, 40, 41, 1000, 1031, 1032, 39000]
        assert decode(encode(sample)) == sample, "codec"
        tasks = [(9, 101, 0.05, deadline, 30), (9, 102, 0.06, deadline, 30),
                 (9, 103, 0.07, deadline, 30), (9, 104, 0.04, deadline, 30),
                 (9, 105, 0.06, deadline, 30), (9, 106, 0.08, deadline, 30),
                 (9, 107, 0.05, deadline, 30), (8, 108, 0.07, deadline, 20)][:workers]
        import multiprocessing as mp
        pool = mp.get_context("spawn").Pool(workers)
        try:
            res = pool.map(run_task, tasks, chunksize=1)
        finally:
            pool.close()
            pool.join()
        out.append("pool%d cpu=%s" % (workers, os.cpu_count()))
        best = None
        for r in res:
            out.append(_short(r))
            if "error" in r or not r["ok"] or r["n"] != 10:
                continue
            if best is None or r["best"] < best["best"]:
                best = r
        tot = sorted((r["layer"] + r["best"], r["n"]) for r in res if "error" not in r)
        tail.append("TOTALS %s" % tot)
        if best is not None and best["best"] < 845:
            tail.append("NEW N10 total=%d (844 + %d = %d+%d; record 1689) seed %d CODE %s" % (
                844 + best["best"], best["best"], best["split"][0], best["split"][1],
                best["seed"], best["code"]))
        else:
            tail.append("N10 best %s, nothing below 1689" % (
                844 + best["best"] if best else None))
    except Exception:  # noqa: BLE001
        tail.append("ERROR " + traceback.format_exc()[-700:])
    tail.append("total %.0fs" % (time.time() - t0))
    return " || ".join(out + tail)
