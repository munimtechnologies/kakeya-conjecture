# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Report-only probe around the 422-point Kakeya set of F_3^8 (k6/sets/kakeya8_422.json).

1. REM: full lines, single-line directions and deletable points of the 422-point set.
2. ST: how its two free layers A (x0 = 1, 96 points) and B (x0 = 2, 115 points) sit
   relative to its layer K7 (x0 = 0, the 211-point set of F_3^7); the same for the
   211-point set (layers 105 / 56 / 50).
3. More layered searches in F_3^8 on the 211-point layer (fresh starts, and from the 422
   set's own layers).
4. Layered sets in F_3^9 with the 422-point set K8 fixed in the layer x0 = 0
   (classical 1023, Lev's pattern 853).  Starts:
   * "lift": L1 = F_3 x A, L2 = F_3 x B.  It is already a Kakeya set: for v = (c, w)
     take x = (0, y) with y in K7, y + w in A, y - w in B (the F_3^8 layer property);
     then x + v = (c, y + w) is in L1 and x - v = (-c, y - w) is in L2.  1055 points.
   * "prod": L1 = L2 = K8, i.e. F_3 x K8, 1266 points.
Point lists are delta-coded as in k8/probe8.py; for F_3^9 an internal point p stands for
the point with index 3^8 + p."""
import collections
import itertools
import json
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

SETS = os.path.join(_ROOT, "k6", "sets")


def load_full(name):
    with open(os.path.join(SETS, name), encoding="utf-8") as f:
        d = json.load(f)
    if "points" in d:
        return sorted(d["points"])
    with open(os.path.join(_ROOT, d["layer0_file"]), encoding="utf-8") as f:
        base = json.load(f)["points"]
    assert decode(d["points_code"]) == d["free_points"]
    return sorted(base + [d["offset"] + p for p in d["free_points"]])


def split(n, idx):
    """Layers (as integer indices in F_3^(n-1)) of a set of F_3^n."""
    N = 3 ** (n - 1)
    L = [set(), set(), set()]
    for i in idx:
        L[i // N].add(i % N)
    return L


def neg(m, i):
    return check_index(m, tuple((-c) % 3 for c in check.decode(m, i)))


def check_index(m, x):
    return sum(c * 3 ** (m - 1 - k) for k, c in enumerate(x))


def add(m, i, j, k=1):
    a, b = check.decode(m, i), check.decode(m, j)
    return check_index(m, tuple((x + k * y) % 3 for x, y in zip(a, b)))


def structure(n, name):
    m = n - 1
    K, A, B = split(n, load_full(name))
    nB = set(neg(m, b) for b in B)
    nA = set(neg(m, a) for a in A)
    U = A | B
    ukak = check.is_kakeya(m, sorted(U))
    # best overlap of A u B with a translate of K and of -K
    tK = collections.Counter(add(m, u, k, 2) for u in U for k in K)   # t = u - k
    nK = [neg(m, k) for k in K]
    tnK = collections.Counter(add(m, u, k, 2) for u in U for k in nK)
    tAB = collections.Counter(add(m, a, b, 2) for a in A for b in B)  # B + t = A
    tAnB = collections.Counter(add(m, a, b, 2) for a in A for b in nB)
    return ("ST%d K%d A%d B%d AB%d AK%d BK%d A-B%d A-A%d uKak=%s tK%d t-K%d tBA%d t-BA%d" % (
        n, len(K), len(A), len(B), len(A & B), len(A & K), len(B & K), len(A & nB),
        len(A & nA), ukak, max(tK.values()), max(tnK.values()), max(tAB.values()),
        max(tAnB.values())))


def rem422(deadline):
    out = {"kind": "REM"}
    idx = load_full("kakeya8_422.json")
    lines, single, rem, unc = removable_of(8, idx)
    out.update(lines=lines, single=single, removable=rem, unc=unc)
    cur, dropped, rng = list(idx), [], random.Random(9)
    while rem and not unc and time.time() < deadline - 150:
        p = rng.choice(rem)
        cur.remove(p)
        dropped.append(p)
        _, _, rem, unc = removable_of(8, cur)
    out["dropped"] = dropped
    out["final"] = len(cur)
    out["final_ok"] = check.is_kakeya(8, cur) if dropped else None
    out["final_code"] = encode(cur) if dropped else ""
    try:
        out["st"] = structure(8, "kakeya8_422.json") + " | " + structure(7, "kakeya7_211.json")
    except Exception:  # noqa: BLE001
        out["st"] = "ST ERROR " + traceback.format_exc()[-300:]
    try:
        out["search"] = l8_search(16, 0.10, None, deadline)
    except Exception:  # noqa: BLE001
        out["search"] = {"kind": "L8", "error": traceback.format_exc()[-300:]}
    return out


def l8_search(seed, noise, start_mode, deadline):
    W = layer.Layer(7, load_full("kakeya7_211.json"))
    start = None
    if start_mode == "422":
        start = [f - W.N for f in load_full("kakeya8_422.json") if f >= W.N]
    budget = max(5.0, deadline - time.time() - 12)
    r = layer.search(W, seed=seed, budget=budget, noise=noise, start=start)
    full = W.full_indices(r["points"])
    return {"kind": "L8", "seed": seed, "noise": noise, "start": start_mode or "-",
            "best": r["best"], "total": len(full), "reached": r["reached"][-3:],
            "steps": r["steps"], "verified": check.is_kakeya(8, full),
            "code": encode(r["points"])}


def l9_search(seed, noise, start_mode, deadline):
    K8 = load_full("kakeya8_422.json")
    W = layer.Layer(8, K8)
    N = W.N
    K, A, B = split(8, K8)
    if start_mode == "lift":
        L1 = [c * 2187 + a for c in range(3) for a in A]
        L2 = [c * 2187 + b for c in range(3) for b in B]
    else:  # "prod"
        L1, L2 = list(K8), list(K8)
    start = sorted(set(L1)) + [N + z for z in sorted(set(L2))]
    # the start must already be a Kakeya set (checked by the layer state: nothing uncovered)
    W.load(start)
    start_unc = len(W.unc)
    budget = max(5.0, deadline - time.time() - 75)
    r = layer.search(W, seed=seed, budget=budget, noise=noise, start=start)
    full = W.full_indices(r["points"])
    ok = check.is_kakeya(9, full) if len(full) < 1023 else None
    return {"kind": "L9", "seed": seed, "noise": noise, "start": "%s/u%d" % (start_mode, start_unc),
            "best": r["best"], "total": len(full), "reached": r["reached"][-3:],
            "steps": r["steps"], "verified": ok, "code": encode(r["points"])}


def run_task(t):
    kind, args = t[0], t[1:]
    try:
        if kind == "REM":
            return rem422(*args)
        if kind == "L8":
            return l8_search(*args)
        if kind == "L9":
            return l9_search(*args)
        return {"kind": kind, "error": "unknown task"}
    except Exception:  # noqa: BLE001 - the report is the product
        return {"kind": kind, "error": traceback.format_exc()[-400:]}


def _short(r):
    if "error" in r:
        return "%s ERROR %s" % (r.get("kind"), r["error"])
    return "%s s%d nz%.2f st%s: x%d t%d last %s st%d ok=%s" % (
        r["kind"], r["seed"], r["noise"], r["start"], r["best"], r["total"],
        r["reached"][-2:], r["steps"], r["verified"])


def probe(seconds=470, workers=7):
    t0 = time.time()
    deadline = t0 + seconds - 20
    out, tail = [], []
    try:
        tasks = [("REM", deadline),
                 ("L8", 17, 0.10, None, deadline),
                 ("L8", 18, 0.10, None, deadline),
                 ("L8", 19, 0.07, "422", deadline),
                 ("L9", 31, 0.05, "lift", deadline),
                 ("L9", 32, 0.10, "lift", deadline),
                 ("L9", 33, 0.10, "prod", deadline)][:workers]
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(workers) as pool:
            res = pool.map(run_task, tasks, chunksize=1)
        out.append("pool%d cpu=%s" % (workers, os.cpu_count()))
        best8, best9, rem = None, None, None
        for r in res:
            if r.get("kind") == "REM" and "error" not in r:
                rem = r
                res.append(r["search"])   # handled when the loop reaches it
                continue
            out.append(_short(r))
            if "error" in r or not r["verified"]:
                continue
            if r["kind"] == "L8" and (best8 is None or r["total"] < best8["total"]):
                best8 = r
            if r["kind"] == "L9" and (best9 is None or r["total"] < best9["total"]):
                best9 = r
        if rem is not None:
            tail.append(rem["st"])
            tail.append("REM422 lines=%d single=%d unc=%d removable=%s dropped=%s final=%d ok=%s %s" % (
                rem["lines"], rem["single"], rem["unc"], rem["removable"][:20], rem["dropped"],
                rem["final"], rem["final_ok"],
                ("CODE " + rem["final_code"]) if rem["final_code"] else ""))
        else:
            tail.append("REM422 failed")
        if best8 is not None and best8["total"] < 422:
            tail.append("NEW N8 total=%d (211 + %d) seed %d CODE %s" % (
                best8["total"], best8["best"], best8["seed"], best8["code"]))
        else:
            tail.append("N8 best %s, nothing below 422" % (best8["total"] if best8 else None))
        if best9 is not None:
            tail.append("N9 total=%d (422 + %d; classical 1023, pattern 853) seed %d CODE %s" % (
                best9["total"], best9["best"], best9["seed"], best9["code"]))
        else:
            tail.append("N9 no verified set below 1023")
    except Exception:  # noqa: BLE001
        tail.append("ERROR " + traceback.format_exc()[-700:])
    tail.append("total %.0fs" % (time.time() - t0))
    return " || ".join(out + tail)
