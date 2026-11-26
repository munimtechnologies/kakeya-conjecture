# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Report-only probe around the 427-point Kakeya set of F_3^8 (k6/sets/kakeya8_427.json).

1. Structure of the 427-point set: number of full lines, directions covered by one line
   only, and points whose deletion keeps it Kakeya (greedy deletion if there are any).
2. Layered searches in F_3^8 (k8/layer.py) with the 211-point set of F_3^7 fixed in the
   layer x0 = 0: from the 427-point set's own layers and from fresh starts.  Lev's pattern
   gives 427; 211 + 2 * 105 = 421 is what "record + 2 k_(n-2)" would give.
3. A first layered set in F_3^9 on top of the 427-point layer (classical 2^10 - 1 = 1023,
   Lev's pattern 853).

Point lists in the report are delta-coded (encode/decode below): sorted integers, each gap
to the previous one (starting from -1) written in base 32, high digits from C32, the last
digit from T32.  Only the free layers are listed: for F_3^n with fixed layer L0 in F_3^m,
an internal point p stands for the point with index 3^m + p of F_3^n."""
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

SETS = os.path.join(_ROOT, "k6", "sets")
T32 = "abcdefghijklmnopqrstuvwxyz012345"
C32 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ6789+/"


def encode(pts):
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


def load(name):
    with open(os.path.join(SETS, name), encoding="utf-8") as f:
        return json.load(f)["points"]


def index(n, x):
    return sum(c * 3 ** (n - 1 - k) for k, c in enumerate(x))


def _plus(x, d):
    return tuple((a + b) % 3 for a, b in zip(x, d))


def line_stats(n, idx):
    """(full lines, single-line directions, essential points, uncovered directions)."""
    pts = [check.decode(n, i) for i in idx]
    S = set(pts)
    lines = single = unc = 0
    ess = set()
    for d in check.directions(n):
        inter, c = None, 0
        for x in pts:
            y = _plus(x, d)
            if y in S:
                z = _plus(y, d)
                if z in S:
                    c += 1
                    t = {x, y, z}
                    inter = t if inter is None else (inter & t)
        if c == 0:
            unc += 1
            continue
        lines += c // 3          # each line is found from each of its 3 points
        if c == 3:
            single += 1
        ess |= inter
    return lines, single, ess, unc


def removable_of(n, idx):
    lines, single, ess, unc = line_stats(n, idx)
    rem = sorted(i for i in idx if check.decode(n, i) not in ess)
    return lines, single, rem, unc


# ---- tasks -------------------------------------------------------------------------------
def l8_search(seed, noise, start_mode, deadline):
    t0 = time.time()
    W = layer.Layer(7, load("kakeya7_211.json"))
    start = None
    if start_mode == "427":
        start = [f - W.N for f in load("kakeya8_427.json") if f >= W.N]
    budget = max(5.0, deadline - time.time() - 12)
    r = layer.search(W, seed=seed, budget=budget, noise=noise, start=start)
    full = W.full_indices(r["points"])
    return {"kind": "L8", "seed": seed, "noise": noise, "start": start_mode or "-",
            "best": r["best"], "total": len(full), "reached": r["reached"][-3:],
            "steps": r["steps"], "verified": check.is_kakeya(8, full),
            "code": encode(r["points"]), "setup": round(time.time() - t0 - r["seconds"], 1)}


def l9_search(seed, noise, start_mode, deadline):
    t0 = time.time()
    W = layer.Layer(8, load("kakeya8_427.json"))
    setup = round(time.time() - t0, 1)
    start = None
    if start_mode == "cl":   # the classical layers {0,1}^8 and {0,2}^8 as a starting point
        A = [index(8, p) for p in itertools.product((0, 1), repeat=8)]
        B = [index(8, p) for p in itertools.product((0, 2), repeat=8)]
        start = A + [W.N + b for b in B]
    budget = max(5.0, deadline - time.time() - 50)
    r = layer.search(W, seed=seed, budget=budget, noise=noise, start=start)
    full = W.full_indices(r["points"])
    ok = check.is_kakeya(9, full) if len(full) < 1023 else None
    return {"kind": "L9", "seed": seed, "noise": noise, "start": start_mode or "-",
            "best": r["best"], "total": len(full), "reached": r["reached"][-3:],
            "steps": r["steps"], "verified": ok, "code": encode(r["points"]), "setup": setup}


def rem_then_search(deadline):
    out = {"kind": "REM"}
    idx = load("kakeya8_427.json")
    lines, single, rem, unc = removable_of(8, idx)
    out.update(lines=lines, single=single, removable=rem, unc=unc)
    cur, dropped, rng = list(idx), [], random.Random(7)
    while rem and not unc and time.time() < deadline - 120:
        p = rng.choice(rem)
        cur.remove(p)
        dropped.append(p)
        _, _, rem, unc = removable_of(8, cur)
    out["dropped"] = dropped
    out["final"] = len(cur)
    out["final_ok"] = check.is_kakeya(8, cur) if dropped else None
    out["final_code"] = encode(cur) if dropped else ""
    try:
        out["search"] = l8_search(15, 0.10, None, deadline)
    except Exception:  # noqa: BLE001
        out["search"] = {"error": traceback.format_exc()[-300:]}
    return out


def run_task(t):
    kind, args = t[0], t[1:]
    try:
        if kind == "REM":
            return rem_then_search(*args)
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
    return "%s s%d nz%.2f st%s: extra %d total %d last %s steps %d ok=%s" % (
        r["kind"], r["seed"], r["noise"], r["start"], r["best"], r["total"],
        r["reached"][-2:], r["steps"], r["verified"])


def probe(seconds=440, workers=7):
    t0 = time.time()
    deadline = t0 + seconds - 20
    out, tail = [], []
    try:
        # self-check of the code before anything is reported with it
        sample = [0, 1, 2, 40, 41, 1000, 1031, 1032, 20000]
        assert decode(encode(sample)) == sample, "codec"
        tasks = [("REM", deadline),
                 ("L8", 11, 0.05, "427", deadline),
                 ("L8", 12, 0.10, "427", deadline),
                 ("L8", 13, 0.10, None, deadline),
                 ("L8", 14, 0.07, None, deadline),
                 ("L9", 21, 0.10, None, deadline),
                 ("L9", 22, 0.10, "cl", deadline)][:workers]
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(workers) as pool:
            res = pool.map(run_task, tasks, chunksize=1)
        out.append("mode=pool%d cpu=%s" % (workers, os.cpu_count()))
        best8, best9, rem = None, None, None
        for r in res:
            if r.get("kind") == "REM" and "error" not in r:
                rem = r
                s = r["search"]
                if "error" in s:
                    out.append("L8 s15 ERROR " + s["error"])
                else:
                    res.append(s)   # handled (once) when the loop reaches it
                continue
            out.append(_short(r))
            if "error" in r or not r["verified"]:
                continue
            if r["kind"] == "L8" and (best8 is None or r["total"] < best8["total"]):
                best8 = r
            if r["kind"] == "L9" and (best9 is None or r["total"] < best9["total"]):
                best9 = r
        if best9 is not None:
            tail.append("N9 total=%d (427 + %d; classical 1023, pattern 853) seed %d CODE %s" % (
                best9["total"], best9["best"], best9["seed"], best9["code"]))
        else:
            tail.append("N9 no verified set below 1023")
        if best8 is not None and best8["total"] < 427:
            tail.append("NEW N8 total=%d (211 + %d) seed %d CODE %s" % (
                best8["total"], best8["best"], best8["seed"], best8["code"]))
        else:
            tail.append("N8 best layered total %s, nothing below 427" % (
                best8["total"] if best8 else None))
        if rem is not None:
            tail.append("REM427 lines=%d single-line dirs=%d uncovered=%d removable=%s "
                        "dropped=%s final=%d ok=%s %s" % (
                            rem["lines"], rem["single"], rem["unc"], rem["removable"][:20],
                            rem["dropped"], rem["final"], rem["final_ok"],
                            ("CODE " + rem["final_code"]) if rem["final_code"] else ""))
        else:
            tail.append("REM427 failed")
    except Exception:  # noqa: BLE001
        tail.append("ERROR " + traceback.format_exc()[-700:])
    tail.append("total %.0fs" % (time.time() - t0))
    return " || ".join(out + tail)
