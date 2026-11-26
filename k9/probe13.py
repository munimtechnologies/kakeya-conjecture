# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Report-only probe for k9/split.py (splits, the doubling lemma and the lifting lemma).

1. Direct split searches on the record sets 105 (F_3^6), 211, 422 and 844 (several seeds,
   in parallel).  Each result is re-checked by split.verify_split (the definition of C(A, B)).
2. CHAIN: the best split of the 105-point set lifted (Lemma 2) through the layered records
   211 -> 422 -> 844 -> 1693, switching to a direct split when one has a smaller excess;
   every split is verified.
3. DOUBLE: Lemma 1 on each record, and the iteration K -> K+ from the 844-point set up to
   F_3^13 as far as time allows; each step's split is verified, which proves the next set
   Kakeya.  The set doubled from the 105-point set (F_3^7) is also checked by
   k8.check.is_kakeya, and the one doubled from the 844-point set (F_3^10) by a pair check.
Codes of the splits (states over sorted K, see k9/split.py) come last in the report."""
import os
import sys
import time
import traceback

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from k9 import split  # noqa: E402

LEVELS = [(6, "kakeya6_105.json"), (7, "kakeya7_211.json"), (8, "kakeya8_422.json"),
          (9, "kakeya9_844.json"), (10, "kakeya10_1693.json")]


def run_task(t):
    name, m, seed, until = t
    try:
        n, K = split.load_set(name)
        assert n == m
        r = split.find_split(m, K, seed=seed, max_steps=10 ** 9, deadline=until)
        A, B = split.split_from_states(K, r["states"])
        ok, e, miss = split.verify_split(m, K, A, B)
        r["ok"], r["miss"], r["name"], r["seed"] = ok, miss, name, seed
        r["e_check"] = e
        del r["states"]
        return r
    except Exception:  # noqa: BLE001
        return {"name": name, "seed": seed, "error": traceback.format_exc()[-300:]}


def pair_kakeya(n, idx):
    """Pair check from the definition (three distinct points are collinear iff they sum to 0)."""
    P = sorted(set(split.digits(i, n) for i in idx))
    S = set(P)
    cov = set()
    for i, x in enumerate(P):
        for y in P[i + 1:]:
            if tuple((-a - b) % 3 for a, b in zip(x, y)) in S:
                d = tuple((b - a) % 3 for a, b in zip(x, y))
                if [c for c in d if c][0] == 2:
                    d = tuple((2 * c) % 3 for c in d)
                cov.add(d)
    return len(cov) == (3 ** n - 1) // 2


def probe(seconds=560, workers=8):
    t0 = time.time()
    deadline = t0 + seconds - 15
    out, tail, codes = [], [], []
    try:
        st = [1, 2, 3, 3, 1]
        assert split.decode_split(split.encode_split(st), 5) == st, "codec"
        tasks = [("kakeya6_105.json", 6, 11, t0 + 50), ("kakeya6_105.json", 6, 12, t0 + 50),
                 ("kakeya6_105.json", 6, 13, t0 + 50), ("kakeya7_211.json", 7, 14, t0 + 130),
                 ("kakeya7_211.json", 7, 15, t0 + 130), ("kakeya8_422.json", 8, 16, t0 + 200),
                 ("kakeya9_844.json", 9, 17, t0 + 250), ("kakeya9_844.json", 9, 18, t0 + 250)]
        import multiprocessing as mp
        pool = mp.get_context("spawn").Pool(workers)
        try:
            res = pool.map(run_task, tasks[:workers], chunksize=1)
        finally:
            pool.close()
            pool.join()
        out.append("pool%d cpu=%s %.0fs" % (workers, os.cpu_count(), time.time() - t0))
        direct = {}
        for r in res:
            if "error" in r:
                out.append("D %s s%d ERROR %s" % (r["name"], r["seed"], r["error"]))
                continue
            out.append("D%d s%d e=%d ok=%s miss=%d st%d lines%d single%d setup%.0f t%.0f" % (
                r["size"], r["seed"], r["excess"], r["ok"], r["miss"], r["steps"], r["lines"],
                r["single"], r["setup_s"], r["search_s"]))
            if r["ok"] and r["e_check"] == r["excess"]:
                b = direct.get(r["name"])
                if b is None or r["excess"] < b["excess"]:
                    direct[r["name"]] = r
        # 2. chain
        chain, A, B, K, m = [], None, None, None, None
        base9 = None
        for lev, name in LEVELS:
            n, Kn = split.load_set(name)
            src = None
            if A is not None:
                try:
                    A, B = split.lift(m, Kn, A, B)
                    src = "L"
                except ValueError:
                    A = None
            d = direct.get(name)
            if d is not None:
                st = split.decode_split(d["code"], len(Kn))
                Ad, Bd = split.split_from_states(Kn, st)
                if A is None or d["excess"] < len(A & B):
                    A, B, src = Ad, Bd, "D"
                    codes.append("CODE%d e=%d %s" % (len(Kn), d["excess"], d["code"]))
            K, m = Kn, n
            if A is None:
                chain.append("%d:none" % len(Kn))
                continue
            ok, e, miss = split.verify_split(m, K, A, B)
            chain.append("%d:e%d%s ok=%s miss=%d dbl=%d" % (len(K), e, src, ok, miss,
                                                             2 * len(K) + e))
            if not ok:
                A = None
            if lev == 6 and ok:
                from k8 import check
                D7, _, _ = split.double(6, K, A, B)
                chain.append("DEF7 %d %s" % (len(D7), check.is_kakeya(7, D7)))
            if lev == 9 and ok:
                base9 = (K, set(A), set(B))
        tail.append("CHAIN " + " | ".join(chain))
        # 3. iteration from the 844-point set
        if base9 is not None:
            K, A, B = base9
            m = 9
            it = []
            D, A, B = split.double(m, K, A, B)
            m = 10
            it.append("n10=%d pairs=%s" % (len(D), pair_kakeya(10, D)))
            while m <= 12:
                need = {10: 20, 11: 60, 12: 200}[m]
                if deadline - time.time() < need:
                    it.append("stop m=%d (time)" % m)
                    break
                ok, e, miss = split.verify_split(m, D, A, B)
                if not ok:
                    it.append("m=%d split FAILED miss=%d" % (m, miss))
                    break
                D, A, B = split.double(m, D, A, B)
                m += 1
                it.append("n%d=%d proved(e%d)" % (m, len(D), e))
            tail.append("ITER " + " ".join(it) + " t%.0f" % (time.time() - t0))
        # codes of the best direct splits of 105 and 211 always
        for name in ("kakeya6_105.json", "kakeya7_211.json"):
            d = direct.get(name)
            if d is not None:
                c = "CODE%d e=%d %s" % (d["size"], d["excess"], d["code"])
                if c not in codes:
                    codes.insert(0, c)
    except Exception:  # noqa: BLE001
        tail.append("ERROR " + traceback.format_exc()[-700:])
    tail.append("total %.0fs" % (time.time() - t0))
    return " || ".join(out + tail + codes)
