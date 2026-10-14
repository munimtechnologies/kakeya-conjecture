"""Search for small rainbow triples in F_3^n and for Kakeya sets built layer by layer.

Parts (results in out/triples.json, log in out/summary.txt; every set found is re-checked in plain
Python before it is reported):
  sigma    least |A| + |M| + |C| of a rainbow triple, n = 1..4 (sigma(n) <= K(n+1, 3)).
  maxsize  least max(|A|, |M|, |C|), n = 1..4: a triple with max s gives K(N, 3) = O(s^(N/n)).
  sizes    for n = 2, 3, size patterns a <= b <= c in increasing order of abc, while abc < 8^n
           (base below 2): first feasible pattern, every earlier one proved infeasible or not.
  doubling least |A| + |C| with M inside A, n = 1..4, plus a solver-free exhaustive check for
           n = 2 (the doubling lemma's constant |A| + |C| over 2^n; Smyth's chain has 2).
  chain    Kakeya sets of F_3^(n+1) built on a fixed Kakeya set B of F_3^n: least
           |A| + |C| + |M minus B|, then B' = A x {0} + (B u M) x {1} + C x {2}; from the 13-point
           set of F_3^3 up to F_3^6 and from the 27-point set of F_3^4 up to F_3^6.

    python search_triples.py --out-dir /out [--workers 4] [--scale 1.0]

By MV2 of Munim, Inc. Made with AI agents working in MV2's name. Apache-2.0.
"""

import argparse
import json
import math
import os
import time
import traceback
from itertools import combinations

from triples import (G1, K33, add, is_kakeya, is_rainbow, lift, product_triple, scale, sizes,
                     vectors)

# The 27-point Kakeya set of F_3^4 found by kakeya-min-v1 (points encoded as sum c_i 3^i).
V1_43 = [0, 1, 2, 8, 9, 10, 12, 13, 15, 17, 18, 21, 25, 34, 37, 39, 43, 45, 50, 56, 57, 60, 64,
         65, 70, 75, 77]

# Best upper bounds on K(n, 3) before this job (kakeya-conjecture README section 7, Smyth).
PREVIOUS = {3: 13, 4: 27, 5: 56, 6: 127, 7: 255}


def decode(p, n):
    return tuple((p // 3 ** j) % 3 for j in range(n))


class Log:
    def __init__(self, out_dir):
        self.dir = out_dir
        self.results = {}

    def line(self, text):
        print(text, flush=True)
        with open(os.path.join(self.dir, "summary.txt"), "a") as f:
            f.write(text + "\n")

    def put(self, key, value):
        self.results.setdefault(key, []).append(value)
        with open(os.path.join(self.dir, "triples.json"), "w") as f:
            json.dump(self.results, f, indent=1)


def model(n, origin=True, middle_in_first=False, caps=None, base=None):
    """CP-SAT model of a rainbow triple (A, M, C) in F_3^n; returns (model, vars, points)."""
    from ortools.sat.python import cp_model
    pts = vectors(n)
    idx = {p: i for i, p in enumerate(pts)}
    m = cp_model.CpModel()
    X = [[m.NewBoolVar("") for _ in pts] for _ in range(3)]
    for d in pts:
        d2 = scale(2, d)
        ys = []
        for a in pts:
            y = m.NewBoolVar("")
            m.AddImplication(y, X[0][idx[a]])
            m.AddImplication(y, X[1][idx[add(a, d)]])
            m.AddImplication(y, X[2][idx[add(a, d2)]])
            ys.append(y)
        m.AddBoolOr(ys)
    if origin:
        # The d = 0 witness lies in all three sets; translating it to the origin is no loss.
        z = idx[(0,) * n]
        for k in range(3):
            m.Add(X[k][z] == 1)
    if middle_in_first:
        for i in range(len(pts)):
            m.AddImplication(X[1][i], X[0][i])
    if caps:
        for k in range(3):
            m.Add(sum(X[k]) <= caps[k])
    return m, X, pts, idx


def run_solver(m, seconds, workers):
    from ortools.sat.python import cp_model
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = float(seconds)
    s.parameters.num_workers = workers
    t0 = time.time()
    st = s.Solve(m)
    found = st in (cp_model.OPTIMAL, cp_model.FEASIBLE)
    return s, st, s.StatusName(st), round(time.time() - t0, 1), found


def extract(s, X, pts):
    return tuple({pts[i] for i in range(len(pts)) if s.Value(X[k][i])} for k in range(3))


def record(T, n):
    a, b, c = sizes(T)
    P = a * b * c
    return {"sizes": [a, b, c], "sum": a + b + c, "product": P,
            "base_from_product": round(P ** (1.0 / (3 * n)), 5),
            "base_from_max": round(max(a, b, c) ** (1.0 / n), 5),
            "verified_rainbow": is_rainbow(*T, n),
            "triple": [sorted(list(p) for p in S) for S in T]}


def part_sigma(log, workers, sc):
    for n, sec in ((1, 10), (2, 60), (3, 600), (4, 900)):
        m, X, pts, _ = model(n)
        for k in range(2):
            m.Add(sum(X[k]) <= sum(X[k + 1]))
        m.Minimize(sum(X[0]) + sum(X[1]) + sum(X[2]))
        s, st, name, secs, found = run_solver(m, sec * sc, workers)
        out = {"n": n, "status": name, "seconds": secs}
        if found:
            out.update(record(extract(s, X, pts), n))
            out["lower_bound"] = int(math.ceil(s.BestObjectiveBound() - 1e-6))
        log.put("sigma", out)
        log.line("sigma n=%d: %s sum %s (lower %s) sizes %s, %.0fs"
                 % (n, name, out.get("sum"), out.get("lower_bound"), out.get("sizes"), secs))


def part_maxsize(log, workers, sc):
    for n, sec in ((1, 10), (2, 60), (3, 600), (4, 900)):
        m, X, pts, _ = model(n)
        S = m.NewIntVar(1, 3 ** n, "S")
        for k in range(3):
            m.Add(sum(X[k]) <= S)
        m.Minimize(S)
        s, st, name, secs, found = run_solver(m, sec * sc, workers)
        out = {"n": n, "status": name, "seconds": secs}
        if found:
            out.update(record(extract(s, X, pts), n))
            out["max"] = max(out["sizes"])
            out["lower_bound"] = int(math.ceil(s.BestObjectiveBound() - 1e-6))
        log.put("maxsize", out)
        log.line("maxsize n=%d: %s max %s (lower %s) sizes %s base %s, %.0fs"
                 % (n, name, out.get("max"), out.get("lower_bound"), out.get("sizes"),
                    out.get("base_from_max"), secs))


def size_patterns(n):
    """Patterns a <= b <= c with pairwise products >= 3^n and abc < 8^n, by increasing abc."""
    N = 3 ** n
    pats = []
    for a in range(1, N + 1):
        for b in range(a, N + 1):
            if a * b < N:
                continue
            for c in range(b, N + 1):
                if a * b * c >= 8 ** n:
                    break
                pats.append((a * b * c, a + b + c, a, b, c))
    return sorted(pats)


def part_sizes(log, workers, sc):
    for n, sec, budget in ((2, 30, 1200), (3, 60, 3000)):
        spent = 0.0
        for P, _, a, b, c in size_patterns(n):
            if spent > budget * sc:
                log.line("sizes n=%d: budget used, stopped before %s" % (n, (a, b, c)))
                break
            m, X, pts, _ = model(n, caps=(a, b, c))
            s, st, name, secs, found = run_solver(m, sec * sc, workers)
            spent += secs
            out = {"n": n, "pattern": [a, b, c], "product": P, "status": name, "seconds": secs,
                   "base_if_feasible": round(P ** (1.0 / (3 * n)), 5)}
            if found:
                out.update(record(extract(s, X, pts), n))
            log.put("sizes", out)
            log.line("sizes n=%d %s (abc %d, base %.4f): %s, %.0fs"
                     % (n, (a, b, c), P, out["base_if_feasible"], name, secs))
            if found:
                break


def exhaustive_doubling_2():
    """Least |A| + |C| over rainbow (A, A, C) in F_3^2 by plain enumeration (sum <= 8)."""
    pts = vectors(2)
    best = None
    for total in range(2, 9):
        for a in range(1, total):
            for A in combinations(pts, a):
                SA = set(A)
                for C in combinations(pts, total - a):
                    if is_rainbow(SA, SA, set(C), 2):
                        return total, (sorted(SA), sorted(C))
    return best


def part_doubling(log, workers, sc):
    t0 = time.time()
    try:
        total, ex = exhaustive_doubling_2()
        log.put("doubling_exhaustive", {"n": 2, "least_sum": total, "example": ex,
                                        "seconds": round(time.time() - t0, 1)})
        log.line("doubling n=2 exhaustive: least |A|+|C| = %d (%.0fs)" % (total, time.time() - t0))
    except Exception:
        log.line("doubling exhaustive failed:\n%s" % traceback.format_exc())
    for n, sec in ((1, 10), (2, 60), (3, 600), (4, 900)):
        m, X, pts, _ = model(n, middle_in_first=True)
        m.Minimize(sum(X[0]) + sum(X[2]))
        s, st, name, secs, found = run_solver(m, sec * sc, workers)
        out = {"n": n, "status": name, "seconds": secs}
        if found:
            T = extract(s, X, pts)
            out.update(record(T, n))
            out["A_plus_C"] = len(T[0]) + len(T[2])
            out["constant"] = round(out["A_plus_C"] / 2 ** n, 5)
            out["middle_in_first"] = T[1] <= T[0]
            out["lower_bound"] = int(math.ceil(s.BestObjectiveBound() - 1e-6))
        log.put("doubling", out)
        log.line("doubling n=%d: %s |A|+|C| %s (lower %s), constant %s (Smyth 2), %.0fs"
                 % (n, name, out.get("A_plus_C"), out.get("lower_bound"), out.get("constant"),
                    secs))


def chain_step(B, n, seconds, workers, hint=None):
    """Least |A| + |C| + |M - B| over rainbow (A, M, C) in F_3^n; B is a Kakeya set of F_3^n."""
    m, X, pts, idx = model(n, origin=False)
    inB = [p in B for p in pts]
    m.Minimize(sum(X[0]) + sum(X[2]) + sum(X[1][i] for i in range(len(pts)) if not inB[i]))
    if hint:
        for k in range(3):
            for i, p in enumerate(pts):
                m.AddHint(X[k][i], 1 if p in hint[k] else 0)
    s, st, name, secs, found = run_solver(m, seconds, workers)
    if not found:
        return {"status": name, "seconds": secs}, None, None
    T = extract(s, X, pts)
    B2 = lift(T[0], set(B) | T[1], T[2])
    out = {"status": name, "seconds": secs, "triple_sizes": list(sizes(T)),
           "lower_bound_on_cost": int(math.ceil(s.BestObjectiveBound() - 1e-6)),
           "cost": len(B2) - len(B), "rainbow": is_rainbow(*T, n)}
    return out, T, B2


def part_chain(log, workers, sc):
    starts = [("K33", set(K33), 3, (300, 900, 1800)),
              ("v1(4,3)", {decode(p, 4) for p in V1_43}, 4, (900, 1800))]
    for label, B, n, limits in starts:
        if not is_kakeya(B, n):
            log.line("chain %s: start is not Kakeya, skipped" % label)
            continue
        T_prev = None
        for sec in limits:
            hint = product_triple(T_prev, G1) if T_prev else None
            try:
                out, T, B2 = chain_step(B, n, sec * sc, workers, hint)
            except Exception:
                log.line("chain %s n=%d failed:\n%s" % (label, n, traceback.format_exc()))
                break
            out.update({"start": label, "from_dim": n, "to_dim": n + 1, "from_size": len(B)})
            if B2 is None:
                log.put("chain", out)
                log.line("chain %s %d->%d: %s, no set" % (label, n, n + 1, out["status"]))
                break
            out["size"] = len(B2)
            out["verified_kakeya"] = is_kakeya(B2, n + 1)
            out["previous_best"] = PREVIOUS.get(n + 1)
            out["improves"] = bool(out["verified_kakeya"] and PREVIOUS.get(n + 1)
                                   and len(B2) < PREVIOUS[n + 1])
            out["set"] = sorted(list(p) for p in B2)
            log.put("chain", out)
            log.line("chain %s %d->%d: |B'| = %d (cost %d, lower %d, %s), Kakeya %s, previous %s%s"
                     % (label, n, n + 1, len(B2), out["cost"], out["lower_bound_on_cost"],
                        out["status"], out["verified_kakeya"], out["previous_best"],
                        " IMPROVES" if out["improves"] else ""))
            if not out["verified_kakeya"]:
                break
            B, n, T_prev = B2, n + 1, T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--scale", type=float, default=1.0, help="multiply every time limit")
    ap.add_argument("--parts", default="maxsize,sizes,chain,doubling,sigma")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    log = Log(a.out_dir)
    parts = {"sigma": part_sigma, "maxsize": part_maxsize, "doubling": part_doubling,
             "sizes": part_sizes, "chain": part_chain}
    for name in a.parts.split(","):
        try:
            parts[name](log, a.workers, a.scale)
        except Exception:
            log.line("part %s failed:\n%s" % (name, traceback.format_exc()))
    log.line("done")


if __name__ == "__main__":
    main()
