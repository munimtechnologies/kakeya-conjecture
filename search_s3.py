"""Targeted search for small rainbow triples in F_3^3 (see TRIPLES.md).

The energy lemma proves s(3) >= 7, and G1^3 gives s(3) <= 8. This job asks:
  - (7, 7, 7): is there a rainbow triple of F_3^3 with three 7-point sets? It would give explicit
    Kakeya sets with base 7^(1/3) = 1.913, better than the two-symbol construction's 1.9626.
  - (6, 7, 7): product 294 < 6.75^3, so such a triple would give base 294^(1/9) ~ 1.880, below the
    ~1.890 of Kopparty-Lev-Saraf-Sudan's random construction.
  - every other size pattern with product below 512 = 8^3 that the energy lemma leaves, in
    increasing order of product, until one is feasible (that product bounds tau(3) from above;
    every pattern proved INFEASIBLE before it bounds tau(3) from below).
Each attempt runs OR-Tools CP-SAT. A symmetric attempt also asks the three sets to be invariant
under the cyclic shift of coordinates ("cyclic") or under x -> -x ("negation"); those searches are
smaller, so they can only find, never rule out. Every triple found is re-checked in plain Python.
Results: out/s3.json and out/summary.txt.

    python search_s3.py --out-dir /out [--workers 2] [--scale 1.0] [--budget 1200]

By MV2 of Munim, Inc. Made with AI agents working in MV2's name. Apache-2.0.
"""

import argparse
import json
import os
import time
import traceback

from triples import add, is_rainbow, scale, sizes, surviving_patterns, vectors

N_DIM = 3


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
        with open(os.path.join(self.dir, "s3.json"), "w") as f:
            json.dump(self.results, f, indent=1)


def symmetry_map(pts, mode):
    idx = {p: i for i, p in enumerate(pts)}
    if mode == "cyclic":
        return [idx[p[1:] + p[:1]] for p in pts]
    if mode == "negation":
        return [idx[scale(2, p)] for p in pts]
    raise ValueError(mode)


def build(n, caps, symmetric):
    from ortools.sat.python import cp_model
    pts = vectors(n)
    idx = {p: i for i, p in enumerate(pts)}
    m = cp_model.CpModel()
    X = [[m.NewBoolVar("x%d_%d" % (k, i)) for i in range(len(pts))] for k in range(3)]
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
    for k in range(3):
        m.Add(sum(X[k]) <= caps[k])
    if symmetric:
        # Translating to the origin could break the symmetry, so the origin is not fixed here.
        sh = symmetry_map(pts, symmetric)
        for k in range(3):
            for i in range(len(pts)):
                m.Add(X[k][i] == X[k][sh[i]])
    else:
        # The d = 0 witness lies in all three sets; translating it to the origin is no loss.
        z = idx[(0,) * n]
        for k in range(3):
            m.Add(X[k][z] == 1)
    return m, X, pts


def attempt(caps, symmetric, seconds, workers):
    from ortools.sat.python import cp_model
    m, X, pts = build(N_DIM, caps, symmetric)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = float(seconds)
    s.parameters.num_workers = workers
    t0 = time.time()
    st = s.Solve(m)
    found = st in (cp_model.OPTIMAL, cp_model.FEASIBLE)
    a, b, c = caps
    P = a * b * c
    rec = {"caps": list(caps), "product": P, "symmetric": symmetric or False,
           "status": s.StatusName(st), "seconds": round(time.time() - t0, 1), "found": found,
           "base_if_feasible": round(P ** (1.0 / (3 * N_DIM)), 5)}
    if found:
        T = tuple({pts[i] for i in range(len(pts)) if s.Value(X[k][i])} for k in range(3))
        rec["sizes"] = list(sizes(T))
        rec["verified_rainbow"] = is_rainbow(*T, N_DIM)
        rec["triple"] = [sorted(list(p) for p in S) for S in T]
    return rec


def describe(rec):
    return ("caps %s (abc %d, base %.4f)%s: %s, %.0fs%s"
            % (tuple(rec["caps"]), rec["product"], rec["base_if_feasible"],
               " " + rec["symmetric"] if rec["symmetric"] else "", rec["status"], rec["seconds"],
               " sizes %s verified %s" % (rec["sizes"], rec["verified_rainbow"])
               if rec["found"] else ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--scale", type=float, default=1.0, help="multiply every time limit")
    ap.add_argument("--budget", type=float, default=1200.0, help="seconds for the pattern sweep")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    log = Log(a.out_dir)

    pats = surviving_patterns(N_DIM, 512)
    log.put("surviving", {"count": len(pats), "patterns": [list(p) for p in pats]})
    log.line("energy lemma: %d size patterns with abc < 512 survive in F_3^3; smallest %s"
             % (len(pats), pats[:5]))

    tried = set()
    found_caps = set()
    plan = (((7, 7, 7), "cyclic", 300), ((7, 7, 7), "negation", 300),
            ((6, 7, 7), "cyclic", 150), ((7, 7, 7), None, 1500), ((6, 7, 7), None, 600))
    for caps, symmetric, sec in plan:
        if caps in found_caps:
            continue
        try:
            rec = attempt(caps, symmetric, sec * a.scale, a.workers)
        except Exception:
            log.line("attempt %s failed:\n%s" % (caps, traceback.format_exc()))
            continue
        log.put("attempts", rec)
        log.line(describe(rec))
        if rec["found"]:
            found_caps.add(caps)
        if not symmetric or rec["found"]:
            tried.add(caps)

    spent = 0.0
    for p in pats:
        if p in tried:
            continue
        if spent > a.budget * a.scale:
            log.line("pattern sweep: budget used, stopped before %s" % (p,))
            break
        try:
            rec = attempt(p, None, 60 * a.scale, a.workers)
        except Exception:
            log.line("pattern %s failed:\n%s" % (p, traceback.format_exc()))
            break
        spent += rec["seconds"]
        log.put("patterns", rec)
        log.line(describe(rec))
        if rec["found"]:
            break
    log.line("done")


if __name__ == "__main__":
    main()
