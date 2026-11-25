# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""ksym: smallest Kakeya sets in F_3^n that are invariant under groups from the field F_{3^n}.

For every group G(r, f) = H_r F_f (see gf.py: H_r = multiplicative subgroup of order r,
F_f = Galois subgroup of order f) this program finds the least size of a G-invariant Kakeya
set in F_3^n, n = 2..8, as a 0/1 program over G-orbits solved by CP-SAT:

    minimise sum_o |o| x_o   subject to: for every representative direction d,
    some line of direction d lies in the union of the chosen orbits.

A G-invariant set is Kakeya iff it covers one direction of each G-orbit of directions.  Every
set found is re-checked point by point with check.is_kakeya (which uses none of this code).
When CP-SAT reports OPTIMAL on a small orbit problem, bb.improve re-proves that no cheaper
union of orbits exists by an independent exhaustive search.

Known: k_2 = 7, k_3 = 13 (Lev), k_4 = 27, 45 <= k_5 <= 53 (Chung, Du, Wesley 53),
k_6 <= 105 and k_7 <= 213 (our walks), classical 2^(n+1) - 1, Lev pattern 107 / 213 / 427.

Phase 2 seeds the incremental WalkSAT search of munimtechnologies/kakeya-conjecture
(k6/fastwalk.py, downloaded at run time) with the best symmetric sets for n = 5, 6, 7.
"""
import argparse
import json
import math
import os
import sys
import time
import traceback
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gf  # noqa: E402
import check  # noqa: E402
import bb  # noqa: E402

RECORD = {2: 7, 3: 13, 4: 27, 5: 53, 6: 105, 7: 213, 8: 427}
LOWER = {2: 7, 3: 13, 4: 27, 5: 45, 6: 57, 7: 103, 8: 186}
CAP = {2: 400, 3: 400, 4: 400, 5: 400, 6: 400, 7: 400, 8: 260}


def divisors(m):
    return [d for d in range(1, m + 1) if m % d == 0]


_FIELDS = {}


def field(n):
    if n not in _FIELDS:
        _FIELDS[n] = gf.Field(n)
    return _FIELDS[n]


def plan_n(n):
    F = field(n)
    out = []
    for r in divisors(F.q - 1):
        for f in divisors(n):
            gens = F.generators(r, f)
            orb, members = gf.orbits(F.q, gens)
            if len(members) <= CAP[n]:
                out.append((len(members), r, f))
    out.sort()
    return out


def build(n, r, f):
    F = field(n)
    gens = F.generators(r, f)
    orb, members = gf.orbits(F.q, gens)
    reps, dsizes, ndirs = gf.direction_orbits(F, gens)
    masks = gf.line_masks(F, orb, reps)
    w = [len(m) for m in members]
    hist = {}
    for x in w:
        hist[x] = hist.get(x, 0) + 1
    return {"F": F, "orb": orb, "members": members, "w": w, "reps": reps, "masks": masks,
            "hist": hist, "ndirs": ndirs, "dsizes": dsizes}


def solve_cpsat(w, masks, seconds, workers=4, seed=1):
    from ortools.sat.python import cp_model
    m = cp_model.CpModel()
    x = [m.NewBoolVar("x%d" % o) for o in range(len(w))]
    ys = {}
    for opts in masks:
        lits = []
        for L in opts:
            bits = gf.mask_bits(L)
            if len(bits) == 1:
                lits.append(x[bits[0]])
                continue
            y = ys.get(L)
            if y is None:
                y = m.NewBoolVar("")
                for b in bits:
                    m.AddImplication(y, x[b])
                ys[L] = y
            lits.append(y)
        m.AddBoolOr(lits)
    m.Minimize(cp_model.LinearExpr.WeightedSum(x, w))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = float(seconds)
    s.parameters.num_workers = workers
    s.parameters.random_seed = seed
    t0 = time.time()
    st = s.Solve(m)
    name = s.StatusName(st)
    res = {"status": name, "seconds": round(time.time() - t0, 1), "aux_vars": len(ys)}
    if name in ("OPTIMAL", "FEASIBLE"):
        res["value"] = int(round(s.ObjectiveValue()))
        res["bound"] = int(math.ceil(s.BestObjectiveBound() - 1e-6))
        res["chosen"] = [o for o in range(len(w)) if s.Value(x[o])]
    return res


def time_for(n, norb):
    if n == 8:
        return 60 if norb <= 60 else 300
    if n == 7:
        return 60 if norb <= 60 else 900
    if norb <= 30:
        return 20
    if norb <= 60:
        return 60
    if norb <= 150:
        return 180
    return 420


def run_instance(n, r, f, seconds, workers):
    t0 = time.time()
    B = build(n, r, f)
    rec = {"n": n, "r": r, "f": f, "group_order": r * f, "orbits": len(B["w"]),
           "orbit_sizes": {str(k): v for k, v in sorted(B["hist"].items())},
           "direction_orbits": len(B["reps"]), "directions": B["ndirs"],
           "masks": sum(len(o) for o in B["masks"]), "build_s": round(time.time() - t0, 1),
           "time_limit": seconds, "poly": B["F"].poly}
    cp = solve_cpsat(B["w"], B["masks"], seconds, workers=workers, seed=n * 1000 + r * 10 + f)
    rec["cpsat"] = {k: v for k, v in cp.items() if k != "chosen"}
    if "chosen" in cp:
        pts = sorted(p for o in cp["chosen"] for p in B["members"][o])
        rec["size"] = len(pts)
        rec["verified"] = (len(pts) == cp["value"] and
                           check.is_kakeya(n, check.to_tuples(n, pts)))
        rec["points"] = pts
        if cp["status"] == "OPTIMAL" and len(B["w"]) <= 70:
            chk = bb.improve(B["w"], B["masks"], cp["value"], seconds=60)
            rec["bb_check"] = {"cheaper_found": chk["found"] is not None,
                               "cheaper_value": chk["value"], "complete": chk["complete"],
                               "nodes": chk["nodes"], "seconds": chk["seconds"]}
    rec["total_s"] = round(time.time() - t0, 1)
    return rec


# ---- phase 2: WalkSAT from the best symmetric sets --------------------------------------
RAW = "https://raw.githubusercontent.com/munimtechnologies/kakeya-conjecture/main/"


def fetch_walk(root="/tmp/kc"):
    for sub in ("k5", "k6"):
        os.makedirs(os.path.join(root, sub), exist_ok=True)
        with open(os.path.join(root, sub, "__init__.py"), "w") as fh:
            fh.write("")
    for path in ("k5/walk.py", "k6/fastwalk.py"):
        with urllib.request.urlopen(RAW + path, timeout=60) as resp:
            data = resp.read()
        with open(os.path.join(root, path), "wb") as fh:
            fh.write(data)
    if root not in sys.path:
        sys.path.insert(0, root)
    from k6 import fastwalk
    return fastwalk


def to_fw_index(n, pts):
    """Our integers (digit i = coefficient of t^i) -> fastwalk indices (first coordinate most
    significant).  This is a coordinate permutation, so Kakeya sets go to Kakeya sets."""
    out = []
    for t in check.to_tuples(n, pts):
        out.append(sum(c * 3 ** (n - 1 - k) for k, c in enumerate(t)))
    return sorted(out)


def phase2(results, seconds, workers, log):
    fw = fetch_walk()
    best = {}
    for rec in results:
        if rec.get("verified") and rec["n"] in (5, 6, 7):
            best.setdefault(rec["n"], []).append(rec)
    tasks = []
    picks = {5: 1, 6: 2, 7: 1}
    seed = 700
    for n in (6, 7, 5):
        cands = sorted(best.get(n, []), key=lambda r: (r["size"], -r["orbits"]))
        used = set()
        for rec in cands:
            if len([t for t in tasks if t[0] == n]) >= picks[n]:
                break
            if rec["size"] in used:
                continue
            used.add(rec["size"])
            seed += 1
            label = "sym-n%d-r%d-f%d-size%d" % (n, rec["r"], rec["f"], rec["size"])
            tasks.append((n, label, to_fw_index(n, rec["points"]), (), seed, seconds, 0.15, 0))
    tasks = tasks[:workers]
    log("phase2 tasks: %s" % [t[1] for t in tasks])
    if not tasks:
        return {"mode": "none", "results": []}
    res, mode = fw.run_pool(tasks, len(tasks))
    out = []
    for r in res:
        if "error" in r:
            out.append({"label": r.get("label"), "error": r["error"]})
            continue
        out.append({"label": r["label"], "n": int(r["label"].split("-")[1][1:]),
                    "best": r["best"], "verified": r["verified"], "reached": r["reached"][-6:],
                    "steps": r["steps"], "points_fw_index": r["points"]})
    return {"mode": mode, "results": out}


def summarise(results):
    best = {}
    for rec in results:
        if rec.get("verified"):
            n = rec["n"]
            if n not in best or rec["size"] < best[n]["size"]:
                best[n] = {"size": rec["size"], "r": rec["r"], "f": rec["f"],
                           "status": rec["cpsat"]["status"], "record": RECORD[n]}
    return best


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/out/ksym.json")
    ap.add_argument("--minutes", type=float, default=290)
    ap.add_argument("--phase2-minutes", type=float, default=40)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--dims", default="2,3,4,5,6,7,8")
    a = ap.parse_args(argv)
    t0 = time.time()
    deadline = t0 + a.minutes * 60

    def log(msg):
        print("[%8.1fs] %s" % (time.time() - t0, msg), flush=True)

    out = {"contract": "ksym-v1", "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "encoding": "point x = sum c_i 3^i, c_i = coefficient of t^i in F_3[t]/(poly); "
                       "poly lists f_0..f_(n-1) of the monic primitive t^n + ... + f_0",
           "instances": [], "skipped": [], "phase2": None, "errors": []}

    def dump():
        out["best"] = summarise(out["instances"])
        out["seconds"] = round(time.time() - t0, 1)
        tmp = a.out + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(out, fh)
        os.replace(tmp, a.out)

    dims = [int(x) for x in a.dims.split(",")]

    def run_dim(n, reserve):
        try:
            todo = plan_n(n)
        except Exception:  # noqa: BLE001
            out["errors"].append({"n": n, "plan": traceback.format_exc()[-800:]})
            return
        log("n=%d: %d instances" % (n, len(todo)))
        for norb, r, f in todo:
            if n >= 2 and (r, f) == (1, 1) and n >= 6:
                continue
            sec = time_for(n, norb)
            if (n, r, f) == (5, 1, 1):
                sec = 600
            if time.time() + sec + reserve > deadline:
                out["skipped"].append([n, r, f, norb])
                continue
            try:
                rec = run_instance(n, r, f, sec, a.workers)
                out["instances"].append(rec)
                log("n=%d r=%d f=%d orbits=%d dir-orbits=%d: %s size=%s bound=%s verified=%s "
                    "bb=%s (%ss)" % (n, r, f, rec["orbits"], rec["direction_orbits"],
                                     rec["cpsat"]["status"], rec.get("size"),
                                     rec["cpsat"].get("bound"), rec.get("verified"),
                                     rec.get("bb_check"), rec["total_s"]))
                if rec.get("verified") and rec["size"] < RECORD[n]:
                    log("NEW n=%d size=%d below %d (r=%d f=%d)" % (n, rec["size"], RECORD[n], r, f))
            except Exception:  # noqa: BLE001
                out["errors"].append({"n": n, "r": r, "f": f, "trace": traceback.format_exc()[-800:]})
                log("ERROR n=%d r=%d f=%d" % (n, r, f))
            dump()

    p2 = a.phase2_minutes * 60
    for n in [d for d in dims if d <= 7]:
        run_dim(n, reserve=p2 + 120)
    if p2 > 0 and time.time() + p2 + 120 < deadline:
        try:
            out["phase2"] = phase2(out["instances"], p2, a.workers, log)
            for r in out["phase2"]["results"]:
                log("phase2 %s" % {k: v for k, v in r.items() if k != "points_fw_index"})
        except Exception:  # noqa: BLE001
            out["errors"].append({"phase2": traceback.format_exc()[-800:]})
            log("phase2 ERROR")
        dump()
    if 8 in dims:
        run_dim(8, reserve=60)
    dump()
    log("best per n: %s" % json.dumps(out["best"]))
    log("done")


if __name__ == "__main__":
    main()
