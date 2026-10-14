"""Exact minimum Kakeya set sizes in AG(n, q) by integer programming, with two independent solvers.

Model: x_p = 1 if point p is in the set, y_L = 1 if line L is the chosen line of its direction.
Exactly one chosen line per direction, every point of a chosen line is in the set, minimise the
number of points. Solved with OR-Tools CP-SAT and with HiGHS (scipy.optimize.milp); every
returned set is re-checked with kakeya.is_kakeya. Optional symmetry fixing uses
kakeya.symmetry_lines; the small cases are also solved without it as a check.

    python solve.py --out-dir /out [--quick]

By MV2 of Munim, Inc. Made with AI agents working in MV2's name. Apache-2.0.
"""

import argparse
import json
import math
import os
import time

from kakeya import AffineSpace, bukh_chao_bound, dvir_bound, is_kakeya, plane_minimum, symmetry_lines

# (n, q, CP-SAT seconds, HiGHS seconds, also solve without symmetry fixing)
CASES = [
    (2, 2, 60, 60, True), (2, 3, 60, 60, True), (2, 4, 60, 60, True), (2, 5, 60, 60, True),
    (2, 7, 60, 60, True), (2, 8, 60, 60, True), (2, 9, 60, 60, True),
    (3, 2, 120, 120, True), (4, 2, 120, 120, True), (5, 2, 300, 300, False),
    (3, 3, 600, 600, True),
    (3, 4, 1200, 600, False),
    (4, 3, 1500, 600, False),
    (3, 5, 1500, 600, False),
    (6, 2, 300, 300, False),
    (4, 4, 900, 600, False),
    (3, 7, 1200, 600, False),
    (5, 3, 1200, 600, False),
]
QUICK = [(2, 3, 20, 20, True), (3, 2, 20, 20, True), (3, 3, 60, 60, True)]


def solve_cpsat(space, seconds, symmetry, workers):
    from ortools.sat.python import cp_model
    m = cp_model.CpModel()
    x = [m.NewBoolVar("x%d" % i) for i in range(space.size)]
    allowed = symmetry_lines(space) if symmetry else {}
    for di, lines in enumerate(space.lines):
        row = []
        for L in lines:
            y = m.NewBoolVar("")
            if di in allowed and L not in allowed[di]:
                m.Add(y == 0)
            for p in L:
                m.AddImplication(y, x[p])
            row.append(y)
        m.AddExactlyOne(row)
    m.Minimize(sum(x))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = float(seconds)
    s.parameters.num_workers = workers
    t0 = time.time()
    st = s.Solve(m)
    out = {"solver": "cp-sat", "symmetry_fixing": symmetry, "status": s.StatusName(st),
           "seconds": round(time.time() - t0, 2), "time_limit": seconds}
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        pts = [i for i in range(space.size) if s.Value(x[i])]
        out["value"] = len(pts)
        out["lower_bound"] = int(math.ceil(s.BestObjectiveBound() - 1e-6))
        out["proved_optimal"] = st == cp_model.OPTIMAL
        out["verified_kakeya"] = is_kakeya(pts, space)
        out["points"] = pts
    return out


def solve_highs(space, seconds, symmetry):
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_matrix
    N = space.size
    flat = [(di, L) for di, lines in enumerate(space.lines) for L in lines]
    M = len(flat)
    allowed = symmetry_lines(space) if symmetry else {}
    rows, cols, vals, lb, ub = [], [], [], [], []
    r = 0
    k = 0
    for lines in space.lines:
        for _ in lines:
            rows.append(r)
            cols.append(N + k)
            vals.append(1.0)
            k += 1
        lb.append(1.0)
        ub.append(1.0)
        r += 1
    for k, (di, L) in enumerate(flat):
        for p in L:
            rows += [r, r]
            cols += [p, N + k]
            vals += [1.0, -1.0]
            lb.append(0.0)
            ub.append(np.inf)
            r += 1
    A = coo_matrix((vals, (rows, cols)), shape=(r, N + M)).tocsr()
    upper = np.ones(N + M)
    for k, (di, L) in enumerate(flat):
        if di in allowed and L not in allowed[di]:
            upper[N + k] = 0.0
    c = np.concatenate([np.ones(N), np.zeros(M)])
    t0 = time.time()
    res = milp(c, constraints=LinearConstraint(A, lb, ub), integrality=np.ones(N + M),
               bounds=Bounds(np.zeros(N + M), upper),
               options={"time_limit": float(seconds), "disp": False})
    out = {"solver": "highs", "symmetry_fixing": symmetry, "status": int(res.status),
           "message": str(res.message), "seconds": round(time.time() - t0, 2),
           "time_limit": seconds}
    if res.x is not None:
        pts = [i for i in range(N) if res.x[i] > 0.5]
        out["value"] = len(pts)
        bound = getattr(res, "mip_dual_bound", None)
        if bound is not None and np.isfinite(bound):
            out["lower_bound"] = int(math.ceil(float(bound) - 1e-6))
        out["proved_optimal"] = res.status == 0
        out["verified_kakeya"] = is_kakeya(pts, space)
        out["points"] = pts
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    results = []
    for n, q, t_cp, t_hi, nosym in (QUICK if a.quick else CASES):
        space = AffineSpace(n, q)
        case = {"n": n, "q": q, "points": space.size, "directions": len(space.directions),
                "lines": sum(len(x) for x in space.lines),
                "dvir_bound": dvir_bound(n, q), "bukh_chao_bound": bukh_chao_bound(n, q),
                "plane_known_minimum": plane_minimum(q) if n == 2 else None, "runs": []}
        case["runs"].append(solve_cpsat(space, t_cp, True, a.workers))
        if nosym:
            case["runs"].append(solve_cpsat(space, t_cp, False, a.workers))
        case["runs"].append(solve_highs(space, t_hi, True))
        proved = [r["value"] for r in case["runs"] if r.get("proved_optimal") and r.get("verified_kakeya")]
        found = [r["value"] for r in case["runs"] if r.get("verified_kakeya")]
        lows = [r["lower_bound"] for r in case["runs"] if "lower_bound" in r]
        case["best_found"] = min(found) if found else None
        case["best_lower_bound"] = max(lows + [case["bukh_chao_bound"]])
        case["proved_minimum"] = proved[0] if proved and len(set(proved)) == 1 else None
        case["solvers_agree"] = len(set(proved)) <= 1
        best_run = min((r for r in case["runs"] if r.get("verified_kakeya")),
                       key=lambda r: r["value"], default=None)
        if best_run:
            case["certificate"] = [list(space.coords[i]) for i in best_run["points"]]
        results.append(case)
        line = "AG(%d,%d): best %s, lower %s, proved %s, agree %s, runs %s" % (
            n, q, case["best_found"], case["best_lower_bound"], case["proved_minimum"],
            case["solvers_agree"],
            [(r["solver"], r["symmetry_fixing"], r["status"], r.get("value"), r["seconds"])
             for r in case["runs"]])
        print(line, flush=True)
        with open(os.path.join(a.out_dir, "results.json"), "w") as f:
            json.dump(results, f, indent=1)
        with open(os.path.join(a.out_dir, "results.txt"), "a") as f:
            f.write(line + "\n")


if __name__ == "__main__":
    main()
