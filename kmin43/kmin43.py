# kmin43.py -- is there a Kakeya set of 26 points in F_3^4?  Known: 18 <= K(4,3) <= 27.
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0.
#
# Reduction (kmin43/prereg-kmin43.md in munimtechnologies/kakeya-conjecture).
# Write F_3^4 = F_3^3 x F_3 with points (y, s) and slices K_s = {y : (y, s) in K}.
#  * A direction (w, 1), w in F_3^3 (27 "transversal" directions), has a line in K iff
#    some a in K_0 has a + w in K_1 and a + 2w in K_2.
#  * A direction (w, 0), w in PG(2, 3) (13 "horizontal" directions), has a line in K iff some
#    slice K_s contains a line of F_3^3 in direction w.
#  * Let K be a Kakeya set with |K| <= 26 and let H be an affine hyperplane with |K n H| = n
#    minimal (n <= 8 because the three parallel sections sum to |K|; n >= 1 because a
#    transversal line meets every parallel hyperplane). Move H to {s = 0} by an affine map
#    of F_3^4, then move K_0 to its canonical representative R under AGL(3, 3) acting on y
#    alone ((y, s) -> (A y + b, s) is affine on F_3^4). The image is again a Kakeya set, so
#    (K_1, K_2) solves MODEL(R, 26): slice 0 equal to R, transversal and horizontal
#    conditions, every one of the 120 hyperplane sections >= |R|, and |K| <= 26.
#  * So if MODEL(R, 26) is infeasible for every orbit representative R with 1 <= |R| <= 8,
#    then K(4, 3) >= 27, and with the known 27-point set K(4, 3) = 27.
# Each model has 54 point variables. It is solved by OR-Tools CP-SAT and, independently, by
# HiGHS (scipy.optimize.milp). Every solution found is re-checked by a brute-force Kakeya test.
#
# Self-checks: AGL(3,3) has 303,264 distinct permutations; orbit sizes of the representatives
# add up to C(27, k) for every k; the known 27-point set (kakeya-min-v1) is Kakeya, and after
# normalisation both solvers accept it at bound 27 and reject the same fixed set at bound 26;
# the unfixed model for its representative is feasible at bound 27.

import argparse
import itertools
import json
import math
import time

import numpy as np

PTS3 = list(itertools.product(range(3), repeat=3))
IDX3 = {p: i for i, p in enumerate(PTS3)}
GROUP_ORDER = 303264  # |AGL(3, 3)| = 27 * 11232


def add(p, q):
    return tuple((a + b) % 3 for a, b in zip(p, q))


def smul(c, p):
    return tuple((c * a) % 3 for a in p)


def dot(u, v):
    return sum(a * b for a, b in zip(u, v)) % 3


def projective_points(n):
    """Directions of F_3^n: non-zero vectors whose first non-zero coordinate is 1."""
    out = []
    for v in itertools.product(range(3), repeat=n):
        nz = [x for x in v if x]
        if nz and nz[0] == 1:
            out.append(v)
    return out


DIR3 = projective_points(3)  # 13
DIR4 = projective_points(4)  # 40


def lines3(w):
    """The 9 lines of F_3^3 in direction w, each as a sorted tuple of 3 point indices."""
    seen = set()
    for a in PTS3:
        seen.add(tuple(sorted(IDX3[add(a, smul(t, w))] for t in range(3))))
    return sorted(seen)


LINES3 = {w: lines3(w) for w in DIR3}


def is_kakeya(K, n):
    """Brute force: every direction of F_3^n has a full line inside K."""
    K = set(tuple(x) for x in K)
    for v in projective_points(n):
        ok = False
        for x in K:
            if add(x, v) in K and add(x, smul(2, v)) in K:
                ok = True
                break
        if not ok:
            return False
    return True


def slice_kakeya(K):
    """The transversal + horizontal test on the slices by the LAST coordinate (for tests)."""
    S = {s: set(tuple(x[:3]) for x in K if x[3] == s) for s in range(3)}
    for w in PTS3:
        if not any(add(a, w) in S[1] and add(a, smul(2, w)) in S[2] for a in S[0]):
            return False
    for w in DIR3:
        ok = False
        for s in range(3):
            for L in LINES3[w]:
                if all(PTS3[i] in S[s] for i in L):
                    ok = True
        if not ok:
            return False
    return True


def sections(K):
    """All 120 affine hyperplane sections of K in F_3^4: (size, normal, c)."""
    out = []
    for nrm in DIR4:
        for c in range(3):
            out.append((sum(1 for x in K if dot(nrm, x) == c), nrm, c))
    return out


# ---------------------------------------------------------------- the group AGL(3, 3)

def det3(M):
    (a, b, c), (d, e, f), (g, h, i) = M
    return (a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)) % 3


def agl3_perms():
    """Array (303264, 27): row g maps point index i to the index of g(PTS3[i])."""
    lin = []
    for ent in itertools.product(range(3), repeat=9):
        M = (ent[0:3], ent[3:6], ent[6:9])
        if det3(M):
            lin.append([IDX3[tuple(dot(M[r], p) for r in range(3))] for p in PTS3])
    lin = np.array(lin, dtype=np.int16)
    T = np.array([[IDX3[add(PTS3[i], b)] for i in range(27)] for b in PTS3], dtype=np.int16)
    perms = T[:, lin]  # [b, m, i] = T[b, lin[m, i]]: first the linear map, then + b
    return perms.reshape(-1, 27)


def codes_of(S, perms):
    S = sorted(S)
    k = len(S)
    img = np.sort(perms[:, S].astype(np.int64), axis=1)
    w = np.array([27 ** (k - 1 - i) for i in range(k)], dtype=np.int64)
    return img @ w


def encode(S):
    S = sorted(S)
    k = len(S)
    return sum(int(x) * 27 ** (k - 1 - i) for i, x in enumerate(S))


def decode(code, k):
    out = []
    for _ in range(k):
        out.append(code % 27)
        code //= 27
    return out[::-1]


def canon(S, perms):
    """(canonical code, index of a group element g with sorted g(S) = decode(code))."""
    c = codes_of(S, perms)
    g = int(np.argmin(c))
    return int(c[g]), g


def stabiliser(S, perms):
    return int((codes_of(S, perms) == encode(S)).sum())


def orbit_reps(kmax, perms, log=None):
    """Canonical codes of the AGL(3,3)-orbits of k-subsets of F_3^3, k = 1..kmax."""
    reps = {1: [canon([0], perms)[0]]}
    for k in range(2, kmax + 1):
        seen = set()
        for code in reps[k - 1]:
            S = decode(code, k - 1)
            for p in range(27):
                if p not in S:
                    seen.add(canon(S + [p], perms)[0])
        reps[k] = sorted(seen)
        if log:
            log(f"k={k}: {len(reps[k])} orbits")
    return reps


def orbit_check(reps, perms):
    """For each k: (sum of orbit sizes, C(27, k)) and the stabiliser orders."""
    out = {}
    for k, codes in reps.items():
        stabs = [stabiliser(decode(c, k), perms) for c in codes]
        total = sum(GROUP_ORDER // s for s in stabs)
        divides = all(GROUP_ORDER % s == 0 for s in stabs)
        out[k] = {"orbits": len(codes), "sum_orbit_sizes": total, "binom": math.comb(27, k),
                  "ok": divides and total == math.comb(27, k), "stabilisers": stabs}
    return out


# ---------------------------------------------------------------- the model

class Lin:
    """A 0/1 feasibility model: variables 0..26 are K_1, 27..53 are K_2."""

    def __init__(self):
        self.n = 54
        self.rows = []  # (list of (var, coef), lo, hi)

    def var(self):
        self.n += 1
        return self.n - 1

    def add(self, terms, lo, hi):
        self.rows.append((list(terms), lo, hi))


def X(s, i):
    return (s - 1) * 27 + i


BIG = 10 ** 6


def build_model(R, bound, fix=None):
    """MODEL(R, bound): slice 0 = R (point indices); variables for slices 1 and 2."""
    R = sorted(R)
    n0 = len(R)
    Rset = set(R)
    m = Lin()
    # transversal directions (w, 1)
    for w in PTS3:
        ts = []
        for a in R:
            t = m.var()
            m.add([(t, 1), (X(1, IDX3[add(PTS3[a], w)]), -1)], -BIG, 0)
            m.add([(t, 1), (X(2, IDX3[add(PTS3[a], smul(2, w))]), -1)], -BIG, 0)
            ts.append(t)
        m.add([(t, 1) for t in ts], 1, BIG)  # empty when R is empty: infeasible row 0 >= 1
    # horizontal directions (w, 0)
    for w in DIR3:
        if any(all(i in Rset for i in L) for L in LINES3[w]):
            continue
        hs = []
        for s in (1, 2):
            for L in LINES3[w]:
                h = m.var()
                for i in L:
                    m.add([(h, 1), (X(s, i), -1)], -BIG, 0)
                hs.append(h)
        m.add([(h, 1) for h in hs], 1, BIG)
    # every hyperplane section has at least n0 points (R is a smallest section)
    for nrm in DIR4:
        alpha, beta = nrm[:3], nrm[3]
        for c in range(3):
            cnt0 = sum(1 for a in R if dot(alpha, PTS3[a]) == c)
            terms = [(X(s, i), 1) for s in (1, 2) for i in range(27)
                     if (dot(alpha, PTS3[i]) + beta * s) % 3 == c]
            if n0 - cnt0 > 0:
                m.add(terms, n0 - cnt0, BIG)
    # size
    m.add([(X(s, i), 1) for s in (1, 2) for i in range(27)], -BIG, bound - n0)
    if fix is not None:
        for s in (1, 2):
            for i in range(27):
                v = 1 if i in fix[s] else 0
                m.add([(X(s, i), 1)], v, v)
    return m


def solve_cpsat(m, seconds, workers):
    from ortools.sat.python import cp_model
    model = cp_model.CpModel()
    v = [model.NewBoolVar(f"v{i}") for i in range(m.n)]
    for terms, lo, hi in m.rows:
        if not terms:
            if lo > 0 or hi < 0:
                return "INFEASIBLE", None, 0.0
            continue
        model.AddLinearConstraint(sum(c * v[j] for j, c in terms), lo, hi)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    try:
        solver.parameters.num_workers = workers
    except Exception:
        solver.parameters.num_search_workers = workers
    t0 = time.time()
    st = solver.Solve(model)
    dt = time.time() - t0
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sol = [int(solver.Value(v[i])) for i in range(54)]
        return "FEASIBLE", sol, dt
    if st == cp_model.INFEASIBLE:
        return "INFEASIBLE", None, dt
    return "UNKNOWN", None, dt


def solve_highs(m, seconds):
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_matrix
    rows, cols, vals, lo, hi = [], [], [], [], []
    r = 0
    for terms, l, h in m.rows:
        if not terms:
            if l > 0 or h < 0:
                return "INFEASIBLE", None, 0.0
            continue
        for j, c in terms:
            rows.append(r)
            cols.append(j)
            vals.append(c)
        lo.append(-np.inf if l <= -BIG else l)
        hi.append(np.inf if h >= BIG else h)
        r += 1
    A = coo_matrix((vals, (rows, cols)), shape=(r, m.n)).tocsr()
    t0 = time.time()
    res = milp(c=np.zeros(m.n), constraints=[LinearConstraint(A, lo, hi)],
               integrality=np.ones(m.n), bounds=Bounds(0, 1),
               options={"time_limit": float(seconds), "disp": False})
    dt = time.time() - t0
    if res.x is not None and res.status in (0, 1):
        sol = [int(round(x)) for x in res.x[:54]]
        return "FEASIBLE", sol, dt
    if res.status == 2:
        return "INFEASIBLE", None, dt
    return "UNKNOWN", None, dt


def assemble(R, sol):
    K1 = [i for i in range(27) if sol[X(1, i)]]
    K2 = [i for i in range(27) if sol[X(2, i)]]
    K = ([PTS3[i] + (0,) for i in R] + [PTS3[i] + (1,) for i in K1]
         + [PTS3[i] + (2,) for i in K2])
    return K, K1, K2


def check_config(R, sol, bound):
    """Independent re-check of a solution: brute-force Kakeya, size, smallest section = |R|."""
    K, _, _ = assemble(R, sol)
    ok_k = is_kakeya(K, 4)
    smin = min(s[0] for s in sections(K))
    ok = ok_k and len(K) <= bound and smin == len(R)
    return ok, {"size": len(K), "kakeya": ok_k, "min_section": smin,
                "points": [list(x) for x in sorted(K)]}


# ---------------------------------------------------------------- normalising a known set

def normalise(K, perms):
    """Affine image of K whose slice s = 0 is a smallest section, moved to its canonical
    representative. Returns (n, code, {s: sorted indices})."""
    n, nrm, c = min(sections(K))
    j = next(i for i in range(4) if nrm[i])
    rows = [tuple(1 if r == i else 0 for r in range(4)) for i in range(4) if i != j] + [nrm]
    K2 = []
    for x in K:
        y = tuple(dot(rw, x) for rw in rows)
        K2.append((y[0], y[1], y[2], (y[3] - c) % 3))
    sl = {s: sorted(IDX3[x[:3]] for x in K2 if x[3] == s) for s in range(3)}
    code, g = canon(sl[0], perms)
    perm = perms[g]
    new = {s: sorted(int(perm[i]) for i in sl[s]) for s in range(3)}
    assert new[0] == decode(code, n), "normalisation failed"
    assert len(K2) == len(set(K2)) == len(K)
    return n, code, new


# The 27-point set of kakeya-min-v1 (munimtechnologies/kakeya-conjecture,
# data/kakeya-min-v1/results.json, entry n = 4, q = 3); re-checked below, and the downloaded
# copy is preferred when it is available.
KNOWN27 = [(0, 0, 0, 0), (1, 0, 0, 0), (2, 0, 0, 0), (2, 2, 0, 0), (0, 0, 1, 0), (1, 0, 1, 0),
           (0, 1, 1, 0), (1, 1, 1, 0), (0, 2, 1, 0), (2, 2, 1, 0), (0, 0, 2, 0), (0, 1, 2, 0),
           (1, 2, 2, 0), (1, 2, 0, 1), (1, 0, 1, 1), (0, 1, 1, 1), (1, 2, 1, 1), (0, 0, 2, 1),
           (2, 1, 2, 1), (2, 0, 0, 2), (0, 1, 0, 2), (0, 2, 0, 2), (1, 0, 1, 2), (2, 0, 1, 2),
           (1, 2, 1, 2), (0, 1, 2, 2), (2, 1, 2, 2)]
KNOWN_URL = ("https://raw.githubusercontent.com/munimtechnologies/kakeya-conjecture/main/"
             "data/kakeya-min-v1/results.json")


def classical(n):
    """The classical Kakeya set in F_3^n (2^(n+1) - 1 points), first coordinate = t."""
    K = set((t,) for t in range(3))
    for m in range(1, n):
        new = set((0,) + y for y in K)
        for d in itertools.product(range(3), repeat=m):
            p = tuple(x * x % 3 for x in d)
            for t in range(3):
                new.add((t,) + tuple((a + t * b) % 3 for a, b in zip(p, d)))
        K = new
    return sorted(K)


def fetch_known(log):
    try:
        import urllib.request
        with urllib.request.urlopen(KNOWN_URL, timeout=60) as fh:
            data = json.load(fh)
        for e in data:
            if e.get("n") == 4 and e.get("q") == 3:
                return [tuple(p) for p in e["certificate"]]
    except Exception as exc:  # network is optional
        log(f"download failed: {exc}")
    return None


def known_check(K, perms, workers, seconds):
    out = {"size": len(K), "kakeya": is_kakeya(K, 4)}
    n, code, sl = normalise(K, perms)
    out.update({"min_section": n, "rep": decode(code, n),
                "slice_sizes": [len(sl[s]) for s in range(3)]})
    fix = {1: set(sl[1]), 2: set(sl[2])}
    for b in (27, 26):
        m = build_model(sl[0], b, fix=fix)
        out[f"cpsat_fixed_{b}"] = solve_cpsat(m, seconds, workers)[0]
        out[f"highs_fixed_{b}"] = solve_highs(m, seconds)[0]
    m = build_model(sl[0], 27)
    st, sol, _ = solve_cpsat(m, seconds, workers)
    out["cpsat_free_27"] = st
    if sol is not None:
        out["cpsat_free_27_check"] = check_config(sl[0], sol, 27)[0]
    out["ok"] = (out["kakeya"] and out["size"] == 27
                 and out["cpsat_fixed_27"] == "FEASIBLE" and out["highs_fixed_27"] == "FEASIBLE"
                 and out["cpsat_fixed_26"] == "INFEASIBLE"
                 and out["highs_fixed_26"] == "INFEASIBLE"
                 and out["cpsat_free_27"] == "FEASIBLE" and out.get("cpsat_free_27_check", False))
    return out


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="kmin43.json")
    ap.add_argument("--bound", type=int, default=26)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--t-cpsat", type=float, default=900)
    ap.add_argument("--t-highs", type=float, default=600)
    ap.add_argument("--t27", type=float, default=60)
    args = ap.parse_args()
    t_start = time.time()
    res = {"job": "kmin43-v1", "question": "is there a Kakeya set of <= 26 points in F_3^4?",
           "made_with": "AI agents working in MV2's name; no human has checked this",
           "bound": args.bound, "bugs": []}
    bugs = res["bugs"]

    def log(msg):
        print(f"[{time.time() - t_start:8.1f}s] {msg}", flush=True)

    def dump():
        res["seconds"] = round(time.time() - t_start, 1)
        with open(args.out, "w") as fh:
            json.dump(res, fh, indent=1)

    perms = agl3_perms()
    nuniq = len(np.unique(perms, axis=0))
    res["group"] = {"rows": int(perms.shape[0]), "distinct": nuniq}
    if perms.shape[0] != GROUP_ORDER or nuniq != GROUP_ORDER:
        bugs.append("AGL(3,3) is not 303,264 distinct permutations")
    kmax = max(9, args.bound // 3)
    reps = orbit_reps(kmax, perms, log)
    oc = orbit_check(reps, perms)
    res["orbits"] = {str(k): {kk: vv for kk, vv in v.items() if kk != "stabilisers"}
                     for k, v in oc.items()}
    for k, v in oc.items():
        if not v["ok"]:
            bugs.append(f"orbit sizes for k={k} add up to {v['sum_orbit_sizes']}, not {v['binom']}")
    dump()

    # known 27-point set
    cands = {"embedded": KNOWN27}
    dl = fetch_known(log)
    if dl is not None:
        cands["downloaded"] = dl
        res["embedded_equals_downloaded"] = sorted(dl) == sorted(KNOWN27)
    res["known27"] = {}
    for name, K in cands.items():
        try:
            res["known27"][name] = known_check(K, perms, args.workers, max(args.t27, 180))
        except Exception as exc:
            res["known27"][name] = {"ok": False, "error": repr(exc)}
        log(f"known27 {name}: {res['known27'][name]}")
    if not any(v["ok"] for v in res["known27"].values()):
        bugs.append("the known 27-point set is not accepted by the model")
    cl = [x[1:] + x[:1] for x in classical(4)]
    res["classical"] = {"size": len(cl), "kakeya": is_kakeya(cl, 4),
                        "slice_test": slice_kakeya(cl)}
    if not (len(cl) == 31 and res["classical"]["kakeya"] and res["classical"]["slice_test"]):
        bugs.append("classical 31-point set check failed")
    dump()

    # main search: bound 26 (default), smallest section n0 = 1 .. bound // 3
    cases = []
    found = []
    for n0 in range(1, args.bound // 3 + 1):
        for code in reps[n0]:
            R = decode(code, n0)
            m = build_model(R, args.bound)
            st1, sol1, dt1 = solve_cpsat(m, args.t_cpsat, args.workers)
            st2, sol2, dt2 = solve_highs(m, args.t_highs)
            case = {"n0": n0, "rep": R, "rep_points": [list(PTS3[i]) for i in R],
                    "cpsat": st1, "cpsat_s": round(dt1, 2), "highs": st2,
                    "highs_s": round(dt2, 2)}
            for name, sol in (("cpsat", sol1), ("highs", sol2)):
                if sol is not None:
                    ok, info = check_config(R, sol, args.bound)
                    case[name + "_check"] = ok
                    if ok:
                        found.append(dict(info, solver=name, n0=n0, rep=R))
                    else:
                        bugs.append(f"{name} solution for rep {R} fails the brute-force check")
            if "INFEASIBLE" in (st1, st2) and found and found[-1]["rep"] == R:
                bugs.append(f"solvers disagree on rep {R}: {st1} / {st2}")
            cases.append(case)
            res["cases"] = cases
            res["found"] = found
            dump()
            log(f"n0={n0} rep={R}: cpsat {st1} {dt1:.1f}s, highs {st2} {dt2:.1f}s")

    # sanity at bound + 1: the solver must find sets there (the known set has 27 points)
    sanity = []
    if args.bound == 26:
        for n0 in range(1, 10):
            for code in reps[n0]:
                R = decode(code, n0)
                st, sol, dt = solve_cpsat(build_model(R, 27), args.t27, args.workers)
                ent = {"n0": n0, "rep": R, "cpsat": st, "s": round(dt, 2)}
                if sol is not None:
                    ok, info = check_config(R, sol, 27)
                    ent["check"] = ok
                    ent["size"] = info["size"]
                    if not ok:
                        bugs.append(f"bound-27 solution for rep {R} fails the brute-force check")
                sanity.append(ent)
        res["sanity27"] = {"feasible_reps": sum(1 for e in sanity if e["cpsat"] == "FEASIBLE"),
                           "infeasible_reps": sum(1 for e in sanity if e["cpsat"] == "INFEASIBLE"),
                           "unknown_reps": sum(1 for e in sanity if e["cpsat"] == "UNKNOWN"),
                           "cases": sanity}
        if res["sanity27"]["feasible_reps"] == 0:
            bugs.append("no 27-point set found at bound 27")

    all_cp = all(c["cpsat"] == "INFEASIBLE" for c in cases)
    all_hi = all(c["highs"] == "INFEASIBLE" for c in cases)
    if bugs:
        verdict = "BUG"
    elif found:
        verdict = "K43-LE-26" if args.bound == 26 else f"FOUND-LE-{args.bound}"
    elif all_cp and all_hi:
        verdict = "K43-EQ-27" if args.bound == 26 else f"NONE-LE-{args.bound}"
    elif all_cp:
        verdict = "K43-EQ-27-CPSAT-ONLY" if args.bound == 26 else "CPSAT-ONLY"
    else:
        verdict = "UNDECIDED"
    res["verdict"] = verdict
    res["cases_total"] = len(cases)
    res["open_cases"] = [c["rep"] for c in cases
                         if c["cpsat"] != "INFEASIBLE" and c["highs"] != "INFEASIBLE"
                         and not any(f["rep"] == c["rep"] for f in found)]
    dump()
    log(f"verdict {verdict}; cases {len(cases)}; found {len(found)}; bugs {bugs}")


if __name__ == "__main__":
    main()
