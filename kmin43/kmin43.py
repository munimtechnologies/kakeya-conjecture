# kmin43.py (v3) -- exact search: is there a Kakeya set of <= 25 (then <= 26) points in F_3^4?
# Known: 25 <= k_4 <= 27 (lower bound: moments/SD.md in munimtechnologies/kakeya-conjecture).
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0.
#
# v3 changes against v2 (the model and the reduction are unchanged):
#  (1) bound 25 is decided first (a "no" proves k_4 >= 26), then bound 26 (a "no" proves k_4 = 27);
#  (2) the smallest-section size n0 only runs up to nmax(bound), from the identity
#        sum over the 120 hyperplanes H of |K n H|^2 = 40 N + 13 N (N - 1)
#      (each point lies in 40 hyperplanes, each pair of points in 13). If every one of the 40
#      parallel classes had all three sections >= m, each class would contribute at most
#      2 m^2 + (N - 2m)^2, so some section is < m as soon as 40 (2 m^2 + (N - 2m)^2) < S.
#      This gives n0 <= 6 for N = 25 and n0 <= 7 for N = 26 (nmax() computes it, max over N' <= bound);
#  (3) symmetry breaking |K_1| <= |K_2| (the linear map (y, s) -> (y, -s) fixes slice 0 and swaps
#      slices 1 and 2), used in the main search only;
#  (4) CP-SAT runs over all cases first; HiGHS re-solves them afterwards as an independent check;
#  (5) the long bound-27 sanity loop of v2 is replaced by the known 27-point set check
#      (fixed slices feasible at 27 and infeasible at 26 in both solvers; free search at 27 finds a
#      verified set with and without symmetry breaking); 2 workers, 1 GB.
#
# Reduction (kmin43/prereg-kmin43.md in munimtechnologies/kakeya-conjecture).
# Write F_3^4 = F_3^3 x F_3 with points (y, s) and slices K_s = {y : (y, s) in K}.
#  * A direction (w, 1), w in F_3^3 (27 "transversal" directions), has a line in K iff
#    some a in K_0 has a + w in K_1 and a + 2w in K_2.
#  * A direction (w, 0), w in PG(2, 3) (13 "horizontal" directions), has a line in K iff some
#    slice K_s contains a line of F_3^3 in direction w.
#  * Let K be a Kakeya set with |K| <= bound and let H be an affine hyperplane with |K n H| = n0
#    minimal. Move H to {s = 0} by an affine map of F_3^4, then move K_0 to its canonical
#    representative R under AGL(3, 3) acting on y alone. The image is again a Kakeya set, so
#    (K_1, K_2) solves MODEL(R, bound): slice 0 equal to R, transversal and horizontal
#    conditions, every one of the 120 hyperplane sections >= |R|, and |K| <= bound.
#  * So if MODEL(R, bound) is infeasible for every orbit representative R with
#    1 <= |R| <= nmax(bound), there is no Kakeya set with <= bound points.
# Every solution found is re-checked by a brute-force Kakeya test.

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


def nmax_for_size(N):
    """Largest m such that a set of N points in F_3^4 can have all 120 sections >= m, by the
    second-moment identity (see the header). 0 if N < 3."""
    S = 40 * N + 13 * N * (N - 1)
    best = 0
    for m in range(1, N // 3 + 1):
        if 40 * (2 * m * m + (N - 2 * m) ** 2) >= S:
            best = m
    return best


def nmax(bound):
    return max(nmax_for_size(N) for N in range(1, bound + 1))


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


def build_model(R, bound, fix=None, sym=False):
    """MODEL(R, bound): slice 0 = R (point indices); variables for slices 1 and 2.
    sym=True adds |K_1| <= |K_2| (valid symmetry breaking, see the header)."""
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
    if sym:
        m.add([(X(1, i), 1) for i in range(27)] + [(X(2, i), -1) for i in range(27)], -BIG, 0)
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
    """The classical Kakeya set in F_3^n (2^(n+1) - 1 points), first coordinate = t.

    K_1 = F_3. Given a Kakeya set K_m in F_3^m (first coordinate t), K_{m+1} is
      {(t, d^2 + t d) : t in F_3, d in F_3^m}  (squares coordinatewise; a line in every
      direction (1, d)), together with {0} x phi(K_m), where phi(t, y) = (t - 1, 2y) is affine.
    """
    K = set((t,) for t in range(3))
    for m in range(1, n):
        new = set((0, (x[0] - 1) % 3) + tuple((2 * a) % 3 for a in x[1:]) for x in K)
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
    for sym in (False, True):
        m = build_model(sl[0], 27, sym=sym)
        st, sol, _ = solve_cpsat(m, seconds, workers)
        key = "cpsat_free_27" + ("_sym" if sym else "")
        out[key] = st
        if sol is not None:
            out[key + "_check"] = check_config(sl[0], sol, 27)[0]
    out["ok"] = (out["kakeya"] and out["size"] == 27
                 and out["cpsat_fixed_27"] == "FEASIBLE" and out["highs_fixed_27"] == "FEASIBLE"
                 and out["cpsat_fixed_26"] == "INFEASIBLE"
                 and out["highs_fixed_26"] == "INFEASIBLE"
                 and out["cpsat_free_27"] == "FEASIBLE" and out.get("cpsat_free_27_check", False)
                 and out["cpsat_free_27_sym"] == "FEASIBLE"
                 and out.get("cpsat_free_27_sym_check", False))
    return out


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="kmin43.json")
    ap.add_argument("--bounds", default="25,26")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--t-cpsat", type=float, default=1200)
    ap.add_argument("--t-highs", type=float, default=900)
    ap.add_argument("--t27", type=float, default=180)
    ap.add_argument("--deadline-min", type=float, default=320)
    args = ap.parse_args()
    bounds = [int(b) for b in args.bounds.split(",")]
    t_start = time.time()
    res = {"job": "kmin43-v3", "question": "is there a Kakeya set of <= 25, or <= 26, points in F_3^4?",
           "made_with": "AI agents working in MV2's name; no human has checked this",
           "bounds": bounds, "nmax": {str(b): nmax(b) for b in bounds}, "bugs": []}
    bugs = res["bugs"]

    def log(msg):
        print(f"[{time.time() - t_start:8.1f}s] {msg}", flush=True)

    def dump():
        res["seconds"] = round(time.time() - t_start, 1)
        with open(args.out, "w") as fh:
            json.dump(res, fh, indent=1)

    def out_of_time():
        return time.time() - t_start > 60 * args.deadline_min

    if res["nmax"].get("25", 6) != 6 or res["nmax"].get("26", 7) != 7:
        bugs.append(f"nmax differs from the hand values 6 (N=25) and 7 (N=26): {res['nmax']}")

    perms = agl3_perms()
    nuniq = len(np.unique(perms, axis=0))
    res["group"] = {"rows": int(perms.shape[0]), "distinct": nuniq}
    if perms.shape[0] != GROUP_ORDER or nuniq != GROUP_ORDER:
        bugs.append("AGL(3,3) is not 303,264 distinct permutations")
    kmax = max(nmax(b) for b in bounds)
    reps = orbit_reps(kmax, perms, log)
    oc = orbit_check(reps, perms)
    res["orbits"] = {str(k): {kk: vv for kk, vv in v.items() if kk != "stabilisers"}
                     for k, v in oc.items()}
    for k, v in oc.items():
        if not v["ok"]:
            bugs.append(f"orbit sizes for k={k} add up to {v['sum_orbit_sizes']}, not {v['binom']}")
    dump()

    # known 27-point set (sanity: the model must accept it)
    cands = {"embedded": KNOWN27}
    dl = fetch_known(log)
    if dl is not None:
        cands["downloaded"] = dl
        res["embedded_equals_downloaded"] = sorted(dl) == sorted(KNOWN27)
    res["known27"] = {}
    for name, K in cands.items():
        try:
            res["known27"][name] = known_check(K, perms, args.workers, args.t27)
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

    # pass 1: CP-SAT on every case, bound by bound
    cases = []
    found = []
    for b in bounds:
        for n0 in range(1, nmax(b) + 1):
            for code in reps[n0]:
                R = decode(code, n0)
                case = {"bound": b, "n0": n0, "rep": R, "rep_points": [list(PTS3[i]) for i in R],
                        "cpsat": "NOT-RUN", "highs": "NOT-RUN"}
                cases.append(case)
                if out_of_time():
                    continue
                st, sol, dt = solve_cpsat(build_model(R, b, sym=True), args.t_cpsat, args.workers)
                case["cpsat"], case["cpsat_s"] = st, round(dt, 2)
                if sol is not None:
                    ok, info = check_config(R, sol, b)
                    case["cpsat_check"] = ok
                    if ok:
                        found.append(dict(info, solver="cpsat", bound=b, n0=n0, rep=R))
                    else:
                        bugs.append(f"cpsat solution for bound {b} rep {R} fails the brute-force check")
                res["cases"], res["found"] = cases, found
                dump()
                log(f"bound={b} n0={n0} rep={R}: cpsat {st} {dt:.1f}s")

    # pass 2: HiGHS re-solves the same cases (independent second solver)
    for case in cases:
        if out_of_time() or case["cpsat"] == "NOT-RUN":
            continue
        R, b = case["rep"], case["bound"]
        st, sol, dt = solve_highs(build_model(R, b, sym=True), args.t_highs)
        case["highs"], case["highs_s"] = st, round(dt, 2)
        if sol is not None:
            ok, info = check_config(R, sol, b)
            case["highs_check"] = ok
            if ok:
                found.append(dict(info, solver="highs", bound=b, n0=case["n0"], rep=R))
            else:
                bugs.append(f"highs solution for bound {b} rep {R} fails the brute-force check")
        if (st == "FEASIBLE" and case["cpsat"] == "INFEASIBLE") or \
                (st == "INFEASIBLE" and case["cpsat"] == "FEASIBLE"):
            bugs.append(f"solvers disagree on bound {b} rep {R}: {case['cpsat']} / {st}")
        res["cases"], res["found"] = cases, found
        dump()
        log(f"bound={b} n0={case['n0']} rep={R}: highs {st} {dt:.1f}s")

    per = {}
    for b in bounds:
        cb = [c for c in cases if c["bound"] == b]
        per[str(b)] = {
            "cases": len(cb),
            "cpsat_all_infeasible": all(c["cpsat"] == "INFEASIBLE" for c in cb),
            "highs_all_infeasible": all(c["highs"] == "INFEASIBLE" for c in cb),
            "cpsat_not_run": sum(1 for c in cb if c["cpsat"] == "NOT-RUN"),
            "highs_not_run": sum(1 for c in cb if c["highs"] == "NOT-RUN"),
            "open_cases": [c["rep"] for c in cb if c["cpsat"] != "INFEASIBLE"
                           and c["highs"] != "INFEASIBLE"],
            "found_sizes": sorted({f["size"] for f in found if f["bound"] == b}),
        }
    res["per_bound"] = per
    if any(f["size"] <= 24 for f in found):
        bugs.append("a verified Kakeya set of <= 24 points contradicts the proof of k_4 >= 25 (SD.md)")
    p25, p26 = per.get("25"), per.get("26")
    if bugs:
        verdict = "BUG"
    elif found:
        verdict = "FOUND-" + str(min(f["size"] for f in found))
    elif p26 and p26["cpsat_all_infeasible"] and p26["highs_all_infeasible"]:
        verdict = "K4-EQ-27"
    elif p26 and p26["cpsat_all_infeasible"]:
        verdict = "K4-EQ-27-CPSAT-ONLY"
    elif p25 and p25["cpsat_all_infeasible"] and p25["highs_all_infeasible"]:
        verdict = "K4-GE-26"
    elif p25 and p25["cpsat_all_infeasible"]:
        verdict = "K4-GE-26-CPSAT-ONLY"
    else:
        verdict = "UNDECIDED"
    res["verdict"] = verdict
    res["cases"], res["found"] = cases, found
    dump()
    log(f"verdict {verdict}; per bound {json.dumps(per)}; bugs {bugs}")


if __name__ == "__main__":
    main()
