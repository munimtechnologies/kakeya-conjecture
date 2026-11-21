# kmin34.py -- the smallest Kakeya set in F_4^3.  Known: 21 <= K(3,4) <= 28.
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0.
#
# F_4 = {0, 1, 2, 3} with 2 = w, 3 = w + 1, w^2 = w + 1 (the encoding of kakeya.py in
# munimtechnologies/kakeya-conjecture): addition is XOR, 2*2 = 3, 2*3 = 1, 3*3 = 2.
#
# Reduction (kmin34/prereg-kmin34.md in munimtechnologies/kakeya-conjecture).
# Write F_4^3 = F_4^2 x F_4 with points (y, s) and slices K_s = {y : (y, s) in K}.
#  * A direction (w, 1), w in F_4^2 (16 "transversal" directions), has a line in K iff some
#    a in K_0 has a + t w in K_t for t = 1, 2, 3.
#  * A direction (w, 0), w in PG(1, 4) (5 "horizontal" directions), has a line in K iff some
#    slice K_s contains a line of F_4^2 in direction w.
#  * Let K be a Kakeya set with |K| <= b and let H be an affine plane with |K n H| = n0
#    minimal (n0 <= b // 4 because four parallel planes partition F_4^3; n0 >= 1 because a
#    transversal line meets every parallel plane). Move H to {s = 0} by an affine map of
#    F_4^3, then move K_0 to its canonical representative R under AGL(2, 4) acting on y alone
#    ((y, s) -> (A y + c, s) is affine on F_4^3). The image is again a Kakeya set, so
#    (K_1, K_2, K_3) solves MODEL(R, b): slice 0 equal to R, transversal and horizontal
#    conditions, every one of the 84 plane sections >= |R|, and |K| <= b.
#  * So if MODEL(R, b) is infeasible for every orbit representative R with 1 <= |R| <= b // 4,
#    then K(3, 4) > b.
# Each model has 48 point variables. It is solved by OR-Tools CP-SAT and, independently, by
# HiGHS (scipy.optimize.milp). Every solution found is re-checked by a brute-force Kakeya test.
# If a set of size m <= b is found, the search is repeated at bound m - 1 (a representative
# proved infeasible by both solvers at a larger bound stays infeasible: the model only gets
# tighter as the bound drops).
#
# Self-checks: AGL(2,4) has 2,880 distinct permutations; orbit sizes add up to C(16, k); the
# known 28-point set (kakeya-min-v1) is Kakeya, and after normalisation both solvers accept it
# at bound 28 and reject the same fixed set at bound 27; the free model for its
# representative is feasible at bound 28.

import argparse
import itertools
import json
import math
import time

import numpy as np

Q = 4
ADD = [[a ^ b for b in range(Q)] for a in range(Q)]
_EXP = [1, 2, 3]  # w^0, w^1, w^2
_LOG = {1: 0, 2: 1, 3: 2}
MUL = [[0 if a == 0 or b == 0 else _EXP[(_LOG[a] + _LOG[b]) % 3] for b in range(Q)]
       for a in range(Q)]

PTS2 = list(itertools.product(range(Q), repeat=2))
IDX2 = {p: i for i, p in enumerate(PTS2)}
GROUP_ORDER = 2880  # |AGL(2, 4)| = 16 * 180


def add(p, q):
    return tuple(a ^ b for a, b in zip(p, q))


def smul(c, p):
    return tuple(MUL[c][a] for a in p)


def dot(u, v):
    s = 0
    for a, b in zip(u, v):
        s ^= MUL[a][b]
    return s


def projective_points(n):
    """Directions of F_4^n: non-zero vectors whose first non-zero coordinate is 1."""
    out = []
    for v in itertools.product(range(Q), repeat=n):
        nz = [x for x in v if x]
        if nz and nz[0] == 1:
            out.append(v)
    return out


DIR2 = projective_points(2)  # 5
DIR3 = projective_points(3)  # 21


def lines2(w):
    """The 4 lines of F_4^2 in direction w, each as a sorted tuple of 4 point indices."""
    seen = set()
    for a in PTS2:
        seen.add(tuple(sorted(IDX2[add(a, smul(t, w))] for t in range(Q))))
    return sorted(seen)


LINES2 = {w: lines2(w) for w in DIR2}


def is_kakeya(K, n):
    """Brute force: every direction of F_4^n has a full line inside K."""
    K = set(tuple(x) for x in K)
    for v in projective_points(n):
        ok = False
        for x in K:
            if all(add(x, smul(t, v)) in K for t in range(1, Q)):
                ok = True
                break
        if not ok:
            return False
    return True


def slice_kakeya(K):
    """The transversal + horizontal test on the slices by the LAST coordinate (for tests)."""
    S = {s: set(tuple(x[:2]) for x in K if x[2] == s) for s in range(Q)}
    for w in PTS2:
        if not any(all(add(a, smul(t, w)) in S[t] for t in range(1, Q)) for a in S[0]):
            return False
    for w in DIR2:
        ok = False
        for s in range(Q):
            for L in LINES2[w]:
                if all(PTS2[i] in S[s] for i in L):
                    ok = True
        if not ok:
            return False
    return True


def sections(K):
    """All 84 affine plane sections of K in F_4^3: (size, normal, c)."""
    out = []
    for nrm in DIR3:
        for c in range(Q):
            out.append((sum(1 for x in K if dot(nrm, x) == c), nrm, c))
    return out


# ---------------------------------------------------------------- the group AGL(2, 4)

def det2(M):
    (a, b), (c, d) = M
    return MUL[a][d] ^ MUL[b][c]  # characteristic 2: minus is plus


def agl2_perms():
    """Array (2880, 16): row g maps point index i to the index of g(PTS2[i])."""
    lin = []
    for ent in itertools.product(range(Q), repeat=4):
        M = (ent[0:2], ent[2:4])
        if det2(M):
            lin.append([IDX2[tuple(dot(M[r], p) for r in range(2))] for p in PTS2])
    lin = np.array(lin, dtype=np.int16)
    T = np.array([[IDX2[add(PTS2[i], b)] for i in range(16)] for b in PTS2], dtype=np.int16)
    perms = T[:, lin]  # [b, m, i] = T[b, lin[m, i]]: first the linear map, then + b
    return perms.reshape(-1, 16)


def codes_of(S, perms):
    S = sorted(S)
    k = len(S)
    img = np.sort(perms[:, S].astype(np.int64), axis=1)
    w = np.array([16 ** (k - 1 - i) for i in range(k)], dtype=np.int64)
    return img @ w


def encode(S):
    S = sorted(S)
    k = len(S)
    return sum(int(x) * 16 ** (k - 1 - i) for i, x in enumerate(S))


def decode(code, k):
    out = []
    for _ in range(k):
        out.append(code % 16)
        code //= 16
    return out[::-1]


def canon(S, perms):
    """(canonical code, index of a group element g with sorted g(S) = decode(code))."""
    c = codes_of(S, perms)
    g = int(np.argmin(c))
    return int(c[g]), g


def stabiliser(S, perms):
    return int((codes_of(S, perms) == encode(S)).sum())


def orbit_reps(kmax, perms, log=None):
    """Canonical codes of the AGL(2,4)-orbits of k-subsets of F_4^2, k = 1..kmax."""
    reps = {1: [canon([0], perms)[0]]}
    for k in range(2, kmax + 1):
        seen = set()
        for code in reps[k - 1]:
            S = decode(code, k - 1)
            for p in range(16):
                if p not in S:
                    seen.add(canon(S + [p], perms)[0])
        reps[k] = sorted(seen)
        if log:
            log(f"k={k}: {len(reps[k])} orbits")
    return reps


def orbit_check(reps, perms):
    """For each k: (sum of orbit sizes, C(16, k)) and the stabiliser orders."""
    out = {}
    for k, codes in reps.items():
        stabs = [stabiliser(decode(c, k), perms) for c in codes]
        total = sum(GROUP_ORDER // s for s in stabs)
        divides = all(GROUP_ORDER % s == 0 for s in stabs)
        out[k] = {"orbits": len(codes), "sum_orbit_sizes": total, "binom": math.comb(16, k),
                  "ok": divides and total == math.comb(16, k), "stabilisers": stabs}
    return out


# ---------------------------------------------------------------- the model

NV = 48  # point variables: slices 1, 2, 3 of 16 points each


class Lin:
    """A 0/1 feasibility model: variables 0..47 are K_1, K_2, K_3; then auxiliary ones.
    aux[j] lists the point variables that auxiliary variable j is bounded by (j <= x)."""

    def __init__(self):
        self.n = NV
        self.rows = []  # (list of (var, coef), lo, hi)
        self.aux = {}

    def var(self, needs):
        j = self.n
        self.n += 1
        self.aux[j] = list(needs)
        for x in needs:
            self.rows.append(([(j, 1), (x, -1)], -BIG, 0))
        return j

    def add(self, terms, lo, hi):
        self.rows.append((list(terms), lo, hi))


def X(s, i):
    return (s - 1) * 16 + i


BIG = 10 ** 6


def build_model(R, bound, fix=None):
    """MODEL(R, bound): slice 0 = R (point indices); variables for slices 1, 2, 3."""
    R = sorted(R)
    n0 = len(R)
    Rset = set(R)
    m = Lin()
    # transversal directions (w, 1): a in R, a + t w in K_t for t = 1, 2, 3
    for w in PTS2:
        ts = []
        for a in R:
            needs = [X(t, IDX2[add(PTS2[a], smul(t, w))]) for t in range(1, Q)]
            ts.append(m.var(needs))
        m.add([(t, 1) for t in ts], 1, BIG)  # empty when R is empty: infeasible row 0 >= 1
    # horizontal directions (w, 0)
    for w in DIR2:
        if any(all(i in Rset for i in L) for L in LINES2[w]):
            continue
        hs = []
        for s in range(1, Q):
            for L in LINES2[w]:
                hs.append(m.var([X(s, i) for i in L]))
        m.add([(h, 1) for h in hs], 1, BIG)
    # every plane section has at least n0 points (R is a smallest section)
    for nrm in DIR3:
        alpha, beta = nrm[:2], nrm[2]
        for c in range(Q):
            cnt0 = sum(1 for a in R if dot(alpha, PTS2[a]) == c)
            terms = [(X(s, i), 1) for s in range(1, Q) for i in range(16)
                     if dot(alpha, PTS2[i]) ^ MUL[beta][s] == c]
            if n0 - cnt0 > 0:
                m.add(terms, n0 - cnt0, BIG)
    # size
    m.add([(X(s, i), 1) for s in range(1, Q) for i in range(16)], -BIG, bound - n0)
    if fix is not None:
        for s in range(1, Q):
            for i in range(16):
                v = 1 if i in fix[s] else 0
                m.add([(X(s, i), 1)], v, v)
    return m


def evaluate(m, point_vals):
    """Solver-free check: with the 48 point values given, set every auxiliary variable to 1
    exactly when all its point variables are 1 (its largest allowed value) and report
    whether every row holds. A 0/1 point assignment extends to a solution of the model iff
    this returns True (auxiliaries only appear with coefficient +1 in >= rows)."""
    v = list(point_vals) + [0] * (m.n - NV)
    for j, needs in m.aux.items():
        v[j] = 1 if all(v[x] for x in needs) else 0
    for terms, lo, hi in m.rows:
        tot = sum(c * v[j] for j, c in terms)
        if tot < lo or tot > hi:
            return False
    return True


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
        sol = [int(solver.Value(v[i])) for i in range(NV)]
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
        sol = [int(round(x)) for x in res.x[:NV]]
        return "FEASIBLE", sol, dt
    if res.status == 2:
        return "INFEASIBLE", None, dt
    return "UNKNOWN", None, dt


def assemble(R, sol):
    sl = {s: [i for i in range(16) if sol[X(s, i)]] for s in range(1, Q)}
    K = [PTS2[i] + (0,) for i in R]
    for s in range(1, Q):
        K += [PTS2[i] + (s,) for i in sl[s]]
    return K, sl


def check_config(R, sol, bound):
    """Independent re-check of a solution: brute-force Kakeya, size, smallest section = |R|."""
    K, _ = assemble(R, sol)
    ok_k = is_kakeya(K, 3)
    smin = min(s[0] for s in sections(K))
    ok = ok_k and len(K) <= bound and smin == len(R) and len(set(K)) == len(K)
    return ok, {"size": len(K), "kakeya": ok_k, "min_section": smin,
                "points": [list(x) for x in sorted(K)]}


def point_vals_of(sl):
    vals = [0] * NV
    for s in range(1, Q):
        for i in sl[s]:
            vals[X(s, i)] = 1
    return vals


# ---------------------------------------------------------------- normalising a known set

def normalise(K, perms):
    """Affine image of K whose slice s = 0 is a smallest section, moved to its canonical
    representative. Returns (n, code, {s: sorted indices})."""
    n, nrm, c = min(sections(K))
    j = next(i for i in range(3) if nrm[i])
    rows = [tuple(1 if r == i else 0 for r in range(3)) for i in range(3) if i != j] + [nrm]
    K2 = []
    for x in K:
        y = tuple(dot(rw, x) for rw in rows)
        K2.append((y[0], y[1], y[2] ^ c))
    sl = {s: sorted(IDX2[x[:2]] for x in K2 if x[2] == s) for s in range(Q)}
    code, g = canon(sl[0], perms)
    perm = perms[g]
    new = {s: sorted(int(perm[i]) for i in sl[s]) for s in range(Q)}
    assert new[0] == decode(code, n), "normalisation failed"
    assert len(K2) == len(set(K2)) == len(K)
    return n, code, new


def from_slices(sl):
    K = []
    for s in range(Q):
        K += [PTS2[i] + (s,) for i in sl[s]]
    return K


# The 28-point set of kakeya-min-v1 (munimtechnologies/kakeya-conjecture,
# data/kakeya-min-v1/results.json, entry n = 3, q = 4); checked by hand to be Kakeya, re-checked
# below, and the downloaded copy is preferred when it is available.
KNOWN28 = [(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0), (2, 2, 0), (3, 2, 0),
           (0, 0, 1), (3, 0, 1), (0, 1, 1), (1, 1, 1), (0, 2, 1), (2, 2, 1), (0, 3, 1),
           (1, 3, 1), (2, 3, 1), (3, 3, 1),
           (0, 1, 2), (1, 1, 2), (2, 1, 2), (3, 1, 2), (1, 2, 2), (2, 2, 2),
           (0, 2, 3), (1, 2, 3), (2, 2, 3), (3, 2, 3), (1, 3, 3), (3, 3, 3)]
KNOWN_URL = ("https://raw.githubusercontent.com/munimtechnologies/kakeya-conjecture/main/"
             "data/kakeya-min-v1/results.json")


def fetch_known(log):
    try:
        import urllib.request
        with urllib.request.urlopen(KNOWN_URL, timeout=60) as fh:
            data = json.load(fh)
        for e in data:
            if e.get("n") == 3 and e.get("q") == 4:
                return [tuple(p) for p in e["certificate"]]
    except Exception as exc:  # network is optional
        log(f"download failed: {exc}")
    return None


def known_check(K, perms, workers, seconds):
    out = {"size": len(K), "kakeya": is_kakeya(K, 3)}
    n, code, sl = normalise(K, perms)
    out.update({"min_section": n, "rep": decode(code, n),
                "slice_sizes": [len(sl[s]) for s in range(Q)]})
    fix = {s: set(sl[s]) for s in range(1, Q)}
    vals = point_vals_of(sl)
    out["evaluate_28"] = evaluate(build_model(sl[0], 28), vals)
    out["evaluate_27"] = evaluate(build_model(sl[0], 27), vals)
    for b in (28, 27):
        m = build_model(sl[0], b, fix=fix)
        out[f"cpsat_fixed_{b}"] = solve_cpsat(m, seconds, workers)[0]
        out[f"highs_fixed_{b}"] = solve_highs(m, seconds)[0]
    m = build_model(sl[0], 28)
    st, sol, _ = solve_cpsat(m, seconds, workers)
    out["cpsat_free_28"] = st
    if sol is not None:
        out["cpsat_free_28_check"] = check_config(sl[0], sol, 28)[0]
    out["ok"] = (out["kakeya"] and out["size"] == 28
                 and out["evaluate_28"] and not out["evaluate_27"]
                 and out["cpsat_fixed_28"] == "FEASIBLE" and out["highs_fixed_28"] == "FEASIBLE"
                 and out["cpsat_fixed_27"] == "INFEASIBLE"
                 and out["highs_fixed_27"] == "INFEASIBLE"
                 and out["cpsat_free_28"] == "FEASIBLE" and out.get("cpsat_free_28_check", False))
    return out


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="kmin34.json")
    ap.add_argument("--bound", type=int, default=27)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--t-cpsat", type=float, default=900)
    ap.add_argument("--t-highs", type=float, default=300)
    ap.add_argument("--t-sanity", type=float, default=30)
    ap.add_argument("--deadline-min", type=float, default=330)
    args = ap.parse_args()
    t_start = time.time()
    res = {"job": "kmin34-v1", "question": "what is the smallest Kakeya set in F_4^3? (21..28)",
           "made_with": "AI agents working in MV2's name; no human has checked this",
           "start_bound": args.bound, "bugs": []}
    bugs = res["bugs"]

    def log(msg):
        print(f"[{time.time() - t_start:8.1f}s] {msg}", flush=True)

    def dump():
        res["seconds"] = round(time.time() - t_start, 1)
        with open(args.out, "w") as fh:
            json.dump(res, fh, indent=1)

    def out_of_time():
        return time.time() - t_start > 60 * args.deadline_min

    perms = agl2_perms()
    nuniq = len(np.unique(perms, axis=0))
    res["group"] = {"rows": int(perms.shape[0]), "distinct": nuniq}
    if perms.shape[0] != GROUP_ORDER or nuniq != GROUP_ORDER:
        bugs.append("AGL(2,4) is not 2,880 distinct permutations")
    kmax = max(7, args.bound // 4)
    reps = orbit_reps(kmax, perms, log)
    oc = orbit_check(reps, perms)
    res["orbits"] = {str(k): {kk: vv for kk, vv in v.items() if kk != "stabilisers"}
                     for k, v in oc.items()}
    for k, v in oc.items():
        if not v["ok"]:
            bugs.append(f"orbit sizes for k={k} add up to {v['sum_orbit_sizes']}, not {v['binom']}")
    dump()

    # known 28-point set
    cands = {"embedded": KNOWN28}
    dl = fetch_known(log)
    if dl is not None:
        cands["downloaded"] = dl
        res["embedded_equals_downloaded"] = sorted(dl) == sorted(KNOWN28)
    res["known28"] = {}
    for name, K in cands.items():
        try:
            res["known28"][name] = known_check(K, perms, args.workers, 180)
        except Exception as exc:
            res["known28"][name] = {"ok": False, "error": repr(exc)}
        log(f"known28 {name}: {res['known28'][name]}")
    if not any(v["ok"] for v in res["known28"].values()):
        bugs.append("the known 28-point set is not accepted by the model")
    dump()

    # main search: bound b, smallest section n0 = 1 .. b // 4; lower b after every find
    infeasible_both = {}  # rep code -> bound at which both solvers proved infeasibility
    infeasible_cp = {}
    passes = []
    found = []
    bound = args.bound
    timed_out = False
    while bound >= 1 and not timed_out:
        cases = []
        hit = None
        for n0 in range(1, bound // 4 + 1):
            for code in reps[n0]:
                R = decode(code, n0)
                key = f"{n0}:{code}"
                if key in infeasible_both:
                    cases.append({"n0": n0, "rep": R, "cpsat": "INFEASIBLE", "highs": "INFEASIBLE",
                                  "inherited_from_bound": infeasible_both[key]})
                    continue
                if out_of_time():
                    timed_out = True
                    cases.append({"n0": n0, "rep": R, "cpsat": "NOT-RUN", "highs": "NOT-RUN"})
                    continue
                m = build_model(R, bound)
                st1, sol1, dt1 = solve_cpsat(m, args.t_cpsat, args.workers)
                st2, sol2, dt2 = solve_highs(m, args.t_highs)
                case = {"n0": n0, "rep": R, "rep_points": [list(PTS2[i]) for i in R],
                        "cpsat": st1, "cpsat_s": round(dt1, 2), "highs": st2,
                        "highs_s": round(dt2, 2)}
                good = []
                for name, sol in (("cpsat", sol1), ("highs", sol2)):
                    if sol is not None:
                        ok, info = check_config(R, sol, bound)
                        case[name + "_check"] = ok
                        if ok and evaluate(m, sol):
                            good.append(dict(info, solver=name, n0=n0, rep=R, bound=bound))
                        else:
                            bugs.append(f"{name} solution for rep {R} at bound {bound} fails "
                                        "the brute-force or row check")
                if good and "INFEASIBLE" in (st1, st2):
                    bugs.append(f"solvers disagree on rep {R} at bound {bound}: {st1} / {st2}")
                if st1 == "INFEASIBLE" and st2 == "INFEASIBLE":
                    infeasible_both[key] = bound
                if st1 == "INFEASIBLE":
                    infeasible_cp.setdefault(key, bound)
                found += good
                cases.append(case)
                log(f"bound {bound} n0={n0} rep={R}: cpsat {st1} {dt1:.1f}s, "
                    f"highs {st2} {dt2:.1f}s")
                if good:
                    hit = min(g["size"] for g in good)
                    break
            if hit is not None:
                break
        passes.append({"bound": bound, "cases": cases, "hit_size": hit})
        res["passes"] = passes
        res["found"] = found
        dump()
        if hit is None:
            break
        bound = hit - 1

    last = passes[-1] if passes else {"bound": args.bound, "cases": [], "hit_size": None}
    all_both = last["hit_size"] is None and all(
        c["cpsat"] == "INFEASIBLE" and c["highs"] == "INFEASIBLE" for c in last["cases"])
    all_cp = last["hit_size"] is None and all(c["cpsat"] == "INFEASIBLE" for c in last["cases"])
    best = min((f["size"] for f in found), default=28)
    res["best_size"] = best
    res["last_bound"] = last["bound"]

    # sanity at bound 28 (only from the start bound): some representative must admit a set
    sanity = []
    if args.bound == 27 and not out_of_time():
        for n0 in range(1, 8):
            for code in reps[n0]:
                R = decode(code, n0)
                st, sol, dt = solve_cpsat(build_model(R, 28), args.t_sanity, args.workers)
                ent = {"n0": n0, "rep": R, "cpsat": st, "s": round(dt, 2)}
                if sol is not None:
                    ok, info = check_config(R, sol, 28)
                    ent["check"] = ok
                    ent["size"] = info["size"]
                    if not ok:
                        bugs.append(f"bound-28 solution for rep {R} fails the brute-force check")
                sanity.append(ent)
        res["sanity28"] = {"feasible_reps": sum(1 for e in sanity if e["cpsat"] == "FEASIBLE"),
                           "infeasible_reps": sum(1 for e in sanity if e["cpsat"] == "INFEASIBLE"),
                           "unknown_reps": sum(1 for e in sanity if e["cpsat"] == "UNKNOWN"),
                           "cases": sanity}
        if res["sanity28"]["feasible_reps"] == 0:
            bugs.append("no 28-point set found at bound 28")

    if bugs:
        verdict = "BUG"
    elif all_both:
        verdict = f"K34-EQ-{best}"
    elif all_cp:
        verdict = f"K34-EQ-{best}-CPSAT-ONLY"
    elif found:
        verdict = f"K34-LE-{best}"
    else:
        verdict = "UNDECIDED"
    res["verdict"] = verdict
    res["cases_total"] = sum(len(p["cases"]) for p in passes)
    res["open_cases"] = [{"bound": last["bound"], "rep": c["rep"]} for c in last["cases"]
                         if c["cpsat"] != "INFEASIBLE" and c["highs"] != "INFEASIBLE"]
    dump()
    log(f"verdict {verdict}; best {best}; passes {len(passes)}; bugs {bugs}")


if __name__ == "__main__":
    main()
