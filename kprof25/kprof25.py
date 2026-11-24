# kprof25.py -- exact search for a 25-point Kakeya set in F_3^4, narrowed by the cut lemma of
# moments/BAD25.md (munimtechnologies/kakeya-conjecture).
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0.
#
# Reduction (BAD25.md). If a Kakeya set with <= 25 points exists, one with exactly 25 exists (add
# points). For a 25-point Kakeya set, sum over the 40 classes of delta = 63 L - 2420 >= 100 with
# delta = 7 s3 - 51 s2 + 3302, and among the 14 cuts allowed by (T) and (SD) only (6,6,13) (+76),
# (5,7,13) (+60) and (6,7,12) (+25) have delta > 0. Hence, after an affine map, one of:
#   A: last class cut (6,6,13)                                      slice 0 = canonical 6-set
#   B: last class cut (5,7,13), no class (6,6,13)                   slice 0 = canonical 5-set
#   C: last class cut (6,7,12), no class (6,6,13) or (5,7,13)       slice 0 = canonical 6-set
# with |K_1| = middle size and |K_2| = largest size (s -> -s fixes slice 0 and swaps 1 and 2).
#
# Model, for each case and each AGL(3,3)-orbit representative R of slice 0 (54 Boolean point
# variables for slices 1 and 2):
#   * a Boolean z_l for each of the 1080 lines of F_3^4 that can lie in K, z_l = AND of its points;
#     every one of the 40 directions has a line (Kakeya); L = number of lines;
#   * the three section sizes of each of the 40 classes form an allowed cut (table, all orders),
#     with that cut's s2 and s3 as two more table columns;
#   * sum of s2 = 13 C(25,2) and sum of s3 = 4 C(25,3) + 9 L (identities true for every set);
#   * |K_1| = b, |K_2| = c.
# All of these hold for every 25-point Kakeya set in the normalised position, so if every model is
# infeasible there is no Kakeya set of 25 points, i.e. k_4 >= 26. Every solution found is
# re-checked by brute force. CP-SAT decides every case; HiGHS (scipy.optimize.milp) re-solves them
# as an independent second solver with a different (0/1 profile-choice) encoding of the table.

import argparse
import itertools
import json
import math
import time
from math import comb

import numpy as np

PTS3 = list(itertools.product(range(3), repeat=3))
IDX3 = {p: i for i, p in enumerate(PTS3)}
PTS4 = list(itertools.product(range(3), repeat=4))
GROUP_ORDER = 303264  # |AGL(3, 3)| = 27 * 11232


def add(p, q):
    return tuple((a + b) % 3 for a, b in zip(p, q))


def smul(c, p):
    return tuple((c * a) % 3 for a in p)


def dot(u, v):
    return sum(a * b for a, b in zip(u, v)) % 3


def projective_points(n):
    out = []
    for v in itertools.product(range(3), repeat=n):
        nz = [x for x in v if x]
        if nz and nz[0] == 1:
            out.append(v)
    return out


DIR3 = projective_points(3)  # 13
DIR4 = projective_points(4)  # 40
DIRIDX4 = {v: i for i, v in enumerate(DIR4)}
LAST = (0, 0, 0, 1)


def X(s, i):
    """Variable index of point (PTS3[i], s) for s = 1, 2."""
    return (s - 1) * 27 + i


def lines4():
    seen = {}
    for v in DIR4:
        for x in PTS4:
            L = tuple(sorted(add(x, smul(t, v)) for t in range(3)))
            if L not in seen:
                seen[L] = v
    return sorted(seen.items())


LINES4 = lines4()  # 1080 (line, direction)


def is_kakeya(K, n=4):
    K = set(tuple(x) for x in K)
    for v in projective_points(n):
        if not any(add(x, v) in K and add(x, smul(2, v)) in K for x in K):
            return False
    return True


def class_profile(K, nrm):
    cnt = [0, 0, 0]
    for x in K:
        cnt[dot(nrm, x)] += 1
    return tuple(cnt)


def count_lines(K):
    Ks = set(tuple(x) for x in K)
    return sum(1 for L, _ in LINES4 if all(p in Ks for p in L))


# ---------------------------------------------------------------- allowed cuts ((T) and (SD))

def energy_cap(x, y, T=27):
    m = min(x, y)
    total = x * y
    if total < T or total > T * m:
        return None
    e = total - T
    if e == 0:
        return T
    q, rem = divmod(e, m - 1)
    ones = T - q - (1 if rem else 0)
    return ones + q * m * m + ((1 + rem) ** 2 if rem else 0)


def sd_floor(x, y, T=27):
    E = energy_cap(x, y, T)
    if E is None:
        return None
    return -(-(T * T) // E)


def allowed_profiles(N):
    out = []
    for a in range(0, N // 3 + 1):
        for b in range(a, (N - a) // 2 + 1):
            c = N - a - b
            if c < b:
                continue
            p = (a, b, c)
            if all(sd_floor(p[i], p[j]) is not None and p[k] >= sd_floor(p[i], p[j])
                   for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0))):
                out.append(p)
    return out


def moments_of(p):
    return sum(comb(c, 2) for c in p), sum(comb(c, 3) for c in p)


def delta25(p):
    s2, s3 = moments_of(p)
    return 7 * s3 - 51 * s2 + 3302


HAND25 = [(4, 9, 12), (4, 10, 11), (5, 7, 13), (5, 8, 12), (5, 9, 11), (5, 10, 10), (6, 6, 13),
          (6, 7, 12), (6, 8, 11), (6, 9, 10), (7, 7, 11), (7, 8, 10), (7, 9, 9), (8, 8, 9)]
CASES = [  # (name, cut of the last class, forbidden cuts); C first: it is the most constrained
    ("C", (6, 7, 12), [(6, 6, 13), (5, 7, 13)]),
    ("B", (5, 7, 13), [(6, 6, 13)]),
    ("A", (6, 6, 13), []),
]


def table(N, forbid=()):
    forbid = set(tuple(sorted(f)) for f in forbid)
    rows = set()
    for p in allowed_profiles(N):
        if p in forbid:
            continue
        s2, s3 = moments_of(p)
        for q in itertools.permutations(p):
            rows.add((q[0], q[1], q[2], s2, s3))
    return sorted(rows)


# ---------------------------------------------------------------- the group AGL(3, 3)

def det3(M):
    (a, b, c), (d, e, f), (g, h, i) = M
    return (a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)) % 3


def agl3_perms():
    lin = []
    for ent in itertools.product(range(3), repeat=9):
        M = (ent[0:3], ent[3:6], ent[6:9])
        if det3(M):
            lin.append([IDX3[tuple(dot(M[r], p) for r in range(3))] for p in PTS3])
    lin = np.array(lin, dtype=np.int16)
    T = np.array([[IDX3[add(PTS3[i], b)] for i in range(27)] for b in PTS3], dtype=np.int16)
    perms = T[:, lin]
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
    c = codes_of(S, perms)
    g = int(np.argmin(c))
    return int(c[g]), g


def stabiliser(S, perms):
    return int((codes_of(S, perms) == encode(S)).sum())


def orbit_reps(kmax, perms):
    reps = {1: [canon([0], perms)[0]]}
    for k in range(2, kmax + 1):
        seen = set()
        for code in reps[k - 1]:
            S = decode(code, k - 1)
            for p in range(27):
                if p not in S:
                    seen.add(canon(S + [p], perms)[0])
        reps[k] = sorted(seen)
    return reps


def orbit_check(reps, perms):
    out = {}
    for k, codes in reps.items():
        stabs = [stabiliser(decode(c, k), perms) for c in codes]
        total = sum(GROUP_ORDER // s for s in stabs)
        divides = all(GROUP_ORDER % s == 0 for s in stabs)
        out[k] = {"orbits": len(codes), "sum_orbit_sizes": total, "binom": math.comb(27, k),
                  "ok": divides and total == math.comb(27, k)}
    return out


def normalise_class(K, nrm, perms):
    """Affine image of K in which class nrm becomes the last-coordinate class, with sections of
    nondecreasing size at s = 0, 1, 2 and slice 0 canonical under AGL(3,3). Returns {s: indices}."""
    sizes = class_profile(K, nrm)
    order = sorted(range(3), key=lambda cc: (sizes[cc], cc))
    c0, mid = order[0], order[1]
    sign = 1 if (mid - c0) % 3 == 1 else 2
    j = next(i for i in range(4) if nrm[i])
    keep = [i for i in range(4) if i != j]
    sl = {0: [], 1: [], 2: []}
    for x in K:
        y = tuple(x[i] for i in keep)
        s = (sign * (dot(nrm, x) - c0)) % 3
        sl[s].append(IDX3[y])
    code, g = canon(sl[0], perms)
    perm = perms[g]
    new = {s: sorted(int(perm[i]) for i in sl[s]) for s in range(3)}
    assert new[0] == decode(code, len(sl[0])), "canonical form mismatch"
    assert [len(new[s]) for s in range(3)] == sorted(sizes), "slice sizes not sorted"
    return new


def assemble(R, K1, K2):
    return ([PTS3[i] + (0,) for i in R] + [PTS3[i] + (1,) for i in K1]
            + [PTS3[i] + (2,) for i in K2])


# ---------------------------------------------------------------- the model

def make_spec(N, R, b, c, forbid=(), fix=None):
    R = sorted(int(i) for i in R)
    Rset = set(R)
    a = len(R)
    if a + b + c != N:
        raise ValueError("sizes do not add up to N")
    lines = []  # (direction index, list of variables); [] = line inside slice 0 = R
    for L, v in LINES4:
        vs = []
        ok = True
        for p in L:
            i, s = IDX3[p[:3]], p[3]
            if s == 0:
                if i not in Rset:
                    ok = False
                    break
            else:
                vs.append(X(s, i))
        if ok:
            lines.append((DIRIDX4[v], vs))
    secs = []
    for nrm in DIR4:
        trip = []
        for cc in range(3):
            const = sum(1 for i in R if dot(nrm[:3], PTS3[i]) == cc)
            vs = [X(s, i) for s in (1, 2) for i in range(27)
                  if (dot(nrm[:3], PTS3[i]) + nrm[3] * s) % 3 == cc]
            trip.append((const, vs))
        secs.append(trip)
    return {"N": N, "R": R, "a": a, "b": b, "c": c, "forbid": [tuple(sorted(f)) for f in forbid],
            "lines": lines, "secs": secs, "table": table(N, forbid), "fix": fix,
            "s2_total": 13 * comb(N, 2), "s3_base": 4 * comb(N, 3)}


def solve_cpsat(spec, seconds, workers):
    from ortools.sat.python import cp_model
    M = cp_model.CpModel()
    x = [M.NewBoolVar(f"x{i}") for i in range(54)]
    M.Add(sum(x[X(1, i)] for i in range(27)) == spec["b"])
    M.Add(sum(x[X(2, i)] for i in range(27)) == spec["c"])
    if spec["fix"] is not None:
        for s in (1, 2):
            for i in range(27):
                M.Add(x[X(s, i)] == (1 if i in spec["fix"][s] else 0))
    by_dir = [[] for _ in DIR4]
    always = [False] * len(DIR4)
    zs = []
    const_lines = 0
    for d, vs in spec["lines"]:
        if not vs:
            const_lines += 1
            always[d] = True
            continue
        z = M.NewBoolVar("")
        for j in vs:
            M.AddImplication(z, x[j])
        M.AddBoolOr([z] + [x[j].Not() for j in vs])
        zs.append(z)
        by_dir[d].append(z)
    for d in range(len(DIR4)):
        if always[d]:
            continue
        if not by_dir[d]:
            return "INFEASIBLE", None, 0.0
        M.AddBoolOr(by_dir[d])
    Lv = M.NewIntVar(len(DIR4), len(LINES4), "L")
    M.Add(Lv == const_lines + sum(zs))
    s2s, s3s = [], []
    for trip in spec["secs"]:
        sz = []
        for const, vs in trip:
            v = M.NewIntVar(0, 81, "")
            M.Add(v == const + sum(x[j] for j in vs))
            sz.append(v)
        s2 = M.NewIntVar(0, 100000, "")
        s3 = M.NewIntVar(0, 1000000, "")
        M.AddAllowedAssignments(sz + [s2, s3], spec["table"])
        s2s.append(s2)
        s3s.append(s3)
    M.Add(sum(s2s) == spec["s2_total"])
    M.Add(sum(s3s) == spec["s3_base"] + 9 * Lv)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    try:
        solver.parameters.num_workers = workers
    except Exception:
        solver.parameters.num_search_workers = workers
    t0 = time.time()
    st = solver.Solve(M)
    dt = time.time() - t0
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return "FEASIBLE", [int(solver.Value(x[i])) for i in range(54)], dt
    if st == cp_model.INFEASIBLE:
        return "INFEASIBLE", None, dt
    return "UNKNOWN", None, dt


def solve_highs(spec, seconds):
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_matrix
    INF = np.inf
    rows = []  # (terms, lo, hi)
    nvar = [54]

    def new():
        nvar[0] += 1
        return nvar[0] - 1

    rows.append(([(X(1, i), 1) for i in range(27)], spec["b"], spec["b"]))
    rows.append(([(X(2, i), 1) for i in range(27)], spec["c"], spec["c"]))
    if spec["fix"] is not None:
        for s in (1, 2):
            for i in range(27):
                v = 1 if i in spec["fix"][s] else 0
                rows.append(([(X(s, i), 1)], v, v))
    by_dir = [[] for _ in DIR4]
    always = [False] * len(DIR4)
    zs = []
    const_lines = 0
    for d, vs in spec["lines"]:
        if not vs:
            const_lines += 1
            always[d] = True
            continue
        z = new()
        for j in vs:
            rows.append(([(z, 1), (j, -1)], -INF, 0))
        rows.append(([(z, 1)] + [(j, -1) for j in vs], -(len(vs) - 1), INF))
        zs.append(z)
        by_dir[d].append(z)
    for d in range(len(DIR4)):
        if always[d]:
            continue
        if not by_dir[d]:
            return "INFEASIBLE", None, 0.0
        rows.append(([(z, 1) for z in by_dir[d]], 1, INF))
    Lv = new()
    rows.append(([(Lv, 1)] + [(z, -1) for z in zs], const_lines, const_lines))
    s2terms, s3terms = [], []
    for trip in spec["secs"]:
        ys = [(new(), t) for t in spec["table"]]
        rows.append(([(y, 1) for y, _ in ys], 1, 1))
        for h in range(3):
            const, vs = trip[h]
            rows.append(([(y, t[h]) for y, t in ys if t[h]] + [(j, -1) for j in vs], const, const))
        s2terms += [(y, t[3]) for y, t in ys]
        s3terms += [(y, t[4]) for y, t in ys]
    rows.append((s2terms, spec["s2_total"], spec["s2_total"]))
    rows.append((s3terms + [(Lv, -9)], spec["s3_base"], spec["s3_base"]))
    n = nvar[0]
    ri, ci, va, lo, hi = [], [], [], [], []
    for r, (terms, l, h) in enumerate(rows):
        for j, cf in terms:
            ri.append(r)
            ci.append(j)
            va.append(cf)
        lo.append(l)
        hi.append(h)
    A = coo_matrix((va, (ri, ci)), shape=(len(rows), n)).tocsr()
    lb = np.zeros(n)
    ub = np.ones(n)
    lb[Lv], ub[Lv] = len(DIR4), len(LINES4)
    t0 = time.time()
    res = milp(c=np.zeros(n), constraints=[LinearConstraint(A, np.array(lo, dtype=float),
                                                            np.array(hi, dtype=float))],
               integrality=np.ones(n), bounds=Bounds(lb, ub),
               options={"time_limit": float(seconds), "disp": False})
    dt = time.time() - t0
    if res.x is not None and res.status in (0, 1):
        return "FEASIBLE", [int(round(v)) for v in res.x[:54]], dt
    if res.status == 2:
        return "INFEASIBLE", None, dt
    return "UNKNOWN", None, dt


def check_solution(spec, sol):
    """Brute-force re-check: Kakeya, size N, last class = (a, b, c), every cut allowed and not
    forbidden."""
    K1 = [i for i in range(27) if sol[X(1, i)]]
    K2 = [i for i in range(27) if sol[X(2, i)]]
    K = assemble(spec["R"], K1, K2)
    allowed = set(allowed_profiles(spec["N"])) - set(spec["forbid"])
    cuts = [tuple(sorted(class_profile(K, nrm))) for nrm in DIR4]
    info = {"size": len(set(K)), "kakeya": is_kakeya(K), "last": list(class_profile(K, LAST)),
            "cuts_ok": all(p in allowed for p in cuts), "lines": count_lines(K),
            "points": [list(p) for p in sorted(K)]}
    ok = (info["kakeya"] and info["size"] == spec["N"] == len(K) and info["cuts_ok"]
          and tuple(info["last"]) == (spec["a"], spec["b"], spec["c"]))
    return ok, info


# ---------------------------------------------------------------- sanity on the known 27-point set

KNOWN27 = [(0, 0, 0, 0), (1, 0, 0, 0), (2, 0, 0, 0), (2, 2, 0, 0), (0, 0, 1, 0), (1, 0, 1, 0),
           (0, 1, 1, 0), (1, 1, 1, 0), (0, 2, 1, 0), (2, 2, 1, 0), (0, 0, 2, 0), (0, 1, 2, 0),
           (1, 2, 2, 0), (1, 2, 0, 1), (1, 0, 1, 1), (0, 1, 1, 1), (1, 2, 1, 1), (0, 0, 2, 1),
           (2, 1, 2, 1), (2, 0, 0, 2), (0, 1, 0, 2), (0, 2, 0, 2), (1, 0, 1, 2), (2, 0, 1, 2),
           (1, 2, 1, 2), (0, 1, 2, 2), (2, 1, 2, 2)]


def known_sanity(K, perms, workers, seconds, highs=True):
    """The machinery must accept a real Kakeya set: fixed slices feasible, a forbidden cut that K
    really has makes it infeasible, and a free search with K's slice 0 and slice sizes finds a set
    that passes the brute-force check."""
    out = {"size": len(set(K)), "kakeya": is_kakeya(K)}
    profs = {u: tuple(sorted(class_profile(K, u))) for u in DIR4}
    allowed = set(allowed_profiles(len(K)))
    out["cuts_allowed"] = all(p in allowed for p in profs.values())
    L = count_lines(K)
    out["lines"] = L
    out["identities"] = (sum(moments_of(p)[0] for p in profs.values()) == 13 * comb(len(K), 2)
                         and sum(moments_of(p)[1] for p in profs.values())
                         == 4 * comb(len(K), 3) + 9 * L)
    u0 = max(DIR4, key=lambda u: (moments_of(profs[u])[1], u))
    sl = normalise_class(K, u0, perms)
    a, b, c = (len(sl[s]) for s in range(3))
    out["chosen_cut"] = [a, b, c]
    out["renormalised_kakeya"] = is_kakeya(assemble(sl[0], sl[1], sl[2]))
    fix = {1: set(sl[1]), 2: set(sl[2])}
    N = len(K)
    out["cpsat_fixed"] = solve_cpsat(make_spec(N, sl[0], b, c, fix=fix), seconds, workers)[0]
    others = sorted(set(p for p in profs.values() if p != (a, b, c)))
    if others:
        neg = make_spec(N, sl[0], b, c, forbid=[others[0]], fix=fix)
        out["forbidden_cut"] = list(others[0])
        out["cpsat_forbid"] = solve_cpsat(neg, seconds, workers)[0]
    free = make_spec(N, sl[0], b, c)
    st, sol, _ = solve_cpsat(free, seconds, workers)
    out["cpsat_free"] = st
    out["cpsat_free_check"] = check_solution(free, sol)[0] if sol is not None else False
    if highs:
        out["highs_fixed"] = solve_highs(make_spec(N, sl[0], b, c, fix=fix), seconds)[0]
        if others:
            out["highs_forbid"] = solve_highs(neg, seconds)[0]
    out["ok"] = (out["kakeya"] and out["size"] == 27 and out["cuts_allowed"] and out["identities"]
                 and out["renormalised_kakeya"] and out["cpsat_fixed"] == "FEASIBLE"
                 and out.get("cpsat_forbid", "INFEASIBLE") == "INFEASIBLE"
                 and out["cpsat_free"] == "FEASIBLE" and out["cpsat_free_check"]
                 and (not highs or (out["highs_fixed"] == "FEASIBLE"
                                    and out.get("highs_forbid", "INFEASIBLE") == "INFEASIBLE")))
    return out


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="kprof25.json")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--t1", type=float, default=240, help="first CP-SAT pass, seconds per case")
    ap.add_argument("--t2", type=float, default=2400, help="second CP-SAT pass on UNKNOWN cases")
    ap.add_argument("--t-highs", type=float, default=600)
    ap.add_argument("--t-sanity", type=float, default=180)
    ap.add_argument("--deadline-min", type=float, default=320)
    ap.add_argument("--only", default="", help="comma list of case names, e.g. C,B")
    args = ap.parse_args()
    t_start = time.time()
    res = {"job": "kprof25-v1",
           "question": "is there a Kakeya set of 25 points in F_3^4? (cases A, B, C of BAD25.md)",
           "made_with": "AI agents working in MV2's name; no human has checked this", "bugs": []}
    bugs = res["bugs"]

    def log(msg):
        print(f"[{time.time() - t_start:8.1f}s] {msg}", flush=True)

    def dump():
        res["seconds"] = round(time.time() - t_start, 1)
        with open(args.out, "w") as fh:
            json.dump(res, fh, indent=1)

    def left():
        return 60 * args.deadline_min - (time.time() - t_start)

    # static checks
    if allowed_profiles(25) != HAND25:
        bugs.append(f"allowed_profiles(25) differs from the hand list: {allowed_profiles(25)}")
    pos = sorted(p for p in HAND25 if delta25(p) > 0)
    if pos != sorted(c[1] for c in CASES):
        bugs.append(f"positive-delta cuts differ: {pos}")
    if len(LINES4) != 1080:
        bugs.append(f"{len(LINES4)} lines, not 1080")

    perms = agl3_perms()
    nuniq = len(np.unique(perms, axis=0))
    res["group"] = {"rows": int(perms.shape[0]), "distinct": nuniq}
    if perms.shape[0] != GROUP_ORDER or nuniq != GROUP_ORDER:
        bugs.append("AGL(3,3) is not 303,264 distinct permutations")
    reps = orbit_reps(6, perms)
    oc = orbit_check(reps, perms)
    res["orbits"] = {str(k): v for k, v in oc.items()}
    for k, v in oc.items():
        if not v["ok"]:
            bugs.append(f"orbit sizes for k={k} add up to {v['sum_orbit_sizes']}, not {v['binom']}")
    log(f"orbits: { {k: v['orbits'] for k, v in oc.items()} }")
    dump()

    try:
        res["known27"] = known_sanity(KNOWN27, perms, args.workers, args.t_sanity, highs=True)
    except Exception as exc:
        res["known27"] = {"ok": False, "error": repr(exc)}
    log(f"known27: {res['known27']}")
    if not res["known27"].get("ok"):
        bugs.append("the known 27-point set is not handled correctly by the model")
    dump()

    only = set(s for s in args.only.split(",") if s)
    cases = []
    for name, cut, forbid in CASES:
        if only and name not in only:
            continue
        a, b, c = cut
        for code in reps[a]:
            R = decode(code, a)
            cases.append({"case": name, "cut": list(cut), "forbid": [list(f) for f in forbid],
                          "rep": R, "rep_points": [list(PTS3[i]) for i in R],
                          "cpsat": "NOT-RUN", "highs": "NOT-RUN"})
    res["cases"] = cases
    found = []
    res["found"] = found

    def spec_of(case):
        a, b, c = case["cut"]
        return make_spec(25, case["rep"], b, c, forbid=[tuple(f) for f in case["forbid"]])

    def record(case, solver, st, sol, dt):
        case[solver], case[solver + "_s"] = st, round(case.get(solver + "_s", 0) + dt, 2)
        if sol is not None:
            ok, info = check_solution(spec_of(case), sol)
            case[solver + "_check"] = ok
            if ok:
                found.append(dict(info, solver=solver, case=case["case"], rep=case["rep"]))
            else:
                bugs.append(f"{solver} solution for case {case['case']} rep {case['rep']} "
                            f"fails the brute-force check: {info}")

    # pass 1 (short) and pass 1b (long, UNKNOWN only): CP-SAT
    for tlim, label in ((args.t1, "1"), (args.t2, "1b")):
        for case in cases:
            if case["cpsat"] in ("INFEASIBLE", "FEASIBLE"):
                continue
            if left() < 120:
                break
            st, sol, dt = solve_cpsat(spec_of(case), min(tlim, max(60, left() - 60)), args.workers)
            record(case, "cpsat", st, sol, dt)
            dump()
            log(f"pass {label} case {case['case']} rep {case['rep']}: cpsat {st} {dt:.1f}s")
    # pass 2: HiGHS on every case CP-SAT decided as infeasible (independent check)
    for case in cases:
        if case["cpsat"] != "INFEASIBLE" or left() < 120:
            continue
        st, sol, dt = solve_highs(spec_of(case), min(args.t_highs, max(60, left() - 60)))
        record(case, "highs", st, sol, dt)
        if st == "FEASIBLE":
            bugs.append(f"solvers disagree on case {case['case']} rep {case['rep']}")
        dump()
        log(f"pass 2 case {case['case']} rep {case['rep']}: highs {st} {dt:.1f}s")

    per = {}
    for name, _, _ in CASES:
        cb = [c for c in cases if c["case"] == name]
        per[name] = {"reps": len(cb),
                     "cpsat_infeasible": sum(1 for c in cb if c["cpsat"] == "INFEASIBLE"),
                     "highs_infeasible": sum(1 for c in cb if c["highs"] == "INFEASIBLE"),
                     "open": [c["rep"] for c in cb if c["cpsat"] != "INFEASIBLE"]}
    res["per_case"] = per
    all_cp = bool(cases) and not only and all(c["cpsat"] == "INFEASIBLE" for c in cases)
    all_hi = all_cp and all(c["highs"] == "INFEASIBLE" for c in cases)
    if bugs:
        verdict = "BUG"
    elif found:
        verdict = "FOUND-25"
    elif all_hi:
        verdict = "K4-GE-26"
    elif all_cp:
        verdict = "K4-GE-26-CPSAT-ONLY"
    else:
        verdict = "UNDECIDED"
    res["verdict"] = verdict
    dump()
    log(f"verdict {verdict}; per case {json.dumps(per)}; bugs {bugs}")


if __name__ == "__main__":
    main()
