# kprof26.py -- exact search for a 26-point Kakeya set in F_3^4. If every case below is
# infeasible there is no Kakeya set of at most 26 points (a superset of a Kakeya set is Kakeya),
# so k_4 = 27 (the 27-point set of Lev 2009 exists): Lev's conjecture.
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0.
# Repository: munimtechnologies/kakeya-conjecture (kprof26/). Imports ../kprof25/kprof25.py.
#
# Reduction (moments/BAD25.md, N = 26). For any 26-point set K in F_3^4, over the 40 classes of
# parallel hyperplanes, sum s2 = 13 C(26,2) and sum s3 = 4 C(26,3) + 9 L (L = lines inside K).
# With delta26(cut) = s3 - 8 s2 + 576 = g(n1) + g(n2) + g(n3) + 576, g(n) = n(n-1)(n-26)/6,
# sum over the classes of delta26 = 9 (L - 40). Among the cuts allowed by (T) and (SD):
#   delta26 > 0 only for (6,6,14) +12, (5,7,14) +9, (6,7,13) +5, (7,7,12) +2;
#   delta26 = 0 exactly for (5,8,13), (6,8,12), (7,8,11), (8,8,10), (8,9,9)
#   (g(n) + g(18 - n) = -408, so every cut with an 8-point section has delta26 = 0).
# Kakeya => L >= 40. Hence a 26-point Kakeya set falls in one of these cases (after an affine map
# sending the chosen class to the last coordinate, slices sorted by size, slice 0 canonical):
#   P_i: some class has the i-th positive cut, no class has an earlier positive cut;
#   Z_j: no positive cut at all. Then every delta26 <= 0 and the sum is 9 (L - 40) >= 0, so every
#        class has a zero cut and L = 40, i.e. EXACTLY one line in each of the 40 directions;
#        some class has the j-th zero cut and no class has an earlier zero cut.
# Each case is decided by CP-SAT for every AGL(3,3)-orbit representative of slice 0, and the
# infeasible ones are re-solved by HiGHS with a different encoding (kprof25.solve_highs; the
# L = 40 / one-line-per-direction constraints of Z_j are implied by its table and identities).

import argparse
import itertools
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (HERE, os.path.join(HERE, "..", "kprof25"), "/tmp"):
    if _p not in sys.path:
        sys.path.append(_p)

import kprof25 as kp  # noqa: E402

N = 26


def g26(n):
    return n * (n - 1) * (n - 26) // 6


def delta26(p):
    s2, s3 = kp.moments_of(p)
    return s3 - 8 * s2 + 576


HAND_POS26 = [(6, 6, 14), (5, 7, 14), (6, 7, 13), (7, 7, 12)]  # by delta26: 12, 9, 5, 2
HAND_ZERO26 = [(5, 8, 13), (6, 8, 12), (7, 8, 11), (8, 8, 10), (8, 9, 9)]


def pos_cuts():
    return sorted((p for p in kp.allowed_profiles(N) if delta26(p) > 0),
                  key=lambda p: (-delta26(p), p))


def zero_cuts():
    return [p for p in kp.allowed_profiles(N) if delta26(p) == 0]


def build_cases():
    out = []
    for i, p in enumerate(HAND_POS26):
        out.append(("P%d" % (i + 1), p, list(HAND_POS26[:i]), False))
    for j, p in enumerate(HAND_ZERO26):
        out.append(("Z%d" % (j + 1), p, list(HAND_ZERO26[:j]), True))
    return out


CASES26 = build_cases()


def zero_table(forbid=()):
    forbid = set(tuple(sorted(f)) for f in forbid)
    rows = set()
    for p in HAND_ZERO26:
        if p in forbid:
            continue
        s2, s3 = kp.moments_of(p)
        for q in itertools.permutations(p):
            rows.add((q[0], q[1], q[2], s2, s3))
    return sorted(rows)


def spec_for(cut, forbid, zero_only, R, fix=None):
    a, b, c = cut
    spec = kp.make_spec(N, R, b, c, forbid=[tuple(f) for f in forbid], fix=fix)
    if zero_only:
        spec["table"] = zero_table(forbid)
        spec["L_exact"] = 40
        spec["cap"] = 1
    else:
        spec["L_exact"] = None
        spec["cap"] = None
    return spec


# ---------------------------------------------------------------- CP-SAT model (kprof25 + extras)

def solve_cpsat26(spec, seconds, workers):
    """kprof25.solve_cpsat plus two optional constraints: L == spec['L_exact'] and at most
    spec['cap'] lines of K in each direction (lines inside slice 0 included)."""
    from ortools.sat.python import cp_model
    M = cp_model.CpModel()
    x = [M.NewBoolVar(f"x{i}") for i in range(54)]
    M.Add(sum(x[kp.X(1, i)] for i in range(27)) == spec["b"])
    M.Add(sum(x[kp.X(2, i)] for i in range(27)) == spec["c"])
    if spec.get("fix") is not None:
        for s in (1, 2):
            for i in range(27):
                M.Add(x[kp.X(s, i)] == (1 if i in spec["fix"][s] else 0))
    ndir = len(kp.DIR4)
    by_dir = [[] for _ in range(ndir)]
    const_dir = [0] * ndir
    zs = []
    for d, vs in spec["lines"]:
        if not vs:
            const_dir[d] += 1
            continue
        z = M.NewBoolVar("")
        for j in vs:
            M.AddImplication(z, x[j])
        M.AddBoolOr([z] + [x[j].Not() for j in vs])
        zs.append(z)
        by_dir[d].append(z)
    cap = spec.get("cap")
    for d in range(ndir):
        if cap is not None:
            room = cap - const_dir[d]
            if room < 0:
                return "INFEASIBLE", None, 0.0
            if by_dir[d]:
                M.Add(sum(by_dir[d]) <= room)
        if const_dir[d] == 0:
            if not by_dir[d]:
                return "INFEASIBLE", None, 0.0
            M.AddBoolOr(by_dir[d])
    Lv = M.NewIntVar(ndir, len(kp.LINES4), "L")
    M.Add(Lv == sum(const_dir) + sum(zs))
    if spec.get("L_exact") is not None:
        M.Add(Lv == int(spec["L_exact"]))
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


# ---------------------------------------------------------------- brute-force predicates

def lines_per_direction(K):
    Ks = set(tuple(p) for p in K)
    cnt = [0] * len(kp.DIR4)
    for Lpts, v in kp.LINES4:
        if all(p in Ks for p in Lpts):
            cnt[kp.DIRIDX4[v]] += 1
    return cnt


def brute_case_ok(K, cut, forbid, zero_only):
    """Does the 26-point set K (already normalised: last class = cut, slices sorted) satisfy every
    condition the model of this case imposes?"""
    if len(set(K)) != N or not kp.is_kakeya(K):
        return False
    cuts = [tuple(sorted(kp.class_profile(K, u))) for u in kp.DIR4]
    allowed = set(HAND_ZERO26 if zero_only else kp.allowed_profiles(N))
    allowed -= set(tuple(sorted(f)) for f in forbid)
    if not all(p in allowed for p in cuts):
        return False
    if tuple(kp.class_profile(K, kp.LAST)) != tuple(cut):
        return False
    if zero_only and lines_per_direction(K) != [1] * len(kp.DIR4):
        return False
    return True


# ---------------------------------------------------------------- sanity checks

def extras_sanity(perms, workers, seconds):
    """The two extra constraints on the real 27-point set: the true L and the true maximum number
    of lines per direction must be feasible, one less (or one more for L) infeasible."""
    K = kp.KNOWN27
    cnt = lines_per_direction(K)
    L, mx = sum(cnt), max(cnt)
    sl = kp.normalise_class(K, kp.DIR4[0], perms)
    a, b, c = (len(sl[s]) for s in range(3))
    fix = {1: set(sl[1]), 2: set(sl[2])}

    def run(L_exact, cap):
        spec = kp.make_spec(27, sl[0], b, c, fix=fix)
        spec["L_exact"], spec["cap"] = L_exact, cap
        return solve_cpsat26(spec, seconds, workers)[0]

    out = {"L": L, "max_per_dir": mx, "min_per_dir": min(cnt),
           "L_true": run(L, None), "L_plus1": run(L + 1, None), "L_minus1": run(L - 1, None),
           "cap_true": run(None, mx), "cap_minus1": run(None, mx - 1), "both_true": run(L, mx)}
    out["ok"] = (min(cnt) >= 1 and out["L_true"] == "FEASIBLE" and out["L_plus1"] == "INFEASIBLE"
                 and out["L_minus1"] == "INFEASIBLE" and out["cap_true"] == "FEASIBLE"
                 and out["cap_minus1"] == "INFEASIBLE" and out["both_true"] == "FEASIBLE")
    return out


def removal_agreement(perms, workers, seconds, which=None):
    """For K = (known 27-point set minus one point), normalised along the class of largest s3,
    the fixed-slice CP-SAT model (all cuts allowed, and the Z1-type model) must agree with the
    brute-force predicate. Any removal that leaves a Kakeya set would itself give k_4 <= 26."""
    rows = []
    idx = range(27) if which is None else which
    for r in idx:
        K = [p for i, p in enumerate(kp.KNOWN27) if i != r]
        u0 = max(kp.DIR4, key=lambda u: (kp.moments_of(tuple(sorted(kp.class_profile(K, u))))[1], u))
        sl = kp.normalise_class(K, u0, perms)
        cut = tuple(len(sl[s]) for s in range(3))
        Kn = kp.assemble(sl[0], sl[1], sl[2])
        fix = {1: set(sl[1]), 2: set(sl[2])}
        row = {"removed": list(kp.KNOWN27[r]), "cut": list(cut), "kakeya": kp.is_kakeya(K)}
        for label, zero_only in (("all", False), ("zero", True)):
            if zero_only and cut not in HAND_ZERO26:
                row[label] = "SKIP"
                continue
            if not zero_only and cut not in set(kp.allowed_profiles(N)):
                row[label] = "SKIP"
                continue
            spec = spec_for(cut, [], zero_only, sl[0], fix=fix)
            st = solve_cpsat26(spec, seconds, workers)[0]
            want = "FEASIBLE" if brute_case_ok(Kn, cut, [], zero_only) else "INFEASIBLE"
            row[label] = st
            row[label + "_agree"] = (st == want)
        rows.append(row)
    ok = all(r.get(k, True) for r in rows for k in ("all_agree", "zero_agree"))
    return {"rows": rows, "ok": ok, "any_kakeya_26": any(r["kakeya"] for r in rows)}


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="kprof26.json")
    ap.add_argument("--job", default="kprof26")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--t1", type=float, default=90, help="first CP-SAT pass, seconds per rep")
    ap.add_argument("--cp-frac", type=float, default=0.72, help="share of the deadline for CP-SAT")
    ap.add_argument("--t-highs", type=float, default=600)
    ap.add_argument("--t-sanity", type=float, default=120)
    ap.add_argument("--deadline-min", type=float, default=330)
    ap.add_argument("--only", default="", help="comma list of case names, e.g. P1,P2")
    args = ap.parse_args()
    t_start = time.time()
    res = {"job": args.job,
           "question": "is there a Kakeya set of 26 points in F_3^4? (cases P1-P4, Z1-Z5)",
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

    def left_cp():
        return 60 * args.deadline_min * args.cp_frac - (time.time() - t_start)

    # static checks
    if pos_cuts() != HAND_POS26:
        bugs.append(f"positive cuts differ from the hand list: {pos_cuts()}")
    if zero_cuts() != HAND_ZERO26:
        bugs.append(f"zero cuts differ from the hand list: {zero_cuts()}")
    if any(g26(n) + g26(18 - n) != -408 for n in range(19)):
        bugs.append("g(n) + g(18-n) = -408 fails")
    if any(delta26(p) != g26(p[0]) + g26(p[1]) + g26(p[2]) + 576 for p in kp.allowed_profiles(N)):
        bugs.append("delta26 != sum g + 576")
    res["allowed26"] = [[list(p), delta26(p)] for p in kp.allowed_profiles(N)]

    only = set(s for s in args.only.split(",") if s)
    sel = [c for c in CASES26 if not only or c[0] in only]
    kmax = max(c[1][0] for c in sel)

    perms = kp.agl3_perms()
    reps = kp.orbit_reps(kmax, perms)
    oc = kp.orbit_check(reps, perms)
    res["orbits"] = {str(k): v for k, v in oc.items()}
    for k, v in oc.items():
        if not v["ok"]:
            bugs.append(f"orbit sizes for k={k} add up to {v['sum_orbit_sizes']}, not {v['binom']}")
    log(f"orbits: { {k: v['orbits'] for k, v in oc.items()} }")
    dump()

    for name, fn in (("known27", lambda: kp.known_sanity(kp.KNOWN27, perms, args.workers,
                                                         args.t_sanity, highs=True)),
                     ("extras", lambda: extras_sanity(perms, args.workers, args.t_sanity)),
                     ("removals", lambda: removal_agreement(perms, args.workers, 60))):
        try:
            res[name] = fn()
        except Exception as exc:
            res[name] = {"ok": False, "error": repr(exc)}
        log(f"{name}: ok={res[name].get('ok')}")
        if not res[name].get("ok"):
            bugs.append(f"sanity check {name} failed")
        dump()
    if res.get("removals", {}).get("any_kakeya_26"):
        res["note"] = "a 26-point Kakeya set obtained by deleting a point of the known 27-point set"

    cases = []
    for name, cut, forbid, zero_only in sel:
        for code in reps[cut[0]]:
            R = kp.decode(code, cut[0])
            cases.append({"case": name, "cut": list(cut), "forbid": [list(f) for f in forbid],
                          "zero_only": zero_only, "rep": R,
                          "cpsat": "NOT-RUN", "highs": "NOT-RUN"})
    res["cases"] = cases
    found = []
    res["found"] = found

    def spec_of(case):
        return spec_for(tuple(case["cut"]), case["forbid"], case["zero_only"], case["rep"])

    def record(case, solver, st, sol, dt):
        case[solver], case[solver + "_s"] = st, round(case.get(solver + "_s", 0) + dt, 2)
        if sol is not None:
            ok, info = kp.check_solution(spec_of(case), sol)
            ok = ok and kp.is_kakeya([tuple(p) for p in info["points"]])
            case[solver + "_check"] = ok
            if ok:
                found.append(dict(info, solver=solver, case=case["case"], rep=case["rep"]))
            else:
                bugs.append(f"{solver} solution for case {case['case']} rep {case['rep']} "
                            f"fails the brute-force check: {info}")

    # pass 1: short CP-SAT on every rep
    for case in cases:
        if left_cp() < 60:
            break
        st, sol, dt = solve_cpsat26(spec_of(case), min(args.t1, left_cp()), args.workers)
        record(case, "cpsat", st, sol, dt)
        dump()
        log(f"pass 1 {case['case']} rep {case['rep']}: cpsat {st} {dt:.1f}s")
    # pass 2: the rest of the CP-SAT budget shared equally among the open reps
    open_ = [c for c in cases if c["cpsat"] not in ("INFEASIBLE", "FEASIBLE")]
    for k, case in enumerate(open_):
        if left_cp() < 60:
            break
        share = max(60.0, left_cp() / (len(open_) - k))
        st, sol, dt = solve_cpsat26(spec_of(case), min(share, left_cp()), args.workers)
        record(case, "cpsat", st, sol, dt)
        dump()
        log(f"pass 2 {case['case']} rep {case['rep']}: cpsat {st} {dt:.1f}s")
    # pass 3: HiGHS on every rep CP-SAT decided as infeasible (independent second solver)
    for case in cases:
        if case["cpsat"] != "INFEASIBLE" or left() < 120:
            continue
        st, sol, dt = kp.solve_highs(spec_of(case), min(args.t_highs, max(60, left() - 60)))
        record(case, "highs", st, sol, dt)
        if st == "FEASIBLE":
            bugs.append(f"solvers disagree on case {case['case']} rep {case['rep']}")
        dump()
        log(f"pass 3 {case['case']} rep {case['rep']}: highs {st} {dt:.1f}s")

    per = {}
    for name, cut, _, _ in sel:
        cb = [c for c in cases if c["case"] == name]
        per[name] = {"cut": list(cut), "reps": len(cb),
                     "cpsat_infeasible": sum(1 for c in cb if c["cpsat"] == "INFEASIBLE"),
                     "highs_infeasible": sum(1 for c in cb if c["highs"] == "INFEASIBLE"),
                     "cpsat_seconds": round(sum(c.get("cpsat_s", 0) for c in cb), 1),
                     "open": [c["rep"] for c in cb if c["cpsat"] != "INFEASIBLE"]}
    res["per_case"] = per
    all_cp = bool(cases) and all(c["cpsat"] == "INFEASIBLE" for c in cases)
    all_hi = all_cp and all(c["highs"] == "INFEASIBLE" for c in cases)
    whole = not only or only >= set(c[0] for c in CASES26)
    if bugs:
        verdict = "BUG"
    elif found:
        verdict = "FOUND-26"
    elif all_hi:
        verdict = "K4-EQ-27" if whole else "PART-CLOSED"
    elif all_cp:
        verdict = ("K4-EQ-27" if whole else "PART-CLOSED") + "-CPSAT-ONLY"
    else:
        verdict = "UNDECIDED"
    res["verdict"] = verdict
    dump()
    log(f"verdict {verdict}; per case {json.dumps(per)}; bugs {bugs}")


if __name__ == "__main__":
    main()
