# kdrat26.py -- independent SAT re-check of the kprof26 cases (is there a 26-point Kakeya set in
# F_3^4?). Each case/representative of kprof26 is written as a CNF, solved by CaDiCaL, and every
# UNSAT answer is checked by drat-trim on CaDiCaL's DRAT proof. If all six live cases (P1-P4,
# Z1, Z2; Z3-Z5 are impossible by the arithmetic remark in kprof26/README.md) are UNSAT with
# verified proofs, then no 26-point Kakeya set exists and k_4 = 27, with a checked certificate
# for the search part (the reduction to the cases is the hand argument in kprof26/README.md).
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0.
# Repository: munimtechnologies/kakeya-conjecture. Uses kprof25.py (geometry, orbit
# representatives, cut tables, model spec) and kprof26.py (cases, zero table) unchanged; the
# solver and the proof checker are independent of kprof25/26's CP-SAT and HiGHS.
#
# CNF encoding of a spec (variables 1..54 are the points of slices 1 and 2, var = kp.X(s,i) + 1):
#   - slice sizes b, c: exact totalizers, unit clauses;
#   - each line not killed by slice 0: z <-> AND(points) (or the point itself if only one is free);
#   - each direction: OR of its lines (unless a line lies inside slice 0); in Z cases at most one;
#   - each of the 40 classes: exact totalizers on two of its sections; every pair of section
#     sizes (q0, q1) that is not in the allowed table is forbidden by one clause (q2 = N - q0 - q1).
# The moment identities are not encoded (they hold for every point set, so leaving them out keeps
# the CNF a correct model; it only makes the solver's job harder).

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (HERE, os.path.join(HERE, "..", "kprof25"), os.path.join(HERE, "..", "kprof26"), "/tmp"):
    if _p not in sys.path:
        sys.path.append(_p)

import kprof25 as kp  # noqa: E402
import kprof26 as k26  # noqa: E402


# ---------------------------------------------------------------- CNF building blocks

class CNF:
    def __init__(self, nvars=0):
        self.n = nvars
        self.clauses = []
        self.trivial_unsat = False

    def new(self):
        self.n += 1
        return self.n

    def add(self, cl):
        cl = list(cl)
        if not cl:
            self.trivial_unsat = True
        self.clauses.append(cl)

    def dimacs(self):
        out = ["p cnf %d %d" % (self.n, len(self.clauses))]
        out.extend(" ".join(str(l) for l in c) + " 0" if c else "0" for c in self.clauses)
        return "\n".join(out) + "\n"


def totalizer(cnf, lits):
    """Exact unary counter: returns U with U[k-1] <=> (number of true lits >= k), k = 1..len."""
    lits = list(lits)
    if not lits:
        return []
    if len(lits) == 1:
        return [lits[0]]
    h = len(lits) // 2
    A = totalizer(cnf, lits[:h])
    B = totalizer(cnf, lits[h:])
    p, q = len(A), len(B)
    R = [cnf.new() for _ in range(p + q)]
    for i in range(p + 1):
        for j in range(q + 1):
            if i + j >= 1:  # (A >= i) and (B >= j) -> R >= i + j
                cl = [R[i + j - 1]]
                if i >= 1:
                    cl.append(-A[i - 1])
                if j >= 1:
                    cl.append(-B[j - 1])
                cnf.add(cl)
            if i + j + 1 <= p + q:  # (A <= i) and (B <= j) -> not (R >= i + j + 1)
                cl = [-R[i + j]]
                if i + 1 <= p:
                    cl.append(A[i])
                if j + 1 <= q:
                    cl.append(B[j])
                cnf.add(cl)
    return R


def neq_lits(U, t):
    """Literals whose OR says 'count != t' (U exact unary outputs). None if count can never be t."""
    m = len(U)
    if t < 0 or t > m:
        return None
    lits = []
    if t >= 1:
        lits.append(-U[t - 1])
    if t + 1 <= m:
        lits.append(U[t])
    return lits


def equal_units(cnf, U, t):
    m = len(U)
    if t < 0 or t > m:
        cnf.add([])
        return
    if t >= 1:
        cnf.add([U[t - 1]])
    if t + 1 <= m:
        cnf.add([-U[t]])


def build_cnf(spec):
    cnf = CNF(54)

    def xv(j):
        return j + 1

    for s, want in ((1, spec["b"]), (2, spec["c"])):
        equal_units(cnf, totalizer(cnf, [xv(kp.X(s, i)) for i in range(27)]), want)
    if spec.get("fix") is not None:
        for s in (1, 2):
            for i in range(27):
                v = xv(kp.X(s, i))
                cnf.add([v] if i in spec["fix"][s] else [-v])
    ndir = len(kp.DIR4)
    by_dir = [[] for _ in range(ndir)]
    const_dir = [0] * ndir
    for d, vs in spec["lines"]:
        if not vs:
            const_dir[d] += 1
            continue
        lits = [xv(j) for j in vs]
        if len(lits) == 1:
            z = lits[0]
        else:
            z = cnf.new()
            for l in lits:
                cnf.add([-z, l])
            cnf.add([z] + [-l for l in lits])
        by_dir[d].append(z)
    cap = spec.get("cap")
    for d in range(ndir):
        if const_dir[d] == 0:
            cnf.add(by_dir[d])
        if cap is not None:
            room = cap - const_dir[d]
            if room < 0:
                cnf.add([])
            elif room == 0:
                for z in by_dir[d]:
                    cnf.add([-z])
            elif room == 1:
                zs = by_dir[d]
                for i in range(len(zs)):
                    for j in range(i + 1, len(zs)):
                        cnf.add([-zs[i], -zs[j]])
            else:
                U = totalizer(cnf, by_dir[d])
                if room < len(U):
                    cnf.add([-U[room]])
    pairs = set((r[0], r[1]) for r in spec["table"])
    for trip in spec["secs"]:
        (c0, v0), (c1, v1) = trip[0], trip[1]
        U0 = totalizer(cnf, [xv(j) for j in v0])
        U1 = totalizer(cnf, [xv(j) for j in v1])
        for t0 in range(len(U0) + 1):
            for t1 in range(len(U1) + 1):
                if (c0 + t0, c1 + t1) in pairs:
                    continue
                cnf.add(neq_lits(U0, t0) + neq_lits(U1, t1))
    return cnf


def unit_propagate(clauses, nvars):
    """Plain unit propagation. Returns 'CONFLICT', 'SAT' (every clause satisfied) or 'OPEN'."""
    val = [0] * (nvars + 1)
    changed = True
    while changed:
        changed = False
        for cl in clauses:
            last, nun, sat = 0, 0, False
            for l in cl:
                v = val[abs(l)]
                if v == 0:
                    nun += 1
                    last = l
                elif (v > 0) == (l > 0):
                    sat = True
                    break
            if sat:
                continue
            if nun == 0:
                return "CONFLICT", val
            if nun == 1:
                val[abs(last)] = 1 if last > 0 else -1
                changed = True
    for cl in clauses:
        if not any(val[abs(l)] != 0 and (val[abs(l)] > 0) == (l > 0) for l in cl):
            return "OPEN", val
    return "SAT", val


def php_cnf(p, h):
    """Pigeonhole: p pigeons into h holes, at most one per hole (UNSAT when p > h)."""
    cnf = CNF(p * h)

    def v(i, j):
        return i * h + j + 1

    for i in range(p):
        cnf.add([v(i, j) for j in range(h)])
    for j in range(h):
        for i1 in range(p):
            for i2 in range(i1 + 1, p):
                cnf.add([-v(i1, j), -v(i2, j)])
    return cnf


# ---------------------------------------------------------------- specs and brute force

def slices_by_last(K):
    sl = {0: [], 1: [], 2: []}
    for p in K:
        sl[p[3]].append(kp.IDX3[tuple(p[:3])])
    return {s: sorted(v) for s, v in sl.items()}


def fixed_spec(K, N, zero_only, forbid=()):
    sl = slices_by_last(K)
    spec = kp.make_spec(N, sl[0], len(sl[1]), len(sl[2]), forbid=list(forbid),
                        fix={1: set(sl[1]), 2: set(sl[2])})
    if zero_only:
        spec["table"] = k26.zero_table(forbid)
        spec["cap"] = 1
    else:
        spec["cap"] = None
    return spec


def brute_ok(K, N, zero_only, forbid=()):
    """Does K satisfy exactly what the CNF of a fixed-slice spec demands?"""
    if len(set(K)) != N or len(K) != N or not kp.is_kakeya(K):
        return False
    allowed = set(k26.HAND_ZERO26 if zero_only else kp.allowed_profiles(N))
    allowed -= set(tuple(sorted(f)) for f in forbid)
    if any(tuple(sorted(kp.class_profile(K, u))) not in allowed for u in kp.DIR4):
        return False
    if zero_only and k26.lines_per_direction(K) != [1] * len(kp.DIR4):
        return False
    return True


def solution_points(spec, model):
    K1 = [i for i in range(27) if model.get(kp.X(1, i) + 1, False)]
    K2 = [i for i in range(27) if model.get(kp.X(2, i) + 1, False)]
    return kp.assemble(spec["R"], K1, K2), [1 if model.get(j + 1, False) else 0 for j in range(54)]


# ---------------------------------------------------------------- external tools

def write_cnf(cnf, path):
    text = cnf.dimacs()
    with open(path, "w") as fh:
        fh.write(text)
    return hashlib.sha256(text.encode()).hexdigest()


def run_cadical(cadical, cnf_path, proof_path, seconds):
    cmd = [cadical, "-q", "-t", str(max(1, int(seconds))), cnf_path]
    if proof_path:
        cmd.append(proof_path)
    t0 = time.time()
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                           timeout=seconds + 180)
    except subprocess.TimeoutExpired:
        return "UNKNOWN", None, time.time() - t0, "killed after timeout"
    dt = time.time() - t0
    if p.returncode == 10:
        model = {}
        for line in p.stdout.splitlines():
            if line.startswith("v "):
                for tok in line[2:].split():
                    v = int(tok)
                    if v != 0:
                        model[abs(v)] = v > 0
        return "SAT", model, dt, ""
    if p.returncode == 20:
        return "UNSAT", None, dt, ""
    if p.returncode == 0:
        return "UNKNOWN", None, dt, ""
    return "ERROR", None, dt, ("exit %d: " % p.returncode) + (p.stderr or p.stdout)[-400:]


def run_drat(drat, cnf_path, proof_path, seconds):
    t0 = time.time()
    try:
        p = subprocess.run([drat, cnf_path, proof_path, "-t", str(max(10, int(seconds)))],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                           timeout=seconds + 180)
    except subprocess.TimeoutExpired:
        return "TIMEOUT", time.time() - t0, ""
    dt = time.time() - t0
    out = p.stdout + "\n" + p.stderr
    s_lines = [l.strip() for l in out.splitlines() if l.startswith("s ")]
    if any(l in ("s VERIFIED", "s TRIVIAL UNSAT") for l in s_lines):
        return "VERIFIED", dt, " | ".join(s_lines)
    if any("NOT VERIFIED" in l for l in s_lines):
        return "NOT-VERIFIED", dt, " | ".join(s_lines)
    return "NO-ANSWER", dt, out[-400:]


def solve_spec(spec, tag, args, seconds, with_proof):
    """Write the CNF, run CaDiCaL (and drat-trim on UNSAT). Returns a result dict."""
    cnf = build_cnf(spec)
    res = {"vars": cnf.n, "clauses": len(cnf.clauses)}
    if cnf.trivial_unsat:
        res.update(status="UNSAT", proof="TRIVIAL-EMPTY-CLAUSE", seconds=0.0)
        return res, None
    os.makedirs(args.tmp, exist_ok=True)
    cpath = os.path.join(args.tmp, tag + ".cnf")
    ppath = os.path.join(args.tmp, tag + ".drat") if with_proof else None
    res["cnf_sha256"] = write_cnf(cnf, cpath)
    st, model, dt, err = run_cadical(args.cadical, cpath, ppath, seconds)
    res.update(status=st, seconds=round(dt, 2))
    if err:
        res["error"] = err
    if st == "UNSAT" and with_proof:
        size = os.path.getsize(ppath) if os.path.exists(ppath) else -1
        res["proof_bytes"] = size
        if size < 0:
            res["proof"] = "MISSING"
        elif size > args.max_proof_mb * 1e6:
            res["proof"] = "TOO-BIG"
        else:
            budget = min(max(600.0, 4 * dt), max(60.0, args.left() - 60))
            pst, pdt, msg = run_drat(args.drat, cpath, ppath, budget)
            res.update(proof=pst, check_seconds=round(pdt, 2), check_msg=msg[-300:])
    for path in (cpath, ppath):
        if path and os.path.exists(path):
            os.remove(path)
    return res, model


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="kdrat26.json")
    ap.add_argument("--cadical", default="cadical")
    ap.add_argument("--drat", default="drat-trim")
    ap.add_argument("--tmp", default="/tmp/kdrat26")
    ap.add_argument("--only", default="Z1,Z2,P1,P2,P3,P4")
    ap.add_argument("--t1", type=float, default=300, help="pass 1 CaDiCaL seconds per rep (proof on)")
    ap.add_argument("--max-proof-mb", type=float, default=300)
    ap.add_argument("--t-free", type=float, default=600)
    ap.add_argument("--deadline-min", type=float, default=340)
    ap.add_argument("--p1-frac", type=float, default=0.6)
    args = ap.parse_args()
    t_start = time.time()
    args.left = lambda: 60 * args.deadline_min - (time.time() - t_start)
    res = {"job": "kdrat26",
           "question": "independent SAT + DRAT re-check of kprof26: is there a 26-point Kakeya set "
                       "in F_3^4? (cases P1-P4, Z1, Z2; Z3-Z5 impossible by arithmetic)",
           "made_with": "AI agents working in MV2's name; no human has checked this", "bugs": []}
    bugs = res["bugs"]

    def log(msg):
        print(f"[{time.time() - t_start:8.1f}s] {msg}", flush=True)

    def dump():
        res["seconds"] = round(time.time() - t_start, 1)
        with open(args.out, "w") as fh:
            json.dump(res, fh, indent=1)

    # tool versions
    for name, cmd in (("cadical", [args.cadical, "--version"]),):
        try:
            res[name + "_version"] = subprocess.run(cmd, stdout=subprocess.PIPE, text=True,
                                                    timeout=60).stdout.strip()[:80]
        except Exception as exc:
            bugs.append(f"{name} not runnable: {exc!r}")
    if k26.pos_cuts() != k26.HAND_POS26 or k26.zero_cuts() != k26.HAND_ZERO26:
        bugs.append("kprof26 cut lists differ from the hand lists")
    dump()

    # sanity 1: the known 27-point set with fixed slices (SAT), and with its own last cut
    # forbidden (UNSAT, proof checked)
    K = [tuple(p) for p in kp.KNOWN27]
    s1 = {}
    try:
        spec = fixed_spec(K, 27, False)
        r, model = solve_spec(spec, "known27", args, 120, False)
        pts = sorted(solution_points(spec, model)[0]) if model else None
        s1["fixed"] = r["status"]
        s1["fixed_same_set"] = pts == sorted(K)
        cut = tuple(sorted(kp.class_profile(K, kp.LAST)))
        spec = fixed_spec(K, 27, False, forbid=[cut])
        r, _ = solve_spec(spec, "known27_forbid", args, 120, True)
        s1["forbid_own_cut"] = r["status"]
        s1["forbid_proof"] = r.get("proof")
        s1["forbid_msg"] = r.get("check_msg")
        # proof pipeline on an instance that needs real search: pigeonhole 8 -> 7 (UNSAT)
        os.makedirs(args.tmp, exist_ok=True)
        cpath, ppath = os.path.join(args.tmp, "php.cnf"), os.path.join(args.tmp, "php.drat")
        write_cnf(php_cnf(8, 7), cpath)
        st, _, dt, err = run_cadical(args.cadical, cpath, ppath, 300)
        pst, pdt, msg = run_drat(args.drat, cpath, ppath, 600) if st == "UNSAT" else ("NOT-RUN", 0, err)
        s1["php"] = {"status": st, "seconds": round(dt, 2), "proof": pst, "msg": msg[-200:]}
        s1["ok"] = (s1["fixed"] == "SAT" and s1["fixed_same_set"] and s1["forbid_own_cut"] == "UNSAT"
                    and st == "UNSAT" and pst == "VERIFIED")
    except Exception as exc:
        s1 = {"ok": False, "error": repr(exc)}
    res["known27"] = s1
    if not s1.get("ok"):
        bugs.append("sanity known27 failed")
    log(f"known27: {s1}")
    dump()

    # sanity 2: known set minus one point, fixed slices, both models, vs brute force
    rows = []
    try:
        for r_i in range(27):
            Kr = [p for i, p in enumerate(K) if i != r_i]
            row = {"removed": list(K[r_i])}
            for label, zero_only in (("all", False), ("zero", True)):
                want = "SAT" if brute_ok(Kr, 26, zero_only) else "UNSAT"
                spec = fixed_spec(Kr, 26, zero_only)
                r, _ = solve_spec(spec, f"rm{r_i}_{label}", args, 60, False)
                row[label] = r["status"]
                row[label + "_agree"] = r["status"] == want
            rows.append(row)
        s2 = {"ok": all(r["all_agree"] and r["zero_agree"] for r in rows),
              "any_sat": any(r["all"] == "SAT" for r in rows), "rows": rows}
    except Exception as exc:
        s2 = {"ok": False, "error": repr(exc), "rows": rows}
    res["removals"] = s2
    if not s2.get("ok"):
        bugs.append("sanity removals failed")
    log(f"removals ok={s2.get('ok')} any_sat={s2.get('any_sat')}")
    dump()

    # orbit representatives (kprof25), checked
    only = [s for s in args.only.split(",") if s]
    sel = [c for c in k26.CASES26 if c[0] in only]
    sel.sort(key=lambda c: only.index(c[0]))
    perms = kp.agl3_perms()
    kmax = max(c[1][0] for c in sel)
    reps = kp.orbit_reps(kmax, perms)
    oc = kp.orbit_check(reps, perms)
    res["orbits"] = {str(k): v for k, v in oc.items()}
    for k, v in oc.items():
        if not v["ok"]:
            bugs.append(f"orbit check failed for k={k}")

    # sanity 3: free search for 27 points with the known set's slice 0 (should be SAT)
    try:
        sl = kp.normalise_class(K, kp.DIR4[0], perms)
        spec = kp.make_spec(27, sl[0], len(sl[1]), len(sl[2]))
        spec["cap"] = None
        r, model = solve_spec(spec, "free27", args, args.t_free, False)
        s3 = {"status": r["status"], "seconds": r["seconds"]}
        if model:
            _, sol = solution_points(spec, model)
            ok, info = kp.check_solution(spec, sol)
            s3["check"] = ok
            if not ok:
                bugs.append(f"free27 solution fails the brute-force check: {info}")
    except Exception as exc:
        s3 = {"status": "ERROR", "error": repr(exc)}
        bugs.append("sanity free27 crashed")
    res["free27"] = s3
    log(f"free27: {s3}")
    dump()

    cases = []
    for name, cut, forbid, zero_only in sel:
        for code in reps[cut[0]]:
            cases.append({"case": name, "cut": list(cut), "forbid": [list(f) for f in forbid],
                          "zero_only": zero_only, "rep": kp.decode(code, cut[0])})
    res["cases"] = cases
    found = []
    res["found"] = found

    def spec_of(case):
        sp = k26.spec_for(tuple(case["cut"]), [tuple(f) for f in case["forbid"]],
                          case["zero_only"], case["rep"])
        return sp

    def record(case, r, model):
        if model is None:
            return
        spec = spec_of(case)
        pts, sol = solution_points(spec, model)
        ok, info = kp.check_solution(spec, sol)
        Kn = [tuple(p) for p in pts]
        ok = ok and k26.brute_case_ok(Kn, tuple(case["cut"]), [tuple(f) for f in case["forbid"]],
                                      case["zero_only"])
        case["sat_check"] = ok
        if ok:
            found.append({"case": case["case"], "rep": case["rep"], "points": info["points"]})
        else:
            bugs.append(f"CaDiCaL model for {case['case']} rep {case['rep']} fails the check")

    # pass 1: CaDiCaL with DRAT proof, then drat-trim
    p1_end = 60 * args.deadline_min * args.p1_frac
    for k, case in enumerate(cases):
        left_p1 = p1_end - (time.time() - t_start)
        if left_p1 < 60:
            break
        r, model = solve_spec(spec_of(case), f"c{k}", args, min(args.t1, left_p1), True)
        case["pass1"] = r
        record(case, r, model)
        dump()
        log(f"pass 1 {case['case']} rep {case['rep']}: {r.get('status')} {r.get('seconds')}s "
            f"proof {r.get('proof')} {r.get('proof_bytes')}")
    # pass 2: undecided reps without proof, the rest of the time shared
    open_ = [c for c in cases if c.get("pass1", {}).get("status") not in ("UNSAT", "SAT")]
    for k, case in enumerate(open_):
        if args.left() < 300:
            break
        share = max(120.0, (args.left() - 240) / (len(open_) - k))
        r, model = solve_spec(spec_of(case), f"o{k}", args, share, False)
        case["pass2"] = r
        record(case, r, model)
        dump()
        log(f"pass 2 {case['case']} rep {case['rep']}: {r.get('status')} {r.get('seconds')}s")

    def unsat(c):
        return (c.get("pass1", {}).get("status") == "UNSAT"
                or c.get("pass2", {}).get("status") == "UNSAT")

    def certified(c):
        p1 = c.get("pass1", {})
        return p1.get("status") == "UNSAT" and p1.get("proof") in ("VERIFIED", "TRIVIAL-EMPTY-CLAUSE")

    per = {}
    for name, cut, _, _ in sel:
        cb = [c for c in cases if c["case"] == name]
        per[name] = {"cut": list(cut), "reps": len(cb),
                     "unsat": sum(1 for c in cb if unsat(c)),
                     "certified": sum(1 for c in cb if certified(c)),
                     "open": [c["rep"] for c in cb if not unsat(c)]}
    res["per_case"] = per
    all_unsat = bool(cases) and all(unsat(c) for c in cases)
    all_cert = all_unsat and all(certified(c) for c in cases)
    whole = set(only) >= {"P1", "P2", "P3", "P4", "Z1", "Z2"}
    if bugs:
        verdict = "BUG"
    elif found:
        verdict = "FOUND-26"
    elif all_cert:
        verdict = "CERT-K4-EQ-27" if whole else "CERT-PART-CLOSED"
    elif all_unsat:
        verdict = "UNSAT-PROOF-INCOMPLETE"
    else:
        verdict = "UNDECIDED"
    res["verdict"] = verdict
    dump()
    log(f"verdict {verdict}; per case {json.dumps(per)}; bugs {bugs}")


if __name__ == "__main__":
    main()
