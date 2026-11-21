# kmin34cert.py -- a certificate-grade re-check that no Kakeya set in F_4^3 has 27 points,
# i.e. K(3, 4) >= 28 (and so K(3, 4) = 28, since 28-point Kakeya sets exist).
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0.
#
# Why: kmin34-v1 (kmin34/kmin34.py) reduced the question to 22 small 0/1 models and found all
# 22 infeasible with the MIP solver HiGHS, but CP-SAT proved only 16 of them within 900 s.
# This job re-checks the same reduction with a different method that leaves a checkable
# certificate: every case is written as a CNF formula by NEW code (nothing is imported from
# kmin34.py), solved by the SAT solver kissat with a DRAT proof, and every proof is checked
# by drat-trim. The orbit representatives are copied from kmin34-v1, but their coverage is
# re-verified here with an independently built group.
#
# The reduction (same as kmin34/prereg-kmin34.md, restated):
#  Points of F_4^3 are (y, s), y in F_4^2, s in F_4; slices K_s = {y : (y, s) in K}.
#  * Direction (w, 1), w in F_4^2 (16 directions): the line through (a, 0) is
#    {(a + t w, t) : t in F_4}, and every line in this direction meets s = 0 once, so K has
#    one iff some a in K_0 has a + t w in K_t for t = 1, 2, 3.
#  * Direction (d, 0), d in PG(1, 4) (5 directions): lines lie inside one slice, so K has one
#    iff some K_s contains a line of F_4^2 in direction d.
#  * If K is Kakeya with |K| <= 27, take a plane H with |K n H| = n0 minimal. A Kakeya set meets
#    every plane (a line not parallel to H meets it), so n0 >= 1, and four parallel planes
#    partition F_4^3, so n0 <= 27 // 4 = 6. An affine map of F_4^3 sends H to {s = 0}; then
#    (y, s) -> (g y, s) with g in AGL(2, 4) is affine on F_4^3 and sends K_0 to any chosen
#    member R of its AGL(2, 4)-orbit. The image is a Kakeya set with slice 0 = R, all 84 plane
#    sections >= |R|, and at most 27 points.
#  * So: if for every listed representative R (whose orbits cover all k-subsets, k = 1..6)
#    the formula CNF(R, 27) is unsatisfiable, no Kakeya set of <= 27 points exists.
#
# Encoding (48 point variables x(s, i), s = 1, 2, 3, i = 4 y0 + y1):
#  * transversal: for each w, a clause OR_a y(w, a), with y(w, a) -> x(t, a + t w), t = 1, 2, 3;
#  * horizontal: for each d not already covered by a line inside R, a clause OR z(s, L) over
#    s = 1, 2, 3 and the 4 lines L in direction d, with z(s, L) -> x(s, i) for i in L;
#  * every plane section >= |R| and the total <= 27 - |R|: Sinz sequential counters.
#  The auxiliary variables only appear negatively except in the OR clauses, so a 0/1 point
#  assignment extends to a model iff it satisfies the conditions; the self-checks below test
#  this against brute force on many sets.

import argparse
import itertools
import json
import math
import os
import random
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor

Q = 4
# F_4 = {0, 1, 2, 3}, 2 = w, 3 = w + 1 = w^2 (the encoding of kakeya.py and kmin34.py).
MUL = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]


def add(p, q):
    return tuple(a ^ b for a, b in zip(p, q))


def smul(c, p):
    return tuple(MUL[c][a] for a in p)


def dot(u, v):
    s = 0
    for a, b in zip(u, v):
        s ^= MUL[a][b]
    return s


PTS2 = [(a, b) for a in range(Q) for b in range(Q)]  # index 4 a + b (as in kmin34.py)
IDX2 = {p: i for i, p in enumerate(PTS2)}
PTS3 = [(a, b, c) for a in range(Q) for b in range(Q) for c in range(Q)]


def proj(n):
    """Directions: non-zero vectors of F_4^n whose first non-zero coordinate is 1."""
    out = []
    for v in itertools.product(range(Q), repeat=n):
        nz = [x for x in v if x]
        if nz and nz[0] == 1:
            out.append(v)
    return out


DIR2 = proj(2)  # 5
DIR3 = proj(3)  # 21


def lines2(d):
    """The 4 lines of F_4^2 in direction d, as sorted tuples of point indices."""
    return sorted({tuple(sorted(IDX2[add(p, smul(t, d))] for t in range(Q))) for p in PTS2})


LINES2 = {d: lines2(d) for d in DIR2}

# Orbit representatives of k-subsets of F_4^2 under AGL(2, 4), k = 1..6, copied from the
# kmin34-v1 output (results/kmin34-v1/kmin34.json); coverage is re-checked by coverage().
REPS = {
    1: [[0]],
    2: [[0, 1]],
    3: [[0, 1, 2], [0, 1, 4]],
    4: [[0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 4, 5], [0, 1, 4, 6]],
    5: [[0, 1, 2, 3, 4], [0, 1, 2, 4, 5], [0, 1, 2, 4, 8], [0, 1, 2, 4, 10], [0, 1, 4, 6, 9]],
    6: [[0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 8], [0, 1, 2, 4, 5, 6], [0, 1, 2, 4, 5, 8],
        [0, 1, 2, 4, 5, 9], [0, 1, 2, 4, 5, 13], [0, 1, 2, 4, 8, 15], [0, 1, 2, 4, 10, 13],
        [0, 1, 4, 6, 9, 10]],
}

# The 28-point Kakeya set of kakeya-min-v1 (also embedded in kmin34.py).
KNOWN28 = [(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0), (2, 2, 0), (3, 2, 0),
           (0, 0, 1), (3, 0, 1), (0, 1, 1), (1, 1, 1), (0, 2, 1), (2, 2, 1), (0, 3, 1),
           (1, 3, 1), (2, 3, 1), (3, 3, 1),
           (0, 1, 2), (1, 1, 2), (2, 1, 2), (3, 1, 2), (1, 2, 2), (2, 2, 2),
           (0, 2, 3), (1, 2, 3), (2, 2, 3), (3, 2, 3), (1, 3, 3), (3, 3, 3)]


# ---------------------------------------------------------------- brute force

def field_ok():
    R = range(Q)
    for a in R:
        for b in R:
            if MUL[a][b] != MUL[b][a]:
                return False
            for c in R:
                if MUL[a][MUL[b][c]] != MUL[MUL[a][b]][c]:
                    return False
                if MUL[a][b ^ c] != MUL[a][b] ^ MUL[a][c]:
                    return False
        if a and sorted(MUL[a]) != [0, 1, 2, 3]:
            return False
    return MUL[1] == [0, 1, 2, 3] and MUL[2][2] == 3


def is_kakeya3(K):
    S = set(tuple(x) for x in K)
    for v in DIR3:
        if not any(all(add(x, smul(t, v)) in S for t in (1, 2, 3)) for x in S):
            return False
    return True


def min_section(K):
    return min(sum(1 for x in K if dot(n, x) == c) for n in DIR3 for c in range(Q))


def from_slices(sl):
    return [PTS2[i] + (s,) for s in range(Q) for i in sl[s]]


# ---------------------------------------------------------------- AGL(2, 4)

def agl2():
    """All 2,880 maps y -> M y + t as tuples: g[i] = index of the image of point i."""
    mats = [((a, b), (c, d)) for a, b, c, d in itertools.product(range(Q), repeat=4)
            if MUL[a][d] ^ MUL[b][c]]
    perms = []
    for M in mats:
        for t in PTS2:
            perms.append(tuple(IDX2[add((dot(M[0], p), dot(M[1], p)), t)] for p in PTS2))
    return perms


def group_ok(perms):
    lines = {L for d in DIR2 for L in LINES2[d]}
    ok = len(perms) == 2880 and len(set(perms)) == 2880 and len(lines) == 20
    ok = ok and all(sorted(g) == list(range(16)) for g in perms)
    for g in perms:
        for L in lines:
            if tuple(sorted(g[i] for i in L)) not in lines:
                return False
    return ok


def coverage(perms):
    out = {}
    for k in range(1, 7):
        seen = set()
        for R in REPS[k]:
            for g in perms:
                seen.add(frozenset(g[i] for i in R))
        out[k] = {"reps": len(REPS[k]), "covered": len(seen), "binom": math.comb(16, k),
                  "ok": len(seen) == math.comb(16, k)}
    return out


def normalise(K, perms):
    """Affine image of K with a smallest plane section as slice 0, moved onto a listed
    representative. Returns (n0, {s: sorted indices}) or (n0, None) if n0 is not 1..6."""
    best = None
    for nrm in DIR3:
        for c in range(Q):
            m = sum(1 for x in K if dot(nrm, x) == c)
            if best is None or m < best[0]:
                best = (m, nrm, c)
    m, nrm, c = best
    if m not in REPS:
        return m, None
    j = next(i for i in range(3) if nrm[i])
    rows = [tuple(1 if r == i else 0 for r in range(3)) for i in range(3) if i != j] + [nrm]
    K2 = [(dot(rows[0], x), dot(rows[1], x), dot(rows[2], x) ^ c) for x in K]
    assert len(set(K2)) == len(K)
    sl = {s: [IDX2[(x[0], x[1])] for x in K2 if x[2] == s] for s in range(Q)}
    targets = {frozenset(R) for R in REPS[m]}
    for g in perms:
        if frozenset(g[i] for i in sl[0]) in targets:
            return m, {s: sorted(g[i] for i in sl[s]) for s in range(Q)}
    raise RuntimeError("no representative found (coverage would have failed)")


# ---------------------------------------------------------------- CNF

class CNF:
    def __init__(self, nfixed=48):
        self.n = nfixed
        self.clauses = []
        self.ands = []  # (aux var, point vars it implies): used only by canonical()
        self.ctrs = []  # (literals, counter matrix s): used only by canonical()

    def new(self):
        self.n += 1
        return self.n

    def add(self, cl):
        self.clauses.append(list(cl))

    def false(self):
        v = self.new()
        self.add([v])
        self.add([-v])

    def atmost(self, lits, k):
        """Sinz sequential counter: at most k of lits are true."""
        n = len(lits)
        if k >= n:
            return
        if k < 0:
            self.false()
            return
        if k == 0:
            for l in lits:
                self.add([-l])
            return
        s = [[self.new() for _ in range(k)] for _ in range(n - 1)]
        self.ctrs.append((list(lits), s))
        self.add([-lits[0], s[0][0]])
        for j in range(1, k):
            self.add([-s[0][j]])
        for i in range(1, n - 1):
            self.add([-lits[i], s[i][0]])
            self.add([-s[i - 1][0], s[i][0]])
            for j in range(1, k):
                self.add([-lits[i], -s[i - 1][j - 1], s[i][j]])
                self.add([-s[i - 1][j], s[i][j]])
            self.add([-lits[i], -s[i - 1][k - 1]])
        self.add([-lits[n - 1], -s[n - 2][k - 1]])

    def atleast(self, lits, m):
        self.atmost([-l for l in lits], len(lits) - m)

    def dimacs(self):
        out = [f"p cnf {self.n} {len(self.clauses)}"]
        out += [" ".join(str(x) for x in c) + " 0" for c in self.clauses]
        return "\n".join(out) + "\n"


def X(s, i):
    return (s - 1) * 16 + i + 1


def build_cnf(R, bound, fix=None):
    R = sorted(R)
    assert R, "slice 0 is never empty"
    Rs = set(R)
    n0 = len(R)
    f = CNF()
    for w in PTS2:  # transversal directions (w, 1)
        ys = []
        for a in R:
            y = f.new()
            ys.append(y)
            needs = [X(t, IDX2[add(PTS2[a], smul(t, w))]) for t in (1, 2, 3)]
            f.ands.append((y, needs))
            for x in needs:
                f.add([-y, x])
        f.add(ys)
    for d in DIR2:  # horizontal directions (d, 0)
        if any(set(L) <= Rs for L in LINES2[d]):
            continue
        zs = []
        for s in (1, 2, 3):
            for L in LINES2[d]:
                z = f.new()
                zs.append(z)
                f.ands.append((z, [X(s, i) for i in L]))
                for i in L:
                    f.add([-z, X(s, i)])
        f.add(zs)
    for nrm in DIR3:  # every plane section >= n0
        for c in range(Q):
            cnt0 = sum(1 for a in R if dot(nrm, PTS2[a] + (0,)) == c)
            lits = [X(s, i) for s in (1, 2, 3) for i in range(16)
                    if dot(nrm, PTS2[i] + (s,)) == c]
            if n0 - cnt0 > 0:
                f.atleast(lits, n0 - cnt0)
    f.atmost([X(s, i) for s in (1, 2, 3) for i in range(16)], bound - n0)
    if fix is not None:
        for s in (1, 2, 3):
            for i in range(16):
                f.add([X(s, i)] if i in fix[s] else [-X(s, i)])
    return f


def slices_from_model(R, true_vars):
    sl = {0: sorted(R)}
    for s in (1, 2, 3):
        sl[s] = [i for i in range(16) if X(s, i) in true_vars]
    return sl


def canonical_ok(f, base):
    """Solver-free completeness check. base: {var: bool} for the variables 1..f's fixed ones.
    Sets every AND-auxiliary to the AND of its point variables and every counter variable
    s[i][j] to 'at least j + 1 of the first i + 1 literals are true', then reports whether
    every clause holds. If the point values satisfy the conditions, this assignment is a model."""
    val = [False] * (f.n + 1)
    for v, b in base.items():
        val[v] = bool(b)

    def lit(l):
        return val[l] if l > 0 else not val[-l]

    for y, needs in f.ands:
        val[y] = all(val[x] for x in needs)
    for lits, s in f.ctrs:
        cnt = 0
        for i in range(len(lits) - 1):
            cnt += 1 if lit(lits[i]) else 0
            for j in range(len(s[i])):
                val[s[i][j]] = cnt >= j + 1
    return all(any(lit(l) for l in cl) for cl in f.clauses)


def up_conflict(f):
    """Solver-free soundness check: unit propagation from the formula's own unit clauses.
    With every point variable fixed by units, an assignment that breaks a condition gives a
    conflict here (all ANDs of an OR clause are forced false, or a counter overflows)."""
    val = {}
    changed = True
    while changed:
        changed = False
        for cl in f.clauses:
            free, sat = [], False
            for l in cl:
                v = val.get(abs(l))
                if v is None:
                    free.append(l)
                elif v == (l > 0):
                    sat = True
                    break
            if sat:
                continue
            if not free:
                return True
            if len(free) == 1:
                val[abs(free[0])] = free[0] > 0
                changed = True
    return False


# ---------------------------------------------------------------- external tools

def run_kissat(kissat, cnf_path, proof_path, timeout):
    """Returns (status, set of true variables or None, seconds, return code)."""
    if proof_path and os.path.exists(proof_path):
        os.remove(proof_path)
    cmd = [kissat, cnf_path] + ([proof_path] if proof_path else [])
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return "UNKNOWN", None, time.time() - t0, None
    dt = time.time() - t0
    status, vals = "UNKNOWN", []
    for line in p.stdout.splitlines():
        if line.startswith("s "):
            word = line[2:].strip()
            if word == "SATISFIABLE":
                status = "SAT"
            elif word == "UNSATISFIABLE":
                status = "UNSAT"
        elif line.startswith("v "):
            vals += [int(x) for x in line[2:].split()]
    if (status == "SAT" and p.returncode != 10) or (status == "UNSAT" and p.returncode != 20):
        status = "ERROR"
    if status == "UNKNOWN" and p.returncode not in (0, 10, 20):
        status = "ERROR"
    true_vars = {v for v in vals if v > 0} if status == "SAT" else None
    return status, true_vars, dt, p.returncode


VERIFIED_LINES = ("s VERIFIED", "s TRIVIAL UNSAT")


def run_drat(drat, cnf_path, proof_path, timeout):
    """Returns (verdict, seconds, last 's' line). verdict is VERIFIED only if drat-trim says so."""
    t0 = time.time()
    try:
        p = subprocess.run([drat, cnf_path, proof_path, "-t", str(int(timeout))],
                           capture_output=True, text=True, timeout=timeout + 120)
    except subprocess.TimeoutExpired:
        return "TIMEOUT", time.time() - t0, ""
    s_lines = [l.strip() for l in p.stdout.splitlines() if l.startswith("s ")]
    last = s_lines[-1] if s_lines else ""
    ok = any(l in VERIFIED_LINES for l in s_lines) and "s NOT VERIFIED" not in s_lines
    return ("VERIFIED" if ok else "NOT-VERIFIED"), time.time() - t0, last


def solve_cnf(kissat, f, path, timeout=120):
    with open(path, "w") as fh:
        fh.write(f.dimacs())
    return run_kissat(kissat, path, None, timeout)


def tool_info(cmd):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        return (p.stdout + p.stderr).strip()[:300]
    except Exception as exc:
        return f"unavailable: {exc!r}"


# ---------------------------------------------------------------- self-checks

def php_cnf(holes):
    """Pigeonhole principle, holes + 1 pigeons: unsatisfiable."""
    f = CNF(nfixed=0)
    var = {}
    for p in range(holes + 1):
        for h in range(holes):
            var[p, h] = f.new()
    for p in range(holes + 1):
        f.add([var[p, h] for h in range(holes)])
    for h in range(holes):
        for p1 in range(holes + 1):
            for p2 in range(p1 + 1, holes + 1):
                f.add([-var[p1, h], -var[p2, h]])
    return f


def self_checks(args, perms, log):
    chk, bugs = {}, []
    tmp = args.tmp
    chk["field_ok"] = field_ok()
    chk["group_ok"] = group_ok(perms)
    cov = coverage(perms)
    chk["coverage"] = {str(k): v for k, v in cov.items()}
    chk["known28_kakeya"] = is_kakeya3(KNOWN28) and len(set(KNOWN28)) == 28
    chk["full_space_kakeya"] = is_kakeya3(PTS3)
    rng = random.Random(2108)  # Bukh-Chao: a Kakeya set in F_4^3 has >= 21 points
    chk["random20_not_kakeya"] = not any(is_kakeya3(rng.sample(PTS3, 20)) for _ in range(200))
    for key in ("field_ok", "group_ok", "known28_kakeya", "full_space_kakeya",
                "random20_not_kakeya"):
        if not chk[key]:
            bugs.append(f"self-check {key} failed")
    for k, v in cov.items():
        if not v["ok"]:
            bugs.append(f"representatives of size {k} cover {v['covered']} of {v['binom']}")
    log(f"pure checks: {json.dumps({k: v for k, v in chk.items() if k != 'coverage'})}")

    # the SAT/DRAT pipeline on known formulas
    p = os.path.join(tmp, "php.cnf")
    with open(p, "w") as fh:
        fh.write(php_cnf(4).dimacs())
    st, _, _, _ = run_kissat(args.kissat, p, p + ".drat", 300)
    dv = run_drat(args.drat, p, p + ".drat", 600)[0] if st == "UNSAT" else "SKIPPED"
    chk["php4_unsat_verified"] = (st, dv)
    if st != "UNSAT" or dv != "VERIFIED":
        bugs.append(f"pigeonhole 5->4: kissat {st}, drat-trim {dv}")
    f = CNF(nfixed=2)
    f.add([1, 2])
    f.add([-1, 2])
    f.add([1, -2])
    st = solve_cnf(args.kissat, f, os.path.join(tmp, "sat.cnf"))[0]
    with open(os.path.join(tmp, "empty.drat"), "w") as fh:
        fh.write("")
    dv = run_drat(args.drat, os.path.join(tmp, "sat.cnf"), os.path.join(tmp, "empty.drat"), 60)[0]
    chk["sat_formula"] = (st, dv)
    if st != "SAT" or dv == "VERIFIED":
        bugs.append(f"satisfiable formula: kissat {st}, drat-trim on empty proof {dv}")

    # the counters against popcount, exhaustively for 5 literals
    bad = 0
    for kind in ("atmost", "atleast"):
        for k in range(-1, 7):
            for bits in range(32):
                f = CNF(nfixed=5)
                getattr(f, kind)([1, 2, 3, 4, 5], k)
                for v in range(5):
                    f.add([v + 1] if bits >> v & 1 else [-(v + 1)])
                pc = bin(bits).count("1")
                want = pc <= k if kind == "atmost" else pc >= k
                st = solve_cnf(args.kissat, f, os.path.join(tmp, "ctr.cnf"))[0]
                if st != ("SAT" if want else "UNSAT"):
                    bad += 1
    chk["counter_mismatches"] = bad
    if bad:
        bugs.append(f"{bad} counter mismatches against popcount")
    log(f"pipeline/counter checks: php {chk['php4_unsat_verified']}, sat {chk['sat_formula']}, "
        f"counter mismatches {bad}")

    # the known 28-point set: normalised, accepted at 28 and refused at 27 (slices fixed)
    n0, sl = normalise(KNOWN28, perms)
    Kn = from_slices(sl)
    k28 = {"n0": n0, "rep": sl[0], "slice_sizes": [len(sl[s]) for s in range(Q)],
           "normalised_kakeya": is_kakeya3(Kn), "normalised_min_section": min_section(Kn)}
    fix = {s: set(sl[s]) for s in (1, 2, 3)}
    k28["fixed_28"] = solve_cnf(args.kissat, build_cnf(sl[0], 28, fix),
                                os.path.join(tmp, "k28.cnf"))[0]
    p = os.path.join(tmp, "k27.cnf")
    with open(p, "w") as fh:
        fh.write(build_cnf(sl[0], 27, fix).dimacs())
    st = run_kissat(args.kissat, p, p + ".drat", 300)[0]
    k28["fixed_27"] = st
    k28["fixed_27_drat"] = run_drat(args.drat, p, p + ".drat", 600)[0] if st == "UNSAT" else ""
    st, tv, _, _ = solve_cnf(args.kissat, build_cnf(sl[0], 28), os.path.join(tmp, "f28.cnf"), 600)
    k28["free_28"] = st
    if tv is not None:
        K = from_slices(slices_from_model(sl[0], tv))
        k28["free_28_witness_ok"] = (is_kakeya3(K) and len(set(K)) == len(K) <= 28
                                     and min_section(K) == n0)
    k28["ok"] = (k28["normalised_kakeya"] and k28["normalised_min_section"] == n0
                 and k28["fixed_28"] == "SAT" and k28["fixed_27"] == "UNSAT"
                 and k28["fixed_27_drat"] == "VERIFIED" and k28["free_28"] == "SAT"
                 and k28.get("free_28_witness_ok", False))
    chk["known28"] = k28
    if not k28["ok"]:
        bugs.append(f"known 28-point set check failed: {k28}")
    log(f"known28: {k28}")

    # CNF + fixed points against brute force on perturbed and random sets
    rng = random.Random(34)
    tally = {"agree_sat": 0, "agree_unsat": 0, "disagree": 0, "skipped": 0}
    examples = []
    for trial in range(args.trials):
        K = set(KNOWN28)
        if trial % 3 == 2:
            p_keep = rng.choice([0.6, 0.8, 0.9])
            K = {x for x in PTS3 if rng.random() < p_keep}
        else:
            for _ in range(rng.randint(0, 3)):
                K.discard(rng.choice(sorted(K)))
            for _ in range(rng.randint(0, 6)):
                K.add(rng.choice(PTS3))
        K = sorted(K)
        n0, sl = normalise(K, perms)
        if sl is None:
            tally["skipped"] += 1
            continue
        b = rng.choice([len(K), len(K) - 1, 64])
        want = is_kakeya3(K) and len(K) <= b
        fix = {s: set(sl[s]) for s in (1, 2, 3)}
        f = build_cnf(sl[0], b, fix)
        st = solve_cnf(args.kissat, f, os.path.join(tmp, "rnd.cnf"))[0]
        base = {X(s, i): (i in fix[s]) for s in (1, 2, 3) for i in range(16)}
        c_ok = canonical_ok(f, base)
        u_cf = up_conflict(f)
        if st == ("SAT" if want else "UNSAT") and c_ok == want and u_cf == (not want):
            tally["agree_sat" if want else "agree_unsat"] += 1
        else:
            tally["disagree"] += 1
            examples.append({"K": [list(x) for x in K], "bound": b, "want": want, "got": st,
                             "canonical_ok": c_ok, "up_conflict": u_cf})
    chk["agreement"] = tally
    chk["disagreements"] = examples[:5]
    if tally["disagree"]:
        bugs.append(f"{tally['disagree']} disagreements between the CNF and brute force")
    if tally["agree_sat"] < 3 or tally["agree_unsat"] < 3:
        bugs.append(f"agreement test too one-sided: {tally}")
    log(f"agreement: {tally}")
    return chk, bugs


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="kmin34cert.json")
    ap.add_argument("--cnf-dir", default="")
    ap.add_argument("--tmp", default="/tmp/kmin34cert")
    ap.add_argument("--kissat", default=os.environ.get("KISSAT", "/opt/kissat/build/kissat"))
    ap.add_argument("--drat", default=os.environ.get("DRAT_TRIM", "/opt/drat-trim/drat-trim"))
    ap.add_argument("--bound", type=int, default=27)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--t-kissat", type=float, default=2400)
    ap.add_argument("--t-drat", type=float, default=3600)
    ap.add_argument("--trials", type=int, default=90)
    ap.add_argument("--deadline-min", type=float, default=320)
    ap.add_argument("--only", default="", help="comma list of case tags, e.g. 5.0,5.3")
    args = ap.parse_args()
    os.makedirs(args.tmp, exist_ok=True)
    if args.cnf_dir:
        os.makedirs(args.cnf_dir, exist_ok=True)
    t_start = time.time()
    lock = threading.Lock()
    res = {"job": "kmin34cert-v1",
           "question": "is there a Kakeya set of at most 27 points in F_4^3? (kmin34-v1: HiGHS "
                       "says no in all 22 cases, CP-SAT in 16)",
           "made_with": "AI agents working in MV2's name; no human has checked this",
           "bound": args.bound, "bugs": []}

    def log(msg):
        print(f"[{time.time() - t_start:8.1f}s] {msg}", flush=True)

    def dump():
        with lock:
            res["seconds"] = round(time.time() - t_start, 1)
            with open(args.out, "w") as fh:
                json.dump(res, fh, indent=1)

    res["tools"] = {
        "kissat_version": tool_info([args.kissat, "--version"]),
        "kissat_commit": tool_info(["git", "-C", os.path.dirname(os.path.dirname(args.kissat)),
                                    "rev-parse", "HEAD"]),
        "drat_trim_commit": tool_info(["git", "-C", os.path.dirname(args.drat),
                                       "rev-parse", "HEAD"]),
    }
    perms = agl2()
    chk, bugs = self_checks(args, perms, log)
    res["self_checks"] = chk
    res["bugs"] = list(bugs)
    dump()

    cases = []
    for n0 in range(1, args.bound // 4 + 1):
        for j, R in enumerate(REPS[n0]):
            tag = f"{n0}.{j}"
            if args.only and tag not in args.only.split(","):
                continue
            cases.append((tag, R))
    res["cases"] = {}

    def work(item):
        tag, R = item
        case = {"n0": len(R), "rep": R, "rep_points": [list(PTS2[i]) for i in R]}
        if time.time() - t_start > 60 * args.deadline_min:
            case["kissat"] = "NOT-RUN"
        else:
            f = build_cnf(R, args.bound)
            text = f.dimacs()
            case["vars"], case["clauses"] = f.n, len(f.clauses)
            path = os.path.join(args.tmp, f"case-{tag}.cnf")
            with open(path, "w") as fh:
                fh.write(text)
            if args.cnf_dir:
                with open(os.path.join(args.cnf_dir, f"case-{tag}.cnf"), "w") as fh:
                    fh.write(text)
            proof = path + ".drat"
            st, tv, dt, rc = run_kissat(args.kissat, path, proof, args.t_kissat)
            case.update({"kissat": st, "kissat_s": round(dt, 1), "kissat_rc": rc})
            if st == "UNSAT":
                case["proof_bytes"] = os.path.getsize(proof) if os.path.exists(proof) else -1
                dv, dt2, last = run_drat(args.drat, path, proof, args.t_drat)
                case.update({"drat": dv, "drat_s": round(dt2, 1), "drat_line": last})
            elif st == "SAT":
                sl = slices_from_model(R, tv)
                K = from_slices(sl)
                case["witness"] = [list(x) for x in K]
                case["witness_ok"] = (is_kakeya3(K) and len(set(K)) == len(K) <= args.bound
                                      and min_section(K) == len(R))
            if os.path.exists(proof):
                os.remove(proof)
        with lock:
            res["cases"][tag] = case
        log(f"case {tag} rep={R}: kissat {case.get('kissat')} {case.get('kissat_s', 0)}s, "
            f"drat {case.get('drat', '-')} {case.get('drat_s', 0)}s")
        dump()

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        list(ex.map(work, cases))

    cs = [res["cases"][t] for t, _ in cases]
    found = [c for c in cs if c.get("kissat") == "SAT" and c.get("witness_ok")]
    badwit = [c for c in cs if c.get("kissat") == "SAT" and not c.get("witness_ok")]
    if badwit:
        res["bugs"].append("a SAT witness fails the brute-force check")
    errors = [c for c in cs if c.get("kissat") == "ERROR"]
    if errors:
        res["bugs"].append(f"{len(errors)} kissat runs ended in an error")
    certified = [c for c in cs if c.get("kissat") == "UNSAT" and c.get("drat") == "VERIFIED"]
    res["certified"] = len(certified)
    res["cases_total"] = len(cs)
    res["open_cases"] = [t for t, _ in cases
                         if not (res["cases"][t].get("kissat") == "UNSAT"
                                 and res["cases"][t].get("drat") == "VERIFIED")]
    if res["bugs"]:
        verdict = "BUG"
    elif found:
        verdict = "FOUND-27"
    elif not args.only and len(certified) == len(cs) == 22:
        verdict = "CERT-EQ-28"
    else:
        verdict = "PARTIAL"
    res["verdict"] = verdict
    dump()
    log(f"verdict {verdict}; certified {len(certified)}/{len(cs)}; bugs {res['bugs']}")


if __name__ == "__main__":
    main()
