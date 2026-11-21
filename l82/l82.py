# l82.py -- is there a 25-point 1-saturating set in PG(7,2)?  Known: 25 <= l(8,2) <= 26.
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0.
#
# l(r,2) = shortest binary linear code of codimension r and covering radius 2 = smallest S in
# F_2^r \ {0} such that every nonzero vector is in S or is a sum of two distinct elements of S
# (S is a 1-saturating set of PG(r-1,2)). Also K(r,2) = l(r,2) + 1 for the Kakeya number over F_2.
# Lower bound l(8,2) >= 25: Struik (IEEE Trans. Inf. Theory 40, 1994). Upper bound 26.
#
# Reduction (l82/prereg-l82sat.md in munimtechnologies/kakeya-conjecture). Let S be 1-saturating
# in F_2^n with |S| = N.
#  L1. S spans F_2^n (otherwise S and all sums lie in a proper subspace).
#  L2. For a nonzero functional g, a_g = |{s in S : g(s) = 1}| satisfies a_g (N + 1 - a_g) >= 2^(n-1):
#      the 2^(n-1) vectors with g = 1 are in S_1 or are s_0 + s_1. So a_lo <= a_g <= N + 1 - a_lo.
#  L3. The a_g are the weights of the [N, n] code {(g(s))_s}, so min_g a_g <= its minimum
#      distance; also min_g a_g <= N 2^(n-1) / (2^n - 1) (average). For n = 8, N = 25:
#      a_lo = 7, and a binary [25,8] code has d <= 10 (Griesmer: d = 11 needs length >= 26) and in
#      fact d <= 9 (Ytrehus-Helleseth 1990: no [25,8,10] code). So min a_g is 7, 8, 9 (or 10).
#  Pick g = f of minimum weight a and coordinates with f = last bit. Then S = T u {(u,1) : u in U},
#  T subset of F_2^m \ {0} (m = n - 1), |T| = N - a, U subset of F_2^m, |U| = a, and
#   (C1) U + (T u {0}) = F_2^m                      (covers the vectors with last bit 1)
#   (C2) T u (T (+) T) u (U (+) U) contains F_2^m \ {0}  ((+): sums of two distinct elements)
#   (C3) every functional (h, e) != (0, 0) has weight |{t : h.t = 1}| + |{u : h.u + e = 1}| in
#        [a, N + 1 - a_lo].
#  The maps (y, s) -> (A y + s c, s) (A in GL(m,2)) fix f and send (T, U) to (A T, A U + c), so U
#  may be replaced by any set in its AGL(m,2)-orbit: U = {0, e_1..e_r} u Y with r = affine rank,
#  Y a set of vectors of weight >= 2 in <e_1..e_r>, taken up to permutations of e_1..e_r.
#  For each such U the existence of T is one CP-SAT model with 2^m - 1 point variables.
# Every set a solver returns is re-checked by brute force on the full 2^n - 1 points.

import argparse
import itertools
import json
import time

import numpy as np


def popcount(x):
    return bin(x).count("1")


def parity(x):
    return popcount(x) & 1


# ---------------------------------------------------------------- brute-force side

def is_saturating(S, n):
    """Every nonzero vector of F_2^n is in S or a sum of two distinct elements of S."""
    S = sorted(set(S))
    if any(s <= 0 or s >= (1 << n) for s in S):
        return False
    cov = set(S)
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            cov.add(S[i] ^ S[j])
    return len(cov) == (1 << n) - 1


def dual_weights(S, n):
    """{g: |{s in S : g.s = 1}|} for every nonzero functional g of F_2^n."""
    return {g: sum(parity(g & s) for s in S) for g in range(1, 1 << n)}


def split(S, n):
    """(T, U): points with last bit 0, and the low parts of points with last bit 1."""
    m = n - 1
    T = sorted(s for s in S if not (s >> m) & 1)
    U = sorted(s & ((1 << m) - 1) for s in S if (s >> m) & 1)
    return T, U


def join(T, U, n):
    m = n - 1
    return sorted(list(T) + [u | (1 << m) for u in U])


def reduced_ok(T, U, m):
    """(C1) and (C2) checked directly."""
    M = 1 << m
    T0 = set(T) | {0}
    c1 = {u ^ t for u in U for t in T0}
    if len(c1) != M:
        return False
    c2 = set(T)
    for x, y in itertools.combinations(sorted(T), 2):
        c2.add(x ^ y)
    for x, y in itertools.combinations(sorted(U), 2):
        c2.add(x ^ y)
    return all(z in c2 for z in range(1, M))


def move_functional_last(S, g, n):
    """Image of S under an invertible linear map whose last output coordinate is g."""
    j = next(i for i in range(n) if (g >> i) & 1)
    rows = [1 << i for i in range(n) if i != j] + [g]
    return sorted(sum(parity(r & x) << k for k, r in enumerate(rows)) for x in S)


def a_range(n, N):
    """(a_lo, wmax, a_avg): L2 bounds and the averaging bound for the minimum weight."""
    half = 1 << (n - 1)
    a_lo = next((a for a in range(0, N + 1) if a * (N + 1 - a) >= half), None)
    if a_lo is None:
        return None, None, None
    return a_lo, N + 1 - a_lo, (N * half) // ((1 << n) - 1)


# ---------------------------------------------------------------- normal forms for U

def permute_bits(y, p):
    out = 0
    for i, pi in enumerate(p):
        if (y >> i) & 1:
            out |= 1 << pi
    return out


def candidates(a, m):
    """[(r, U)]: one U per orbit of the coordinate permutations of e_1..e_r on the normal forms
    {0, e_1..e_r} u Y (Y: |Y| = a - r - 1 vectors of weight >= 2 in <e_1..e_r>). Every a-subset
    of F_2^m is affinely equivalent to one of these."""
    out = []
    for r in range(0, m + 1):
        k = a - r - 1
        if k < 0:
            continue
        pool = [y for y in range(1 << r) if popcount(y) >= 2]
        if len(pool) < k:
            continue
        base = [0] + [1 << i for i in range(r)]
        if k == 0:
            out.append((r, sorted(base)))
            continue
        perms = list(itertools.permutations(range(r)))
        P = np.array([[permute_bits(y, p) for y in range(1 << r)] for p in perms], dtype=np.int64)
        w = np.array([(1 << r) ** (k - 1 - i) for i in range(k)], dtype=np.int64)
        seen = set()
        for Y in itertools.combinations(pool, k):
            img = np.sort(P[:, list(Y)], axis=1)
            key = int((img @ w).min())
            if key in seen:
                continue
            seen.add(key)
            out.append((r, sorted(base + list(Y))))
    return out


# ---------------------------------------------------------------- the CP-SAT model

def build(n, N, a, U, fixT=None):
    from ortools.sat.python import cp_model
    m = n - 1
    M = 1 << m
    a_lo, wmax, _ = a_range(n, N)
    if wmax is None:  # N too small for any saturating set; keep only the trivial bound
        wmax = N
    model = cp_model.CpModel()
    t = [None] + [model.NewBoolVar(f"t{z}") for z in range(1, M)]
    model.Add(sum(t[1:]) == N - a)
    if not U:  # (C1) cannot hold: make the model infeasible explicitly
        dummy = model.NewBoolVar("empty_U")
        model.Add(dummy == 1)
        model.Add(dummy == 0)
        return model, t
    Uset = set(U)
    # (C1): every y in F_2^m is u + t' with u in U, t' in T u {0}
    for y in range(M):
        if y in Uset:
            continue
        model.AddBoolOr([t[y ^ u] for u in U])
    # (C2): every nonzero z is in T, in T (+) T, or in U (+) U
    UU = {x ^ y for x, y in itertools.combinations(U, 2)}
    for z in range(1, M):
        if z in UU:
            continue
        lits = [t[z]]
        for v in range(1, M):
            v2 = v ^ z
            if v < v2:
                p = model.NewBoolVar("")
                model.AddImplication(p, t[v])
                model.AddImplication(p, t[v2])
                lits.append(p)
        model.AddBoolOr(lits)
    # (C3): weights of all other functionals in [a, wmax]
    for e in (0, 1):
        for h in range(1, M):
            const = sum(1 for u in U if parity(h & u) ^ e)
            terms = [t[z] for z in range(1, M) if parity(h & z)]
            model.Add(sum(terms) >= a - const)
            model.Add(sum(terms) <= wmax - const)
    if fixT is not None:
        fs = set(fixT)
        for z in range(1, M):
            model.Add(t[z] == (1 if z in fs else 0))
    return model, t


def solve(model, t, seconds, workers):
    from ortools.sat.python import cp_model
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
        return "FEASIBLE", [z for z in range(1, len(t)) if solver.Value(t[z])], dt
    if st == cp_model.INFEASIBLE:
        return "INFEASIBLE", None, dt
    return "UNKNOWN", None, dt


# ---------------------------------------------------------------- the search

def search(n, N, a_values, seconds, workers, log, deadline=None, stop_on_found=True,
           retry_seconds=None, progress=None):
    """Run every (a, U) case. Returns dict with cases, found sets, counts. progress(cases,
    found) is called after every case (used to write partial results)."""
    m = n - 1
    cases = []
    found = []
    bugs = []
    for a in a_values:
        if a > N or a > (1 << m):
            continue
        for r, U in candidates(a, m):
            case = {"a": a, "r": r, "U": U}
            if deadline is not None and time.time() > deadline:
                case["status"] = "NOT-RUN"
                cases.append(case)
                continue
            model, t = build(n, N, a, U)
            st, T, dt = solve(model, t, seconds, workers)
            case.update({"status": st, "s": round(dt, 2)})
            if T is not None:
                S = join(T, U, n)
                ok = len(S) == N and len(set(S)) == N and is_saturating(S, n)
                case["check"] = ok
                if ok:
                    found.append({"a": a, "U": U, "S": S})
                else:
                    bugs.append(f"n={n} N={N} a={a} U={U}: solver set fails the brute-force check")
            cases.append(case)
            if log:
                log(f"n={n} N={N} a={a} r={r} U={U}: {st} {dt:.1f}s")
            if progress:
                progress(cases, found)
            if found and stop_on_found:
                return {"cases": cases, "found": found, "bugs": bugs}
    if retry_seconds:
        for case in cases:
            if case["status"] != "UNKNOWN":
                continue
            if deadline is not None and time.time() > deadline:
                break
            model, t = build(n, N, case["a"], case["U"])
            st, T, dt = solve(model, t, retry_seconds, workers)
            case.update({"status": st, "s_retry": round(dt, 2)})
            if T is not None:
                S = join(T, case["U"], n)
                ok = len(S) == N and len(set(S)) == N and is_saturating(S, n)
                case["check"] = ok
                if ok:
                    found.append({"a": case["a"], "U": case["U"], "S": S})
                else:
                    bugs.append(f"n={n} N={N} a={case['a']} U={case['U']}: retry set fails check")
            if log:
                log(f"retry n={n} N={N} a={case['a']} U={case['U']}: {st} {dt:.1f}s")
            if progress:
                progress(cases, found)
            if found and stop_on_found:
                break
    return {"cases": cases, "found": found, "bugs": bugs}


def counts(cases):
    out = {}
    for c in cases:
        d = out.setdefault(str(c["a"]), {"cases": 0, "INFEASIBLE": 0, "FEASIBLE": 0,
                                         "UNKNOWN": 0, "NOT-RUN": 0})
        d["cases"] += 1
        d[c["status"]] += 1
    return out


KNOWN = {4: 5, 5: 9, 6: 13, 7: 19}  # exact l(n,2): n = 4, 5 classical; 6, 7 Struik 1994


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="l82.json")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--t-case", type=float, default=30)
    ap.add_argument("--t-retry", type=float, default=300)
    ap.add_argument("--t-sanity", type=float, default=60)
    ap.add_argument("--deadline-min", type=float, default=250)
    ap.add_argument("--sanity-min", type=float, default=50)
    args = ap.parse_args()
    t_start = time.time()
    deadline = t_start + 60 * args.deadline_min
    res = {"job": "l82sat-v1", "question": "is there a 25-point 1-saturating set in PG(7,2)?",
           "made_with": "AI agents working in MV2's name; no human has checked this", "bugs": []}
    bugs = res["bugs"]

    def log(msg):
        print(f"[{time.time() - t_start:8.1f}s] {msg}", flush=True)

    def dump():
        res["seconds"] = round(time.time() - t_start, 1)
        with open(args.out, "w") as fh:
            json.dump(res, fh, indent=1)

    # sanity: reproduce the known exact values l(5,2) = 9, l(6,2) = 13, l(7,2) = 19
    # (at most --sanity-min minutes; cases not reached count as NOT-RUN -> SANITY-INCOMPLETE)
    sanity = {}
    sanity_deadline = t_start + 60 * args.sanity_min
    for n in (5, 6, 7):
        L = KNOWN[n]
        ent = {}
        for N, want in ((L - 1, "none"), (L, "found")):
            a_lo, wmax, a_avg = a_range(n, N)
            avals = list(range(a_lo, a_avg + 1)) if a_lo is not None else []
            out = search(n, N, avals, args.t_sanity, args.workers, None,
                         deadline=sanity_deadline, stop_on_found=True)
            bugs.extend(out["bugs"])
            got = {"a_values": avals, "counts": counts(out["cases"]),
                   "found": len(out["found"]),
                   "example": out["found"][0]["S"] if out["found"] else None}
            if want == "none":
                if out["found"]:
                    bugs.append(f"sanity: a {N}-point saturating set in F_2^{n} contradicts "
                                f"l({n},2) = {L}")
                got["ok"] = (not out["found"]
                             and all(c["status"] == "INFEASIBLE" for c in out["cases"]))
            else:
                if not out["found"] and all(c["status"] == "INFEASIBLE" for c in out["cases"]):
                    bugs.append(f"sanity: every case infeasible for N={N} in F_2^{n}, but "
                                f"l({n},2) = {L}")
                got["ok"] = bool(out["found"])
            ent[str(N)] = got
            log(f"sanity n={n} N={N}: {got['counts']} found={got['found']} ok={got['ok']}")
        sanity[str(n)] = ent
    res["sanity"] = sanity
    res["sanity_all_ok"] = all(g["ok"] for e in sanity.values() for g in e.values())
    dump()

    # main: n = 8, N = 25, a = 7, 8, 9 (Ytrehus-Helleseth), then a = 10 (Griesmer only)
    n, N = 8, 25
    a_lo, wmax, a_avg = a_range(n, N)
    res["main"] = {"n": n, "N": N, "a_lo": a_lo, "wmax": wmax, "a_avg": a_avg,
                   "a_values": [7, 8, 9, 10]}
    assert (a_lo, wmax) == (7, 19)
    cand_counts = {str(a): len(candidates(a, n - 1)) for a in (7, 8, 9, 10)}
    res["main"]["candidates"] = cand_counts
    log(f"candidates per a: {cand_counts}")
    dump()
    allcases = []
    found = []
    for block in ([7, 8, 9], [10]):  # first pass over the block, then retries of UNKNOWN cases
        done = list(allcases)

        def progress(cases, fnd, done=done):
            res["main"]["cases"] = done + cases
            res["main"]["counts"] = counts(done + cases)
            res["found"] = found + fnd
            dump()

        out = search(n, N, block, args.t_case, args.workers, log, deadline=deadline,
                     stop_on_found=True, retry_seconds=args.t_retry, progress=progress)
        bugs.extend(out["bugs"])
        allcases += out["cases"]
        found += out["found"]
        res["main"]["counts"] = counts(allcases)
        res["main"]["cases"] = allcases
        res["found"] = found
        dump()
        if found:
            break

    c = counts(allcases)

    def closed(a):
        d = c.get(str(a))
        return d is not None and d["INFEASIBLE"] == d["cases"] == cand_counts[str(a)]

    if bugs:
        verdict = "BUG"
    elif found:  # a brute-force-verified set is a certificate on its own
        verdict = "L82-EQ-25"
    elif not res["sanity_all_ok"]:
        verdict = "SANITY-INCOMPLETE"
    elif all(closed(a) for a in (7, 8, 9, 10)):
        verdict = "L82-EQ-26"
    elif all(closed(a) for a in (7, 8, 9)):
        verdict = "L82-EQ-26-YH"
    else:
        verdict = "UNDECIDED"
    res["verdict"] = verdict
    res["open_cases"] = [{"a": x["a"], "U": x["U"], "status": x["status"]} for x in allcases
                         if x["status"] in ("UNKNOWN", "NOT-RUN")]
    dump()
    log(f"verdict {verdict}; counts {c}; found {len(found)}; bugs {bugs}")


if __name__ == "__main__":
    main()
