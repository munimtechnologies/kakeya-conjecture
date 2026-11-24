# l82b.py -- l82sat v2: the cases that l82sat-v1 left open, now with symmetry breaking.
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0.
#
# Question: is there a 25-point 1-saturating set in PG(7,2), i.e. is l(8,2) = 25 (Struik's 1994
# lower bound) or 26 (best known construction)?  The reduction, the case split (a, U) and the
# CP-SAT model are those of l82.py (munimtechnologies/kakeya-conjecture, l82/), imported unchanged.
#
# New in v2: after U is fixed, the remaining symmetry is the group G_U of maps
# (y, s) -> (A y + s c, s) with A U + c = U, acting on T by T -> A T. The model of l82.build is
# G_U-invariant (each of (C1), (C2), (C3) and |T| is preserved; proof in README-l82b.md). For any
# subset H of G_U we add lexicographic-leader constraints  t >=_lex t o A  (A = linear part of
# h in H, variables ordered z = 1 .. 2^m - 1). The lex-largest member of every G_U-orbit satisfies
# all of them, so a feasible model stays feasible: the constraints remove only symmetric copies.
# H = (i) the transvections e_j -> e_j + e_i (j outside the span of U), which fix U pointwise,
#     (ii) up to K affine maps that permute U (found by depth-first search, fewest moved points).
# The tests check this on every case of l(5,2) and l(6,2): with and without the constraints the
# models agree on feasibility.

import argparse
import json
import time

import l82
from l82 import a_range, candidates, counts, is_saturating, join, KNOWN

# ---------------------------------------------------------------- linear algebra over F_2


def lin_image(x, cols, c=0):
    """c + A x where A e_i = cols[i]."""
    y = c
    i = 0
    while x:
        if x & 1:
            y ^= cols[i]
        x >>= 1
        i += 1
    return y


def rank(vecs, m):
    piv = [0] * m
    r = 0
    for v in vecs:
        for i in reversed(range(m)):
            if not (v >> i) & 1:
                continue
            if piv[i]:
                v ^= piv[i]
            else:
                piv[i] = v
                r += 1
                break
    return r


def affine_rank(U, m):
    u0 = U[0]
    return rank([u ^ u0 for u in U[1:]], m)


# ---------------------------------------------------------------- the symmetry group G_U


def stabilizer_maps(U, m, node_cap=400000):
    """Non-identity affine maps (cols, c) of F_2^m with A U + c = U and A e_j = e_j outside the
    span of U. U must be in the normal form of l82.candidates: 0 and e_1..e_r in U, U inside
    <e_1..e_r>. Returns (maps, complete)."""
    U = sorted(U)
    Uset = set(U)
    r = affine_rank(U, m)
    assert 0 in Uset and all((1 << i) in Uset for i in range(r)) and max(U) < (1 << r), U
    by_high = [[u for u in U if (1 << k) <= u < (1 << (k + 1))] for k in range(r)]
    maps = []
    state = {"nodes": 0, "complete": True}

    def rec(k, c, cols, span):
        if state["nodes"] >= node_cap:
            state["complete"] = False
            return
        state["nodes"] += 1
        if k == r:
            if c == 0 and all(cols[i] == 1 << i for i in range(r)):
                return
            maps.append((tuple(cols + [1 << j for j in range(r, m)]), c))
            return
        for w in U:
            col = w ^ c
            if col in span:
                continue
            ncols = cols + [col]
            if all(lin_image(u, ncols, c) in Uset for u in by_high[k]):
                rec(k + 1, c, ncols, span | {s ^ col for s in span})

    for c in U:
        rec(0, c, [], {0})
    return maps, state["complete"]


def transvections(r, m):
    """A e_j = e_j + e_i for j >= r (outside span U), i != j; these fix span U pointwise."""
    out = []
    for j in range(r, m):
        for i in range(m):
            if i != j:
                cols = [1 << k for k in range(m)]
                cols[j] ^= 1 << i
                out.append((tuple(cols), 0))
    return out


def perm_of(cols, m):
    return [lin_image(z, cols) for z in range(1 << m)]


def symmetry_perms(U, m, K=40, node_cap=400000):
    """Permutations P of 0..2^m-1 (linear parts of elements of G_U) used for lex-leader."""
    r = affine_rank(sorted(U), m)
    maps, complete = stabilizer_maps(U, m, node_cap)
    perms = []
    seen = set()
    for cols, _ in transvections(r, m):
        P = perm_of(cols, m)
        key = tuple(P)
        if key not in seen and any(P[z] != z for z in range(1 << m)):
            seen.add(key)
            perms.append(P)
    scored = []
    for cols, _ in maps:
        P = perm_of(cols, m)
        moved = sum(1 for z in range(1 << m) if P[z] != z)
        if moved:
            scored.append((moved, P))
    scored.sort(key=lambda x: (x[0], x[1]))
    k = 0
    for _, P in scored:
        if k >= K:
            break
        key = tuple(P)
        if key in seen:
            continue
        seen.add(key)
        perms.append(P)
        k += 1
    return perms, {"stab_maps": len(maps), "stab_complete": complete, "transvections":
                   len(transvections(r, m)), "used": len(perms)}


def add_lex_geq(model, t, P):
    """t >=_lex (t o P) over z = 1 .. len(t)-1 (fixed points of P skipped)."""
    prev = None  # None stands for 'true' (empty prefix is equal)
    for z in range(1, len(t)):
        pz = P[z]
        if pz == z:
            continue
        x, y = t[z], t[pz]
        if prev is None:
            model.AddImplication(y, x)
        else:
            model.AddBoolOr([prev.Not(), x, y.Not()])
        nxt = model.NewBoolVar("")
        if prev is None:
            model.AddBoolOr([x.Not(), y.Not(), nxt])
            model.AddBoolOr([x, y, nxt])
        else:
            model.AddBoolOr([prev.Not(), x.Not(), y.Not(), nxt])
            model.AddBoolOr([prev.Not(), x, y, nxt])
        prev = nxt


def build(n, N, a, U, sym=True, K=40, fixT=None):
    model, t = l82.build(n, N, a, U, fixT=fixT)
    info = {"used": 0}
    if sym and U:
        perms, info = symmetry_perms(U, n - 1, K)
        for P in perms:
            add_lex_geq(model, t, P)
    return model, t, info


# ---------------------------------------------------------------- running cases


def check_found(T, U, n, N):
    S = join(T, U, n)
    return S, (len(S) == N and len(set(S)) == N and is_saturating(S, n))


def run_cases(n, N, cases, schedule, workers, deadline, K, log=None, progress=None,
              stop_on_found=True):
    """cases: dicts with a, r, U, status. schedule: list of per-case limits in seconds; the
    value None means 'share the remaining time equally among the open cases'."""
    found, bugs = [], []
    for p, tl in enumerate(schedule):
        todo = [c for c in cases if c["status"] in ("UNKNOWN", "NOT-RUN")]
        if not todo:
            break
        if tl is None:
            if deadline is None:
                break
            tl = max(60.0, (deadline - time.time()) / len(todo))
        for c in todo:
            if deadline is not None and time.time() + 5 > deadline:
                break
            lim = tl if deadline is None else min(tl, max(5.0, deadline - time.time()))
            model, t, info = build(n, N, c["a"], c["U"], sym=True, K=K)
            st, T, dt = l82.solve(model, t, lim, workers)
            c["status"] = st
            c[f"s{p + 1}"] = round(dt, 2)
            c["sym"] = info
            if T is not None:
                S, ok = check_found(T, c["U"], n, N)
                c["check"] = ok
                if ok:
                    found.append({"a": c["a"], "U": c["U"], "S": S})
                else:
                    bugs.append(f"n={n} N={N} a={c['a']} U={c['U']}: solver set fails check")
            if log:
                log(f"pass {p + 1} n={n} N={N} a={c['a']} r={c['r']} U={c['U']}: {st} "
                    f"{dt:.1f}s (limit {lim:.0f}s, sym {info.get('used')})")
            if progress:
                progress()
            if found and stop_on_found:
                return found, bugs
    return found, bugs


# ---------------------------------------------------------------- which cases v1 left open

# l82sat-v1 (compete-results/l82sat-v1/out/l82.json, 2026-11-24): CP-SAT closed these; every
# other case listed by select_cases was UNKNOWN (a = 7, 8, 9) or NOT-RUN (a = 10).
V1_OPEN_A7 = [[0, 1, 2, 4, 7, 8, 16]]
V1_CLOSED_A8 = [[0, 1, 2, 3, 4, 5, 8, 15], [0, 1, 2, 3, 4, 8, 13, 14], [0, 1, 2, 3, 4, 8, 13, 15]]
V1_OPEN_A9_LOW_R = [[0, 1, 2, 4, 7, 8, 15, 16, 25]]  # the one open a = 9 case with r <= 5


def select_cases(a_list, part=(0, 1), m=7):
    out = []
    for a in a_list:
        for r, U in candidates(a, m):
            if a == 7:
                keep = U in V1_OPEN_A7
            elif a == 8:
                keep = U not in V1_CLOSED_A8
            elif a == 9:
                keep = r >= 6 or U in V1_OPEN_A9_LOW_R
            else:
                keep = True
            if keep:
                out.append({"a": a, "r": r, "U": U, "status": "NOT-RUN"})
    i, k = part
    return [c for j, c in enumerate(out) if j % k == i]


# ---------------------------------------------------------------- main


def sanity(ns, minutes, workers, K, log):
    """With the symmetry constraints on, reproduce l(n,2) for n in ns: every case infeasible at
    l - 1, a verified set found at l."""
    res, bugs = {}, []
    dl = time.time() + 60 * minutes
    for n in ns:
        L = KNOWN[n]
        ent = {}
        for N, want in ((L - 1, "none"), (L, "found")):
            a_lo, _, a_avg = a_range(n, N)
            cases = [{"a": a, "r": r, "U": U, "status": "NOT-RUN"}
                     for a in range(a_lo, a_avg + 1) for r, U in candidates(a, n - 1)]
            found, b = run_cases(n, N, cases, [60, None], workers, dl, K, None, None, True)
            bugs += b
            c = counts(cases)
            got = {"counts": c, "found": len(found), "example": found[0]["S"] if found else None}
            if want == "none":
                if found:
                    bugs.append(f"sanity: {N}-point saturating set in F_2^{n} contradicts l = {L}")
                got["ok"] = not found and all(x["status"] == "INFEASIBLE" for x in cases)
            else:
                if not found and all(x["status"] == "INFEASIBLE" for x in cases):
                    bugs.append(f"sanity: all cases infeasible at N = l({n},2) = {L}: the "
                                f"symmetry constraints cut off every solution")
                got["ok"] = bool(found)
            ent[str(N)] = got
            log(f"sanity n={n} N={N}: {c} found={len(found)} ok={got['ok']}")
        res[str(n)] = ent
    return res, bugs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", default="l82b")
    ap.add_argument("--out", default="l82b.json")
    ap.add_argument("--a", default="7,8")
    ap.add_argument("--part", default="0/1")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--K", type=int, default=40)
    ap.add_argument("--t1", type=float, default=120)
    ap.add_argument("--t2", type=float, default=900)
    ap.add_argument("--deadline-min", type=float, default=330)
    ap.add_argument("--sanity-n", default="5,6,7")
    ap.add_argument("--sanity-min", type=float, default=40)
    args = ap.parse_args()
    t0 = time.time()
    deadline = t0 + 60 * args.deadline_min
    res = {"job": args.job, "question": "is there a 25-point 1-saturating set in PG(7,2)?",
           "made_with": "AI agents working in MV2's name; no human has checked this",
           "args": vars(args), "bugs": []}

    def log(msg):
        print(f"[{time.time() - t0:8.1f}s] {msg}", flush=True)

    def dump():
        res["seconds"] = round(time.time() - t0, 1)
        with open(args.out, "w") as fh:
            json.dump(res, fh, indent=1)

    ns = [int(x) for x in args.sanity_n.split(",") if x]
    san, b = sanity(ns, args.sanity_min, args.workers, args.K, log)
    res["sanity"] = san
    res["sanity_all_ok"] = all(g["ok"] for e in san.values() for g in e.values())
    res["bugs"] += b
    dump()

    a_list = [int(x) for x in args.a.split(",")]
    i, k = (int(x) for x in args.part.split("/"))
    cases = select_cases(a_list, (i, k))
    res["main"] = {"n": 8, "N": 25, "a": a_list, "part": args.part, "n_cases": len(cases),
                   "cases": cases}
    log(f"main: {len(cases)} cases for a = {a_list}, part {args.part}")
    dump()
    found, b = run_cases(8, 25, cases, [args.t1, args.t2, None], args.workers, deadline,
                         args.K, log, dump, True)
    res["bugs"] += b
    res["found"] = found
    res["main"]["counts"] = counts(cases)
    if res["bugs"]:
        verdict = "BUG"
    elif found:
        verdict = "FOUND-25"
    elif not res["sanity_all_ok"]:
        verdict = "SANITY-INCOMPLETE"
    elif all(c["status"] == "INFEASIBLE" for c in cases):
        verdict = "PART-CLOSED"
    else:
        verdict = "PART-OPEN"
    res["verdict"] = verdict
    res["open_cases"] = [{"a": c["a"], "U": c["U"], "status": c["status"]} for c in cases
                         if c["status"] != "INFEASIBLE"]
    dump()
    log(f"verdict {verdict}; counts {res['main']['counts']}; found {len(found)}; "
        f"bugs {res['bugs']}")


if __name__ == "__main__":
    main()
