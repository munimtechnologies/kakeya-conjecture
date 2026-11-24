#!/usr/bin/env python3
# kcutc.py - kcut.py plus a canonical-B symmetry reduction, run on the seven N = 26 cuts that decide
# Lev's conjecture k_4 = 27 for Kakeya sets in F_3^4.
# MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name;
# no human mathematician has checked it yet. Licence: Apache-2.0. Python standard library only.
"""
Why 26 is enough: a Kakeya set plus any extra point is still a Kakeya set, so if no 26-point Kakeya
set exists in F_3^4 then none of 25 or fewer points exists either, and k_4 = 27 (the 27-point set
KNOWN27 exists). Which cuts must be searched for N = 26 comes from the moment identity
    sum over the 40 hyperplane classes of (n1 - 8)(n2 - 8)(n3 - 8) = 18 (L - 40)
(see moments/BAD25.md and moments/PRODUCT.md) together with (T), (SD) and the min-section lemma:
some class has cut (6,6,14), (5,7,14), (6,7,13) or (7,7,12), or else every class has an 8-point
section and one class has cut (5,8,13), (6,8,12) or (7,8,11).

The new reduction. Fix the smallest section A at s = 0 (an AGL(3,3)-orbit representative). For every
g(x) = M x + t in the stabiliser of A and every v, the affine bijection (x, s) -> (M x + t + s v, s)
of F_3^4 fixes A, sends B to M B + (t + v) and C to M C + (t + 2 v). As v is free, the search for C
gives the same answer for B and for every M (B - b) with M in Lin(Stab A) and b in B (all of these
contain 0). kcut.py tried every B containing 0; here B is kept only if its bitmask is the smallest
of these images. Each orbit keeps exactly one representative (its minimum), so the search is still
complete, and it does about |B| * |Lin(Stab A)| times fewer C-searches.
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kcut as K  # noqa: E402


def stabiliser(S, lin, perms, Am):
    """(number of affine maps fixing the set Am, list of their distinct linear parts).
    perms[i] is the affine map with linear part lin[i // P] (see kcut.affine_perms)."""
    P = S.P
    A = K.bits(Am)
    n_aff = 0
    idxs = set()
    for i, p in enumerate(perms):
        mm = 0
        for a in A:
            mm |= 1 << p[a]
        if mm == Am:
            n_aff += 1
            idxs.add(i // P)
    return n_aff, [lin[j] for j in sorted(idxs)]


def is_canonical(S, Bm, Ms):
    """True iff Bm is the smallest bitmask among M (B - b), M in Ms (linear maps), b in B."""
    pts = K.bits(Bm)
    sub = S.sub
    for M in Ms:
        img = [M[x] for x in pts]
        for mb in img:
            mm = 0
            for y in img:
                mm |= 1 << sub[y][mb]
            if mm < Bm:
                return False
    return True


def run_canon(cs, b, Ms, stop_at_first=True):
    """Like CutSearch.run(b) (B contains 0), but solves for C only for canonical B."""
    S = cs.S
    P = S.P
    FULL = S.FULL
    md = cs.md
    amax = len(cs.A)
    stats = {"B_cover": 0, "B_canon": 0}
    found = []

    def rec(start, left, cov, Bm):
        if left == 0:
            if cov == FULL:
                stats["B_cover"] += 1
                if is_canonical(S, Bm, Ms):
                    stats["B_canon"] += 1
                    C = cs.solve(Bm)
                    if C is not None:
                        found.append((Bm, C))
            return
        if cov.bit_count() + left * amax < P:
            return
        for x in range(start, P - left + 1):
            rec(x + 1, left - 1, cov | md[x], Bm | (1 << x))
            if found and stop_at_first:
                return

    rec(1, b - 1, md[0], 1)
    return {"B_cover": stats["B_cover"], "B_canon": stats["B_canon"],
            "found": found[0] if found else None, "n_found": len(found)}


def k3_control_canon(N):
    """End-to-end control in F_3^3 with the reduction: N = 12 must give nothing, N = 13 a set."""
    S = K.Space(2)
    lin, perms = K.affine_perms(S)
    found = None
    cuts = 0
    for a in range(1, S.P + 1):
        reps = None
        for b in range(a, S.P + 1):
            c = N - a - b
            if c < b or c > S.P:
                continue
            cuts += 1
            if reps is None:
                reps = K.orbit_reps(S, perms, a)
            for A, _ in reps:
                Am = K.mask_of(A)
                _, Ms = stabiliser(S, lin, perms, Am)
                r = run_canon(K.CutSearch(S, A, c), b, Ms)
                if r["found"] is not None and found is None:
                    Bm, Cm = r["found"]
                    KK = K.assemble(S, Am, Bm, Cm)
                    found = {"cut": [a, b, c], "K": KK, "is_kakeya": K.is_kakeya(KK, 3), "size": len(KK)}
    return {"N": N, "cuts": cuts, "found": found}


# ------------------------------------------------------------------ the job

_S = None


def _init():
    global _S
    _S = K.Space(3)


def _task(t):
    tag, cut, A, Ms = t
    a, b, c = cut
    t0 = time.time()
    r = run_canon(K.CutSearch(_S, A, c), b, Ms)
    out = {"tag": tag, "cut": list(cut), "A": A, "B_cover": r["B_cover"], "B_canon": r["B_canon"],
           "lin_stab": len(Ms), "seconds": round(time.time() - t0, 1), "found": None}
    if r["found"] is not None:
        Bm, Cm = r["found"]
        KK = K.assemble(_S, K.mask_of(A), Bm, Cm)
        out["found"] = KK
        out["found_size"] = len(KK)
        out["found_is_kakeya"] = K.is_kakeya(KK, 4)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="kcutc.json")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--cuts", default="26", help="'26' (the seven N = 26 cuts) or e.g. '7,8,11;6,8,12'")
    ap.add_argument("--no-control24", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    res = {"contract": "kcutc-v1", "bugs": [], "controls": {}, "per_cut": {}, "tasks": []}

    def log(msg):
        print("[%8.1fs] %s" % (time.time() - t0, msg), flush=True)

    def save():
        tmp = args.out + ".tmp"
        with open(tmp, "w") as f:
            json.dump(res, f, indent=1)
        os.replace(tmp, args.out)

    if args.cuts == "26":
        cuts = list(K.CUTS26)
    else:
        cuts = [tuple(int(x) for x in s.split(",")) for s in args.cuts.split(";")]

    for N in (12, 13):
        c = k3_control_canon(N)
        res["controls"]["k3_canon_N%d" % N] = c
        log("k3 control (canonical B) N=%d: %s" % (N, "found" if c["found"] else "none"))
    if res["controls"]["k3_canon_N12"]["found"] is not None:
        res["bugs"].append("k3 control: a 12-point Kakeya set in F_3^3 was reported")
    f13 = res["controls"]["k3_canon_N13"]["found"]
    if f13 is None or not f13["is_kakeya"] or f13["size"] > 13:
        res["bugs"].append("k3 control: no verified 13-point Kakeya set in F_3^3")
    save()

    S = K.Space(3)
    lin, perms = K.affine_perms(S)
    if len(lin) != 11232 or len(perms) != 303264:
        res["bugs"].append("|GL(3,3)| = %d, |AGL(3,3)| = %d" % (len(lin), len(perms)))
    sizes = sorted({cut[0] for cut in cuts} | {6})
    reps = {}
    stabs = {}
    for k in sizes:
        reps[k] = K.orbit_reps(S, perms, k)
        for A, orb in reps[k]:
            n_aff, Ms = stabiliser(S, lin, perms, K.mask_of(A))
            if n_aff * orb != len(perms):
                res["bugs"].append("orbit-stabiliser fails for A=%s: %d * %d" % (A, n_aff, orb))
            stabs[tuple(A)] = Ms
        log("%d-subsets of AG(3,3): %d orbits, linear stabiliser orders %s"
            % (k, len(reps[k]), [len(stabs[tuple(A)]) for A, _ in reps[k]]))
    res["orbits"] = {str(k): [[o, len(stabs[tuple(A)])] for A, o in reps[k]] for k in sizes}
    # control: the known 27-point set has a class with cut (6, 8, 13); A = its 6-point section is not an
    # orbit representative, which also tests the reduction for an arbitrary A.
    sec = K.sections(K.KNOWN27, 4)
    _, Ms27 = stabiliser(S, lin, perms, K.mask_of(sec[1]))
    del perms
    save()

    tasks = [("control-known27", (6, 8, 13), sec[1], Ms27)]
    for cut in cuts:
        key = "%d,%d,%d" % cut
        res["per_cut"][key] = {"N": sum(cut), "reps": len(reps[cut[0]]), "done": 0, "B_cover": 0,
                               "B_canon": 0, "found": None, "seconds": 0.0}
        for A, _ in reps[cut[0]]:
            tasks.append(("N%d" % sum(cut), cut, A, stabs[tuple(A)]))
    if not args.no_control24:
        for A, _ in reps[6]:
            tasks.append(("control-24", (6, 6, 12), A, stabs[tuple(A)]))
    log("%d tasks" % len(tasks))
    save()

    from multiprocessing import Pool
    with Pool(args.workers, initializer=_init) as pool:
        for out in pool.imap_unordered(_task, tasks):
            res["tasks"].append({k: out[k] for k in ("tag", "cut", "A", "B_cover", "B_canon", "lin_stab",
                                                     "seconds")} | {"found": out["found"] is not None})
            if out["tag"] == "control-known27":
                ok = out["found"] is not None and out.get("found_is_kakeya") and out["found_size"] <= 27
                res["controls"]["known27_cut_6_8_13"] = {"found": out["found"], "ok": bool(ok)}
                if not ok:
                    res["bugs"].append("control: the cut (6,8,13) of the known 27-point set was not found")
            elif out["tag"] == "control-24":
                res["controls"].setdefault("control24_done", 0)
                res["controls"]["control24_done"] += 1
                if out["found"] is not None:
                    res["bugs"].append("control: a set of <= 24 points with cut (6,6,12) was reported: %s"
                                       % out["found"])
            else:
                pc = res["per_cut"]["%d,%d,%d" % tuple(out["cut"])]
                pc["done"] += 1
                pc["B_cover"] += out["B_cover"]
                pc["B_canon"] += out["B_canon"]
                pc["seconds"] = round(pc["seconds"] + out["seconds"], 1)
                if out["found"] is not None:
                    if not out.get("found_is_kakeya"):
                        res["bugs"].append("reported set fails is_kakeya: %s" % out["found"])
                    elif pc["found"] is None:
                        pc["found"] = out["found"]
            log("%s cut=%s A=%s B_cover=%d B_canon=%d %s %.1fs" % (
                out["tag"], out["cut"], out["A"], out["B_cover"], out["B_canon"],
                "FOUND" if out["found"] is not None else "none", out["seconds"]))
            save()

    pcs = list(res["per_cut"].values())
    if any(pc["found"] is not None for pc in pcs):
        v = "FOUND"
    elif all(pc["done"] == pc["reps"] for pc in pcs):
        v = "CLOSED"
    else:
        v = "INCOMPLETE"
    if args.cuts == "26":
        v += "-26"
    res["verdict26"] = v
    res["verdict"] = "BUG" if res["bugs"] else v
    res["seconds"] = round(time.time() - t0, 1)
    save()
    log("verdict: %s" % res["verdict"])


if __name__ == "__main__":
    main()
