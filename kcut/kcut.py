#!/usr/bin/env python3
# kcut.py - solver-free exhaustive search for small Kakeya sets in F_3^n through one hyperplane cut.
# MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name;
# no human mathematician has checked it yet. Licence: Apache-2.0. Python standard library only.
"""
Write F_3^n = F_3^m x F_3 (m = n - 1). The parallel class {s = 0}, {s = 1}, {s = 2} splits a set K
into three sections A (s = 0), B (s = 1), C (s = 2), all subsets of F_3^m. K is a Kakeya set iff

  (1) for every d in F_3^m there is a in A with a + d in B and a + 2d in C
      (this is a line of K in the direction (d, 1); these are all directions not parallel to the class),
  (2) every direction of F_3^m has a full line inside A, inside B or inside C
      (the directions parallel to the class).

Symmetries used (all are affine bijections of F_3^n, so they map Kakeya sets to Kakeya sets):
  * (x, s) -> (x, +-s + t) permutes the three positions, so any section may be put at s = 0, 1, 2;
  * (x, s) -> (g x, s) with g in AGL(m, 3), so A may be taken as a fixed orbit representative;
  * (x, s) -> (x - s v, s) (a shear) fixes A and translates B by -v, so B may be taken to contain 0.

So: if a Kakeya set has a class with section sizes (a, b, c), then for some orbit representative A of
a-subsets there is a b-set B containing 0 and a set C with |C| <= c satisfying (1) and (2).
`CutSearch(S, A, c).run(b)` decides this exactly by enumerating every such B and running an exact
set-cover search for C (no solver, no heuristics, no time limits).

For d in F_3^m the points c that serve d are W[d] = {-(a + b) : a in A, b in B, b - a = d}; the line
is (a, a + d, a + 2d) and a + 2d = -(a + b) over F_3. Condition (1) says C meets every W[d].
"""
import argparse
import itertools
import json
import os
import sys
import time
from math import comb


# ------------------------------------------------------------------ the space F_3^m

class Space:
    def __init__(self, m):
        self.m = m
        P = 3 ** m
        self.P = P
        self.FULL = (1 << P) - 1
        self.pts = [tuple((i // 3 ** k) % 3 for k in range(m)) for i in range(P)]
        idx = {p: i for i, p in enumerate(self.pts)}
        self.idx = idx
        self.add = [[idx[tuple((x + y) % 3 for x, y in zip(self.pts[i], self.pts[j]))]
                     for j in range(P)] for i in range(P)]
        self.neg = [idx[tuple((-x) % 3 for x in p)] for p in self.pts]
        self.sub = [[self.add[i][self.neg[j]] for j in range(P)] for i in range(P)]  # i - j
        self.negsum = [[self.neg[self.add[i][j]] for j in range(P)] for i in range(P)]  # -(i + j)
        self.scal = [[0] * P, list(range(P)), self.neg[:]]  # scal[x][v] = x * v, x in {0, 1, 2}
        dirs, dir_of = [], {}
        for v in range(1, P):
            if v in dir_of:
                continue
            dir_of[v] = len(dirs)
            dir_of[self.neg[v]] = len(dirs)
            dirs.append(v)
        self.dirs = dirs
        self.dir_of = dir_of
        self.ndirs = len(dirs)
        self.ALLDIRS = (1 << self.ndirs) - 1
        self.lines = []
        for v in dirs:
            ls = set()
            for p in range(P):
                q = self.add[p][v]
                r = self.add[q][v]
                ls.add((1 << p) | (1 << q) | (1 << r))
            self.lines.append(sorted(ls))


def mask_of(it):
    m = 0
    for i in it:
        m |= 1 << i
    return m


def bits(mask):
    out = []
    while mask:
        low = mask & -mask
        out.append(low.bit_length() - 1)
        mask ^= low
    return out


def line_dirs(S, mask):
    """Bitmask of the directions k such that some line of direction k lies inside the set `mask`."""
    pts = bits(mask)
    out = 0
    for i in range(len(pts)):
        x = pts[i]
        for j in range(i + 1, len(pts)):
            y = pts[j]
            z = S.negsum[x][y]
            if (mask >> z) & 1:
                out |= 1 << S.dir_of[S.sub[y][x]]
    return out


# ------------------------------------------------------------------ exact search for C

def solve_C(S, W, H, need, mC):
    """Find C (bitmask), |C| <= mC, meeting every W[d] and containing a line in every direction of
    `need`. H[c] is the set of d with c in W[d]. Returns a mask or None. Exact (complete) search."""
    P = S.P
    FULL = S.FULL
    lines = S.lines
    needk = [k for k in range(S.ndirs) if (need >> k) & 1]

    def phase2(Cm, size, forb):
        target = -1
        for k in needk:
            for L in lines[k]:
                if (L & Cm) == L:
                    break
            else:
                target = k
                break
        if target < 0:
            return Cm
        for L in lines[target]:
            if L & forb:
                continue
            cost = (L & ~Cm).bit_count()
            if size + cost <= mC:
                r = phase2(Cm | L, size + cost, forb)
                if r is not None:
                    return r
        return None

    def phase1(Cm, covd, size, forb):
        if covd == FULL:
            return phase2(Cm, size, forb)
        if size >= mC:
            return None
        items = []
        for d in range(P):
            if not ((covd >> d) & 1):
                w = W[d] & ~forb
                if not w:
                    return None
                items.append((w.bit_count(), w))
        items.sort()
        U = 0
        lb = 0
        for _, w in items:
            if not (w & U):
                lb += 1
                U |= w
        if size + lb > mC:
            return None
        w = items[0][1]
        f = forb
        while w:
            low = w & -w
            c = low.bit_length() - 1
            r = phase1(Cm | low, covd | H[c], size + 1, f)
            if r is not None:
                return r
            f |= low  # later branches: c is not in C (disjoint branches)
            w ^= low
        return None

    return phase1(0, 0, 0, 0)


class CutSearch:
    def __init__(self, S, A, mC):
        self.S = S
        self.A = sorted(A)
        self.Am = mask_of(self.A)
        self.mC = mC
        P = S.P
        self.md = [0] * P  # md[x] = set of d with x in A + d
        for x in range(P):
            mm = 0
            for a in self.A:
                mm |= 1 << S.sub[x][a]
            self.md[x] = mm
        self.dirsA = line_dirs(S, self.Am)

    def tables(self, Blist):
        S = self.S
        W = [0] * S.P
        H = [0] * S.P
        sub = S.sub
        negsum = S.negsum
        for a in self.A:
            for b in Blist:
                d = sub[b][a]
                c = negsum[a][b]
                W[d] |= 1 << c
                H[c] |= 1 << d
        return W, H

    def solve(self, Bm):
        """A mask C with |C| <= mC such that (A, B, C) satisfies (1) and (2), or None."""
        S = self.S
        W, H = self.tables(bits(Bm))
        for w in W:
            if not w:
                return None
        need = S.ALLDIRS & ~self.dirsA & ~line_dirs(S, Bm)
        return solve_C(S, W, H, need, self.mC)

    def run(self, b, stop_at_first=True, with_origin=True):
        """Enumerate every b-set B (containing 0 if with_origin) with B - A = F_3^m and solve for C."""
        S = self.S
        P = S.P
        FULL = S.FULL
        md = self.md
        amax = len(self.A)
        stats = {"B_cover": 0}
        found = []

        def rec(start, left, cov, Bm):
            if left == 0:
                if cov == FULL:
                    stats["B_cover"] += 1
                    C = self.solve(Bm)
                    if C is not None:
                        found.append((Bm, C))
                return
            if cov.bit_count() + left * amax < P:
                return
            for x in range(start, P - left + 1):
                rec(x + 1, left - 1, cov | md[x], Bm | (1 << x))
                if found and stop_at_first:
                    return

        if with_origin:
            rec(1, b - 1, md[0], 1)
        else:
            rec(0, b, 0, 0)
        return {"B_cover": stats["B_cover"], "found": found[0] if found else None,
                "n_found": len(found)}


# ------------------------------------------------------------------ group and orbits

def affine_perms(S):
    """All elements of AGL(m, 3) as permutations (bytes) of the P points; also the linear maps."""
    P, m = S.P, S.m
    lin = []
    for cols in itertools.product(range(P), repeat=m):
        img = []
        for i in range(P):
            v = 0
            for k, x in enumerate(S.pts[i]):
                v = S.add[v][S.scal[x][cols[k]]]
            img.append(v)
        if len(set(img)) == P:
            lin.append(img)
    perms = []
    for img in lin:
        for t in range(P):
            perms.append(bytes(S.add[t][j] for j in img))
    return lin, perms


def orbit_reps(S, perms, k):
    """One representative (sorted list) of every AGL(m,3)-orbit of k-subsets, with orbit sizes."""
    seen = set()
    reps = []
    for cmb in itertools.combinations(range(S.P), k):
        mk = mask_of(cmb)
        if mk in seen:
            continue
        orb = set()
        for p in perms:
            mm = 0
            for i in cmb:
                mm |= 1 << p[i]
            orb.add(mm)
        seen |= orb
        reps.append((list(cmb), len(orb)))
    total = sum(o for _, o in reps)
    if total != comb(S.P, k):
        raise AssertionError("orbit sizes add up to %d, not C(%d,%d)" % (total, S.P, k))
    return reps


# ------------------------------------------------------------------ independent checks

def is_kakeya(K, n):
    """Plain definition: for every direction v of F_3^n some p in K has p + v and p + 2v in K."""
    Ks = set(tuple(p) for p in K)
    for v in itertools.product(range(3), repeat=n):
        nz = [x for x in v if x]
        if not nz or nz[0] != 1:
            continue
        ok = False
        for p in Ks:
            q = tuple((a + b) % 3 for a, b in zip(p, v))
            r = tuple((a + 2 * b) % 3 for a, b in zip(p, v))
            if q in Ks and r in Ks:
                ok = True
                break
        if not ok:
            return False
    return True


def assemble(S, Am, Bm, Cm):
    """The set K in F_3^(m+1) with sections A (s=0), B (s=1), C (s=2); masks in, tuples out."""
    K = []
    for s, mk in enumerate((Am, Bm, Cm)):
        for i in bits(mk):
            K.append(tuple(S.pts[i]) + (s,))
    return K


def sections(K, n):
    """Sections of K by the last coordinate, as index lists of F_3^(n-1)."""
    S = Space(n - 1)
    out = [[], [], []]
    for p in K:
        out[p[-1]].append(S.idx[tuple(p[:-1])])
    return [sorted(x) for x in out]


KNOWN27 = [(0, 0, 0, 0), (1, 0, 0, 0), (2, 0, 0, 0), (2, 2, 0, 0), (0, 0, 1, 0), (1, 0, 1, 0),
           (0, 1, 1, 0), (1, 1, 1, 0), (0, 2, 1, 0), (2, 2, 1, 0), (0, 0, 2, 0), (0, 1, 2, 0),
           (1, 2, 2, 0), (1, 2, 0, 1), (1, 0, 1, 1), (0, 1, 1, 1), (1, 2, 1, 1), (0, 0, 2, 1),
           (2, 1, 2, 1), (2, 0, 0, 2), (0, 1, 0, 2), (0, 2, 0, 2), (1, 0, 1, 2), (2, 0, 1, 2),
           (1, 2, 1, 2), (0, 1, 2, 2), (2, 1, 2, 2)]

# Cuts that must occur (moments/BAD25.md, moments/min_section.py):
#  N = 25: some class has cut (6,6,13), (5,7,13) or (6,7,12).
#  N = 26: some class has cut (6,6,14), (5,7,14), (6,7,13) or (7,7,12); otherwise every class has an
#          8-point section and some hyperplane has <= 7 points, so some class has cut (5,8,13),
#          (6,8,12) or (7,8,11).
CUTS25 = [(6, 6, 13), (5, 7, 13), (6, 7, 12)]
CUTS26 = [(6, 6, 14), (5, 7, 14), (6, 7, 13), (5, 8, 13), (7, 7, 12), (6, 8, 12), (7, 8, 11)]


def k3_control(N):
    """End-to-end control in F_3^3 (m = 2): search every cut a <= b <= c of N points.
    Lev's k_3 = 13 means N = 12 must give nothing and N = 13 must give a set."""
    S = Space(2)
    _, perms = affine_perms(S)
    reps = {}
    found = None
    cuts = 0
    for a in range(1, S.P + 1):
        for b in range(a, S.P + 1):
            c = N - a - b
            if c < b or c > S.P:
                continue
            cuts += 1
            if a not in reps:
                reps[a] = orbit_reps(S, perms, a)
            for A, _ in reps[a]:
                r = CutSearch(S, A, c).run(b)
                if r["found"] is not None:
                    Bm, Cm = r["found"]
                    K = assemble(S, mask_of(A), Bm, Cm)
                    if found is None:
                        found = {"cut": [a, b, c], "K": K, "is_kakeya": is_kakeya(K, 3), "size": len(K)}
    return {"N": N, "cuts": cuts, "found": found}


# ------------------------------------------------------------------ the job

_S = None


def _init():
    global _S
    _S = Space(3)


def _task(t):
    tag, cut, A = t
    a, b, c = cut
    t0 = time.time()
    r = CutSearch(_S, A, c).run(b)
    out = {"tag": tag, "cut": list(cut), "A": A, "B_cover": r["B_cover"],
           "seconds": round(time.time() - t0, 1), "found": None}
    if r["found"] is not None:
        Bm, Cm = r["found"]
        K = assemble(_S, mask_of(A), Bm, Cm)
        out["found"] = K
        out["found_size"] = len(K)
        out["found_is_kakeya"] = is_kakeya(K, 4)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="kcut.json")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--skip26", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    res = {"contract": "kcut-v1", "bugs": [], "controls": {}, "per_cut": {}, "tasks": []}

    def log(msg):
        print("[%8.1fs] %s" % (time.time() - t0, msg), flush=True)

    def save():
        tmp = args.out + ".tmp"
        with open(tmp, "w") as f:
            json.dump(res, f, indent=1)
        os.replace(tmp, args.out)

    # controls in F_3^3: k_3 = 13
    for N in (12, 13):
        c = k3_control(N)
        res["controls"]["k3_N%d" % N] = c
        log("k3 control N=%d: %s" % (N, "found" if c["found"] else "none"))
    if res["controls"]["k3_N12"]["found"] is not None:
        res["bugs"].append("k3 control: a 12-point Kakeya set in F_3^3 was reported")
    f13 = res["controls"]["k3_N13"]["found"]
    if f13 is None or not f13["is_kakeya"] or f13["size"] > 13:
        res["bugs"].append("k3 control: no verified 13-point Kakeya set in F_3^3")
    if not is_kakeya(KNOWN27, 4):
        res["bugs"].append("KNOWN27 is not Kakeya")
    save()

    S = Space(3)
    lin, perms = affine_perms(S)
    res["controls"]["GL33_order"] = len(lin)
    if len(lin) != 11232:
        res["bugs"].append("|GL(3,3)| = %d, expected 11232" % len(lin))
    sizes = sorted({a for a, _, _ in CUTS25 + CUTS26})
    reps = {}
    for k in sizes:
        reps[k] = orbit_reps(S, perms, k)
        log("orbit representatives of %d-subsets of AG(3,3): %d" % (k, len(reps[k])))
    res["orbits"] = {str(k): [o for _, o in reps[k]] for k in sizes}
    del perms
    save()

    sec = sections(KNOWN27, 4)  # sizes 13, 6, 8 at s = 0, 1, 2
    tasks = [("control-known27", (6, 8, 13), sec[1])]
    for A, _ in reps[6]:
        tasks.append(("control-24", (6, 6, 12), A))
    cuts = CUTS25 + ([] if args.skip26 else CUTS26)
    for cut in cuts:
        key = "%d,%d,%d" % cut
        res["per_cut"][key] = {"N": sum(cut), "reps": len(reps[cut[0]]), "done": 0, "B_cover": 0,
                               "found": None, "seconds": 0.0}
        for A, _ in reps[cut[0]]:
            tasks.append(("N%d" % sum(cut), cut, A))
    log("%d tasks" % len(tasks))

    from multiprocessing import Pool
    with Pool(args.workers, initializer=_init) as pool:
        for out in pool.imap_unordered(_task, tasks):
            res["tasks"].append({k: out[k] for k in ("tag", "cut", "A", "B_cover", "seconds")}
                                | {"found": out["found"] is not None})
            if out["tag"] == "control-known27":
                ok = out["found"] is not None and out.get("found_is_kakeya") and out["found_size"] <= 27
                res["controls"]["known27_cut_6_8_13"] = {"found": out["found"], "ok": bool(ok)}
                if not ok:
                    res["bugs"].append("control: the cut (6,8,13) of the known 27-point set was not found")
            elif out["tag"] == "control-24":
                if out["found"] is not None:
                    res["bugs"].append("control: a set of <= 24 points with cut (6,6,12) was reported: %s"
                                       % out["found"])
            else:
                pc = res["per_cut"]["%d,%d,%d" % tuple(out["cut"])]
                pc["done"] += 1
                pc["B_cover"] += out["B_cover"]
                pc["seconds"] = round(pc["seconds"] + out["seconds"], 1)
                if out["found"] is not None:
                    if not out.get("found_is_kakeya"):
                        res["bugs"].append("reported set fails is_kakeya: %s" % out["found"])
                    elif pc["found"] is None:
                        pc["found"] = out["found"]
            log("%s cut=%s A=%s B_cover=%d %s %.1fs" % (out["tag"], out["cut"], out["A"], out["B_cover"],
                "FOUND" if out["found"] is not None else "none", out["seconds"]))
            save()

    def verdict(cutlist, N):
        pcs = [res["per_cut"].get("%d,%d,%d" % c) for c in cutlist]
        if any(pc is None for pc in pcs):
            return "NOT-RUN"
        if any(pc["found"] is not None for pc in pcs):
            return "FOUND-%d" % N
        if all(pc["done"] == pc["reps"] for pc in pcs):
            return "CLOSED-%d" % N
        return "INCOMPLETE"

    res["verdict25"] = verdict(CUTS25, 25)
    res["verdict26"] = verdict(CUTS26, 26)
    res["verdict"] = "BUG" if res["bugs"] else ("%s / %s" % (res["verdict25"], res["verdict26"]))
    res["seconds"] = round(time.time() - t0, 1)
    save()
    log("verdict: %s" % res["verdict"])


if __name__ == "__main__":
    main()
