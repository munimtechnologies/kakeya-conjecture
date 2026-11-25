# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Focused local search for small Kakeya sets in F_3^n (standard library only).

The smallest Kakeya set known in F_3^5 has 53 points (Chung, Du and Wesley, arXiv:2608.23691),
matching the size that V. Lev's 2009 recurrence k_n = k_(n-1) + 2 k_(n-2) predicts.  This
module looks for smaller ones.  A set of 52 points would break that prediction; failing to
find one is only weak evidence, never a proof.

Search (WalkSAT style, at a fixed size budget m):
  * state: a set K, optionally a union of orbits of a small group G of linear maps;
  * u = number of directions (of the (3^n - 1)/2) with no line inside K;
  * step: pick an uncovered direction at random, complete one of its lines that needs the
    fewest new points (ties broken by how many uncovered directions the addition covers),
    then remove orbits of least cost (cost = directions that would lose their last line)
    until |K| <= m, never removing what was added in the last `tenure` steps;
  * when u = 0 the set is recorded and m drops by one;
  * after `stall` steps without progress: restart from the smallest set found on the current
    trajectory minus a random orbit, and every `fresh`-th restart from a new greedy set.
Every reported set is re-checked by is_kakeya, which uses no table of this module.
"""
import itertools
import random
import time


# ---- tables ------------------------------------------------------------------------------
def setup(n):
    pts = list(itertools.product(range(3), repeat=n))
    idx = {p: i for i, p in enumerate(pts)}
    dirs = [v for v in pts if any(v) and [c for c in v if c][0] == 1]
    P, D = len(pts), len(dirs)
    LA = [0] * (P * D)
    LB = [0] * (P * D)
    for i, p in enumerate(pts):
        b = i * D
        for j, d in enumerate(dirs):
            LA[b + j] = idx[tuple((x + y) % 3 for x, y in zip(p, d))]
            LB[b + j] = idx[tuple((x + 2 * y) % 3 for x, y in zip(p, d))]
    lines = []
    for j in range(D):
        seen = set()
        L = []
        for i in range(P):
            t = tuple(sorted((i, LA[i * D + j], LB[i * D + j])))
            if t not in seen:
                seen.add(t)
                L.append(t)
        lines.append(L)
    return {"n": n, "pts": pts, "idx": idx, "dirs": dirs, "P": P, "D": D,
            "LA": LA, "LB": LB, "lines": lines}


# ---- symmetry groups (generators; linear bijections of F_3^n, n >= 3) --------------------
def _neg(p):
    return tuple(-c for c in p)


def _refl1(p):
    return (-p[0],) + tuple(p[1:])


def _shear(p):
    return (p[0], p[1] + p[0]) + tuple(p[2:])


def _cyc3(p):
    return (p[1], p[2], p[0]) + tuple(p[3:])


GROUPS = {"trivial": [], "neg": [_neg], "refl1": [_refl1], "shear": [_shear],
          "cyc3": [_cyc3]}


def orbits(T, gens):
    idx, pts, P = T["idx"], T["pts"], T["P"]
    perms = [[idx[tuple(c % 3 for c in f(p))] for p in pts] for f in gens]
    orb = [-1] * P
    members = []
    for i in range(P):
        if orb[i] >= 0:
            continue
        o = len(members)
        orb[i] = o
        mem = [i]
        stack = [i]
        while stack:
            x = stack.pop()
            for pm in perms:
                y = pm[x]
                if orb[y] < 0:
                    orb[y] = o
                    mem.append(y)
                    stack.append(y)
        members.append(mem)
    return orb, members


# ---- independent check -------------------------------------------------------------------
def is_kakeya(n, K):
    """Does K (iterable of point tuples) contain a line in every direction of F_3^n?"""
    S = set(tuple(x) for x in K)
    for d in itertools.product(range(3), repeat=n):
        nz = [c for c in d if c]
        if not nz or nz[0] != 1:
            continue
        ok = False
        for x in S:
            if (tuple((a + b) % 3 for a, b in zip(x, d)) in S and
                    tuple((a + 2 * b) % 3 for a, b in zip(x, d)) in S):
                ok = True
                break
        if not ok:
            return False
    return True


# ---- search state ------------------------------------------------------------------------
class Walker:
    def __init__(self, T, gens, rng):
        self.T = T
        self.rng = rng
        self.orb, self.mem = orbits(T, gens)
        self.smax = max(len(m) for m in self.mem)
        self.reset()

    def reset(self):
        self.inK = bytearray(self.T["P"])
        self.c = [0] * self.T["D"]
        self.size = 0
        self.Kset = set()

    def _add_point(self, q):
        D, LA, LB, inK, c = self.T["D"], self.T["LA"], self.T["LB"], self.inK, self.c
        b = q * D
        for j in range(D):
            if inK[LA[b + j]] and inK[LB[b + j]]:
                c[j] += 1
        inK[q] = 1

    def _remove_point(self, p):
        D, LA, LB, inK, c = self.T["D"], self.T["LA"], self.T["LB"], self.inK, self.c
        inK[p] = 0
        b = p * D
        for j in range(D):
            if inK[LA[b + j]] and inK[LB[b + j]]:
                c[j] -= 1

    def add_orbit(self, o):
        if o in self.Kset:
            return
        for q in self.mem[o]:
            self._add_point(q)
        self.Kset.add(o)
        self.size += len(self.mem[o])

    def remove_orbit(self, o):
        if o not in self.Kset:
            return
        for p in self.mem[o]:
            self._remove_point(p)
        self.Kset.discard(o)
        self.size -= len(self.mem[o])

    def uncovered(self):
        return self.c.count(0)

    def points(self):
        return sorted(p for o in self.Kset for p in self.mem[o])

    def load(self, pts):
        self.reset()
        for o in sorted(set(self.orb[p] for p in pts)):
            self.add_orbit(o)

    def costs(self):
        """cost[o] = number of directions that lose their last line if orbit o is removed."""
        D, LA, LB, inK, c = self.T["D"], self.T["LA"], self.T["LB"], self.inK, self.c
        orb, mem, smax = self.orb, self.mem, self.smax
        Kpts = [p for o in self.Kset for p in mem[o]]
        cost = dict.fromkeys(self.Kset, 0)
        for j in range(D):
            cj = c[j]
            if cj == 0 or cj > smax:
                continue
            found = []
            for p in Kpts:
                b = p * D + j
                a = LA[b]
                e = LB[b]
                if p < a and p < e and inK[a] and inK[e]:
                    found.append((p, a, e))
                    if len(found) == cj:
                        break
            for o in set(orb[x] for x in found[0]):
                if all(any(orb[x] == o for x in ln) for ln in found):
                    cost[o] += 1
        return cost

    def greedy(self):
        """Random greedy Kakeya set: complete a cheapest line of each uncovered direction."""
        rng, orb, mem = self.rng, self.orb, self.mem
        order = list(range(self.T["D"]))
        rng.shuffle(order)
        for j in order:
            if self.c[j]:
                continue
            best, bl = None, []
            for ln in self.T["lines"][j]:
                os_ = set(orb[x] for x in ln if not self.inK[x])
                add = sum(len(mem[o]) for o in os_)
                if best is None or add < best:
                    best, bl = add, [os_]
                elif add == best:
                    bl.append(os_)
            for o in rng.choice(bl):
                self.add_orbit(o)

    def prune(self):
        """Remove orbits whose removal keeps the set Kakeya, in random order, until none."""
        changed = True
        while changed:
            changed = False
            order = list(self.Kset)
            self.rng.shuffle(order)
            for o in order:
                self.remove_orbit(o)
                if self.uncovered():
                    self.add_orbit(o)
                else:
                    changed = True


def run(n, group="trivial", seed=1, budget=60.0, noise=0.1, tenure=5, stall=4000,
        fresh=4, ncand=6, T=None, start=None):
    """Search for small Kakeya sets in F_3^n for `budget` seconds.  Returns a dict with the
    smallest set found (as point indices, row-major base 3), when each size was first
    reached, and counters."""
    t0 = time.time()
    if T is None:
        T = setup(n)
    rng = random.Random(seed)
    W = Walker(T, GROUPS[group], rng)
    orb, mem, D = W.orb, W.mem, T["D"]
    if start is not None:
        W.load(start)
    else:
        W.greedy()
    W.prune()
    best_pts = W.points()
    reached = [(len(best_pts), round(time.time() - t0, 1))]
    traj_best = list(best_pts)
    added_at = {}
    it = 0
    last_imp = 0
    best_u_here = 10 ** 9
    restarts = 0
    greedy_starts = 1 if start is None else 0
    target = W.size - 1

    def shrink(protect):
        while W.size > target:
            cost = W.costs()
            allowed = [o for o in W.Kset if o not in protect and
                       it - added_at.get(o, -10 ** 9) > tenure]
            if not allowed:
                allowed = [o for o in W.Kset if o not in protect] or list(W.Kset)
            if rng.random() < noise:
                o = rng.choice(allowed)
            else:
                mc = min(cost[o] for o in allowed)
                o = rng.choice([o for o in allowed if cost[o] == mc])
            W.remove_orbit(o)

    shrink(set())
    while time.time() - t0 < budget:
        it += 1
        u = W.uncovered()
        if u == 0:
            pts = W.points()
            if len(pts) < len(traj_best):
                traj_best = pts
            if len(pts) < len(best_pts):
                best_pts = pts
                reached.append((len(pts), round(time.time() - t0, 1)))
            target = W.size - 1
            shrink(set())
            last_imp = it
            best_u_here = W.uncovered()
            continue
        if u < best_u_here:
            best_u_here = u
            last_imp = it
        if it - last_imp > stall:
            restarts += 1
            if fresh and restarts % fresh == 0:
                W.reset()
                W.greedy()
                W.prune()
                greedy_starts += 1
                traj_best = W.points()
                if len(traj_best) < len(best_pts):
                    best_pts = list(traj_best)
                    reached.append((len(best_pts), round(time.time() - t0, 1)))
            else:
                W.load(traj_best)
            target = W.size - 1
            added_at = {}
            W.remove_orbit(rng.choice(list(W.Kset)))
            shrink(set())
            last_imp = it
            best_u_here = W.uncovered()
            continue
        unc = [j for j in range(D) if W.c[j] == 0]
        j = rng.choice(unc)
        cands = []
        for ln in T["lines"][j]:
            os_ = frozenset(orb[x] for x in ln if not W.inK[x])
            add = sum(len(mem[o]) for o in os_)
            cands.append((add, os_))
        mn = min(a for a, _ in cands)
        if rng.random() < noise:
            pool = [x for x in cands if x[0] <= mn + 1]
            pick = rng.choice(pool)[1]
        else:
            pool = [x[1] for x in cands if x[0] == mn]
            rng.shuffle(pool)
            pool = pool[:ncand]
            bestu, bl = None, []
            for os_ in pool:
                for o in os_:
                    W.add_orbit(o)
                uu = W.uncovered()
                for o in os_:
                    W.remove_orbit(o)
                if bestu is None or uu < bestu:
                    bestu, bl = uu, [os_]
                elif uu == bestu:
                    bl.append(os_)
            pick = rng.choice(bl)
        for o in pick:
            W.add_orbit(o)
            added_at[o] = it
        shrink(set(pick))
    return {"n": n, "group": group, "seed": seed, "best": len(best_pts), "reached": reached,
            "iterations": it, "restarts": restarts, "greedy_starts": greedy_starts,
            "seconds": round(time.time() - t0, 1), "points": best_pts,
            "verified": is_kakeya(n, [T["pts"][p] for p in best_pts])}


def suite(plan):
    """plan: list of (n, group, seed, budget[, noise]).  Returns the list of run results."""
    tables = {}
    out = []
    for item in plan:
        n, group, seed, budget = item[:4]
        noise = item[4] if len(item) > 4 else 0.1
        if n not in tables:
            tables[n] = setup(n)
        out.append(run(n, group, seed, budget, noise=noise, T=tables[n]))
    return out


def report(results, show_points_at_most=53):
    parts = []
    for r in results:
        s = ("n=%d group=%s seed=%d: best %d (verified %s) reached %s; %d steps, %d restarts, "
             "%d greedy starts, %.0fs" % (r["n"], r["group"], r["seed"], r["best"],
                                         r["verified"], r["reached"], r["iterations"],
                                         r["restarts"], r["greedy_starts"], r["seconds"]))
        if r["n"] == 5 and r["best"] <= show_points_at_most:
            s += " POINTS " + ",".join(str(p) for p in r["points"])
        parts.append(s)
    return " || ".join(parts)
