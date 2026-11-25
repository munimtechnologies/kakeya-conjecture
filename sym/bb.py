# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Exact branch and bound for the orbit covering problem (independent of the CP-SAT solver).

Input: weights w[o] of the orbits, and for each representative direction a list of bit masks
(sets of orbits whose union contains a line of that direction).  Find a set M of orbits of
least total weight such that every direction has a mask inside M.

Branching: take the uncovered direction with the fewest options that still fit under the
current upper bound; branch over its masks.  Bound: current cost + the largest, over the
uncovered directions, of the cheapest way to cover that direction alone.  Visited sets are
memoised.  `improve(..., ub)` searches for a set strictly cheaper than ub; when the search
completes without one, ub is proved optimal for this orbit problem.
"""
import time


def _cost(w, L):
    c = 0
    while L:
        low = L & -L
        c += w[low.bit_length() - 1]
        L ^= low
    return c


def improve(w, reps, ub, seconds=60.0, memo_cap=3_000_000):
    t0 = time.time()
    st = {"ub": ub, "best": None, "nodes": 0, "timeout": False}
    seen = set()

    def rec(M, c):
        if st["timeout"]:
            return
        st["nodes"] += 1
        if st["nodes"] % 4096 == 0 and time.time() - t0 > seconds:
            st["timeout"] = True
            return
        if M in seen:
            return
        if len(seen) < memo_cap:
            seen.add(M)
        ubc = st["ub"]
        pick, pick_n, lb = None, None, 0
        for opts in reps:
            covered = False
            fit = []
            mn = None
            for L in opts:
                rest = L & ~M
                if rest == 0:
                    covered = True
                    break
                a = _cost(w, rest)
                if mn is None or a < mn:
                    mn = a
                if c + a < ubc:
                    fit.append((a, rest))
            if covered:
                continue
            if mn is None or c + mn >= ubc:
                return
            if mn > lb:
                lb = mn
            if pick is None or len(fit) < pick_n:
                pick, pick_n = fit, len(fit)
        if pick is None:
            st["ub"] = c
            st["best"] = M
            return
        if c + lb >= ubc:
            return
        pick.sort()
        done = set()
        for a, rest in pick:
            if rest in done:
                continue
            done.add(rest)
            if c + a < st["ub"]:
                rec(M | rest, c + a)

    rec(0, 0)
    return {"found": st["best"], "value": st["ub"] if st["best"] is not None else None,
            "complete": not st["timeout"], "nodes": st["nodes"],
            "seconds": round(time.time() - t0, 1)}


def minimum(w, reps, seconds=60.0):
    """Exact minimum (or the best found before the time limit)."""
    ub = sum(w) + 1
    r = improve(w, reps, ub, seconds)
    return r
