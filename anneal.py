"""Simulated annealing on the witness function of a rainbow triple (see TRIPLES.md).

A triple is rainbow exactly when some function f on F_3^m has f(d) in A, f(d) + d in M and
f(d) + 2d in C for every d; the smallest such sets are the three images. This searches f directly.
The cost is how far the three image sizes are above the target caps; a move changes one value
f(d), usually to a value that lands f(d), f(d) + d or f(d) + 2d on a point already in use.
Restarts cycle through the product G1^m (images of 2^m points each), a random function and the
best function so far. Every triple found is re-checked with triples.is_rainbow.

    python anneal.py --out-dir /out --plan 3:7,7,7:600 3:6,7,7:300 4:14,15,15:300

Each plan item is dim:caps:seconds. Results: out/anneal.json, and lines in out/anneal.txt.

By MV2 of Munim, Inc. Made with AI agents working in MV2's name. Apache-2.0.
"""

import argparse
import json
import math
import os
import random
import time

from triples import is_rainbow, sizes, vectors

STEPS_PER_RESTART = 400000


def tables(m):
    pts = vectors(m)
    idx = {p: i for i, p in enumerate(pts)}
    ADD = [[idx[tuple((x + y) % 3 for x, y in zip(p, q))] for q in pts] for p in pts]
    TWICE = [idx[tuple((2 * x) % 3 for x in p)] for p in pts]
    return pts, idx, ADD, TWICE


def g1_power(pts, idx):
    g = (0, 0, 1)
    return [idx[tuple(g[x] for x in p)] for p in pts]


def image_counts(f, ADD, TWICE):
    N = len(f)
    cnt = [[0] * N for _ in range(3)]
    for d, u in enumerate(f):
        cnt[0][u] += 1
        cnt[1][ADD[u][d]] += 1
        cnt[2][ADD[u][TWICE[d]]] += 1
    return cnt


def excess(size, caps):
    return sum(max(0, s - c) for s, c in zip(size, caps))


def run_one(f, caps, ADD, TWICE, rng, steps, t_hot=1.5, t_cold=0.05):
    """One annealing run from f. Returns (best cost, best f)."""
    N = len(f)
    f = list(f)
    A, M, C = image_counts(f, ADD, TWICE)
    size = [sum(1 for c in row if c) for row in (A, M, C)]
    c0, c1, c2 = caps
    cost = excess(size, caps)
    best_cost, best_f = cost, list(f)
    T = t_hot
    for step in range(steps):
        if step % 2000 == 0:
            T = t_hot * (t_cold / t_hot) ** (step / steps)
        d = rng.randrange(N)
        e = rng.randrange(N)
        r = rng.random()
        if r < 0.3:
            v = f[e]
        elif r < 0.6:
            v = ADD[ADD[f[e]][e]][TWICE[d]]        # v + d = f(e) + e
        elif r < 0.9:
            v = ADD[ADD[f[e]][TWICE[e]]][d]        # v + 2d = f(e) + 2e
        else:
            v = rng.randrange(N)
        u = f[d]
        if v == u:
            continue
        um, uc = ADD[u][d], ADD[u][TWICE[d]]
        vm, vc = ADD[v][d], ADD[v][TWICE[d]]
        s0 = size[0] - (A[u] == 1) + (A[v] == 0)
        s1 = size[1] - (M[um] == 1) + (M[vm] == 0)
        s2 = size[2] - (C[uc] == 1) + (C[vc] == 0)
        new = max(0, s0 - c0) + max(0, s1 - c1) + max(0, s2 - c2)
        delta = new - cost
        if delta <= 0 or rng.random() < math.exp(-delta / T):
            A[u] -= 1
            A[v] += 1
            M[um] -= 1
            M[vm] += 1
            C[uc] -= 1
            C[vc] += 1
            size = [s0, s1, s2]
            f[d] = v
            cost = new
            if cost < best_cost:
                best_cost, best_f = cost, list(f)
                if cost == 0:
                    break
    return best_cost, best_f


def search(m, caps, seconds, rng):
    pts, idx, ADD, TWICE = tables(m)
    N = len(pts)
    deadline = time.time() + seconds
    best_cost, best_f, restarts = None, None, 0
    while time.time() < deadline:
        kind = restarts % 3
        if kind == 0 or best_f is None:
            start = g1_power(pts, idx)
        elif kind == 1:
            start = [rng.randrange(N) for _ in range(N)]
        else:
            start = best_f
        restarts += 1
        cost, f = run_one(start, caps, ADD, TWICE, rng, STEPS_PER_RESTART)
        if best_cost is None or cost < best_cost:
            best_cost, best_f = cost, f
        if best_cost == 0:
            break
    T = ({pts[u] for u in best_f},
         {pts[ADD[u][d]] for d, u in enumerate(best_f)},
         {pts[ADD[u][TWICE[d]]] for d, u in enumerate(best_f)})
    rec = {"dim": m, "caps": list(caps), "seconds": seconds, "restarts": restarts,
           "best_cost": best_cost, "best_sizes": list(sizes(T)), "found": best_cost == 0}
    if best_cost == 0:
        rec["verified_rainbow"] = is_rainbow(*T, m)
        rec["triple"] = [sorted(list(p) for p in S) for S in T]
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--plan", nargs="+", default=["3:7,7,7:600", "3:6,7,7:300"])
    ap.add_argument("--seed", type=int, default=20261014)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rng = random.Random(a.seed)
    results = []
    for item in a.plan:
        dim, caps, sec = item.split(":")
        caps = tuple(int(x) for x in caps.split(","))
        rec = search(int(dim), caps, float(sec), rng)
        P = caps[0] * caps[1] * caps[2]
        rec["base_if_found"] = round(P ** (1.0 / (3 * int(dim))), 5)
        results.append(rec)
        line = ("anneal F_3^%s caps %s (base %.4f): %s, best cost %s, best sizes %s, %d restarts%s"
                % (dim, caps, rec["base_if_found"], "FOUND" if rec["found"] else "not found",
                   rec["best_cost"], tuple(rec["best_sizes"]), rec["restarts"],
                   ", verified %s" % rec["verified_rainbow"] if rec["found"] else ""))
        print(line, flush=True)
        with open(os.path.join(a.out_dir, "anneal.txt"), "a") as f:
            f.write(line + "\n")
        with open(os.path.join(a.out_dir, "anneal.json"), "w") as f:
            json.dump(results, f, indent=1)


if __name__ == "__main__":
    main()
