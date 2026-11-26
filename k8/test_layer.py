# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Unit tests for k8/layer.py: tables against tuple arithmetic, the incremental counters
against a recount from scratch, the layer criterion against the stand-alone check of
k8/check.py, and small searches whose results are re-checked by k8/check.py."""
import itertools
import random
import unittest

try:
    from . import check, layer
except ImportError:  # pragma: no cover
    import check
    import layer


def tup(m, i):
    return tuple((i // 3 ** (m - 1 - k)) % 3 for k in range(m))


def idx(t):
    z = 0
    for a in t:
        z = 3 * z + a
    return z


class LayerTests(unittest.TestCase):
    def test_tables(self):
        m = 3
        L0 = layer.classical(m)
        W = layer.Layer(m, L0)
        self.assertEqual(W.L, 15)
        for v in range(W.N):
            tv = tup(m, v)
            for i, x in enumerate(W.X):
                tx = tup(m, x)
                y = idx(tuple((a + b) % 3 for a, b in zip(tx, tv)))
                z = idx(tuple((a + 2 * b) % 3 for a, b in zip(tx, tv)))
                self.assertEqual(W.YV[v * W.L + i], y)
                self.assertEqual(W.ZV[v * W.L + i], z)
                self.assertEqual(W.V1[y * W.L + i], v)
                self.assertEqual(W.V2[z * W.L + i], v)

    def recount(self, W):
        N, L = W.N, W.L
        cnt, cost = [0] * N, [0] * (2 * N)
        for v in range(N):
            ls = [i for i in range(L)
                  if W.inS[W.YV[v * L + i]] and W.inS[N + W.ZV[v * L + i]]]
            cnt[v] = len(ls)
            if len(ls) == 1:
                cost[W.YV[v * L + ls[0]]] += 1
                cost[N + W.ZV[v * L + ls[0]]] += 1
        return cnt, cost

    def test_incremental_state(self):
        m = 3
        W = layer.Layer(m, layer.classical(m))
        rng = random.Random(7)
        for step in range(600):
            p = rng.randrange(2 * W.N)
            if W.inS[p] and rng.random() < 0.45:
                W.remove(p)
            else:
                W.add(p)
            if step % 50 == 49:
                cnt, cost = self.recount(W)
                self.assertEqual(W.cnt, cnt)
                self.assertEqual(W.cost, cost)
                self.assertEqual(sorted(W.unc), [v for v in range(W.N) if cnt[v] == 0])
                self.assertEqual(sorted(W.S), [p for p in range(2 * W.N) if W.inS[p]])

    def test_layer_criterion_matches_direct_check(self):
        m = 2
        W = layer.Layer(m, layer.K7_PLANE)
        self.assertTrue(check.is_kakeya(2, layer.K7_PLANE))
        rng = random.Random(3)
        seen = {True: 0, False: 0}
        for trial in range(300):
            k = rng.randrange(4, 12)
            S = rng.sample(range(2 * W.N), k)
            W.load(S)
            ok = check.is_kakeya(3, W.full_indices())
            self.assertEqual(ok, not W.unc, (S,))
            seen[ok] += 1
        # make sure both outcomes were exercised: complete the last set greedily
        W.greedy(random.Random(1))
        self.assertTrue(check.is_kakeya(3, W.full_indices()))
        self.assertGreater(seen[False], 0)

    def test_small_searches(self):
        # F_3^3 from the 7-point plane set, F_3^4 from the 15-point classical set of F_3^3
        for m, L0, cap in ((2, layer.K7_PLANE, 15), (3, layer.classical(3), 31)):
            W = layer.Layer(m, L0)
            r = layer.search(W, seed=1, budget=40.0, max_steps=4000)
            full = W.full_indices(r["points"])
            self.assertTrue(check.is_kakeya(m + 1, full), m)
            self.assertEqual(len(full), len(L0) + r["best"])
            self.assertLessEqual(len(full), cap)
            print("layer search m=%d |L0|=%d extra=%d total=%d steps=%d" % (
                m, len(L0), r["best"], len(full), r["steps"]))

    def test_check_module(self):
        self.assertEqual(sum(1 for _ in check.directions(8)), 3280)
        C = layer.classical(4)
        self.assertTrue(check.is_kakeya(4, C))
        self.assertFalse(check.is_kakeya(4, [c for c in C if c != 0]))


if __name__ == "__main__":
    unittest.main()
