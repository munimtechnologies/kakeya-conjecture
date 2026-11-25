# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Offline checks of k6/fastwalk.py: tables against direct coordinate arithmetic, the
incremental line counts and costs against brute force after random edits, the classical and
product constructions, and the search itself in F_3^3 (k_3 = 13, Lev 2009) and F_3^4
(k_4 = 27, PROOF26.md), where it must reach the known optimum and never go below it."""
import itertools
import random
import unittest

try:
    from . import fastwalk as fw
except ImportError:  # pragma: no cover
    import fastwalk as fw

from k5 import walk as k5walk


def coords(n, p):
    return tuple((p // 3 ** (n - 1 - k)) % 3 for k in range(n))


def index(v):
    s = 0
    for c in v:
        s = 3 * s + c % 3
    return s


def brute_uncovered(T, K):
    S = set(K)
    D, LA, LB = T["D"], T["LA"], T["LB"]
    return {j for j in range(D)
            if not any(x in S and LA[x * D + j] in S and LB[x * D + j] in S for x in S)}


class Tables(unittest.TestCase):
    def test_tables_match_coordinates(self):
        for n in (1, 2, 3, 4):
            T = fw.setup(n)
            self.assertEqual(T["D"], (3 ** n - 1) // 2)
            for p in range(T["P"]):
                v = coords(n, p)
                for j, dj in enumerate(T["dirs"]):
                    d = coords(n, dj)
                    self.assertEqual(T["LA"][p * T["D"] + j], index([a + b for a, b in zip(v, d)]))
                    self.assertEqual(T["LB"][p * T["D"] + j],
                                     index([a + 2 * b for a, b in zip(v, d)]))
            for j in range(T["D"]):
                self.assertEqual(len(T["anchors"][j]), 3 ** (n - 1))
            for dj in T["dirs"]:
                d = coords(n, dj)
                self.assertEqual([c for c in d if c][0], 1)


class Incremental(unittest.TestCase):
    def test_random_edits(self):
        rng = random.Random(7)
        for n in (3, 4):
            T = fw.setup(n)
            frozen = [0, 1] if n == 3 else []
            W = fw.FW(T, frozen)
            for step in range(250):
                p = rng.randrange(T["P"])
                if W.inK[p] and not W.fz[p]:
                    W.remove(p)
                else:
                    W.add(p)
                K = sorted(W.K)
                self.assertEqual(K, [q for q in range(T["P"]) if W.inK[q]])
                unc = brute_uncovered(T, K)
                self.assertEqual(set(W.unc), unc)
                if step % 10 == 0:
                    for q in K:
                        rest = [x for x in K if x != q]
                        self.assertEqual(W.cost[q], len(brute_uncovered(T, rest) - unc),
                                         "n=%d step=%d q=%d" % (n, step, q))
            for q in frozen:
                self.assertTrue(W.inK[q])


class Constructions(unittest.TestCase):
    def test_classical(self):
        for n in range(1, 6):
            K = fw.classical(n)
            self.assertEqual(len(K), 2 ** (n + 1) - 1)
            self.assertTrue(fw.independent_check(n, K))

    def test_plane_and_products(self):
        self.assertTrue(fw.independent_check(2, fw.K7_PLANE))
        self.assertEqual(len(fw.K7_PLANE), 7)
        K = fw.product(fw.classical(2), 2, fw.K7_PLANE, 2)
        self.assertEqual(len(K), 49)
        self.assertTrue(fw.independent_check(4, K))
        self.assertTrue(fw.is_kakeya_idx(fw.setup(4), K))
        K3 = fw.product(fw.full(1), 1, fw.K7_PLANE, 2)
        self.assertTrue(fw.independent_check(3, K3))

    def test_independent_check_rejects(self):
        K = fw.classical(3)
        T = fw.setup(3)
        for p in K:
            rest = [q for q in K if q != p]
            self.assertEqual(fw.independent_check(3, rest), fw.is_kakeya_idx(T, rest))
        self.assertFalse(fw.independent_check(3, list(range(12))))


class Search(unittest.TestCase):
    def test_f3_3_reaches_13(self):
        T = fw.setup(3)
        r = fw.search(T, seed=1, budget=20, stop_at=13)
        self.assertEqual(r["best"], 13)
        self.assertTrue(fw.independent_check(3, r["points"]))

    def test_f3_4_reaches_27(self):
        T = fw.setup(4)
        r = fw.search(T, seed=1, budget=90, stop_at=27)
        self.assertEqual(r["best"], 27, r["reached"])
        self.assertTrue(fw.independent_check(4, r["points"]))

    def test_frozen_kept(self):
        T = fw.setup(3)
        frozen = fw.classical(2)  # inside the layer x0 = 0
        r = fw.search(T, seed=3, budget=5, frozen=frozen)
        self.assertTrue(set(frozen) <= set(r["points"]))
        self.assertTrue(fw.independent_check(3, r["points"]))

    def test_start_from_classical(self):
        T = fw.setup(4)
        r = fw.search(T, seed=2, budget=20, start=fw.classical(4), stop_at=27)
        self.assertLessEqual(r["best"], 31)
        self.assertGreaterEqual(r["best"], 27)
        self.assertTrue(fw.independent_check(4, r["points"]))


if __name__ == "__main__":
    unittest.main()
