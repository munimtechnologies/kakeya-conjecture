# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""AAB pairs of F_3^3 (aab/AAB.md): there is no pair with (a, b, e) equal to (7, 7, 1),
(6, 8, 1) or (8, 6, 1).  These are the only size triples allowed by k_3 = 13 whose
products a^t + b^t - e^t would beat the records 105 (F_3^6) or 844 (F_3^9), so products of
F_3^3 pairs cannot improve those records.

Method: every a-point set A of F_3^3 containing 0 (any pair can be translated so) whose
pairs cover all 13 directions, and for each an exact branch-and-bound search for the least
B with one point inside A (padding allowed).  Positive controls check that the same code
finds the classical pair ({0,1}^3, {0,2}^3) and the product lemma on small cases."""
import os
import sys
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

try:
    from . import core
except ImportError:  # pragma: no cover
    import core


class TestControls(unittest.TestCase):
    def test_pair_of_f3_and_products(self):
        self.assertTrue(core.check_pair(1, [0, 1], [0, 2]))
        self.assertFalse(core.check_pair(1, [0, 1], [2]))
        m, A, B = 1, [0, 1], [0, 2]
        for k in range(2, 5):
            m, A, B = core.product(m, A, B, 1, [0, 1], [0, 2])
            self.assertEqual(m, k)
            self.assertTrue(core.check_pair(m, A, B))
            K = set(A) | set(B)
            self.assertEqual(len(K), 2 ** (k + 1) - 1)
            self.assertTrue(core.kakeya_by_pairs(m, K))

    def test_search_finds_the_cube_pair(self):
        S = core.Space(3)
        _, A, B = core.product(2, [0, 1, 3, 4], [0, 2, 6, 8], 1, [0, 1], [0, 2])
        self.assertEqual(sorted(A), [0, 1, 3, 4, 9, 10, 12, 13])
        self.assertTrue(core.check_pair(3, A, B))
        r = core.analyse(S, A, 1, 7)
        self.assertIsNotNone(r)
        f, Bf = r
        self.assertLessEqual(f, 7)
        self.assertTrue(core.check_pair(3, A, Bf))
        self.assertEqual(len(set(A) & set(Bf)), 1)
        self.assertIsNone(core.analyse(S, A, 1, f - 1))

    def test_not_a_pair(self):
        # A = {0, 1, 2} x {0} x {0}: one direction only, so no B works
        S = core.Space(3)
        self.assertIsNone(core.analyse(S, [0, 9, 18], 3, 27))


class TestNoSmallPairInF3Cubed(unittest.TestCase):
    def test_no_7_7_1_no_6_8_1_no_8_6_1(self):
        # (a, share, nshares, Es, limit): limit = b - 1 points of B outside A
        tasks = [(6, 0, 1, (1,), 7), (7, 0, 2, (1,), 6), (7, 1, 2, (1,), 6)]
        tasks += [(8, s, 5, (1,), 5) for s in range(5)]
        res = core.run_pool(core.scan_task, tasks, workers=8)
        seen = {}
        for r in res:
            self.assertTrue(r["done"], r)
            seen[r["a"]] = seen.get(r["a"], 0) + r["seen"]
            hit = r["res"]["1"]
            self.assertIsNone(hit["best"], "pair found: a=%d %r" % (r["a"], hit["example"]))
        # every size has sets A whose pairs cover all 13 directions, so the search ran
        for a in (6, 7, 8):
            self.assertGreater(seen.get(a, 0), 0)


if __name__ == "__main__":
    unittest.main()
