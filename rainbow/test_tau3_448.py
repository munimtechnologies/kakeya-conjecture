# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""tau(3) = 448 (rainbow/TAU3.md): every rainbow triple (A, M, C) of F_3^3 has
|A| |M| |C| >= 448, and the (7, 8, 8) triple found by CP-SAT in October 2026 attains it.
In particular there is no rainbow triple of F_3^3 with three sets of at most 7 points, and
the base of dimension 3 is 448^(1/9) = 1.9705...

Method (SCAN.md): normal form as in test_tau3_308.py.  A product below 448 forces |A| <= 7
(8^3 = 512) and, by lemma H (9 tau(2) = 576 > 448), A spanning F_3^3, so |A| is 4 to 7.
Every class, every allowed size of M and every M in normal form is scanned, with the least
C found exactly; none gives a product below 448."""
import math
import os
import sys
import time
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

try:
    from . import rscan
except ImportError:  # pragma: no cover
    import rscan

A788 = [(0, 0, 0), (0, 0, 1), (0, 2, 1), (1, 0, 0), (1, 2, 1), (2, 0, 1), (2, 2, 2)]
M788 = [(0, 0, 0), (0, 1, 2), (0, 2, 0), (1, 0, 2), (1, 1, 2), (2, 0, 0), (2, 0, 2),
        (2, 1, 0)]
C788 = [(0, 0, 0), (0, 0, 2), (0, 1, 0), (0, 1, 2), (0, 2, 0), (1, 0, 2), (2, 0, 0),
        (2, 1, 2)]


def ix(pts):
    return sorted(rscan.to_idx(p) for p in pts)


class TestTau3Is448(unittest.TestCase):
    def test_a_upper_bound_448(self):
        A, M, C = ix(A788), ix(M788), ix(C788)
        self.assertEqual(len(A) * len(M) * len(C), 448)
        self.assertTrue(rscan.check_triple(3, A, M, C))

    def test_b_class_list_of_7_sets_is_complete(self):
        cl = rscan.classes(3, 7)
        full = [st for k, _, st in cl if k == 3]
        plan = [st for k, _, st in cl if k == 2]
        self.assertEqual(len(full) + len(plan), len(cl))
        self.assertEqual(sum(303264 // st for st in full),
                         math.comb(27, 7) - 39 * math.comb(9, 7))
        self.assertEqual(sum(432 // st for st in plan), math.comb(9, 7))

    def test_c_no_rainbow_triple_of_f3_cubed_below_448(self):
        self.assertEqual(rscan.sizes_for(3, 4, 448), [7, 8, 9, 10])
        self.assertEqual(rscan.sizes_for(3, 5, 448), [6, 7, 8, 9])
        self.assertEqual(rscan.sizes_for(3, 6, 448), [6, 7, 8])
        self.assertEqual(rscan.sizes_for(3, 7, 448), [7])
        self.assertEqual(rscan.sizes_for(3, 8, 448), [])
        deadline = time.time() + 420
        res = rscan.sweep(3, 448, (4, 5, 6, 7), deadline, workers=8, full_only=True)
        by = rscan.summarise(res)
        msg = rscan.describe(by)
        hits = [r for r in res if r["best"] is not None]
        self.assertFalse(hits, "TRIPLE BELOW 448 " + repr(
            [(r["best"], r["A"], r["M"], r["C"]) for r in hits[:3]]) + " || " + msg)
        self.assertTrue(all(r["complete"] for r in res), "INCOMPLETE || " + msg)
        self.assertEqual(sorted(by), [4, 5, 6, 7], msg)
        for x in (4, 5, 6, 7):
            nfull = sum(1 for k, _, _ in rscan.classes(3, x) if k == 3)
            self.assertEqual(by[x]["reps"], nfull, msg)
            self.assertGreater(by[x]["nodes"], 0, msg)


if __name__ == "__main__":
    unittest.main()
