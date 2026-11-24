# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Checks of SD2.md: with the refined sum-difference bound, a 44-point Kakeya set in F_3^5
survives the moment LP only through hyperplane cuts of profile (10, 11, 23)."""
import unittest
from math import comb

try:
    from . import sd2
except ImportError:  # run from inside the folder
    import sd2


class TestEmax(unittest.TestCase):
    def test_values_checked_by_hand(self):
        self.assertEqual(sd2.emax(10, 10, 81), 282)   # columns 10, 10, 2, then 78 ones
        self.assertEqual(sd2.emax(10, 11, 81), 386)   # 10, 10, 10, 3, then 77 ones
        self.assertEqual(sd2.emax(9, 11, 81), 249)
        self.assertEqual(sd2.emax(9, 12, 81), 336)
        self.assertEqual(sd2.emax(11, 11, 81), 561)
        self.assertEqual(sd2.emax(9, 9, 81), 81)
        self.assertIsNone(sd2.emax(8, 10, 81))        # 80 < 81: (T) fails

    def test_emax_is_a_true_maximum_small_cases(self):
        # brute force over all column vectors for tiny parameters
        import itertools
        for T in (3, 4, 5):
            for x in range(1, 5):
                for y in range(1, 5):
                    m = min(x, y)
                    best = None
                    for r in itertools.product(range(1, m + 1), repeat=T):
                        if sum(r) == x * y:
                            v = sum(t * t for t in r)
                            best = v if best is None or v > best else best
                    self.assertEqual(sd2.emax(x, y, T), best, (x, y, T))

    def test_min_third(self):
        self.assertEqual(sd2.min_third(10, 10, 81, refined=False), 24)  # 6561/282 = 23.3
        self.assertEqual(sd2.min_third(10, 10, 81), 25)                 # 6561/263 = 24.9
        self.assertEqual(sd2.min_third(10, 11, 81), 19)                 # 6561/357 = 18.4
        self.assertEqual(sd2.min_third(9, 11, 81), 29)
        self.assertEqual(sd2.min_third(9, 12, 81), 22)


class TestN44(unittest.TestCase):
    def setUp(self):
        self.T, self.H, self.a1, self.a2 = sd2.constants(5)
        self.base = (14, 15, 15)

    def test_constants(self):
        self.assertEqual((self.T, self.H, self.a1, self.a2), (81, 121, 40, 13))
        S2 = self.a1 * comb(44, 2)
        S3 = self.a2 * comb(44, 3) + (self.a1 - self.a2) * self.H
        self.assertEqual(S2 - self.H * sd2.s2(self.base), 1419)
        self.assertEqual(S3 - self.H * sd2.s3(self.base), 21285)
        self.assertEqual(21285, 15 * 1419)

    def test_profile_lists(self):
        P1 = sd2.profiles(44, self.T, "sd1")
        P2 = sd2.profiles(44, self.T, "sd2")
        self.assertIn((10, 10, 24), P1)
        self.assertNotIn((10, 10, 24), P2)
        for p in (self.base, (10, 11, 23), (11, 11, 22), (9, 12, 23)):
            self.assertIn(p, P2)
        for p in ((9, 11, 24), (9, 10, 25), (8, 11, 25), (9, 9, 26), (8, 13, 23)):
            self.assertNotIn(p, P2)
        self.assertTrue(set(P2) <= set(P1))

    def test_base_is_the_unique_smallest_s2(self):
        P2 = sd2.profiles(44, self.T, "sd2")
        b = sd2.s2(self.base)
        self.assertTrue(all(sd2.s2(p) > b for p in P2 if p != self.base))

    def test_only_10_11_23_reaches_ratio_15(self):
        P2 = sd2.profiles(44, self.T, "sd2")
        b2, b3 = sd2.s2(self.base), sd2.s3(self.base)
        threats = [p for p in P2 if p != self.base
                   and sd2.s3(p) - b3 >= 15 * (sd2.s2(p) - b2)]
        self.assertEqual(threats, [(10, 11, 23)])

    def test_lp(self):
        self.assertTrue(sd2.lp_feasible(5, 44, "sd1"))
        self.assertTrue(sd2.lp_feasible(5, 44, "sd2"))
        self.assertFalse(sd2.lp_feasible(5, 44, "sd2", exclude=((10, 11, 23),)))
        # with SD alone, removing (10,11,23) is not enough: (10,10,24) still carries the LP
        self.assertTrue(sd2.lp_feasible(5, 44, "sd1", exclude=((10, 11, 23),)))

    def test_at_least_20_cuts(self):
        # k cuts of profile (10,11,23); every other allowed profile has ratio <= 14.9
        P2 = sd2.profiles(44, self.T, "sd2")
        b2, b3 = sd2.s2(self.base), sd2.s3(self.base)
        others = [p for p in P2 if p not in (self.base, (10, 11, 23))]
        self.assertTrue(all(10 * (sd2.s3(p) - b3) <= 149 * (sd2.s2(p) - b2) for p in others))
        d2 = sd2.s2((10, 11, 23)) - b2
        d3 = sd2.s3((10, 11, 23)) - b3
        self.assertEqual((d2, d3), (52, 782))
        # 782 k + 14.9 (1419 - 52 k) >= 21285  <=>  72 k >= 1419
        ks = [k for k in range(0, 28) if 10 * 782 * k + 149 * (1419 - 52 * k) >= 10 * 21285]
        self.assertEqual(min(ks), 20)


if __name__ == "__main__":
    unittest.main()
