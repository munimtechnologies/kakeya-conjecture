# test_sd_bound.py -- checks for sd_bound.py and SD.md (k_4 >= 25, k_5 >= 44).
# MV2 of Munim, Inc. Made with AI agents working in MV2's name. No human mathematician has checked
# this. Licence: Apache-2.0. Standard library only.

import unittest
from math import comb

from moments import sd_bound as S


class EnergyCap(unittest.TestCase):
    def test_hand_values(self):
        # values used by hand in SD.md
        self.assertEqual(S.energy_cap(5, 6, 27), 42)
        self.assertEqual(S.energy_cap(4, 7, 27), 30)
        self.assertEqual(S.energy_cap(3, 10, 27), 38)
        self.assertEqual(S.energy_cap(4, 8, 27), 50)
        self.assertEqual(S.energy_cap(3, 9, 27), 27)
        self.assertEqual(S.energy_cap(9, 9, 81), 81)
        self.assertEqual(S.energy_cap(9, 10, 81), 164)
        self.assertEqual(S.energy_cap(9, 11, 81), 249)
        self.assertEqual(S.energy_cap(10, 10, 81), 282)
        self.assertEqual(S.energy_cap(3, 4, 9), 20)
        self.assertEqual(S.sd_floor(5, 6, 27), 18)
        self.assertEqual(S.sd_floor(10, 10, 81), 24)
        self.assertEqual(S.sd_floor(9, 11, 81), 27)
        self.assertEqual(S.sd_floor(9, 10, 81), 41)
        self.assertIsNone(S.energy_cap(2, 4, 9))

    def test_cap_is_a_true_maximum(self):
        # brute force over all admissible r-vectors for small cases (T = 9 and T = 3)
        def brute(x, y, T):
            m, total = min(x, y), x * y
            best = None

            # nonincreasing r-vectors only (the sum of squares does not depend on the order)
            def rec(i, left, acc, top):
                nonlocal best
                if i == T:
                    if left == 0 and (best is None or acc > best):
                        best = acc
                    return
                rest = T - i - 1
                for r in range(min(m, top), 0, -1):
                    if rest <= left - r <= rest * r:
                        rec(i + 1, left - r, acc + r * r, r)

            rec(0, total, 0, m)
            return best

        for T in (3, 9):
            for x in range(1, T + 1):
                for y in range(x, T + 1):
                    if x * y >= T:
                        self.assertEqual(S.energy_cap(x, y, T), brute(x, y, T), (x, y, T))


class AG33(unittest.TestCase):
    def test_f_values(self):
        f, table = S.d3_table(5)
        self.assertEqual(f, [3, 5, 6, 7, 9])
        self.assertEqual([table[m] for m in range(3, 9)], [1, 1, 2, 3, 4, 4])
        self.assertEqual(S.D3(8), 4)
        self.assertEqual(S.D3(9), 13)


class KnownValues(unittest.TestCase):
    """The method must reproduce the known exact values k_2 = 7 and k_3 = 13 and must not exclude
    sizes that are known to occur (27 in F_3^4, 53 in F_3^5)."""

    def test_k2_k3(self):
        self.assertEqual(S.threshold(2), 7)
        self.assertEqual(S.threshold(3), 13)

    def test_known_sets_not_excluded(self):
        self.assertTrue(S.feasible(4, 27, sd=True, hz=True))
        self.assertTrue(S.feasible(5, 53, sd=True))

    def test_old_bounds_reproduced(self):
        # moments + (T) alone: k_4 >= 24 (moments/README.md); for n = 5 that method already gives 42
        self.assertEqual(S.threshold(4, sd=False), 24)
        self.assertEqual(S.threshold(5, sd=False), 42)


class NewBounds(unittest.TestCase):
    def test_k4_at_least_25(self):
        self.assertEqual(S.threshold(4, sd=True, hz=True), 25)
        # both new conditions are needed at N = 24
        self.assertTrue(S.feasible(4, 24, sd=True, hz=False))
        self.assertTrue(S.feasible(4, 24, sd=False, hz=True))
        self.assertFalse(S.feasible(4, 24, sd=True, hz=True))

    def test_k4_hand_certificate(self):
        # 11 s3 - 85 s2 <= -5300 on every allowed profile of N = 24, tight at (7,8,9) and (6,6,12)
        self.assertEqual(S.certificate_max(4, 24, 11, 85, sd=True, hz=True), -5300)
        self.assertTrue(S.certificate_refutes(4, 24, 11, 85, sd=True, hz=True))
        allowed = S.profiles(4, 24, sd=True, hz=True)
        self.assertIn((7, 8, 9), allowed)
        self.assertIn((6, 6, 12), allowed)
        for p in [(8, 8, 8), (5, 6, 13), (4, 7, 13), (4, 8, 12), (3, 9, 12), (3, 10, 11)]:
            self.assertNotIn(p, allowed)
        self.assertEqual(len(allowed), 11)
        # sum: 11*(4*C(24,3) + 9L) - 85*13*C(24,2) <= 40*(-5300)  <=>  99 L <= 3924, but L >= 40
        self.assertEqual(85 * 13 * comb(24, 2) - 11 * 4 * comb(24, 3) - 40 * 5300, 3924)
        self.assertGreater(99 * 40, 3924)

    def test_k5_at_least_44(self):
        self.assertEqual(S.threshold(5, sd=True), 44)

    def test_k5_hand_certificates(self):
        self.assertEqual(S.certificate_max(5, 41, 13, 185, sd=False), -34918)
        self.assertTrue(S.certificate_refutes(5, 41, 13, 185, sd=False))
        self.assertEqual(S.certificate_max(5, 42, 3, 43, sd=True), -8463)
        self.assertTrue(S.certificate_refutes(5, 42, 3, 43, sd=True))
        self.assertEqual(S.certificate_max(5, 43, 22, 321, sd=True), -66101)
        self.assertTrue(S.certificate_refutes(5, 43, 22, 321, sd=True))
        a43 = S.profiles(5, 43, sd=True)
        self.assertIn((14, 14, 15), a43)
        self.assertIn((10, 11, 22), a43)
        for p in [(9, 9, 25), (9, 10, 24), (9, 11, 23), (10, 10, 23), (8, 11, 24)]:
            self.assertNotIn(p, a43)
        self.assertIn((10, 10, 24), S.profiles(5, 44, sd=True))


if __name__ == "__main__":
    unittest.main()
