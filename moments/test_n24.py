# Tests for N24.md. MV2 of Munim, Inc. Made with AI agents. Apache-2.0. Standard library only.
import os
import sys
import unittest
from math import comb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ag33_lines as ag  # noqa: E402
import kakeya_moments as km  # noqa: E402
import n24  # noqa: E402


class TestGeometry(unittest.TestCase):
    def test_counts(self):
        pts, _, dirs, lines = ag.geometry()
        self.assertEqual(len(pts), 27)
        self.assertEqual(len(dirs), 13)
        self.assertEqual([len(f) for f in lines], [9] * 13)
        for fam in lines:
            # the 9 parallel lines of one direction partition the 27 points
            u = 0
            for m in fam:
                self.assertEqual(bin(m).count("1"), 3)
                self.assertEqual(u & m, 0)
                u |= m
            self.assertEqual(u, (1 << 27) - 1)


class TestPartialKakeyaAG33(unittest.TestCase):
    def test_small_unions(self):
        # f(k) = fewest points carrying lines in k distinct directions: f(1..5) = 3, 5, 6, 7, 9
        for k, f in [(1, 3), (2, 5), (3, 6), (4, 7), (5, 9)]:
            self.assertTrue(ag.fits(k, f), k)
            self.assertFalse(ag.fits(k, f - 1), k)

    def test_d_upper_is_sharp_up_to_8(self):
        for m in range(3, 9):
            self.assertEqual(ag.max_directions(m), n24.D_UPPER[m], m)

    def test_paraboloid_has_six_directions(self):
        # the 9 points (x, y, xy) carry lines in exactly 6 directions, so D(9) >= 6 and f(6) = 9
        mask = ag.mask_of([(x, y, x * y) for x in range(3) for y in range(3)])
        self.assertEqual(bin(mask).count("1"), 9)
        self.assertEqual(ag.directions_inside(mask), 6)
        self.assertTrue(ag.fits(6, 9))


class TestN24(unittest.TestCase):
    def test_profiles_and_coverage(self):
        all24 = km.profiles(3, 4, 24)
        self.assertEqual(len(all24), 17)
        self.assertEqual([P for P in all24 if not n24.covers(P)], [(8, 8, 8)])

    def test_every_class_has_a_section_of_9(self):
        # three sections of at most 8 points carry at most 4 + 4 + 4 = 12 < 13 directions
        for a in range(1, 9):
            for b in range(a, 9):
                for c in range(b, 9):
                    self.assertFalse(n24.covers((a, b, c)))

    def test_line_bounds(self):
        all24 = km.profiles(3, 4, 24)
        ok24 = [P for P in all24 if n24.covers(P)]
        self.assertEqual(n24.max_lines(24, all24), 44)
        self.assertEqual(n24.max_lines(24, ok24), 41)

    def test_deficits(self):
        table = {(7, 8, 9): 0, (5, 6, 13): 0, (7, 7, 10): 12, (6, 6, 12): 21, (6, 7, 11): 36,
                 (6, 8, 10): 45, (6, 9, 9): 48, (5, 7, 12): 72, (4, 7, 13): 120, (5, 8, 11): 120,
                 (5, 9, 10): 144, (4, 8, 12): 225, (4, 9, 11): 288, (4, 10, 10): 309,
                 (3, 9, 12): 480, (3, 10, 11): 540, (8, 8, 8): -15}
        self.assertEqual({P: n24.deficit(P) for P in km.profiles(3, 4, 24)}, table)
        # sum over the 40 classes of the deficits = 6780 - 162 L (L = lines of K)
        for L in (40, 41, 42):
            total = -8835 * 40 - (18 * (4 * comb(24, 3) + 9 * L) - 141 * 13 * comb(24, 2))
            self.assertEqual(total, 6780 - 162 * L)
        self.assertEqual(6780 - 162 * 40, 300)
        self.assertEqual(6780 - 162 * 41, 138)

    def test_identity_arithmetic(self):
        # L = 40, N = 24: sum_H n_H lambda_H = 4920, sum_H lambda_H (C(n_H,2) - lambda_H) = 18428
        N, sr = 24, 120
        self.assertEqual((40 + 13 * (N - 1)) * sr, 40680)
        self.assertEqual((40680 - 27 * 40 * N) % 3, 0)
        self.assertEqual((40680 - 27 * 40 * N) // 3, 4920)
        per_point = 13 * (N - 1) + 4 * comb(N - 1, 2)  # plus 9 r_p
        self.assertEqual(per_point, 1311)
        self.assertEqual((per_point * sr - 27 * 13 * comb(N, 2)) % 3, 0)
        self.assertEqual((per_point * sr - 27 * 13 * comb(N, 2)) // 3, 20148)
        # sum_H (27 + 3 lambda_H)^2 = 27 sum r^2 + 13 (sum r)^2 with sum_H lambda_H = 13 * 40
        const = 120 * 27 ** 2 + 2 * 27 * 3 * 13 * 40
        self.assertEqual(const, 171720)
        self.assertEqual((13 * sr ** 2 - const) % 9, 0)
        self.assertEqual((13 * sr ** 2 - const) // 9, 1720)
        self.assertEqual(20148 - 1720, 18428)


if __name__ == "__main__":
    unittest.main()
