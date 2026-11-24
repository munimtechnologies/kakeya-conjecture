# test_plane_lp.py -- checks for plane_lp.py (PLANE.md): the identities and per-array conditions
# hold on a real 27-point Kakeya set, the enumeration finds its arrays, the exact simplex works,
# and the relaxation behaves as it must at N = 24 (impossible) and N = 27 (a set exists).
# MV2 of Munim, Inc. Made with AI agents working in MV2's name. No human mathematician has checked
# this. Licence: Apache-2.0. Standard library only.

import unittest
from math import comb

from moments import plane_lp as P


class Simplex(unittest.TestCase):
    def test_small_systems(self):
        # x1 + x2 = 2, x1 + 2 x2 = 3: x = (1, 1)
        self.assertTrue(P.phase1_feasible([(1, 1), (1, 2)], [2, 3]))
        # x1 + x2 = 1 and x1 + x2 = 2: impossible
        self.assertFalse(P.phase1_feasible([(1, 1), (1, 1)], [1, 2]))
        # x1 - x2 = -1 needs x2 = x1 + 1; with x1 + x2 = 1: x = (0, 1)
        self.assertTrue(P.phase1_feasible([(1, 1), (1, -1)], [1, -1]))
        # 2 x1 = 1 has the rational solution 1/2
        self.assertTrue(P.phase1_feasible([(2,)], [1]))
        # x1 + x2 + x3 = 1, x1 - x3 = 2: needs x1 >= 2 > 1
        self.assertFalse(P.phase1_feasible([(1, 1), (1, 0), (1, -1)], [1, 2]))


class Geometry(unittest.TestCase):
    def test_lines_of_ag23(self):
        lines = [l for fam in P.LINES2 for l in fam]
        self.assertEqual(len(set(lines)), 12)
        for fam in P.LINES2:
            self.assertEqual(sorted(c for l in fam for c in l), list(range(9)))
        # any two cells lie on exactly one line
        for a in range(9):
            for b in range(a + 1, 9):
                self.assertEqual(sum(1 for l in lines if a in l and b in l), 1)

    def test_plane_classes(self):
        pc = P.plane_classes()
        self.assertEqual(len(pc), 130)

    def test_small_tables(self):
        # 9 pairs with distinct differences have 9 distinct sums, but a 3-point third section
        # takes only 3 of them (the sum-difference condition with T = 9)
        self.assertEqual(P.cross_cap(3, 3, 3), 8)
        self.assertEqual(P.cross_cap(3, 3, 9), 9)
        self.assertEqual(P.cross_cap(1, 9, 9), 9)
        self.assertEqual(P.cross_cap(2, 4, 9), 8)
        self.assertEqual(P.cross_cap(2, 2, 9), 4)


class RealSet(unittest.TestCase):
    """The identities and every per-array condition must hold for a real Kakeya set."""

    @classmethod
    def setUpClass(cls):
        cls.K = P.K27
        cls.L = P.count_lines(cls.K)
        cls.arrays = [P.array_of(cls.K, u, v) for u, v in P.plane_classes()]

    def test_is_kakeya(self):
        self.assertEqual(len(set(self.K)), 27)
        self.assertTrue(P.is_kakeya(self.K))
        self.assertGreaterEqual(self.L, 40)

    def test_arrays_valid(self):
        for M in self.arrays:
            self.assertEqual(sum(M), 27)
            self.assertTrue(P.array_ok(M), M)

    def test_identities(self):
        N, L = 27, self.L
        C = [comb(N, k) for k in range(P.KMAX + 1)]
        tot_c = [0] * (P.KMAX + 1)
        tot_h = [0] * (P.KMAX + 1)
        tot_lc = tot_xc = 0
        for M in self.arrays:
            c, h, lc, xc = P.raw_moments(M)
            for k in range(P.KMAX + 1):
                tot_c[k] += c[k]
                tot_h[k] += h[k]
            tot_lc += lc
            tot_xc += xc
        self.assertEqual(tot_c[2], 13 * C[2])
        self.assertEqual(tot_c[3], C[3] + 12 * L)
        self.assertEqual(tot_h[2], 13 * 13 * C[2])
        self.assertEqual(tot_h[3], 13 * (4 * C[3] + 9 * L))
        self.assertEqual(tot_h[4] - 39 * tot_c[4], 13 * C[4])
        for k in range(5, P.KMAX + 1):
            self.assertLessEqual(tot_h[k] - 39 * tot_c[k], 13 * C[k], k)
        self.assertGreaterEqual(tot_lc, 13 * L)
        self.assertGreaterEqual(tot_xc, 117 * L)

    def test_real_point_satisfies_lp(self):
        # the 130 real arrays with s = L - 40 and exact slacks solve the system of lp_system
        N, L = 27, self.L
        F = sorted({P.features(M) for M in self.arrays})
        cols, b = P.lp_system(N, 40, F)
        weight = {f: 0 for f in F}
        for M in self.arrays:
            weight[P.features(M)] += 1
        x = [weight[f] for f in F] + [L - 40]
        nk = P.KMAX - 4
        lhs = [sum(xj * col[i] for xj, col in zip(x, cols)) for i in range(len(b))]
        for k in range(nk):
            self.assertGreaterEqual(b[5 + k] - lhs[5 + k], 0)
        for r in (5 + nk, 6 + nk):
            self.assertGreaterEqual(lhs[r] - b[r], 0)
        for r in range(5):
            self.assertEqual(lhs[r], b[r], r)
        self.assertTrue(P.phase1_feasible(cols, b))

    def test_enumeration_finds_real_arrays(self):
        found = set(P.valid_arrays(27))
        for M in self.arrays:
            S = P.sort_rows(M)
            self.assertEqual(P.features(S), P.features(M))
            self.assertIn(S, found, M)


class Relaxation(unittest.TestCase):
    def test_known_cases(self):
        self.assertFalse(P.plane_feasible(24))   # already impossible by SD.md
        self.assertTrue(P.plane_feasible(27))    # a 27-point Kakeya set exists

    def test_sample_array(self):
        # an array with line profiles (6,6,13) and three times (8,8,9)
        M = (5, 4, 4, 2, 2, 2, 2, 2, 2)
        self.assertTrue(P.array_ok(M))
        c, h, lc, xc = P.raw_moments(M)
        self.assertEqual((c[2], c[3], c[4]), (28, 18, 7))
        self.assertEqual((h[2], h[3], h[4]), (384, 914, 1543))


if __name__ == "__main__":
    unittest.main()
