# test_plane_lp2.py -- checks for plane_lp2.py (PLANE2.md): the exact triple tables are symmetric,
# monotone, below the old capacities and match hand values; a real 27-point Kakeya set satisfies every
# new condition with its true number of lines; the sharper relaxation still admits N = 27.
# MV2 of Munim, Inc. Made with AI agents working in MV2's name. No human mathematician has checked
# this. Licence: Apache-2.0. Standard library only.

import unittest
from itertools import permutations

from moments import plane_lp as P
from moments import plane_lp2 as Q


class Tables(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.G, cls.T3 = Q.tables()

    def test_arithmetic(self):
        for i in range(9):
            for j in range(9):
                self.assertEqual(Q.ADD[i][Q.SUB[i][j]], j)
                self.assertEqual(Q.ADD[i][j], Q.ADD[j][i])

    def test_symmetric_in_all_three_sets(self):
        G, T3 = self.G, self.T3
        for a in range(10):
            for b in range(10):
                for c in range(10):
                    for p, q, r in permutations((a, b, c)):
                        self.assertEqual(G[a][b][c], G[p][q][r], (a, b, c))
                        self.assertEqual(T3[a][b][c], T3[p][q][r], (a, b, c))

    def test_monotone_and_below_old_caps(self):
        G, T3 = self.G, self.T3
        for a in range(10):
            for b in range(10):
                for c in range(10):
                    self.assertLessEqual(G[a][b][c], T3[a][b][c])
                    self.assertLessEqual(T3[a][b][c], min(a * b, a * c, b * c))
                    self.assertLessEqual(G[a][b][c], P.cross_cap(a, b, c), (a, b, c))
                    if c < 9:
                        self.assertLessEqual(G[a][b][c], G[a][b][c + 1])
                        self.assertLessEqual(T3[a][b][c], T3[a][b][c + 1])

    def test_hand_values(self):
        G, T3 = self.G, self.T3
        self.assertEqual(T3[1][1][1], 1)
        self.assertEqual(T3[2][2][2], 3)
        self.assertEqual(G[2][2][2], 3)
        self.assertEqual(T3[3][3][3], 9)
        self.assertEqual(G[1][9][9], 9)
        self.assertEqual(T3[9][9][9], 81)
        self.assertEqual(G[9][9][9], 9)
        self.assertEqual(T3[0][5][5], 0)


class RealSet(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.K = P.K27
        cls.L = P.count_lines(cls.K)
        cls.arrays = [P.array_of(cls.K, u, v) for u, v in P.plane_classes()]

    def test_new_conditions_hold(self):
        self.assertTrue(P.is_kakeya(self.K))
        tot_lc = tot_xc = 0
        for M in self.arrays:
            self.assertTrue(Q.array_ok2(M, self.L), M)
            lc, xc = Q.caps2(M)
            tot_lc += lc
            tot_xc += xc
        self.assertGreaterEqual(tot_lc, 13 * self.L)
        self.assertGreaterEqual(tot_xc, 117 * self.L)

    def test_real_point_satisfies_sharper_lp(self):
        N, L = 27, self.L
        F = sorted({Q.features2(M) for M in self.arrays})
        cols, b = P.lp_system(N, 40, F)
        weight = {f: 0 for f in F}
        for M in self.arrays:
            weight[Q.features2(M)] += 1
        x = [weight[f] for f in F] + [L - 40]
        nk = P.KMAX - 4
        lhs = [sum(xj * col[i] for xj, col in zip(x, cols)) for i in range(len(b))]
        for r in range(5):
            self.assertEqual(lhs[r], b[r], r)
        for k in range(nk):
            self.assertGreaterEqual(b[5 + k] - lhs[5 + k], 0)
        for r in (5 + nk, 6 + nk):
            self.assertGreaterEqual(lhs[r] - b[r], 0)
        found = set(Q.valid_arrays2(27))
        for M in self.arrays:
            self.assertIn(P.sort_rows(M), found, M)

    def test_sharper_relaxation_admits_27(self):
        self.assertTrue(Q.plane2_feasible(27))


if __name__ == "__main__":
    unittest.main()
