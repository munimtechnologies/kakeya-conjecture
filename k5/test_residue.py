# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Checks of RESIDUE.md: the residue lemma for hyperplane sections over F_3, and the three
explicit points that the profile method cannot exclude (26 points in F_3^4, 45 in F_3^5,
78 in F_3^6)."""
import itertools
import random
import unittest

try:
    from . import quad, residue
except ImportError:  # run from inside the folder
    import quad
    import residue


class TestResidueLemma(unittest.TestCase):
    def check_set(self, K, n):
        N = len(K)
        for xi in quad.proj_points(n):
            secs = quad.sections(K, xi)
            l, R = residue.moments(K, xi)
            self.assertEqual(tuple(t % 3 for t in secs), residue.predicted_sections(N, l, R))
            self.assertEqual(quad.QK(K, xi), (N * R - l * l) % 3)
            g = residue.group_of(secs, N)
            if N % 3 == 0 and l != 0:
                self.assertEqual(g, "P")
            elif N % 3 == 0 or l == 0:
                self.assertEqual(g, R)

    def test_all_subsets_f3_2(self):
        space = list(itertools.product(range(3), repeat=2))
        for mask in range(1, 2 ** 9):
            K = [space[i] for i in range(9) if mask >> i & 1]
            self.check_set(K, 2)

    def test_random_sets(self):
        rng = random.Random(20261125)
        for n, N in ((3, 13), (3, 12), (4, 26), (4, 27), (4, 30), (5, 44), (5, 45), (5, 31)):
            space = list(itertools.product(range(3), repeat=n))
            for _ in range(3):
                self.check_set(rng.sample(space, N), n)

    def test_rank_at_most_one_when_n_divisible_by_3(self):
        rng = random.Random(3)
        for n, N in ((3, 9), (3, 12), (4, 27), (5, 45)):
            space = list(itertools.product(range(3), repeat=n))
            H = (3 ** n - 1) // 2
            for _ in range(3):
                K = rng.sample(space, N)
                z = [0, 0, 0]
                for xi in quad.proj_points(n):
                    z[quad.QK(K, xi)] += 1
                self.assertIn(tuple(z), [(H, 0, 0), ((3 ** (n - 1) - 1) // 2, 0, 3 ** (n - 1))])

    def test_lev_set(self):
        lev = [
            (1, 2, 0, 0), (1, 0, 1, 0), (0, 1, 1, 0), (1, 2, 1, 0), (0, 0, 2, 0), (2, 1, 2, 0),
            (0, 0, 0, 1), (1, 0, 0, 1), (2, 0, 0, 1), (0, 1, 0, 1), (2, 0, 1, 1), (2, 1, 1, 1),
            (1, 2, 1, 1), (0, 0, 2, 1),
            (0, 0, 0, 2), (2, 0, 0, 2), (1, 1, 0, 2), (2, 1, 0, 2), (0, 2, 0, 2), (2, 2, 0, 2),
            (0, 0, 1, 2), (1, 1, 1, 2), (0, 2, 1, 2), (0, 0, 2, 2), (2, 0, 2, 2), (0, 1, 2, 2),
            (1, 1, 2, 2)]
        self.assertTrue(quad.is_kakeya(lev, 4))
        self.check_set(lev, 4)


class TestWhereTheProfileMethodStops(unittest.TestCase):
    """Explicit integer points meeting every profile-level constraint."""

    def test_points(self):
        for name, pt in residue.FEASIBLE.items():
            self.assertEqual(residue.check_point(pt), [], name)

    def test_quad_lp_agrees(self):
        # weak duality: quad.py's LP cannot exclude the quadric type of a feasible point
        for name, pt in residue.FEASIBLE.items():
            z = residue.q_type(pt)
            res = quad.certify(pt["n"], pt["N"], z)
            self.assertFalse(res["excluded"], (name, z, res))

    def test_profiles_sit_in_their_quad_groups(self):
        for name, pt in residue.FEASIBLE.items():
            N = pt["N"]
            qvals = set()
            for xi in quad.proj_points(pt["n"]):
                l, R = pt["l"](xi), pt["R"](xi)
                qvals.add((N * R - l * l) % 3)
            for p in pt["counts"]:
                self.assertIn(quad.residue(tuple(sorted(p)), N), qvals, (name, p))

    def test_n45_divisible_type_needs_17_to_23_heavy_classes(self):
        # l = 0, R = 0 at N = 45 in F_3^5: every size divisible by 3.  With D2 = s2 - 315 and
        # P = s3 - 1365 - 14 D2: sum D2 = 1485, sum P = 1782 + 27 (L - 121) >= 1782.
        sd2 = residue.sd2
        heavy = (9, 12, 24)
        group = [p for p in sd2.profiles(45, 81, "sd2") if all(t % 3 == 0 for t in p)]
        self.assertIn(heavy, group)
        S2, need = quad.totals(5, 45)
        self.assertEqual(S2 - 121 * 315, 1485)
        self.assertEqual(need - 121 * 1365 - 14 * 1485, 1782)
        for p in group:
            D2 = sd2.s2(p) - 315
            P = sd2.s3(p) - 1365 - 14 * D2
            self.assertGreaterEqual(D2, 0)
            if p == heavy:
                self.assertEqual((D2, P), (63, 81))
            else:
                self.assertLessEqual(P, D2, p)   # every other profile: P <= D2
        # x heavy classes: 1782 <= 81 x + (1485 - 63 x) and 63 x <= 1485
        ok = [x for x in range(0, 122) if 81 * x + (1485 - 63 * x) >= 1782 and 63 * x <= 1485]
        self.assertEqual((min(ok), max(ok)), (17, 23))

    def test_line_counts_are_integers_at_least_h(self):
        expect = {"f3_4_n26": (40, 53), "f3_5_n45": (121, 125), "f3_6_n78": (364, 364)}
        for name, pt in residue.FEASIBLE.items():
            H = (3 ** pt["n"] - 1) // 2
            self.assertEqual((H, pt["L"]), expect[name])


if __name__ == "__main__":
    unittest.main()
