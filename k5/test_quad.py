# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Checks of QUAD45.md: the mod-3 quadric constraint on hyperplane sections, and the proof
that no Kakeya set in F_3^5 has 44 points (so k_5 >= 45)."""
import itertools
import random
import unittest
from math import comb

try:
    from . import sd2, quad
except ImportError:  # run from inside the folder
    import sd2
    import quad

# Lev's 27-point Kakeya set in F_3^4 (as in k5/k5cut.py)
LEV27 = [
    (1, 2, 0, 0), (1, 0, 1, 0), (0, 1, 1, 0), (1, 2, 1, 0), (0, 0, 2, 0), (2, 1, 2, 0),
    (0, 0, 0, 1), (1, 0, 0, 1), (2, 0, 0, 1), (0, 1, 0, 1), (2, 0, 1, 1), (2, 1, 1, 1),
    (1, 2, 1, 1), (0, 0, 2, 1),
    (0, 0, 0, 2), (2, 0, 0, 2), (1, 1, 0, 2), (2, 1, 0, 2), (0, 2, 0, 2), (2, 2, 0, 2),
    (0, 0, 1, 2), (1, 1, 1, 2), (0, 2, 1, 2), (0, 0, 2, 2), (2, 0, 2, 2), (0, 1, 2, 2),
    (1, 1, 2, 2)]

EXPECTED_N5 = sorted([(121, 0, 0), (40, 81, 0), (40, 0, 81), (13, 54, 54), (67, 27, 27),
                      (40, 54, 27), (40, 27, 54), (49, 36, 36), (31, 45, 45), (40, 45, 36),
                      (40, 36, 45)])


def all_form_counts(n):
    pts = quad.proj_points(n)
    mons = quad.monomials(n)
    prods = [[v[i] * v[j] for (i, j) in mons] for v in pts]
    seen = set()
    for coefs in itertools.product(range(3), repeat=len(mons)):
        z = [0, 0, 0]
        for pr in prods:
            z[sum(c * m for c, m in zip(coefs, pr)) % 3] += 1
        seen.add(tuple(z))
    return seen


class TestQuadricTypes(unittest.TestCase):
    def test_point_counts(self):
        for n in range(1, 7):
            self.assertEqual(len(quad.proj_points(n)), (3 ** n - 1) // 2)

    def test_n5_table(self):
        got = sorted(z for _, z in quad.quadric_types(5))
        self.assertEqual(got, EXPECTED_N5)

    def test_classification_exhaustive_n3_n4(self):
        for n in (3, 4):
            reps = set(z for _, z in quad.quadric_types(n))
            self.assertEqual(all_form_counts(n), reps, n)

    def test_classification_random_n5(self):
        rng = random.Random(20261125)
        reps = set(z for _, z in quad.quadric_types(5))
        pts = quad.proj_points(5)
        mons = quad.monomials(5)
        for _ in range(1500):
            coefs = [rng.randrange(3) for _ in mons]
            self.assertIn(quad.value_counts(coefs, pts, mons), reps)


class TestModThreeIdentity(unittest.TestCase):
    def test_random_sets(self):
        rng = random.Random(7)
        for n, N in ((3, 13), (3, 9), (4, 26), (4, 27), (5, 44), (5, 30)):
            space = list(itertools.product(range(3), repeat=n))
            for _ in range(2):
                K = rng.sample(space, N)
                for xi in quad.proj_points(n):
                    self.assertEqual((quad.s2_of(K, xi) - comb(N, 2) + quad.QK(K, xi)) % 3, 0)

    def test_lev_set(self):
        self.assertTrue(quad.is_kakeya(LEV27, 4))
        T = 27
        for xi in quad.proj_points(4):
            p = tuple(sorted(quad.sections(LEV27, xi)))
            self.assertTrue(sd2.allowed(p, T, "sd2"), p)          # SD2 holds on a real set
            self.assertEqual(quad.residue(p, 27), quad.QK(LEV27, xi))


class TestControls(unittest.TestCase):
    """The method must not exclude sizes where Kakeya sets exist."""

    def test_f3_3_size_13_not_excluded(self):
        self.assertFalse(quad.excluded(3, 13)[0])

    def test_f3_4_size_27_not_excluded(self):
        self.assertFalse(quad.excluded(4, 27)[0])

    def test_lev_profiles_feasible_for_its_own_type(self):
        # Lev's set itself is a feasible point of the LP for its own quadric type.
        pts = quad.proj_points(4)
        z = [0, 0, 0]
        for xi in pts:
            z[quad.QK(LEV27, xi)] += 1
        res = quad.certify(4, 27, z)
        self.assertFalse(res["excluded"])


class TestNoKakeya44(unittest.TestCase):
    """k_5 >= 45: every quadric type is excluded at N = 44."""

    def test_all_types_excluded_with_checked_certificates(self):
        ok, per = quad.excluded(5, 44)
        self.assertTrue(ok)
        for name, v in per.items():
            self.assertTrue(v["excluded"], name)
            if v["reason"] == "lp":
                self.assertTrue(quad.check_certificate(5, 44, v["z"], v["mu"]), name)

    def test_hand_proof_facts(self):
        T = 81
        P = sd2.profiles(44, T, "sd2")
        base, heavy = (14, 15, 15), (10, 11, 23)

        def d2(p):
            return sd2.s2(p) - 301

        def f(p):
            return (sd2.s3(p) - 1274) - 15 * d2(p)

        for p in P:
            prod = (p[0] - 16) * (p[1] - 16) * (p[2] - 16)
            self.assertEqual(prod % 2, 0)
            self.assertEqual(f(p), 1 - 2 * d2(p) + prod // 2)          # closed form
            self.assertEqual(f(p) % 3, 2 if d2(p) % 3 == 1 else 0)      # f mod 3
            self.assertEqual(quad.residue(p, 44), (-d2(p)) % 3)         # Q_K = -Delta2
            if p not in (base, heavy):
                self.assertLessEqual(f(p), -1)
                self.assertGreaterEqual(d2(p), 1)
                if d2(p) % 3 != 1:
                    self.assertLessEqual(f(p), -3)
        self.assertEqual((f(base), d2(base), f(heavy), d2(heavy)), (0, 0, 2, 52))
        # identities: sum Delta2 = 1419, sum f = 27 (L - 121)
        S2, need = quad.totals(5, 44)
        self.assertEqual(S2 - 121 * 301, 1419)
        self.assertEqual(need - 121 * 1274 - 15 * 1419, 0)
        # quadric types: z1 <= 18 only for the zero form and 2 l^2 (z0 = 40)
        small = sorted(z for z in EXPECTED_N5 if z[1] <= 18)
        self.assertEqual(small, [(40, 0, 81), (121, 0, 0)])
        # those two: k heavy classes, d defects, b base classes <= z0
        for z in small:
            for k in range(0, 28):
                for d in range(0, 122):
                    b = 121 - k - d
                    if b < 0 or b > z[0] or k > z[2]:
                        continue
                    # defects: |f| >= 1 each, total |f| <= 2k; Delta2 >= 1 each, total 1419 - 52k
                    feasible = (d <= 2 * k) and (d <= 1419 - 52 * k) and (k > 0 or d == 0)
                    if k == 0:
                        feasible = feasible and (1419 == 0)
                    self.assertFalse(feasible, (z, k, d))


if __name__ == "__main__":
    unittest.main()
