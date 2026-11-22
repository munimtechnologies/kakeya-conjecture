# Tests for kakeya_moments.py. MV2 of Munim, Inc. Made with AI agents. Apache-2.0.
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import kakeya_moments as km  # noqa: E402


class TestCounts(unittest.TestCase):
    def test_flat_counts(self):
        # F_3^4: 40 hyperplane classes, 13 hyperplanes through a line, 4 through a plane
        self.assertEqual((km.gauss1(3, 4), km.gauss1(3, 3), km.gauss1(3, 2)), (40, 13, 4))
        # F_3^3: 13 plane classes, 4 planes through a line, 1 through three non-collinear points
        self.assertEqual((km.gauss1(3, 3), km.gauss1(3, 2), km.gauss1(3, 1)), (13, 4, 1))
        self.assertEqual(km.gauss1(3, 0), 0)

    def test_targets_f34(self):
        # N = 23 in F_3^4: sum C(n_H,2) = 13*253, sum C(n_H,3) >= 4*1771 + 9*40
        self.assertEqual(km.targets(3, 4, 23), (40, 3289, 7444))

    def test_profile_lists_f34(self):
        # hand-enumerated lists in README.md
        self.assertEqual(km.profiles(3, 4, 21), [
            (3, 9, 9), (4, 7, 10), (4, 8, 9), (5, 6, 10), (5, 7, 9), (5, 8, 8),
            (6, 6, 9), (6, 7, 8), (7, 7, 7)])
        self.assertEqual(len(km.profiles(3, 4, 22)), 11)
        self.assertEqual(km.profiles(3, 4, 23), [
            (3, 9, 11), (3, 10, 10), (4, 7, 12), (4, 8, 11), (4, 9, 10), (5, 6, 12),
            (5, 7, 11), (5, 8, 10), (5, 9, 9), (6, 6, 11), (6, 7, 10), (6, 8, 9),
            (7, 7, 9), (7, 8, 8)])
        self.assertEqual(len(km.profiles(3, 4, 24)), 17)
        self.assertEqual(km.profiles(3, 4, 16), [])
        self.assertEqual(km.profiles(3, 4, 17), [(5, 6, 6)])

    def test_moments(self):
        self.assertEqual(km.moments((7, 8, 8)), (77, 147))
        self.assertEqual(km.moments((5, 6, 13)), (103, 316))


class TestKnownValues(unittest.TestCase):
    def test_plane(self):
        # K(2,3) = 7 (Blokhuis-Mazzocca: q(q+1)/2 + (q-1)/2 for odd q)
        self.assertEqual(km.lower_bound(3, 2), 7)

    def test_f33(self):
        # K(3,3) = 13 is known by exhaustive search; the moment bound reaches it exactly
        self.assertEqual(km.lower_bound(3, 3), 13)
        c = km.certificate(3, 3, 12, 16, 5)
        self.assertEqual(c["max_profile_value"], -226)
        self.assertEqual(c["collinear_triples_at_most"], 12)
        self.assertEqual(c["collinear_triples_needed"], 13)


class TestF34(unittest.TestCase):
    def test_small_N_by_second_moment(self):
        # N = 17, 18, 19: even the largest possible sum of C(n_H,2) is too small
        for N in (17, 18, 19):
            C, S2, _ = km.targets(3, 4, N)
            top = max(km.moments(P)[0] for P in km.profiles(3, 4, N))
            self.assertLess(C * top, S2, N)

    def test_certificates(self):
        # hand certificates of README.md: N, c, max profile value, lines at most
        for N, cn, cd, top, lmax in [(20, 13, 2, -561, 30), (21, 7, 1, -336, 38),
                                     (22, 7, 1, -364, 33), (23, 37, 5, -2114, 38)]:
            c = km.certificate(3, 4, N, cn, cd)
            self.assertEqual(c["max_profile_value"], top, N)
            self.assertEqual(c["collinear_triples_at_most"], lmax, N)
            self.assertLess(lmax, 40)

    def test_bound_24(self):
        for N in range(1, 24):
            self.assertFalse(km.feasible(3, 4, N), N)
        self.assertTrue(km.feasible(3, 4, 24))
        self.assertEqual(km.lower_bound(3, 4), 24)


if __name__ == "__main__":
    unittest.main()
