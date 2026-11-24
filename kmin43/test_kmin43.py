# Tests for kmin43.py (v3). MV2 of Munim, Inc. Made with AI agents working in MV2's name.
# Licence: Apache-2.0.
import itertools
import math
import random
import unittest

import numpy as np

import kmin43 as km

ALL4 = [tuple(p) for p in np.ndindex(3, 3, 3, 3)]


def is_kakeya_q3(K, n):
    """A second brute-force Kakeya test (all 3^n - 1 nonzero directions, no projective
    normalisation), independent of km.is_kakeya."""
    K = set(tuple(x) for x in K)
    for v in itertools.product(range(3), repeat=n):
        if not any(v):
            continue
        if not any(tuple((a + b) % 3 for a, b in zip(x, v)) in K
                   and tuple((a + 2 * b) % 3 for a, b in zip(x, v)) in K for x in K):
            return False
    return True


class TestBasics(unittest.TestCase):
    def test_directions(self):
        self.assertEqual(len(km.DIR3), 13)
        self.assertEqual(len(km.DIR4), 40)
        for w in km.DIR3:
            self.assertEqual(len(km.LINES3[w]), 9)
            pts = sorted(i for L in km.LINES3[w] for i in L)
            self.assertEqual(pts, list(range(27)))

    def test_classical(self):
        for n, size in ((1, 3), (2, 7), (3, 15), (4, 31)):
            cl = km.classical(n)
            self.assertEqual(len(cl), size, n)
            self.assertTrue(km.is_kakeya(cl, n), n)
            self.assertTrue(is_kakeya_q3(cl, n), n)

    def test_slice_decomposition_matches_brute_force(self):
        rng = random.Random(7)
        cl = [x[1:] + x[:1] for x in km.classical(4)]
        self.assertTrue(km.slice_kakeya(cl))
        seen = {True: 0, False: 0}
        for trial in range(60):
            K = set(cl)
            for _ in range(rng.randint(0, 3)):
                K.discard(rng.choice(sorted(K)))
            for _ in range(rng.randint(0, 6)):
                K.add(rng.choice(ALL4))
            a = km.is_kakeya(K, 4)
            self.assertEqual(a, km.slice_kakeya(K))
            self.assertEqual(a, is_kakeya_q3(K, 4))
            seen[a] += 1
        self.assertGreater(seen[True], 0)
        self.assertGreater(seen[False], 0)

    def test_second_moment_identity(self):
        # sum over the 120 sections of |K n H|^2 = 40 N + 13 N (N - 1), for any set
        rng = random.Random(5)
        for N in (1, 2, 7, 25, 26, 40, 81):
            for _ in range(5):
                K = rng.sample(ALL4, N)
                s2 = sum(n * n for n, _, _ in km.sections(K))
                self.assertEqual(s2, 40 * N + 13 * N * (N - 1))

    def test_nmax_hand_values(self):
        self.assertEqual(km.nmax_for_size(25), 6)
        self.assertEqual(km.nmax_for_size(26), 7)
        self.assertEqual(km.nmax_for_size(27), 7)
        self.assertEqual(km.nmax(25), 6)
        self.assertEqual(km.nmax(26), 7)

    def test_nmax_is_respected_by_random_sets(self):
        rng = random.Random(9)
        for N in (20, 25, 26, 27):
            for _ in range(200):
                K = rng.sample(ALL4, N)
                self.assertLessEqual(min(n for n, _, _ in km.sections(K)), km.nmax_for_size(N))

    def test_known_set_respects_nmax(self):
        self.assertLessEqual(min(n for n, _, _ in km.sections(km.KNOWN27)), km.nmax_for_size(27))
        cl = km.classical(4)
        self.assertLessEqual(min(n for n, _, _ in km.sections(cl)), km.nmax_for_size(31))


class TestGroup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.perms = km.agl3_perms()

    def test_group_is_agl33(self):
        P = self.perms
        self.assertEqual(P.shape, (km.GROUP_ORDER, 27))
        self.assertTrue((np.sort(P, axis=1) == np.arange(27)).all())
        self.assertEqual(len(np.unique(P, axis=0)), km.GROUP_ORDER)

    def test_group_preserves_lines(self):
        rng = random.Random(3)
        lines = set(L for w in km.DIR3 for L in km.LINES3[w])
        for _ in range(200):
            g = self.perms[rng.randrange(km.GROUP_ORDER)]
            L = rng.choice(sorted(lines))
            self.assertIn(tuple(sorted(int(g[i]) for i in L)), lines)

    def test_orbits_small(self):
        reps = km.orbit_reps(4, self.perms)
        self.assertEqual(len(reps[1]), 1)
        self.assertEqual(len(reps[2]), 1)
        self.assertEqual(len(reps[3]), 2)  # collinear or not
        oc = km.orbit_check(reps, self.perms)
        for k in range(1, 5):
            self.assertTrue(oc[k]["ok"], (k, oc[k]["sum_orbit_sizes"]))

    def test_canon_invariant(self):
        rng = random.Random(11)
        for _ in range(5):
            S = rng.sample(range(27), 6)
            g = self.perms[rng.randrange(km.GROUP_ORDER)]
            T = [int(g[i]) for i in S]
            self.assertEqual(km.canon(S, self.perms)[0], km.canon(T, self.perms)[0])

    def test_known_set_accepted(self):
        K = km.fetch_known(print) or km.KNOWN27
        self.assertEqual(len(set(K)), 27)
        self.assertTrue(km.is_kakeya(K, 4))
        self.assertTrue(is_kakeya_q3(K, 4))
        out = km.known_check(K, self.perms, 2, 180)
        self.assertTrue(out["ok"], out)

    def test_model_solution_is_checked(self):
        cl = [x[1:] + x[:1] for x in km.classical(4)]
        n, code, sl = km.normalise(cl, self.perms)
        for sym in (False, True):
            m = km.build_model(sl[0], 31, sym=sym)
            st, sol, _ = km.solve_cpsat(m, 120, 2)
            self.assertEqual(st, "FEASIBLE")
            ok, info = km.check_config(sl[0], sol, 31)
            self.assertTrue(ok, info)
        st2, sol2, _ = km.solve_highs(km.build_model(sl[0], 31, sym=True), 180)
        self.assertEqual(st2, "FEASIBLE")
        self.assertTrue(km.check_config(sl[0], sol2, 31)[0])

    def test_tiny_sections_infeasible(self):
        # |K_0| = 1 forces |K_1| = |K_2| = 27 (all 27 differences from one point), so 26 fails.
        R = km.decode(km.orbit_reps(1, self.perms)[1][0], 1)
        self.assertEqual(km.solve_cpsat(km.build_model(R, 26, sym=True), 60, 2)[0], "INFEASIBLE")
        self.assertEqual(km.solve_highs(km.build_model(R, 26, sym=True), 60)[0], "INFEASIBLE")


if __name__ == "__main__":
    unittest.main()
