# Tests for kmin43.py. MV2 of Munim, Inc. Made with AI agents working in MV2's name.
# Licence: Apache-2.0.
import math
import random
import unittest

import numpy as np

import kmin43 as km


class TestBasics(unittest.TestCase):
    def test_directions(self):
        self.assertEqual(len(km.DIR3), 13)
        self.assertEqual(len(km.DIR4), 40)
        for w in km.DIR3:
            self.assertEqual(len(km.LINES3[w]), 9)
            pts = sorted(i for L in km.LINES3[w] for i in L)
            self.assertEqual(pts, list(range(27)))

    def test_classical(self):
        cl = km.classical(4)
        self.assertEqual(len(cl), 31)
        self.assertTrue(km.is_kakeya(cl, 4))
        self.assertTrue(km.is_kakeya(km.classical(3), 3))
        self.assertEqual(len(km.classical(3)), 15)

    def test_slice_decomposition_matches_brute_force(self):
        rng = random.Random(7)
        cl = [x[1:] + x[:1] for x in km.classical(4)]
        self.assertTrue(km.slice_kakeya(cl))
        allpts = [tuple(p) for p in np.ndindex(3, 3, 3, 3)]
        seen = {True: 0, False: 0}
        for trial in range(60):
            K = set(cl)
            for _ in range(rng.randint(0, 3)):
                K.discard(rng.choice(sorted(K)))
            for _ in range(rng.randint(0, 6)):
                K.add(rng.choice(allpts))
            a = km.is_kakeya(K, 4)
            self.assertEqual(a, km.slice_kakeya(K))
            seen[a] += 1
        self.assertGreater(seen[True], 0)
        self.assertGreater(seen[False], 0)

    def test_sections_partition(self):
        cl = km.classical(4)
        secs = km.sections(cl)
        self.assertEqual(len(secs), 120)
        by = {}
        for n, nrm, c in secs:
            by.setdefault(nrm, 0)
            by[nrm] += n
        self.assertTrue(all(v == 31 for v in by.values()))


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

    def test_encode_decode(self):
        for S in ([0], [3, 7, 26], [0, 1, 2, 5, 8, 13, 21, 26, 25]):
            self.assertEqual(km.decode(km.encode(S), len(S)), sorted(S))

    def test_orbits_small(self):
        reps = km.orbit_reps(4, self.perms)
        self.assertEqual(len(reps[1]), 1)
        self.assertEqual(len(reps[2]), 1)
        self.assertEqual(len(reps[3]), 2)  # collinear or not
        oc = km.orbit_check(reps, self.perms)
        for k in range(1, 5):
            self.assertTrue(oc[k]["ok"], (k, oc[k]["sum_orbit_sizes"]))
            self.assertEqual(oc[k]["sum_orbit_sizes"], math.comb(27, k))

    def test_canon_invariant(self):
        rng = random.Random(11)
        for _ in range(5):
            S = rng.sample(range(27), 6)
            g = self.perms[rng.randrange(km.GROUP_ORDER)]
            T = [int(g[i]) for i in S]
            self.assertEqual(km.canon(S, self.perms)[0], km.canon(T, self.perms)[0])

    def test_known_set_accepted(self):
        K = km.fetch_known(print) or km.KNOWN27  # prefer the published certificate
        self.assertEqual(len(set(K)), 27)
        self.assertTrue(km.is_kakeya(K, 4))
        n, code, sl = km.normalise(K, self.perms)
        self.assertLessEqual(n, 9)
        Kn = ([km.PTS3[i] + (0,) for i in sl[0]] + [km.PTS3[i] + (1,) for i in sl[1]]
              + [km.PTS3[i] + (2,) for i in sl[2]])
        self.assertEqual(len(Kn), 27)
        self.assertTrue(km.is_kakeya(Kn, 4))
        self.assertEqual(min(s[0] for s in km.sections(Kn)), n)
        out = km.known_check(K, self.perms, 4, 180)
        self.assertTrue(out["ok"], out)

    def test_model_solution_is_checked(self):
        cl = [x[1:] + x[:1] for x in km.classical(4)]
        n, code, sl = km.normalise(cl, self.perms)
        m = km.build_model(sl[0], 31)
        st, sol, _ = km.solve_cpsat(m, 60, 2)
        self.assertEqual(st, "FEASIBLE")
        ok, info = km.check_config(sl[0], sol, 31)
        self.assertTrue(ok, info)
        st2, sol2, _ = km.solve_highs(m, 60)
        self.assertEqual(st2, "FEASIBLE")
        self.assertTrue(km.check_config(sl[0], sol2, 31)[0])

    def test_tiny_sections_infeasible(self):
        # |K_0| = 1 forces |K_1| = |K_2| = 27 (all 27 differences from one point), so 26 fails.
        R = km.decode(km.orbit_reps(1, self.perms)[1][0], 1)
        self.assertEqual(km.solve_cpsat(km.build_model(R, 26), 60, 2)[0], "INFEASIBLE")
        self.assertEqual(km.solve_highs(km.build_model(R, 26), 60)[0], "INFEASIBLE")


if __name__ == "__main__":
    unittest.main()
