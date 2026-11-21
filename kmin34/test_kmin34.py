# Tests for kmin34.py. MV2 of Munim, Inc. Made with AI agents working in MV2's name.
# Licence: Apache-2.0.
import itertools
import math
import random
import unittest

import numpy as np

import kmin34 as km

ALL3 = [tuple(p) for p in itertools.product(range(4), repeat=3)]


def perturb(rng, base, drop, extra):
    K = set(base)
    for _ in range(rng.randint(0, drop)):
        K.discard(rng.choice(sorted(K)))
    for _ in range(rng.randint(0, extra)):
        K.add(rng.choice(ALL3))
    return K


class TestField(unittest.TestCase):
    def test_axioms(self):
        A, M = km.ADD, km.MUL
        R = range(4)
        for a, b, c in itertools.product(R, R, R):
            self.assertEqual(A[A[a][b]][c], A[a][A[b][c]])
            self.assertEqual(M[M[a][b]][c], M[a][M[b][c]])
            self.assertEqual(M[a][A[b][c]], A[M[a][b]][M[a][c]])
        for a in R:
            self.assertEqual(A[a][0], a)
            self.assertEqual(M[a][1], a)
            if a:
                self.assertEqual(sum(1 for b in R if M[a][b] == 1), 1)
        self.assertEqual(M[2][2], 3)  # w^2 = w + 1, as in kakeya.py
        self.assertEqual(M[2][3], 1)
        self.assertEqual(M[3][3], 2)

    def test_directions_and_lines(self):
        self.assertEqual(len(km.DIR2), 5)
        self.assertEqual(len(km.DIR3), 21)
        for w in km.DIR2:
            self.assertEqual(len(km.LINES2[w]), 4)
            pts = sorted(i for L in km.LINES2[w] for i in L)
            self.assertEqual(pts, list(range(16)))


class TestKakeya(unittest.TestCase):
    def test_known28_is_kakeya(self):
        K = km.KNOWN28
        self.assertEqual(len(set(K)), 28)
        self.assertTrue(km.is_kakeya(K, 3))
        self.assertTrue(km.slice_kakeya(K))
        self.assertTrue(km.is_kakeya(ALL3, 3))
        no_bottom = [x for x in ALL3 if x[2] != 0]  # transversal lines all meet z = 0
        self.assertFalse(km.is_kakeya(no_bottom, 3))
        self.assertFalse(km.slice_kakeya(no_bottom))

    def test_slice_decomposition_matches_brute_force(self):
        rng = random.Random(7)
        seen = {True: 0, False: 0}
        for _ in range(80):
            K = perturb(rng, km.KNOWN28, 3, 6)
            a = km.is_kakeya(K, 3)
            self.assertEqual(a, km.slice_kakeya(K))
            seen[a] += 1
        self.assertGreater(seen[True], 0)
        self.assertGreater(seen[False], 0)

    def test_sections_partition(self):
        secs = km.sections(km.KNOWN28)
        self.assertEqual(len(secs), 84)
        by = {}
        for n, nrm, c in secs:
            by[nrm] = by.get(nrm, 0) + n
        self.assertEqual(len(by), 21)
        self.assertTrue(all(v == 28 for v in by.values()))


class TestGroup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.perms = km.agl2_perms()

    def test_group_is_agl24(self):
        P = self.perms
        self.assertEqual(P.shape, (km.GROUP_ORDER, 16))
        self.assertTrue((np.sort(P, axis=1) == np.arange(16)).all())
        self.assertEqual(len(np.unique(P, axis=0)), km.GROUP_ORDER)

    def test_group_preserves_lines(self):
        lines = set(L for w in km.DIR2 for L in km.LINES2[w])
        self.assertEqual(len(lines), 20)
        for g in self.perms[::7]:
            for L in lines:
                self.assertIn(tuple(sorted(int(g[i]) for i in L)), lines)

    def test_encode_decode(self):
        for S in ([0], [3, 7, 15], [0, 1, 2, 5, 8, 13, 14]):
            self.assertEqual(km.decode(km.encode(S), len(S)), sorted(S))

    def test_orbits(self):
        reps = km.orbit_reps(7, self.perms)
        self.assertEqual(len(reps[1]), 1)
        self.assertEqual(len(reps[2]), 1)
        self.assertEqual(len(reps[3]), 2)  # collinear or not
        oc = km.orbit_check(reps, self.perms)
        for k in range(1, 8):
            self.assertTrue(oc[k]["ok"], (k, oc[k]["sum_orbit_sizes"]))
            self.assertEqual(oc[k]["sum_orbit_sizes"], math.comb(16, k))

    def test_canon_invariant(self):
        rng = random.Random(11)
        for _ in range(20):
            S = rng.sample(range(16), rng.randint(1, 7))
            g = self.perms[rng.randrange(km.GROUP_ORDER)]
            T = [int(g[i]) for i in S]
            self.assertEqual(km.canon(S, self.perms)[0], km.canon(T, self.perms)[0])

    def test_normalise_and_evaluate(self):
        # Solver-free: a normalised Kakeya set satisfies MODEL(R, |K|) and violates
        # MODEL(R, |K| - 1); a normalised non-Kakeya set violates MODEL(R, 64).
        rng = random.Random(5)
        kinds = {True: 0, False: 0}
        for _ in range(40):
            K = sorted(perturb(rng, km.KNOWN28, 4, 8))
            n, code, sl = km.normalise(K, self.perms)
            if n == 0:
                continue
            Kn = km.from_slices(sl)
            self.assertEqual(len(Kn), len(K))
            kak = km.is_kakeya(K, 3)
            self.assertEqual(km.is_kakeya(Kn, 3), kak)
            self.assertEqual(min(s[0] for s in km.sections(Kn)), n)
            self.assertEqual(sl[0], km.decode(code, n))
            vals = km.point_vals_of(sl)
            if kak:
                self.assertTrue(km.evaluate(km.build_model(sl[0], len(K)), vals))
                self.assertFalse(km.evaluate(km.build_model(sl[0], len(K) - 1), vals))
            else:
                self.assertFalse(km.evaluate(km.build_model(sl[0], 64), vals))
            kinds[kak] += 1
        self.assertGreater(kinds[True], 0)
        self.assertGreater(kinds[False], 0)

    def test_known_set_accepted_by_solvers(self):
        out = km.known_check(km.KNOWN28, self.perms, 4, 180)
        self.assertTrue(out["ok"], out)

    def test_single_point_section_infeasible(self):
        # |K_0| = 1 forces K_1 = K_2 = K_3 = F_4^2 (a + t w covers the plane), so 49 > 27.
        R = km.decode(km.orbit_reps(1, self.perms)[1][0], 1)
        self.assertEqual(km.solve_cpsat(km.build_model(R, 27), 60, 2)[0], "INFEASIBLE")
        self.assertEqual(km.solve_highs(km.build_model(R, 27), 60)[0], "INFEASIBLE")
        self.assertEqual(km.solve_cpsat(km.build_model(R, 49), 60, 2)[0], "FEASIBLE")


if __name__ == "__main__":
    unittest.main()
