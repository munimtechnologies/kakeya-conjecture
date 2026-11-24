# Tests for kprof25.py. MV2 of Munim, Inc. Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0.
import itertools
import random
import unittest
from math import comb

import numpy as np

import kprof25 as kp

ALL4 = [tuple(p) for p in itertools.product(range(3), repeat=4)]


def is_kakeya_q3(K):
    """Second brute-force test over all 80 non-zero vectors (no projective normalisation)."""
    K = set(tuple(x) for x in K)
    for v in ALL4:
        if not any(v):
            continue
        if not any(tuple((a + b) % 3 for a, b in zip(x, v)) in K
                   and tuple((a + 2 * b) % 3 for a, b in zip(x, v)) in K for x in K):
            return False
    return True


class TestStatic(unittest.TestCase):
    def test_lines(self):
        self.assertEqual(len(kp.LINES4), 1080)
        per = {}
        for L, v in kp.LINES4:
            per[v] = per.get(v, 0) + 1
            self.assertEqual(len(set(L)), 3)
        self.assertEqual(sorted(per.values()), [27] * 40)

    def test_profiles_and_table(self):
        self.assertEqual(kp.allowed_profiles(25), kp.HAND25)
        self.assertEqual(sorted(p for p in kp.HAND25 if kp.delta25(p) > 0),
                         sorted(c[1] for c in kp.CASES))
        t = kp.table(25)
        self.assertEqual(len(t), 69)
        for row in t:
            self.assertEqual(sum(row[:3]), 25)
            self.assertEqual(kp.moments_of(row[:3]), (row[3], row[4]))
        # (6,6,13) has 3 orders and (5,7,13) has 6
        self.assertEqual(len(kp.table(25, [(6, 6, 13), (5, 7, 13)])), 60)

    def test_delta_identity_random(self):
        rng = random.Random(1)
        for _ in range(20):
            K = rng.sample(ALL4, 25)
            cuts = [kp.class_profile(K, u) for u in kp.DIR4]
            L = kp.count_lines(K)
            self.assertEqual(sum(kp.delta25(c) for c in cuts), 63 * L - 2420)

    def test_known_set(self):
        self.assertTrue(kp.is_kakeya(kp.KNOWN27))
        self.assertTrue(is_kakeya_q3(kp.KNOWN27))
        allowed = set(kp.allowed_profiles(27))
        for u in kp.DIR4:
            self.assertIn(tuple(sorted(kp.class_profile(kp.KNOWN27, u))), allowed)

    def test_spec_matches_brute_force(self):
        rng = random.Random(5)
        for trial in range(25):
            a = rng.choice([5, 6, 9])
            b = rng.randint(a, 12)
            c = rng.randint(b, 20)
            R = sorted(rng.sample(range(27), a))
            K1 = rng.sample(range(27), b)
            K2 = rng.sample(range(27), c)
            if trial % 5 == 0:  # include a full line of F_3^4 crossing the slices
                y = rng.randrange(27)
                w = kp.PTS3[rng.randrange(27)]
                p0 = kp.PTS3[y]
                R = sorted(set(R[:-1]) | {y})
                i1 = kp.IDX3[kp.add(p0, w)]
                i2 = kp.IDX3[kp.add(p0, kp.smul(2, w))]
                K1 = list(set(K1[:-1]) | {i1})
                K2 = list(set(K2[:-1]) | {i2})
            spec = kp.make_spec(len(R) + len(K1) + len(K2), R, len(K1), len(K2))
            sol = [0] * 54
            for i in K1:
                sol[kp.X(1, i)] = 1
            for i in K2:
                sol[kp.X(2, i)] = 1
            K = kp.assemble(R, K1, K2)
            n_lines = sum(1 for _, vs in spec["lines"] if all(sol[j] for j in vs))
            self.assertEqual(n_lines, kp.count_lines(K))
            for u, trip in zip(kp.DIR4, spec["secs"]):
                sizes = tuple(const + sum(sol[j] for j in vs) for const, vs in trip)
                self.assertEqual(sizes, kp.class_profile(K, u))
            self.assertEqual(kp.class_profile(K, kp.LAST), (len(R), len(K1), len(K2)))


class TestGroupAndSolvers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.perms = kp.agl3_perms()

    def test_group(self):
        P = self.perms
        self.assertEqual(P.shape, (kp.GROUP_ORDER, 27))
        self.assertTrue((np.sort(P, axis=1) == np.arange(27)).all())

    def test_orbits(self):
        reps = kp.orbit_reps(4, self.perms)
        self.assertEqual([len(reps[k]) for k in (1, 2, 3)], [1, 1, 2])
        oc = kp.orbit_check(reps, self.perms)
        self.assertTrue(all(v["ok"] for v in oc.values()), oc)

    def test_normalise_every_class(self):
        for u in kp.DIR4:
            sl = kp.normalise_class(kp.KNOWN27, u, self.perms)
            K = kp.assemble(sl[0], sl[1], sl[2])
            self.assertEqual(len(set(K)), 27)
            self.assertTrue(is_kakeya_q3(K))
            self.assertEqual(kp.class_profile(K, kp.LAST),
                             tuple(sorted(kp.class_profile(kp.KNOWN27, u))))
            self.assertEqual(sorted(tuple(sorted(kp.class_profile(K, v))) for v in kp.DIR4),
                             sorted(tuple(sorted(kp.class_profile(kp.KNOWN27, v))) for v in kp.DIR4))

    def test_known_sanity(self):
        out = kp.known_sanity(kp.KNOWN27, self.perms, 2, 150, highs=True)
        self.assertTrue(out["ok"], out)

    def test_impossible_sizes_infeasible(self):
        # a 3-point slice 0 with sizes (3, 9, 13) violates (SD): not in the table, so infeasible
        R = [0, 1, 2]
        spec = kp.make_spec(25, R, 9, 13)
        self.assertEqual(kp.solve_cpsat(spec, 60, 2)[0], "INFEASIBLE")
        self.assertEqual(kp.solve_highs(spec, 60)[0], "INFEASIBLE")


if __name__ == "__main__":
    unittest.main()
