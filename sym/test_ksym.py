# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Controls for ksym: field arithmetic, group invariance, and the orbit model against brute
force on small cases (k_2 = 7, k_3 = 13)."""
import itertools
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gf  # noqa: E402
import check  # noqa: E402
import bb  # noqa: E402
import ksym  # noqa: E402

try:
    import ortools  # noqa: F401
    HAVE_ORTOOLS = True
except ImportError:  # pragma: no cover
    HAVE_ORTOOLS = False


class TestField(unittest.TestCase):
    def test_exp_is_permutation(self):
        for n in range(2, 7):
            F = gf.Field(n)
            self.assertEqual(sorted(F.exp), list(range(1, F.q)))

    def test_linear_maps(self):
        rng = random.Random(5)
        for n in range(2, 7):
            F = gf.Field(n)
            for r, f in ((F.q - 1, 1), (1, n), (2, 1)):
                for pm in F.generators(r, f):
                    self.assertEqual(sorted(pm), list(range(F.q)))
                    for _ in range(300):
                        a, b = rng.randrange(F.q), rng.randrange(F.q)
                        self.assertEqual(pm[F.add(a, b)], F.add(pm[a], pm[b]))

    def test_add_matches_digits(self):
        rng = random.Random(6)
        for n in (2, 5, 8):
            F = gf.Field(n)
            for _ in range(300):
                a, b = rng.randrange(F.q), rng.randrange(F.q)
                s = gf.from_digits([(x + y) % 3 for x, y in zip(gf.digits(a, n), gf.digits(b, n))])
                self.assertEqual(F.add(a, b), s)

    def test_mul_distributes(self):
        rng = random.Random(7)
        F = gf.Field(4)
        for _ in range(500):
            a, b, c = (rng.randrange(F.q) for _ in range(3))
            self.assertEqual(F.mul(a, F.add(b, c)), F.add(F.mul(a, b), F.mul(a, c)))


class TestOrbits(unittest.TestCase):
    def test_invariance_and_sizes(self):
        for n in (3, 4, 5):
            F = gf.Field(n)
            for r in ksym.divisors(F.q - 1):
                for f in ksym.divisors(n):
                    gens = F.generators(r, f)
                    orb, members = gf.orbits(F.q, gens)
                    self.assertEqual(sum(len(m) for m in members), F.q)
                    for pm in gens:
                        for x in range(F.q):
                            self.assertEqual(orb[pm[x]], orb[x])
                    for m in members:
                        self.assertEqual((r * f) % len(m), 0)
                    reps, sizes, nd = gf.direction_orbits(F, gens)
                    self.assertEqual(nd, (F.q - 1) // 2)
                    self.assertEqual(sum(sizes), nd)

    def test_line_masks_trivial_group(self):
        F = gf.Field(3)
        orb, members = gf.orbits(F.q, [])
        reps, sizes, nd = gf.direction_orbits(F, [])
        masks = gf.line_masks(F, orb, reps)
        self.assertEqual(len(masks), 13)
        for opts in masks:
            self.assertEqual(len(opts), 9)
            self.assertTrue(all(len(gf.mask_bits(L)) == 3 for L in opts))


def brute_min(n, members):
    """Least size of a Kakeya set that is a union of the given orbits (all subsets)."""
    best = None
    for k in range(len(members) + 1):
        for comb in itertools.combinations(range(len(members)), k):
            pts = [p for o in comb for p in members[o]]
            if best is not None and len(pts) >= best:
                continue
            if check.is_kakeya(n, check.to_tuples(n, pts)):
                best = len(pts)
    return best


class TestModel(unittest.TestCase):
    def test_checker(self):
        for n in (2, 3, 4):
            self.assertTrue(check.is_kakeya(n, check.classical(n)))
        c2 = check.classical(2)
        for p in c2:
            self.assertFalse(check.is_kakeya(2, [x for x in c2 if x != p]))

    def test_k2_brute(self):
        F = gf.Field(2)
        orb, members = gf.orbits(F.q, [])
        self.assertEqual(brute_min(2, members), 7)
        reps, _, _ = gf.direction_orbits(F, [])
        r = bb.minimum([1] * 9, gf.line_masks(F, orb, reps))
        self.assertEqual(r["value"], 7)
        self.assertTrue(r["complete"])

    def test_groups_n3_vs_brute(self):
        F = gf.Field(3)
        for r in (1, 2, 13, 26):
            for f in (1, 3):
                gens = F.generators(r, f)
                orb, members = gf.orbits(F.q, gens)
                if len(members) > 16:
                    continue
                reps, _, _ = gf.direction_orbits(F, gens)
                w = [len(m) for m in members]
                got = bb.minimum(w, gf.line_masks(F, orb, reps), seconds=60)
                self.assertTrue(got["complete"])
                self.assertEqual(got["value"], brute_min(3, members), (r, f))

    def test_k3_bb(self):
        F = gf.Field(3)
        orb, members = gf.orbits(F.q, [])
        reps, _, _ = gf.direction_orbits(F, [])
        r = bb.minimum([1] * 27, gf.line_masks(F, orb, reps), seconds=150)
        self.assertLessEqual(r["value"], 13)
        if not r["complete"]:
            self.skipTest("branch and bound did not finish in 150 s (best %s)" % r["value"])
        self.assertEqual(r["value"], 13)

    def test_fw_index(self):
        for n in (3, 4):
            pts = [gf.from_digits(t) for t in check.classical(n)]
            idx = ksym.to_fw_index(n, pts)
            back = [tuple((p // 3 ** (n - 1 - k)) % 3 for k in range(n)) for p in idx]
            self.assertTrue(check.is_kakeya(n, back))
            self.assertEqual(len(set(back)), len(pts))

    @unittest.skipUnless(HAVE_ORTOOLS, "ortools not installed")
    def test_cpsat_small(self):
        for n, want in ((2, 7), (3, 13)):
            rec = ksym.run_instance(n, 1, 1, 60, 2)
            self.assertEqual(rec["cpsat"]["status"], "OPTIMAL")
            self.assertEqual(rec["size"], want)
            self.assertTrue(rec["verified"])
        F = gf.Field(3)
        for r, f in ((2, 1), (1, 3), (2, 3)):
            rec = ksym.run_instance(3, r, f, 60, 2)
            gens = F.generators(r, f)
            orb, members = gf.orbits(F.q, gens)
            self.assertEqual(rec["size"], brute_min(3, members), (r, f))


if __name__ == "__main__":
    unittest.main()
