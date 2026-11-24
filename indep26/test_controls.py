# indep26/test_controls.py - controls for the independent 26-point search in indep26/kk.py.
# MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name;
# no human mathematician has checked it yet. Licence: Apache-2.0. Python standard library only.
import itertools
import random
import unittest
from math import comb

from indep26 import kk


def _cov_brute(S, A, B):
    """True if B - A is all of F_3^m (computed directly)."""
    diffs = set()
    for y in kk.bits(B):
        for x in kk.bits(A):
            diffs.add(S.sub[y][x])
    return len(diffs) == S.N


class TestSpaceAndGroup(unittest.TestCase):
    def test_space(self):
        for m, ndirs in ((2, 4), (3, 13)):
            S = kk.space(m)
            self.assertEqual(len(S.dirs), ndirs)
            for L in S.lines_by_dir:
                self.assertEqual(len(L), S.N // 3)
                self.assertEqual(sum(l.bit_count() for l in L), S.N)    # a partition
                u = 0
                for l in L:
                    u |= l
                self.assertEqual(u, (1 << S.N) - 1)

    def test_group_orders(self):
        S2, S3 = kk.space(2), kk.space(3)
        self.assertEqual(kk.group_order(kk.linear_generators(S2), 9), 48)          # |GL(2,3)|
        self.assertEqual(kk.group_order(kk.linear_generators(S3), 27), 11232)      # |GL(3,3)|
        self.assertEqual(kk.group_order(kk.linear_generators(S2) + kk.translation_generators(S2), 9),
                         432)                                                       # |AGL(2,3)|

    def test_orbits(self):
        S2, S3 = kk.space(2), kk.space(3)
        for a in range(10):
            reps = kk.orbit_reps(S2, a)
            self.assertEqual(sum(o for _, o in reps), comb(9, a))
            for _, o in reps:
                self.assertEqual(432 % o, 0)
        self.assertEqual(len(kk.orbit_reps(S3, 3)), 2)       # a line, or three non-collinear points
        for a in (1, 2, 3, 4, 5, 6, 7):
            reps = kk.orbit_reps(S3, a)
            self.assertEqual(sum(o for _, o in reps), comb(27, a))
            for R, o in reps:
                self.assertEqual(R.bit_count(), a)
                self.assertEqual(303264 % o, 0)
        for R, o in kk.orbit_reps(S3, 5):
            self.assertEqual(len(kk.orbit_of(S3, R)), o)


class TestEnumB(unittest.TestCase):
    def test_enum_B_dim2_all(self):
        S = kk.space(2)
        for a in range(1, 6):
            for A, _ in kk.orbit_reps(S, a):
                for b in range(1, 8):
                    brute = []
                    for rest in itertools.combinations(range(1, 9), b - 1):
                        B = 1 | kk.mask_of(rest)
                        if _cov_brute(S, A, B):
                            brute.append(B)
                    self.assertEqual(sorted(kk.enum_B(S, A, b)), sorted(brute), (a, A, b))

    def test_enum_B_dim3_one_case(self):
        S = kk.space(3)
        A = kk.orbit_reps(S, 5)[0][0]
        cov = [kk.mask_of(S.sub[y][x] for x in kk.bits(A)) for y in range(27)]
        full = (1 << 27) - 1
        brute = []
        for rest in itertools.combinations(range(1, 27), 6):
            u = cov[0]
            for y in rest:
                u |= cov[y]
            if u == full:
                brute.append(1 | kk.mask_of(rest))
        got = kk.enum_B(S, A, 7)
        self.assertEqual(sorted(got), sorted(brute))
        for B in got[:200]:
            self.assertTrue(_cov_brute(S, A, B))


class TestCSearch(unittest.TestCase):
    def test_c_search_matches_brute_force_minimum_dim2(self):
        S = kk.space(2)
        checked = 0
        for a in (2, 3, 4):
            for A, _ in kk.orbit_reps(S, a):
                for rest_size in range(0, 5):
                    for rest in itertools.combinations(range(1, 9), rest_size):
                        B = 1 | kk.mask_of(rest)
                        if not _cov_brute(S, A, B):
                            continue
                        # independent brute force over all 512 third sections
                        Q = []
                        for v in range(9):
                            q = 0
                            for x in kk.bits(A):
                                y = S.add[x][v]
                                if (B >> y) & 1:
                                    q |= 1 << S.add[y][v]
                            Q.append(q)
                        best = None
                        for C in range(512):
                            if best is not None and C.bit_count() >= best:
                                continue
                            if all(q & C for q in Q) and kk.valid_ABC(S, A, B, C):
                                best = C.bit_count()
                        found_min = None
                        for c in range(10):
                            C = kk.solve_AB(S, A, B, c)
                            if C is not None:
                                self.assertLessEqual(C.bit_count(), c)
                                self.assertTrue(kk.valid_ABC(S, A, B, C))
                                found_min = c
                                break
                        self.assertEqual(found_min, best, (A, B))
                        checked += 1
        self.assertGreater(checked, 20)


class TestIdentitiesAndCuts(unittest.TestCase):
    def test_product_and_pair_identities(self):
        rng = random.Random(20261124)
        cases = [(3, n) for n in (8, 12, 13, 15)] + [(4, n) for n in (20, 25, 26, 27, 30)]
        sets = [(4, kk.KNOWN27)]
        for n, N in cases:
            pts = list(itertools.product(range(3), repeat=n))
            sets.append((n, rng.sample(pts, N)))
        for n, K in sets:
            N = len(K)
            L = kk.count_lines(K)
            sizes = kk.class_sizes(K, n)
            self.assertEqual(len(sizes), (3 ** n - 1) // 2)
            self.assertEqual(sum(x * y * z for x, y, z in sizes), 3 ** (n - 2) * (comb(N, 3) + 2 * L))
            self.assertEqual(sum(x * y + x * z + y * z for x, y, z in sizes), 3 ** (n - 1) * comb(N, 2))
        self.assertTrue(kk.count_lines(kk.KNOWN27) >= 40)

    def test_cut_list(self):
        self.assertEqual(sorted(kk.cuts26_from_identity()), sorted(kk.CUTS26))
        self.assertEqual(len(kk.CUTS26), 10)
        # N = 26, n = 4: sum over the 40 classes of (n1-8)(n2-8)(n3-8) = 18 (L - 40).
        # Per class (n1-8)(n2-8)(n3-8) = e3 - 8 e2 + 64*26 - 512; e1 = 26.
        for a, b, c in kk.cuts_T(26, 3):
            e2 = a * b + a * c + b * c
            self.assertEqual((a - 8) * (b - 8) * (c - 8), a * b * c - 8 * e2 + 1152)
        self.assertEqual(9 * comb(26, 3) - 8 * 27 * comb(26, 2) + 40 * 1152, -720)
        # zero case: sum over classes of n1^2+n2^2+n3^2 = 40*26^2 - 2*27*C(26,2) = 9490, while a
        # class with every section >= 8 and one section = 8 has at most 228.
        self.assertEqual(40 * 26 ** 2 - 2 * 27 * comb(26, 2), 9490)
        big = [t for t in kk.cuts_T(26, 3) if min(t) >= 8 and 8 in t]
        self.assertEqual(sorted(big), [(8, 8, 10), (8, 9, 9)])
        self.assertEqual(max(sum(x * x for x in t) for t in big), 228)
        self.assertLess(40 * 228, 9490)


class TestPipeline(unittest.TestCase):
    def test_quick_controls(self):
        self.assertTrue(kk.quick_controls())

    def test_full_pipeline_finds_known27_in_its_orbit(self):
        # KNOWN27, class x_3 = t: sections of sizes 13, 6, 8. In the cyclic order (6, 8, 13) the
        # pipeline (orbit representative of the 6-section, B through 0, C <= 13) must find a set.
        S = kk.space(3)
        secs = kk.coordinate_sections(kk.KNOWN27, 3)
        self.assertEqual([s.bit_count() for s in secs], [13, 6, 8])
        reps = kk.orbit_reps(S, 6)
        orb = kk.orbit_of(S, secs[1])
        idx = [i for i, (R, _) in enumerate(reps) if R in orb]
        self.assertEqual(len(idx), 1)
        R, o = reps[idx[0]]
        self.assertEqual(o, len(orb))
        found = None
        for B in kk.enum_B(S, R, 8):
            C = kk.solve_AB(S, R, B, 13)
            if C is not None:
                found = (B, C)
                break
        self.assertIsNotNone(found)
        B, C = found
        self.assertTrue(kk.valid_ABC(S, R, B, C))
        K = kk.assemble(S, R, B, C)
        self.assertTrue(kk.is_kakeya(K, 4))
        self.assertLessEqual(len(K), 27)


if __name__ == "__main__":
    unittest.main()
