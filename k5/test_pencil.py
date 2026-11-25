# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Checks of the pencil LP formulation (pencil.py) on actual Kakeya sets, standard library only.

Every constraint of build_lp is checked to hold exactly for the true vector of actual Kakeya
sets (Lev's 27-point set in F_3^4, random minimal sets in F_3^3 and F_3^4, a random set in
F_3^5); the symmetry-reduced enumeration is compared with brute force; and the exact Farkas
checker is tested on a toy system.
"""
import itertools
import unittest
from collections import Counter
from math import comb

try:
    from . import pencil
except ImportError:
    import pencil

LEV = [(1, 2, 0, 0), (1, 0, 1, 0), (0, 1, 1, 0), (1, 2, 1, 0), (0, 0, 2, 0), (2, 1, 2, 0),
       (0, 0, 0, 1), (1, 0, 0, 1), (2, 0, 0, 1), (0, 1, 0, 1), (2, 0, 1, 1), (2, 1, 1, 1),
       (1, 2, 1, 1), (0, 0, 2, 1),
       (0, 0, 0, 2), (2, 0, 0, 2), (1, 1, 0, 2), (2, 1, 0, 2), (0, 2, 0, 2), (2, 2, 0, 2),
       (0, 0, 1, 2), (1, 1, 1, 2), (0, 2, 1, 2), (0, 0, 2, 2), (2, 0, 2, 2), (0, 1, 2, 2),
       (1, 1, 2, 2)]


class TestCounting(unittest.TestCase):
    def test_constants(self):
        self.assertEqual([pencil.gauss2(m) for m in range(1, 6)], [0, 1, 13, 130, 1210])
        for n in range(3, 7):
            k = pencil.consts(n)
            self.assertEqual(4 * k["P"], k["a1"] * k["H"])
            self.assertEqual(k["H"] - k["a1"], k["T"])

    def test_pencils_and_pairs(self):
        for n in (3, 4):
            k = pencil.consts(n)
            bases = pencil.pencil_bases(n)
            self.assertEqual(len(bases), k["P"])
            per_class, per_pair = Counter(), Counter()
            for x, y in bases:
                cl = [pencil.norm(tuple((a * s + b * t) % 3 for s, t in zip(x, y)))
                      for a, b in pencil.DIRS]
                self.assertEqual(len(set(cl)), 4)
                for c in cl:
                    per_class[c] += 1
                for c1, c2 in itertools.combinations(sorted(cl), 2):
                    per_pair[(c1, c2)] += 1
            self.assertEqual(set(per_class.values()), {k["a1"]})
            self.assertEqual(len(per_pair), comb(k["H"], 2))
            self.assertEqual(set(per_pair.values()), {1})

    def test_plane_tables(self):
        pts = list(itertools.product(range(3), repeat=2))
        lines = set()
        for p in pts:
            for d in pencil.proj_points(2):
                lines.add((frozenset(tuple((a + t * b) % 3 for a, b in zip(p, d))
                                     for t in range(3)), d))
        self.assertEqual(len(lines), 12)
        best_l, best_d = [0] * 10, [0] * 10
        for mask in range(512):
            S = frozenset(pts[i] for i in range(9) if mask >> i & 1)
            inside = [(L, d) for L, d in lines if L <= S]
            best_l[len(S)] = max(best_l[len(S)], len(inside))
            best_d[len(S)] = max(best_d[len(S)], len(set(d for _, d in inside)))
        self.assertEqual(best_l, pencil.MAXLINES2)
        self.assertEqual(best_d, pencil.MAXDIRS2)

    def test_line_formula_bound(self):
        for m in range(28):
            self.assertLessEqual(pencil.max_dirs(m, 3), 13)
            self.assertLessEqual(pencil.max_dirs(m, 3), pencil.max_lines(m, 3))
        self.assertEqual(pencil.max_dirs(12, 3), 12)
        self.assertEqual(pencil.max_dirs(13, 3), 13)


class TestEnumeration(unittest.TestCase):
    def brute(self, n, N):
        D = pencil.consts(n)["D"]
        allowed = set(pencil.profile_list(n, N))
        out = set()
        for m in itertools.product(range(D + 1), repeat=8):
            last = N - sum(m)
            if 0 <= last <= D:
                f = pencil.array_features(m + (last,), n, allowed)
                if f is not None:
                    out.add(f)
        return out

    def test_symmetry_reduction_matches_brute_force(self):
        for N in (12, 13, 15):
            self.assertEqual(set(pencil.enumerate_features(3, N)), self.brute(3, N), N)

    def test_prune_keeps_a_dominator(self):
        E = pencil.enumerate_features(4, 27)
        kept = pencil.prune(E)
        self.assertTrue(set(kept) <= set(E))
        by_key = {}
        for g in kept:
            by_key.setdefault((g[0], g[1]), []).append(g)
        for f in E:
            self.assertTrue(any(g[2] >= f[2] and g[3] >= f[3]
                                for g in by_key[(f[0], f[1])]))


class TestRealSets(unittest.TestCase):
    def check_real(self, K, n, enum=True):
        self.assertTrue(pencil.is_kakeya(K, n))
        k = pencil.consts(n)
        data = pencil.real_data(K, n)
        N = data["N"]
        self.assertNotIn(None, data["feats"])
        self.assertEqual(len(data["arrays"]), k["P"])
        self.assertEqual(sum(sum(comb(v, 2) for v in m) for m in data["arrays"]),
                         k["B1"] * comb(N, 2))
        self.assertEqual(sum(sum(comb(v, 3) for v in m) for m in data["arrays"]),
                         k["B2"] * comb(N, 3) + (k["B1"] - k["B2"]) * data["L"])
        profs = pencil.profile_list(n, N)
        for p in data["classes"]:
            self.assertIn(p, profs)
        if enum:
            E = pencil.enumerate_features(n, N)
            for f in data["feats"]:
                self.assertIn(f, E)
            featlist = sorted(E)
        else:
            featlist = sorted(data["feats"])
        z = [0, 0, 0]
        for p, c in data["classes"].items():
            z[pencil.residue(p, N)] += c
        self.assertIn(tuple(z), pencil.quadric_types(n))
        U = pencil.unbalanced_set(profs, N, 1)
        u = sum(c for p, c in data["classes"].items() if p in U)
        eq, ub, bounds = pencil.build_lp(n, N, profs, featlist, tuple(z), U, u)
        v = pencil.real_vector(n, data, profs, featlist)
        self.assertEqual(pencil.check_point(eq, ub, v), [])
        for x, b in zip(v, bounds):
            self.assertTrue(0 <= x <= b)
        return N

    def test_lev(self):
        self.assertEqual(self.check_real(LEV, 4), 27)

    def test_random_f3_3(self):
        sizes = set()
        for seed in range(6):
            sizes.add(self.check_real(pencil.greedy_kakeya(3, seed), 3))
        self.assertTrue(min(sizes) >= 13)

    def test_random_f3_4(self):
        self.assertGreaterEqual(self.check_real(pencil.greedy_kakeya(4, 1), 4), 27)

    def test_random_f3_5(self):
        self.check_real(pencil.greedy_kakeya(5, 0), 5, enum=False)


class TestFarkas(unittest.TestCase):
    def test_toy(self):
        eq = [({0: 1, 1: 1}, 1)]
        ub = [({0: -1}, -2)]
        bounds = [5, 5]
        self.assertTrue(pencil.check_certificate(eq, ub, bounds, {0: 1}, {0: 1})[0])
        self.assertFalse(pencil.check_certificate(eq, ub, bounds, {0: -1}, {})[0])
        self.assertFalse(pencil.check_certificate(eq, ub, bounds, {0: 1}, {0: -1})[0])
        # a feasible system admits no certificate: x0 = 1 works for x0 + x1 = 1, x0 >= 1
        eq2, ub2 = [({0: 1, 1: 1}, 1)], [({0: -1}, -1)]
        for lam in (-2, -1, 0, 1, 2):
            for mu in (0, 1, 2):
                self.assertFalse(pencil.check_certificate(eq2, ub2, bounds, {0: lam},
                                                          {0: mu})[0])


if __name__ == "__main__":
    unittest.main()
