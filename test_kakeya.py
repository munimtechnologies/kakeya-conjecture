import contextlib
import io
import os
import tempfile
import unittest
from itertools import combinations

from kakeya import (GF, KNOWN_F2, AffineSpace, brute_force_minimum, bukh_chao_bound,
                    counting_bound, dvir_bound, is_difference_cover_f2, is_kakeya, main,
                    missing_directions, parabola_kakeya, parabola_size, plane_minimum,
                    product_set, smyth_kakeya, smyth_size, symmetry_lines)


class FieldTests(unittest.TestCase):
    def test_tables_are_fields(self):
        for q in (2, 3, 4, 5, 7, 8, 9, 16):
            F = GF(q)
            R = range(q)
            for a in R:
                self.assertEqual(F.add[a][0], a)
                self.assertEqual(F.mul[a][1], a)
                self.assertEqual(F.add[a][F.neg[a]], 0)
                self.assertEqual(sorted(F.add[a]), list(R))
                if a:
                    self.assertEqual(sorted(F.mul[a]), list(R), (q, a))
                for b in R:
                    self.assertEqual(F.add[a][b], F.add[b][a])
                    self.assertEqual(F.mul[a][b], F.mul[b][a])

    def test_associative_and_distributive(self):
        for q in (4, 8, 9):
            F = GF(q)
            for a in range(q):
                for b in range(q):
                    for c in range(q):
                        self.assertEqual(F.mul[F.mul[a][b]][c], F.mul[a][F.mul[b][c]])
                        self.assertEqual(F.add[F.add[a][b]][c], F.add[a][F.add[b][c]])
                        self.assertEqual(F.mul[a][F.add[b][c]],
                                         F.add[F.mul[a][b]][F.mul[a][c]])

    def test_characteristic(self):
        for q, p in ((4, 2), (8, 2), (9, 3), (16, 2)):
            F = GF(q)
            for a in range(q):
                s = 0
                for _ in range(p):
                    s = F.add[s][a]
                self.assertEqual(s, 0)

    def test_unsupported(self):
        with self.assertRaises(ValueError):
            GF(6)


class GeometryTests(unittest.TestCase):
    def test_counts(self):
        for n, q in ((2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (2, 5), (4, 3), (3, 4)):
            S = AffineSpace(n, q)
            self.assertEqual(len(S.directions), (q ** n - 1) // (q - 1))
            for lines in S.lines:
                self.assertEqual(len(lines), q ** (n - 1))
                covered = set()
                for L in lines:
                    self.assertEqual(len(L), q)
                    self.assertFalse(covered & L)
                    covered |= L
                self.assertEqual(covered, set(range(q ** n)))

    def test_two_points_lie_on_one_line(self):
        S = AffineSpace(2, 4)
        all_lines = [L for lines in S.lines for L in lines]
        self.assertEqual(len(all_lines), 20)
        for a in range(16):
            for b in range(a + 1, 16):
                self.assertEqual(sum(1 for L in all_lines if a in L and b in L), 1)

    def test_is_kakeya(self):
        S = AffineSpace(2, 3)
        self.assertTrue(is_kakeya(range(9), S))
        self.assertFalse(is_kakeya([], S))
        self.assertTrue(is_kakeya(set(range(9)) - {4}, S))
        self.assertFalse(is_kakeya(S.lines[0][0], S))
        self.assertEqual(len(missing_directions([], S)), 4)
        self.assertEqual(missing_directions(range(9), S), [])

    def test_symmetry_lines(self):
        for n, q in ((2, 3), (3, 2), (3, 3), (4, 3)):
            S = AffineSpace(n, q)
            fixed = symmetry_lines(S)
            self.assertEqual(len(fixed), 2)
            for i, lines in fixed.items():
                self.assertEqual(len(lines), 1 if (n == 2 or i == 0) else 2)
                for L in lines:
                    self.assertIn(L, S.lines[i])
            self.assertIn(0, fixed[0][0])


class MinimumTests(unittest.TestCase):
    def test_plane_formula(self):
        self.assertEqual([plane_minimum(q) for q in (2, 3, 4, 5, 7, 8, 9)],
                         [3, 7, 10, 17, 31, 36, 49])

    def test_plane_brute_force_matches_known_minimum(self):
        for q in (2, 3, 4, 5):
            S = AffineSpace(2, q)
            best, pts = brute_force_minimum(S)
            self.assertEqual(best, plane_minimum(q), q)
            self.assertEqual(len(pts), best)
            self.assertTrue(is_kakeya(pts, S))
            fixed_best, fixed_pts = brute_force_minimum(S, symmetry_lines(S))
            self.assertEqual(fixed_best, best, q)
            self.assertTrue(is_kakeya(fixed_pts, S))

    def test_binary_cube(self):
        # Over F_2 a line is a pair {x, x + d}, so a Kakeya set must realise all 7 nonzero
        # differences: k points give k(k-1)/2 >= 7 pairs, so k >= 5, and
        # {0, e1, e2, e3, e1+e2+e3} works.
        S = AffineSpace(3, 2)
        best, pts = brute_force_minimum(S)
        self.assertEqual(best, 5)
        self.assertEqual(best, KNOWN_F2[3])
        self.assertTrue(is_kakeya(pts, S))
        self.assertTrue(is_kakeya([0, 1, 2, 4, 7], S))
        self.assertFalse(is_kakeya([0, 1, 2, 4], S))
        self.assertEqual(brute_force_minimum(S, symmetry_lines(S))[0], 5)
        self.assertEqual(brute_force_minimum(AffineSpace(2, 2))[0], KNOWN_F2[2])

    def test_bounds(self):
        self.assertEqual(dvir_bound(3, 3), 10)
        self.assertEqual(dvir_bound(2, 5), 15)
        self.assertEqual(bukh_chao_bound(3, 3), 10)
        self.assertEqual(bukh_chao_bound(3, 5), 39)
        self.assertEqual(bukh_chao_bound(2, 3), 6)
        self.assertEqual(bukh_chao_bound(1, 7), 7)
        # 2187/125, 19683/625, 1024/49, 16384/343, 16807/169 rounded up.
        self.assertEqual([bukh_chao_bound(4, 3), bukh_chao_bound(5, 3), bukh_chao_bound(3, 4),
                          bukh_chao_bound(4, 4), bukh_chao_bound(3, 7)], [18, 32, 21, 48, 100])
        self.assertEqual([dvir_bound(4, 3), dvir_bound(5, 3), dvir_bound(3, 4), dvir_bound(4, 4),
                          dvir_bound(3, 5), dvir_bound(3, 7)], [15, 21, 20, 35, 35, 84])
        for q in (2, 3, 4, 5, 7, 8, 9):
            self.assertLessEqual(bukh_chao_bound(2, q), plane_minimum(q))
            self.assertLessEqual(counting_bound(2, q), plane_minimum(q))
            self.assertLessEqual(dvir_bound(2, q), plane_minimum(q))

    def test_dvir_is_exact_in_even_planes(self):
        # binom(q + 1, 2) = q(q + 1)/2, the even-q plane minimum.
        for q in (2, 4, 8, 16):
            self.assertEqual(dvir_bound(2, q), plane_minimum(q))
        self.assertEqual([dvir_bound(2, q) for q in (2, 4, 8, 16)], [3, 10, 36, 136])

    def test_counting_bound(self):
        self.assertEqual([counting_bound(n, 2) for n in range(1, 7)], [2, 3, 5, 6, 9, 12])
        for n, v in KNOWN_F2.items():
            self.assertLessEqual(counting_bound(n, 2), v)
        self.assertEqual(counting_bound(2, 3), 6)
        self.assertEqual(counting_bound(3, 3), 10)
        self.assertEqual(counting_bound(2, 5), 12)


class BinaryTests(unittest.TestCase):
    def test_difference_cover_equals_kakeya_in_the_cube(self):
        S = AffineSpace(3, 2)
        for mask in range(256):
            pts = [i for i in range(8) if mask >> i & 1]
            self.assertEqual(is_kakeya(pts, S), is_difference_cover_f2(pts, 3), pts)

    def test_four_dimensional_binary_minimum_is_six(self):
        S = AffineSpace(4, 2)
        for pts in combinations(range(16), 5):
            self.assertFalse(is_kakeya(pts, S))
            self.assertFalse(is_difference_cover_f2(pts, 4))
        found = 0
        for pts in combinations(range(16), 6):
            k = is_kakeya(pts, S)
            self.assertEqual(k, is_difference_cover_f2(pts, 4))
            found += k
        self.assertGreater(found, 0)
        # {0, e1, e2, e3, e4, e1+e2+e3+e4}: its 15 pairwise sums are the 15 nonzero vectors.
        sidon = [0, 1, 2, 4, 8, 15]
        self.assertTrue(is_kakeya(sidon, S))
        self.assertEqual(len({a ^ b for a, b in combinations(sidon, 2)}), 15)
        self.assertEqual(KNOWN_F2[4], 6)


class ConstructionTests(unittest.TestCase):
    def test_parabola_sets(self):
        for n, q in ((1, 3), (2, 3), (2, 5), (2, 7), (2, 9), (3, 3), (3, 5), (4, 3), (3, 7)):
            S, pts = parabola_kakeya(n, q)
            self.assertTrue(is_kakeya(pts, S), (n, q))
            self.assertEqual(len(pts), parabola_size(n, q), (n, q))
            if n == 2:
                self.assertEqual(len(pts), plane_minimum(q), q)
        self.assertEqual([parabola_size(3, 3), parabola_size(3, 5), parabola_size(4, 3),
                          parabola_size(3, 7)], [17, 61, 43, 145])

    def test_smyth_sets(self):
        for n, q in ((1, 3), (2, 2), (2, 3), (2, 4), (3, 3), (4, 3), (3, 4), (3, 5), (5, 3),
                     (4, 4)):
            S, pts = smyth_kakeya(n, q)
            self.assertTrue(is_kakeya(pts, S), (n, q))
            self.assertEqual(len(pts), smyth_size(n, q), (n, q))
        self.assertEqual(smyth_size(2, 3), plane_minimum(3))
        self.assertEqual([smyth_size(3, 3), smyth_size(4, 3), smyth_size(5, 3), smyth_size(3, 4),
                          smyth_size(4, 4)], [15, 31, 63, 34, 96])
        # Removing any 0/1 vector other than the all-ones vector loses a direction.
        S, pts = smyth_kakeya(3, 3)
        for v in (0, 1, 3, 4, 9, 10, 12):
            self.assertFalse(is_kakeya([p for p in pts if p != v], S), v)

    def test_parabola_needs_odd_q(self):
        with self.assertRaises(ValueError):
            parabola_kakeya(2, 4)

    def test_products(self):
        _, plane3 = parabola_kakeya(2, 3)
        prod = product_set(plane3, 2, plane3, 3)
        self.assertEqual(len(prod), 49)
        self.assertTrue(is_kakeya(prod, AffineSpace(4, 3)))
        _, plane2 = brute_force_minimum(AffineSpace(2, 2))
        prod = product_set(plane2, 2, [0, 1], 2)
        self.assertEqual(len(prod), 6)
        self.assertTrue(is_kakeya(prod, AffineSpace(3, 2)))
        prod = product_set([0, 1, 2, 4, 7], 3, [0, 1], 2)
        self.assertEqual(len(prod), 10)
        self.assertTrue(is_kakeya(prod, AffineSpace(4, 2)))


class CommandLineTests(unittest.TestCase):
    def run_main(self, argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = main(argv)
        return code, out.getvalue()

    def test_bounds_and_parabola(self):
        code, text = self.run_main(["bounds", "3", "3"])
        self.assertEqual(code, 0)
        self.assertIn("Dvir bound       10", text)
        self.assertIn("parabola set     17", text)
        code, text = self.run_main(["parabola", "2", "3"])
        self.assertEqual(code, 0)
        self.assertEqual(len(text.strip().splitlines()), 7)

    def test_check(self):
        S, pts = parabola_kakeya(2, 5)
        with tempfile.TemporaryDirectory() as d:
            good = os.path.join(d, "good.txt")
            bad = os.path.join(d, "bad.txt")
            with open(good, "w") as f:
                f.write("# parabola set in AG(2,5)\n")
                for p in pts:
                    f.write("%d %d\n" % S.coords[p])
            with open(bad, "w") as f:
                for p in pts[1:]:
                    f.write("%d %d\n" % S.coords[p])
            self.assertEqual(self.run_main(["check", "2", "5", good])[0], 0)
            code, text = self.run_main(["check", "2", "5", bad])
            self.assertEqual(code, 1)
            self.assertIn("NOT a Kakeya set", text)

    def test_usage(self):
        self.assertEqual(self.run_main([])[0], 2)


if __name__ == "__main__":
    unittest.main()
