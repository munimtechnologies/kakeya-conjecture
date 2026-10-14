import random
import unittest
from itertools import combinations, permutations, product

from triples import (G1, K33, PLANE_PATTERNS, TWO_SYMBOL_BASE, add, directions, doubling_chain,
                     energy_excludes, energy_lower, energy_upper, in_b, in_w, is_kakeya,
                     is_kakeya_by_layers, is_rainbow, layered_step, lift, missing_differences,
                     product_triple, rainbow_with_sizes, scale, sizes, slices, smyth, spread,
                     tripled, tripled_point, vectors, wwb_base, wwb_kakeya_size, wwb_product,
                     wwb_sizes, wwb_triple, wwb_witness)


def subsets(points):
    for r in range(len(points) + 1):
        for c in combinations(points, r):
            yield set(c)


class BasicTests(unittest.TestCase):
    def test_vectors_and_directions(self):
        for n in range(5):
            self.assertEqual(len(vectors(n)), 3 ** n)
            self.assertEqual(len(directions(n)), (3 ** n - 1) // 2)

    def test_g1_is_rainbow(self):
        self.assertTrue(is_rainbow(*G1, 1))
        A, M, C = G1
        self.assertTrue(M <= A and M <= {(0,), (1,)})

    def test_missing_differences(self):
        self.assertEqual(missing_differences({(0,)}, {(0,)}, {(0,)}, 1), [(1,), (2,)])
        self.assertEqual(missing_differences(*G1, 1), [])


class SymmetryTests(unittest.TestCase):
    def test_all_six_orders_agree_in_f3(self):
        subs = list(subsets(vectors(1)))
        for T in product(subs, repeat=3):
            values = {is_rainbow(*(T[i] for i in p), 1) for p in permutations(range(3))}
            self.assertEqual(len(values), 1, T)

    def test_all_six_orders_agree_in_f3_squared(self):
        rng = random.Random(20261014)
        pts = vectors(2)
        cases = [slices(K33), product_triple(G1, G1)]
        for _ in range(400):
            cases.append(tuple(set(rng.sample(pts, rng.randint(2, 7))) for _ in range(3)))
        for T in cases:
            values = {is_rainbow(*(T[i] for i in p), 2) for p in permutations(range(3))}
            self.assertEqual(len(values), 1, T)

    def test_translation_invariance(self):
        T = slices(K33)
        for c in vectors(2):
            moved = [{tuple((x + y) % 3 for x, y in zip(p, c)) for p in S} for S in T]
            self.assertTrue(is_rainbow(*moved, 2))


class LayerLemmaTests(unittest.TestCase):
    def test_exhaustive_in_the_plane(self):
        for K in subsets(vectors(2)):
            self.assertEqual(is_kakeya(K, 2), is_kakeya_by_layers(K, 1), K)

    def test_random_sets_in_f3_cubed(self):
        rng = random.Random(7)
        pts = vectors(3)
        cases = [set(K33)] + [set(K33) - {p} for p in K33]
        for _ in range(300):
            cases.append(set(rng.sample(pts, rng.randint(10, 18))))
        kakeya = 0
        for K in cases:
            a = is_kakeya(K, 3)
            self.assertEqual(a, is_kakeya_by_layers(K, 2))
            kakeya += a
        self.assertGreaterEqual(kakeya, 1)

    def test_k33_slices(self):
        S = slices(K33)
        self.assertEqual(sizes(S), (4, 6, 3))
        self.assertTrue(is_rainbow(*S, 2))
        self.assertEqual(lift(*S), set(K33))
        for p in K33:
            smaller = slices(set(K33) - {p})
            self.assertFalse(is_kakeya(set(K33) - {p}, 3))
            self.assertFalse(is_rainbow(*smaller, 2) and all(
                any(_line(smaller[t], d) for t in range(3)) for d in directions(2)))


def _line(S, d):
    return any(tuple((a[i] + d[i]) % 3 for i in range(len(d))) in S
               and tuple((a[i] + 2 * d[i]) % 3 for i in range(len(d))) in S for a in S)


class ConstructionTests(unittest.TestCase):
    def test_doubling_chain_from_a_point(self):
        for k in range(6):
            K = smyth(k)
            self.assertEqual(len(K), 2 ** (k + 1) - 1)
            self.assertTrue(is_kakeya(K, k))

    def test_doubling_chain_from_k33(self):
        T = (set(K33), set(K33), set(K33))
        for k in range(3):
            B, T2 = doubling_chain(K33, T, 3, k)
            self.assertEqual(len(B), 13 + (2 ** k - 1) * 26)
            self.assertTrue(is_kakeya(B, 3 + k))
            self.assertTrue(is_rainbow(*T2, 3 + k))
            self.assertTrue(T2[1] <= T2[0] and T2[1] <= B)

    def test_doubling_chain_rejects_bad_middle(self):
        with self.assertRaises(ValueError):
            doubling_chain(K33, ({(0, 0, 0)}, {(1, 1, 1)}, {(0, 0, 0)}), 3, 1)

    def test_layered_step(self):
        S = slices(K33)
        B = smyth(2)
        K = layered_step(B, S)
        self.assertTrue(is_kakeya(K, 3))
        self.assertEqual(len(K), len(S[0]) + len(B | S[1]) + len(S[2]))

    def test_products_and_cyclic_product(self):
        T = product_triple(G1, G1)
        self.assertEqual(sizes(T), (4, 4, 4))
        self.assertTrue(is_rainbow(*T, 2))
        A, M, C = G1
        star = product_triple(product_triple(G1, (M, C, A)), (C, A, M))
        self.assertEqual(sizes(star), (8, 8, 8))
        self.assertTrue(is_rainbow(*star, 3))
        mixed = product_triple(slices(K33), G1)
        self.assertTrue(is_rainbow(*mixed, 3))


class SmallestTriplesTests(unittest.TestCase):
    """Exhaustive over all 512 triples of subsets of F_3."""

    def setUp(self):
        subs = list(subsets(vectors(1)))
        self.rainbow = [T for T in product(subs, repeat=3) if is_rainbow(*T, 1)]

    def test_sigma_1(self):
        self.assertEqual(min(sum(sizes(T)) for T in self.rainbow), 6)

    def test_max_and_product_1(self):
        self.assertEqual(min(max(sizes(T)) for T in self.rainbow), 2)
        a, b, c = min((sizes(T) for T in self.rainbow), key=lambda s: s[0] * s[1] * s[2])
        self.assertEqual(a * b * c, 8)

    def test_pairwise_products(self):
        for T in self.rainbow:
            a, b, c = sizes(T)
            self.assertTrue(a * b >= 3 and b * c >= 3 and a * c >= 3)

    def test_doubling_constant_1(self):
        best = min(len(A) + len(C) for A, M, C in self.rainbow if M <= A)
        self.assertEqual(best, 4)

    def test_exhaustive_search_matches_brute_force(self):
        feasible = {sizes(T) for T in self.rainbow}
        for p in product(range(1, 4), repeat=3):
            T = rainbow_with_sizes(*p, 1)
            self.assertEqual(T is not None, p in feasible, p)
            if T is not None:
                self.assertEqual(sizes(T), p)
                self.assertTrue(is_rainbow(*T, 1))

    def test_energy_lemma_is_sound_in_f3(self):
        feasible = {sizes(T) for T in self.rainbow}
        for p in product(range(1, 4), repeat=3):
            if energy_excludes(*p, 1):
                self.assertNotIn(p, feasible, p)
        for p in feasible:
            self.assertFalse(energy_excludes(*p, 1), p)


class EnergyLemmaTests(unittest.TestCase):
    def test_spread(self):
        self.assertEqual(spread(9, 4, 3), 21)
        self.assertEqual(spread(9, 5, 3), 17)
        self.assertEqual(spread(27, 6, 6), 123)
        self.assertEqual(spread(0, 0, 3), 0)
        self.assertIsNone(spread(1, 0, 3))
        self.assertIsNone(spread(13, 4, 3))

    def test_hand_values(self):
        self.assertEqual(energy_upper(3, 4, 9), 20)
        self.assertEqual(energy_lower(3, 4, 4, 9), 24)
        self.assertEqual(energy_lower(3, 4, 5, 9), 20)
        self.assertEqual(energy_upper(2, 5, 9), 12)
        self.assertEqual(energy_lower(2, 5, 5, 9), 18)
        self.assertEqual(energy_upper(6, 6, 27), 86)
        self.assertEqual(energy_lower(6, 6, 6, 27), 132)
        self.assertIsNone(energy_upper(2, 4, 9))
        self.assertIsNone(energy_lower(2, 4, 5, 9))

    def test_energy_matches_real_triples(self):
        cases = [(product_triple(G1, G1), 2), (slices(K33), 2),
                 (product_triple(slices(K33), G1), 3)]
        for T, n in cases:
            N = 3 ** n
            A, M, C = T
            for X, Y, Z in ((A, M, C), (A, C, M), (M, C, A)):
                counts = {}
                for x in X:
                    for y in Y:
                        v = tuple((p + q) % 3 for p, q in zip(x, y))
                        counts[v] = counts.get(v, 0) + 1
                E = sum(c * c for c in counts.values())
                self.assertLessEqual(energy_lower(len(X), len(Y), len(Z), N), E)
                self.assertLessEqual(E, energy_upper(len(X), len(Y), N))

    def test_plane_patterns_by_energy(self):
        excluded = {p: energy_excludes(*p, 2) for p in PLANE_PATTERNS}
        self.assertEqual([p for p in PLANE_PATTERNS if not excluded[p]], [(3, 4, 5)])
        self.assertFalse(energy_excludes(4, 4, 4, 2))

    def test_s3_at_least_seven(self):
        for a in range(1, 7):
            for b in range(a, 7):
                for c in range(b, 7):
                    self.assertTrue(energy_excludes(a, b, c, 3), (a, b, c))
        self.assertFalse(energy_excludes(7, 7, 7, 3))
        self.assertFalse(energy_excludes(6, 7, 7, 3))


class PlaneTests(unittest.TestCase):
    """tau(2) = 64, sigma(2) = 12, s(2) = 4: no rainbow triple of F_3^2 beats base 2."""

    def test_pattern_list_is_complete(self):
        pats = [(a, b, c) for a in range(1, 10) for b in range(a, 10) for c in range(b, 10)
                if a * b * c < 64 and a * b >= 9]
        self.assertEqual(pats, PLANE_PATTERNS)
        small_sum = [(a, b, c) for a in range(1, 10) for b in range(a, 10) for c in range(b, 10)
                     if a + b + c <= 11 and a * b >= 9]
        self.assertTrue(set(small_sum) <= set(PLANE_PATTERNS))
        self.assertIn((3, 3, 3), PLANE_PATTERNS)

    def test_no_triple_below_base_two(self):
        for p in PLANE_PATTERNS:
            self.assertIsNone(rainbow_with_sizes(*p, 2), p)

    def test_four_four_four_exists(self):
        T = rainbow_with_sizes(4, 4, 4, 2)
        self.assertIsNotNone(T)
        self.assertEqual(sizes(T), (4, 4, 4))
        self.assertTrue(is_rainbow(*T, 2))


class TwoSymbolTests(unittest.TestCase):
    """(W, W, B) is rainbow; tripled, it is a Kakeya set; its sizes give base 432^(1/9) < 2."""

    def test_witness_on_random_vectors(self):
        rng = random.Random(432)
        for m in (1, 2, 3, 5, 8, 13, 30, 61):
            for _ in range(300):
                d = tuple(rng.randrange(3) for _ in range(m))
                a = wwb_witness(d)
                self.assertTrue(in_w(a), (d, a))
                self.assertTrue(in_w(add(a, d)), (d, a))
                self.assertTrue(in_b(add(a, scale(2, d))), (d, a))

    def test_rainbow_exhaustive_small(self):
        for m in range(1, 6):
            T = wwb_triple(m)
            self.assertEqual(sizes(T), wwb_sizes(m))
            self.assertTrue(is_rainbow(*T, m))
            for d in vectors(m):
                a = wwb_witness(d)
                self.assertIn(a, T[0])
                self.assertIn(add(a, d), T[1])
                self.assertIn(add(a, scale(2, d)), T[2])

    def test_tripling_lemma_small(self):
        K = tripled(G1)
        self.assertEqual(len(K), 19)
        self.assertTrue(is_kakeya(K, 3))
        self.assertTrue(is_kakeya(tripled(slices(K33)), 6))
        for m in (1, 2):
            K = tripled(wwb_triple(m))
            self.assertEqual(len(K), wwb_kakeya_size(m))
            self.assertTrue(is_kakeya(K, 3 * m))

    def test_tripled_set_size_formula(self):
        for m in (1, 2, 3):
            self.assertEqual(len(tripled(wwb_triple(m))), wwb_kakeya_size(m))

    def test_tripling_point_exhaustive_in_f3_9(self):
        m = 3
        K = tripled(wwb_triple(m))
        for d in vectors(3 * m):
            x = tripled_point(d, m, wwb_witness)
            self.assertIn(x, K)
            self.assertIn(add(x, d), K)
            self.assertIn(add(x, scale(2, d)), K)

    def test_tripling_point_large(self):
        rng = random.Random(9)
        for m in (4, 7, 20):
            for _ in range(200):
                d = tuple(rng.randrange(3) for _ in range(3 * m))
                x = tripled_point(d, m, wwb_witness)
                for t, roles in ((0, "WWB"), (1, "WBW"), (2, "BWW")):
                    p = add(x, scale(t, d))
                    for k, r in enumerate(roles):
                        block = p[k * m:(k + 1) * m]
                        self.assertTrue(in_w(block) if r == "W" else in_b(block), (d, t, k))

    def test_counts_and_thresholds(self):
        self.assertEqual(wwb_sizes(1), (3, 3, 1))
        self.assertEqual(wwb_sizes(14), (49149, 49149, 1471))
        self.assertEqual(wwb_product(14), 3553383199671)
        self.assertEqual([m for m in range(1, 15) if wwb_product(m) < 8 ** m], [14])
        self.assertLess(wwb_kakeya_size(17), 2 ** 52 - 1)
        self.assertEqual(min(m for m in range(1, 40)
                             if wwb_kakeya_size(m) < 2 ** (3 * m + 1) - 1), 17)
        self.assertLess(wwb_base(14), 2)
        self.assertAlmostEqual(TWO_SYMBOL_BASE, 2 ** (4 / 9) * 3 ** (1 / 3), places=12)
        self.assertLess(abs(wwb_base(3000) - TWO_SYMBOL_BASE), 0.001)

    def test_clean_bound(self):
        # |W|^2 |B| <= 9 * 4^m * 2^(m h(1/3)) = 9 * 432^(m/3), so K(3m, 3) <= 27 * 432^(m/3).
        for m in range(1, 200):
            self.assertLessEqual(wwb_product(m) ** 3, 729 * 432 ** m)
            self.assertLessEqual(wwb_kakeya_size(m) ** 3, 3 ** 9 * 432 ** m)


if __name__ == "__main__":
    unittest.main()
