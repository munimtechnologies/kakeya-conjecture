# test_product.py - checks for moments/PRODUCT.md (standard library only).
# MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name;
# no human mathematician has checked it yet. Apache-2.0.
import itertools
import random
import unittest
from math import comb

CUTS26 = [(6, 6, 14), (5, 7, 14), (6, 7, 13), (5, 8, 13), (7, 7, 12), (6, 8, 12), (7, 8, 11)]


def normals(n):
    """One normal vector per parallel class of hyperplanes of F_3^n ((3^n - 1) / 2 of them)."""
    out = []
    for v in itertools.product(range(3), repeat=n):
        nz = [x for x in v if x]
        if nz and nz[0] == 1:
            out.append(v)
    return out


def section_sizes(K, xi):
    s = [0, 0, 0]
    for p in K:
        s[sum(a * b for a, b in zip(p, xi)) % 3] += 1
    return s


def n_lines(K):
    Ks = set(K)
    cnt = 0
    for x in K:
        for y in K:
            if x != y and tuple((-a - b) % 3 for a, b in zip(x, y)) in Ks:
                cnt += 1
    assert cnt % 6 == 0
    return cnt // 6


def emax(x, y, T=27):
    """Largest sum of r_d^2 over T integers 1 <= r_d <= min(x, y) with sum x*y (greedy = optimal,
    since a convex function on a box-and-sum polytope peaks at a vertex)."""
    m = min(x, y)
    extra = x * y - T
    assert 0 <= extra <= T * (m - 1)
    if m == 1:
        return T
    k, rem = divmod(extra, m - 1)
    total = k * m * m + (T - k) * 1
    if rem:
        total += (1 + rem) ** 2 - 1
    return total


def emax_brute(x, y, T):
    m = min(x, y)
    best = -1

    def rec(i, left, acc):
        nonlocal best
        if i == T:
            if left == 0:
                best = max(best, acc)
            return
        for r in range(1, m + 1):
            if left - r < T - i - 1 or left - r > (T - i - 1) * m:
                continue
            rec(i + 1, left - r, acc + r * r)

    rec(0, x * y, 0)
    return best


def allowed(cut, T=27):
    """(T): every two sections have product >= T; (SD): every section >= T^2 / Emax(other two)."""
    for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
        a, b, c = cut[i], cut[j], cut[k]
        if a * b < T:
            return False
        if c * emax(a, b, T) < T * T:
            return False
    return True


class TestProductIdentity(unittest.TestCase):
    def test_identities_on_random_sets(self):
        rng = random.Random(11)
        for n in (2, 3, 4):
            pts = list(itertools.product(range(3), repeat=n))
            xis = normals(n)
            self.assertEqual(len(xis), (3 ** n - 1) // 2)
            for _ in range(12 if n < 4 else 8):
                N = rng.randint(3, len(pts))
                K = rng.sample(pts, N)
                L = n_lines(K)
                cuts = [section_sizes(K, xi) for xi in xis]
                prod = sum(a * b * c for a, b, c in cuts)
                e2 = sum(a * b + a * c + b * c for a, b, c in cuts)
                s2 = sum(comb(a, 2) + comb(b, 2) + comb(c, 2) for a, b, c in cuts)
                self.assertEqual(prod, 3 ** (n - 2) * (comb(N, 3) + 2 * L) if n >= 2 else prod)
                self.assertEqual(e2, 3 ** (n - 1) * comb(N, 2))
                self.assertEqual(s2, (3 ** (n - 1) - 1) // 2 * comb(N, 2))

    def test_known27(self):
        K27 = [(0, 0, 0, 0), (1, 0, 0, 0), (2, 0, 0, 0), (2, 2, 0, 0), (0, 0, 1, 0), (1, 0, 1, 0),
               (0, 1, 1, 0), (1, 1, 1, 0), (0, 2, 1, 0), (2, 2, 1, 0), (0, 0, 2, 0), (0, 1, 2, 0),
               (1, 2, 2, 0), (1, 2, 0, 1), (1, 0, 1, 1), (0, 1, 1, 1), (1, 2, 1, 1), (0, 0, 2, 1),
               (2, 1, 2, 1), (2, 0, 0, 2), (0, 1, 0, 2), (0, 2, 0, 2), (1, 0, 1, 2), (2, 0, 1, 2),
               (1, 2, 1, 2), (0, 1, 2, 2), (2, 1, 2, 2)]
        L = n_lines(K27)
        cuts = [section_sizes(K27, xi) for xi in normals(4)]
        self.assertEqual(sum(a * b * c for a, b, c in cuts), 9 * (comb(27, 3) + 2 * L))
        # N = 27 form: sum of (n1 - 9)(n2 - 9)(n3 - 9) = 18 (L - 36)
        self.assertEqual(sum((a - 9) * (b - 9) * (c - 9) for a, b, c in cuts), 18 * (L - 36))
        self.assertGreaterEqual(L, 40)

    def test_shifted_forms(self):
        """sum (n_i - 8) products: N = 26 gives 18 (L - 40); N = 25 gives 18 L - 580."""
        rng = random.Random(12)
        pts = list(itertools.product(range(3), repeat=4))
        xis = normals(4)
        for N, const in ((26, -720), (25, -580)):
            for _ in range(4):
                K = rng.sample(pts, N)
                L = n_lines(K)
                cuts = [section_sizes(K, xi) for xi in xis]
                self.assertEqual(sum((a - 8) * (b - 8) * (c - 8) for a, b, c in cuts), 18 * L + const)


class TestCutLists(unittest.TestCase):
    def test_emax_greedy_is_optimal(self):
        for T in (5, 6):
            for x in range(1, 6):
                for y in range(x, 7):
                    if T <= x * y <= T * x:
                        self.assertEqual(emax(x, y, T), emax_brute(x, y, T), (x, y, T))

    def test_sd_values(self):
        self.assertEqual(emax(5, 6), 42)
        self.assertEqual(emax(4, 8), 50)   # so (4, 8) forces a third section of at least 15
        self.assertFalse(allowed((4, 8, 14)))
        self.assertFalse(allowed((4, 7, 15)))
        self.assertFalse(allowed((5, 6, 15)))

    def test_cuts26_complete(self):
        """Every cut a 26-point Kakeya set could need, by the product identity, is in CUTS26."""
        pos, zero_small = set(), set()
        for a in range(27):
            for b in range(a, 27):
                c = 26 - a - b
                if c < b:
                    continue
                cut = (a, b, c)
                if not allowed(cut):
                    continue
                p = (a - 8) * (b - 8) * (c - 8)
                if p > 0:
                    pos.add(cut)
                elif p == 0 and a <= 7:
                    zero_small.add(cut)
        self.assertTrue(pos <= set(CUTS26), pos - set(CUTS26))
        self.assertTrue(zero_small <= set(CUTS26), zero_small - set(CUTS26))
        self.assertEqual(pos, {(6, 6, 14), (5, 7, 14), (6, 7, 13), (7, 7, 12)})
        self.assertEqual(zero_small, {(5, 8, 13), (6, 8, 12), (7, 8, 11)})

    def test_min_section_26(self):
        """Summed over the 40 classes, n1^2 + n2^2 + n3^2 = 40 N + 13 N (N - 1). For N = 26 that is
        9490, but if every section had at least 8 points each class would be (8,8,10) or (8,9,9),
        with at most 228, and 40 * 228 = 9120 < 9490. So some hyperplane has at most 7 points."""
        N = 26
        total = 40 * N + 13 * N * (N - 1)
        self.assertEqual(total, 9490)
        big = [(a, b, N - a - b) for a in range(8, N) for b in range(a, N) if N - a - b >= b]
        self.assertEqual(sorted(big), [(8, 8, 10), (8, 9, 9)])
        self.assertLess(40 * max(a * a + b * b + c * c for a, b, c in big), total)


if __name__ == "__main__":
    unittest.main()
