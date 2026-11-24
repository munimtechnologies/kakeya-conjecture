# test_bad25.py -- checks every claim of BAD25.md (bad25.py) in exact integer arithmetic.
# MV2 of Munim, Inc. Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0. Standard library only.

import random
import unittest
from itertools import combinations_with_replacement
from math import comb

from moments import bad25 as B

DELTA = {(4, 9, 12): -50, (4, 10, 11): -81, (5, 7, 13): 60, (5, 8, 12): 0, (5, 9, 11): -36,
         (5, 10, 10): -48, (6, 6, 13): 76, (6, 7, 12): 25, (6, 8, 11): -9, (6, 9, 10): -26,
         (7, 7, 11): 0, (7, 8, 10): -15, (7, 9, 9): -20, (8, 8, 9): -18}


class Bad25(unittest.TestCase):
    def test_profile_list(self):
        self.assertEqual(B.allowed_profiles(25), B.HAND25)
        # (T) alone would also allow these six; (SD) removes them
        for p in [(3, 9, 13), (3, 10, 12), (3, 11, 11), (4, 7, 14), (4, 8, 13), (5, 6, 14)]:
            self.assertGreaterEqual(min(p[0] * p[1], p[0] * p[2], p[1] * p[2]), 27)
            self.assertNotIn(p, B.allowed_profiles(25))
        # agreement with SD.md at N = 24: 12 cuts pass (T) and (SD), (8,8,8) among them
        a24 = B.allowed_profiles(24)
        self.assertEqual(len(a24), 12)
        self.assertIn((8, 8, 8), a24)
        self.assertNotIn((5, 6, 13), a24)

    def test_delta_values(self):
        for p in B.HAND25:
            self.assertEqual(B.delta(p), DELTA[p], p)
        pos = sorted(p for p in B.HAND25 if B.delta(p) > 0)
        self.assertEqual(pos, sorted(B.BAD25))
        self.assertEqual(sorted(p for p in B.HAND25 if B.delta(p) == 0), [(5, 8, 12), (7, 7, 11)])

    def test_identity_algebra(self):
        for L in range(0, 1081):
            lhs = 7 * (4 * comb(25, 3) + 9 * L) - 51 * 13 * comb(25, 2) + 3302 * 40
            self.assertEqual(lhs, 63 * L - 2420)

    def test_identity_on_random_sets(self):
        rng = random.Random(25)
        P = B.points()
        self.assertEqual(len(B.directions()), 40)
        for trial in range(40):
            K = rng.sample(P, 25)
            if trial % 4 == 0:  # add structure: force some lines in
                K = list(set(K[:19]) | {(0, 0, 0, t) for t in range(3)} | {(t, 1, 2, 0) for t in range(3)})
                extra = [p for p in P if p not in K]
                K += rng.sample(extra, 25 - len(K))
            self.assertEqual(len(set(K)), 25)
            cs = B.cuts(K)
            L = B.count_lines(K)
            s2 = sum(B.moments_of(c)[0] for c in cs)
            s3 = sum(B.moments_of(c)[1] for c in cs)
            self.assertEqual(s2, 3900)
            self.assertEqual(s3, 9200 + 9 * L)
            self.assertEqual(sum(B.delta(c) for c in cs), 63 * L - 2420)

    def test_line_count_brute(self):
        P = B.points()
        self.assertEqual(B.count_lines(P), 1080)
        self.assertEqual(B.count_lines([(0, 0, 0, 0), (1, 1, 0, 0), (2, 2, 0, 0)]), 1)

    def test_consequences(self):
        # a Kakeya set has L >= 40, so the delta sum is >= 100; one bad class gives at most 76
        self.assertEqual(63 * 40 - 2420, 100)
        self.assertLess(max(B.delta(p) for p in B.BAD25), 100)
        # two classes are enough only for these pairs
        ok_pairs = sorted(pq for pq in combinations_with_replacement(sorted(B.BAD25), 2)
                          if B.delta(pq[0]) + B.delta(pq[1]) >= 100)
        self.assertEqual(ok_pairs, [((5, 7, 13), (5, 7, 13)), ((5, 7, 13), (6, 6, 13)),
                                    ((6, 6, 13), (6, 6, 13)), ((6, 6, 13), (6, 7, 12))])
        # without a 13-point section only (6,7,12) is positive: at least four such classes
        self.assertLess(3 * B.delta((6, 7, 12)), 100)
        self.assertGreaterEqual(4 * B.delta((6, 7, 12)), 100)
        # 13 is the largest section; 13-point sections occur only in (5,7,13) and (6,6,13)
        self.assertEqual(max(p[2] for p in B.HAND25), 13)
        self.assertEqual(sorted(p for p in B.HAND25 if p[2] == 13), [(5, 7, 13), (6, 6, 13)])
        # with all other classes at delta <= 0, the line count is at most (2420 + positive part)/63
        self.assertEqual((2420 + 40 * 76) // 63, 86)

    def test_n26(self):
        self.assertEqual(B.allowed_profiles(26), B.HAND26)
        d26 = {p: B.delta26(p) for p in B.HAND26}
        self.assertEqual(sorted(p for p in B.HAND26 if d26[p] > 0), sorted(B.BAD26))
        self.assertEqual([d26[p] for p in B.BAD26], [2, 5, 9, 12])
        zero = sorted(p for p in B.HAND26 if d26[p] == 0)
        self.assertEqual(zero, [(5, 8, 13), (6, 8, 12), (7, 8, 11), (8, 8, 10), (8, 9, 9)])
        # the zero cuts are exactly the allowed cuts with an 8-point section
        self.assertEqual(zero, sorted(p for p in B.HAND26 if 8 in p))
        for L in range(0, 1081):
            lhs = (4 * comb(26, 3) + 9 * L) - 8 * 13 * comb(26, 2) + 576 * 40
            self.assertEqual(lhs, 9 * (L - 40))
        rng = random.Random(26)
        P = B.points()
        for _ in range(20):
            K = rng.sample(P, 26)
            self.assertEqual(sum(B.delta26(c) for c in B.cuts(K)), 9 * (B.count_lines(K) - 40))

    def test_known_27_set(self):
        K = B.KNOWN27
        self.assertEqual(len(set(K)), 27)
        self.assertTrue(B.is_kakeya(K))
        allowed27 = set(B.allowed_profiles(27))
        cs = B.cuts(K)
        for c in cs:
            self.assertIn(tuple(sorted(c)), allowed27)
        L = B.count_lines(K)
        self.assertGreaterEqual(L, 40)
        self.assertEqual(sum(B.moments_of(c)[1] for c in cs), 4 * comb(27, 3) + 9 * L)
        self.assertEqual(sum(B.moments_of(c)[0] for c in cs), 13 * comb(27, 2))


if __name__ == "__main__":
    unittest.main()
