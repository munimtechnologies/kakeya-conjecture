# test_min_section.py -- checks the second-moment identity and the small-section lemma
# (min_section.py). MV2 of Munim, Inc. Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0. Standard library only.

import random
import unittest

from moments import min_section as M


class MinSection(unittest.TestCase):
    def test_counts(self):
        self.assertEqual(len(M.points()), 81)
        self.assertEqual(len(M.directions()), 40)

    def test_identity_on_random_sets(self):
        rng = random.Random(2026)
        P = M.points()
        for N in (1, 2, 3, 13, 25, 26, 27, 53, 81):
            for _ in range(4):
                K = rng.sample(P, N)
                secs = M.section_sizes(K)
                self.assertTrue(all(sum(t) == N for t in secs))
                self.assertEqual(sum(a * a for t in secs for a in t), M.second_moment(N))

    def test_pairs_lie_in_13_hyperplanes(self):
        rng = random.Random(7)
        P = M.points()
        D = M.directions()
        for _ in range(30):
            x, y = rng.sample(P, 2)
            n = sum(1 for u in D
                    if sum(a * b for a, b in zip(u, x)) % 3 == sum(a * b for a, b in zip(u, y)) % 3)
            self.assertEqual(n, 13)

    def test_hand_values(self):
        self.assertEqual(M.second_moment(25), 8800)
        self.assertEqual(M.second_moment(26), 9490)
        self.assertEqual(M.nmax_for_size(25), 6)
        self.assertEqual(M.nmax_for_size(26), 7)

    def test_lemma_on_random_sets(self):
        rng = random.Random(11)
        P = M.points()
        for N in (24, 25, 26, 27):
            for _ in range(100):
                K = rng.sample(P, N)
                self.assertLessEqual(min(min(t) for t in M.section_sizes(K)), M.nmax_for_size(N))


if __name__ == "__main__":
    unittest.main()
