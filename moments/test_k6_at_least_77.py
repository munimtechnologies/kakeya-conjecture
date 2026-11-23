# test_k6_at_least_77.py -- if this file is in the repository, its test passed: the relaxation of
# sd_bound.py (hyperplane moments + (T) + (SD), see SD.md) rules out every N < 77 in F_3^6, so every
# Kakeya set in F_3^6 has at least 77 points (k_6 >= 77). Bukh-Chao give k_6 >= 57.
# MV2 of Munim, Inc. Made with AI agents working in MV2's name. No human mathematician has checked
# this. Licence: Apache-2.0. Standard library only.

import unittest

from moments import sd_bound as S

X = 77


class K6(unittest.TestCase):
    def test_main_results_hold(self):
        self.assertEqual(S.threshold(3), 13)
        self.assertEqual(S.threshold(4, sd=True, hz=True), 25)
        self.assertEqual(S.threshold(5, sd=True), 44)

    def test_k6_at_least_X(self):
        for N in range(1, X):
            self.assertFalse(S.feasible(6, N, sd=True), N)


if __name__ == "__main__":
    unittest.main()
