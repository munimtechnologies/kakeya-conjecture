# test_k7_at_least_132.py -- if this file is in the repository, its test passed: the relaxation of
# sd_bound.py (hyperplane moments + (T) + (SD), see SD.md) rules out every N < 132 in F_3^7, so every
# Kakeya set in F_3^7 has at least 132 points (k_7 >= 132). Bukh-Chao give k_7 >= 103.
# MV2 of Munim, Inc. Made with AI agents working in MV2's name. No human mathematician has checked
# this. Licence: Apache-2.0. Standard library only.

import unittest

from moments import sd_bound as S

X = 132


class K7(unittest.TestCase):
    def test_main_results_hold(self):
        self.assertEqual(S.threshold(3), 13)
        self.assertEqual(S.threshold(5, sd=True), 44)

    def test_k7_at_least_X(self):
        for N in range(1, X):
            self.assertFalse(S.feasible(7, N, sd=True), N)


if __name__ == "__main__":
    unittest.main()
