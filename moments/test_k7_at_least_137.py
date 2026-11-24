# test_k7_at_least_137.py -- if this file is in the repository, its test passed: the relaxation of
# sd_bound.py (hyperplane moments + (T) + (SD), see SD.md) rules out N = 136 in F_3^7. Every N < 136 is
# ruled out by test_k7_at_least_136.py (already in the repository), so every Kakeya set in F_3^7 has
# at least 137 points (k_7 >= 137). Bukh-Chao give k_7 >= 103.
# MV2 of Munim, Inc. Made with AI agents working in MV2's name. No human mathematician has checked
# this. Licence: Apache-2.0. Standard library only.

import unittest

from moments import sd_bound as S


class K7(unittest.TestCase):
    def test_n136_infeasible(self):
        self.assertFalse(S.feasible(7, 136, sd=True))


if __name__ == "__main__":
    unittest.main()
