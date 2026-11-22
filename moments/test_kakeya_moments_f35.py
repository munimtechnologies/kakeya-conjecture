# F_3^5 check for kakeya_moments.py. MV2 of Munim, Inc. Made with AI agents. Apache-2.0.
# This file is committed only if it passes: then every N <= 40 is excluded by the exact
# integer test in kakeya_moments.feasible(), i.e. every Kakeya set in F_3^5 has >= 41 points.
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import kakeya_moments as km  # noqa: E402


class TestF35(unittest.TestCase):
    def test_counts(self):
        # 121 hyperplane classes, 40 hyperplanes through a line, 13 through a plane
        self.assertEqual((km.gauss1(3, 5), km.gauss1(3, 4), km.gauss1(3, 3)), (121, 40, 13))
        self.assertEqual(km.targets(3, 5, 40), (121, 31200, 128440 + 27 * 121))

    def test_certificate_40(self):
        # with c = 13.9 the best profile is (13, 13, 14): 10*936 - 139*247 = -24973
        c = km.certificate(3, 5, 40, 139, 10)
        self.assertEqual(c["max_profile_value"], -24973)
        self.assertLess(c["collinear_triples_at_most"], 121)

    def test_bound_at_least_41(self):
        for N in range(1, 41):
            self.assertFalse(km.feasible(3, 5, N), N)
        self.assertGreaterEqual(km.lower_bound(3, 5), 41)


if __name__ == "__main__":
    unittest.main()
