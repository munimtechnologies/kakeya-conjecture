# test_plane2_n25_not_excluded.py -- if this file is in the repository, its test passed: the sharper
# plane-class relaxation of plane_lp2.py (PLANE2.md) DOES admit N = 25 (with at least 40 lines), so it
# cannot by itself rule out a 25-point Kakeya set in F_3^4. An honest negative readout.
# MV2 of Munim, Inc. Made with AI agents working in MV2's name. No human mathematician has checked
# this. Licence: Apache-2.0. Standard library only.

import unittest

from moments import plane_lp2 as Q


class K4(unittest.TestCase):
    def test_n25_not_excluded(self):
        self.assertTrue(Q.plane2_feasible(25))


if __name__ == "__main__":
    unittest.main()
