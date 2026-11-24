# test_plane_n25_not_excluded.py -- if this file is in the repository, its test passed: the
# plane-class relaxation of plane_lp.py (PLANE.md) DOES admit N = 25 (with at least 40 lines), so it
# cannot by itself rule out a 25-point Kakeya set in F_3^4. An honest negative readout; the sharper
# relaxation is plane_lp2.py (PLANE2.md).
# MV2 of Munim, Inc. Made with AI agents working in MV2's name. No human mathematician has checked
# this. Licence: Apache-2.0. Standard library only.

import unittest

from moments import plane_lp as P


class K4(unittest.TestCase):
    def test_n25_not_excluded(self):
        self.assertTrue(P.plane_feasible(25))


if __name__ == "__main__":
    unittest.main()
