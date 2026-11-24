# indep26/test_cut_6_6_14.py - independent check: no 26-point Kakeya set in F_3^4 has a hyperplane class with cut (6, 6, 14).
# MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name;
# no human mathematician has checked it yet. Licence: Apache-2.0. Python standard library only.
# This file is committed only if the exhaustive search below finds nothing (see indep26/README.md).
import unittest
from math import comb

from indep26 import kk

CUT = (6, 6, 14)
PART, PARTS = 0, 1


class TestCut(unittest.TestCase):
    def test_no_26_point_kakeya_set_with_this_cut(self):
        self.assertTrue(kk.quick_controls())
        a, b, c = CUT
        self.assertEqual(a + b + c, 26)
        self.assertIn(CUT, kk.CUTS26)
        S = kk.space(3)
        reps = kk.orbit_reps(S, a)
        self.assertEqual(sum(o for _, o in reps), comb(27, a))
        rows = kk.search_cut(S, CUT, PART, PARTS)
        self.assertEqual([r["rep"] for r in rows], [i for i in range(len(reps)) if i % PARTS == PART])
        hits = [r for r in rows if r["found"] is not None]
        self.assertEqual(hits, [], "candidate found for cut %r: %r" % (CUT, hits))


if __name__ == "__main__":
    unittest.main()
