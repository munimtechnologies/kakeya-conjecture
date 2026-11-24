# gate26/test_cut_5_8_13.py - no 26-point Kakeya set in F_3^4 has a hyperplane class with cut (5, 8, 13).
# MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name;
# no human mathematician has checked it yet. Licence: Apache-2.0. Python standard library only.
# This file is committed only if the exhaustive search below finds nothing (see gate26/README.md).
import unittest
from math import comb

from gate26 import engine as E

CUT = (5, 8, 13)


class TestCut(unittest.TestCase):
    def test_no_26_point_kakeya_set_with_this_cut(self):
        a, b, c = CUT
        self.assertEqual(a + b + c, 26)
        self.assertIn(CUT, E.K.CUTS26)
        reps = E.reps(a)
        self.assertEqual(sum(o for _, o in reps), comb(27, a))
        rows = E.close_cut(CUT)
        self.assertEqual(len(rows), len(reps))
        hits = [r for r in rows if r["found"] is not None]
        self.assertEqual(hits, [], "candidate found for cut %r: %r" % (CUT, hits))


if __name__ == "__main__":
    unittest.main()
