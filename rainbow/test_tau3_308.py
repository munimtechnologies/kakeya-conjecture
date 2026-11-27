# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Every rainbow triple (A, M, C) of F_3^3 has |A| |M| |C| >= 308 (rainbow/RAINBOW3.md).

Since (3 / 2^(2/3))^9 = 307.55..., no rainbow triple of F_3^3 gives, through the tripling
lemma of TRIPLES.md, a Kakeya growth base over F_3 below the 2010 bound 3 / 2^(2/3) =
1.88988... of Kopparty, Lev, Saraf and Sudan.

Method (SCAN.md): put the triple in normal form (|A| <= |M| <= |C|, A a canonical affine
class representative, 0 in M and M least among its translates).  A product below 308 forces
|A| <= 6 (7^3 = 343), and by the hyperplane lemma H a smallest set inside a plane forces a
product >= 9 tau(2) = 576, so A spans F_3^3 and has 4, 5 or 6 points.  For every such class,
every allowed size of M and every M in normal form, the least C is found exactly; none
gives a product below 308."""
import os
import sys
import time
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

try:
    from . import rscan
except ImportError:  # pragma: no cover
    import rscan


class TestTau3AtLeast308(unittest.TestCase):
    def test_no_rainbow_triple_of_f3_cubed_below_308(self):
        # the sizes the normal form allows below 308
        self.assertEqual(rscan.sizes_for(3, 4, 308), [7, 8])
        self.assertEqual(rscan.sizes_for(3, 5, 308), [6, 7])
        self.assertEqual(rscan.sizes_for(3, 6, 308), [6, 7])
        self.assertEqual(rscan.sizes_for(3, 7, 308), [])
        deadline = time.time() + 420
        res = rscan.sweep(3, 308, (4, 5, 6), deadline, workers=8, full_only=True)
        by = rscan.summarise(res)
        msg = rscan.describe(by)
        hits = [r for r in res if r["best"] is not None]
        self.assertFalse(hits, "TRIPLE FOUND " + repr(
            [(r["best"], r["A"], r["M"], r["C"]) for r in hits[:3]]) + " || " + msg)
        self.assertTrue(all(r["complete"] for r in res), "INCOMPLETE || " + msg)
        self.assertEqual(sorted(by), [4, 5, 6], msg)
        full = {x: sum(1 for k, _, _ in rscan.classes(3, x) if k == 3) for x in (4, 5, 6)}
        for x in (4, 5, 6):
            self.assertEqual(by[x]["reps"], full[x], msg)
            self.assertGreater(by[x]["nodes"], 0, msg)


if __name__ == "__main__":
    unittest.main()
