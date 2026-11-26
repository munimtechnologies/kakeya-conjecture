# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Tests for k8/census.py.

The F_3^3 tests check exhaustively the hand proof in k8/CENSUS.md: k_2 = 7, every 7-point
subset of F_3^2 is Kakeya, and no layered Kakeya set of F_3^3 has 13 points or fewer (so the
13-point Kakeya sets of F_3^3 have no Kakeya plane section)."""
import itertools
import unittest

try:
    from . import census, check, layer
except ImportError:  # pragma: no cover
    import census
    import check
    import layer


class CensusTests(unittest.TestCase):
    def test_plane_kakeya_sets(self):
        pts = range(9)
        self.assertFalse(any(check.is_kakeya(2, list(S))
                             for S in itertools.combinations(pts, 6)))
        self.assertTrue(all(check.is_kakeya(2, list(S))
                            for S in itertools.combinations(pts, 7)))

    def test_layered_f33_needs_14(self):
        # |L0| = 7 needs extra >= 7, |L0| = 8 needs extra >= 6, |L0| = 9 needs extra >= 5
        for size in (7, 8, 9):
            for L0 in itertools.combinations(range(9), size):
                self.assertIsNone(census.exact_min_extra(2, list(L0), 13 - size), L0)

    def test_layered_f33_exact(self):
        ex = census.exact_min_extra(2, layer.K7_PLANE, 8)
        self.assertIn(ex, (7, 8))
        self.assertEqual(census.exact_min_extra(2, layer.classical(2), 8), ex)
        W = layer.Layer(2, layer.K7_PLANE)
        r = layer.search(W, seed=1, budget=20.0, max_steps=3000)
        self.assertTrue(check.is_kakeya(3, W.full_indices(r["points"])))
        self.assertGreaterEqual(r["best"], ex)
        print("F_3^3: least layered Kakeya set has %d points (search found %d)" % (
            7 + ex, 7 + r["best"]))

    def test_census_phase_small(self):
        ph = census.census_phase(3, seconds=10.0, per=0.5, max_add=1, record=27, seed0=1,
                                 walk_budget=0.3)
        self.assertEqual(ph["unverified"], 0)
        self.assertGreaterEqual(ph["distinct_L0"], 1)
        self.assertTrue(ph["stats"])
        self.assertGreaterEqual(ph["best"][0], 27)   # k_4 = 27 (k4 proof in this repo)
        self.assertEqual(ph["new"], [])
        print("census F_3^4: " + census._fmt_stats(ph["stats"]))


if __name__ == "__main__":
    unittest.main()
