# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Smoke test of walk.py: in F_3^3 the search reaches a 13-point Kakeya set (13 is optimal,
Lev 2009), the independent checker accepts it, and rejects a set with one point removed."""
import unittest

try:
    from . import walk
except ImportError:
    import walk


class WalkN3(unittest.TestCase):
    def test_tables(self):
        T = walk.setup(3)
        self.assertEqual(T["D"], 13)
        for j in range(T["D"]):
            self.assertEqual(len(T["lines"][j]), 9)
            for ln in T["lines"][j]:
                s = [sum(T["pts"][p][i] for p in ln) % 3 for i in range(3)]
                self.assertEqual(s, [0, 0, 0])

    def test_reaches_13(self):
        r = walk.run(3, "trivial", seed=1, budget=30.0)
        self.assertEqual(r["best"], 13, walk.report([r]))
        self.assertTrue(r["verified"])
        T = walk.setup(3)
        K = [T["pts"][p] for p in r["points"]]
        self.assertTrue(walk.is_kakeya(3, K))
        for i in range(len(K)):
            self.assertFalse(walk.is_kakeya(3, K[:i] + K[i + 1:]))


if __name__ == "__main__":
    unittest.main()
