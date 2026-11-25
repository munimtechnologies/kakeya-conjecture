# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""The numbers quoted in SD2.md 'Result (k5cut-v1)' are the ones in results/k5cut-v1.json."""
import json
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))


class K5cutResult(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(HERE, "results", "k5cut-v1.json"), encoding="utf-8") as f:
            self.r = json.load(f)

    def test_verdict(self):
        m = self.r["main"]
        self.assertEqual(self.r["verdict"], "CLOSED")
        self.assertEqual(m["verdict"], "CLOSED")
        self.assertEqual((m["x"], m["y"], m["z"], m["d"]), (10, 11, 23, 4))
        self.assertEqual(m["n_solutions"], 0)
        self.assertEqual(m["solutions"], [])
        self.assertEqual(m["Emin"], 315)

    def test_table(self):
        p = self.r["main"]["per_k1"]
        want = {"8": (11, 75, 197100, 129367), "9": (10, 162, 11664, 3159), "10": (10, 20, 20, 20)}
        self.assertEqual(sorted(p), sorted(want))
        for k, (xm, reps, xa, kept) in want.items():
            row = p[k]
            self.assertEqual((row["Xmax"], row["reps"], row["XA"], row["XA_kept"]), (xm, reps, xa, kept))
            self.assertEqual(row["B"], 0)
            self.assertEqual(row["solutions"], 0)

    def test_controls(self):
        self.assertTrue(all(self.r["orbit_sum_checks"][k] for k in ("1", "2", "3", "4", "5")))
        c = self.r["control_lev"]
        self.assertEqual(c["sizes"], [6, 8, 13])
        self.assertTrue(c["star"] and c["maxcol_allowed"] and c["I_excess_ok"])
        self.assertEqual((c["E"], c["Emin"], c["maxcol"], c["I_excess"], c["min_cover"]),
                         (110, 78, 5, 2, 13))


if __name__ == "__main__":
    unittest.main()
