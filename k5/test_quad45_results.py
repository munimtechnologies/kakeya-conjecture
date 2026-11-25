# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Every number in k5/results/quad45.json, recomputed exactly from sd2.py and quad.py."""
import json
import os
import unittest
from math import comb

try:
    from . import quad, sd2
except ImportError:
    import quad
    import sd2

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    with open(os.path.join(HERE, "results", "quad45.json")) as f:
        return json.load(f)


class Quad45Results(unittest.TestCase):
    def setUp(self):
        self.R = load()

    def test_constants_and_identities(self):
        R = self.R
        T, H, a1, a2 = sd2.constants(5)
        self.assertEqual((R["T"], R["classes"], R["a1"], R["a2"]), (T, H, a1, a2))
        N = R["N"]
        self.assertEqual(R["pair_total"], a1 * comb(N, 2))
        base, heavy = tuple(R["base_profile"]), tuple(R["heavy_profile"])
        self.assertEqual(R["base_s2"], sd2.s2(base))
        self.assertEqual(R["base_s3"], sd2.s3(base))
        self.assertEqual(R["sum_delta2"], a1 * comb(N, 2) - H * sd2.s2(base))
        # sum f = sum Delta3 - 15 sum Delta2 = 27 (L - 121): constant part vanishes
        self.assertEqual(a2 * comb(N, 3) + (a1 - a2) * H - H * sd2.s3(base)
                         - 15 * R["sum_delta2"], 0)
        self.assertEqual(a1 - a2, 27)
        d2 = sd2.s2(heavy) - sd2.s2(base)
        f = sd2.s3(heavy) - sd2.s3(base) - 15 * d2
        self.assertEqual((R["heavy_delta2"], R["heavy_f"]), (d2, f))
        self.assertEqual(R["max_heavy_classes"], R["sum_delta2"] // d2)

    def test_profile_classes(self):
        R = self.R
        N, T = R["N"], R["T"]
        base, heavy = tuple(R["base_profile"]), tuple(R["heavy_profile"])
        profs = sd2.profiles(N, T, "sd2")
        self.assertIn(base, profs)
        self.assertIn(heavy, profs)
        for p in profs:
            d2 = sd2.s2(p) - sd2.s2(base)
            f = sd2.s3(p) - sd2.s3(base) - 15 * d2
            if p == base:
                self.assertEqual((d2, f), (0, 0))
            elif p == heavy:
                self.assertEqual((d2, f), (52, 2))
            else:
                self.assertGreaterEqual(d2, 1, p)
                self.assertLessEqual(f, -1, p)
                if d2 % 3 != 1:
                    self.assertLessEqual(f, -3, p)
            # Delta2 = -Q_K (mod 3) because C(44,2) - 301 = 645 is divisible by 3
            self.assertEqual((-d2) % 3, quad.residue(p, N))

    def test_quadric_table_and_case_split(self):
        R = self.R
        table = {name: list(z) for name, z in quad.quadric_types(5)}
        self.assertEqual(table, R["quadric_types"])
        small = sorted(k for k, z in table.items() if z[1] <= R["max_z1"])
        self.assertEqual(small, sorted(R["types_with_z1_at_most_18"]))
        self.assertEqual(R["max_z1"], (2 * R["max_heavy_classes"]) // 3)
        # rank1-a2: base classes only where Q_K = 0
        z0 = table["rank1-a2"][0]
        k = R["max_heavy_classes"]
        self.assertGreaterEqual(2 * k, R["classes"] - k - z0)          # forces k = 27
        self.assertLess(2 * (k - 1), R["classes"] - (k - 1) - z0)
        self.assertLess(R["sum_delta2"] - R["heavy_delta2"] * k, R["classes"] - k - z0)

    def test_lp_at_44_and_45(self):
        R = self.R
        ok44, per44 = quad.excluded(5, 44)
        self.assertTrue(ok44)
        self.assertEqual(sum(v["excluded"] for v in per44.values()), R["lp_types_excluded_at_44"])
        ok45, per45 = quad.excluded(5, 45)
        self.assertFalse(ok45)
        self.assertEqual(sum(v["excluded"] for v in per45.values()), R["lp_types_excluded_at_45"])
        self.assertEqual(sorted(k for k, v in per45.items() if not v["excluded"]),
                         sorted(R["lp_types_surviving_at_45"]))
        for name, v in per44.items():
            if v["mu"] is not None:
                self.assertTrue(quad.check_certificate(5, 44, v["z"], v["mu"]), name)

    def test_bounds(self):
        b = self.R["bounds"]
        self.assertEqual(b["k5_lower_this_result"], self.R["N"] + 1)
        self.assertLess(b["k5_lower_published_general"], b["k5_lower_this_result"])
        self.assertLess(b["k5_lower_this_result"], b["k5_upper_chung_du_wesley"])


if __name__ == "__main__":
    unittest.main()
