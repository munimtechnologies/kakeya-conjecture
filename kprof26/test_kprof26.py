# Tests for kprof26.py. MV2 of Munim, Inc. Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0.
import random
import unittest

import kprof26 as k6
from kprof26 import kp


class TestStatic(unittest.TestCase):
    def test_g_identity(self):
        for n in range(19):
            self.assertEqual(k6.g26(n) + k6.g26(18 - n), -408)
        self.assertEqual([k6.g26(n) for n in (6, 7, 8, 13, 14)], [-100, -133, -168, -338, -364])

    def test_hand_lists(self):
        self.assertEqual(k6.pos_cuts(), k6.HAND_POS26)
        self.assertEqual(k6.zero_cuts(), k6.HAND_ZERO26)
        self.assertEqual([k6.delta26(p) for p in k6.HAND_POS26], [12, 9, 5, 2])
        for p in k6.HAND_ZERO26:
            self.assertIn(8, p)
        for p in kp.allowed_profiles(26):
            self.assertEqual(k6.delta26(p), sum(k6.g26(n) for n in p) + 576)
            if p not in k6.HAND_POS26 and p not in k6.HAND_ZERO26:
                self.assertLess(k6.delta26(p), 0)
            if 8 in p:
                self.assertEqual(k6.delta26(p), 0)

    def test_delta_identity_random(self):
        rng = random.Random(7)
        for _ in range(20):
            K = rng.sample(kp.PTS4, 26)
            L = kp.count_lines(K)
            s = sum(k6.delta26(kp.class_profile(K, u)) for u in kp.DIR4)
            self.assertEqual(s, 9 * L - 360)

    def test_cases_and_tables(self):
        names = [c[0] for c in k6.CASES26]
        self.assertEqual(names, ["P1", "P2", "P3", "P4", "Z1", "Z2", "Z3", "Z4", "Z5"])
        self.assertEqual(k6.CASES26[3][2], [(6, 6, 14), (5, 7, 14), (6, 7, 13)])
        self.assertEqual(k6.CASES26[8][2], k6.HAND_ZERO26[:4])
        self.assertEqual(len(k6.zero_table()), 24)
        self.assertEqual(len(k6.zero_table([(5, 8, 13)])), 18)
        for row in k6.zero_table():
            self.assertEqual(sum(row[:3]), 26)
            self.assertEqual(kp.moments_of(row[:3]), (row[3], row[4]))
        for name, cut, forbid, zero_only in k6.CASES26:
            self.assertEqual(sum(cut), 26)
            self.assertEqual(list(cut), sorted(cut))
            self.assertNotIn(cut, forbid)

    def test_lines_per_direction(self):
        cnt = k6.lines_per_direction(kp.KNOWN27)
        self.assertEqual(sum(cnt), kp.count_lines(kp.KNOWN27))
        self.assertGreaterEqual(min(cnt), 1)


class TestSolver(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.perms = kp.agl3_perms()

    def test_extras_sanity(self):
        out = k6.extras_sanity(self.perms, 2, 60)
        self.assertTrue(out["ok"], out)

    def test_removals_agree(self):
        out = k6.removal_agreement(self.perms, 2, 30, which=[0, 9, 13, 26])
        self.assertTrue(out["ok"], out)

    def test_impossible_sizes(self):
        # sizes (3, 9, 14) of slice 0 / 1 / 2 are not an allowed cut for 26: infeasible
        spec = k6.spec_for((3, 9, 14), [], False, [0, 1, 2])
        self.assertEqual(k6.solve_cpsat26(spec, 60, 2)[0], "INFEASIBLE")


if __name__ == "__main__":
    unittest.main()
