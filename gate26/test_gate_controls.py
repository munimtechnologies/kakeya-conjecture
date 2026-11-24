# gate26/test_gate_controls.py - controls that every gate26 commit runs before its cut check.
# MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name;
# no human mathematician has checked it yet. Licence: Apache-2.0. Python standard library only.
import itertools
import unittest
from math import comb

from gate26 import engine as E

K = E.K


class TestGateControls(unittest.TestCase):
    def test_group_order(self):
        n_lin, perms = E.group()
        self.assertEqual(n_lin, 11232)          # |GL(3,3)| = 26 * 24 * 18
        self.assertEqual(len(perms), 11232 * 27)

    def test_orbit_sizes_of_5_sets(self):
        self.assertEqual(sum(o for _, o in E.reps(5)), comb(27, 5))

    def test_k3_end_to_end(self):
        # Lev 2009: the smallest Kakeya set in F_3^3 has 13 points. The same engine, run on every cut
        # in F_3^3, must find nothing at 12 points and a verified set at 13.
        self.assertIsNone(K.k3_control(12)["found"])
        f = K.k3_control(13)["found"]
        self.assertIsNotNone(f)
        self.assertTrue(f["is_kakeya"])
        self.assertLessEqual(f["size"], 13)

    def test_known27_is_recovered(self):
        # Positive control in F_3^4: from the 6- and 8-point sections of the known 27-point set
        # (put at s = 0 and s = 1 by the affine map s -> s - 1), the exact C-search must find a C of
        # at most 13 points that completes them to a Kakeya set.
        S = E.space()
        self.assertTrue(K.is_kakeya(K.KNOWN27, 4))
        self.assertEqual(len(set(K.KNOWN27)), 27)
        sec = K.sections(K.KNOWN27, 4)
        self.assertEqual([len(x) for x in sec], [13, 6, 8])
        cs = K.CutSearch(S, sec[1], 13)
        Cm = cs.solve(K.mask_of(sec[2]))
        self.assertIsNotNone(Cm)
        pts = K.assemble(S, K.mask_of(sec[1]), K.mask_of(sec[2]), Cm)
        self.assertTrue(K.is_kakeya(pts, 4))
        self.assertLessEqual(len(pts), 27)

    def test_B_enumeration_matches_brute_force(self):
        # For every 5-point representative A, kcut's pruned enumeration must count exactly the 7-sets
        # B containing 0 with B - A = F_3^3 that a plain loop over all C(26, 6) candidates counts.
        # (With |C| <= 0 the C-search never succeeds, so run() visits and counts every such B.)
        S = E.space()
        pts = [tuple((i // 3 ** k) % 3 for k in range(3)) for i in range(27)]
        self.assertEqual([tuple(p) for p in S.pts], pts)
        idx = {p: i for i, p in enumerate(pts)}
        full = (1 << 27) - 1
        total = 0
        for A, _ in E.reps(5):
            md = []
            for x in range(27):
                m = 0
                for a in A:
                    m |= 1 << idx[tuple((u - v) % 3 for u, v in zip(pts[x], pts[a]))]
                md.append(m)
            brute = 0
            for rest in itertools.combinations(range(1, 27), 6):
                m = md[0]
                for x in rest:
                    m |= md[x]
                if m == full:
                    brute += 1
            got = K.CutSearch(S, list(A), 0).run(7)["B_cover"]
            self.assertEqual(got, brute, "A=%r" % (A,))
            total += brute
        self.assertGreaterEqual(total, 0)


if __name__ == "__main__":
    unittest.main()
