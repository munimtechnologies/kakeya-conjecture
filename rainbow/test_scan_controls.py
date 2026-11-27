# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Controls for rainbow/rscan.py (rainbow/SCAN.md), plus the first case of the F_3^3 bound:
no rainbow triple of F_3^3 whose smallest set has 4 points has |A| |M| |C| < 308.

Controls: the checker on the classical triple ({0,1}, {0,1}, {0,2}) and its products; the
exact hitting-set routine against brute force on 300 random instances; the affine class
lists against orbit-stabiliser counts; tau(1) = 8 and tau(2) = 64 re-found with nothing
below them (TRIPLES.md proves both by hand); lemma H one dimension down (A on a line of the
plane gives nothing below 9 tau(1) = 72); and the (7, 8, 8) triple found by CP-SAT in
October 2026 is rainbow and is re-found by the scan."""
import math
import os
import random
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

# copied from compete-results/kakeya-s3-v2/out/s3.json (CP-SAT, 2026-10-15, caps (7, 8, 8))
A788 = [(0, 0, 0), (0, 0, 1), (0, 2, 1), (1, 0, 0), (1, 2, 1), (2, 0, 1), (2, 2, 2)]
M788 = [(0, 0, 0), (0, 1, 2), (0, 2, 0), (1, 0, 2), (1, 1, 2), (2, 0, 0), (2, 0, 2),
        (2, 1, 0)]
C788 = [(0, 0, 0), (0, 0, 2), (0, 1, 0), (0, 1, 2), (0, 2, 0), (1, 0, 2), (2, 0, 0),
        (2, 1, 2)]


def ix(pts):
    return sorted(rscan.to_idx(p) for p in pts)


def translation_least_form(m, M, C):
    """The translate (M - t, C + t) with t in M whose M has the least bit mask."""
    N, add, neg = rscan.tables(m)
    best = None
    for t in M:
        Mt = sorted(add[p][neg[t]] for p in M)
        Ct = sorted(add[p][t] for p in C)  # C - 2t = C + t
        key = sum(1 << p for p in Mt)
        if best is None or key < best[0]:
            best = (key, Mt, Ct)
    return best[1], best[2]


class TestControls(unittest.TestCase):
    def test_a_checker_and_products(self):
        G1 = ([0, 1], [0, 1], [0, 2])
        self.assertTrue(rscan.check_triple(1, *G1))
        self.assertFalse(rscan.check_triple(1, [0, 1], [0, 1], [0, 1]))
        self.assertFalse(rscan.check_triple(1, [0], [0, 1], [0, 2]))
        T, m = G1, 1
        for k in (2, 3):
            T = rscan.product_triple(m, T, 1, G1)
            m = k
            self.assertTrue(rscan.check_triple(m, *T))
            self.assertEqual([len(X) for X in T], [2 ** k] * 3)
        base = 3 / 2 ** (2 / 3)
        self.assertLess(307, base ** 9)
        self.assertGreater(308, base ** 9)
        self.assertGreater(308 ** (1 / 9), base)

    def test_b_hitting_set_matches_brute_force(self):
        rng = random.Random(20261127)
        for _ in range(300):
            n = rng.randint(3, 9)
            sets = set()
            for _ in range(rng.randint(1, 12)):
                s = 0
                for _ in range(rng.randint(1, 3)):
                    s |= 1 << rng.randrange(n)
                sets.add(s)
            want = rscan.brute_min_hit(sets, n)
            got = rscan.min_hit(sets, n)
            self.assertIsNotNone(got)
            self.assertEqual(got.bit_count(), want, (sets, n))
            self.assertTrue(all(s & got for s in sets))
            if want > 0:
                self.assertIsNone(rscan.min_hit(sets, want - 1))

    def test_c_classes_are_complete(self):
        self.assertEqual(rscan.agl_order(2), 432)
        self.assertEqual(rscan.agl_order(3), 303264)
        self.assertEqual(sorted(k for k, _, _ in rscan.classes(3, 3)), [1, 2])
        self.assertEqual([k for k, _, _ in rscan.classes(3, 4)].count(3), 1)
        for x in (4, 5, 6):
            cl = rscan.classes(3, x)
            full = [st for k, _, st in cl if k == 3]
            plan = [st for k, _, st in cl if k == 2]
            self.assertEqual(len(full) + len(plan), len(cl))
            for st in full:
                self.assertEqual(303264 % st, 0)
            for st in plan:
                self.assertEqual(432 % st, 0)
            # 39 planes; an x-set (x >= 4) spanning a plane lies in exactly one of them
            self.assertEqual(sum(303264 // st for st in full),
                             math.comb(27, x) - 39 * math.comb(9, x), x)
            self.assertEqual(sum(432 // st for st in plan), math.comb(9, x), x)
        # the plane: orbit of a class spanning the plane 432/|Stab|; of a class inside a
        # line 12 lines x 6/|Stab| (AGL(1, 3) has order 6); of a point 9
        orbit = {2: lambda st: 432 // st, 1: lambda st: 72 // st, 0: lambda st: 9 // st}
        for x in range(1, 10):
            self.assertEqual(sum(orbit[k](st) for k, _, st in rscan.classes(2, x)),
                             math.comb(9, x), x)
        self.assertEqual([st for k, _, st in rscan.classes(3, 4) if k == 3], [24])
        S = [0, 1, 5, 13, 22]
        k, form, _ = rscan.canon(3, S)
        img = sorted(rscan.to_idx(((p[0] + p[1]) % 3, (p[1] + 2) % 3, (p[2] + p[0]) % 3))
                     for p in (rscan.to_pt(i, 3) for i in S))
        self.assertEqual(rscan.canon(3, img)[1], form)
        self.assertEqual(rscan.canon(3, list(form))[1], form)

    def test_d_tau1_tau2(self):
        dl = time.time() + 120
        for m, tau in ((1, 8), (2, 64)):
            xs = tuple(range(1, 3 ** m + 1))
            res = rscan.sweep(m, tau + 1, xs, dl, workers=1)
            self.assertTrue(all(r["complete"] for r in res))
            found = [r for r in res if r["best"] is not None]
            self.assertTrue(found, m)
            self.assertEqual(min(r["best"] for r in found), tau, m)
            for r in found:
                self.assertTrue(rscan.check_triple(m, r["A"], r["M"], r["C"]))
                self.assertEqual(len(r["A"]) * len(r["M"]) * len(r["C"]), r["best"])
            res = rscan.sweep(m, tau, xs, dl, workers=1)
            self.assertTrue(all(r["complete"] for r in res))
            self.assertFalse([r for r in res if r["best"] is not None], m)

    def test_e_lemma_h_in_the_plane(self):
        # A inside a line of F_3^2 gives nothing below 9 tau(1) = 72
        dl = time.time() + 120
        n = 0
        for x in (1, 2, 3):
            for k, A, _ in rscan.classes(2, x):
                if k > 1:
                    continue
                for s in rscan.sizes_for(2, x, 72):
                    for pre in rscan.prefixes(9, s):
                        r = rscan.scan_task((2, A, s, pre, 72, dl, False))
                        self.assertTrue(r["complete"])
                        self.assertIsNone(r["best"], r)
                        n += 1
        self.assertGreater(n, 0)

    def test_f_known_788_is_refound(self):
        A, M, C = ix(A788), ix(M788), ix(C788)
        self.assertEqual((len(A), len(M), len(C)), (7, 8, 8))
        self.assertTrue(rscan.check_triple(3, A, M, C))
        self.assertEqual(rscan.affine_dim(3, A), 3)
        Ms, Cs = translation_least_form(3, M, C)
        self.assertEqual(Ms[0], 0)
        self.assertTrue(rscan.check_triple(3, A, Ms, Cs))
        dl = time.time() + 90
        r = rscan.scan_task((3, A, 8, tuple(Ms[1:]), 449, dl, False))
        self.assertTrue(r["complete"], r)
        self.assertEqual(r["covers"], 1)
        self.assertEqual(r["least"], 1)
        self.assertIsNotNone(r["best"], r)
        self.assertLessEqual(r["best"], 448)
        self.assertTrue(rscan.check_triple(3, A, r["M"], r["C"]))
        r = rscan.scan_task((3, A, 8, tuple(Ms[1:3]), 449, dl, True))
        self.assertIsNotNone(r["best"], r)
        self.assertLessEqual(r["best"], 448)
        self.assertTrue(rscan.check_triple(3, A, r["M"], r["C"]))


class TestFourPointCase(unittest.TestCase):
    def test_no_triple_below_308_with_a_4_point_smallest_set(self):
        # by lemma H the smallest set spans F_3^3; 4 spanning points form one affine class
        deadline = time.time() + 300
        res = rscan.sweep(3, 308, (4,), deadline, workers=8, full_only=True)
        by = rscan.summarise(res)
        msg = rscan.describe(by)
        hits = [r for r in res if r["best"] is not None]
        self.assertFalse(hits, "TRIPLE FOUND " + repr(
            [(r["best"], r["A"], r["M"], r["C"]) for r in hits[:3]]) + " || " + msg)
        self.assertTrue(all(r["complete"] for r in res), "INCOMPLETE || " + msg)
        self.assertEqual(sorted(by), [4], msg)
        self.assertEqual(by[4]["reps"], 1, msg)
        self.assertGreater(by[4]["nodes"], 0, msg)


if __name__ == "__main__":
    unittest.main()
