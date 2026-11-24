# test_kcut.py - unit tests for kcut.py (standard library only).
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name. Apache-2.0.
import itertools
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kcut as K  # noqa: E402


def brute_ok(S, Am, Bm, Cm):
    """Conditions (1) and (2) straight from the definitions (no W/H tables)."""
    for d in range(S.P):
        ok = False
        for a in range(S.P):
            if (Am >> a) & 1:
                b = S.add[a][d]
                c = S.add[b][d]
                if (Bm >> b) & 1 and (Cm >> c) & 1:
                    ok = True
                    break
        if not ok:
            return False
    for k in range(S.ndirs):
        if not any((L & X) == L for L in S.lines[k] for X in (Am, Bm, Cm)):
            return False
    return True


def rand_mask(rng, P, k):
    return K.mask_of(rng.sample(range(P), k))


class TestSpace(unittest.TestCase):
    def test_counts(self):
        for m, nd in ((2, 4), (3, 13)):
            S = K.Space(m)
            self.assertEqual(S.ndirs, nd)
            all_lines = set()
            for k in range(nd):
                self.assertEqual(len(S.lines[k]), S.P // 3)
                for L in S.lines[k]:
                    self.assertEqual(L.bit_count(), 3)
                    self.assertEqual(K.line_dirs(S, L), 1 << k)
                    all_lines.add(L)
            self.assertEqual(len(all_lines), nd * S.P // 3)

    def test_tables(self):
        S = K.Space(3)
        rng = random.Random(1)
        for _ in range(200):
            i, j = rng.randrange(27), rng.randrange(27)
            pi, pj = S.pts[i], S.pts[j]
            self.assertEqual(S.pts[S.sub[i][j]], tuple((x - y) % 3 for x, y in zip(pi, pj)))
            self.assertEqual(S.pts[S.negsum[i][j]], tuple((-x - y) % 3 for x, y in zip(pi, pj)))

    def test_group_orders(self):
        S2 = K.Space(2)
        lin, perms = K.affine_perms(S2)
        self.assertEqual(len(lin), 48)
        self.assertEqual(len(perms), 432)
        for k in range(1, 6):
            K.orbit_reps(S2, perms, k)  # raises if orbit sizes do not add up
        S3 = K.Space(3)
        lin, perms = K.affine_perms(S3)
        self.assertEqual(len(lin), 11232)
        reps = K.orbit_reps(S3, perms, 3)
        self.assertEqual(sorted(o for _, o in reps), [117, 2808])


class TestSearch(unittest.TestCase):
    def test_brute_matches_definition(self):
        S = K.Space(2)
        rng = random.Random(2)
        for _ in range(300):
            Am, Bm, Cm = (rand_mask(rng, 9, rng.randint(1, 7)) for _ in range(3))
            self.assertEqual(brute_ok(S, Am, Bm, Cm), K.is_kakeya(K.assemble(S, Am, Bm, Cm), 3))

    def test_solve_matches_brute(self):
        S = K.Space(2)
        rng = random.Random(3)
        allC = list(range(1 << 9))
        n_yes = 0
        for _ in range(400):
            A = sorted(rng.sample(range(9), rng.randint(1, 6)))
            Bm = rand_mask(rng, 9, rng.randint(1, 7))
            mC = rng.randint(1, 9)
            r = K.CutSearch(S, A, mC).solve(Bm)
            Am = K.mask_of(A)
            exists = any(c.bit_count() <= mC and brute_ok(S, Am, Bm, c) for c in allC)
            self.assertEqual(r is not None, exists, (A, Bm, mC))
            if r is not None:
                n_yes += 1
                self.assertLessEqual(r.bit_count(), mC)
                self.assertTrue(brute_ok(S, Am, Bm, r))
        self.assertGreater(n_yes, 20)

    def test_run_matches_brute(self):
        S = K.Space(2)
        rng = random.Random(4)
        allC = [c for c in range(1 << 9)]
        for _ in range(12):
            A = sorted(rng.sample(range(9), rng.randint(2, 4)))
            Am = K.mask_of(A)
            b = rng.randint(2, 5)
            mC = rng.randint(3, 7)
            r = K.CutSearch(S, A, mC).run(b, stop_at_first=False)
            cover = 0
            exists = False
            for rest in itertools.combinations(range(1, 9), b - 1):
                Bm = 1 | K.mask_of(rest)
                covers = all(any((Am >> a) & 1 and (Bm >> S.add[a][d]) & 1 for a in range(9))
                             for d in range(9))
                cover += covers
                if covers and not exists and any(c.bit_count() <= mC and brute_ok(S, Am, Bm, c) for c in allC):
                    exists = True
            self.assertEqual(r["B_cover"], cover)
            self.assertEqual(r["n_found"] > 0, exists)

    def test_k3_is_13(self):
        self.assertIsNone(K.k3_control(12)["found"])
        f = K.k3_control(13)["found"]
        self.assertIsNotNone(f)
        self.assertTrue(f["is_kakeya"])
        self.assertLessEqual(f["size"], 13)


class TestKnown27(unittest.TestCase):
    def test_known27(self):
        self.assertEqual(len(set(K.KNOWN27)), 27)
        self.assertTrue(K.is_kakeya(K.KNOWN27, 4))
        for p in K.KNOWN27:
            self.assertFalse(K.is_kakeya([q for q in K.KNOWN27 if q != p], 4))

    def test_known27_cut(self):
        S = K.Space(3)
        sec = K.sections(K.KNOWN27, 4)
        self.assertEqual([len(x) for x in sec], [13, 6, 8])
        m0, m1, m2 = (K.mask_of(x) for x in sec)
        self.assertTrue(brute_ok(S, m0, m1, m2))
        self.assertTrue(brute_ok(S, m1, m2, m0))  # positions shifted by s -> s - 1
        r = K.CutSearch(S, sec[1], 13).solve(m2)
        self.assertIsNotNone(r)
        self.assertTrue(K.is_kakeya(K.assemble(S, m1, m2, r), 4))
        r = K.CutSearch(S, sec[0], 8).solve(m1)
        self.assertIsNotNone(r)
        self.assertTrue(K.is_kakeya(K.assemble(S, m0, m1, r), 4))
        # shear: B translated to contain 0, C translated the opposite way, still Kakeya
        b0 = sec[2][0]
        B2 = K.mask_of(S.sub[x][b0] for x in sec[2])
        C2 = K.mask_of(S.add[x][b0] for x in sec[0])
        self.assertTrue(B2 & 1)
        self.assertTrue(brute_ok(S, m1, B2, C2))
        self.assertIsNotNone(K.CutSearch(S, sec[1], 13).solve(B2))


if __name__ == "__main__":
    unittest.main()
