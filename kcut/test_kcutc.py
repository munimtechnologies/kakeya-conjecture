# test_kcutc.py - unit tests for the canonical-B reduction in kcutc.py (standard library only).
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name. Apache-2.0.
import itertools
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kcut as K  # noqa: E402
import kcutc as KC  # noqa: E402

S2 = K.Space(2)
LIN2, PERMS2 = K.affine_perms(S2)


def covering_Bs(S, A, b):
    """Every b-set B containing 0 with B - A = F_3^m (brute force)."""
    out = []
    for rest in itertools.combinations(range(1, S.P), b - 1):
        Bm = 1 | K.mask_of(rest)
        if all(any((Bm >> S.add[a][d]) & 1 for a in A) for d in range(S.P)):
            out.append(Bm)
    return out


def image(S, M, w, Bm):
    """Bitmask of M B + w."""
    return K.mask_of(S.add[w][M[x]] for x in K.bits(Bm))


class TestStabiliser(unittest.TestCase):
    def test_orbit_stabiliser_and_group(self):
        ident = list(range(S2.P))
        for k in range(1, 6):
            for A, orb in K.orbit_reps(S2, PERMS2, k):
                n_aff, Ms = KC.stabiliser(S2, LIN2, PERMS2, K.mask_of(A))
                self.assertEqual(n_aff * orb, len(PERMS2), A)
                self.assertIn(ident, Ms)
                keys = set(tuple(M) for M in Ms)
                for M1 in Ms:
                    for M2 in Ms:
                        self.assertIn(tuple(M1[M2[x]] for x in range(S2.P)), keys)
                # every M is the linear part of some affine map fixing A
                Am = K.mask_of(A)
                for M in Ms:
                    self.assertTrue(any(image(S2, M, t, Am) == Am for t in range(S2.P)))

    def test_canonical_is_orbit_minimum(self):
        rng = random.Random(5)
        for _ in range(25):
            A = sorted(rng.sample(range(9), rng.randint(2, 4)))
            b = rng.randint(3, 5)
            _, Ms = KC.stabiliser(S2, LIN2, PERMS2, K.mask_of(A))
            cov = covering_Bs(S2, A, b)
            canon = set(B for B in cov if KC.is_canonical(S2, B, Ms))
            minima = set()
            for B in cov:
                imgs = [image(S2, M, w, B) for M in Ms for w in range(S2.P)]
                minima.add(min(x for x in imgs if x & 1))
            self.assertEqual(canon, minima, (A, b))


class TestReduction(unittest.TestCase):
    def test_symmetry_claim(self):
        """solve(B) and solve(M B + w) agree for M in Lin(Stab A), any w (checked with the solver,
        which test_kcut.py checks against brute force)."""
        rng = random.Random(6)
        yes = no = 0
        for _ in range(150):
            A = sorted(rng.sample(range(9), rng.randint(2, 4)))
            _, Ms = KC.stabiliser(S2, LIN2, PERMS2, K.mask_of(A))
            Bm = K.mask_of(rng.sample(range(9), rng.randint(2, 6)))
            mC = rng.randint(2, 7)
            cs = K.CutSearch(S2, A, mC)
            r0 = cs.solve(Bm) is not None
            yes += r0
            no += not r0
            M = rng.choice(Ms)
            w = rng.randrange(9)
            self.assertEqual(cs.solve(image(S2, M, w, Bm)) is not None, r0, (A, Bm, M, w, mC))
        self.assertGreater(yes, 2)
        self.assertGreater(no, 2)

    def test_canon_matches_full(self):
        rng = random.Random(7)
        yes = no = 0
        for _ in range(60):
            A = sorted(rng.sample(range(9), rng.randint(2, 4)))
            b = rng.randint(2, 5)
            mC = rng.randint(3, 7)
            _, Ms = KC.stabiliser(S2, LIN2, PERMS2, K.mask_of(A))
            cs = K.CutSearch(S2, A, mC)
            full = cs.run(b, stop_at_first=False)
            red = KC.run_canon(cs, b, Ms, stop_at_first=False)
            self.assertEqual(red["B_cover"], full["B_cover"])
            self.assertLessEqual(red["B_canon"], red["B_cover"])
            self.assertEqual(red["found"] is not None, full["n_found"] > 0, (A, b, mC))
            if red["found"] is not None:
                Bm, Cm = red["found"]
                self.assertTrue(K.is_kakeya(K.assemble(S2, K.mask_of(A), Bm, Cm), 3))
                yes += 1
            else:
                no += 1
        self.assertGreater(yes, 2)
        self.assertGreater(no, 2)

    def test_k3_is_13_with_reduction(self):
        self.assertIsNone(KC.k3_control_canon(12)["found"])
        f = KC.k3_control_canon(13)["found"]
        self.assertIsNotNone(f)
        self.assertTrue(f["is_kakeya"])
        self.assertLessEqual(f["size"], 13)


class TestCuts(unittest.TestCase):
    def test_delta26_is_a_product(self):
        """s3 - 8 s2 + 576 = (n1 - 8)(n2 - 8)(n3 - 8) / 2 for every cut of 26 points (BAD25.md)."""
        from math import comb
        for n1 in range(27):
            for n2 in range(27 - n1):
                n3 = 26 - n1 - n2
                s2 = sum(comb(x, 2) for x in (n1, n2, n3))
                s3 = sum(comb(x, 3) for x in (n1, n2, n3))
                self.assertEqual(2 * (s3 - 8 * s2 + 576), (n1 - 8) * (n2 - 8) * (n3 - 8))


if __name__ == "__main__":
    unittest.main()
