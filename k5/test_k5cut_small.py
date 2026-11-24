# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Controls for k5cut: orbit enumeration, pruning conditions and the full pipeline,
checked against brute force in small dimensions."""
import itertools
import random
import unittest
from math import comb

try:
    from . import k5cut
except ImportError:
    import k5cut


class TestOrbits(unittest.TestCase):
    def test_orbit_sums_full(self):
        for d, maxk in ((2, 9), (3, 5), (4, 3)):
            sp = k5cut.Space(d)
            lev = sp.orbit_levels(maxk, {k: 10 ** 9 for k in range(1, maxk + 1)})
            for k in range(1, maxk + 1):
                self.assertEqual(sum(s for _, s in lev[k]), comb(sp.n, k), (d, k))

    def test_known_counts(self):
        # AG(3,3): 5 orbits of 5-sets (as in the kcutc run); AG(2,3): 2 orbits of 3-sets
        sp3 = k5cut.Space(3)
        lev = sp3.orbit_levels(5, {k: 10 ** 9 for k in range(1, 6)})
        self.assertEqual(len(lev[5]), 5)
        sp2 = k5cut.Space(2)
        lev2 = sp2.orbit_levels(3, {k: 10 ** 9 for k in range(1, 4)})
        self.assertEqual(len(lev2[3]), 2)

    def test_filtered_levels_against_brute_force(self):
        sp = k5cut.Space(3)
        for b in (0, 1, 2):
            lev = sp.orbit_levels(5, {k: b for k in range(1, 6)})
            for k in range(2, 6):
                brute = sum(1 for S in itertools.combinations(range(27), k) if sp.excess(S) <= b)
                self.assertEqual(sum(s for _, s in lev[k]), brute, (b, k))

    def test_canon_is_invariant(self):
        rng = random.Random(7)
        sp = k5cut.Space(4)
        for _ in range(20):
            S = rng.sample(range(81), rng.randint(2, 9))
            key = sp.canon(S)[0]
            # random affine map: random invertible matrix and translation
            while True:
                M = [[rng.randrange(3) for _ in range(4)] for _ in range(4)]
                if _det3(M) % 3:
                    break
            t = rng.randrange(81)
            img = []
            for p in S:
                c = sp.coords[p]
                w = [sum(M[i][j] * c[j] for j in range(4)) % 3 for i in range(4)]
                q = sum(w[i] * 3 ** i for i in range(4))
                img.append(sp.add[q][t])
            self.assertEqual(sp.canon(img)[0], key)


def _det3(M):
    # determinant of a 4x4 integer matrix (exact), reduced mod 3 by the caller
    n = len(M)
    if n == 1:
        return M[0][0]
    return sum((-1) ** j * M[0][j] * _det3([row[:j] + row[j + 1:] for row in M[1:]]) for j in range(n))


def brute_Z(sp, x, y):
    """Least |C| over all A, B of sizes x, y in F_3^d (d small) satisfying (*); None if never."""
    best = None
    for A in itertools.combinations(range(sp.n), x):
        for B in itertools.combinations(range(sp.n), y):
            r, rows = k5cut.rows_and_columns(sp, A, B)
            if min(r) == 0:
                continue
            lim = (best - 1) if best else sp.n
            while lim >= 1 and k5cut.cover(rows, sp.full, lim) is not None:
                best = lim
                lim -= 1
    return best


class TestPipelineDim2(unittest.TestCase):
    def test_pipeline_matches_brute_force(self):
        sp = k5cut.Space(2)
        checked = 0
        for x in range(1, 10):
            for y in range(x, 10):
                if x * y < 9:
                    continue
                Zb = brute_Z(sp, x, y)
                if Zb is None:
                    res = k5cut.search_profile(sp, x, y, 9, stop_at_first=True)
                    self.assertNotEqual(res["verdict"], "FEASIBLE", (x, y))
                    continue
                res = k5cut.search_profile(sp, x, y, Zb, stop_at_first=True)
                self.assertEqual(res["verdict"], "FEASIBLE", (x, y, Zb))
                if Zb > 1:
                    res = k5cut.search_profile(sp, x, y, Zb - 1)
                    self.assertIn(res["verdict"], ("CLOSED", "CLOSED-by-bounds"), (x, y, Zb))
                checked += 1
        self.assertGreater(checked, 20)


class TestNecessaryConditions(unittest.TestCase):
    def _check(self, sp, A, B, z):
        Emin, need = k5cut.plan(sp, len(A), len(B), z)
        r, rows = k5cut.rows_and_columns(sp, A, B)
        self.assertGreaterEqual(sum(t * t for t in r), Emin)
        k1 = max(r)
        self.assertIn(k1, need)
        v0 = r.index(k1)
        I = sorted(set(A) & {sp.sub[b][v0] for b in B})
        self.assertEqual(len(I), k1)
        self.assertLessEqual(sp.excess(I), need[k1])

    def test_lev_cut(self):
        sp = k5cut.Space(3)
        A, B, C = k5cut.lev_sections(sp)
        self.assertEqual((len(A), len(B), len(C)), (6, 8, 13))
        self.assertTrue(k5cut.check_star(sp, A, B, C))
        self.assertLessEqual(k5cut.min_cover_size(sp, A, B), 13)
        self._check(sp, A, B, 13)

    def test_random_instances_dim3(self):
        sp = k5cut.Space(3)
        rng = random.Random(2026)
        found = 0
        tries = 0
        while found < 25 and tries < 200000:
            tries += 1
            x, y = rng.choice(((6, 8), (5, 7), (6, 7), (7, 7)))
            A = rng.sample(range(27), x)
            B = rng.sample(range(27), y)
            r, rows = k5cut.rows_and_columns(sp, A, B)
            if min(r) == 0:
                continue
            z = k5cut.min_cover_size(sp, A, B)
            self._check(sp, A, B, z)
            found += 1
        self.assertEqual(found, 25)


if __name__ == "__main__":
    unittest.main()
