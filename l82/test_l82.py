# Tests for l82.py. MV2 of Munim, Inc. Made with AI agents working in MV2's name.
# Licence: Apache-2.0.
import itertools
import random
import unittest

import l82


def brute_min(n):
    pts = list(range(1, 1 << n))
    for k in range(1, len(pts) + 1):
        for S in itertools.combinations(pts, k):
            if l82.is_saturating(S, n):
                return k
    return None


def gl_cols(m):
    """All invertible m x m binary matrices, as tuples of column images."""
    out = []
    for cols in itertools.product(range(1, 1 << m), repeat=m):
        span = {0}
        ok = True
        for c in cols:
            if c in span:
                ok = False
                break
            span |= {s ^ c for s in span}
        if ok:
            out.append(cols)
    return out


def apply(cols, x):
    y = 0
    for i, c in enumerate(cols):
        if (x >> i) & 1:
            y ^= c
    return y


class TestBasics(unittest.TestCase):
    def test_brute_force_small(self):
        self.assertEqual(brute_min(2), 2)
        self.assertEqual(brute_min(3), 4)
        self.assertEqual(brute_min(4), 5)
        self.assertTrue(l82.is_saturating([1, 2, 4, 8, 15], 4))
        self.assertFalse(l82.is_saturating([1, 2, 4, 8], 4))
        self.assertFalse(l82.is_saturating([0, 1, 2, 4, 8, 15], 4))

    def test_a_range(self):
        self.assertEqual(l82.a_range(8, 25), (7, 19, 12))
        self.assertEqual(l82.a_range(4, 4), (None, None, None))
        self.assertEqual(l82.a_range(4, 5), (2, 4, 2))
        self.assertEqual(l82.a_range(7, 18), (5, 14, 9))

    def test_split_join(self):
        rng = random.Random(1)
        for _ in range(50):
            S = sorted(rng.sample(range(1, 64), rng.randint(1, 20)))
            T, U = l82.split(S, 6)
            self.assertEqual(l82.join(T, U, 6), S)

    def test_reduction_equivalence(self):
        rng = random.Random(2)
        seen = {True: 0, False: 0}
        for _ in range(300):
            n = rng.choice((4, 5, 6))
            S = rng.sample(range(1, 1 << n), rng.randint(3, 2 * n + 2))
            T, U = l82.split(S, n)
            a = l82.is_saturating(S, n)
            self.assertEqual(a, l82.reduced_ok(T, U, n - 1))
            seen[a] += 1
            # weights of (h, e) from the split agree with the direct dual weights
            W = l82.dual_weights(S, n)
            m = n - 1
            for h in range(1 << m):
                for e in (0, 1):
                    if h == 0 and e == 0:
                        continue
                    w = (sum(l82.parity(h & t) for t in T)
                         + sum(1 for u in U if l82.parity(h & u) ^ e))
                    self.assertEqual(w, W[h | (e << m)])
        self.assertGreater(seen[True], 10)
        self.assertGreater(seen[False], 10)

    def test_move_functional(self):
        rng = random.Random(3)
        for _ in range(100):
            n = rng.choice((4, 5, 6))
            S = rng.sample(range(1, 1 << n), rng.randint(3, 15))
            g = rng.randrange(1, 1 << n)
            S2 = l82.move_functional_last(S, g, n)
            self.assertEqual(len(set(S2)), len(S))
            self.assertNotIn(0, S2)
            self.assertEqual(l82.dual_weights(S2, n)[1 << (n - 1)], l82.dual_weights(S, n)[g])
            self.assertEqual(l82.is_saturating(S2, n), l82.is_saturating(S, n))
            self.assertEqual(sorted(l82.dual_weights(S2, n).values()),
                             sorted(l82.dual_weights(S, n).values()))

    def test_candidates_complete_m3(self):
        G = gl_cols(3)
        self.assertEqual(len(G), 168)
        AG = [(cols, c) for cols in G for c in range(8)]

        def key(X):
            return min(tuple(sorted(apply(cols, x) ^ c for x in X)) for cols, c in AG)

        for a in range(1, 9):
            want = {key(X) for X in itertools.combinations(range(8), a)}
            cands = l82.candidates(a, 3)
            for r, U in cands:
                self.assertEqual(len(set(U)), a)
                self.assertTrue(all(0 <= u < 8 for u in U))
            got = {key(U) for r, U in cands}
            self.assertEqual(got, want, a)

    def test_candidates_counts_m7(self):
        # raw normal forms before the coordinate-permutation reduction bound the counts
        for a in (7, 8, 9, 10):
            c = l82.candidates(a, 7)
            self.assertGreater(len(c), 0)
            for r, U in c:
                self.assertEqual(len(set(U)), a)
                self.assertTrue(all(0 <= u < 128 for u in U))
                base = {0} | {1 << i for i in range(r)}
                self.assertTrue(base <= set(U))
                for y in set(U) - base:
                    self.assertGreaterEqual(l82.popcount(y), 2)
                    self.assertLess(y, 1 << r)


class TestModel(unittest.TestCase):
    def test_model_matches_reduced(self):
        rng = random.Random(4)
        seen = {True: 0, False: 0}
        n = 5
        for _ in range(40):
            N = rng.randint(8, 13)
            S = rng.sample(range(1, 32), N)
            W = l82.dual_weights(S, n)
            g = min(W, key=lambda k: (W[k], k))
            S2 = l82.move_functional_last(S, g, n)
            T, U = l82.split(S2, n)
            model, t = l82.build(n, N, len(U), U, fixT=T)
            st, sol, _ = l82.solve(model, t, 30, 2)
            want = l82.is_saturating(S2, n)
            self.assertEqual(st == "FEASIBLE", want, (S2, st))
            self.assertIn(st, ("FEASIBLE", "INFEASIBLE"))
            seen[want] += 1
        self.assertGreater(seen[True], 0)
        self.assertGreater(seen[False], 0)

    def test_pipeline_reproduces_known_values(self):
        for n in (4, 5, 6):
            L = l82.KNOWN[n]
            for N, expect in ((L - 1, False), (L, True)):
                a_lo, wmax, a_avg = l82.a_range(n, N)
                avals = list(range(a_lo, a_avg + 1)) if a_lo is not None else []
                out = l82.search(n, N, avals, 60, 2, None, stop_on_found=True)
                self.assertEqual(out["bugs"], [])
                self.assertEqual(bool(out["found"]), expect, (n, N))
                if expect:
                    S = out["found"][0]["S"]
                    self.assertEqual(len(S), N)
                    self.assertTrue(l82.is_saturating(S, n))
                else:
                    self.assertTrue(all(c["status"] == "INFEASIBLE" for c in out["cases"]))


if __name__ == "__main__":
    unittest.main()
