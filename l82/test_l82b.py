# Tests for l82b.py (symmetry breaking for the l82sat cases). Apache-2.0.
# MV2 of Munim, Inc.; made with AI agents working in MV2's name.

import random
import unittest

import l82
import l82b


class TestGroup(unittest.TestCase):
    def test_rank(self):
        self.assertEqual(l82b.rank([1, 2, 3], 3), 2)
        self.assertEqual(l82b.rank([1, 2, 4], 3), 3)
        self.assertEqual(l82b.affine_rank([0, 1, 2, 3], 3), 2)
        self.assertEqual(l82b.affine_rank([0, 1, 2, 4, 7], 3), 3)

    def test_affine_basis_group_is_symmetric_group(self):
        # an affine basis of AG(4,2) has 5 points; every permutation extends uniquely
        maps, complete = l82b.stabilizer_maps([0, 1, 2, 4, 8], 4)
        self.assertTrue(complete)
        self.assertEqual(len(maps), 119)
        self.assertEqual(len(set(maps)), 119)

    def test_maps_fix_U_and_are_bijective(self):
        for m, a in ((4, 5), (5, 6), (7, 8)):
            cands = l82.candidates(a, m)
            random.seed(a)
            for r, U in random.sample(cands, min(6, len(cands))):
                maps, _ = l82b.stabilizer_maps(U, m, node_cap=50000)
                for cols, c in maps[:300]:
                    img = sorted(l82b.lin_image(u, cols, c) for u in U)
                    self.assertEqual(img, sorted(U))
                    P = l82b.perm_of(cols, m)
                    self.assertEqual(sorted(P), list(range(1 << m)))
                for cols, c in l82b.transvections(r, m):
                    self.assertEqual(c, 0)
                    for u in U:
                        self.assertEqual(l82b.lin_image(u, cols), u)
                    self.assertEqual(sorted(l82b.perm_of(cols, m)), list(range(1 << m)))

    def test_model_is_invariant(self):
        # a solution mapped by a group element is again a solution (fixT model feasible)
        n, N = 5, 9
        m = n - 1
        a_lo, _, a_avg = l82.a_range(n, N)
        hits = 0
        for a in range(a_lo, a_avg + 1):
            for r, U in l82.candidates(a, m):
                model, t = l82.build(n, N, a, U)
                st, T, _ = l82.solve(model, t, 20, 2)
                if T is None:
                    continue
                maps, _ = l82b.stabilizer_maps(U, m)
                for cols, c in (maps[:20] + l82b.transvections(r, m)[:10]):
                    T2 = sorted(l82b.lin_image(z, cols) for z in T)
                    self.assertTrue(l82.reduced_ok(T2, U, m), (U, cols, c))
                    hits += 1
        self.assertGreater(hits, 0)

    def test_lex_leader_keeps_feasibility(self):
        # every (a, U) case of l(5,2) = 9 and l(6,2) = 13, at N = l - 1 and N = l:
        # feasible with symmetry constraints <=> feasible without
        compared = 0
        for n in (5, 6):
            L = l82.KNOWN[n]
            for N in (L - 1, L):
                a_lo, _, a_avg = l82.a_range(n, N)
                for a in range(a_lo, a_avg + 1):
                    for r, U in l82.candidates(a, n - 1):
                        m0, t0 = l82.build(n, N, a, U)
                        s0, _, _ = l82.solve(m0, t0, 30, 2)
                        m1, t1, info = l82b.build(n, N, a, U, sym=True, K=40)
                        s1, T1, _ = l82.solve(m1, t1, 30, 2)
                        self.assertIn(s0, ("FEASIBLE", "INFEASIBLE"), (n, N, U))
                        self.assertEqual(s0, s1, (n, N, a, U, info))
                        if T1 is not None:
                            S, ok = l82b.check_found(T1, U, n, N)
                            self.assertTrue(ok)
                        compared += 1
        self.assertGreater(compared, 10)

    def test_lex_constraint_semantics(self):
        # brute force on 4 variables: t >=_lex t o P exactly
        from ortools.sat.python import cp_model
        P = [0, 2, 3, 1, 4]  # rotate positions 1, 2, 3
        sols = set()
        for mask in range(16):
            model = cp_model.CpModel()
            t = [None] + [model.NewBoolVar(f"t{i}") for i in range(1, 5)]
            for i in range(1, 5):
                model.Add(t[i] == ((mask >> (i - 1)) & 1))
            l82b.add_lex_geq(model, t, P)
            st = cp_model.CpSolver().Solve(model)
            x = [None] + [(mask >> (i - 1)) & 1 for i in range(1, 5)]
            y = [None] + [x[P[i]] for i in range(1, 5)]
            want = x[1:] >= y[1:]
            self.assertEqual(st in (cp_model.OPTIMAL, cp_model.FEASIBLE), want, mask)

    def test_select_cases_matches_v1(self):
        self.assertEqual(len(l82b.select_cases([7])), 1)
        self.assertEqual(len(l82b.select_cases([8])), 35)
        self.assertEqual(len(l82b.select_cases([9])), 34)
        self.assertEqual(len(l82b.select_cases([10])), 499)
        p0 = l82b.select_cases([7, 8], (0, 2))
        p1 = l82b.select_cases([7, 8], (1, 2))
        self.assertEqual(len(p0) + len(p1), 36)
        for U in l82b.V1_CLOSED_A8 + l82b.V1_OPEN_A7 + l82b.V1_OPEN_A9_LOW_R:
            a = len(U)
            self.assertIn(U, [x[1] for x in l82.candidates(a, 7)])


if __name__ == "__main__":
    unittest.main()
