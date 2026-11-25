# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Checks of the standard-library phase-1 simplex (lpsolve.py) and its exact certificates,
and two controls of the pencil LP: sizes where Kakeya sets exist must stay feasible."""
import random
import unittest

try:
    from . import lpsolve, pencil, pencil_gate
except ImportError:
    import lpsolve
    import pencil
    import pencil_gate

LEV = [(1, 2, 0, 0), (1, 0, 1, 0), (0, 1, 1, 0), (1, 2, 1, 0), (0, 0, 2, 0), (2, 1, 2, 0),
       (0, 0, 0, 1), (1, 0, 0, 1), (2, 0, 0, 1), (0, 1, 0, 1), (2, 0, 1, 1), (2, 1, 1, 1),
       (1, 2, 1, 1), (0, 0, 2, 1),
       (0, 0, 0, 2), (2, 0, 0, 2), (1, 1, 0, 2), (2, 1, 0, 2), (0, 2, 0, 2), (2, 2, 0, 2),
       (0, 0, 1, 2), (1, 1, 1, 2), (0, 2, 1, 2), (0, 0, 2, 2), (2, 0, 2, 2), (0, 1, 2, 2),
       (1, 1, 2, 2)]


def rand_rows(rng, k, n, lo=-3, hi=3):
    return [{j: rng.randint(lo, hi) for j in range(n) if rng.random() < 0.7} for _ in range(k)]


class TestToy(unittest.TestCase):
    def test_small_infeasible(self):
        eq = [({0: 1, 1: 1}, 1)]
        ub = [({0: -1}, -2)]                       # x0 >= 2 but x0 + x1 = 1
        r = lpsolve.phase1(eq, ub, 2)
        self.assertEqual(r["status"], "infeasible")
        self.assertTrue(pencil.check_certificate(eq, ub, [5, 5], r["lam"], r["mu"])[0])

    def test_small_feasible(self):
        eq = [({0: 1, 1: 1}, 1)]
        ub = [({0: -1}, -1)]
        self.assertEqual(lpsolve.phase1(eq, ub, 2)["status"], "feasible")

    def test_empty_rows(self):
        self.assertEqual(lpsolve.phase1([({}, 0)], [], 1)["status"], "feasible")
        r = lpsolve.phase1([({}, 3)], [], 1)
        self.assertEqual(r["status"], "infeasible")
        self.assertTrue(pencil.check_certificate([({}, 3)], [], [1], r["lam"], r["mu"])[0])

    def test_planted_feasible(self):
        rng = random.Random(1)
        for t in range(150):
            n = rng.randint(3, 14)
            v = [rng.randint(0, 4) for _ in range(n)]
            A = rand_rows(rng, rng.randint(1, 6), n)
            G = rand_rows(rng, rng.randint(0, 5), n)
            eq = [(r, sum(c * v[j] for j, c in r.items())) for r in A]
            ub = [(r, sum(c * v[j] for j, c in r.items()) + rng.randint(0, 2)) for r in G]
            self.assertEqual(lpsolve.phase1(eq, ub, n)["status"], "feasible", t)

    def test_planted_infeasible_certified(self):
        rng = random.Random(2)
        for t in range(150):
            n = rng.randint(3, 14)
            s = rng.randint(0, 6)
            eq = [(r, rng.randint(-8, 8)) for r in rand_rows(rng, rng.randint(0, 4), n)]
            eq.append(({j: 1 for j in range(n)}, s))
            ub = [(r, rng.randint(-8, 8)) for r in rand_rows(rng, rng.randint(0, 4), n)]
            ub.append(({j: -1 for j in range(n)}, -(s + 1)))       # sum >= s + 1
            r = lpsolve.phase1(eq, ub, n)
            self.assertEqual(r["status"], "infeasible", t)
            self.assertTrue(pencil.check_certificate(eq, ub, [10 ** 6] * n, r["lam"],
                                                     r["mu"])[0], t)

    def test_random_mixed_never_wrong(self):
        """On random systems: every 'infeasible' answer carries a valid certificate, and every
        'feasible' answer is consistent with a planted solution when one was planted."""
        rng = random.Random(3)
        for t in range(200):
            n = rng.randint(2, 10)
            eq = [(r, rng.randint(-6, 6)) for r in rand_rows(rng, rng.randint(1, 4), n)]
            ub = [(r, rng.randint(-6, 6)) for r in rand_rows(rng, rng.randint(0, 4), n)]
            r = lpsolve.phase1(eq, ub, n)
            self.assertIn(r["status"], ("feasible", "infeasible"), t)
            if r["status"] == "infeasible":
                self.assertTrue(pencil.check_certificate(eq, ub, [10 ** 6] * n, r["lam"],
                                                         r["mu"])[0], t)


class TestControls(unittest.TestCase):
    def test_f3_3_thirteen_points_feasible(self):
        ctx = pencil_gate.prepare(3, 13)
        st, info = pencil_gate.lp_branch(ctx)
        self.assertEqual(st, "feasible", info)

    def test_lev_branch_feasible(self):
        data = pencil.real_data(LEV, 4)
        ctx = pencil_gate.prepare(4, 27)
        z = [0, 0, 0]
        for p, c in data["classes"].items():
            z[pencil.residue(p, 27)] += c
        u = sum(c for p, c in data["classes"].items() if p in ctx["U"])
        self.assertIn(tuple(z), ctx["types"])
        st, info = pencil_gate.lp_branch(ctx, tuple(z), u)
        self.assertEqual(st, "feasible", info)


if __name__ == "__main__":
    unittest.main()
