# Tests for kdrat26.py (solver-free: the CNF encoding is evaluated by unit propagation and
# compared with brute force). MV2 of Munim, Inc.; made with AI agents. Licence: Apache-2.0.

import itertools
import unittest

import kdrat26 as kd
import kprof25 as kp
import kprof26 as k26


def assign_units(cnf, values):
    """values: dict var -> bool; returns the clause list with the unit clauses appended."""
    return cnf.clauses + [[v if b else -v] for v, b in values.items()]


class Totalizer(unittest.TestCase):
    def test_exact_by_unit_propagation(self):
        for m in range(1, 8):
            cnf = kd.CNF(m)
            U = kd.totalizer(cnf, list(range(1, m + 1)))
            self.assertEqual(len(U), m)
            for bits in itertools.product((False, True), repeat=m):
                cl = assign_units(cnf, {i + 1: bits[i] for i in range(m)})
                st, val = kd.unit_propagate(cl, cnf.n)
                self.assertEqual(st, "SAT", (m, bits))
                cnt = sum(bits)
                for k in range(1, m + 1):
                    lit = U[k - 1]
                    got = (val[abs(lit)] > 0) == (lit > 0)
                    self.assertEqual(got, cnt >= k, (m, bits, k))

    def test_neq_and_equal(self):
        m = 5
        cnf = kd.CNF(m)
        U = kd.totalizer(cnf, list(range(1, m + 1)))
        self.assertIsNone(kd.neq_lits(U, -1))
        self.assertIsNone(kd.neq_lits(U, m + 1))
        for bits in itertools.product((False, True), repeat=m):
            st, val = kd.unit_propagate(assign_units(cnf, {i + 1: bits[i] for i in range(m)}), cnf.n)
            cnt = sum(bits)
            for t in range(m + 1):
                lits = kd.neq_lits(U, t)
                sat = any((val[abs(l)] > 0) == (l > 0) for l in lits)
                self.assertEqual(sat, cnt != t, (bits, t))
                c2 = kd.CNF(cnf.n)
                c2.clauses = list(cnf.clauses)
                kd.equal_units(c2, U, t)
                st2, _ = kd.unit_propagate(assign_units(c2, {i + 1: bits[i] for i in range(m)}), c2.n)
                self.assertEqual(st2 == "SAT", cnt == t, (bits, t))

    def test_php_shape(self):
        c = kd.php_cnf(4, 3)
        self.assertEqual(c.n, 12)
        self.assertEqual(len(c.clauses), 4 + 3 * 6)
        self.assertIn("p cnf 12 22", c.dimacs())


class Encoding(unittest.TestCase):
    K = [tuple(p) for p in kp.KNOWN27]

    def up(self, spec):
        cnf = kd.build_cnf(spec)
        if cnf.trivial_unsat:
            return "CONFLICT"
        st, _ = kd.unit_propagate(cnf.clauses, cnf.n)
        self.assertIn(st, ("SAT", "CONFLICT"))
        return st

    def test_known27_fixed(self):
        self.assertTrue(kp.is_kakeya(self.K))
        self.assertEqual(self.up(kd.fixed_spec(self.K, 27, False)), "SAT")
        cut = tuple(sorted(kp.class_profile(self.K, kp.LAST)))
        self.assertEqual(self.up(kd.fixed_spec(self.K, 27, False, forbid=[cut])), "CONFLICT")

    def test_removals_agree_with_brute_force(self):
        for r in range(0, 27, 3):
            Kr = [p for i, p in enumerate(self.K) if i != r]
            for zero_only in (False, True):
                want = "SAT" if kd.brute_ok(Kr, 26, zero_only) else "CONFLICT"
                self.assertEqual(self.up(kd.fixed_spec(Kr, 26, zero_only)), want, (r, zero_only))

    def test_case_specs_build(self):
        for name, cut, forbid, zero_only in k26.CASES26:
            R = list(range(cut[0]))  # any slice 0 of the right size, just to build the CNF
            spec = k26.spec_for(cut, forbid, zero_only, R)
            cnf = kd.build_cnf(spec)
            self.assertLess(cnf.n, 300000)
            self.assertLess(len(cnf.clauses), 1000000)
            if zero_only:
                self.assertEqual(spec["cap"], 1)

    def test_solution_points_roundtrip(self):
        sl = kd.slices_by_last(self.K)
        spec = kd.fixed_spec(self.K, 27, False)
        model = {}
        for s in (1, 2):
            for i in range(27):
                model[kp.X(s, i) + 1] = i in sl[s]
        pts, sol = kd.solution_points(spec, model)
        self.assertEqual(sorted(pts), sorted(self.K))
        ok, info = kp.check_solution(spec, sol)
        self.assertTrue(ok, info)


if __name__ == "__main__":
    unittest.main()
