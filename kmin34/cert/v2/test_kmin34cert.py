# Tests for kmin34cert.py. MV2 of Munim, Inc. Made with AI agents working in MV2's name.
# Licence: Apache-2.0. The solver-free tests need only the standard library; the last class
# needs the kissat and drat-trim binaries (built in the job's setup) and is skipped without them.
import math
import os
import random
import unittest

import kmin34cert as kc


def popcount(b):
    return bin(b).count("1")


class TestBasics(unittest.TestCase):
    def test_field(self):
        self.assertTrue(kc.field_ok())
        self.assertEqual(kc.MUL[2][3], 1)
        self.assertEqual(kc.MUL[3][3], 2)

    def test_directions_and_lines(self):
        self.assertEqual(len(kc.DIR2), 5)
        self.assertEqual(len(kc.DIR3), 21)
        for d in kc.DIR2:
            self.assertEqual(len(kc.LINES2[d]), 4)
            self.assertEqual(sorted(i for L in kc.LINES2[d] for i in L), list(range(16)))
        # hand example: direction (0, 1) gives the columns {4a, 4a+1, 4a+2, 4a+3}
        self.assertEqual(kc.LINES2[(0, 1)], [(0, 1, 2, 3), (4, 5, 6, 7), (8, 9, 10, 11),
                                             (12, 13, 14, 15)])

    def test_kakeya_brute_force(self):
        self.assertTrue(kc.is_kakeya3(kc.KNOWN28))
        self.assertEqual(len(set(kc.KNOWN28)), 28)
        self.assertTrue(kc.is_kakeya3(kc.PTS3))
        # every line of direction (w, 1) meets s = 0, so removing that plane kills them all
        no_bottom = [x for x in kc.PTS3 if x[2] != 0]
        self.assertFalse(kc.is_kakeya3(no_bottom))
        self.assertEqual(kc.min_section(kc.PTS3), 16)

    def test_job_tag(self):
        self.assertEqual(kc.JOB, "kmin34cert-v2")


class TestGroup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.perms = kc.agl2()

    def test_group(self):
        self.assertTrue(kc.group_ok(self.perms))

    def test_coverage(self):
        cov = kc.coverage(self.perms)
        for k in range(1, 7):
            self.assertTrue(cov[k]["ok"], cov[k])
            self.assertEqual(cov[k]["covered"], math.comb(16, k))
        self.assertEqual(sum(len(kc.REPS[k]) for k in range(1, 7)), 22)

    def test_normalise_known28(self):
        n0, sl = kc.normalise(kc.KNOWN28, self.perms)
        self.assertEqual(n0, 4)
        self.assertIn(sl[0], kc.REPS[4])
        self.assertEqual(sum(len(sl[s]) for s in range(4)), 28)
        Kn = kc.from_slices(sl)
        self.assertTrue(kc.is_kakeya3(Kn))
        self.assertEqual(kc.min_section(Kn), 4)


class TestCounters(unittest.TestCase):
    def test_atmost_atleast_solver_free(self):
        for n in range(1, 6):
            lits = list(range(1, n + 1))
            for kind in ("atmost", "atleast"):
                for k in range(-1, n + 2):
                    for bits in range(2 ** n):
                        f = kc.CNF(nfixed=n)
                        getattr(f, kind)(lits, k)
                        base = {}
                        for v in range(n):
                            on = bool(bits >> v & 1)
                            base[v + 1] = on
                            f.add([v + 1] if on else [-(v + 1)])
                        pc = popcount(bits)
                        want = pc <= k if kind == "atmost" else pc >= k
                        self.assertEqual(kc.canonical_ok(f, base), want, (kind, n, k, bits))
                        self.assertEqual(kc.up_conflict(f), not want, (kind, n, k, bits))

    def test_dimacs_header(self):
        f = kc.CNF(nfixed=3)
        f.atmost([1, 2, 3], 1)
        lines = f.dimacs().strip().split("\n")
        self.assertEqual(lines[0], f"p cnf {f.n} {len(f.clauses)}")
        self.assertEqual(len(lines) - 1, len(f.clauses))
        self.assertTrue(all(l.endswith(" 0") for l in lines[1:]))


class TestModel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.perms = kc.agl2()

    def check(self, K, bound):
        n0, sl = kc.normalise(sorted(K), self.perms)
        if sl is None:
            return None
        fix = {s: set(sl[s]) for s in (1, 2, 3)}
        f = kc.build_cnf(sl[0], bound, fix)
        base = {kc.X(s, i): (i in fix[s]) for s in (1, 2, 3) for i in range(16)}
        want = kc.is_kakeya3(K) and len(K) <= bound
        self.assertEqual(kc.canonical_ok(f, base), want)
        self.assertEqual(kc.up_conflict(f), not want)
        return want

    def test_known28(self):
        self.assertTrue(self.check(kc.KNOWN28, 28))
        self.assertFalse(self.check(kc.KNOWN28, 27))

    def test_perturbed_sets(self):
        rng = random.Random(3)
        seen = {True: 0, False: 0}
        for _ in range(24):
            K = set(kc.KNOWN28)
            for _ in range(rng.randint(0, 2)):
                K.discard(rng.choice(sorted(K)))
            for _ in range(rng.randint(0, 5)):
                K.add(rng.choice(kc.PTS3))
            r = self.check(sorted(K), rng.choice([len(K), 64]))
            if r is not None:
                seen[r] += 1
        self.assertGreater(seen[True], 0)
        self.assertGreater(seen[False], 0)

    def test_all_case_formulas_build(self):
        for n0 in range(1, 7):
            for R in kc.REPS[n0]:
                f = kc.build_cnf(R, 27)
                self.assertGreater(f.n, 48)
                self.assertTrue(all(1 <= abs(l) <= f.n for cl in f.clauses for l in cl))
                self.assertTrue(all(cl for cl in f.clauses))

    def test_single_point_slice_needs_49(self):
        # |K_0| = 1 forces K_1 = K_2 = K_3 = F_4^2 (the 16 lines a + t w cover each slice).
        full = {s: set(range(16)) for s in (1, 2, 3)}
        base = {kc.X(s, i): True for s in (1, 2, 3) for i in range(16)}
        self.assertTrue(kc.canonical_ok(kc.build_cnf([0], 49, full), base))
        self.assertTrue(kc.up_conflict(kc.build_cnf([0], 48, full)))


KISSAT = os.environ.get("KISSAT", "/opt/kissat/build/kissat")
DRAT = os.environ.get("DRAT_TRIM", "/opt/drat-trim/drat-trim")


@unittest.skipUnless(os.path.exists(KISSAT) and os.path.exists(DRAT), "kissat/drat-trim absent")
class TestTools(unittest.TestCase):
    def test_pigeonhole_certified(self):
        tmp = "/tmp/kmin34cert-test"
        os.makedirs(tmp, exist_ok=True)
        p = os.path.join(tmp, "php.cnf")
        with open(p, "w") as fh:
            fh.write(kc.php_cnf(3).dimacs())
        st = kc.run_kissat(KISSAT, p, p + ".drat", 120)[0]
        self.assertEqual(st, "UNSAT")
        self.assertEqual(kc.run_drat(DRAT, p, p + ".drat", 120)[0], "VERIFIED")

    def test_single_point_slice_unsat(self):
        tmp = "/tmp/kmin34cert-test"
        os.makedirs(tmp, exist_ok=True)
        p = os.path.join(tmp, "one.cnf")
        with open(p, "w") as fh:
            fh.write(kc.build_cnf([0], 27).dimacs())
        self.assertEqual(kc.run_kissat(KISSAT, p, p + ".drat", 120)[0], "UNSAT")
        self.assertEqual(kc.run_drat(DRAT, p, p + ".drat", 300)[0], "VERIFIED")
        st, tv, _, _ = kc.solve_cnf(KISSAT, kc.build_cnf([0], 49), p, 120)
        self.assertEqual(st, "SAT")
        K = kc.from_slices(kc.slices_from_model([0], tv))
        self.assertTrue(kc.is_kakeya3(K))


if __name__ == "__main__":
    unittest.main()
