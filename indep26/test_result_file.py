# indep26/test_result_file.py - consistency checks for indep26/k4_result.json (the k_4 = 27 summary).
# MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name;
# no human mathematician has checked it yet. Licence: Apache-2.0. Python standard library only.
import json
import os
import unittest
from math import comb

from indep26 import kk

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    with open(os.path.join(HERE, "k4_result.json"), encoding="utf-8") as f:
        return json.load(f)


class TestResultFile(unittest.TestCase):
    def test_cut_list_matches_the_search(self):
        R = load()
        cuts = [tuple(c["cut"]) for c in R["cuts"]]
        self.assertEqual(sorted(cuts), sorted(kk.CUTS26))
        self.assertEqual(sorted(cuts), sorted(kk.cuts26_from_identity()))
        for c in R["cuts"]:
            a, b, cc = c["cut"]
            self.assertEqual(a + b + cc, 26)
            self.assertTrue(a <= b <= cc and a * b >= 27)
            positive = (a - 8) * (b - 8) * (cc - 8) > 0
            self.assertEqual(c["case"], "positive" if positive else "zero")
            # every cut is searched by indep26, and has a second argument (gate26 search or lemma)
            self.assertTrue(c["indep26_search"])
            self.assertTrue((c["gate26_search"] is None) != (c["other_argument"] is None))

    def test_counts(self):
        R = load()
        n = R["counts"]
        cuts = R["cuts"]
        self.assertEqual(n["cuts_needed"], len(cuts))
        self.assertEqual(n["cuts_searched_by_indep26"], sum(1 for c in cuts if c["indep26_search"]))
        self.assertEqual(n["cuts_searched_by_gate26"], sum(1 for c in cuts if c["gate26_search"]))
        self.assertEqual(n["cuts_searched_twice"],
                         sum(1 for c in cuts if c["indep26_search"] and c["gate26_search"]))
        self.assertEqual(n["cuts_one_search_plus_lemma"], sum(1 for c in cuts if c["other_argument"]))
        self.assertEqual(n["indep26_commit_checks"], len({s for c in cuts for s in c["indep26_search"]}))
        self.assertEqual(n["gate26_commit_checks"], len({c["gate26_search"] for c in cuts if c["gate26_search"]}))

    def test_identity_numbers(self):
        I = load()["identities"]
        N = 26
        # sum over 40 classes of (n1-8)(n2-8)(n3-8) = 9 C(N,3) + 18 L - 8 * 27 C(N,2) + 40 (64 N - 512)
        const = 9 * comb(N, 3) - 8 * 27 * comb(N, 2) + 40 * (64 * N - 512)
        self.assertEqual(const, I["N26_constant"])
        self.assertEqual(const, -18 * 40)
        self.assertEqual(40 * N * N - 2 * 27 * comb(N, 2), I["zero_case_sum_of_squares"])
        self.assertEqual(max(8 * 8 + 8 * 8 + 10 * 10, 8 * 8 + 9 * 9 + 9 * 9), I["zero_case_max_if_all_sections_at_least_8"])
        self.assertEqual(40 * I["zero_case_max_if_all_sections_at_least_8"], I["zero_case_bound_40_times_228"])
        self.assertLess(I["zero_case_bound_40_times_228"], I["zero_case_sum_of_squares"])

    def test_known_27_point_set(self):
        R = load()
        self.assertEqual(R["conclusion"]["k4_lower"], 27)
        self.assertEqual(R["conclusion"]["k4_upper"], len(set(kk.KNOWN27)))
        self.assertTrue(kk.is_kakeya(kk.KNOWN27, 4))
        hist = [h["value"] for h in R["lower_bound_history"]]
        self.assertEqual(hist, sorted(hist))
        self.assertEqual(hist[-1], 27)


if __name__ == "__main__":
    unittest.main()
