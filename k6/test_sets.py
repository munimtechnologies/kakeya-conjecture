# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Stand-alone check of the record Kakeya sets in k6/sets/ (standard library only).

It shares no code with the search: it decodes each stored integer into a coordinate tuple
of F_3^n and checks, for every direction d (nonzero, first nonzero entry 1), that some
point x of the set has x + d and x + 2d in the set as well.  It also checks the numbers
quoted in k6/results/records.json (sizes, layer sizes, the classical construction, the
Bukh-Chao bound and the values of Lev's pattern)."""
import itertools
import json
import math
import os
import unittest
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
# current records
SETS = {6: ("kakeya6_105.json", 105), 7: ("kakeya7_213.json", 213)}
# earlier, superseded sets that stay in the repository and are still checked
EARLIER = [(6, "kakeya6_106.json", 106)]


def load(name):
    with open(os.path.join(HERE, "sets", name), encoding="utf-8") as f:
        return json.load(f)


def decode(n, i):
    return tuple((i // 3 ** (n - 1 - k)) % 3 for k in range(n))


def directions(n):
    for d in itertools.product(range(3), repeat=n):
        nz = [c for c in d if c]
        if nz and nz[0] == 1:
            yield d


def uncovered_directions(n, pts):
    S = set(pts)
    bad = []
    for d in directions(n):
        if not any(tuple((a + b) % 3 for a, b in zip(x, d)) in S and
                   tuple((a + 2 * b) % 3 for a, b in zip(x, d)) in S for x in S):
            bad.append(d)
    return bad


def classical(n):
    S = set()
    for pat in itertools.product((0, 1), repeat=n):
        for v in (1, 2):
            S.add(tuple(v * c for c in pat))
    return S


def check_file(tc, n, name, size):
    data = load(name)
    tc.assertEqual(data["n"], n)
    tc.assertEqual(data["size"], size)
    idx = data["points"]
    tc.assertEqual(len(idx), size)
    tc.assertEqual(len(set(idx)), size)
    tc.assertTrue(all(0 <= i < 3 ** n for i in idx))
    pts = [decode(n, i) for i in idx]
    tc.assertEqual(len(set(pts)), size)
    tc.assertEqual(uncovered_directions(n, pts), [], "%s n=%d" % (name, n))
    return pts


class RecordSets(unittest.TestCase):
    def test_checker_on_known_cases(self):
        # the classical set is Kakeya with 2^(n+1) - 1 points; dropping its origin is not
        for n in (2, 3, 4):
            C = classical(n)
            self.assertEqual(len(C), 2 ** (n + 1) - 1)
            self.assertEqual(uncovered_directions(n, C), [])
            self.assertNotEqual(uncovered_directions(n, C - {(0,) * n}), [])
        self.assertEqual(sum(1 for _ in directions(6)), 364)
        self.assertEqual(sum(1 for _ in directions(7)), 1093)

    def test_record_sets(self):
        for n, (name, size) in SETS.items():
            check_file(self, n, name, size)

    def test_earlier_sets(self):
        for n, name, size in EARLIER:
            check_file(self, n, name, size)

    def test_quoted_numbers(self):
        with open(os.path.join(HERE, "results", "records.json"), encoding="utf-8") as f:
            R = json.load(f)
        lev = {1: 3, 2: 7}
        for n in range(3, 8):
            lev[n] = lev[n - 1] + 2 * lev[n - 2]
        self.assertEqual([lev[n] for n in range(1, 8)], [3, 7, 13, 27, 53, 107, 213])
        for n in range(1, 8):  # the same values as Lev's 'pattern' 2k+1, 2k-1 alternately
            self.assertEqual(lev[n], (5 * 2 ** n + (-1) ** n) // 3)
        for n, (name, size) in SETS.items():
            r = R["records"][str(n)]
            self.assertEqual(r["size"], size)
            self.assertEqual(r["set_file"], "k6/sets/" + name)
            self.assertEqual(r["directions"], (3 ** n - 1) // 2)
            self.assertEqual(r["classical_construction"], 2 ** (n + 1) - 1)
            self.assertEqual(r["lev_pattern_value"], lev[n])
            self.assertEqual(r["bukh_chao_lower_bound"],
                             math.ceil(Fraction(3 ** n) / Fraction(5, 3) ** (n - 1)))
            pts = [decode(n, i) for i in load(name)["points"]]
            layers = [sum(1 for x in pts if x[0] == c) for c in range(3)]
            self.assertEqual(r["layers_x0"], layers)
        for n, name, size in EARLIER:
            e = [x for x in R["earlier"] if x["set_file"] == "k6/sets/" + name]
            self.assertEqual(len(e), 1)
            self.assertEqual(e[0]["size"], size)
            self.assertEqual(e[0]["n"], n)
        self.assertEqual(R["lev_pattern"]["values"], {str(n): lev[n] for n in range(1, 8)})
        # the record in F_3^6 lies below the pattern's value 107 and below 2 * 53 = 106
        self.assertLess(SETS[6][1], lev[6])
        self.assertEqual(SETS[6][1], 2 * lev[5] - 1)
        # the recurrence on true values would ask k_7 <= k_6 + 2 k_5 <= 105 + 2 * 53 = 211
        self.assertEqual(R["recurrence_on_true_values"]["k7_at_most"], SETS[6][1] + 2 * lev[5])


if __name__ == "__main__":
    unittest.main()
