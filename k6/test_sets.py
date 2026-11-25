# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Stand-alone check of the record Kakeya sets in k6/sets/ (standard library only).

It shares no code with the search: it decodes each stored integer into a coordinate tuple
of F_3^n and checks, for every direction d (nonzero, first nonzero entry 1), that some
point x of the set has x + d and x + 2d in the set as well.  It also checks the numbers
quoted in k6/results/records.json (sizes, the classical construction, the Bukh-Chao bound
and the values of Lev's pattern)."""
import itertools
import json
import math
import os
import unittest
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
SETS = {6: ("kakeya6_106.json", 106), 7: ("kakeya7_213.json", 213)}


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
            data = load(name)
            self.assertEqual(data["n"], n)
            self.assertEqual(data["size"], size)
            idx = data["points"]
            self.assertEqual(len(idx), size)
            self.assertEqual(len(set(idx)), size)
            self.assertTrue(all(0 <= i < 3 ** n for i in idx))
            pts = [decode(n, i) for i in idx]
            self.assertEqual(len(set(pts)), size)
            self.assertEqual(uncovered_directions(n, pts), [], "n=%d" % n)

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
        self.assertEqual(R["lev_pattern"]["values"], {str(n): lev[n] for n in range(1, 8)})
        # the record in F_3^6 lies below the pattern's value 107
        self.assertLess(SETS[6][1], lev[6])


if __name__ == "__main__":
    unittest.main()
