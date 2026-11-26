# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Stand-alone check of the record Kakeya sets in k6/sets/ (standard library only).

It shares no code with the search: it decodes each stored integer into a coordinate tuple
of F_3^n and checks, for every direction d (nonzero, first nonzero entry 1), that some
point x of the set has x + d and x + 2d in the set as well.  For n >= 10 the same property
is checked by running over pairs instead: three distinct points x, y, z of F_3^n lie on a
line iff x + y + z = 0, so the directions covered by K are exactly the directions y - x of
the pairs x != y in K with -(x + y) in K.  Both methods are compared on small cases.
It also checks the numbers quoted in k6/results/records.json (sizes, layer sizes, the
classical construction, the Bukh-Chao bound and the values of Lev's pattern).

Most files list their points directly ("points").  kakeya8_422.json, kakeya9_846.json,
kakeya9_844.json and kakeya10_1702.json are stored in two parts: their layer x0 = 0 is the
file named in "layer0_file" (itself two-part for the 846-, 844- and 1702-point sets), the
rest is offset + the free points.  points_of() rebuilds them.  kakeya8_422.json lists its
free points and also keeps the compact "points_code" printed by the search (checked to
decode to the same list); the other two-part files keep only the code, decoded here."""
import itertools
import json
import math
import os
import unittest
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
# current records
SETS = {6: ("kakeya6_105.json", 105), 7: ("kakeya7_211.json", 211),
        8: ("kakeya8_422.json", 422), 9: ("kakeya9_844.json", 844),
        10: ("kakeya10_1702.json", 1702)}
# records that are also checked to be minimal under deleting one point
MINIMAL = (6, 7)
# earlier, superseded sets that stay in the repository and are still checked
EARLIER = [(6, "kakeya6_106.json", 106), (7, "kakeya7_213.json", 213),
           (8, "kakeya8_427.json", 427), (9, "kakeya9_846.json", 846)]

T32 = "abcdefghijklmnopqrstuvwxyz012345"
C32 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ6789+/"


def decode_code(s):
    """Inverse of the search's delta code (written again here, not imported)."""
    out, prev, val = [], -1, 0
    for ch in s:
        if ch in C32:
            val = val * 32 + C32.index(ch)
            continue
        val = val * 32 + T32.index(ch)
        prev += val
        out.append(prev)
        val = 0
    return out


def load(name):
    with open(os.path.join(HERE, "sets", name), encoding="utf-8") as f:
        return json.load(f)


def points_of(data):
    if "points" in data:
        return list(data["points"])
    with open(os.path.join(ROOT, data["layer0_file"]), encoding="utf-8") as f:
        base = points_of(json.load(f))      # the layer file may be two-part itself
    free = decode_code(data["points_code"])
    if "free_points" in data:
        assert free == data["free_points"], "points_code != free_points"
    if "free_count" in data:
        assert len(free) == data["free_count"], "free_count"
    off = data["offset"]
    assert off == 3 ** (data["n"] - 1)
    assert all(0 <= b < off for b in base)
    assert all(0 <= p < 2 * off for p in free)
    return sorted(list(base) + [off + p for p in free])


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


def uncovered_by_pairs(n, pts):
    """The same answer as uncovered_directions, from the pairs of points of the set."""
    S = set(pts)
    P = sorted(S)
    cov = set()
    for i, x in enumerate(P):
        for y in P[i + 1:]:
            if tuple((-a - b) % 3 for a, b in zip(x, y)) in S:
                d = tuple((b - a) % 3 for a, b in zip(x, y))
                if [c for c in d if c][0] == 2:
                    d = tuple((2 * c) % 3 for c in d)
                cov.add(d)
    return [d for d in directions(n) if d not in cov]


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
    idx = points_of(data)
    tc.assertEqual(len(idx), size)
    tc.assertEqual(len(set(idx)), size)
    tc.assertTrue(all(0 <= i < 3 ** n for i in idx))
    pts = [decode(n, i) for i in idx]
    tc.assertEqual(len(set(pts)), size)
    if n >= 10:
        tc.assertEqual(uncovered_by_pairs(n, pts), [], "%s n=%d" % (name, n))
    else:
        tc.assertEqual(uncovered_directions(n, pts), [], "%s n=%d" % (name, n))
    return pts


class RecordSets(unittest.TestCase):
    def test_checker_on_known_cases(self):
        # the classical set is Kakeya with 2^(n+1) - 1 points; dropping its origin is not
        for n in (2, 3, 4):
            C = classical(n)
            self.assertEqual(len(C), 2 ** (n + 1) - 1)
            self.assertEqual(uncovered_directions(n, C), [])
            self.assertEqual(uncovered_by_pairs(n, C), [])
            self.assertNotEqual(uncovered_directions(n, C - {(0,) * n}), [])
            self.assertEqual(uncovered_by_pairs(n, C - {(0,) * n}),
                             uncovered_directions(n, C - {(0,) * n}))
        # the two methods agree on the 105-point record with a point removed
        pts = [decode(6, i) for i in points_of(load("kakeya6_105.json"))]
        T = set(pts[1:])
        self.assertNotEqual(uncovered_directions(6, T), [])
        self.assertEqual(uncovered_by_pairs(6, T), uncovered_directions(6, T))
        self.assertEqual(uncovered_by_pairs(6, pts), [])
        self.assertEqual(sum(1 for _ in directions(6)), 364)
        self.assertEqual(sum(1 for _ in directions(7)), 1093)
        self.assertEqual(sum(1 for _ in directions(8)), 3280)
        self.assertEqual(sum(1 for _ in directions(9)), 9841)
        self.assertEqual(sum(1 for _ in directions(10)), 29524)
        # the decoder on a hand-made example: gaps 11, 18, 6*32+24 = 216
        self.assertEqual(decode_code("lsGy"), [10, 28, 244])

    def test_record_sets(self):
        for n, (name, size) in SETS.items():
            check_file(self, n, name, size)

    def test_record_sets_are_minimal_under_deletion(self):
        # no single point of a record set can be deleted (it would be a smaller record)
        for n in MINIMAL:
            name, size = SETS[n]
            pts = [decode(n, i) for i in points_of(load(name))]
            S = set(pts)
            dirs = list(directions(n))
            for p in pts:
                T = S - {p}
                # directions whose lines in S pass through p are the only candidates
                cand = [d for d in dirs
                        if any(tuple((a - k * b) % 3 for a, b in zip(p, d)) in S and
                               tuple((a + (1 - k) * b) % 3 for a, b in zip(p, d)) in S and
                               tuple((a + (2 - k) * b) % 3 for a, b in zip(p, d)) in S
                               for k in (0, 1, 2))]
                lost = [d for d in cand
                        if not any(tuple((a + b) % 3 for a, b in zip(x, d)) in T and
                                   tuple((a + 2 * b) % 3 for a, b in zip(x, d)) in T
                                   for x in T)]
                self.assertNotEqual(lost, [], (name, p))

    def test_earlier_sets(self):
        for n, name, size in EARLIER:
            check_file(self, n, name, size)

    def test_quoted_numbers(self):
        with open(os.path.join(HERE, "results", "records.json"), encoding="utf-8") as f:
            R = json.load(f)
        lev = {1: 3, 2: 7}
        for n in range(3, 11):
            lev[n] = lev[n - 1] + 2 * lev[n - 2]
        self.assertEqual([lev[n] for n in range(1, 11)],
                         [3, 7, 13, 27, 53, 107, 213, 427, 853, 1707])
        for n in range(1, 11):  # the same values as Lev's 'pattern' 2k+1, 2k-1 alternately
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
            pts = [decode(n, i) for i in points_of(load(name))]
            layers = [sum(1 for x in pts if x[0] == c) for c in range(3)]
            if "layers_x0" in r:
                self.assertEqual(r["layers_x0"], layers)
            else:   # only the layer x0 = 0 and the total of the other two are quoted
                self.assertEqual(r["layer0_x0"], layers[0])
                self.assertEqual(r["free_layers_total"], layers[1] + layers[2])
        for n, name, size in EARLIER:
            e = [x for x in R["earlier"] if x["set_file"] == "k6/sets/" + name]
            self.assertEqual(len(e), 1)
            self.assertEqual(e[0]["size"], size)
            self.assertEqual(e[0]["n"], n)
        self.assertEqual(R["lev_pattern"]["values"], {str(n): lev[n] for n in range(1, 11)})
        # the record in F_3^6 lies below the pattern's value 107 and below 2 * 53 = 106
        self.assertLess(SETS[6][1], lev[6])
        self.assertEqual(SETS[6][1], 2 * lev[5] - 1)
        # the record in F_3^7 lies below the pattern's value 213 and equals 105 + 2 * 53
        self.assertLess(SETS[7][1], lev[7])
        self.assertEqual(R["recurrence_on_true_values"]["k7_at_most"], SETS[6][1] + 2 * lev[5])
        self.assertEqual(SETS[7][1], SETS[6][1] + 2 * lev[5])
        # the record in F_3^8 lies below the classical 511 and below the pattern's value 427;
        # it is the 211-point record plus 211 more points, one more than 211 + 2 * 105 = 421
        self.assertLess(SETS[8][1], 2 ** 9 - 1)
        self.assertLess(SETS[8][1], lev[8])
        self.assertEqual(SETS[8][1], 2 * SETS[7][1])
        self.assertEqual(R["recurrence_on_true_values"]["k8_if_recurrence_bound"],
                         SETS[7][1] + 2 * SETS[6][1])
        self.assertEqual(SETS[8][1], R["recurrence_on_true_values"]["k8_if_recurrence_bound"] + 1)
        # the earlier 427-point set had 216 points outside its 211-point layer
        self.assertEqual(427 - SETS[7][1], 2 * SETS[6][1] + 6)
        # the record in F_3^9 lies below the classical 1023 and the pattern's value 853;
        # it is the 422-point record plus 422 more points, exactly 2 * 422, which is also
        # 422 + 2 * 211, the value the recurrence on our own records gives
        self.assertLess(SETS[9][1], 2 ** 10 - 1)
        self.assertLess(SETS[9][1], lev[9])
        self.assertEqual(lev[9] - SETS[9][1], 9)
        self.assertEqual(SETS[9][1], 2 * SETS[8][1])
        self.assertEqual(R["recurrence_on_true_values"]["k9_if_recurrence_bound"],
                         SETS[8][1] + 2 * SETS[7][1])
        self.assertEqual(SETS[9][1], R["recurrence_on_true_values"]["k9_if_recurrence_bound"])
        # the earlier 846-point set was 2 * 422 + 2
        self.assertEqual(846, 2 * SETS[8][1] + 2)
        # the record in F_3^10 lies below the classical 2047 and the pattern's value 1707;
        # it is the earlier 846-point set of F_3^9 plus 856 more points (846 + 10)
        self.assertLess(SETS[10][1], 2 ** 11 - 1)
        self.assertLess(SETS[10][1], lev[10])
        self.assertEqual(lev[10] - SETS[10][1], 5)
        self.assertEqual(SETS[10][1], 2 * 846 + 10)
        self.assertEqual(R["records"]["10"]["layers_x0"][0], 846)
        self.assertEqual(R["recurrence_on_true_values"]["k10_if_recurrence_bound"],
                         SETS[9][1] + 2 * SETS[8][1])


if __name__ == "__main__":
    unittest.main()
