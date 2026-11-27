# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Direct checks of the doubled Kakeya sets in F_3^11 and F_3^12, and of the numbers quoted
in k9/results/doubling.json.

The sets are rebuilt from the stored split of the 105-point set (k9/splits/), lifted through
the layered records 211 -> 422 -> 844 (Lemma 2 of k9/SPLIT.md) and doubled (Lemma 1):
844 -> 1689 (F_3^10) -> 3379 (F_3^11) -> 6759 (F_3^12) -> 13519 (F_3^13).

* The 3379-point set of F_3^11 is checked over all 88573 directions by pairs (three distinct
  points are collinear iff they sum to 0), without using the lemmas.
* The 6759-point set of F_3^12 is checked by the layer criterion: its layer x0 = 0 is the
  3379-point set (just checked), and every v in F_3^11 (177147 vectors) has some x in that
  layer with x + v in the layer x0 = 1 and x + 2v = -(x + (x + v)) in the layer x0 = 2.
* The split of the 6759-point set is verified from the definition (531441 vectors), which by
  Lemma 1 (hand proof only) proves the 13519-point set of F_3^13.
The checkers here are this file's own loops; only the construction uses k9/split.py."""
import json
import math
import os
import sys
import unittest
from fractions import Fraction

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from k9 import split  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def tup(i, n):
    return tuple((i // 3 ** (n - 1 - k)) % 3 for k in range(n))


def covered_directions_by_pairs(n, idx):
    P = sorted(set(tup(i, n) for i in idx))
    S = set(P)
    cov = set()
    for i, x in enumerate(P):
        for y in P[i + 1:]:
            if tuple((-a - b) % 3 for a, b in zip(x, y)) in S:
                d = tuple((b - a) % 3 for a, b in zip(x, y))
                if [c for c in d if c][0] == 2:
                    d = tuple((2 * c) % 3 for c in d)
                cov.add(d)
    return len(cov)


def chain_to_10():
    with open(os.path.join(HERE, "splits", "kakeya6_105_split.json"), encoding="utf-8") as f:
        d = json.load(f)
    m, K = split.load_set(d["set_file"])
    A, B = split.split_from_states(K, split.decode_split(d["states_code"], len(K)))
    for name in ("kakeya7_211.json", "kakeya8_422.json", "kakeya9_844.json"):
        n, Kn = split.load_set(name)
        A, B = split.lift(m, Kn, A, B)
        K, m = Kn, n
    return split.double(9, K, A, B)


class DoubledBig(unittest.TestCase):
    def test_f3_11_f3_12_and_the_split_for_f3_13(self):
        D10, A10, B10 = chain_to_10()
        self.assertEqual(len(D10), 1689)
        D11, A11, B11 = split.double(10, D10, A10, B10)
        self.assertEqual(len(D11), 3379)
        self.assertEqual(len(set(D11)), 3379)
        self.assertTrue(all(0 <= i < 3 ** 11 for i in D11))
        # F_3^11: all 88573 directions, from pairs
        self.assertEqual(covered_directions_by_pairs(11, D11), (3 ** 11 - 1) // 2)
        # F_3^12: layer criterion over the 3379-point layer
        D12, A12, B12 = split.double(11, D11, A11, B11)
        self.assertEqual(len(D12), 6759)
        self.assertEqual(len(set(D12)), 6759)
        N = 3 ** 11
        self.assertEqual(sorted(i for i in D12 if i < N), sorted(D11))
        L0 = [tup(i, 11) for i in D12 if i < N]
        L1 = [tup(i - N, 11) for i in D12 if N <= i < 2 * N]
        L2 = set(tup(i - 2 * N, 11) for i in D12 if 2 * N <= i < 3 * N)
        self.assertEqual(len(L0) + len(L1) + len(L2), 6759)
        cov = set()
        for x in L0:
            for y in L1:
                if tuple((-a - b) % 3 for a, b in zip(x, y)) in L2:
                    cov.add(tuple((b - a) % 3 for a, b in zip(x, y)))
        self.assertEqual(len(cov), 3 ** 11)
        # the split of the 6759-point set (Corollary), for Lemma 1 into F_3^13
        T = set(tup(k, 12) for k in D12)
        At = [tup(a, 12) for a in sorted(A12)]
        Bt = [tup(b, 12) for b in sorted(B12)]
        self.assertEqual(set(A12) | set(B12), set(D12))
        self.assertEqual(len(set(A12) & set(B12)), 1)
        cv = set()
        for a in At:
            for b in Bt:
                if tuple((-x - y) % 3 for x, y in zip(a, b)) in T:
                    cv.add(tuple((y - x) % 3 for x, y in zip(a, b)))
        self.assertEqual(len(cv), 3 ** 12)
        self.assertEqual(2 * len(D12) + 1, 13519)


class QuotedNumbers(unittest.TestCase):
    def test_doubling_json(self):
        with open(os.path.join(HERE, "results", "doubling.json"), encoding="utf-8") as f:
            R = json.load(f)
        lev = {n: (5 * 2 ** n + (-1) ** n) // 3 for n in range(1, 16)}
        self.assertEqual([lev[n] for n in range(1, 6)], [3, 7, 13, 27, 53])
        for n in range(3, 16):
            self.assertEqual(lev[n], lev[n - 1] + 2 * lev[n - 2])
        ours = {6: 105, 7: 211, 8: 422, 9: 844}
        for n in range(10, 16):
            ours[n] = 845 * 2 ** (n - 9) - 1
        for n in range(10, 15):
            self.assertEqual(ours[n + 1], 2 * ours[n] + 1)
        self.assertEqual(ours[10], 2 * 844 + 1)
        for key, r in R["table"].items():
            n = int(key)
            self.assertEqual(r["ours"], ours[n], n)
            self.assertEqual(r["lev_pattern_value"], lev[n], n)
            self.assertEqual(r["below_lev_by"], lev[n] - ours[n], n)
            self.assertEqual(r["classical_construction"], 2 ** (n + 1) - 1, n)
            self.assertEqual(r["directions"], (3 ** n - 1) // 2, n)
            self.assertEqual(r["bukh_chao_lower_bound"],
                             math.ceil(Fraction(3 ** n) / Fraction(5, 3) ** (n - 1)), n)
        self.assertEqual(sorted(int(k) for k in R["table"]), list(range(6, 14)))
        prev = {10: 1693}
        for n, p in prev.items():
            self.assertEqual(R["table"][str(n)]["previous_record"], p)
        # the gap to the pattern: 25 * 2^n / 1536 + 1 + (-1)^n / 3 for n >= 9
        for n in range(9, 16):
            self.assertEqual(Fraction(lev[n] - ours[n]),
                             Fraction(25 * 2 ** n, 1536) + 1 + Fraction((-1) ** n, 3))
        F = R["family"]
        self.assertEqual(F["base_size"], 844)
        self.assertEqual(F["excess"], 1)
        self.assertEqual(F["ratio_numerator"], 845)
        self.assertEqual(F["ratio_denominator"], 512)
        self.assertLess(Fraction(845, 512), Fraction(5, 3))


if __name__ == "__main__":
    unittest.main()
