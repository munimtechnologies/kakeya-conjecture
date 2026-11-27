# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Tests of the stored splits k9/splits/*.json and the record they give in F_3^10.

* The split codes are decoded twice (k9/split.py decode_split and this file's own loop) and
  must agree; every state must be 1, 2 or 3 and exactly one point must be in both A and B
  (excess 1).
* The split of the 105-point Kakeya set of F_3^6 and the direct split of the 211-point set
  of F_3^7 are checked from the definition C(A, B) = F_3^m by this file's own loop and by
  split.verify_split.  Lemma 1 on them gives sets of 211 (F_3^7) and 423 (F_3^8) points,
  checked here as Kakeya sets directly.
* The split of the 105-point set is lifted (Lemma 2) through the layered records
  211 -> 422 -> 844 -> 1693; each lifted split is verified from the definition, with
  excess 1.
* Lemma 1 on the 844-point set of F_3^9 gives a Kakeya set of 2 * 844 + 1 = 1689 points in
  F_3^10 (the previous record was 1693).  It is checked here by this file's own pair checker
  over all 29524 directions (three distinct points of F_3^n are collinear iff they sum to 0),
  which shares no code with k9/split.py and does not use the lemma.  Its own split (the
  Corollary of k9/SPLIT.md) is verified too, which by Lemma 1 proves the 3379-point set of
  F_3^11; k9/test_doubled_big.py checks that set directly."""
import itertools
import json
import os
import sys
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from k9 import split  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_"


def read_split(name):
    with open(os.path.join(HERE, "splits", name), encoding="utf-8") as f:
        return json.load(f)


def states_of_code(code, count):
    """This file's own decoder: character ALPHA[16 s0 + 4 s1 + s2] -> s0, s1, s2."""
    out = []
    for ch in code:
        v = ALPHA.index(ch)
        out.extend((v // 16, (v // 4) % 4, v % 4))
    return out[:count], out[count:]


def tup(i, n):
    return tuple((i // 3 ** (n - 1 - k)) % 3 for k in range(n))


def dirs(n):
    for d in itertools.product(range(3), repeat=n):
        nz = [c for c in d if c]
        if nz and nz[0] == 1:
            yield d


def kakeya_by_directions(n, idx):
    S = set(tup(i, n) for i in idx)
    for d in dirs(n):
        if not any(tuple((a + b) % 3 for a, b in zip(x, d)) in S and
                   tuple((a + 2 * b) % 3 for a, b in zip(x, d)) in S for x in S):
            return False
    return True


def kakeya_by_pairs(n, idx):
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
    return len(cov) == (3 ** n - 1) // 2


def split_covers(n, K, A, B):
    """C(A, B) = {b - a : a in A, b in B, -(a + b) in K} is all of F_3^n (own loop)."""
    T = set(tup(k, n) for k in K)
    At = [tup(a, n) for a in sorted(A)]
    Bt = [tup(b, n) for b in sorted(B)]
    cov = set()
    for a in At:
        for b in Bt:
            if tuple((-x - y) % 3 for x, y in zip(a, b)) in T:
                cov.add(tuple((y - x) % 3 for x, y in zip(a, b)))
    return len(cov) == 3 ** n


class StoredSplits(unittest.TestCase):
    def load(self, fname):
        d = read_split(fname)
        n, K = split.load_set(d["set_file"])
        self.assertEqual(n, d["n"])
        self.assertEqual(len(K), d["size"])
        self.assertEqual(K, sorted(K))
        st, pad = states_of_code(d["states_code"], len(K))
        self.assertEqual(len(st), len(K))
        self.assertTrue(len(pad) < 3 and all(s == 0 for s in pad))
        self.assertEqual(st, split.decode_split(d["states_code"], len(K)))
        self.assertTrue(all(s in (1, 2, 3) for s in st))
        self.assertEqual(st.count(3), d["excess"])
        A, B = split.split_from_states(K, st)
        self.assertEqual(A | B, set(K))
        self.assertEqual(len(A & B), d["excess"])
        return d, n, K, A, B

    def test_split_of_the_105_point_set(self):
        d, n, K, A, B = self.load("kakeya6_105_split.json")
        self.assertEqual((n, len(K), d["excess"]), (6, 105, 1))
        self.assertTrue(kakeya_by_directions(6, K))
        self.assertTrue(split_covers(6, K, A, B))
        self.assertEqual(split.verify_split(6, K, A, B), (True, 1, 0))
        # Lemma 1: a Kakeya set of 2 * 105 + 1 = 211 points in F_3^7
        D7, A7, B7 = split.double(6, K, A, B)
        self.assertEqual(len(D7), 211)
        self.assertTrue(kakeya_by_directions(7, D7))
        self.assertTrue(split_covers(7, D7, A7, B7))

    def test_split_of_the_211_point_set(self):
        d, n, K, A, B = self.load("kakeya7_211_split.json")
        self.assertEqual((n, len(K), d["excess"]), (7, 211, 1))
        self.assertTrue(split_covers(7, K, A, B))
        self.assertEqual(split.verify_split(7, K, A, B), (True, 1, 0))
        D8, _, _ = split.double(7, K, A, B)
        self.assertEqual(len(D8), 423)
        self.assertTrue(kakeya_by_pairs(8, D8))

    def test_verifier_detects_failure(self):
        d, n, K, A, B = self.load("kakeya6_105_split.json")
        # without the common point the zero vector is not covered
        both = A & B
        self.assertEqual(len(both), 1)
        self.assertFalse(split_covers(6, K, A - both, B))
        # moving every point of B into A only leaves B empty, nothing is covered
        self.assertFalse(split_covers(6, K, set(K), both))


class ChainTo10(unittest.TestCase):
    def test_lift_through_the_records_and_double(self):
        d = read_split("kakeya6_105_split.json")
        m, K = split.load_set(d["set_file"])
        st = split.decode_split(d["states_code"], len(K))
        A, B = split.split_from_states(K, st)
        self.assertEqual(split.verify_split(m, K, A, B), (True, 1, 0))
        base9 = None
        for name, size in (("kakeya7_211.json", 211), ("kakeya8_422.json", 422),
                           ("kakeya9_844.json", 844), ("kakeya10_1693.json", 1693)):
            n, Kn = split.load_set(name)
            self.assertEqual((n, len(Kn)), (m + 1, size))
            # the layer x0 = 0 of each record is exactly the previous record
            self.assertEqual(set(i for i in Kn if i < 3 ** m), set(K))
            A, B = split.lift(m, Kn, A, B)
            self.assertEqual(split.verify_split(n, Kn, A, B), (True, 1, 0))
            if n == 9:
                base9 = (Kn, set(A), set(B))
            K, m = Kn, n
        K9, A9, B9 = base9
        # Lemma 1 on the 844-point set: 2 * 844 + 1 = 1689 points in F_3^10
        D10, A10, B10 = split.double(9, K9, A9, B9)
        self.assertEqual(len(D10), 1689)
        self.assertEqual(len(set(D10)), 1689)
        self.assertTrue(all(0 <= i < 3 ** 10 for i in D10))
        self.assertLess(len(D10), 1693)
        self.assertEqual(sorted(i for i in D10 if i < 3 ** 9), sorted(K9))
        # direct check of all 29524 directions, independent of the lemmas
        self.assertTrue(kakeya_by_pairs(10, D10))
        # its split has excess 1 again, which proves 2 * 1689 + 1 = 3379 points in F_3^11
        self.assertEqual(split.verify_split(10, D10, A10, B10), (True, 1, 0))
        D11, _, _ = split.double(10, D10, A10, B10)
        self.assertEqual(len(D11), 3379)


if __name__ == "__main__":
    unittest.main()
