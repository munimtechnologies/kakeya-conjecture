# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Tests of k9/split.py: the doubling lemma (Lemma 1), the lifting lemma (Lemma 2) and a
split of the 105-point Kakeya set of F_3^6 found by the search, lifted through the layered
records 211 -> 422 -> 844 and doubled into F_3^10 and F_3^11.

The doubled sets are checked by this file's own Kakeya checkers (by directions for small n,
by pairs x + y + z = 0 for n = 10), which share no code with k9/split.py."""
import itertools
import os
import sys
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from k9 import split  # noqa: E402


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


class Codec(unittest.TestCase):
    def test_roundtrip(self):
        st = [1, 2, 3, 3, 1, 2, 2]
        code = split.encode_split(st)
        self.assertEqual(len(code), 3)
        self.assertEqual(code[0], "b")          # 16 * 1 + 4 * 2 + 3 = 27 -> 'b'
        self.assertEqual(split.decode_split(code, 7), st)


class ClassicalChain(unittest.TestCase):
    def test_lemmas_rebuild_the_classical_sizes(self):
        # K = F_3, A = {0, 1}, B = {0, 2} is a split with excess 1; doubling gives
        # 2^(n+1) - 1 points in F_3^n, the size of the classical construction
        K, A, B = [0, 1, 2], {0, 1}, {0, 2}
        ok, e, miss = split.verify_split(1, K, A, B)
        self.assertEqual((ok, e, miss), (True, 1, 0))
        m = 1
        for n in range(2, 6):
            K, A, B = split.double(m, K, A, B)
            m = n
            self.assertEqual(len(K), 2 ** (n + 1) - 1)
            self.assertTrue(kakeya_by_directions(n, K))
            self.assertEqual(split.verify_split(n, K, A, B), (True, 1, 0))

    def test_verify_detects_failure(self):
        # without the common point 0 the zero vector is missing
        self.assertEqual(split.verify_split(1, [0, 1, 2], {0, 1}, {2})[0], False)
        # A u B must be K
        self.assertEqual(split.verify_split(1, [0, 1, 2], {0}, {0, 2})[0], False)


class RecordChain(unittest.TestCase):
    def test_split_of_105_lifted_and_doubled(self):
        n6, K6 = split.load_set("kakeya6_105.json")
        self.assertEqual((n6, len(K6)), (6, 105))
        self.assertTrue(kakeya_by_directions(6, K6))
        best = None
        for seed in range(1, 6):
            r = split.find_split(6, K6, seed=seed, max_steps=100000)
            if best is None or r["excess"] < best["excess"]:
                best = r
            if best["excess"] == 1:
                break
        st = split.decode_split(best["code"], 105)
        A, B = split.split_from_states(K6, st)
        ok, e, miss = split.verify_split(6, K6, A, B)
        self.assertTrue(ok, (e, miss))
        self.assertEqual(e, best["excess"])
        self.assertLessEqual(e, 4, "excess %d: the doubled F_3^10 set would not beat 1693" % e)
        # Lemma 1 on the 105-point set itself: a Kakeya set of 210 + e points in F_3^7
        D7, _, _ = split.double(6, K6, A, B)
        self.assertEqual(len(D7), 210 + e)
        self.assertTrue(kakeya_by_directions(7, D7))
        # Lemma 2 through the layered records 211 -> 422 -> 844
        K, m = K6, 6
        for name, size in (("kakeya7_211.json", 211), ("kakeya8_422.json", 422),
                           ("kakeya9_844.json", 844)):
            n, Kn = split.load_set(name)
            self.assertEqual((n, len(Kn)), (m + 1, size))
            self.assertEqual(set(i for i in Kn if i < 3 ** m), set(K))
            A, B = split.lift(m, Kn, A, B)
            self.assertEqual(split.verify_split(n, Kn, A, B), (True, e, 0))
            K, m = Kn, n
        # Lemma 1 on the 844-point set: a Kakeya set of 1688 + e points in F_3^10
        D10, A10, B10 = split.double(9, K, A, B)
        self.assertEqual(len(D10), 1688 + e)
        self.assertLess(len(D10), 1693)
        self.assertTrue(kakeya_by_pairs(10, D10))
        # and once more: its split (Corollary) gives 2 (1688 + e) + e points in F_3^11
        self.assertEqual(split.verify_split(10, D10, A10, B10), (True, e, 0))
        D11, _, _ = split.double(10, D10, A10, B10)
        self.assertEqual(len(D11), 2 * (1688 + e) + e)


if __name__ == "__main__":
    unittest.main()
