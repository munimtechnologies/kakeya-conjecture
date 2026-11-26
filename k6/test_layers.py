# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""A second, independent proof that the layered record sets are Kakeya sets:
k6/sets/kakeya7_211.json in F_3^7, k6/sets/kakeya8_422.json and the earlier
k6/sets/kakeya8_427.json in F_3^8, k6/sets/kakeya9_844.json and the earlier
k6/sets/kakeya9_846.json in F_3^9, and k6/sets/kakeya10_1693.json and the earlier
k6/sets/kakeya10_1695.json and k6/sets/kakeya10_1702.json in F_3^10.

Write a point of F_3^n as (c, x) with c in F_3 and x in F_3^(n-1), and let L_c be the layer
{x : (c, x) in K}.

* A direction (0, d) with d != 0 is covered if L_0 contains a full line in direction d,
  i.e. if L_0 is a Kakeya set of F_3^(n-1).
* A direction (1, v), v in F_3^(n-1) (all 3^(n-1) of them, v = 0 included), is covered iff
  some line {(0, x), (1, x + v), (2, x + 2v)} lies in K, i.e. iff some x in L_0 has x + v in
  L_1 and x + 2v in L_2.  (Every line in such a direction meets each layer exactly once.)
* Every direction of F_3^n is a nonzero multiple of exactly one of these.

For the 211-point set, L_0 is the 105-point set of k6/sets/kakeya6_105.json, checked to be
Kakeya in F_3^6 here by its own loop.  For the 422- and 427-point sets, L_0 is the 211-point
set of F_3^7, which this file proves Kakeya by the same argument one dimension down; for
the 846- and 844-point sets, L_0 is the 422-point set of F_3^8; for the 1702-point set it is
the 846-point set of F_3^9 and for the 1695- and 1693-point sets the 844-point set.  So the
chains 105 -> 211 -> 422 -> 844 -> 1693 (and -> 1695) and 105 -> 211 -> 422 -> 846 -> 1702
are proved here without k6/test_sets.py.

For the sets in F_3^10 the vectors v are collected from pairs instead of tested one by one:
with y = x + v in L_1, the third point x + 2v = 2y - x = -(x + y), so the covered vectors are
exactly {y - x : x in L_0, y in L_1, -(x + y) in L_2}; the test checks that these are all
3^9 = 19683 vectors.  Both methods are compared on the 422-point set.  The file shares no
code with k6/test_sets.py or with any search (the compact point code is decoded by its own
loop)."""
import itertools
import json
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
_HI = "ABCDEFGHIJKLMNOPQRSTUVWXYZ6789+/"   # leading base-32 digits of a gap
_LO = "abcdefghijklmnopqrstuvwxyz012345"   # last base-32 digit of a gap


def unpack(code):
    """Sorted integers from the gap code (first gap counted from -1)."""
    vals, cur, acc = [], -1, 0
    for ch in code:
        k = _LO.find(ch)
        if k < 0:
            acc = 32 * acc + _HI.index(ch)
        else:
            cur += 32 * acc + k
            vals.append(cur)
            acc = 0
    return vals


def read(path):
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    if "points" in data:
        return list(data["points"])
    # two-part file: the layer x0 = 0 is another file (possibly two-part itself),
    # the rest is offset + the free points
    base = read(os.path.join(os.path.dirname(HERE), data["layer0_file"]))
    free = data["free_points"] if "free_points" in data else unpack(data["points_code"])
    return list(base) + [data["offset"] + p for p in free]


def load(name):
    return read(os.path.join(HERE, "sets", name))


def digits(i, n):
    out = []
    for _ in range(n):
        i, r = divmod(i, 3)
        out.append(r)
    return tuple(reversed(out))


def plus(x, v, k=1):
    return tuple((a + k * b) % 3 for a, b in zip(x, v))


def kakeya_6(S):
    for d in itertools.product(range(3), repeat=6):
        if not any(d) or [c for c in d if c][0] != 1:
            continue
        if not any(plus(x, d) in S and plus(x, d, 2) in S for x in S):
            return False
    return True


def layers(name, n):
    K = [digits(i, n) for i in load(name)]
    return K, [set(p[1:] for p in K if p[0] == c) for c in range(3)]


def missing_vectors(L):
    m = len(next(iter(L[0])))
    return [v for v in itertools.product(range(3), repeat=m)
            if not any(plus(x, v) in L[1] and plus(x, v, 2) in L[2] for x in L[0])]


def covered_vectors(L):
    """{y - x : x in L_0, y in L_1, -(x + y) in L_2}, the vectors v with a line."""
    out = set()
    for x in L[0]:
        for y in L[1]:
            if tuple((-a - b) % 3 for a, b in zip(x, y)) in L[2]:
                out.add(tuple((b - a) % 3 for a, b in zip(x, y)))
    return out


class Layers211(unittest.TestCase):
    def test_layer_argument(self):
        K, L = layers("kakeya7_211.json", 7)
        self.assertEqual(len(set(K)), 211)
        self.assertEqual([len(x) for x in L], [105, 56, 50])
        A = set(digits(i, 6) for i in load("kakeya6_105.json"))
        self.assertEqual(L[0], A)
        self.assertTrue(kakeya_6(L[0]))
        self.assertEqual(missing_vectors(L), [])

    def test_argument_detects_failure(self):
        # removing any one point of layer 1 or 2 breaks some condition (both layers are tight)
        K, L = layers("kakeya7_211.json", 7)
        V = list(itertools.product(range(3), repeat=6))
        for c in (1, 2):
            for y in sorted(L[c]):
                L1 = L[1] - {y} if c == 1 else L[1]
                L2 = L[2] - {y} if c == 2 else L[2]
                self.assertTrue(any(not any(plus(x, v) in L1 and plus(x, v, 2) in L2
                                            for x in L[0]) for v in V), (c, y))


class Layers422(unittest.TestCase):
    def test_layer_argument(self):
        K, L = layers("kakeya8_422.json", 8)
        self.assertEqual(len(K), 422)
        self.assertEqual(len(set(K)), 422)
        self.assertEqual([len(x) for x in L], [211, 96, 115])
        # layer 0 is exactly the 211-point set of F_3^7 (proved Kakeya in Layers211)
        A = set(digits(i, 7) for i in load("kakeya7_211.json"))
        self.assertEqual(L[0], A)
        # all 2187 vectors v of F_3^7
        self.assertEqual(missing_vectors(L), [])

    def test_argument_detects_failure(self):
        # with a single point left in layer 2 at most 211 of the 2187 vectors can be covered
        # (one per x in layer 0), so the check must report missing vectors
        K, L = layers("kakeya8_422.json", 8)
        one = next(iter(sorted(L[2])))
        self.assertGreaterEqual(len(missing_vectors([L[0], L[1], {one}])), 2187 - 211)

    def test_pair_method_agrees(self):
        # the pair method used for the sets in F_3^10 gives the same answer here, also
        # after one point of layer 1 is removed
        K, L = layers("kakeya8_422.json", 8)
        allv = set(itertools.product(range(3), repeat=7))
        self.assertEqual(covered_vectors(L), allv)
        y = sorted(L[1])[0]
        Lm = [L[0], L[1] - {y}, L[2]]
        self.assertEqual(sorted(allv - covered_vectors(Lm)), missing_vectors(Lm))
        self.assertNotEqual(missing_vectors(Lm), [])


class Layers427(unittest.TestCase):
    def test_layer_argument(self):
        K, L = layers("kakeya8_427.json", 8)
        self.assertEqual(len(K), 427)
        self.assertEqual(len(set(K)), 427)
        self.assertEqual([len(x) for x in L], [211, 114, 102])
        A = set(digits(i, 7) for i in load("kakeya7_211.json"))
        self.assertEqual(L[0], A)
        self.assertEqual(missing_vectors(L), [])


class Layers846(unittest.TestCase):
    def test_layer_argument(self):
        K, L = layers("kakeya9_846.json", 9)
        self.assertEqual(len(K), 846)
        self.assertEqual(len(set(K)), 846)
        self.assertEqual([len(x) for x in L], [422, 204, 220])
        # |L1| |L2| >= 3^8 is necessary (each pair (y, z) serves one v at most)
        self.assertGreaterEqual(len(L[1]) * len(L[2]), 6561)
        # layer 0 is exactly the 422-point set of F_3^8 (proved Kakeya in Layers422)
        A = set(digits(i, 8) for i in load("kakeya8_422.json"))
        self.assertEqual(L[0], A)
        # all 6561 vectors v of F_3^8
        self.assertEqual(missing_vectors(L), [])

    def test_decoder(self):
        # gaps 11, 18, 6*32+24 = 216, written by hand
        self.assertEqual(unpack("lsGy"), [10, 28, 244])


class Layers844(unittest.TestCase):
    def test_layer_argument(self):
        K, L = layers("kakeya9_844.json", 9)
        self.assertEqual(len(K), 844)
        self.assertEqual(len(set(K)), 844)
        self.assertEqual([len(x) for x in L], [422, 207, 215])
        self.assertGreaterEqual(len(L[1]) * len(L[2]), 6561)
        # layer 0 is exactly the 422-point set of F_3^8 (proved Kakeya in Layers422)
        A = set(digits(i, 8) for i in load("kakeya8_422.json"))
        self.assertEqual(L[0], A)
        # all 6561 vectors v of F_3^8
        self.assertEqual(missing_vectors(L), [])


def check_ten(tc, name, size, sizes, layer0_name):
    K, L = layers(name, 10)
    tc.assertEqual(len(K), size)
    tc.assertEqual(len(set(K)), size)
    tc.assertEqual([len(x) for x in L], sizes)
    tc.assertGreaterEqual(len(L[1]) * len(L[2]), 19683)
    # layer 0 is exactly the named set of F_3^9 (proved Kakeya above)
    A = set(digits(i, 9) for i in load(layer0_name))
    tc.assertEqual(L[0], A)
    # all 19683 vectors v of F_3^9 have a line
    cov = covered_vectors(L)
    tc.assertEqual(len(cov), 3 ** 9)
    tc.assertEqual(cov, set(itertools.product(range(3), repeat=9)))


class Layers1702(unittest.TestCase):
    def test_layer_argument(self):
        check_ten(self, "kakeya10_1702.json", 1702, [846, 443, 413], "kakeya9_846.json")


class Layers1695(unittest.TestCase):
    def test_layer_argument(self):
        check_ten(self, "kakeya10_1695.json", 1695, [844, 398, 453], "kakeya9_844.json")


class Layers1693(unittest.TestCase):
    def test_layer_argument(self):
        check_ten(self, "kakeya10_1693.json", 1693, [844, 430, 419], "kakeya9_844.json")


if __name__ == "__main__":
    unittest.main()
