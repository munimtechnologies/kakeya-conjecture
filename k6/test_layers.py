# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""A second, independent proof that the layered record sets are Kakeya sets:
k6/sets/kakeya7_211.json in F_3^7, and k6/sets/kakeya8_422.json and the earlier
k6/sets/kakeya8_427.json in F_3^8.

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
set of F_3^7, which this file proves Kakeya by the same argument one dimension down.  So the
chain 105 -> 211 -> 422 is proved here without k6/test_sets.py.  It shares no code with
k6/test_sets.py or with any search."""
import itertools
import json
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    with open(os.path.join(HERE, "sets", name), encoding="utf-8") as f:
        data = json.load(f)
    if "points" in data:
        return data["points"]
    # two-part file: the layer x0 = 0 is another file, the rest is offset + free_points
    with open(os.path.join(os.path.dirname(HERE), data["layer0_file"]), encoding="utf-8") as f:
        base = json.load(f)["points"]
    return list(base) + [data["offset"] + p for p in data["free_points"]]


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


class Layers427(unittest.TestCase):
    def test_layer_argument(self):
        K, L = layers("kakeya8_427.json", 8)
        self.assertEqual(len(K), 427)
        self.assertEqual(len(set(K)), 427)
        self.assertEqual([len(x) for x in L], [211, 114, 102])
        A = set(digits(i, 7) for i in load("kakeya7_211.json"))
        self.assertEqual(L[0], A)
        self.assertEqual(missing_vectors(L), [])


if __name__ == "__main__":
    unittest.main()
