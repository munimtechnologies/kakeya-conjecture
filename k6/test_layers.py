# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""A second, independent proof that k6/sets/kakeya7_211.json is a Kakeya set in F_3^7.

Write a point of F_3^7 as (c, x) with c in F_3 and x in F_3^6, and let L_c be the layer
{x : (c, x) in K}.

* A direction (0, d) with d != 0 is covered if L_0 contains a full line in direction d,
  i.e. if L_0 is a Kakeya set of F_3^6.
* A direction (1, v), v in F_3^6 (all 729 of them, v = 0 included), is covered iff some
  line {(0, x), (1, x + v), (2, x + 2v)} lies in K, i.e. iff some x in L_0 has x + v in L_1
  and x + 2v in L_2.  (Every line in such a direction meets each layer exactly once.)
* Every direction of F_3^7 is a nonzero multiple of exactly one of these.

The test checks that L_0 is exactly the 105-point set of k6/sets/kakeya6_105.json, that this
set is Kakeya in F_3^6 (directly, by its own loop), and the 729 conditions on L_1 and L_2.
It shares no code with k6/test_sets.py or with any search."""
import itertools
import json
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    with open(os.path.join(HERE, "sets", name), encoding="utf-8") as f:
        return json.load(f)["points"]


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


class Layers211(unittest.TestCase):
    def test_layer_argument(self):
        K = [digits(i, 7) for i in load("kakeya7_211.json")]
        self.assertEqual(len(set(K)), 211)
        L = [set(p[1:] for p in K if p[0] == c) for c in range(3)]
        self.assertEqual([len(x) for x in L], [105, 56, 50])
        A = set(digits(i, 6) for i in load("kakeya6_105.json"))
        self.assertEqual(L[0], A)
        self.assertTrue(kakeya_6(L[0]))
        missing = [v for v in itertools.product(range(3), repeat=6)
                   if not any(plus(x, v) in L[1] and plus(x, v, 2) in L[2] for x in L[0])]
        self.assertEqual(missing, [])

    def test_argument_detects_failure(self):
        # removing any one point of layer 1 or 2 breaks some condition (both layers are tight)
        K = [digits(i, 7) for i in load("kakeya7_211.json")]
        L = [set(p[1:] for p in K if p[0] == c) for c in range(3)]
        V = list(itertools.product(range(3), repeat=6))
        for c in (1, 2):
            for y in sorted(L[c]):
                L1 = L[1] - {y} if c == 1 else L[1]
                L2 = L[2] - {y} if c == 2 else L[2]
                self.assertTrue(any(not any(plus(x, v) in L1 and plus(x, v, 2) in L2
                                            for x in L[0]) for v in V), (c, y))


if __name__ == "__main__":
    unittest.main()
