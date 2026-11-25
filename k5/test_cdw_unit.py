# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Offline checks of cdw.py: the .npy reader on hand-made files, and the exact neighbourhood
search against brute force in F_3^3 (where k_3 = 13, Lev 2009)."""
import itertools
import struct
import unittest

try:
    from . import cdw, walk
except ImportError:
    import cdw
    import walk


def make_npy(descr, shape, vals, fortran=False):
    hdr = "{'descr': '%s', 'fortran_order': %s, 'shape': %r, }" % (descr, fortran, shape)
    pad = 64 - (10 + len(hdr) + 1) % 64
    hdr = hdr + " " * pad + "\n"
    code = cdw._CODES[descr[1:]]
    body = struct.pack("<" + str(len(vals)) + code, *vals)
    return b"\x93NUMPY\x01\x00" + struct.pack("<H", len(hdr)) + hdr.encode("latin1") + body


class Npy(unittest.TestCase):
    def test_points_rows(self):
        shape, vals, _ = cdw.parse_npy(make_npy("<i8", (2, 3), [0, 1, 2, 2, 2, -1]))
        self.assertEqual(shape, (2, 3))
        self.assertEqual(cdw.to_points(shape, vals, 3), [(0, 1, 2), (2, 2, 2)])

    def test_fortran(self):
        # C-order matrix [[0,1,2],[1,1,0]] stored column by column
        shape, vals, _ = cdw.parse_npy(make_npy("<i4", (2, 3), [0, 1, 1, 1, 2, 0], True))
        self.assertEqual(vals, [0, 1, 2, 1, 1, 0])

    def test_indicator(self):
        flat = [0] * 27
        flat[5] = flat[26] = 1
        shape, vals, _ = cdw.parse_npy(make_npy("|b1", (3, 3, 3), flat))
        self.assertEqual(cdw.to_points(shape, vals, 3), [(0, 1, 2), (2, 2, 2)])


class Neighbourhood(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.T = cdw.tables(3)
        r = walk.run(3, "trivial", 1, 10.0, T=cls.T)
        cls.K13 = r["points"]
        cls.extra = [p for p in range(27) if p not in set(cls.K13)][:2]

    def test_seed_is_optimal(self):
        self.assertEqual(len(self.K13), 13)
        self.assertTrue(walk.is_kakeya(3, [self.T["pts"][p] for p in self.K13]))

    def test_no_smaller_neighbour(self):
        by_r, _ = cdw.neighbourhood(3, self.K13, 3, 10 ** 12, 1)
        for r, a in by_r.items():
            self.assertEqual(a["done"], a["total"])
            self.assertEqual(a["found"], set(), "k_3 = 13 forbids a 12-point Kakeya set")

    def test_recovers_seed(self):
        K14 = sorted(self.K13 + self.extra[:1])
        NB = cdw.Nbhd(self.T, K14)
        target = frozenset(self.K13)
        self.assertEqual(sum(f == target for R in itertools.combinations(K14, 1)
                             for f in NB.search(R, 0)), 1)
        self.assertEqual(sum(f == target for R in itertools.combinations(K14, 2)
                             for f in NB.search(R, 1)), 13)

    def test_matches_brute_force(self):
        T, pts = self.T, self.T["pts"]
        K15 = sorted(self.K13 + self.extra)
        NB = cdw.Nbhd(T, K15)
        for r in (1, 2, 3):
            for R in itertools.combinations(K15, r):
                keep = [p for p in K15 if p not in R]
                found = NB.search(R, r - 1)
                for f in found:
                    self.assertTrue(walk.is_kakeya(3, [pts[p] for p in f]))
                    self.assertTrue(set(keep) <= f and len(f) <= len(keep) + r - 1)
                rest = [p for p in range(27) if p not in set(keep)]
                brute = any(walk.is_kakeya(3, [pts[p] for p in keep + list(A)])
                            for s in range(r) for A in itertools.combinations(rest, s))
                self.assertEqual(bool(found), brute, "R=%r" % (R,))


if __name__ == "__main__":
    unittest.main()
