# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Checks of k5/results/local53.json: internal consistency, and (when the published file can
be downloaded) the 53-point set's invariants and the complete r <= 2 neighbourhood again."""
import json
import os
import time
import unittest

try:
    from . import cdw, walk
except ImportError:  # pragma: no cover
    import cdw
    import walk

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    with open(os.path.join(HERE, "results", "local53.json"), encoding="utf-8") as f:
        return json.load(f)


class Consistency(unittest.TestCase):
    def test_counts(self):
        R = load()
        s = R["set_check"]
        self.assertEqual(sum(s["hyperplane_profiles"].values()), 121)
        for k in s["hyperplane_profiles"]:
            self.assertEqual(sum(map(int, k.split(","))), 53)
        self.assertEqual(sum(s["point_degree_histogram"].values()), 53)
        self.assertEqual(s["directions_with_1_line"] + s["directions_with_2_lines"], 121)
        self.assertEqual(s["directions_with_1_line"] + 2 * s["directions_with_2_lines"],
                         s["lines_inside"])
        self.assertEqual(sum(int(d) * c for d, c in s["point_degree_histogram"].items()),
                         3 * s["lines_inside"])
        from math import comb
        for row in R["exact_neighbourhood"]:
            self.assertEqual(row["deletions_total"], comb(53, row["r"]))
            self.assertEqual(row["kakeya_sets_found"], 0)


class Recheck(unittest.TestCase):
    def test_download_and_r2(self):
        try:
            P53, note = cdw.load_set(5, 53)
        except Exception as e:  # noqa: BLE001
            self.skipTest("download failed: %s" % e)
        R = load()["set_check"]
        self.assertEqual(len(P53), 53)
        self.assertTrue(walk.is_kakeya(5, P53))
        T = cdw.tables(5)
        K = cdw.indices(T, P53)
        inv = cdw.invariants(T, K)
        self.assertEqual(inv["lines"], R["lines_inside"])
        self.assertEqual(inv["removable"], 0)
        self.assertEqual(inv["centres"], 0)
        self.assertEqual({",".join(map(str, k)): v for k, v in inv["profiles"].items()},
                         R["hyperplane_profiles"])
        by_r, _ = cdw.neighbourhood(5, K, 2, time.time() + 10 ** 6, 1)
        for r, a in by_r.items():
            self.assertEqual(a["done"], a["total"])
            self.assertEqual(a["found"], set())


if __name__ == "__main__":
    unittest.main()
