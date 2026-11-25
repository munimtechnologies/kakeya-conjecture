# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Recompute k5/results/pencil26.json: the pencil LP excludes none of the 9 quadric types at
26 points in F_3^4 (a negative result about the method; see k5/PENCIL26.md)."""
import json
import os
import time
import unittest

try:
    from . import pencil_gate
except ImportError:
    import pencil_gate

HERE = os.path.dirname(os.path.abspath(__file__))


class Pencil26Results(unittest.TestCase):
    def test_no_type_excluded(self):
        with open(os.path.join(HERE, "results", "pencil26.json"), encoding="utf-8") as fh:
            data = json.load(fh)
        t0 = time.time()
        ctx = pencil_gate.prepare(data["n"], data["N"], data["rule"])
        self.assertEqual(len(ctx["profs"]), data["profiles"])
        self.assertEqual(len(ctx["feats"]), data["features"])
        self.assertEqual(len(ctx["types"]), len(data["types"]))
        for rec in data["types"]:
            z = ctx["types"][rec["k"]]
            self.assertEqual(list(z), rec["z"])
            res = pencil_gate.run_type(ctx, z, deadline=t0 + 540)
            self.assertNotIn("timeout_at_u", res, pencil_gate.summary(ctx, res))
            self.assertFalse(res["excluded"], pencil_gate.summary(ctx, res))
            self.assertEqual(res["z_level"], rec["z_level"], pencil_gate.summary(ctx, res))
            self.assertEqual(res["feasible_u"], rec["feasible_u"], pencil_gate.summary(ctx, res))
            self.assertEqual(res["unsettled_u"], [], pencil_gate.summary(ctx, res))


if __name__ == "__main__":
    unittest.main()
