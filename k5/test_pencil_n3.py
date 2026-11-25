# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Calibration: the pencil LP with exact certificates excludes 12 points in F_3^3
(re-deriving Lev's k_3 >= 13 without a search).  Passes only if every branch is certified."""
import time
import unittest

try:
    from . import pencil_gate
except ImportError:
    import pencil_gate


class PencilN3(unittest.TestCase):
    def test_twelve_points_excluded(self):
        t0 = time.time()
        ctx = pencil_gate.prepare(3, 12)
        results = [pencil_gate.run_type(ctx, z, deadline=t0 + 500) for z in ctx["types"]]
        self.assertTrue(all(r["excluded"] for r in results),
                        pencil_gate.summary(ctx, results))


if __name__ == "__main__":
    unittest.main()
