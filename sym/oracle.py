# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Standard-library pass of ksym: for every field group G(r, f) with at most `max_orbits`
orbits in F_3^n, n = 3..8, the least size of a G-invariant Kakeya set by the exact branch and
bound of bb.py (no solver).  Every set found is re-checked with check.is_kakeya.
Returns one line per instance: n r f |G| orbits dir-orbits -> min (complete or not)."""
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gf  # noqa: E402
import check  # noqa: E402
import bb  # noqa: E402
import ksym  # noqa: E402


def run(dims=(3, 4, 5, 6, 7, 8), max_orbits=60, per=15.0, total=380.0):
    t0 = time.time()
    lines = []
    best = {}
    for n in dims:
        try:
            todo = [t for t in ksym.plan_n(n) if t[0] <= max_orbits]
        except Exception:  # noqa: BLE001
            lines.append("n=%d plan ERROR %s" % (n, traceback.format_exc()[-300:]))
            continue
        for norb, r, f in todo:
            if time.time() - t0 > total:
                lines.append("n=%d r=%d f=%d SKIPPED (time)" % (n, r, f))
                continue
            try:
                B = ksym.build(n, r, f)
                res = bb.minimum(B["w"], B["masks"], seconds=per)
                if res["found"] is None:
                    lines.append("n=%d r=%d f=%d orbits=%d: none found (%s nodes)" % (
                        n, r, f, norb, res["nodes"]))
                    continue
                chosen = gf.mask_bits(res["found"])
                pts = sorted(p for o in chosen for p in B["members"][o])
                ok = len(pts) == res["value"] and check.is_kakeya(n, check.to_tuples(n, pts))
                tag = "min" if res["complete"] else "best"
                lines.append("n=%d r=%d f=%d |G|=%d orbits=%d dirorbits=%d: %s %d verified=%s %.0fs" % (
                    n, r, f, r * f, norb, len(B["reps"]), tag, res["value"], ok, res["seconds"]))
                if ok and (n not in best or res["value"] < best[n][0]):
                    best[n] = (res["value"], r, f, tag)
                    if res["value"] < ksym.RECORD[n]:
                        lines.append("NEW n=%d size=%d POINTS %s" % (n, res["value"], ",".join(map(str, pts))))
            except Exception:  # noqa: BLE001
                lines.append("n=%d r=%d f=%d ERROR %s" % (n, r, f, traceback.format_exc()[-300:]))
    lines.append("BEST %s" % best)
    lines.append("total %.0fs" % (time.time() - t0))
    return " || ".join(lines)
