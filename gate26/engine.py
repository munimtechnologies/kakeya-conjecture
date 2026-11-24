# gate26/engine.py - shared helpers for the commit-time checks in gate26/.
# MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name;
# no human mathematician has checked it yet. Licence: Apache-2.0. Python standard library only.
"""
Loads the solver-free search kcut/kcut.py by file path (kcut/ is not a Python package, so the
commit-time test runner, `python -m unittest discover`, never ran the tests inside kcut/) and caches
the group data that every check needs.
"""
import functools
import importlib.util
import os
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
KCUT_PATH = os.path.join(os.path.dirname(_HERE), "kcut", "kcut.py")


def _load():
    spec = importlib.util.spec_from_file_location("kcut_for_gate26", KCUT_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


K = _load()


@functools.lru_cache(maxsize=None)
def space():
    """F_3^3 as used by kcut (points 0..26)."""
    return K.Space(3)


@functools.lru_cache(maxsize=None)
def group():
    """(|GL(3,3)|, list of all 303264 elements of AGL(3,3) as point permutations)."""
    lin, perms = K.affine_perms(space())
    return len(lin), perms


@functools.lru_cache(maxsize=None)
def reps(k):
    """One representative of every AGL(3,3)-orbit of k-subsets of F_3^3, with its orbit size.
    kcut.orbit_reps raises if the orbit sizes do not add up to C(27, k)."""
    return tuple((tuple(A), orb) for A, orb in K.orbit_reps(space(), group()[1], k))


def close_cut(cut):
    """Run kcut's exact search for one cut (a, b, c) over every orbit representative A of a-sets:
    every b-set B containing 0 with B - A = F_3^3, and an exact search for C with |C| <= c.
    Returns one row per representative; row["found"] is None unless a Kakeya set was found."""
    a, b, c = cut
    S = space()
    rows = []
    for A, orb in reps(a):
        t0 = time.time()
        r = K.CutSearch(S, list(A), c).run(b)
        found = None
        if r["found"] is not None:
            Bm, Cm = r["found"]
            pts = K.assemble(S, K.mask_of(A), Bm, Cm)
            found = {"K": pts, "size": len(pts), "is_kakeya": K.is_kakeya(pts, 4)}
        rows.append({"A": list(A), "orbit": orb, "B_cover": r["B_cover"], "found": found,
                     "seconds": round(time.time() - t0, 1)})
    return rows
