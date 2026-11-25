# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Run the pencil LP of pencil.py (k5/PENCIL.md) with the standard-library solver lpsolve.py.

Every branch that the solver finds infeasible is backed by integer Farkas multipliers that
pencil.check_certificate verifies exactly (with the column bounds of build_lp).  A branch
counts as excluded only if that exact check passes.

Branching (as in pencil_lp.py): for each quadric type z of QUAD45.md, first the LP with the
three residue rows; if it is not excluded, every u = 0..H (number of classes with an
unbalanced profile, with the pair-incidence row sum_f C(k_U(f), 2) y_f = C(u, 2)).
N points are excluded for type z if every branch of z is excluded.
"""
import time

try:
    from . import lpsolve, pencil
except ImportError:
    import lpsolve
    import pencil


def prepare(n, N, rule="sd2"):
    t0 = time.time()
    profs = pencil.profile_list(n, N, rule)
    feats = pencil.prune(pencil.enumerate_features(n, N, rule))
    return {"n": n, "N": N, "rule": rule, "profs": profs, "feats": feats,
            "U": pencil.unbalanced_set(profs, N, 1),
            "types": sorted(set(pencil.quadric_types(n))),
            "H": pencil.consts(n)["H"], "enum_seconds": round(time.time() - t0, 1)}


def lp_branch(ctx, z=None, u=None):
    """Solve one branch.  Returns (status, info) with status 'excluded' (certified),
    'feasible', 'uncertified' (solver says infeasible, exact check failed) or 'unknown'."""
    U = ctx["U"] if u is not None else None
    eq, ub, bounds = pencil.build_lp(ctx["n"], ctx["N"], ctx["profs"], ctx["feats"], z, U, u)
    ncols = len(bounds)
    Lc = ncols - 1
    res = lpsolve.phase1(eq, ub, ncols, extra_ub=[({Lc: 1}, bounds[Lc])])
    st = res["status"]
    if st == "infeasible":
        ok, S = pencil.check_certificate(eq, ub, bounds, res["lam"], res["mu"])
        return ("excluded" if ok else "uncertified"), {"iters": res["iters"],
                                                        "obj": res.get("obj")}
    return st, {"iters": res.get("iters"), "obj": res.get("obj"), "why": res.get("why")}


def run_type(ctx, z, deadline=None):
    """All branches of quadric type z.  Returns a result dict; 'excluded' is True only if
    every branch was excluded with an exactly checked certificate."""
    t0 = time.time()
    out = {"z": z, "lps": 0}
    st, info = lp_branch(ctx, z)
    out["lps"] += 1
    out["z_level"] = st
    out["z_iters"] = info.get("iters")
    if st == "excluded":
        out["excluded"] = True
        out["seconds"] = round(time.time() - t0, 1)
        return out
    feas, other, done = [], [], []
    for u in range(ctx["H"] + 1):
        if deadline is not None and time.time() > deadline:
            out["timeout_at_u"] = u
            break
        s, inf = lp_branch(ctx, z, u)
        out["lps"] += 1
        done.append(u)
        if s == "feasible":
            feas.append(u)
        elif s != "excluded":
            other.append((u, s, inf.get("why")))
    out["feasible_u"] = feas
    out["unsettled_u"] = other
    out["u_done"] = len(done)
    out["excluded"] = (not feas and not other and len(done) == ctx["H"] + 1)
    out["seconds"] = round(time.time() - t0, 1)
    return out


def summary(ctx, res):
    return ("n=%d N=%d rule=%s: %d profiles, %d features (enumeration %.1fs); result %r"
            % (ctx["n"], ctx["N"], ctx["rule"], len(ctx["profs"]), len(ctx["feats"]),
               ctx["enum_seconds"], res))


def check_type(test, n, N, k, budget=500):
    """unittest helper: assert that quadric type number k (in sorted order) is excluded."""
    t0 = time.time()
    ctx = prepare(n, N)
    if k >= len(ctx["types"]):
        test.skipTest("only %d quadric types" % len(ctx["types"]))
    res = run_type(ctx, ctx["types"][k], deadline=t0 + budget)
    test.assertTrue(res["excluded"], summary(ctx, res))
    return ctx, res
