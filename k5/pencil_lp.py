#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Solve the pencil LP of pencil.py with HiGHS (scipy) and extract exact Farkas certificates.

For each (n, N): enumerate pencil features; solve the LP without residue / pair constraints;
if feasible, for each quadric type z; if still feasible, for each u = number of classes with an
unbalanced profile (pair-incidence constraint).  N is excluded if every branch is infeasible.
Every infeasible branch gets a Farkas certificate, checked exactly by pencil.check_certificate.
Controls: N = 13 (n=3), N = 27 (n=4) and N = 53 (n=5) must stay feasible (sets exist), and the
(z, u) branch of a real Kakeya set must be feasible.
"""
import argparse
import json
import time

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, vstack, hstack, identity

try:
    from . import pencil
except ImportError:
    import pencil


def to_csr(rows, ncols):
    data, ri, ci = [], [], []
    for k, (r, _) in enumerate(rows):
        for j, v in r.items():
            data.append(float(v))
            ri.append(k)
            ci.append(j)
    return csr_matrix((data, (ri, ci)), shape=(len(rows), ncols))


def solve(eq, ub, bounds):
    """0 = feasible, 2 = infeasible, other = solver trouble."""
    n = len(bounds)
    A = to_csr(eq, n)
    G = to_csr(ub, n)
    res = linprog(np.zeros(n), A_ub=G, b_ub=np.array([float(h) for _, h in ub]),
                  A_eq=A, b_eq=np.array([float(b) for _, b in eq]),
                  bounds=[(0, float(u)) for u in bounds], method="highs")
    return res.status


def farkas(eq, ub, bounds):
    """Find lam (free), mu >= 0 minimising b.lam + h.mu + sum_j s_j U_j with r_j + s_j >= 0,
    |lam|, mu <= 1000; return an exactly checked certificate or None."""
    n = len(bounds)
    ne, ni = len(eq), len(ub)
    A = to_csr(eq, n)
    G = to_csr(ub, n)
    MT = hstack([A.T, G.T, identity(n)]).tocsr()          # r_j + s_j
    c = np.concatenate([[float(b) for _, b in eq], [float(h) for _, h in ub],
                        [float(u) for u in bounds]])
    bnds = [(-1000, 1000)] * ne + [(0, 1000)] * ni + [(0, None)] * n
    res = linprog(c, A_ub=-MT, b_ub=np.zeros(n), bounds=bnds, method="highs")
    if res.status != 0 or res.fun > -1e-7:
        return None, (res.status, None if res.status != 0 else float(res.fun))
    w = res.x
    lam = {i: repr(float(w[i])) for i in range(ne) if abs(w[i]) > 1e-13}
    mu = {i: repr(max(0.0, float(w[ne + i]))) for i in range(ni) if w[ne + i] > 1e-13}
    ok, S = pencil.check_certificate(eq, ub, bounds, lam, mu)
    if not ok:
        return None, ("unverified", float(res.fun))
    return {"lam": lam, "mu": mu, "S": str(float(S))}, ("ok", float(res.fun))


def branch(n, N, profs, featlist, z, U, u, want_cert=False):
    eq, ub, bounds = pencil.build_lp(n, N, profs, featlist, z, U, u)
    st = solve(eq, ub, bounds)
    out = {"status": int(st)}
    if st == 2 and want_cert:
        cert, info = farkas(eq, ub, bounds)
        out["certificate"] = cert
        out["farkas"] = info
    return out


def run_case(n, N, rule, real=None, log=print):
    t0 = time.time()
    feats = pencil.enumerate_features(n, N, rule)
    featlist = pencil.prune(feats)
    profs = pencil.profile_list(n, N, rule)
    res = {"n": n, "N": N, "rule": rule, "profiles": len(profs), "features_all": len(feats),
           "features": len(featlist),
           "enum_seconds": round(time.time() - t0, 1)}
    log("n=%d N=%d: %d profiles, %d features (%.1fs)" % (n, N, len(profs), len(featlist),
                                                           time.time() - t0))
    U = pencil.unbalanced_set(profs, N, 1)
    res["U"] = sorted(U)
    base = branch(n, N, profs, featlist, None, None, None)
    res["base"] = base
    feasible = []
    infeasible_branches = []          # (z, u) with None meaning "no constraint"
    if base["status"] == 2:
        infeasible_branches.append((None, None))
    else:
        res["types"] = {}
        for z in sorted(set(pencil.quadric_types(n))):
            bz = branch(n, N, profs, featlist, z, None, None)
            entry = {"status": bz["status"]}
            if bz["status"] == 2:
                infeasible_branches.append((z, None))
            else:
                us = []
                for u in range(consts_H(n) + 1):
                    bu = branch(n, N, profs, featlist, z, U, u)
                    if bu["status"] != 2:
                        us.append((u, bu["status"]))
                    else:
                        infeasible_branches.append((z, u))
                entry["feasible_u"] = us
                if us:
                    feasible.append((z, us))
            res["types"][str(z)] = entry
            log("  z=%s status=%s feasible u=%s" % (z, entry["status"],
                                                    entry.get("feasible_u")))
    res["excluded"] = (base["status"] == 2) or not feasible
    certs = {}
    if res["excluded"]:
        for z, u in infeasible_branches:
            b = branch(n, N, profs, featlist, z, U if u is not None else None, u,
                       want_cert=True)
            key = "base" if z is None else ("%s" % (z,) if u is None else "%s|u=%d" % (z, u))
            certs[key] = b.get("certificate")
    res["all_certified"] = bool(certs) and all(c is not None for c in certs.values())
    res["n_certificates"] = len(certs)
    res["n_uncertified"] = sum(1 for c in certs.values() if c is None)
    if real is not None:
        z, u = real
        br = branch(n, N, profs, featlist, z, U, u)
        res["real_branch"] = {"z": z, "u": u, "status": br["status"]}
    res["seconds"] = round(time.time() - t0, 1)
    log("  -> excluded=%s all_certified=%s (%d certs) [%.1fs]" % (
        res["excluded"], res["all_certified"], len(certs), time.time() - t0))
    return res, certs


def consts_H(n):
    return pencil.consts(n)["H"]


def real_branch_of(K, n, rule="sd2"):
    data = pencil.real_data(K, n, rule)
    N = data["N"]
    profs = pencil.profile_list(n, N, rule)
    z = [0, 0, 0]
    for p, c in data["classes"].items():
        z[pencil.residue(p, N)] += c
    U = pencil.unbalanced_set(profs, N, 1)
    u = sum(c for p, c in data["classes"].items() if p in U)
    return tuple(z), u


LEV = [(1, 2, 0, 0), (1, 0, 1, 0), (0, 1, 1, 0), (1, 2, 1, 0), (0, 0, 2, 0), (2, 1, 2, 0),
       (0, 0, 0, 1), (1, 0, 0, 1), (2, 0, 0, 1), (0, 1, 0, 1), (2, 0, 1, 1), (2, 1, 1, 1),
       (1, 2, 1, 1), (0, 0, 2, 1),
       (0, 0, 0, 2), (2, 0, 0, 2), (1, 1, 0, 2), (2, 1, 0, 2), (0, 2, 0, 2), (2, 2, 0, 2),
       (0, 0, 1, 2), (1, 1, 1, 2), (0, 2, 1, 2), (0, 0, 2, 2), (2, 0, 2, 2), (0, 1, 2, 2),
       (1, 1, 2, 2)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="pencil.json")
    ap.add_argument("--certdir", default=".")
    args = ap.parse_args()
    t0 = time.time()
    out = {"contract": "k5pencil-v1", "started": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                               time.gmtime()), "cases": []}

    def save():
        with open(args.out, "w") as f:
            json.dump(out, f, indent=1, default=str)

    # a real 13-point set in F_3^3 for the control, if random reduction finds one
    k13 = None
    for seed in range(400):
        K = pencil.greedy_kakeya(3, seed)
        if len(K) == 13:
            k13 = K
            break
    out["k13_found"] = k13 is not None
    plan = [(3, 12, None), (3, 13, real_branch_of(k13, 3) if k13 else None),
            (4, 25, None), (4, 26, None), (4, 27, real_branch_of(LEV, 4)),
            (5, 45, None)]
    for rule in ("sd2",):
        for n, N, real in plan:
            res, certs = run_case(n, N, rule, real)
            out["cases"].append(res)
            with open("%s/cert_n%d_N%d_%s.json" % (args.certdir, n, N, rule), "w") as f:
                json.dump(certs, f)
            save()
    # n = 5: if 45 is excluded, push upward; 53 is a control (a 53-point set exists)
    if out["cases"][-1]["excluded"]:
        for N in range(46, 53):
            res, certs = run_case(5, N, "sd2")
            out["cases"].append(res)
            with open("%s/cert_n5_N%d_sd2.json" % (args.certdir, N), "w") as f:
                json.dump(certs, f)
            save()
            if not res["excluded"]:
                break
    res, _ = run_case(5, 53, "sd2")
    res["role"] = "control: a 53-point Kakeya set exists, must be feasible"
    out["cases"].append(res)
    # the class-only rule (T) for the main cases, to see how much SD2 matters
    for n, N in ((4, 26), (5, 45)):
        res, certs = run_case(n, N, "T")
        out["cases"].append(res)
        with open("%s/cert_n%d_N%d_T.json" % (args.certdir, n, N), "w") as f:
            json.dump(certs, f)
        save()
    controls = {}
    for c in out["cases"]:
        if (c["n"], c["N"]) in ((3, 13), (4, 27), (5, 53)):
            controls["%d,%d,%s" % (c["n"], c["N"], c["rule"])] = not c["excluded"]
        if "real_branch" in c:
            controls["real %d,%d" % (c["n"], c["N"])] = c["real_branch"]["status"] == 0
    out["controls_ok"] = controls
    out["seconds"] = round(time.time() - t0, 1)
    save()
    print("controls", controls)


if __name__ == "__main__":
    main()
