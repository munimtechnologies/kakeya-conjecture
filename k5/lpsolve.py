# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""A small phase-1 simplex in pure Python (standard library only) with exact certificates.

Decides whether {v >= 0 : A v = b, G v <= h} has a solution, for integer rows given as
(dict column -> int, rhs).  The simplex runs in floating point, but it never *proves*
anything by itself: when it ends with a positive phase-1 objective, the dual vector of the
final basis is recomputed exactly in rational arithmetic (solving B^T y = c_B with Fractions)
and returned as integer Farkas multipliers (lam for the equality rows, mu >= 0 for the
inequality rows).  The caller then checks them exactly (pencil.check_certificate).  A wrong
floating-point decision can therefore only cost a certificate, never produce a false proof.

Why a phase-1 dual is a Farkas certificate: after flipping rows so that every right-hand side
is >= 0, the phase-1 problem is  min 1.a  s.t.  M v + a = b', v, a >= 0.  At an optimal basis
B the duals y = c_B B^-1 satisfy y.M_j <= 0 for every structural and slack column j and
y.b' = optimum > 0.  With w_i = sign_i y_i (multipliers of the original rows), lam = -w on
equality rows and mu = -w on inequality rows: the slack columns give mu >= 0, the structural
columns give A^T lam + G^T mu >= 0, and b.lam + h.mu = -y.b' < 0.  So no v >= 0 exists.
Extra <= rows (for example variable upper bounds) may be passed; their multipliers are
dropped, which is valid when the checker accounts for those bounds itself.
"""
from fractions import Fraction
from math import gcd

PIV_TOL = 1e-9
PRICE_TOL = 1e-9


def _build(eq, ub, extra, ncols):
    rows = list(eq) + list(ub) + list(extra)
    m = len(rows)
    sign = [1 if rhs >= 0 else -1 for _, rhs in rows]
    b = [sign[i] * rows[i][1] for i in range(m)]
    colr = [[] for _ in range(ncols)]
    colv = [[] for _ in range(ncols)]
    for i, (r, _) in enumerate(rows):
        for j, c in r.items():
            if c:
                colr[j].append(i)
                colv[j].append(sign[i] * c)
    ne = len(eq)
    for i in range(ne, m):                 # one slack per <= row
        colr.append([i])
        colv.append([sign[i]])
    nstruct = len(colr)
    for i in range(m):                     # one artificial per row
        colr.append([i])
        colv.append([1])
    return m, sign, b, colr, colv, nstruct


def _invert(M):
    """Float Gauss-Jordan inverse with partial pivoting, or None if (near) singular."""
    m = len(M)
    A = [[float(x) for x in M[i]] + [1.0 if i == j else 0.0 for j in range(m)]
         for i in range(m)]
    for c in range(m):
        p = max(range(c, m), key=lambda r: abs(A[r][c]))
        if abs(A[p][c]) < 1e-12:
            return None
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        rc = [x / pv for x in A[c]]
        A[c] = rc
        for r in range(m):
            if r != c:
                f = A[r][c]
                if f != 0.0:
                    A[r] = [x - f * y for x, y in zip(A[r], rc)]
    return [row[m:] for row in A]


def _exact_dual(basis, colr, colv, m, nstruct):
    """Solve B^T y = c_B exactly (Fractions).  None if B is singular."""
    M = [[Fraction(0)] * m + [Fraction(1 if basis[k] >= nstruct else 0)] for k in range(m)]
    for k, j in enumerate(basis):
        for i, v in zip(colr[j], colv[j]):
            M[k][i] = Fraction(v)
    for c in range(m):
        p = next((r for r in range(c, m) if M[r][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        rc = [x / pv for x in M[c]]
        M[c] = rc
        for r in range(m):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], rc)]
    return [M[i][m] for i in range(m)]


def _lcm(a, b):
    return a // gcd(a, b) * b


def phase1(eq, ub, ncols, extra_ub=(), max_iter=20000, block=1500, refactor=60):
    """Feasibility of {v >= 0 : eq rows, ub rows, extra_ub rows}.

    Returns a dict with "status" in {"feasible", "infeasible", "unknown"}, "iters", "obj";
    for "infeasible" also integer multipliers "lam" {eq row: int} and "mu" {ub row: int >= 0}
    (extra_ub rows get no multiplier).  The certificate must still be checked exactly."""
    m, sign, b, colr, colv, nstruct = _build(eq, ub, extra_ub, ncols)
    ncol = len(colr)
    cols = [list(zip(colr[j], [float(v) for v in colv[j]])) for j in range(ncol)]
    basis = list(range(nstruct, nstruct + m))
    inb = [False] * ncol
    for j in basis:
        inb[j] = True
    Binv = [[1.0 if i == k else 0.0 for i in range(m)] for k in range(m)]
    bf = [float(v) for v in b]
    xB = list(bf)
    scale = 1.0 + max([abs(v) for v in bf] + [0.0])
    tol_obj = 1e-9 * scale
    best_obj, stall, start = float("inf"), 0, 0
    obj = sum(xB)
    for it in range(1, max_iter + 1):
        if it % refactor == 0:
            Bm = [[0.0] * m for _ in range(m)]
            for k, j in enumerate(basis):
                for i, v in cols[j]:
                    Bm[i][k] = v
            inv = _invert(Bm)
            if inv is None:
                return {"status": "unknown", "why": "singular basis", "iters": it}
            Binv = inv
            xB = [sum(a * c for a, c in zip(row, bf)) for row in Binv]
            xB = [0.0 if (v < 0 and v > -1e-9 * scale) else v for v in xB]
        art = [basis[k] >= nstruct for k in range(m)]
        obj = sum(x for x, a in zip(xB, art) if a)
        if obj <= tol_obj:
            return {"status": "feasible", "iters": it, "obj": obj}
        if obj < best_obj - 1e-12 * scale:
            best_obj, stall = obj, 0
        else:
            stall += 1
        bland = stall > 40
        y = [0.0] * m
        for k in range(m):
            if art[k]:
                rk = Binv[k]
                for i in range(m):
                    y[i] += rk[i]
        q, dq = -1, -PRICE_TOL
        if bland:
            for j in range(nstruct):
                if not inb[j]:
                    d = -sum([y[i] * a for i, a in cols[j]])
                    if d < -PRICE_TOL:
                        q = j
                        break
        else:
            j, scanned, seen = start, 0, 0
            while scanned < nstruct:
                if not inb[j]:
                    d = -sum([y[i] * a for i, a in cols[j]])
                    if d < dq:
                        dq, q = d, j
                scanned += 1
                seen += 1
                j += 1
                if j == nstruct:
                    j = 0
                if seen >= block:
                    if q >= 0:
                        break
                    seen = 0
            start = j
        if q < 0:
            y_exact = _exact_dual(basis, colr, colv, m, nstruct)
            if y_exact is None:
                return {"status": "unknown", "why": "exact basis singular", "iters": it,
                        "obj": obj}
            ne, nu = len(eq), len(ub)
            w = [sign[i] * y_exact[i] for i in range(m)]
            lam = {i: -w[i] for i in range(ne) if w[i] != 0}
            mu = {k: -w[ne + k] for k in range(nu) if w[ne + k] < 0}
            den = 1
            for v in list(lam.values()) + list(mu.values()):
                den = _lcm(den, v.denominator)
            lam = {i: int(v * den) for i, v in lam.items()}
            mu = {k: int(v * den) for k, v in mu.items()}
            return {"status": "infeasible", "iters": it, "obj": obj, "lam": lam, "mu": mu}
        cq = cols[q]
        u = [sum([Binv[k][i] * a for i, a in cq]) for k in range(m)]
        tmin = None
        for k in range(m):
            if u[k] > PIV_TOL:
                t = max(xB[k], 0.0) / u[k]
                if tmin is None or t < tmin:
                    tmin = t
        if tmin is None:
            return {"status": "unknown", "why": "unbounded direction", "iters": it, "obj": obj}
        cand = [k for k in range(m)
                if u[k] > PIV_TOL and max(xB[k], 0.0) / u[k] <= tmin + 1e-10 * (1.0 + tmin)]
        if bland:
            r = min(cand, key=lambda k: basis[k])
        else:
            r = max(cand, key=lambda k: (art[k], u[k]))
        piv = u[r]
        Br = [x / piv for x in Binv[r]]
        Binv[r] = Br
        xr = xB[r] / piv
        xB[r] = xr
        for k in range(m):
            if k != r:
                f = u[k]
                if f != 0.0:
                    Binv[k] = [a - f * c for a, c in zip(Binv[k], Br)]
                    xB[k] -= f * xr
        inb[basis[r]] = False
        basis[r] = q
        inb[q] = True
    return {"status": "unknown", "why": "iteration limit", "iters": max_iter, "obj": obj}
